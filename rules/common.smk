import os
import re
import glob
import pandas as pd

# ------------------------------------------------------------------
# Paths inside the workflow itself
#
# Snakemake resolves `script:` and `conda:` paths relative to the FILE THE
# RULE IS DEFINED IN, not relative to the main Snakefile. Once rules live in
# rules/*.smk, a bare "scripts/foo.py" would resolve to rules/scripts/foo.py
# and a bare "environment_RNAseq.yaml" to rules/environment_RNAseq.yaml.
# Resolving both against workflow.basedir (the directory holding the main
# Snakefile) keeps them correct no matter which module a rule sits in, and
# also keeps working if the pipeline is invoked from another directory with
# --snakefile / --directory.
# ------------------------------------------------------------------
WORKFLOW_DIR = workflow.basedir
SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "scripts")


def workflow_path(path):
    """Resolve a workflow-relative path (e.g. a conda env) to an absolute one."""
    return path if os.path.isabs(path) else os.path.join(WORKFLOW_DIR, path)


# ------------------------------------------------------------------
# Sample sheet parsing
# ------------------------------------------------------------------
samples_df = pd.read_csv(config["samples_tsv"], sep="\t", dtype=str).fillna("")
samples_df["read_pair1"] = samples_df["read_pair1"].str.strip()
samples_df["read_pair2"] = samples_df["read_pair2"].str.strip()

# Strip trailing slashes so that preproc_dir(sample) always matches the
# `pre_processing_dir` wildcard value Snakemake derives from an output
# path. Without this, checkpoints.<name>.get(pre_processing_dir=...)
# silently fails to match for any row written with a trailing "/".
samples_df["pre_processing_dir"] = (
    samples_df["pre_processing_dir"].str.strip().str.rstrip("/")
)

# Optional column: samples without it simply skip the bigwig branch.
if "generate_bw" not in samples_df.columns:
    samples_df["generate_bw"] = "no"
samples_df["generate_bw"] = samples_df["generate_bw"].str.strip().str.lower()

samples_df = samples_df.set_index("sample_name", drop=False)

SAMPLES = list(samples_df.index)
CONDA_ENV = workflow_path(config["conda_env"])

# ------------------------------------------------------------------
# BigWig / genome-browser settings
# ------------------------------------------------------------------
BW_CFG = config.get("bigwig", {})
BEDGRAPH_TO_BIGWIG = BW_CFG.get("bedgraph_to_bigwig", "bedGraphToBigWig")
HUB_EMAIL = BW_CFG.get("hub_email", "")
SCALE_FACTOR = float(BW_CFG.get("scale_factor", 1e8))
SCALE_SUFFIX = BW_CFG.get("scale_suffix", "RP10_8")
STRAND_THRESHOLD = float(BW_CFG.get("strand_threshold", 0.75))

# By default Snakemake wildcards match ANY characters, including "/".
# Without this constraint, the `sample` wildcard can accidentally absorb
# extra path components when Snakemake matches on-disk files against rule
# output patterns (e.g. a downloaded fastq sitting inside the download_sra
# checkpoint's own output directory) -- since sample names never contain
# "/", this constraint prevents that.
wildcard_constraints:
    sample="[^/]+"

SRA_PATTERN = re.compile(r"^[SED]RR\d+$")

# Validate species column up front so bad values fail fast, before any
# jobs are dispatched, rather than surfacing as a KeyError mid-run.
VALID_SPECIES = set(config["references"].keys())
for _sample in samples_df.index:
    _sp = samples_df.loc[_sample, "species"].strip().lower()
    if _sp not in VALID_SPECIES:
        raise ValueError(
            f"Sample '{_sample}' has species '{_sp}' which is not one of "
            f"the species defined under 'references' in the config: {sorted(VALID_SPECIES)}"
        )
samples_df["species"] = samples_df["species"].str.strip().str.lower()


def species(sample):
    return samples_df.loc[sample, "species"]


def star_genome_dir(sample):
    return config["references"][species(sample)]["star_genome_dir"]


def salmon_index_dir(sample):
    return config["references"][species(sample)]["salmon_index_dir"]


def transcripts_bed(sample):
    return config["references"][species(sample)]["transcripts_bed"]


def chrom_sizes(sample):
    return config["references"][species(sample)]["chrom_sizes"]


def genome_build(sample):
    return config["references"][species(sample)]["genome_build"]


def wants_bigwig(sample):
    return samples_df.loc[sample, "generate_bw"] in ("yes", "y", "true", "1")


# Samples that opted in to the bigwig / genome-browser branch. This is read
# straight from the sample sheet, so it is known at parse time (no checkpoint
# resolution needed) and can be used to build the final target list.
BW_SAMPLES = [s for s in samples_df.index if wants_bigwig(s)]

# Fail fast if a bigwig sample's species lacks the extra reference keys,
# rather than dying halfway through on a compute node.
for _sample in BW_SAMPLES:
    _refs = config["references"][samples_df.loc[_sample, "species"]]
    for _key in ("transcripts_bed", "chrom_sizes", "genome_build"):
        if _key not in _refs:
            raise ValueError(
                f"Sample '{_sample}' has generate_bw=yes but species "
                f"'{samples_df.loc[_sample, 'species']}' is missing "
                f"'{_key}' under 'references' in the config."
            )


# ------------------------------------------------------------------
# Analyses driven by experimental_design.tsv
#
# Both the DESeq2 (dge) and rMATS/StringTie (splicing) branches read the same
# experimental_design.tsv. It is only required when at least one of them is
# enabled -- with both disabled the pipeline parses and runs without it.
# ------------------------------------------------------------------
DGE_CFG = config.get("dge", {})
DGE_ENABLED = bool(DGE_CFG.get("enabled", False))
DGE_OUTDIR = DGE_CFG.get("output_dir", "DGE")
DGE_CONDA_ENV = workflow_path(DGE_CFG.get("conda_env", CONDA_ENV))
DGE_COUNTS_FROM_ABUNDANCE = DGE_CFG.get("counts_from_abundance", "lengthScaledTPM")

SPL_CFG = config.get("splicing", {})
SPLICING_ENABLED = bool(SPL_CFG.get("enabled", False))
SPLICING_OUTDIR = SPL_CFG.get("output_dir", "splicing")
SPLICING_CONDA_ENV = workflow_path(SPL_CFG.get("conda_env", CONDA_ENV))

DESIGN_REQUIRED = DGE_ENABLED or SPLICING_ENABLED

design_df = None
CONTRASTS = []
DESIGN_GROUPS = {}

if DESIGN_REQUIRED:
    _design_path = (
        DGE_CFG.get("experimental_design")
        or SPL_CFG.get("experimental_design")
        or "experimental_design.tsv"
    )
    if not os.path.exists(_design_path):
        _who = " and ".join(
            n for n, on in (("dge", DGE_ENABLED), ("splicing", SPLICING_ENABLED)) if on
        )
        raise ValueError(
            f"{_who} enabled but the experimental design file "
            f"'{_design_path}' does not exist. Either create it, or set "
            "dge.enabled and splicing.enabled to false in the config."
        )

    design_df = pd.read_csv(_design_path, sep="\t", dtype=str).fillna("")
    for _col in ("sample_name", "group", "contrast_A", "contrast_B"):
        if _col not in design_df.columns:
            raise ValueError(
                f"'{_design_path}' is missing the required column '{_col}'. "
                "Expected columns: sample_name, group, contrast_A, contrast_B."
            )
        design_df[_col] = design_df[_col].str.strip()

    design_df = design_df[design_df["sample_name"] != ""]
    design_df = design_df.set_index("sample_name", drop=False)

    # Every sample in the design must be one the pipeline actually processes,
    # otherwise there would be no quant.sf / BAM to analyse.
    _unknown = [s for s in design_df.index if s not in samples_df.index]
    if _unknown:
        raise ValueError(
            f"These samples appear in '{_design_path}' but not in "
            f"{config['samples_tsv']}: {_unknown}"
        )

    for _sample in design_df.index:
        _grp = design_df.loc[_sample, "group"]
        if not _grp:
            raise ValueError(
                f"Sample '{_sample}' in '{_design_path}' has no group assigned."
            )
        DESIGN_GROUPS.setdefault(_grp, []).append(_sample)

    # Contrast pairs live on whichever rows have both columns filled in;
    # the remaining rows leave them blank.
    for _sample in design_df.index:
        _a = design_df.loc[_sample, "contrast_A"]
        _b = design_df.loc[_sample, "contrast_B"]
        if _a and _b and (_a, _b) not in CONTRASTS:
            CONTRASTS.append((_a, _b))

    for _a, _b in CONTRASTS:
        for _side, _grp in (("contrast_A", _a), ("contrast_B", _b)):
            if _grp not in DESIGN_GROUPS:
                raise ValueError(
                    f"{_side} '{_grp}' in '{_design_path}' does not match any "
                    f"value in the 'group' column. Known groups: "
                    f"{sorted(DESIGN_GROUPS)}"
                )
            if DGE_ENABLED and len(DESIGN_GROUPS[_grp]) < 2:
                raise ValueError(
                    f"Group '{_grp}' has only {len(DESIGN_GROUPS[_grp])} "
                    "sample(s). DESeq2 cannot estimate dispersion without at "
                    "least two replicates per group in a contrast."
                )
        if _a == _b:
            raise ValueError(
                f"Contrast '{_a}' vs '{_b}' compares a group against itself."
            )

    if not CONTRASTS:
        raise ValueError(
            f"'{_design_path}' defines no contrast pairs. Fill in contrast_A "
            "and contrast_B on at least one row."
        )


# ------------------------------------------------------------------
# Splicing (StringTie + rMATS) sample set
#
# Only samples belonging to a group that appears in a contrast take part.
# A sample can sit in the design file without being in any contrast; it is
# then excluded here.
# ------------------------------------------------------------------
SPLICING_GROUPS = []
SPLICING_SAMPLES = []

if SPLICING_ENABLED:
    for _a, _b in CONTRASTS:
        for _grp in (_a, _b):
            if _grp not in SPLICING_GROUPS:
                SPLICING_GROUPS.append(_grp)
    for _grp in SPLICING_GROUPS:
        for _sample in DESIGN_GROUPS[_grp]:
            if _sample not in SPLICING_SAMPLES:
                SPLICING_SAMPLES.append(_sample)

    # rMATS takes a single --gtf, so every sample in the run has to share one
    # annotation. Mixed species would silently compare against the wrong one.
    _spl_species = {samples_df.loc[s, "species"] for s in SPLICING_SAMPLES}
    if len(_spl_species) > 1:
        raise ValueError(
            "Splicing analysis spans more than one species "
            f"({sorted(_spl_species)}). rMATS uses a single merged annotation, "
            "so run each species as a separate pipeline."
        )
    SPLICING_SPECIES = _spl_species.pop()

    for _key in ("gtf",):
        if _key not in config["references"][SPLICING_SPECIES]:
            raise ValueError(
                f"splicing.enabled is true but species '{SPLICING_SPECIES}' is "
                f"missing '{_key}' under 'references' in the config."
            )

    # Used to name the merged assembly, the rMATS prep list/output and the
    # shared rMATS tmp directory.
    SPLICING_TAG = "all"
else:
    SPLICING_SPECIES = None
    SPLICING_TAG = ""


def contrast_samples(contrast_a, contrast_b):
    """Samples belonging to either group of a contrast, reference group first."""
    ref = [s for s in design_df.index if design_df.loc[s, "group"] == contrast_a]
    trt = [s for s in design_df.index if design_df.loc[s, "group"] == contrast_b]
    return ref + trt


def contrast_groups(contrast_a, contrast_b):
    """Group labels aligned 1:1 with contrast_samples()."""
    return [
        design_df.loc[s, "group"]
        for s in contrast_samples(contrast_a, contrast_b)
    ]


def group_samples(group):
    """Samples belonging to a single design group."""
    return list(DESIGN_GROUPS[group])


def is_sra(sample):
    return bool(SRA_PATTERN.match(samples_df.loc[sample, "read_pair1"]))


def raw_dir(sample):
    return samples_df.loc[sample, "raw_data_dir"]


def preproc_dir(sample):
    return samples_df.loc[sample, "pre_processing_dir"]


def sub_dir(sample, name):
    return os.path.join(preproc_dir(sample), name)


def local_is_paired(sample):
    """Pairing for non-SRA samples is known directly from the sample sheet."""
    return samples_df.loc[sample, "read_pair2"] != ""


def get_raw_fastqs(sample):
    """
    Returns a list of raw fastq paths for a sample.
    For SRA samples this resolves the download checkpoint first, so it
    can only be called from within an input function (never at parse time).
    """
    if is_sra(sample):
        ckpt_dir = checkpoints.download_sra.get(
            sample=sample, pre_processing_dir=preproc_dir(sample)
        ).output[0]
        files = sorted(glob.glob(os.path.join(ckpt_dir, "*.fastq.gz")))
        if len(files) not in (1, 2):
            raise ValueError(
                f"Unexpected number of fastqs downloaded for {sample}: {files}"
            )
        return files
    else:
        r1 = os.path.join(raw_dir(sample), samples_df.loc[sample, "read_pair1"])
        r2 = samples_df.loc[sample, "read_pair2"]
        if r2:
            return [r1, os.path.join(raw_dir(sample), r2)]
        return [r1]


def sample_is_paired(sample):
    """
    Resolves pairing for ANY sample (local or SRA). Only call from within
    an input function, since it may need to resolve the SRA checkpoint.
    """
    if is_sra(sample):
        return len(get_raw_fastqs(sample)) == 2
    return local_is_paired(sample)


# ------------------------------------------------------------------
# Final targets
# ------------------------------------------------------------------
def aggregate_final_outputs(wildcards):
    targets = []
    for sample in SAMPLES:
        targets.append(
            os.path.join(sub_dir(sample, "star_alignment"),
                         f"{sample}_Aligned.out.unique.sorted.bam.bai")
        )
        targets.append(
            os.path.join(sub_dir(sample, "salmon_quant"), sample, "quant.sf")
        )
    for sample in BW_SAMPLES:
        targets.append(
            os.path.join(sub_dir(sample, "genome_browser"), sample)
        )
    # CONTRASTS is populated whenever the design file is read, which now
    # happens for splicing too -- so the DESeq2 targets need their own guard.
    if DGE_ENABLED:
        for contrast_a, contrast_b in CONTRASTS:
            targets.append(
                os.path.join(DGE_OUTDIR,
                             f"DESeq2_salmon_diff_{contrast_a}_vs_{contrast_b}.tab")
            )
            targets.append(
                os.path.join(DGE_OUTDIR,
                             f"DESeq2_salmon_PCA_{contrast_a}_vs_{contrast_b}.pdf")
            )
            targets.append(
                os.path.join(DGE_OUTDIR,
                             f"DESeq2_salmon_PCA_{contrast_a}_vs_{contrast_b}_labels.pdf")
            )
    if SPLICING_ENABLED:
        for contrast_a, contrast_b in CONTRASTS:
            targets.append(
                os.path.join(SPLICING_OUTDIR, "rMATs_final_output",
                             f"SE.MATS.JC_sig_IDed_{contrast_a}_vs_{contrast_b}.txt")
            )
    targets.append("QC/multiqc_report.html")
    return targets
