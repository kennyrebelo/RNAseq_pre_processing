# RNA-seq pre-processing, differential gene expression and splicing quantification pipeline

A Snakemake workflow that takes raw RNA-seq reads (local FASTQs or SRA
accessions) through trimming, alignment, transcript quantification, QC, and
optionally through genome-browser track generation, differential gene
expression, and differential splicing with novel-isoform annotation.

Built for a PBS/Torque HPC cluster, but the workflow itself is
scheduler-agnostic — cluster specifics live entirely in a Snakemake profile.

<p align="center">
  <img src="docs/images/pipeline_subway_compact.svg" alt="Pipeline subway map" width="100%">
</p>

<p align="center">
  <sub>
    <a href="docs/images/pipeline_subway_detailed.svg">Detailed view</a> — every rule as its own station ·
    <a href="docs/images/pipeline_rulegraph.svg">Snakemake rulegraph</a>
  </sub>
</p>

---

## What it does

| Stage | Tools | Always runs? |
|---|---|---|
| SRA download | `fastq-dump` | Only for samples given an SRA accession |
| Read trimming + FastQC | Trim Galore, FastQC | Yes |
| Genome alignment | STAR | Yes |
| Unique-read filtering, sorting, indexing | samtools | Yes |
| Transcript quantification | Salmon | Yes |
| Aggregate QC report | MultiQC | Yes |
| Strandedness inference | RSeQC `infer_experiment.py` | When bigwig or splicing is enabled |
| Strand-split bigwigs + UCSC track hub | pysam, bedtools, `bedGraphToBigWig` | Per sample, via `generate_bw` |
| Differential gene expression | tximport, DESeq2 | Via `dge.enabled` |
| Transcript assembly + differential splicing | StringTie, rMATS | Via `splicing.enabled` |
| Cryptic/annotated exon calling | custom scripts | With splicing |

Single-end and paired-end samples are handled automatically, as are mixtures of
human and mouse samples within a run (except in the splicing branch, which
requires one species — see [Configuration](docs/CONFIGURATION.md)).

## Documentation

- **[Installation](docs/INSTALLATION.md)** — dependencies, conda environments, cluster profile
- **[Configuration](docs/CONFIGURATION.md)** — every tunable option, and what is deliberately fixed
- **[Usage](docs/USAGE.md)** — running, testing, resuming, common commands
- **[Outputs](docs/OUTPUTS.md)** — what gets written where
- **[Troubleshooting](docs/TROUBLESHOOTING.md)** — failure modes seen in practice

## Quick start

```bash
# 1. clone
git clone https://github.com/<user>/<repo>.git
cd <repo>

# 2. environment holding Snakemake itself
#    NOTE: the pipeline requires Snakemake 7.18.1. Installing Snakemake
#    unpinned will give you 8.x or later, whose profile format, cluster
#    submission interface and conda handling are incompatible with this
#    workflow. Use the pinned environment file:
conda env create -f environment_snakemake.yaml
conda activate snakemake

# 3. verify before going further
snakemake --version        # must print 7.18.1
python -c "import pandas"  # must not error

# 4. describe your data
#    edit samples.tsv, experimental_design.tsv, config_RNAseq_pipeline.yaml
#    (the versions shipped here are annotated templates with
#     /path/to/... placeholders and one example of every input type)

# 5. install the cluster profile, if running on PBS/Torque
cp -r profiles/pbs-torque ~/.config/snakemake/
chmod +x ~/.config/snakemake/pbs-torque/{pbs-submit.py,pbs-status.py,pbs-jobscript.sh}
#    then edit the conda path at the top of pbs-jobscript.sh -- see
#    profiles/pbs-torque/README.md

# 6. check the plan without running anything
snakemake --profile pbs-torque --use-conda -n

# 7. run
snakemake --profile pbs-torque --use-conda
```

Not on PBS? Drop `--profile` and give `--cores` instead — the workflow itself
contains no scheduler-specific code. See
[profiles/pbs-torque/README.md](profiles/pbs-torque/README.md) for adapting to
SLURM, SGE or LSF.

The three tool environments (`environment_RNAseq.yaml`,
`environment_DESeq2.yaml`, `environment_rMATS.yaml`) are built automatically by
`--use-conda` on first run. Expect 15–30 minutes for that build. They are
separate from `environment_snakemake.yaml`, which you create yourself in step 2
and which holds no bioinformatics tools.

## Directory structure

Everything below lives in one directory, and the workflow is run from it.
The `rules/` and `scripts/` directory names are referenced by `Snakefile` and
`rules/common.smk`, so they cannot be renamed without editing those two files.

```
<project>/
├── Snakefile                          # configfile, rule all, module includes
├── config_RNAseq_pipeline.yaml        # all tunable settings
├── samples.tsv                        # one row per sample
├── experimental_design.tsv            # group + contrast definitions
│
├── LICENSE                            # MIT
├── CITATION.cff                       # how to cite this pipeline
├── CHANGELOG.md                       # what changed between versions
├── .gitignore
│
├── environment_snakemake.yaml         # Snakemake itself — PINNED, see Quick start
├── environment_RNAseq.yaml            # main toolchain
├── environment_DESeq2.yaml            # R / Bioconductor, DGE only
├── environment_rMATS.yaml             # StringTie + rMATS, splicing only
│
├── rules/
│   ├── common.smk                     # config, sample sheet, helpers  (included first)
│   ├── preprocessing.smk              # SRA download, Trim Galore
│   ├── alignment.smk                  # STAR, unique filtering, Salmon
│   ├── qc.smk                         # MultiQC
│   ├── bigwig.smk                     # strandedness, bigwigs, track hub
│   ├── dge.smk                        # DESeq2
│   └── splicing.smk                   # StringTie + rMATS
│
├── profiles/
│   └── pbs-torque/                    # example PBS/Torque profile
│       ├── config.yaml                # cluster command, job limits, resources
│       ├── pbs-submit.py              # builds and issues qsub
│       ├── pbs-status.py              # reports job state to Snakemake
│       ├── pbs-jobscript.sh           # per-job wrapper  ← must be edited
│       └── README.md
│
├── scripts/
│   ├── split_strands.py
│   ├── make_bigwig.py
│   ├── make_genome_browser_hub.py
│   ├── deseq2_contrast.R
│   ├── splicing_library_info.py
│   ├── assign_gene_ids.py
│   ├── filter_se_events.py
│   └── identify_cryptic_exons.py
│
└── docs/
    ├── INSTALLATION.md
    ├── CONFIGURATION.md
    ├── USAGE.md
    ├── OUTPUTS.md
    ├── TROUBLESHOOTING.md
    └── images/
        ├── pipeline_subway_compact.svg     # the map shown at the top of this file
        ├── pipeline_subway_detailed.svg    # every rule as its own station
        └── pipeline_rulegraph.svg          # raw `snakemake --rulegraph` output
```

`samples.tsv`, `experimental_design.tsv` and `config_RNAseq_pipeline.yaml` ship
as annotated templates: every path is a `/path/to/...` placeholder, and the
sample sheet includes one row of each input type (paired-end, single-end, SRA
accession, and a second species) with `generate_bw` set both ways.

`profiles/` is only needed on a cluster. It is copied to
`~/.config/snakemake/` rather than used in place, since that is where Snakemake
looks up a profile by name.

Generated at runtime (not tracked in git):

```
├── processed/                         # or whatever pre_processing_dir points to
│   ├── trimmed_reads/  star_alignment/  salmon_quant/
│   ├── bigwig/  genome_browser/
│   └── logs/<rule>/<sample>.log
├── QC/multiqc_report.html
├── DGE/
├── splicing/
├── logs/cluster/                      # scheduler stdout/stderr
└── .snakemake/                        # Snakemake state and built conda envs
```

## Design notes

**The Snakemake environment is pinned, and the cluster depends on it.**
`environment_snakemake.yaml` pins Snakemake to 7.18.1. Beyond Snakemake itself,
`pandas` is a hard requirement of this environment rather than a convenience:
`rules/common.smk` parses the sample sheet with it, and `pbs-jobscript.sh`
activates this same environment on every compute node, where each job re-reads
the workflow before running its rule. Anything imported at workflow-parse time
must therefore be present here, not only in the tool environments built by
`--use-conda`. The same coupling means a workflow file edited on the login node
must be visible to the compute nodes — a worker that reads a stale copy will
fail at DAG build with an error that looks nothing like a rule failure.

**Three separate conda environments for the tools.** The R/Bioconductor stack
and the StringTie/rMATS stack are isolated from the main toolchain so that
changing one does not force the solver to re-resolve the others. Each pins
`sysroot_linux-64=2.17` to keep the solver on builds that run against older
cluster glibc, and uses `nodefaults` to keep stale `pkgs/main` and `pkgs/r`
builds out of the solve.

**Library type is inferred, not declared.** RSeQC's `infer_experiment.py`
determines strandedness per sample, which then sets StringTie's `--rf`/`--fr`
flag and rMATS's `--libType`. If samples in a splicing run disagree, the
pipeline stops rather than picking one — an inverted library type silently
flips every junction's strand assignment.

**Checkpoints where output is unknown ahead of time.** SRA downloads may
produce one or two FASTQs, and strandedness is only known after alignment. Both
are Snakemake `checkpoint`s, so the DAG is re-evaluated once they complete.

**Failure is loud.** Every rule writes a `log:`, which survives job failure even
though Snakemake deletes failed outputs. Sample sheet and design file problems
are caught at DAG-build time rather than mid-run on a compute node.

## Citation

If you use this pipeline, please cite the underlying tools: Snakemake,
Trim Galore, Cutadapt, FastQC, STAR, Salmon, samtools, MultiQC, RSeQC,
bedtools, DESeq2, tximport, StringTie, and rMATS.

The subway-map figures in this repository were produced with
[slitherFlow](https://github.com/mikemartinez99/slitherFlow) (Martinez, M.),
which is not a dependency of the workflow itself.

## Development note

The workflow files, helper scripts and repository structure in this project were
drafted with the assistance of AI (Anthropic's Claude), then tested, debugged and
validated against real datasets on the target HPC cluster by the author.
