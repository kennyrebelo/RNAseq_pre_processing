"""
Resolve the three run-wide parameters rMATS needs, from outputs the pipeline
has already produced:

    libtype      --libType   fr-firststrand | fr-secondstrand | fr-unstranded
    read_type    -t          paired | single
    read_length  --readLength

rMATS accepts a single value for each across the whole run, so all samples in
the splicing analysis must agree on library type and on paired/single. A
disagreement is an error rather than something to average over -- an inverted
libType silently flips which transcript strand every junction is assigned to.

LIBRARY TYPE
------------
infer_experiment.py reports the fraction of reads consistent with each of two
mutually exclusive configurations, written as "read number, read strand,
transcript strand":

    1++,1--,2+-,2-+   (++,-- single-end)   read 1 SENSE to transcript
    1+-,1-+,2++,2--   (+-,-+ single-end)   read 1 ANTISENSE to transcript

    sense fraction dominant     -> fr-secondstrand  (StringTie --fr)
    antisense fraction dominant -> fr-firststrand   (StringTie --rf, dUTP)
    neither dominant (~0.5/0.5) -> fr-unstranded    (StringTie: no flag)

"Fraction of reads failed to determine" is a confidence check, not a third
category: a few percent is normal, a large fraction means the transcript BED
does not match the genome the BAM was aligned to.

READ LENGTH
-----------
Taken from Trim Galore's FastQC report rather than by scanning BAMs. FastQC's
"Sequence length" field gives the observed range (e.g. "20-151") and the upper
bound is the read length rMATS wants; --variable-read-length then handles the
shorter trimmed reads. Reading one small zip per sample is far cheaper than
touching the alignments.
"""

import os
import re
import sys
import zipfile

_log = open(snakemake.log[0], "w")
_stdout, _stderr = sys.stdout, sys.stderr
sys.stdout = _log
sys.stderr = _log

samples = list(snakemake.params.samples)
paired_flags = list(snakemake.params.paired)
threshold = float(snakemake.params.strand_threshold)

strand_files = list(snakemake.input.strand)
fastqc_files = list(snakemake.input.fastqc)

# The input lists are built in SPLICING_SAMPLES order, so they line up.
assert len(strand_files) == len(samples)
assert len(fastqc_files) == len(samples)


def parse_infer_experiment(path):
    """
    Returns (libtype, sense_fraction, antisense_fraction, failed_fraction).

    RSeQC prints the sense group first and the antisense group second, so the
    two 'explained by' lines are read in order rather than by matching their
    strand notation, which differs between single- and paired-end output.
    """
    fractions = []
    failed = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("Fraction of reads failed to determine"):
                failed = float(line.rsplit(":", 1)[1])
            elif line.startswith('Fraction of reads explained by "'):
                fractions.append(float(line.rsplit(":", 1)[1]))

    if len(fractions) < 2:
        raise ValueError(
            f"Could not parse two strandedness fractions from {path}. Check "
            "that infer_experiment.py ran successfully and that the transcript "
            "BED matches the genome the BAM was aligned to."
        )

    sense, antisense = fractions[0], fractions[1]
    if sense >= threshold:
        libtype = "fr-secondstrand"
    elif antisense >= threshold:
        libtype = "fr-firststrand"
    else:
        libtype = "fr-unstranded"
    return libtype, sense, antisense, failed


def parse_max_read_length(zip_path):
    """Upper bound of FastQC's 'Sequence length' field."""
    with zipfile.ZipFile(zip_path) as zf:
        data_members = [n for n in zf.namelist() if n.endswith("fastqc_data.txt")]
        if not data_members:
            raise ValueError(f"No fastqc_data.txt inside {zip_path}")
        with zf.open(data_members[0]) as fh:
            for raw in fh:
                line = raw.decode("utf-8", errors="replace").strip()
                if line.startswith("Sequence length"):
                    value = line.split("\t", 1)[1].strip()
                    # Either "151" or "20-151"
                    numbers = [int(n) for n in re.findall(r"\d+", value)]
                    if not numbers:
                        raise ValueError(
                            f"Could not read a read length from '{line}' in {zip_path}"
                        )
                    return max(numbers)
    raise ValueError(f"No 'Sequence length' line found in {zip_path}")


libtypes = {}
read_lengths = {}

print("Per-sample inference\n" + "-" * 70)
for sample, strand_path, fastqc_path, is_paired in zip(
    samples, strand_files, fastqc_files, paired_flags
):
    libtype, sense, antisense, failed = parse_infer_experiment(strand_path)
    max_len = parse_max_read_length(fastqc_path)

    libtypes[sample] = libtype
    read_lengths[sample] = max_len

    print(
        f"{sample}\n"
        f"    sense fraction      : {sense}\n"
        f"    antisense fraction  : {antisense}\n"
        f"    failed to determine : {failed}\n"
        f"    -> libType          : {libtype}\n"
        f"    max read length     : {max_len}\n"
        f"    layout              : {'paired' if is_paired else 'single'}"
    )
    if failed is not None and failed > 0.30:
        print(
            f"    WARNING: {failed:.1%} of reads could not be assigned. That "
            "usually means the transcript BED does not match the genome build "
            "used for alignment; treat the libType call as unreliable."
        )

print("-" * 70)

# --- all samples must agree on library type ---
distinct_libtypes = sorted(set(libtypes.values()))
if len(distinct_libtypes) > 1:
    detail = "\n".join(f"    {s}: {t}" for s, t in libtypes.items())
    raise ValueError(
        "Samples in the splicing analysis disagree on library type:\n"
        f"{detail}\n"
        "rMATS takes one --libType for the whole run. Either split the run by "
        "library type, or check whether a borderline sample is being "
        "misclassified (splicing.strand_threshold controls how confident the "
        "call has to be)."
    )
libtype = distinct_libtypes[0]

# --- all samples must agree on layout ---
distinct_layouts = sorted({"paired" if p else "single" for p in paired_flags})
if len(distinct_layouts) > 1:
    detail = "\n".join(
        f"    {s}: {'paired' if p else 'single'}" for s, p in zip(samples, paired_flags)
    )
    raise ValueError(
        "Samples in the splicing analysis mix paired- and single-end data:\n"
        f"{detail}\n"
        "rMATS takes one -t for the whole run, so these must be analysed "
        "separately."
    )
read_type = distinct_layouts[0]

# rMATS wants the longest read in the run; --variable-read-length covers the rest.
read_length = max(read_lengths.values())

stringtie_flag = {
    "fr-firststrand": "--rf",
    "fr-secondstrand": "--fr",
    "fr-unstranded": "",
}[libtype]

print("\nRun-wide values")
print(f"    libtype         = {libtype}")
print(f"    read_type       = {read_type}")
print(f"    read_length     = {read_length}")
print(f"    stringtie_flag  = {stringtie_flag or '(none)'}")

os.makedirs(os.path.dirname(snakemake.output.info) or ".", exist_ok=True)
with open(snakemake.output.info, "w") as fh:
    fh.write(f"libtype={libtype}\n")
    fh.write(f"read_type={read_type}\n")
    fh.write(f"read_length={read_length}\n")
    fh.write(f"stringtie_flag={stringtie_flag}\n")
    fh.write(f"n_samples={len(samples)}\n")
    for sample in samples:
        fh.write(
            f"# {sample}\tlibtype={libtypes[sample]}\t"
            f"max_read_length={read_lengths[sample]}\n"
        )

print(f"\nWrote {snakemake.output.info}")

# Restore the real streams before closing, so anything Snakemake writes after
# this script returns does not hit a closed file.
sys.stdout, sys.stderr = _stdout, _stderr
_log.close()
