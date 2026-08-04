"""
Generate normalised bigwig coverage tracks from BAM files.

Adapted from the lab's makeBW.py. Behaviour preserved from the original:

  * coverage comes from `genomeCoverageBed -ibam ... -bg -split`
  * values are scaled to reads-per-scale_factor:
        scaled = raw_coverage * scale_factor / mapped_reads
    with scale_factor = 1e8 reproducing the original "RP10_8" tracks
  * the reverse-strand track is negated so it plots below the axis

Deliberate differences from the original script:

  * THE F/R INVERSION HAS BEEN REMOVED. The original built the forward-strand
    bedgraph from R_*.bam and the reverse-strand bedgraph from F_*.bam. Each
    track is now built from its own BAM, so the F track really does carry the
    reads split_strands.py assigned to the forward bucket. See split_strands.py
    for the strand-assignment logic this now lines up with.
  * no `track type=bedGraph ...` header line is written into the bedgraph.
    The original wrote one and then stripped it again before conversion,
    because bedGraphToBigWig rejects track lines. Skipping it is equivalent
    and avoids a failure mode if scaling is ever disabled.
  * bedgraphs are sorted with `sort -k1,1 -k2,2n` before conversion.
    bedGraphToBigWig requires this ordering and the original relied on the
    BAM header order happening to match.
  * negation is applied to the scaled float rather than to an int()-cast of
    the raw coverage, which avoids discarding sub-integer coverage values.

Intermediate bedgraphs are removed once the bigwigs exist; the strand BAMs are
marked temp() in the Snakefile and cleaned up by Snakemake.
"""

import os
import subprocess
import sys

# `script:` rules don't redirect output automatically. Route both this
# script's own output and the stderr of every subprocess it spawns to the
# log file Snakemake handed us. Note that subprocesses inherit the OS-level
# file descriptor, not sys.stderr, so each subprocess.run() below has to be
# given the handle explicitly.
LOG = open(snakemake.log[0], "w")
sys.stdout = LOG
sys.stderr = LOG

mode = snakemake.params.mode
chrom_sizes = snakemake.params.chrom_sizes
bg2bw = snakemake.params.bedgraph_to_bigwig
scale_factor = float(snakemake.params.scale_factor)
workdir = snakemake.params.workdir

os.makedirs(workdir, exist_ok=True)


def read_count(path):
    with open(path) as fh:
        value = fh.read().strip().split()[0]
    count = float(value)
    if count <= 0:
        raise ValueError(
            f"Mapped-read count read from {path} is {count}; cannot normalise "
            "coverage. Check that the BAM is not empty."
        )
    return count


def genome_coverage(bam, out_bedgraph):
    """Raw per-base coverage bedgraph from a BAM."""
    with open(out_bedgraph, "w") as fh:
        subprocess.run(
            ["genomeCoverageBed", "-ibam", bam, "-g", chrom_sizes, "-bg", "-split"],
            stdout=fh,
            stderr=LOG,
            check=True,
        )


def scale_bedgraph(in_bedgraph, out_bedgraph, tags, negate=False):
    """Normalise coverage values and optionally flip their sign."""
    factor = scale_factor / tags
    with open(in_bedgraph) as fin, open(out_bedgraph, "w") as fout:
        for line in fin:
            if line.startswith("track") or not line.strip():
                continue
            chrom, start, end, value = line.rstrip("\n").split("\t")[:4]
            scaled = float(value) * factor
            if negate:
                scaled = -scaled
            fout.write(f"{chrom}\t{start}\t{end}\t{scaled}\n")


def bedgraph_to_bigwig(bedgraph, bigwig):
    """Sort the bedgraph, then convert it with the UCSC binary."""
    sorted_bg = bedgraph + ".sorted"
    with open(sorted_bg, "w") as fh:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", bedgraph],
            stdout=fh,
            stderr=LOG,
            check=True,
        )
    subprocess.run(
        [bg2bw, sorted_bg, chrom_sizes, bigwig],
        stdout=LOG,
        stderr=LOG,
        check=True,
    )
    os.remove(sorted_bg)


tmp_files = []
tags = read_count(snakemake.input.counts)
print(f"[make_bigwig] mode={mode}, mapped reads={tags:.0f}")

if mode == "stranded":
    fwd_raw = os.path.join(workdir, f"{snakemake.wildcards.sample}_F.bedgraph")
    rev_raw = os.path.join(workdir, f"{snakemake.wildcards.sample}_R.bedgraph")
    fwd_scaled = os.path.join(workdir, f"{snakemake.wildcards.sample}_F_scaled.bedgraph")
    rev_scaled = os.path.join(workdir, f"{snakemake.wildcards.sample}_R_scaled.bedgraph")

    # Each strand's coverage now comes from its own BAM (the original
    # makeBW.py had these swapped).
    genome_coverage(snakemake.input.fwd_bam, fwd_raw)
    genome_coverage(snakemake.input.rev_bam, rev_raw)

    scale_bedgraph(fwd_raw, fwd_scaled, tags, negate=False)
    scale_bedgraph(rev_raw, rev_scaled, tags, negate=True)

    bedgraph_to_bigwig(fwd_scaled, snakemake.output.fwd_bw)
    bedgraph_to_bigwig(rev_scaled, snakemake.output.rev_bw)

    tmp_files = [fwd_raw, rev_raw, fwd_scaled, rev_scaled]

elif mode == "unstranded":
    raw = os.path.join(workdir, f"{snakemake.wildcards.sample}.bedgraph")
    scaled = os.path.join(workdir, f"{snakemake.wildcards.sample}_scaled.bedgraph")

    genome_coverage(snakemake.input.bam, raw)
    scale_bedgraph(raw, scaled, tags, negate=False)
    bedgraph_to_bigwig(scaled, snakemake.output.bw)

    tmp_files = [raw, scaled]

else:
    raise ValueError(f"Unknown mode: {mode!r} (expected 'stranded' or 'unstranded')")

for path in tmp_files:
    if os.path.exists(path):
        os.remove(path)

print("[make_bigwig] done")

LOG.close()
