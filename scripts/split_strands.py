"""
Split a coordinate-sorted, indexed BAM into forward- and reverse-strand SAM
files.

The read-assignment logic is taken verbatim from the reference snippet used in
the lab's existing genome-browser workflow.

ORIENTATION: this assignment assumes a reverse-stranded (dUTP-style) library,
where read 1 maps antisense to the transcript. A read-1 alignment on the minus
strand therefore comes from a plus-strand transcript and is written to
`forward`; likewise a single-end reverse alignment is written to `forward`.
Read 2, which maps sense, is assigned directly.

The downstream make_bigwig.py step used to invert this again (building the F
bedgraph from R_*.bam). That inversion has been removed, so the F track is now
built from the reads assigned to `forward` here. If your libraries are
forward-stranded rather than dUTP, flip the four branches below -- not the
bedgraph step -- so the two stay consistent.
"""

import sys

import pysam

# `script:` rules don't redirect output automatically the way a shell rule
# does -- the log path is just handed to us. Send everything there.
_log = open(snakemake.log[0], "w")
sys.stdout = _log
sys.stderr = _log

infile = pysam.AlignmentFile(snakemake.input.bam, "rb")
forward = pysam.AlignmentFile(snakemake.output.fwd, "w", template=infile)
reverse = pysam.AlignmentFile(snakemake.output.rev, "w", template=infile)

n_fwd = 0
n_rev = 0

for read in infile.fetch(until_eof=True):
    if read.is_paired:
        if read.is_read1 and read.is_reverse:
            forward.write(read)
            n_fwd += 1
        elif read.is_read1:
            reverse.write(read)
            n_rev += 1
        elif read.is_read2 and read.is_reverse:
            reverse.write(read)
            n_rev += 1
        else:
            forward.write(read)
            n_fwd += 1
    else:
        if read.is_reverse:
            forward.write(read)
            n_fwd += 1
        else:
            reverse.write(read)
            n_rev += 1

forward.close()
reverse.close()
infile.close()

print(f"[split_strands] wrote {n_fwd} reads to {snakemake.output.fwd}")
print(f"[split_strands] wrote {n_rev} reads to {snakemake.output.rev}")

_log.close()
