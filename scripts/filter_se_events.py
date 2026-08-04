"""
Keep only well-supported, significant skipped-exon events from rMATS.

Adapted from Filter_rows_SE.py. An event survives when all three hold:

  * every replicate in BOTH samples has at least `min_reads` supporting reads,
    where support is the mean of inclusion and skipping counts for that
    replicate
  * FDR < `fdr`
  * |IncLevelDifference| >= `min_delta_psi`

Output columns, as in the original: chr, exonStart_0base, exonEnd, GeneID,
geneSymbol, strand, FDR, IncLevel1, IncLevel2, IncLevelDifference.

Differences from the standalone script:

  * paths and thresholds come from Snakemake rather than being hard-coded
  * read counts are parsed as float, not int. rMATS writes non-integer values
    in some columns, and int() on those raised ValueError mid-file.
  * replicates whose counts are "NA" (rMATS emits these when a replicate has
    no coverage at the event) make the event fail the read filter instead of
    crashing the parse
"""

import os
import sys

_log = open(snakemake.log[0], "w")
_stdout, _stderr = sys.stdout, sys.stderr
sys.stdout = _log
sys.stderr = _log

in_path = snakemake.input.se
out_path = snakemake.output.sig

min_reads = float(snakemake.params.min_reads)
max_fdr = float(snakemake.params.fdr)
min_delta_psi = float(snakemake.params.min_delta_psi)

print(f"Input : {in_path}")
print(f"Output: {out_path}")
print(f"Filters: min_reads >= {min_reads}, FDR < {max_fdr}, "
      f"|dPSI| >= {min_delta_psi}\n")


def parse_counts(field):
    """
    Comma-separated per-replicate counts. Returns None when any replicate is
    missing, which makes the event fail the read-support filter.
    """
    values = []
    for token in field.split(","):
        token = token.strip()
        if token in ("", "NA"):
            return None
        try:
            values.append(float(token))
        except ValueError:
            return None
    return values or None


os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

n_total = 0
n_kept = 0
n_skipped_unparsed = 0

with open(in_path) as fin, open(out_path, "w") as fout:
    for line in fin:
        event = line.rstrip("\n\r").split("\t")
        if len(event) < 23:
            continue

        chrom = event[3]
        exon_start = event[5]
        exon_end = event[6]
        gene_id = event[1]
        gene_symbol = event[2]
        strand = event[4]
        fdr = event[19]
        psi_1 = event[20]
        psi_2 = event[21]
        dlt_psi = event[22]

        out_fields = [chrom, exon_start, exon_end, gene_id, gene_symbol,
                      strand, fdr, psi_1, psi_2, dlt_psi]

        if event[0] == "ID":
            fout.write("\t".join(out_fields) + "\n")
            continue

        n_total += 1

        s1_inc = parse_counts(event[12])
        s1_skp = parse_counts(event[13])
        s2_inc = parse_counts(event[14])
        s2_skp = parse_counts(event[15])

        if None in (s1_inc, s1_skp, s2_inc, s2_skp):
            n_skipped_unparsed += 1
            continue

        avg_s1 = [(i + s) / 2 for i, s in zip(s1_inc, s1_skp)]
        avg_s2 = [(i + s) / 2 for i, s in zip(s2_inc, s2_skp)]

        try:
            fdr_value = float(fdr)
            delta_value = abs(float(dlt_psi))
        except ValueError:
            n_skipped_unparsed += 1
            continue

        if (
            min(avg_s1) >= min_reads
            and min(avg_s2) >= min_reads
            and fdr_value < max_fdr
            and delta_value >= min_delta_psi
        ):
            fout.write("\t".join(out_fields) + "\n")
            n_kept += 1

print(f"Events read   : {n_total}")
print(f"Events kept   : {n_kept}")
print(f"Unparseable   : {n_skipped_unparsed} (missing or non-numeric values)")

sys.stdout, sys.stderr = _stdout, _stderr
_log.close()
