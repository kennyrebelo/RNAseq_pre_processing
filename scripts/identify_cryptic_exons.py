"""
Label each filtered skipped-exon event as Annotated or Cryptic.

Adapted from Parse_rMATS_to_identify_cryptic_isoforms.py. An exon is called
Annotated when its coordinates match an exon in the reference GTF exactly, and
Cryptic otherwise.

The +1 on the start coordinate is retained from the original and is the point
of the comparison: rMATS reports exonStart_0base in BED-style half-open
coordinates while GTF is 1-based inclusive, so the start must be shifted before
the two can be compared.

Differences from the standalone script:

  * paths come from Snakemake; the GTF is the reference annotation for the
    sample's species, taken from the config
  * annotated exons are held in a set keyed on (chrom, start, end) rather than
    a list of concatenated strings. The original did a linear scan of every
    annotated exon per event, and concatenating chrom+start+end can collide
    (chr1|100|200 and chr11|00|200 both became "chr1100200")
  * the header row is carried through and given an "exon_status" column, so
    the output is still readable as a table
  * the original required the GTF's header lines to be deleted by hand; lines
    beginning with # are skipped here
"""

import os
import sys

_log = open(snakemake.log[0], "w")
_stdout, _stderr = sys.stdout, sys.stderr
sys.stdout = _log
sys.stderr = _log

events_path = snakemake.input.sig
gtf_path = snakemake.params.ref_gtf
out_path = snakemake.output.ided

print(f"Events     : {events_path}")
print(f"Annotation : {gtf_path}")
print(f"Output     : {out_path}\n")

# ------------------------------------------------------------------
# Annotated exon coordinates
# ------------------------------------------------------------------
annotated = set()

with open(gtf_path) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n\r").split("\t")
        if len(parts) < 5 or parts[2] != "exon":
            continue
        annotated.add((parts[0], parts[3], parts[4]))

print(f"Annotated exons in GTF: {len(annotated)}")

# ------------------------------------------------------------------
# Classify each event
# ------------------------------------------------------------------
os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

n_annotated = 0
n_cryptic = 0

with open(events_path) as fin, open(out_path, "w") as fout:
    for line in fin:
        parsed = line.rstrip("\n\r").split("\t")
        if len(parsed) < 3:
            continue

        # Header row from filter_se_events.py
        if parsed[1] in ("exonStart_0base", "longExonStart_0base"):
            fout.write("\t".join(parsed + ["exon_status"]) + "\n")
            continue

        chrom = parsed[0]
        try:
            # rMATS start is 0-based half-open; GTF is 1-based inclusive.
            start = str(int(parsed[1]) + 1)
        except ValueError:
            # Any other non-numeric first data column is a header variant.
            fout.write("\t".join(parsed + ["exon_status"]) + "\n")
            continue
        end = parsed[2]

        if (chrom, start, end) in annotated:
            status = "Annotated"
            n_annotated += 1
        else:
            status = "Cryptic"
            n_cryptic += 1

        fout.write("\t".join(parsed + [status]) + "\n")

print(f"\nAnnotated events: {n_annotated}")
print(f"Cryptic events  : {n_cryptic}")

sys.stdout, sys.stderr = _stdout, _stderr
_log.close()
