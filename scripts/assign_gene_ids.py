"""
Give StringTie's novel transcripts a gene identity.

Adapted from Parse_GTF_to_fix_gene_id_issue.py. The logic is unchanged:

  * transcripts that already carry a gene_name get a composite gene id of the
    form "<gene_id>-<gene_name>", and their exon coordinates are recorded
  * transcripts with no gene_name (StringTie's MSTRG.* novel entries) are
    matched against those recorded exons; a transcript sharing an exact exon
    with a named gene inherits that gene's composite id
  * novel transcripts that match nothing are dropped, as in the original

Differences from the standalone script:

  * input and output paths come from Snakemake instead of being hard-coded
  * the exon lookup is a set keyed on chrom/start/end rather than a list that
    is re-scanned per novel transcript. The original was O(novel x exons) and
    becomes very slow on a merged assembly spanning tens of samples; this is
    the same comparison, just indexed.
  * coordinate keys are tuples rather than concatenated strings, so a
    chromosome ending in digits cannot collide with a start coordinate
    (e.g. chr1 + 100 + 200 and chr11 + 00 + 200 both gave "chr1100200")
  * header lines are preserved, as in the original
"""

import os
import sys

_log = open(snakemake.log[0], "w")
_stdout, _stderr = sys.stdout, sys.stderr
sys.stdout = _log
sys.stderr = _log

gtf_in = snakemake.input.merged
gtf_out = snakemake.output.assigned

print(f"Reading {gtf_in}")

header_lines = []
rows = []

with open(gtf_in) as fh:
    for line in fh:
        if line.startswith("#"):
            header_lines.append(line.rstrip("\n\r"))
        else:
            rows.append(line.rstrip("\n\r").split("\t"))

print(f"  {len(header_lines)} header lines, {len(rows)} feature lines")


def attr(attribute_string, key):
    """Value of a GTF attribute, or None when absent."""
    token = f'{key} "'
    if token not in attribute_string:
        return None
    return attribute_string.split(token, 1)[1].split('"', 1)[0]


# ------------------------------------------------------------------
# Pass 1: transcripts that already have a gene_name
# ------------------------------------------------------------------
named_rows = []          # rows rewritten with the composite gene id
exon_to_gene = {}        # (chrom, start, end) -> composite gene id
novel_by_gene = {}       # original gene_id -> rows lacking gene_name

for row in rows:
    attrs = row[8]
    gene_id = attr(attrs, "gene_id")
    if gene_id is None:
        continue

    gene_name = attr(attrs, "gene_name")
    if gene_name is not None:
        new_gene_id = f"{gene_id}-{gene_name}"
        # Replace only the gene_id value, leaving the rest of the attribute
        # string untouched.
        new_attrs = 'gene_id "' + new_gene_id + attrs.split(gene_id, 1)[1]
        new_row = row[:-1] + [new_attrs]
        named_rows.append(new_row)

        if row[2] == "exon":
            exon_to_gene.setdefault((row[0], row[3], row[4]), new_gene_id)
    else:
        novel_by_gene.setdefault(gene_id, []).append(row)

print(f"  {len(named_rows)} rows from genes with a gene_name")
print(f"  {len(exon_to_gene)} distinct annotated exons recorded")
print(f"  {len(novel_by_gene)} novel gene ids without a gene_name")

# ------------------------------------------------------------------
# Pass 2: assign novel transcripts by shared exon
# ------------------------------------------------------------------
transcript_to_gene = {}

for novel_gene, novel_rows in novel_by_gene.items():
    for row in novel_rows:
        transcript_id = attr(row[8], "transcript_id")
        if transcript_id is None or transcript_id in transcript_to_gene:
            continue
        hit = exon_to_gene.get((row[0], row[3], row[4]))
        if hit is not None:
            transcript_to_gene[transcript_id] = hit

print(f"  {len(transcript_to_gene)} novel transcripts matched to a named gene")

assigned_novel_rows = []
unmatched = 0

for novel_rows in novel_by_gene.values():
    for row in novel_rows:
        transcript_id = attr(row[8], "transcript_id")
        if transcript_id in transcript_to_gene:
            rest = ";".join(row[8].split(";")[1:])
            new_attrs = f'gene_id "{transcript_to_gene[transcript_id]}";{rest}'
            assigned_novel_rows.append(row[:-1] + [new_attrs])
        else:
            unmatched += 1

print(f"  {len(assigned_novel_rows)} novel rows reassigned")
print(f"  {unmatched} novel rows dropped (no shared exon with a named gene)")

# ------------------------------------------------------------------
# Write out
# ------------------------------------------------------------------
os.makedirs(os.path.dirname(gtf_out) or ".", exist_ok=True)

with open(gtf_out, "w") as out:
    for line in header_lines:
        out.write(line + "\n")
    for row in named_rows:
        out.write("\t".join(row) + "\n")
    for row in assigned_novel_rows:
        out.write("\t".join(row) + "\n")

print(f"\nWrote {gtf_out}")
print(f"  {len(named_rows) + len(assigned_novel_rows)} feature lines total")

sys.stdout, sys.stderr = _stdout, _stderr
_log.close()
