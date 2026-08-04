"""
Build a UCSC track-hub folder for one sample.

Adapted from create_genome_browser_folders.py. Produces:

    <hub_dir>/
        hub.txt
        genomes.txt
        <genome>/
            trackDb.txt
            <sample>_F_<suffix>.bw     (stranded)
            <sample>_R_<suffix>.bw
          or
            <sample>_unstranded_<suffix>.bw   (unstranded)

Deliberate differences from the original script:

  * the sample name and genome build are passed in explicitly rather than
    parsed out of the bigwig filename, so sample names containing "_F_" or
    "_R_" cannot be truncated
  * bigwigs are COPIED into the hub folder rather than moved. Moving would
    delete the rule's own declared inputs, which makes Snakemake regenerate
    them on every subsequent run.
"""

import os
import shutil
import sys

# `script:` rules don't redirect output automatically -- send it to the log.
_log = open(snakemake.log[0], "w")
sys.stdout = _log
sys.stderr = _log

hub_dir = snakemake.output.hub_dir
sample = snakemake.params.sample
genome = snakemake.params.genome
email = snakemake.params.email

bw_files = list(snakemake.input.bw)
stranded = len(bw_files) > 1

genome_dir = os.path.join(hub_dir, genome)
os.makedirs(genome_dir, exist_ok=True)

# ---------------- hub.txt ----------------
with open(os.path.join(hub_dir, "hub.txt"), "w") as fh:
    fh.write(f"hub {sample}\n")
    fh.write(f"shortLabel {sample}\n")
    fh.write(f"longLabel {sample}\n")
    fh.write("genomesFile genomes.txt\n")
    fh.write(f"email {email}\n")

# ---------------- genomes.txt ----------------
with open(os.path.join(hub_dir, "genomes.txt"), "w") as fh:
    fh.write(f"genome {genome}\n")
    fh.write(f"trackDb {genome}/trackDb.txt\n")

# ---------------- trackDb.txt ----------------
with open(os.path.join(genome_dir, "trackDb.txt"), "w") as fh:
    fh.write(f"track {sample}\n")
    fh.write("container multiWig\n")
    fh.write(f"shortLabel {sample}\n")
    fh.write(f"longLabel {sample}\n")
    fh.write("type bigWig 0 60000\n")
    fh.write("viewLimits -10000:10000\n")
    fh.write("visibility full\n")
    fh.write("maxHeightPixels 160:120:11\n")
    fh.write("aggregate solidOverlay\n")
    fh.write("showSubtrackColorOnUi on\n")
    fh.write("windowingFunction maximum\n")
    fh.write("priority 1.4\n")
    fh.write("configurable on\n")
    fh.write("autoScale on\n")

    if stranded:
        fwd_bw = next(p for p in bw_files if "_F_" in os.path.basename(p))
        rev_bw = next(p for p in bw_files if "_R_" in os.path.basename(p))

        fh.write("\n")
        fh.write(f"track {sample}_F\n")
        fh.write(f"bigDataUrl {os.path.basename(fwd_bw)}\n")
        fh.write(f"shortLabel {sample}\n")
        fh.write(f"longLabel {sample}\n")
        fh.write(f"parent {sample}\n")
        fh.write("type bigWig\n")
        fh.write("color 0,0,255\n\n")

        fh.write(f"track {sample}_R\n")
        fh.write(f"bigDataUrl {os.path.basename(rev_bw)}\n")
        fh.write(f"shortLabel {sample}\n")
        fh.write(f"longLabel {sample}\n")
        fh.write(f"parent {sample}\n")
        fh.write("type bigWig\n")
        fh.write("color 255,0,0\n")
    else:
        fh.write(f"bigDataUrl {os.path.basename(bw_files[0])}\n")
        fh.write("type bigWig\n")
        fh.write("color 0,0,0\n")

# ---------------- copy bigwigs in ----------------
for bw in bw_files:
    shutil.copy2(bw, os.path.join(genome_dir, os.path.basename(bw)))

print(f"[genome_browser_hub] built hub for {sample} ({genome}) at {hub_dir}")
print(f"[genome_browser_hub] stranded={stranded}, tracks={[os.path.basename(b) for b in bw_files]}")

_log.close()
