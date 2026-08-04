# ==================================================================
# RNA-seq pre-processing, bigwig, differential expression and
# splicing pipeline.
#
# This file is deliberately short: it loads the configuration and then
# includes one module per stage from rules/. Each module is plain
# Snakemake syntax -- `include:` is a textual include evaluated in
# order, so anything defined in an earlier module (helper functions,
# constants) is available to later ones. That ordering is why
# common.smk comes first.
#
#   rules/common.smk         config, sample sheet, design file, helpers
#   rules/preprocessing.smk  SRA download, Trim Galore
#   rules/alignment.smk      STAR, unique-read filtering, Salmon
#   rules/qc.smk             MultiQC
#   rules/bigwig.smk         strandedness, strand split, bigwigs, track hub
#   rules/dge.smk            DESeq2
#   rules/splicing.smk       StringTie + rMATS
#
# `rule all` sits here rather than in a module because Snakemake uses
# the FIRST rule it encounters as the default target -- putting it
# after the other includes would make download_sra the default.
# ==================================================================

configfile: "config_RNAseq_pipeline.yaml"


include: "rules/common.smk"


rule all:
    input:
        aggregate_final_outputs


include: "rules/preprocessing.smk"
include: "rules/alignment.smk"
include: "rules/qc.smk"
include: "rules/bigwig.smk"
include: "rules/dge.smk"
include: "rules/splicing.smk"
