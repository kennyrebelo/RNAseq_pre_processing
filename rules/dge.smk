# ==================================================================
# 8. Differential gene expression with DESeq2
#    One table + one PCA plot per contrast pair defined in
#    experimental_design.tsv. Skipped entirely when dge.enabled is false.
# ==================================================================

# Contrast group names routinely contain underscores (e.g. "kd_c"), so an
# unconstrained wildcard would split "nt_vs_kd_c" ambiguously. Restrict both
# wildcards to the exact group names parsed from the design file.
_CONTRAST_NAMES = sorted({n for pair in CONTRASTS for n in pair}, key=len, reverse=True)
CONTRAST_REGEX = (
    "|".join(re.escape(n) for n in _CONTRAST_NAMES) if _CONTRAST_NAMES else "[^/]+"
)

wildcard_constraints:
    contrast_a=CONTRAST_REGEX,
    contrast_b=CONTRAST_REGEX,


def contrast_quant_files(contrast_a, contrast_b):
    return [
        os.path.join(sub_dir(s, "salmon_quant"), s, "quant.sf")
        for s in contrast_samples(contrast_a, contrast_b)
    ]


rule deseq2_contrast:
    input:
        quant=lambda wc: contrast_quant_files(wc.contrast_a, wc.contrast_b),
    output:
        table=os.path.join(DGE_OUTDIR,
                           "DESeq2_salmon_diff_{contrast_a}_vs_{contrast_b}.tab"),
        pca=os.path.join(DGE_OUTDIR,
                         "DESeq2_salmon_PCA_{contrast_a}_vs_{contrast_b}.pdf"),
        pca_labels=os.path.join(DGE_OUTDIR,
                                "DESeq2_salmon_PCA_{contrast_a}_vs_{contrast_b}_labels.pdf"),
    params:
        # Aligned 1:1 with input.quant, so the R script never has to infer
        # sample identity or condition from file paths.
        sample_names=lambda wc: contrast_samples(wc.contrast_a, wc.contrast_b),
        sample_groups=lambda wc: contrast_groups(wc.contrast_a, wc.contrast_b),
        reference=lambda wc: wc.contrast_a,
        treatment=lambda wc: wc.contrast_b,
        counts_from_abundance=DGE_COUNTS_FROM_ABUNDANCE,
    log:
        os.path.join(DGE_OUTDIR, "logs",
                     "deseq2_{contrast_a}_vs_{contrast_b}.log"),
    threads: DGE_CFG.get("threads", 4)
    resources:
        mem_mb=DGE_CFG.get("mem_mb", 16000)
    conda:
        DGE_CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "deseq2_contrast.R")
