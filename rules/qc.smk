# ------------------------------------------------------------------
# 6. MultiQC — collects STAR logs + FastQC outputs into QC/ and runs
#    multiqc there
# ------------------------------------------------------------------
def all_qc_inputs(wildcards):
    files = []
    for sample in SAMPLES:
        files.append(
            os.path.join(sub_dir(sample, "star_alignment"), f"{sample}_Log.final.out")
        )
        if sample_is_paired(sample):
            files.append(
                os.path.join(sub_dir(sample, "trimmed_reads"), f"{sample}_R1_val_1_fastqc.zip")
            )
            files.append(
                os.path.join(sub_dir(sample, "trimmed_reads"), f"{sample}_R2_val_2_fastqc.zip")
            )
        else:
            files.append(
                os.path.join(sub_dir(sample, "trimmed_reads"), f"{sample}_trimmed_fastqc.zip")
            )
    return files


rule collect_and_run_multiqc:
    input:
        all_qc_inputs
    output:
        "QC/multiqc_report.html"
    log:
        # This rule aggregates across all samples, so it has no {sample}
        # wildcard to nest under -- its log lives beside the report instead.
        "QC/logs/collect_and_run_multiqc.log"
    conda:
        CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p QC
        cp {input} QC/
        cd QC && multiqc --force ./
        """
