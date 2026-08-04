# ------------------------------------------------------------------
# 3. Alignment with STAR
# ------------------------------------------------------------------
rule star_align:
    input:
        fq=lambda wc: trimmed_fastqs(wc.sample)
    output:
        bam=os.path.join("{pre_processing_dir}", "star_alignment", "{sample}_Aligned.out.bam"),
        log=os.path.join("{pre_processing_dir}", "star_alignment", "{sample}_Log.final.out"),
    params:
        outprefix=lambda wc: os.path.join(sub_dir(wc.sample, "star_alignment"), f"{wc.sample}_"),
        outdir=lambda wc: sub_dir(wc.sample, "star_alignment"),
        genome_dir=lambda wc: star_genome_dir(wc.sample),
    log:
        os.path.join("{pre_processing_dir}", "logs", "star_align", "{sample}.log")
    threads: config["resources"]["alignment"]["threads"]
    resources:
        mem_mb=config["resources"]["alignment"]["mem_mb"]
    conda:
        CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        STAR --runThreadN {threads} --runMode alignReads \
             --genomeDir {params.genome_dir} \
             --readFilesIn {input.fq} \
             --chimSegmentMin 20 \
             --outSAMtype BAM Unsorted \
             --readFilesCommand gunzip -c \
             --outFileNamePrefix {params.outprefix}
        """


# ------------------------------------------------------------------
# 4. Filter uniquely mapped reads, sort, index
# ------------------------------------------------------------------
rule filter_sort_index:
    input:
        bam=os.path.join("{pre_processing_dir}", "star_alignment", "{sample}_Aligned.out.bam"),
    output:
        sorted_bam=protected(
            os.path.join("{pre_processing_dir}", "star_alignment",
                         "{sample}_Aligned.out.unique.sorted.bam")
        ),
        bai=protected(
            os.path.join("{pre_processing_dir}", "star_alignment",
                         "{sample}_Aligned.out.unique.sorted.bam.bai")
        ),
    params:
        outdir=lambda wc: sub_dir(wc.sample, "star_alignment"),
        unique_bam=lambda wc: os.path.join(
            sub_dir(wc.sample, "star_alignment"), f"{wc.sample}_Aligned.out.unique.bam"
        ),
    log:
        os.path.join("{pre_processing_dir}", "logs", "filter_sort_index", "{sample}.log")
    threads: config["resources"]["sorting"]["threads"]
    resources:
        mem_mb=config["resources"]["sorting"]["mem_mb"]
    conda:
        CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        samtools view -q 255 -b {input.bam} > {params.unique_bam}
        samtools sort -@ {threads} -o {output.sorted_bam} {params.unique_bam}
        samtools index {output.sorted_bam}
        rm -f {params.unique_bam} {input.bam}
        """


# ------------------------------------------------------------------
# 5. Quantification with Salmon
# ------------------------------------------------------------------
def salmon_reads_arg(sample):
    """
    Build the salmon read-input flags as a plain string, resolved in
    Python rather than inside the shell template. This avoids referencing
    an input.fq[1] index that doesn't exist for single-end samples --
    Snakemake formats the *entire* shell string (including untaken
    bash branches) when logging a job, so any literal index reference
    must always be valid regardless of pairing.
    """
    files = trimmed_fastqs(sample)
    if sample_is_paired(sample):
        return f"-1 {files[0]} -2 {files[1]}"
    return f"-r {files[0]}"


rule salmon_quant:
    input:
        fq=lambda wc: trimmed_fastqs(wc.sample)
    output:
        quant=os.path.join("{pre_processing_dir}", "salmon_quant", "{sample}", "quant.sf")
    params:
        outdir=lambda wc: os.path.join(sub_dir(wc.sample, "salmon_quant"), wc.sample),
        index=lambda wc: salmon_index_dir(wc.sample),
        reads_arg=lambda wc: salmon_reads_arg(wc.sample),
    log:
        os.path.join("{pre_processing_dir}", "logs", "salmon_quant", "{sample}.log")
    threads: config["resources"]["alignment"]["threads"]
    resources:
        mem_mb=config["resources"]["alignment"]["mem_mb"]
    conda:
        CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        salmon quant -i {params.index} -l A {params.reads_arg} \
            -p {threads} --validateMappings -o {params.outdir}
        """
