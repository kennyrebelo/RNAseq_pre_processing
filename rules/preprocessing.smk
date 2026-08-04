# ------------------------------------------------------------------
# 1. SRA download (only invoked for samples whose read_pair1 is an
#    SRA accession; produces 1 or 2 fastq.gz files depending on
#    whether the run is single- or paired-end). Downloaded directly
#    into a subfolder of the sample's pre_processing_dir — no
#    raw_data_dir needed for SRA samples.
# ------------------------------------------------------------------
checkpoint download_sra:
    output:
        directory(os.path.join("{pre_processing_dir}", "sra_download_{sample}"))
    params:
        accession=lambda wc: samples_df.loc[wc.sample, "read_pair1"]
    log:
        os.path.join("{pre_processing_dir}", "logs", "download_sra", "{sample}.log")
    threads: 2
    conda:
        CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output}
        cd {output}
        fastq-dump --split-files --gzip {params.accession}
        """


# ------------------------------------------------------------------
# 2. Read trimming with Trim Galore (single-end and paired-end)
# ------------------------------------------------------------------
def trim_input(wc):
    return get_raw_fastqs(wc.sample)


def strip_fastq_ext(filepath):
    """
    trim_galore names its outputs from the INPUT file's basename, not
    from any sample name -- e.g. input 'foo_R1_001.fastq.gz' produces
    'foo_R1_001_trimmed.fq.gz', regardless of what Snakemake's {sample}
    wildcard is. This strips the fastq extension from an input path so
    we can predict trim_galore's actual output filename and rename it
    to the sample-based name our rule declares as output.
    """
    base = os.path.basename(filepath)
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return os.path.splitext(base)[0]


rule trim_galore_se:
    input:
        lambda wc: trim_input(wc) if not sample_is_paired(wc.sample) else []
    output:
        fq=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_trimmed.fq.gz"),
        html=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_trimmed_fastqc.html"),
        zip=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_trimmed_fastqc.zip"),
    params:
        outdir=lambda wc: sub_dir(wc.sample, "trimmed_reads"),
        in_base=lambda wc, input: strip_fastq_ext(input[0]),
    log:
        os.path.join("{pre_processing_dir}", "logs", "trim_galore_se", "{sample}.log")
    threads: config["resources"]["trimming"]["threads"]
    resources:
        mem_mb=config["resources"]["trimming"]["mem_mb"]
    conda:
        CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        trim_galore --gzip --fastqc --cores {threads} --output_dir {params.outdir} {input}
        mv {params.outdir}/{params.in_base}_trimmed.fq.gz {output.fq}
        mv {params.outdir}/{params.in_base}_trimmed_fastqc.html {output.html}
        mv {params.outdir}/{params.in_base}_trimmed_fastqc.zip {output.zip}
        """


rule trim_galore_pe:
    input:
        lambda wc: trim_input(wc) if sample_is_paired(wc.sample) else ["dummy_r1", "dummy_r2"]
    output:
        fq1=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_R1_val_1.fq.gz"),
        fq2=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_R2_val_2.fq.gz"),
        html1=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_R1_val_1_fastqc.html"),
        zip1=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_R1_val_1_fastqc.zip"),
        html2=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_R2_val_2_fastqc.html"),
        zip2=os.path.join("{pre_processing_dir}", "trimmed_reads", "{sample}_R2_val_2_fastqc.zip"),
    params:
        outdir=lambda wc: sub_dir(wc.sample, "trimmed_reads"),
        r1_base=lambda wc, input: strip_fastq_ext(input[0]),
        r2_base=lambda wc, input: strip_fastq_ext(input[1]),
    log:
        os.path.join("{pre_processing_dir}", "logs", "trim_galore_pe", "{sample}.log")
    threads: config["resources"]["trimming"]["threads"]
    resources:
        mem_mb=config["resources"]["trimming"]["mem_mb"]
    conda:
        CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        trim_galore --gzip --fastqc --cores {threads} --output_dir {params.outdir} --paired {input[0]} {input[1]}
        mv {params.outdir}/{params.r1_base}_val_1.fq.gz {output.fq1}
        mv {params.outdir}/{params.r2_base}_val_2.fq.gz {output.fq2}
        mv {params.outdir}/{params.r1_base}_val_1_fastqc.html {output.html1}
        mv {params.outdir}/{params.r1_base}_val_1_fastqc.zip {output.zip1}
        mv {params.outdir}/{params.r2_base}_val_2_fastqc.html {output.html2}
        mv {params.outdir}/{params.r2_base}_val_2_fastqc.zip {output.zip2}
        """


def trimmed_fastqs(sample):
    """Path(s) to trimmed fastq(s) for downstream (STAR / Salmon) rules."""
    trim_dir = sub_dir(sample, "trimmed_reads")
    if sample_is_paired(sample):
        return [
            os.path.join(trim_dir, f"{sample}_R1_val_1.fq.gz"),
            os.path.join(trim_dir, f"{sample}_R2_val_2.fq.gz"),
        ]
    return [os.path.join(trim_dir, f"{sample}_trimmed.fq.gz")]
