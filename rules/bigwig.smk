# ==================================================================
# 7. BigWig / genome-browser branch
#    Only runs for samples with generate_bw = yes in samples.tsv.
#
#    infer_strandedness (checkpoint)
#         |
#         +-- stranded -----> split_strands -> strand_sam_to_sorted_bam
#         |                              -> make_bigwig_stranded (F + R)
#         |
#         +-- unstranded ---> make_bigwig_unstranded (single track)
#                                          |
#                                    genome_browser_hub
# ==================================================================

# ------------------------------------------------------------------
# 7a. Strandedness inference with RSeQC
#
# infer_experiment.py reports the fraction of reads whose mapping strand is
# consistent with each of two mutually exclusive configurations. The notation
# reads "read number, read strand, transcript strand", so 1++ means read 1 on
# the + strand of a + strand transcript.
#
#   group 1  paired: 1++,1--,2+-,2-+     single: ++,--
#            -> read 1 is SENSE to the transcript
#   group 2  paired: 1+-,1-+,2++,2--     single: +-,-+
#            -> read 1 is ANTISENSE to the transcript
#
# Mapping onto the library types rMATS (--libType) and StringTie expect:
#
#   group 1 dominant   -> fr-secondstrand   (StringTie --fr)
#                         read 1 sense; ligation / directional mRNA kits
#   group 2 dominant   -> fr-firststrand    (StringTie --rf)
#                         read 1 antisense; dUTP protocols such as TruSeq
#                         Stranded -- the most common current chemistry
#   both around 0.5    -> fr-unstranded     (StringTie: no strand flag)
#                         non-directional library
#
# Equivalents elsewhere: Salmon ISR / ISF / IU (SR / SF / U for single-end),
# HISAT2 --rna-strandness RF / FR, htseq-count reverse / yes / no.
#
# Two caveats. The "fr-" prefix also asserts inward-facing mates, which
# infer_experiment.py does not test -- it measures strandedness only, and
# standard Illumina paired-end libraries satisfy the orientation half. And
# "Fraction of reads failed to determine" is the confidence check: a few
# percent is normal (regions where both strands are transcribed), but a large
# fraction means the transcript BED does not match the genome the BAM was
# aligned to, which invalidates the call rather than softening it.
#
# An intermediate result such as 0.65 / 0.33 is not a valid middle ground; it
# usually signals degraded input, an annotation mismatch, or a mixed library.
# splicing.strand_threshold governs how confident the call has to be.
#
# This checkpoint runs for any sample with generate_bw = yes, and for every
# sample taking part in the splicing analysis -- the latter regardless of
# generate_bw, since rMATS and StringTie cannot be run without knowing the
# library type.
# ------------------------------------------------------------------
checkpoint infer_strandedness:
    input:
        bam=os.path.join("{pre_processing_dir}", "star_alignment",
                         "{sample}_Aligned.out.unique.sorted.bam"),
        bai=os.path.join("{pre_processing_dir}", "star_alignment",
                         "{sample}_Aligned.out.unique.sorted.bam.bai"),
    output:
        txt=os.path.join("{pre_processing_dir}", "bigwig",
                         "{sample}_infer_experiment.txt"),
    params:
        bed=lambda wc: transcripts_bed(wc.sample),
        outdir=lambda wc: sub_dir(wc.sample, "bigwig"),
    log:
        os.path.join("{pre_processing_dir}", "logs", "infer_strandedness", "{sample}.log")
    threads: 1
    resources:
        mem_mb=config["resources"]["sorting"]["mem_mb"]
    conda:
        CONDA_ENV
    shell:
        """
        mkdir -p {params.outdir}
        infer_experiment.py -r {params.bed} -i {input.bam} > {output.txt} 2> {log}
        """


def parse_infer_experiment(path):
    """
    Reads RSeQC infer_experiment.py output and returns True when the library
    looks strand-specific.

    RSeQC reports two 'Fraction of reads explained by ...' values. For an
    unstranded library they are ~0.5 / 0.5; for a stranded library one of
    them dominates. We call the library stranded when the larger fraction
    reaches STRAND_THRESHOLD.
    """
    fractions = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('Fraction of reads explained by "'):
                try:
                    fractions.append(float(line.rsplit(":", 1)[1]))
                except ValueError:
                    continue
    if not fractions:
        raise ValueError(
            f"Could not parse any strandedness fractions from {path}. "
            "Check that infer_experiment.py ran successfully and that the "
            "transcript BED matches the genome the BAM was aligned to."
        )
    return max(fractions) >= STRAND_THRESHOLD


def is_stranded(sample):
    """
    Resolves the infer_strandedness checkpoint for a sample. Only call this
    from inside an input function -- never at parse time.
    """
    ckpt = checkpoints.infer_strandedness.get(
        sample=sample, pre_processing_dir=preproc_dir(sample)
    )
    return parse_infer_experiment(ckpt.output.txt)


def unique_sorted_bam(sample):
    return os.path.join(
        sub_dir(sample, "star_alignment"),
        f"{sample}_Aligned.out.unique.sorted.bam",
    )


# ------------------------------------------------------------------
# 7b. Split alignments into forward / reverse strand SAM files
#     (stranded libraries only). The SAMs are marked temp() so
#     Snakemake removes them automatically once the sorted BAMs exist.
# ------------------------------------------------------------------
rule split_strands:
    input:
        bam=lambda wc: unique_sorted_bam(wc.sample),
        bai=lambda wc: unique_sorted_bam(wc.sample) + ".bai",
        strand_txt=os.path.join("{pre_processing_dir}", "bigwig",
                                "{sample}_infer_experiment.txt"),
    output:
        fwd=temp(os.path.join("{pre_processing_dir}", "bigwig", "F_{sample}.sam")),
        rev=temp(os.path.join("{pre_processing_dir}", "bigwig", "R_{sample}.sam")),
    log:
        os.path.join("{pre_processing_dir}", "logs", "split_strands", "{sample}.log")
    threads: 1
    resources:
        mem_mb=config["resources"]["sorting"]["mem_mb"]
    conda:
        CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "split_strands.py")


# ------------------------------------------------------------------
# 7c. SAM -> BAM -> sorted BAM, using the same thread count as the
#     post-alignment sorting step. Intermediates are removed.
# ------------------------------------------------------------------
rule strand_sam_to_sorted_bam:
    input:
        sam=os.path.join("{pre_processing_dir}", "bigwig", "{strand}_{sample}.sam"),
    output:
        sorted_bam=temp(
            os.path.join("{pre_processing_dir}", "bigwig",
                         "{strand}_{sample}_sorted.bam")
        ),
    params:
        unsorted=lambda wc: os.path.join(
            sub_dir(wc.sample, "bigwig"), f"{wc.strand}_{wc.sample}.bam"
        ),
    wildcard_constraints:
        strand="[FR]"
    log:
        os.path.join("{pre_processing_dir}", "logs", "strand_sam_to_sorted_bam",
                     "{strand}_{sample}.log")
    threads: config["resources"]["sorting"]["threads"]
    resources:
        mem_mb=config["resources"]["sorting"]["mem_mb"]
    conda:
        CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        samtools view -bS {input.sam} > {params.unsorted}
        samtools sort -@ {threads} -o {output.sorted_bam} {params.unsorted}
        rm -f {params.unsorted}
        """


# ------------------------------------------------------------------
# 7d. Mapped-read count, used to normalise coverage
# ------------------------------------------------------------------
rule count_mapped_reads:
    input:
        bam=lambda wc: unique_sorted_bam(wc.sample),
    output:
        counts=os.path.join("{pre_processing_dir}", "bigwig",
                            "{sample}_read_count.txt"),
    log:
        os.path.join("{pre_processing_dir}", "logs", "count_mapped_reads", "{sample}.log")
    threads: 1
    resources:
        mem_mb=config["resources"]["sorting"]["mem_mb"]
    conda:
        CONDA_ENV
    shell:
        "samtools view -c {input.bam} > {output.counts} 2> {log}"


# ------------------------------------------------------------------
# 7e. BigWig generation
# ------------------------------------------------------------------
rule make_bigwig_stranded:
    input:
        fwd_bam=os.path.join("{pre_processing_dir}", "bigwig", "F_{sample}_sorted.bam"),
        rev_bam=os.path.join("{pre_processing_dir}", "bigwig", "R_{sample}_sorted.bam"),
        counts=os.path.join("{pre_processing_dir}", "bigwig", "{sample}_read_count.txt"),
    output:
        fwd_bw=os.path.join("{pre_processing_dir}", "bigwig",
                            "{sample}_F_" + SCALE_SUFFIX + ".bw"),
        rev_bw=os.path.join("{pre_processing_dir}", "bigwig",
                            "{sample}_R_" + SCALE_SUFFIX + ".bw"),
    params:
        mode="stranded",
        chrom_sizes=lambda wc: chrom_sizes(wc.sample),
        bedgraph_to_bigwig=BEDGRAPH_TO_BIGWIG,
        scale_factor=SCALE_FACTOR,
        workdir=lambda wc: sub_dir(wc.sample, "bigwig"),
    log:
        os.path.join("{pre_processing_dir}", "logs", "make_bigwig", "{sample}_stranded.log")
    threads: config["resources"]["bigwig"]["threads"]
    resources:
        mem_mb=config["resources"]["bigwig"]["mem_mb"]
    conda:
        CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "make_bigwig.py")


rule make_bigwig_unstranded:
    input:
        bam=lambda wc: unique_sorted_bam(wc.sample),
        counts=os.path.join("{pre_processing_dir}", "bigwig", "{sample}_read_count.txt"),
        strand_txt=os.path.join("{pre_processing_dir}", "bigwig",
                                "{sample}_infer_experiment.txt"),
    output:
        bw=os.path.join("{pre_processing_dir}", "bigwig",
                        "{sample}_unstranded_" + SCALE_SUFFIX + ".bw"),
    params:
        mode="unstranded",
        chrom_sizes=lambda wc: chrom_sizes(wc.sample),
        bedgraph_to_bigwig=BEDGRAPH_TO_BIGWIG,
        scale_factor=SCALE_FACTOR,
        workdir=lambda wc: sub_dir(wc.sample, "bigwig"),
    log:
        os.path.join("{pre_processing_dir}", "logs", "make_bigwig", "{sample}_unstranded.log")
    threads: config["resources"]["bigwig"]["threads"]
    resources:
        mem_mb=config["resources"]["bigwig"]["mem_mb"]
    conda:
        CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "make_bigwig.py")


def bigwig_files(sample):
    """
    Which bigwig(s) this sample produces -- two (F/R) for a strand-specific
    library, one for an unstranded one. Resolves the strandedness checkpoint,
    so only call from inside an input function.
    """
    bw_dir = sub_dir(sample, "bigwig")
    if is_stranded(sample):
        return [
            os.path.join(bw_dir, f"{sample}_F_{SCALE_SUFFIX}.bw"),
            os.path.join(bw_dir, f"{sample}_R_{SCALE_SUFFIX}.bw"),
        ]
    return [os.path.join(bw_dir, f"{sample}_unstranded_{SCALE_SUFFIX}.bw")]


# ------------------------------------------------------------------
# 7f. UCSC track hub folder
# ------------------------------------------------------------------
rule genome_browser_hub:
    input:
        bw=lambda wc: bigwig_files(wc.sample),
    output:
        hub_dir=directory(os.path.join("{pre_processing_dir}", "genome_browser", "{sample}")),
    params:
        sample=lambda wc: wc.sample,
        genome=lambda wc: genome_build(wc.sample),
        email=HUB_EMAIL,
    log:
        os.path.join("{pre_processing_dir}", "logs", "genome_browser_hub", "{sample}.log")
    threads: 1
    resources:
        mem_mb=config["resources"]["sorting"]["mem_mb"]
    conda:
        CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "make_genome_browser_hub.py")
