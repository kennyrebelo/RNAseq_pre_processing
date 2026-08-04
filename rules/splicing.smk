# ==================================================================
# 9. Splicing quantification and novel-isoform annotation
#    (StringTie assembly + merge -> gene-id assignment -> rMATS
#     prep/post -> significance filter -> cryptic-exon annotation)
#
#    Skipped entirely when splicing.enabled is false. Samples taking
#    part are those whose group appears in a contrast pair.
# ==================================================================

SPLICING_DIR = SPLICING_OUTDIR
STRINGTIE_DIR = os.path.join(SPLICING_DIR, "stringtie")
RMATS_DIR = os.path.join(SPLICING_DIR, "rmats")
RMATS_FINAL_DIR = os.path.join(SPLICING_DIR, "rMATs_final_output")

MERGED_GTF = os.path.join(STRINGTIE_DIR, f"stringtie_merged_{SPLICING_TAG}.gtf")
ASSIGNED_GTF = os.path.join(STRINGTIE_DIR,
                            f"stringtie_merged_{SPLICING_TAG}_assigned.gtf")
RMATS_TMP = os.path.join(RMATS_DIR, f"tmp_merged_{SPLICING_TAG}")
LIBRARY_INFO = os.path.join(SPLICING_DIR, "library_info.txt")


def splicing_gtf():
    """Reference annotation for the (single) species in the splicing run."""
    return config["references"][SPLICING_SPECIES]["gtf"]


def sample_bam(sample):
    return unique_sorted_bam(sample)


def sample_fastqc_zip(sample):
    """
    Trim Galore's FastQC report for a sample, used to read the maximum read
    length. Resolves the SRA checkpoint for downloaded samples, so this may
    only be called from an input function.
    """
    trim_dir = sub_dir(sample, "trimmed_reads")
    if sample_is_paired(sample):
        return os.path.join(trim_dir, f"{sample}_R1_val_1_fastqc.zip")
    return os.path.join(trim_dir, f"{sample}_trimmed_fastqc.zip")


# ------------------------------------------------------------------
# 9a. Library type + read length for the whole splicing run
#
#     rMATS takes a single --libType, -t and --readLength for the entire
#     run, so these are resolved once across all splicing samples rather
#     than per sample. The values are written to a small text file that
#     downstream rules read at runtime -- this keeps the strandedness
#     result out of the DAG-building stage, where reading a checkpoint's
#     contents inside params would be fragile.
# ------------------------------------------------------------------
rule splicing_library_info:
    input:
        strand=lambda wc: [
            os.path.join(sub_dir(s, "bigwig"), f"{s}_infer_experiment.txt")
            for s in SPLICING_SAMPLES
        ],
        fastqc=lambda wc: [sample_fastqc_zip(s) for s in SPLICING_SAMPLES],
    output:
        info=LIBRARY_INFO,
    params:
        samples=SPLICING_SAMPLES,
        paired=lambda wc: [sample_is_paired(s) for s in SPLICING_SAMPLES],
        strand_threshold=SPL_CFG.get("strand_threshold",
                                     BW_CFG.get("strand_threshold", 0.75)),
    log:
        os.path.join(SPLICING_DIR, "logs", "splicing_library_info.log"),
    threads: 1
    resources:
        mem_mb=SPL_CFG.get("library_info", {}).get("mem_mb", 4000)
    conda:
        SPLICING_CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "splicing_library_info.py")


# ------------------------------------------------------------------
# 9b. Per-sample transcript assembly
# ------------------------------------------------------------------
rule stringtie_assemble:
    input:
        bam=lambda wc: sample_bam(wc.sample),
        bai=lambda wc: sample_bam(wc.sample) + ".bai",
        info=LIBRARY_INFO,
    output:
        gtf=temp(os.path.join(STRINGTIE_DIR, "ST_out_{sample}.gtf")),
    params:
        ref_gtf=lambda wc: splicing_gtf(),
        min_junction_coverage=SPL_CFG.get("stringtie_min_junction_coverage", 2),
        outdir=STRINGTIE_DIR,
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SPLICING_SAMPLES) or "[^/]+"
    log:
        os.path.join(SPLICING_DIR, "logs", "stringtie_assemble", "{sample}.log"),
    threads: SPL_CFG.get("stringtie", {}).get("threads", 4)
    resources:
        mem_mb=SPL_CFG.get("stringtie", {}).get("mem_mb", 5000)
    conda:
        SPLICING_CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}

        # StringTie's strand flag follows the inferred library type:
        #   fr-firststrand  -> --rf     fr-secondstrand -> --fr
        #   fr-unstranded   -> no flag
        STRAND_FLAG=$(grep '^stringtie_flag=' {input.info} | cut -d= -f2-)

        stringtie -o {output.gtf} $STRAND_FLAG \
            -G {params.ref_gtf} \
            -j {params.min_junction_coverage} \
            -p {threads} \
            {input.bam}
        """


# ------------------------------------------------------------------
# 9c. Merge the per-sample assemblies into one non-redundant set
# ------------------------------------------------------------------
rule stringtie_merge_list:
    input:
        gtfs=expand(os.path.join(STRINGTIE_DIR, "ST_out_{sample}.gtf"),
                    sample=SPLICING_SAMPLES),
    output:
        listfile=os.path.join(STRINGTIE_DIR,
                              f"get_list_to_merge_{SPLICING_TAG}.txt"),
    log:
        os.path.join(SPLICING_DIR, "logs", "stringtie_merge_list.log"),
    threads: 1
    resources:
        mem_mb=1000
    conda:
        SPLICING_CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        # One assembly per line, as stringtie --merge expects.
        printf '%s\\n' {input.gtfs} > {output.listfile}
        cat {output.listfile}
        """


rule stringtie_merge:
    input:
        listfile=os.path.join(STRINGTIE_DIR,
                              f"get_list_to_merge_{SPLICING_TAG}.txt"),
        gtfs=expand(os.path.join(STRINGTIE_DIR, "ST_out_{sample}.gtf"),
                    sample=SPLICING_SAMPLES),
    output:
        merged=temp(MERGED_GTF),
    params:
        ref_gtf=lambda wc: splicing_gtf(),
    log:
        os.path.join(SPLICING_DIR, "logs", "stringtie_merge.log"),
    threads: SPL_CFG.get("stringtie_merge", {}).get("threads", 24)
    resources:
        mem_mb=SPL_CFG.get("stringtie_merge", {}).get("mem_mb", 20000)
    conda:
        SPLICING_CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        stringtie --merge -G {params.ref_gtf} -p {threads} \
            -o {output.merged} {input.listfile}
        """


# ------------------------------------------------------------------
# 9d. Give the novel isoforms a gene identity
#
#     stringtie --merge emits novel transcripts with MSTRG gene ids and no
#     gene_name. This assigns them to an annotated gene when they share an
#     exon with one. Only the assigned GTF is kept -- the per-sample
#     assemblies and the raw merge are marked temp() above, so Snakemake
#     removes them once this rule succeeds.
# ------------------------------------------------------------------
rule assign_gene_ids:
    input:
        merged=MERGED_GTF,
    output:
        assigned=ASSIGNED_GTF,
    log:
        os.path.join(SPLICING_DIR, "logs", "assign_gene_ids.log"),
    threads: 1
    resources:
        mem_mb=SPL_CFG.get("assign_gene_ids", {}).get("mem_mb", 16000)
    conda:
        SPLICING_CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "assign_gene_ids.py")


# ------------------------------------------------------------------
# 9e. rMATS prep -- run once over every splicing sample
# ------------------------------------------------------------------
rule rmats_prep_bam_list:
    input:
        bams=lambda wc: [sample_bam(s) for s in SPLICING_SAMPLES],
        bais=lambda wc: [sample_bam(s) + ".bai" for s in SPLICING_SAMPLES],
    output:
        listfile=os.path.join(RMATS_DIR, f"prep_samples_{SPLICING_TAG}.txt"),
    log:
        os.path.join(SPLICING_DIR, "logs", "rmats_prep_bam_list.log"),
    threads: 1
    resources:
        mem_mb=1000
    conda:
        SPLICING_CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        # rMATS expects one line of comma-separated BAM paths.
        printf '%s' "$(echo {input.bams} | tr ' ' ',')" > {output.listfile}
        echo
        cat {output.listfile}
        """


rule rmats_prep:
    input:
        listfile=os.path.join(RMATS_DIR, f"prep_samples_{SPLICING_TAG}.txt"),
        gtf=ASSIGNED_GTF,
        info=LIBRARY_INFO,
    output:
        tmpdir=directory(RMATS_TMP),
        done=touch(os.path.join(RMATS_DIR, f"prep_{SPLICING_TAG}.done")),
    params:
        anchor_length=SPL_CFG.get("anchor_length", 4),
        outdir=os.path.join(RMATS_DIR, f"prep_output_{SPLICING_TAG}"),
    log:
        os.path.join(SPLICING_DIR, "logs", "rmats_prep.log"),
    threads: SPL_CFG.get("rmats_prep", {}).get("threads", 20)
    resources:
        mem_mb=SPL_CFG.get("rmats_prep", {}).get("mem_mb", 10000)
    conda:
        SPLICING_CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output.tmpdir} {params.outdir}

        LIBTYPE=$(grep '^libtype=' {input.info} | cut -d= -f2-)
        READTYPE=$(grep '^read_type=' {input.info} | cut -d= -f2-)
        READLENGTH=$(grep '^read_length=' {input.info} | cut -d= -f2-)
        echo "libType=$LIBTYPE  -t=$READTYPE  readLength=$READLENGTH"

        rmats.py --b1 {input.listfile} \
            --gtf {input.gtf} \
            -t $READTYPE \
            --readLength $READLENGTH \
            --variable-read-length \
            --anchorLength {params.anchor_length} \
            --nthread {threads} \
            --od {params.outdir} \
            --tmp {output.tmpdir} \
            --libType $LIBTYPE \
            --task prep
        """


# ------------------------------------------------------------------
# 9f. rMATS post -- once per contrast pair
# ------------------------------------------------------------------
rule rmats_contrast_bam_lists:
    input:
        bams=lambda wc: [
            sample_bam(s)
            for s in group_samples(wc.contrast_a) + group_samples(wc.contrast_b)
        ],
    output:
        b1=os.path.join(RMATS_DIR, "bam_lists", "{contrast_a}_vs_{contrast_b}_b1.txt"),
        b2=os.path.join(RMATS_DIR, "bam_lists", "{contrast_a}_vs_{contrast_b}_b2.txt"),
    params:
        b1_bams=lambda wc: ",".join(sample_bam(s) for s in group_samples(wc.contrast_a)),
        b2_bams=lambda wc: ",".join(sample_bam(s) for s in group_samples(wc.contrast_b)),
    log:
        os.path.join(SPLICING_DIR, "logs", "rmats_bam_lists",
                     "{contrast_a}_vs_{contrast_b}.log"),
    threads: 1
    resources:
        mem_mb=1000
    conda:
        SPLICING_CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        printf '%s' '{params.b1_bams}' > {output.b1}
        printf '%s' '{params.b2_bams}' > {output.b2}
        echo "b1:"; cat {output.b1}; echo
        echo "b2:"; cat {output.b2}; echo
        """


rule rmats_post:
    input:
        b1=os.path.join(RMATS_DIR, "bam_lists", "{contrast_a}_vs_{contrast_b}_b1.txt"),
        b2=os.path.join(RMATS_DIR, "bam_lists", "{contrast_a}_vs_{contrast_b}_b2.txt"),
        gtf=ASSIGNED_GTF,
        info=LIBRARY_INFO,
        tmpdir=RMATS_TMP,
        prep_done=os.path.join(RMATS_DIR, f"prep_{SPLICING_TAG}.done"),
    output:
        se=os.path.join(RMATS_DIR, "rmats_output_{contrast_a}_vs_{contrast_b}",
                        "SE.MATS.JC.txt"),
    params:
        outdir=lambda wc: os.path.join(
            RMATS_DIR, f"rmats_output_{wc.contrast_a}_vs_{wc.contrast_b}"
        ),
        anchor_length=SPL_CFG.get("anchor_length", 4),
        cstat=SPL_CFG.get("cstat", 0.01),
    log:
        os.path.join(SPLICING_DIR, "logs", "rmats_post",
                     "{contrast_a}_vs_{contrast_b}.log"),
    threads: SPL_CFG.get("rmats_post", {}).get("threads", 4)
    resources:
        mem_mb=SPL_CFG.get("rmats_post", {}).get("mem_mb", 10000)
    conda:
        SPLICING_CONDA_ENV
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}

        LIBTYPE=$(grep '^libtype=' {input.info} | cut -d= -f2-)
        READTYPE=$(grep '^read_type=' {input.info} | cut -d= -f2-)
        READLENGTH=$(grep '^read_length=' {input.info} | cut -d= -f2-)

        rmats.py --b1 {input.b1} --b2 {input.b2} \
            --gtf {input.gtf} \
            -t $READTYPE \
            --readLength $READLENGTH \
            --variable-read-length \
            --anchorLength {params.anchor_length} \
            --nthread {threads} \
            --od {params.outdir} \
            --tmp {input.tmpdir} \
            --libType $LIBTYPE \
            --cstat {params.cstat} \
            --task post
        """


# ------------------------------------------------------------------
# 9g. Keep only significant, well-supported skipped-exon events
# ------------------------------------------------------------------
rule filter_se_events:
    input:
        se=os.path.join(RMATS_DIR, "rmats_output_{contrast_a}_vs_{contrast_b}",
                        "SE.MATS.JC.txt"),
    output:
        sig=os.path.join(RMATS_DIR, "rmats_output_{contrast_a}_vs_{contrast_b}",
                         "SE.MATS.JC_sig.txt"),
    params:
        min_reads=SPL_CFG.get("min_reads", 2),
        fdr=SPL_CFG.get("fdr", 0.05),
        min_delta_psi=SPL_CFG.get("min_delta_psi", 0.1),
    log:
        os.path.join(SPLICING_DIR, "logs", "filter_se_events",
                     "{contrast_a}_vs_{contrast_b}.log"),
    threads: 1
    resources:
        mem_mb=4000
    conda:
        SPLICING_CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "filter_se_events.py")


# ------------------------------------------------------------------
# 9h. Flag each surviving exon as annotated or cryptic
# ------------------------------------------------------------------
rule identify_cryptic_exons:
    input:
        sig=os.path.join(RMATS_DIR, "rmats_output_{contrast_a}_vs_{contrast_b}",
                         "SE.MATS.JC_sig.txt"),
    output:
        ided=os.path.join(RMATS_FINAL_DIR,
                          "SE.MATS.JC_sig_IDed_{contrast_a}_vs_{contrast_b}.txt"),
    params:
        ref_gtf=lambda wc: splicing_gtf(),
    log:
        os.path.join(SPLICING_DIR, "logs", "identify_cryptic_exons",
                     "{contrast_a}_vs_{contrast_b}.log"),
    threads: 1
    resources:
        mem_mb=SPL_CFG.get("identify_cryptic", {}).get("mem_mb", 16000)
    conda:
        SPLICING_CONDA_ENV
    script:
        os.path.join(SCRIPTS_DIR, "identify_cryptic_exons.py")
