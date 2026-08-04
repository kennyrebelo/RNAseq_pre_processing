# Configuration

Three files describe a run: `samples.tsv` (what data), `experimental_design.tsv`
(how to compare it), and `config_RNAseq_pipeline.yaml` (everything else). The
Snakefile and the modules under `rules/` should not need editing.

---

## 1. `samples.tsv`

Tab-separated, one row per sample.

| Column | Required | Description |
|---|---|---|
| `sample_name` | yes | Name used throughout the pipeline; replaces the raw filename. Must be unique and must not contain `/` |
| `species` | yes | Must match a key under `references:` in the config (`human` or `mouse` as shipped) |
| `read_pair1` | yes | Read 1 filename, **or** an SRA accession (`SRR…`, `ERR…`, `DRR…`) |
| `read_pair2` | no | Read 2 filename. Blank means single-end |
| `raw_data_dir` | for local samples | Directory holding the FASTQs. Leave blank for SRA samples |
| `pre_processing_dir` | yes | Directory under which this sample's outputs are created |
| `generate_bw` | no | `yes` or `no`. Defaults to `no` if the column is absent entirely |

```tsv
sample_name	species	read_pair1	read_pair2	raw_data_dir	pre_processing_dir	generate_bw
total_nt_rep1	human	NT_rep1_R1.fastq.gz	NT_rep1_R2.fastq.gz	/data/raw	/data/processed	yes
sample_se	mouse	SE_sample.fastq.gz		/data/raw	/data/processed	no
sample_sra	human	SRR398210			/data/processed	no
```

**Inferred, not declared:**

- **Single vs paired-end** — from whether `read_pair2` is filled. For SRA
  samples it comes from how many FASTQs `fastq-dump --split-files` returns.
- **Strandedness** — from RSeQC, per sample.
- **Read length** — from the FastQC report, as the maximum observed length.

Samples may mix species and mix single/paired-end freely, *except* within the
splicing branch (see below).

---

## 2. `experimental_design.tsv`

Required only when `dge.enabled` or `splicing.enabled` is true.

| Column | Description |
|---|---|
| `sample_name` | Must match a `sample_name` in `samples.tsv` |
| `group` | Condition this replicate belongs to |
| `contrast_A` | Reference / baseline group |
| `contrast_B` | Treatment group |

Fill `contrast_A`/`contrast_B` on as many rows as you have comparisons; leave
them blank on the rest. Row position is irrelevant.

```tsv
sample_name	group	contrast_A	contrast_B
total_nt_rep1	nt	nt	kd_c
total_nt_rep2	nt
total_nt_rep3	nt
total_kd_c_rep1	kd_c
total_kd_c_rep2	kd_c
total_kd_c_rep3	kd_c
```

**Direction:** `contrast_A` is the reference. `nt` / `kd_c` reports
log2FoldChange of **kd_c relative to nt** — positive means higher in `kd_c`.

**Which samples participate:** every sample whose `group` appears in some
contrast. A sample listed in the design but not in any contrast is excluded
from both DGE and splicing.

**Validated at DAG-build time**, before any job is submitted:

- every design sample exists in `samples.tsv`
- every contrast group exists in the `group` column
- no group is contrasted against itself
- (DGE only) each group has ≥ 2 replicates — DESeq2 cannot estimate dispersion
  otherwise
- (splicing only) all participating samples share one species

---

## 3. `config_RNAseq_pipeline.yaml`

### References

Per species. Add further species by adding another key.

```yaml
references:
  human:
    star_genome_dir:  "/path/to/STAR_index"
    salmon_index_dir: "/path/to/salmon_index"
    transcripts_bed:  "/path/to/transcripts.bed"      # bigwig + splicing
    chrom_sizes:      "/path/to/hg38.chrom.sizes"     # bigwig
    genome_build:     "hg38"                          # track hub label
    gtf:              "/path/to/gencode.v41.gtf"      # splicing
```

`transcripts_bed`, `chrom_sizes`, `genome_build` and `gtf` are only required
for species used by the branches that need them; missing keys are caught at
parse time with a message naming the sample and the key.

### Resources

```yaml
resources:
  trimming:  {threads: 8,  mem_mb: 10000}
  alignment: {threads: 20, mem_mb: 40000}   # STAR and Salmon
  sorting:   {threads: 8,  mem_mb: 10000}   # also strand-BAM sorting
  bigwig:    {threads: 4,  mem_mb: 16000}
```

`disk_mb` is not set here — it belongs in the cluster profile's
`set-resources`, since it is a scheduler concern rather than a pipeline one.

### BigWig branch

```yaml
bigwig:
  bedgraph_to_bigwig: "/path/to/bedGraphToBigWig"
  hub_email: "you@institution.edu"
  scale_factor: 100000000      # coverage × factor / mapped reads
  scale_suffix: "RP10_8"       # appears in output filenames
  strand_threshold: 0.75
```

### Differential expression

```yaml
dge:
  enabled: true
  experimental_design: "experimental_design.tsv"
  output_dir: "DGE"
  conda_env: "environment_DESeq2.yaml"
  counts_from_abundance: "lengthScaledTPM"    # passed to tximport
  threads: 4
  mem_mb: 16000
```

### Splicing

```yaml
splicing:
  enabled: true
  experimental_design: "experimental_design.tsv"
  output_dir: "splicing"
  conda_env: "environment_rMATS.yaml"

  stringtie_min_junction_coverage: 2    # StringTie -j
  stringtie:       {threads: 4,  mem_mb: 5000}
  stringtie_merge: {threads: 24, mem_mb: 20000}

  anchor_length: 4     # rMATS --anchorLength: minimum nucleotides that must
                       # map to each end of a given splice junction
  cstat: 0.01          # rMATS --cstat: splicing-difference cutoff
  rmats_prep: {threads: 20, mem_mb: 10000}
  rmats_post: {threads: 4,  mem_mb: 10000}

  strand_threshold: 0.75

  min_reads: 2         # per-replicate support, mean of inclusion and skipping
  fdr: 0.05
  min_delta_psi: 0.1
```

### The two master switches

| `dge.enabled` | `splicing.enabled` | Effect |
|---|---|---|
| false | false | Pre-processing + QC only. `experimental_design.tsv` not required or read |
| true | false | Adds DESeq2 |
| false | true | Adds StringTie + rMATS, and forces strandedness inference for participating samples regardless of `generate_bw` |
| true | true | Both |

---

## What is **not** configurable

Deliberately fixed, to keep the pipeline's behaviour predictable. Changing any
of these means editing the relevant `rules/*.smk` file.

### Tool parameters

| Tool | Fixed setting | Rationale |
|---|---|---|
| Trim Galore | `--gzip --fastqc`, automatic adapter detection | Adapter type is detected per sample |
| STAR | `--chimSegmentMin 20`, `--outSAMtype BAM Unsorted`, `--readFilesCommand gunzip -c` | Chimeric detection on; sorting is a separate step |
| samtools | `view -q 255` for unique reads | STAR's MAPQ 255 convention for unique alignments |
| Salmon | `-l A`, `--validateMappings` | Library type auto-detected by Salmon itself |
| RSeQC | default sampling depth | |
| StringTie | `-G <gtf>`, strand flag from inference | |
| rMATS | `--variable-read-length` | Trimming produces variable lengths |
| DESeq2 | design `~group`, Wald test, `vst(blind=FALSE)` for PCA | Single-factor comparison |

### Behaviour

- **Library type, read length and paired/single are inferred**, never set by
  hand. To override, adjust `strand_threshold`.
- **rMATS runs once per species, per run.** It takes one `--libType`, one `-t`
  and one `--readLength`, so all participating samples must agree. Disagreement
  is an error naming the offending samples, not something averaged away.
- **Only SE (skipped exon) events are carried downstream.** rMATS produces
  A3SS, A5SS, MXE and RI too, and they are written to the rMATS output
  directory, but the significance filter and cryptic-exon annotation only
  process `SE.MATS.JC.txt`.
- **JC, not JCEC.** Junction counts only; reads on the exon body are not used.
- **Strand assignment assumes a reverse-stranded (dUTP) library.**
  `scripts/split_strands.py` writes reverse-strand read 1 to the forward
  bucket. For forward-stranded libraries, flip the four branches in that script
  — and only that script, since the historical inversion in the bedgraph step
  has been removed.
- **Gene identity for DESeq2 is parsed from GENCODE FASTA headers** in
  `quant.sf`, not from a separate tx2gene table.
- **Novel isoforms are assigned by exact shared exon.** A StringTie transcript
  with no `gene_name` inherits the gene of any named transcript it shares an
  exon with; transcripts matching nothing are dropped.

### Output naming and layout

Fixed by rule definitions:

```
{sample}_Aligned.out.unique.sorted.bam
{sample}_F_RP10_8.bw                              # suffix from scale_suffix
DESeq2_salmon_diff_{A}_vs_{B}.tab
DESeq2_salmon_PCA_{A}_vs_{B}[_labels].pdf
stringtie_merged_{groups}_assigned.gtf
SE.MATS.JC_sig_IDed_{A}_vs_{B}.txt
```

Only the parent directories are configurable, via `output_dir` and
`pre_processing_dir`.

### File lifetime

- **`protected()`** — `*_Aligned.out.unique.sorted.bam` and its `.bai`. Made
  read-only after creation; must be `chmod u+w` and removed by hand to
  regenerate.
- **`temp()`** — strand SAM/BAM intermediates, per-sample StringTie assemblies,
  the raw merged GTF. Deleted automatically once consumed.
- Bedgraphs and unsorted BAMs are removed by the scripts and shell blocks that
  create them.
