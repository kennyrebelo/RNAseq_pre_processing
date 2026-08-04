# Outputs

Per-sample results are written under that sample's `pre_processing_dir`.
Cross-sample results (`QC/`, `DGE/`, `splicing/`) are written relative to the
directory the workflow runs from.

```
<pre_processing_dir>/
├── trimmed_reads/
│   ├── {sample}_R1_val_1.fq.gz              # paired-end
│   ├── {sample}_R2_val_2.fq.gz
│   ├── {sample}_R{1,2}_val_{1,2}_fastqc.{html,zip}
│   ├── {sample}_trimmed.fq.gz               # single-end
│   └── {sample}_trimmed_fastqc.{html,zip}
│
├── star_alignment/
│   ├── {sample}_Aligned.out.unique.sorted.bam       ← protected
│   ├── {sample}_Aligned.out.unique.sorted.bam.bai   ← protected
│   └── {sample}_Log.final.out
│
├── salmon_quant/{sample}/
│   ├── quant.sf
│   └── logs/, aux_info/, ...
│
├── bigwig/
│   ├── {sample}_infer_experiment.txt
│   ├── {sample}_read_count.txt
│   ├── {sample}_F_RP10_8.bw                 # stranded
│   ├── {sample}_R_RP10_8.bw
│   └── {sample}_unstranded_RP10_8.bw        # unstranded
│
├── genome_browser/{sample}/
│   ├── hub.txt
│   ├── genomes.txt
│   └── {genome_build}/
│       ├── trackDb.txt
│       └── *.bw
│
└── logs/<rule>/{sample}.log
```

Trim Galore names its outputs after the input filename; the pipeline renames
them to `{sample}_*` so downstream rules have predictable inputs. The
`_Aligned.out.bam` and `_Aligned.out.unique.bam` intermediates are deleted once
the sorted BAM exists.

## QC

```
QC/
├── multiqc_report.html
├── multiqc_data/
├── logs/collect_and_run_multiqc.log
└── (copies of the STAR logs and FastQC zips it aggregated)
```

MultiQC runs with `--force`, so re-runs overwrite rather than failing.

## Genome browser tracks

Bigwigs exist in two places. `bigwig/` holds the canonical copies;
`genome_browser/{sample}/{genome_build}/` holds copies alongside the hub files.
They are copied rather than moved — moving would delete the rule's own declared
inputs and cause regeneration on every subsequent run.

Coverage is normalised as `raw_coverage × scale_factor / mapped_reads`, with
`scale_factor` defaulting to 1e8 (the "RP10_8" convention). The reverse-strand
track is negated so it plots below the axis.

Point UCSC at `genome_browser/{sample}/hub.txt`.

## Differential expression

```
DGE/
├── DESeq2_salmon_diff_{A}_vs_{B}.tab
├── DESeq2_salmon_PCA_{A}_vs_{B}.pdf
├── DESeq2_salmon_PCA_{A}_vs_{B}_labels.pdf
└── logs/deseq2_{A}_vs_{B}.log
```

One set per contrast pair.

`*_diff_*.tab` — tab-separated, unquoted, no row names:

| Column | Description |
|---|---|
| `gene_id` | e.g. `ENSG00000001036.14` |
| `gene_name` | e.g. `FUCA2` |
| `baseMean` | Mean normalised count across all samples |
| `log2FoldChange` | **contrast_B relative to contrast_A** |
| `lfcSE` | Standard error of the estimate |
| `stat` | Wald statistic |
| `pvalue` | Unadjusted |
| `padj` | Benjamini–Hochberg adjusted |

The two PCAs are the same plot; `_labels` adds sample names to each point.

The log records the resolved `.libPaths()` and loaded package versions before
anything else, which makes environment problems immediately visible.

## Splicing

```
splicing/
├── library_info.txt
├── stringtie/
│   ├── get_list_to_merge_{groups}.txt
│   └── stringtie_merged_{groups}_assigned.gtf
├── rmats/
│   ├── prep_samples_{groups}.txt
│   ├── tmp_merged_{groups}/
│   ├── prep_output_{groups}/
│   ├── bam_lists/{A}_vs_{B}_b{1,2}.txt
│   └── rmats_output_{A}_vs_{B}/
│       ├── SE.MATS.JC.txt
│       ├── SE.MATS.JC_sig.txt
│       └── A3SS/A5SS/MXE/RI outputs
├── rMATs_final_output/
│   └── SE.MATS.JC_sig_IDed_{A}_vs_{B}.txt
└── logs/<rule>/...
```

`{groups}` is the participating group names joined by `_`, e.g. `nt_kd_c`.

**`library_info.txt`** — the run-wide parameters, resolved once and read by
every downstream rule:

```
libtype=fr-firststrand
read_type=paired
read_length=151
stringtie_flag=--rf
n_samples=6
# total_nt_rep1	libtype=fr-firststrand	max_read_length=151
...
```

**`stringtie_merged_{groups}_assigned.gtf`** — the merged assembly with novel
transcripts given a gene identity of the form `<gene_id>-<gene_name>`. The
per-sample assemblies and the raw merge are `temp()` and removed once this
exists.

**`rMATs_final_output/SE.MATS.JC_sig_IDed_{A}_vs_{B}.txt`** — the primary
result. One file per contrast:

| Column | Description |
|---|---|
| `chr`, `exonStart_0base`, `exonEnd` | Skipped exon coordinates |
| `GeneID`, `geneSymbol`, `strand` | Gene identity |
| `FDR` | rMATS FDR |
| `IncLevel1`, `IncLevel2` | Per-replicate PSI for each group |
| `IncLevelDifference` | PSI difference |
| `exon_status` | `Annotated` or `Cryptic` |

`Cryptic` means the exon's coordinates do not match any exon in the reference
GTF exactly. rMATS reports `exonStart_0base` half-open while GTF is 1-based
inclusive, so +1 is applied to the start before comparison.

Only events passing `min_reads`, `fdr` and `min_delta_psi` reach this file; the
filter's log reports how many were read, kept, and skipped.

## Logs

| Location | Contents |
|---|---|
| `<pre_processing_dir>/logs/<rule>/{sample}.log` | Per-sample tool output |
| `QC/logs/`, `DGE/logs/`, `splicing/logs/` | Cross-sample rules |
| `logs/cluster/` | Scheduler stdout/stderr, if the profile is configured for it |
| `.snakemake/log/` | Snakemake's own run logs |

Rule logs persist through failure. Scheduler logs are the fallback when a job
dies before its rule log is written.
