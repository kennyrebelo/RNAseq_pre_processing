# Installation

## Requirements

### System

| Requirement | Notes |
|---|---|
| Linux x86-64 | glibc 2.17 or newer; the environments pin `sysroot_linux-64=2.17` for this reason |
| conda or miniconda | Any recent version. `mamba` is strongly recommended as the solver |
| Job scheduler (optional) | PBS/Torque, SLURM, SGE — via a Snakemake profile. The workflow also runs on a single machine |
| Disk | Roughly 30–60 GB per sample for intermediates, plus ~15 GB for the three conda environments |
| Memory | 40 GB for the STAR and Salmon steps (configurable). A human STAR index needs ~32 GB resident |
| Network | Required on first run to build conda environments, and for SRA samples |

### Software installed by you

Only two things:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake mamba
conda activate snakemake
```

Everything else is declared in the three environment files and built
automatically by `--use-conda`.

Tested with Snakemake 7.x. The workflow uses `checkpoint`, `temp()`,
`protected()`, `directory()` and `script:` directives, all long-standing
features.

### Software installed automatically

**`environment_RNAseq.yaml`** — the main toolchain:

| Package | Pinned | Purpose |
|---|---|---|
| `trim-galore` | 0.6.7 | Adapter and quality trimming |
| `fastqc` | — | Per-sample read QC |
| `star` | — | Genome alignment |
| `salmon` | 2.3.4 | Transcript quantification |
| `samtools` | — | Filtering, sorting, indexing, read counts |
| `multiqc` | — | Aggregate QC report |
| `sra-tools` | — | `fastq-dump` for SRA samples |
| `rseqc` | — | `infer_experiment.py` strandedness inference |
| `pysam` | — | Strand splitting for bigwigs |
| `bedtools` | — | `genomeCoverageBed` for coverage |
| `ucsc-bedgraphtobigwig` | — | Fallback if the configured binary is unavailable |

`trim-galore` is pinned to 0.6.7 because later bioconda builds have been
compiled against a newer glibc than some cluster nodes provide.

**`environment_DESeq2.yaml`** — differential expression only:

`r-base=4.4`, `bioconductor-deseq2`, `bioconductor-tximport`, `r-ggplot2`,
`r-readr`.

Only R is pinned. Bioconda tags every R package build with the R version it was
compiled against (`r44_*`), so pinning R constrains the whole Bioconductor
stack to one coherent generation without naming each package. Pinning
individual package versions instead tends to produce unsatisfiable solves.

**`environment_rMATS.yaml`** — splicing only:

`python=3.10`, `stringtie=2.2.1`, `rmats=4.3.0`, `samtools`.

rMATS is pinned because `scripts/filter_se_events.py` assumes the 4.3.0
`SE.MATS.JC.txt` column layout.

### Why `nodefaults` matters

All three environment files list:

```yaml
channels:
  - conda-forge
  - bioconda
  - nodefaults
```

Channels in an environment file are **merged** with those in `~/.condarc`, not
substituted for them. Omitting `defaults` is therefore not enough to keep
`pkgs/main` and `pkgs/r` out of the solve — `nodefaults` is required. Those
channels carry very old Bioconductor builds, and a solver free to mix them will
happily pair a current DESeq2 with a MatrixGenerics four releases too old,
producing:

```
object 'colMeans' is not exported by 'namespace:MatrixGenerics'
```

## Reference data

Not distributed with the pipeline. You need, per species:

| File | Used by | Required when |
|---|---|---|
| STAR genome index | `star_align` | Always |
| Salmon index | `salmon_quant` | Always |
| Transcript BED (12-column) | `infer_experiment.py` | bigwig or splicing enabled |
| Chromosome sizes | `bedGraphToBigWig` | bigwig enabled |
| GENCODE GTF | StringTie, rMATS, cryptic-exon calling | splicing enabled |

Paths go in `config_RNAseq_pipeline.yaml` under `references:`.

Two constraints worth knowing before you build indices:

- **The Salmon index must retain full GENCODE FASTA headers.** The DESeq2 step
  parses gene ID and gene symbol out of `quant.sf`'s `Name` column
  (`ENST…|ENSG…|…|SYMBOL|…`). An index built with `salmon index --gencode`
  strips everything after the first `|` and the step will stop with an explicit
  message.
- **The transcript BED must match the genome build the BAMs are aligned to.**
  A mismatch shows up as a large "fraction of reads failed to determine" in the
  strandedness output and makes the library-type call unreliable.

## Cluster profile

The workflow contains no scheduler-specific code. Cluster behaviour comes from
a Snakemake profile, typically at `~/.config/snakemake/<profile-name>/`:

```
<profile-name>/
├── config.yaml          # cluster command, resources, job limits
├── pbs-submit.py        # builds and issues qsub
├── pbs-status.py        # reports job state back to Snakemake
└── pbs-jobscript.sh     # template wrapping each job
```

Settings that matter:

```yaml
jobs: 30                 # max concurrent cluster jobs
latency-wait: 120        # seconds to wait for outputs on a shared filesystem
retries: 2               # retry a failed job before giving up
keep-going: true
use-conda: true
conda-frontend: mamba
rerun-incomplete: true
conda-prefix: /path/to/shared/envs   # optional: reuse envs across directories
```

`latency-wait` and `retries` are worth setting generously. On a busy shared
filesystem the default 5-second wait is not enough for an output written on a
compute node to become visible to the master process, and a single spurious
failure with `retries: 0` takes down an otherwise complete run.

**The compute-node jobscript must activate conda itself.** A job launched by
the scheduler gets a fresh shell that does not inherit your interactive
environment, so `pbs-jobscript.sh` needs:

```sh
#!/bin/sh
# properties = {properties}
source /path/to/conda/etc/profile.d/conda.sh
conda activate snakemake
{exec_job}
```

Without this, jobs fail with `command not found` for tools that are plainly
installed.

## Verifying the install

```bash
conda activate snakemake
snakemake --version

# build the environments without running any analysis
snakemake --profile <your-profile> --use-conda --conda-create-envs-only

# confirm the DAG resolves against your sample sheet
snakemake --profile <your-profile> --use-conda -n
```

Building the environments separately is worth the few extra minutes: a solver
problem surfaces immediately instead of 20 jobs into a run.
