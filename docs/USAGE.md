# Usage

All commands are run from the project directory — the one containing
`Snakefile`, `rules/` and `scripts/`.

## Before every run

```bash
conda activate snakemake

# what would run, and why
snakemake --profile <profile> --use-conda -n

# the same, with the actual shell commands
snakemake --profile <profile> --use-conda -n -p
```

A dry run costs seconds and catches the errors that are most expensive to hit
later: sample sheet typos, missing reference paths, design-file inconsistencies.
All of those are validated at DAG-build time and reported with the offending
sample or key named.

## Running

### On a cluster

```bash
snakemake --profile <profile> --use-conda
```

The master process dispatches one scheduler job per rule instance, so it needs
to stay alive for the duration of the run. Submit it as a job itself, or run it
under `tmux`/`screen`.

Submitting the master as a PBS job, with generous wall time and memory since it
supervises everything:

```bash
submitjob -n skm -m 100 -w 100 snakemake --profile <profile> --use-conda
```

Or wrap it, which avoids re-typing the flags and guarantees the right conda
environment on whichever node the master lands on:

```bash
cat > run_rnaseq.sh << 'EOF'
#!/bin/bash
source /path/to/conda/etc/profile.d/conda.sh
conda activate snakemake
snakemake --profile <profile> --use-conda "$@"
EOF
chmod +x run_rnaseq.sh

submitjob -n skm -m 100 -w 100 ./run_rnaseq.sh
```

### On a single machine

```bash
snakemake --use-conda --cores 24
```

Everything except the scheduler interaction works identically.

## Testing on a subset

Full runs are expensive. To validate configuration changes quickly, subsample
the FASTQs to 100,000 reads:

```bash
mkdir -p test_data
zcat /data/raw/sample_R1.fastq.gz | head -n 400000 | gzip > test_data/sample_R1.fastq.gz
zcat /data/raw/sample_R2.fastq.gz | head -n 400000 | gzip > test_data/sample_R2.fastq.gz
```

400,000 lines = 100,000 reads. `head` closes the pipe once satisfied, so `zcat`
stops early and this is fast even on very large files. R1 and R2 stay in sync
because paired FASTQs are written in matching order.

For SRA samples, cap the download instead:

```bash
fastq-dump --split-files --gzip -X 100000 SRR398210
```

Then point a copy of the sample sheet at `test_data/` and give it a separate
`pre_processing_dir` so test outputs never collide with real ones.

## Resuming and re-running

Snakemake re-runs a job when its output is missing or its inputs, code,
parameters, or software environment have changed. After an interrupted run,
resubmitting the same command picks up where it stopped.

```bash
# unlock after a killed run
snakemake --profile <profile> --use-conda --unlock

# force one rule to re-run
snakemake --profile <profile> --use-conda --forcerun deseq2_contrast

# rebuild one specific file
snakemake --profile <profile> --use-conda DGE/DESeq2_salmon_diff_nt_vs_kd_c.tab

# keep partial outputs from failed jobs, for debugging
snakemake --profile <profile> --use-conda --keep-incomplete
```

### Re-running a sample whose inputs changed

Because the sorted BAM and its index are `protected()`, they must be cleared by
hand:

```bash
chmod u+w processed/star_alignment/<sample>_Aligned.out.unique.sorted.bam*
rm -f    processed/star_alignment/<sample>_Aligned.out.unique.sorted.bam*
rm -f    processed/trimmed_reads/<sample>_*
rm -f    processed/salmon_quant/<sample>/quant.sf
rm -f    processed/bigwig/<sample>_infer_experiment.txt
```

Then dry-run to confirm only that sample is rebuilt.

## Visualising the workflow

```bash
# rule-level: one node per rule, stable across runs
snakemake --profile <profile> --use-conda --forceall --rulegraph > rulegraph.dot
dot -Tsvg rulegraph.dot > rulegraph.svg

# job-level: one node per sample instance
snakemake --profile <profile> --use-conda --dag > dag.dot
dot -Tsvg dag.dot > dag.svg
```

For a more readable rendering, [snakevision](https://github.com/OpenOmics/snakevision)
takes either `.dot` file:

```bash
snakevision -s all -o pipeline.svg rulegraph.dot
```

**A caveat on graphs and checkpoints.** The workflow uses `checkpoint`s for SRA
download and strandedness inference. Until those have actually run, Snakemake
cannot resolve what depends on them, so branches downstream appear collapsed or
missing — most visibly, `collect_and_run_multiqc` shows only one incoming edge.
Generate the graph after a completed run to see the real structure. Note that
`--rulegraph` never resolves checkpoints at all, because it works on rule
definitions rather than concrete jobs.

## Monitoring

```bash
qstat -u $USER                        # queue state
tail -f processed/logs/star_align/<sample>.log
ls -lt logs/cluster/ | head           # scheduler stdout/stderr
```

Every rule writes a `log:` under `<pre_processing_dir>/logs/<rule>/<sample>.log`.
These survive job failure — Snakemake deletes a failed job's *outputs* but never
its logs — so they are the first place to look when something goes wrong.

## Interpreting the strandedness report

Worth reading once per project:

```bash
cat splicing/logs/splicing_library_info.log
```

It reports, per sample, the sense and antisense fractions, the resulting library
type, the maximum read length, and the layout. Replicates from one library prep
should agree closely. A sample sitting near 0.5/0.5 while its siblings are at
0.85 is a signal that something is wrong with that sample — in practice, most
often a sample-sheet error pairing the wrong FASTQs — rather than a genuine
difference in chemistry.
