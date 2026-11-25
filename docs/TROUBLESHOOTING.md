# Troubleshooting

Every entry below is a failure mode encountered while developing and running
this pipeline, not a hypothetical.

## Where to look first

Snakemake's console error is usually a summary. The real message is in the
rule's own log:

```
Error in rule <name>:
    log: path/to/rule.log (check log file(s) for error message)
```

Rule logs survive failure — Snakemake deletes a failed job's *outputs* but never
its logs. If the rule log is empty or absent, the job died before reaching it;
check the scheduler log (`logs/cluster/`, or `tracejob <jobid>` on PBS).

---

## Environment and conda

### `mamba: command not found` during DAG building

The master Snakemake process is not running in an environment containing mamba.
On a cluster this usually means the job that launched Snakemake never activated
the right conda environment — an interactive `conda activate` does not propagate
to a job the scheduler starts in a fresh shell.

Wrap the invocation:

```bash
cat > run.sh << 'EOF'
#!/bin/bash
source /path/to/conda/etc/profile.d/conda.sh
conda activate snakemake
snakemake --profile <profile> --use-conda
EOF
chmod +x run.sh
```

### `command not found` for a tool that is installed

The compute-node jobscript is not activating conda. Add to
`<profile>/pbs-jobscript.sh`, before `{exec_job}`:

```sh
source /path/to/conda/etc/profile.d/conda.sh
conda activate snakemake
```

Use the path from `conda info --base` for the conda installation that owns your
environments — a cluster often has a system-wide install as well, and they are
not interchangeable.

### `ModuleNotFoundError: No module named 'pandas'` at DAG building

The environment holding Snakemake is missing pandas. `rules/common.smk` imports
it to parse the sample sheet, so this fails before any rule runs.

It appears in two places. On the login node it means the environment was built
by hand rather than from `environment_snakemake.yaml`. On a compute node it
means `pbs-jobscript.sh` is activating a *different* environment from the one
the master process uses — every cluster job re-reads the workflow and rebuilds
the DAG, so it needs pandas too.

```bash
conda env update -f environment_snakemake.yaml
```

### `version 'GLIBC_2.xx' not found`

A conda package was built against a newer glibc than the compute nodes provide.
Check the node's version with `ldd --version`.

The environments pin `sysroot_linux-64=2.17` to constrain the solver, but that
only helps when multiple builds exist. When a package ships only one build and
it is too new, pin that package to an older version — this is why `trim-galore`
is pinned to `0.6.7`.

Diagnose precisely:

```bash
ldd /path/to/env/bin/<tool> | grep "not found"
```

### `object 'colMeans' is not exported by 'namespace:MatrixGenerics'`

Two possibilities, in order of likelihood.

**R is loading packages from outside the conda environment.** `R_LIBS_USER` is
commonly exported in HPC shell profiles and points at a personal library built
against a different R version. R prepends it to `.libPaths()`, so DESeq2 loads
from conda while MatrixGenerics loads from your home directory. `Rscript
--vanilla` does *not* prevent this — it ignores `.Renviron` but still honours
`R_LIBS_USER`.

`scripts/deseq2_contrast.R` guards against this by resetting `.libPaths()` from
`CONDA_PREFIX` and printing the result at the top of its log. Check there first.

**The environment genuinely mixed Bioconductor releases.** Confirm with:

```bash
mamba list -p <env> | grep -E "^r-base|deseq2|matrixgenerics|biocgenerics"
```

Versions should share one `r44_*` build generation.

### Solve fails with "excluded by strict repo priority", or picks ancient packages

The `defaults` channels are in the solve. Channels in an environment file are
merged with `~/.condarc`, not substituted — `nodefaults` is required in the
channel list. The Snakemake, DESeq2 and rMATS environment files have it;
`environment_RNAseq.yaml` still lists `defaults`, so it is the one most likely
to hit this.

To remove the risk globally: `conda config --remove channels defaults`.

### `Download error (28) Timeout was reached`

Transient network failure while fetching a channel index. Re-run. If it recurs,
test connectivity from the node that failed:

```bash
curl -sS -o /dev/null -w '%{http_code}\n' https://conda.anaconda.org/bioconda/noarch/repodata.json
```

### Environments rebuild in every new project directory

Snakemake hashes the environment path into the env identity. Share them:

```yaml
# in <profile>/config.yaml
conda-prefix: /path/to/shared/snakemake_envs
```

---

## Cluster and filesystem

### `Missing files after N seconds ... completed successfully, but some output files are missing`

The job succeeded but its output was not yet visible to the master process.
Common on shared filesystems under load, especially for jobs that finish quickly.

```yaml
latency-wait: 120
retries: 2
```

Check the job actually ran (`tracejob <jobid>`) before assuming a code problem.

### `Directory cannot be locked`

Left over from a killed run.

```bash
snakemake --profile <profile> --use-conda --unlock
```

Confirm nothing is still running first: `qstat -u $USER`.

### `qsub: Illegal attribute or resource value`

Conflicting resource specifications in the submit command — for example
`ncpus=1` from the profile's `cluster:` line alongside `nodes=1:ppn=8` computed
by the submit script. PBS accepts one CPU syntax, not both. Remove the
redundant one from `<profile>/config.yaml`.

### Scheduler logs pile up in the working directory

PBS writes `.o`/`.e` files where the job was submitted unless told otherwise.
Either point `-o` at a directory in the profile's `cluster:` line (the path must
end in `/` and the directory must already exist), or patch `pbs-submit.py` to
default to `<cwd>/logs/cluster/`, which is portable across project directories.

---

## Pipeline logic

### `Missing wildcard values for <name>`

A checkpoint's `.get()` call was not given every wildcard in the checkpoint's
output pattern. Supply them all.

### Wildcard captures too much (e.g. `sample=RNAseq_human/SRR398210_1.fastq.gz`)

Snakemake wildcards match `/` by default. The workflow constrains
`sample="[^/]+"` globally; a new wildcard needs the same treatment.

### `IndexError: list index out of range` when formatting a shell command

Snakemake formats the *entire* shell string when logging a job, including
branches of an `if`/`else` that will not execute. A literal `{input.fq[1]}`
therefore fails for single-end samples even inside an `else` branch that never
runs. Build such arguments in Python via `params` instead of branching in the
shell.

### `A report already exists! Use -f to overwrite it`

MultiQC will not overwrite by default. The rule passes `--force`; if you have
edited it, restore that.

### `Samples in the splicing analysis disagree on library type`

rMATS takes one `--libType` for the whole run, so this is a hard stop.

Read `splicing/logs/splicing_library_info.log`. Replicates from one prep should
agree closely. **A single sample near 0.5/0.5 while its siblings sit at 0.85 is
almost never a borderline chemistry difference — check that sample's row in
`samples.tsv` first.** Pairing the wrong FASTQs (two R2 files, or R1 from one
replicate with R2 from another) destroys the read-1/read-2 strand relationship
and produces exactly this signature.

If the libraries genuinely differ, split them into separate runs.

### `Samples ... mix paired- and single-end data`

Same cause, same fix: rMATS takes one `-t`. Run them separately.

### quant.sf headers are bare transcript IDs

The Salmon index was built with `--gencode`, which strips everything after the
first `|`. The DESeq2 step needs the full header to recover gene ID and symbol.
Rebuild the index from a FASTA with intact GENCODE headers.

### A rule keeps re-running when nothing changed

Check the profile's `rerun-triggers`. `code` re-runs a rule when its definition
or script changes; `software-env` when the environment changes. Both are useful,
but if you are iterating on a script and want to skip a long branch, narrow the
list or use `--forcerun`/`--touch` deliberately.

### `protected()` output cannot be overwritten

Sorted BAMs and their indexes are read-only by design:

```bash
chmod u+w <file>
rm -f <file>
```

---

## Interpreting warnings that are not errors

| Message | Meaning |
|---|---|
| `Replace <rule> with dynamic branch` | Normal checkpoint re-evaluation |
| `Downstream jobs will be updated after completion` | Normal; a checkpoint is pending |
| `input: <TBD>` in a dry run | Depends on an unresolved checkpoint; resolves at runtime |
| `Fraction of reads failed to determine: 0.11` | Normal. Above ~0.30 suggests the transcript BED does not match the alignment genome |
| `Using an excessive number of cores ... (you asked for 8 cores)` | Trim Galore's own advice; harmless |
| `--validateMappings has no effect (deprecated)` | Salmon ≥ 1.10 makes selective alignment the default |
| Empty rule log | The tool wrote nothing to stderr. Normal for some rules on a clean run |
