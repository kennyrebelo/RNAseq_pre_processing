# Example PBS / Torque profile

A working Snakemake profile for PBS Pro, OpenPBS and Torque. The workflow
itself contains no scheduler-specific code — everything cluster-related lives
here, so adapting the pipeline to SLURM or SGE means replacing this directory,
not editing rules.

## Contents

| File | Role |
|---|---|
| `config.yaml` | Snakemake settings: cluster command, job limits, resources |
| `pbs-submit.py` | Builds and issues the `qsub` command for each job |
| `pbs-status.py` | Reports job state back to Snakemake |
| `pbs-jobscript.sh` | Template wrapping every job; **must be edited** |

## Install

```bash
mkdir -p ~/.config/snakemake
cp -r profiles/pbs-torque ~/.config/snakemake/

chmod +x ~/.config/snakemake/pbs-torque/pbs-submit.py \
         ~/.config/snakemake/pbs-torque/pbs-status.py \
         ~/.config/snakemake/pbs-torque/pbs-jobscript.sh
```

The execute bit matters. Snakemake invokes `pbs-submit.py` and `pbs-status.py`
as commands, not via `python <script>`, so without it they fail with a
permission error.

Then:

```bash
snakemake --profile pbs-torque --use-conda
```

## Required edit

**`pbs-jobscript.sh`** contains a placeholder:

```sh
source /path/to/your/conda/etc/profile.d/conda.sh
conda activate snakemake
```

Replace the path with the output of `conda info --base` for the conda
installation that owns your environments, and the environment name if yours is
not called `snakemake`.

This is not optional. A job launched by the scheduler starts in a fresh shell
that does not inherit your interactive environment, so without these lines
every job fails with `command not found` for tools that are plainly installed.
Note that a cluster often has a system-wide conda alongside your own, and they
are not interchangeable — using the wrong one gives "environment not found"
even though `conda env list` shows it.

The environment activated here is the one holding **Snakemake itself**. The
per-rule tool environments are activated by Snakemake, not here.

## Optional edits

**Shared conda environments.** Snakemake hashes the environment path into the
environment's identity, so a new project directory rebuilds all of them. To
build once and reuse, uncomment in `config.yaml`:

```yaml
conda-prefix: /path/to/shared/snakemake_envs
```

**Queue, account or node selection.** `pbs-submit.py` already accepts `-q`,
`-A`, `-N` and the rest; add them to the `cluster:` line in `config.yaml`:

```yaml
cluster: "pbs-submit.py --depend \"{dependencies}\" -q longqueue -A my_account -l mem={resources.mem_mb}mb,walltime={resources.runtime} "
```

**Concurrency.** `jobs: 30` caps simultaneous cluster jobs. Lower it if your
site enforces a per-user limit.

## How resources are passed

Deliberately split, to avoid the two halves contradicting each other:

| Setting | Where it lives | Why |
|---|---|---|
| `mem_mb`, `threads`, `runtime` | pipeline's `config_RNAseq_pipeline.yaml` | Properties of the analysis |
| `disk_mb` | this profile's `set-resources` | Scheduler concern; not set by the pipeline |

The `cluster:` line passes `mem` and `walltime` only. CPU count is **not**
passed there: `pbs-submit.py` derives `nodes=1:ppn=<threads>` from each rule's
`threads` itself, and sending both produces a contradictory request that PBS
rejects with:

```
qsub: Illegal attribute or resource value
```

If you add `ncpus=` back to the `cluster:` line, this is what happens.

## Scheduler logs

`pbs-submit.py` sends PBS's own stdout/stderr to `<working directory>/logs/cluster/`
and creates that directory if needed, rather than scattering `.o`/`.e` files
across the directory you submitted from. The path is derived from `os.getcwd()`
at submission time, so it follows the project directory with no per-project
configuration.

To override, set `-o`/`-e` on the `cluster:` line. PBS requires that a
directory path end in `/` and that the directory already exist.

These are distinct from the pipeline's own rule logs under
`<pre_processing_dir>/logs/`, which are usually the more informative of the two.

## Settings worth keeping

```yaml
latency-wait: 120
retries: 2
```

Both address the same class of problem on a shared filesystem. A job's output
is written on a compute node and must become visible to the master process;
under load the default 5-second wait is not enough, and Snakemake reports
`completed successfully, but some output files are missing` for jobs that in
fact worked. With `retries: 0`, one such event ends a run that may be nearly
complete.

## Adapting to another scheduler

Replace this directory. Snakemake maintains profiles for common schedulers at
<https://github.com/snakemake-profiles>, including `slurm`, `sge` and `lsf`.
The only pipeline-aware part of this profile is the `set-resources` block, which
names rules — copy those lines into whichever profile you adopt.
