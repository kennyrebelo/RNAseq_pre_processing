#!/bin/sh
# properties = {properties}

# ---------------------------------------------------------------------------
# EDIT THIS PATH before first use.
#
# A job launched by the scheduler starts in a fresh shell that does NOT
# inherit the environment you had when you submitted. Without these two
# lines, every job fails with "command not found" for tools that are
# plainly installed.
#
# Use the base reported by `conda info --base` for the conda installation
# that owns your environments. A cluster often has a system-wide conda as
# well, and the two are not interchangeable.
#
# The environment activated here is the one containing Snakemake itself,
# not the per-rule tool environments -- those are activated by Snakemake.
# ---------------------------------------------------------------------------
source /path/to/your/conda/etc/profile.d/conda.sh
conda activate snakemake

{exec_job}
