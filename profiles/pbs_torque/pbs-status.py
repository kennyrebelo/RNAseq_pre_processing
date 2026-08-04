#!/usr/bin/env python3
"""
Report a PBS job's state back to Snakemake.

Snakemake calls this with a job id and expects exactly one of
"success", "failed" or "running" on stdout.

Two changes from the minimal version of this script that circulates widely:

  * the job_state field is located by name rather than by a fixed line
    index into `qstat -f -x` output. The field's position varies between
    PBS/OpenPBS/Torque versions and with which attributes are set, so a
    hard-coded index silently reads the wrong line on some systems.

  * qstat failures are retried before being reported. A single unreadable
    reply -- the scheduler being momentarily busy, or the job not yet
    registered right after submission -- otherwise ends the whole run.

    After MAX_ATTEMPTS the script reports "failed" rather than "running".
    That is a deliberate trade: reporting "running" forever would hang the
    workflow if the job has genuinely aged out of the scheduler's history,
    which PBS does after a configurable retention period. With `retries`
    set in the profile, a spurious "failed" costs one resubmission.
"""

import subprocess
import sys
import time

MAX_ATTEMPTS = 3
RETRY_DELAY_SECONDS = 5

jobid = sys.argv[1]


def query(jid):
    """Raw `qstat -f -x` output, or None if the query could not be made."""
    try:
        result = subprocess.run(
            ["qstat", "-f", "-x", jid],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
    except (subprocess.CalledProcessError, OSError):
        return None
    return result.stdout


def field(lines, key):
    """
    Value of a `key = value` attribute in qstat -f output, or None.

    Attributes may be indented and wrapped; this reads the first line whose
    left-hand side matches, which is sufficient for the short scalar fields
    used here.
    """
    for line in lines:
        if "=" not in line:
            continue
        name, _, value = line.partition("=")
        if name.strip() == key:
            return value.strip()
    return None


output = None
for attempt in range(1, MAX_ATTEMPTS + 1):
    output = query(jobid)
    if output is not None:
        break
    if attempt < MAX_ATTEMPTS:
        time.sleep(RETRY_DELAY_SECONDS)

if output is None:
    # Could not reach the scheduler, or the job is no longer known to it.
    print("failed")
    sys.exit(0)

lines = output.splitlines()
job_state = field(lines, "job_state")

if job_state is None:
    # Reply was readable but carried no job_state. Treat as still running so
    # a malformed single reply does not end the run; the next poll will
    # usually succeed.
    print("running")
elif job_state == "F":
    exit_status = field(lines, "Exit_status")
    if exit_status is None:
        # Finished but the exit status has not been written yet.
        print("running")
    elif exit_status == "0":
        print("success")
    else:
        print("failed")
else:
    # Q (queued), R (running), H (held), E (exiting), ...
    print("running")
