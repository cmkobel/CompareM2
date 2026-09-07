"""Stopping a run that is already going.

**Snakemake will not do this for us**, and the reason is worth writing down
because it looks like it should.

Its scheduler reaches `executor.cancel()` from exactly one place — a
`KeyboardInterrupt` inside its own scheduling loop (`job_scheduler.py`) — and it
installs the SIGTERM handler that would get it there inside a
`try/except ValueError` that silently skips when the scheduler is not on the
main thread. `runner.run()` puts Snakemake on a worker thread so the interface
keeps its own loop, which means that handler is never installed and there is
nothing to signal. Reaching into the running workflow for the executor object is
no better: the profile branch goes through `args_to_api()` and hands back a
bool.

So each half is done here, directly, and the two halves are not symmetric:

  - **A queue.** The SLURM executor plugin submits every job with
    `--job-name <run_uuid>` *in order to* make `--name`-based cancellation
    possible — its own comment says so — and announces the id as
    `SLURM run ID: <uuid>`. One `scancel` therefore stops the whole run
    whatever is still queued, and it is the same command the plugin's own
    `cancel_slurm_jobs()` would run.
  - **This machine.** Snakemake's local executor `cancel()` is
    `self.pool.shutdown()`: it stops *scheduling* and waits for what is
    running. It does not kill anything. So the process tree is walked and
    signalled here.

Both leave partial output files behind. That is what `--rerun-incomplete`
handles, which `runner.run()` passes by default, so the next run redoes a rule
that was interrupted rather than trusting half a result. Neither releases
Snakemake's lock on the output directory — nothing here can, because the lock
outlives the process — so a cancelled run is followed by `--unlock`, or `u` in
the interface.
"""

from __future__ import annotations

import getpass
import os
import signal
import subprocess
from time import sleep


def descendants(pid: int) -> list[int] | None:
    """Every live process below `pid`, deepest first — or None if unknowable.

    `None` means `ps` could not be read, which is not the same answer as the
    empty list and must not be reported as one: the caller says "nothing was
    running" for `[]`, and that would be a lie about a tree it failed to see.

    Deepest first because that is the order they have to be signalled in: a
    tool runs under a `conda run` wrapper under a spawned Snakemake, and
    killing the wrapper first leaves the tool running with a reparented
    parent and no one to collect it.

    **A zombie is not a live process and is not in the list.** One that has
    exited but whose parent has not yet reaped it keeps its place in the
    process table, and nothing can be done to it — SIGKILL at a zombie returns
    successfully and changes nothing. Counting them made a cancelled run
    report a SIGTERM that had worked as having "needed SIGKILL", and turned
    the CI Linux runners red on 2026-09-07 while the same suite passed on
    macOS: `pgrep -P`, which this used to walk with, lists a zombie on Linux
    and does not on macOS. `ps` lists it on both, so the state has to be read
    either way and the platform difference stops mattering. They are still
    walked *through*, in case a table lags — though a dying process's children
    are reparented at once, so nothing should be hidden below one.

    One `ps` snapshot rather than one `pgrep` per node: it is a single process
    either way, and the tree cannot shift underneath the walk. psutil would do
    this too and is not worth a dependency of the *pipeline* for one function;
    `ps` is on both platforms the pipeline runs on, and on Linux it comes from
    the same `procps` package `pgrep` did, so nothing new has to be present.
    """
    try:
        out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,stat="],
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        # `subprocess.run` does not raise on a non-zero exit, so without this
        # a `ps` that refused the arguments — a busybox build, a `hidepid`
        # mount — reads as an empty process table, and the caller tells the
        # user nothing was running while the tools carry on. Not knowing and
        # knowing there is nothing are different answers, so they get
        # different return values.
        return None

    children: dict[int, list[int]] = {}
    reaped: set[int] = set()
    for line in out.stdout.splitlines():
        fields = line.split()
        if len(fields) < 3 or not (fields[0].isdigit() and fields[1].isdigit()):
            continue
        child, parent, state = int(fields[0]), int(fields[1]), fields[2]
        if child == parent:
            # A process listed as its own parent would be walked forever.
            continue
        children.setdefault(parent, []).append(child)
        if state.startswith("Z"):
            reaped.add(child)

    found: list[int] = []
    frontier = [pid]
    while frontier:
        for child in children.get(frontier.pop(), ()):
            if child not in reaped:
                found.append(child)
            frontier.append(child)
    found.reverse()
    return found


def stop_local(pid: int | None = None, grace: float = 3.0) -> str:
    """SIGTERM every descendant of this process, then SIGKILL what is left.

    Returns a line to print, because the caller is a TUI that has already given
    the terminal back and the outcome is the last thing the user sees.

    The grace period is not politeness: a tool killed outright leaves its
    output half-written *and* its temporary files behind, and several of the
    fourteen clean up on SIGTERM. Three seconds is enough for that and short
    enough that nobody wonders whether the quit worked.
    """
    root = os.getpid() if pid is None else pid
    targets = descendants(root)
    if targets is None:
        return ("The process table could not be read, so nothing was "
                "signalled. Check with: ps -f -u $USER")
    if not targets:
        return "No job processes were running."

    signalled = set()
    for target in targets:
        try:
            os.kill(target, signal.SIGTERM)
            signalled.add(target)
        except ProcessLookupError:
            # Finished between listing and signalling. Nothing to report: the
            # outcome the user asked for is the one that happened.
            continue
        except PermissionError:
            continue

    sleep(grace)
    # Re-scanned rather than re-checking the list that was signalled.
    # Snakemake is on a daemon thread that only stops when it next calls into
    # the interface — measured at about a second after the app comes down — so
    # it can start one more job in the window between the signal and this
    # wait, and a job that started *during* the cancelling is exactly the
    # survivor the user would notice. Whatever is still below the root now gets
    # SIGKILL, whether it ignored the SIGTERM or never received one.
    killed = set()
    for target in descendants(root) or ():
        try:
            os.kill(target, signal.SIGKILL)
            killed.add(target)
        except OSError:
            continue

    # The union, not the first pass's count. The two scans see different sets —
    # a job Snakemake started during the grace period is in the second and not
    # the first, and one that exited between listing and signalling is in
    # neither — so counting only the SIGTERMs could report "Stopped 0 job
    # processes (2 needed SIGKILL)", which says two contradictory things about
    # the same run.
    stopped = signalled | killed
    plural = "process" if len(stopped) == 1 else "processes"
    note = f"Stopped {len(stopped)} job {plural}"
    if killed:
        note += f" ({len(killed)} needed SIGKILL)"
    return note + "."


def stop_slurm(run_uuid: str, user: str | None = None) -> str:
    """`scancel` every job of this run, by the job name the plugin gave them.

    One call regardless of how many jobs are in the queue, and it covers jobs
    submitted but not yet started — which is the case that matters most, since
    a killed frontend otherwise leaves a queue full of work nothing will
    collect.

    `-u` is included so a name collision with another user's job produces
    nothing rather than a permission error, and because `--me` is not in every
    SLURM old enough to still be running on a cluster.
    """
    who = getpass.getuser() if user is None else user
    try:
        done = subprocess.run(
            ["scancel", "--name", run_uuid, "-u", who],
            capture_output=True, text=True, timeout=60)
    except FileNotFoundError:
        return ("scancel is not on PATH here, so the queue was left alone — "
                f"cancel the run with: scancel --name {run_uuid} -u {who}")
    except subprocess.TimeoutExpired:
        return (f"scancel did not answer within a minute. Check with: "
                f"squeue --name {run_uuid} -u {who}")
    if done.returncode != 0:
        detail = (done.stderr or done.stdout).strip().splitlines()
        first = detail[0] if detail else f"exit code {done.returncode}"
        return f"scancel failed ({first}). Jobs may still be queued."
    # scancel says nothing on success and exits 0 on a name that matches
    # nothing, so this is a statement about what was asked, not about what
    # was in the queue.
    return f"Asked SLURM to cancel run {run_uuid}."
