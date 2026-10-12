"""Run the generated workflow and report progress as structured events.

Snakemake owns execution. This module owns *observation* — it turns Snakemake's
log records into a small, stable event type so the CLI and the TUI can both
show progress without either knowing anything about Snakemake's internals.

Snakemake 9 exposes structured fields on its log records (`jobid`, `rule_name`,
`done`, `total`), so nothing here parses text, and no fork of Snakemake is
needed to get at them.
"""

from __future__ import annotations

import logging
import queue
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator


# What the SLURM executor plugin prints once, at info level, when it starts:
# `SLURM run ID: <uuid>`. Matched as a prefix rather than parsed, and kept here
# rather than inline because it is a string in someone else's package — if a
# plugin release changes the wording, cancellation goes quiet, and this is the
# line a test pins.
SLURM_RUN_ID = "SLURM run ID:"

# The thread `run()` puts Snakemake on, named so that leaving can find it.
THREAD = "comparem2-snakemake"

# How often a reader with nothing to read checks whether it has been told to
# stop. A job can run for an hour without a log line, and a quit cannot wait
# for one.
STOP_POLL = 0.25

# The first Snakemake whose log records say which event they are. Before it
# `snakemake.logging.get_event` does not exist (9.13.7 has none), `_Capture`
# drops every record, and a run under the interface shows no job at all.
EVENTS_SINCE = "9.14.0"


def events_unavailable() -> str | None:
    """Why this Snakemake cannot be followed by the interface, or None."""
    try:
        from snakemake.logging import get_event  # noqa: F401
    except ImportError:
        from importlib.metadata import version

        return (f"--tui needs Snakemake {EVENTS_SINCE} or newer to follow a "
                f"run, and this is {version('snakemake')}. The same command "
                "without --tui works.")
    return None


def stop_scheduling() -> bool:
    """Tell a running Snakemake to start nothing further. True if one was told.

    `JobScheduler.exit_gracefully` is Snakemake's own SIGTERM handler: the jobs
    running carry on, nothing new is scheduled, and when the last one ends
    Snakemake shuts down and takes its lock with it. It is never installed
    here, because the scheduler is not on the main thread (see `cancel.py`),
    so it is called directly. The scheduler is found among live objects
    because neither route into Snakemake hands one back: the profile branch
    gets a bool from `args_to_api`.

    Measured on 2026-10-12: without this, both of the quit dialog's options
    let a new bakta start after the key that promised nothing would.
    """
    import gc

    try:
        from snakemake.scheduling.job_scheduler import JobScheduler
    except ImportError:  # pragma: no cover - checked by events_unavailable
        return False
    told = False
    for obj in gc.get_objects():
        # `type()`, not `isinstance()`: JobScheduler derives from an ABC, whose
        # check reads `__class__`, and some objects in a live process lie
        # about theirs — `issubclass() arg 1 must be a class`.
        if issubclass(type(obj), JobScheduler):
            obj.exit_gracefully()
            told = True
    return told


def snakemake_alive(wait: float = 0.0) -> bool:
    """Whether `run()`'s Snakemake thread is still going, after up to `wait` s."""
    for thread in threading.enumerate():
        if thread.name == THREAD:
            thread.join(wait)
            return thread.is_alive()
    return False


@dataclass(frozen=True)
class Event:
    """One thing that happened during a run."""

    kind: str  # started | job_started | job_finished | job_error | progress
    #          | slurm_run_id | done | error
    rule: str | None = None
    jobid: int | None = None
    done: int | None = None
    total: int | None = None
    message: str = ""


class _Capture(logging.Handler):
    """Turns Snakemake log records into `Event`s on a queue.

    Snakemake does not spell a job the same way twice, verified against 9.26.1
    by dumping the records: `job_info` carries `jobid` and `rule_name`,
    `job_finished` carries `job_id` — with an underscore — and no rule name at
    all, and `job_started` carries only a list of ids. So `job_info` is the one
    place a rule name can be had, and the mapping is remembered here to give
    the finish event its name back.

    Reading `jobid` on a `job_finished` record therefore yielded None, the TUI
    dropped every finish event as unattributable, and a workflow that completed
    6 of 6 steps reported that nothing had run.
    """

    def __init__(self, sink: Callable[[Event], None]) -> None:
        super().__init__()
        self.sink = sink
        self.rules: dict[int, str] = {}

    def emit(self, record: logging.LogRecord) -> None:
        try:
            message = record.getMessage()
        except Exception:  # pragma: no cover - a handler must not raise
            message = ""
        # The name every job of this run is submitted under. It is the only
        # handle on a queue full of them — `scancel --name` stops the run in
        # one call, whatever is still pending — and the plugin sets it *in
        # order to* make that possible. It arrives as a plain info message with
        # no structured event attached, so it has to be read before the
        # get_event() filter below drops it.
        if message.startswith(SLURM_RUN_ID):
            self.sink(Event("slurm_run_id",
                            message=message[len(SLURM_RUN_ID):].strip()))

        try:
            from snakemake.logging import get_event

            event = get_event(record)
        except Exception:  # pragma: no cover - snakemake internals
            return
        if event is None:
            return

        name = getattr(event, "value", str(event))
        data = record.__dict__
        jobid = data.get("jobid")
        if jobid is None:
            jobid = data.get("job_id")
        rule = data.get("rule_name") or self.rules.get(jobid)

        if name == "workflow_started":
            self.sink(Event("started"))
        elif name == "job_info":
            if jobid is not None and rule:
                self.rules[jobid] = rule
            self.sink(Event("job_started", rule=rule, jobid=jobid))
        elif name == "job_finished":
            self.sink(Event("job_finished", rule=rule, jobid=jobid))
        elif name in ("job_error", "error"):
            self.sink(Event("job_error", rule=rule, jobid=jobid,
                            message=str(record.getMessage())[:400]))
        elif name == "progress":
            self.sink(Event("progress", done=data.get("done"), total=data.get("total")))


def _profile_argv(snakefile: Path, profile: str, cores: int | None,
                  workdir: Path | None, dry_run: bool, keep_going: bool,
                  rerun_incomplete: bool, conda_prefix: Path | None,
                  deploy: bool = True) -> list[str]:
    """The command line the profile branch hands to Snakemake's own parser.

    Kept separate so it can be asserted in a test without a Snakemake install,
    which is the only way this gets checked on a machine that cannot run the
    workflow at all.
    """
    argv = ["--snakefile", str(snakefile),
            # Under a TUI, Snakemake's own logger plugin writing to stderr
            # scribbles over the Textual display. The non-profile branch gets
            # the same effect from `OutputSettings(quiet={"all"})`.
            "--quiet", "all",
            "--profile", profile]
    if workdir is not None:
        argv += ["--directory", str(workdir)]
    # Two settings, as in the API branch: whether to deploy, and where. A
    # missing prefix means Snakemake's default location, not no deployment.
    if deploy:
        argv += ["--software-deployment-method", "conda"]
    if conda_prefix is not None:
        argv += ["--conda-prefix", str(conda_prefix)]
    # Omitted unless asked for, so a profile's own `cores:` survives. It does
    # not affect how many jobs reach the queue — see cli.py for the measurement.
    if cores is not None:
        argv += ["--cores", str(cores)]
    if rerun_incomplete:
        argv.append("--rerun-incomplete")
    if keep_going:
        argv.append("--keep-going")
    if dry_run:
        argv.append("--dry-run")
    return argv


# A queued job is retried once unless the profile says otherwise. On GenomeDK
# on 2026-10-09 one bakta job died on one node, where conda itself failed to
# start; 42 jobs on 9 other nodes were fine. Without a retry that one job cost
# panaroo, snp-dists and FastTree for every genome. Local runs get no retry:
# there, a failure is the input's, and a retry only doubles the wait for it.
QUEUE_RETRIES = 1


def queue_retries(parsed) -> None:
    """Give a profile run one retry, unless its profile set its own.

    Applied after Snakemake's parser, which has already read the profile, so a
    profile's `retries: 3` stands. The one case this cannot see is a profile
    that sets `retries: 0` on purpose: it parses the same as no setting.
    """
    if not parsed.retries:
        parsed.retries = QUEUE_RETRIES


def queue_retries_flag(argv: list[str]) -> list[str]:
    """`queue_retries()` for the CLI, which hands Snakemake a command line.

    The same decision, made by the same parser over the same argv, so the two
    paths cannot disagree about what a profile said. Until 2026-10-12 only the
    TUI's path applied it, and a failed SLURM job under a plain `comparem2`
    was not retried while `--help` promised it would be.

    A command line Snakemake's parser rejects gets nothing added, and the
    parser's complaint is swallowed here: the run about to start makes the same
    complaint, once, in Snakemake's own words.
    """
    import contextlib
    import io

    from snakemake.cli import parse_args

    try:
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            _, parsed = parse_args(argv)
    except SystemExit:
        return []
    if parsed.retries:
        return []
    queue_retries(parsed)
    return ["--retries", str(parsed.retries)]


def run(snakefile: Path, cores: int | None, workdir: Path | None = None,
        dry_run: bool = False, keep_going: bool = False,
        rerun_incomplete: bool = True,
        conda_prefix: Path | None = None,
        deploy: bool = True,
        profile: str | None = None,
        stop: threading.Event | None = None) -> Iterator[Event]:
    """Execute the workflow, yielding events as they happen.

    Snakemake runs on a worker thread so the caller — a TUI, usually — keeps
    its own loop responsive.

    **`profile` takes a second route through Snakemake, deliberately.** A
    profile is a construct of Snakemake's *command line*: `config.yaml` is read
    as argparse defaults, which is how one file can set the executor, the
    account, `jobs`, `default-resources` and `set-resources` at once. The API
    used below has no notion of it — `execute_workflow()` takes `executor` and
    `executor_settings`, and mapping a profile onto those by hand would mean
    re-implementing that parser and getting a subset of it right. So with a
    profile, this calls `snakemake.cli.parse_args()` and `args_to_api()`, which
    is Snakemake's own CLI running in this process.

    The cost is that `args_to_api` swallows the exception and reports failure
    as `False`, where the API branch below yields the exception text as an
    `error` event. Both still emit every job event, because those come from the
    log handler and not from either call. Keeping the API branch for the local
    case is what makes the profile branch additive: the path that has actually
    been run stays byte-identical.

    The other cost is that a *command line* can be rejected, and Snakemake
    rejects one by calling `exit()` — which is why `SystemExit` is caught
    below, separately and before `Exception`, since it is not one.

    `workdir` is not optional in practice: Snakemake locks its working
    directory rather than its output paths, so leaving it unset makes every run
    in one checkout share `./.snakemake`, and a killed run leaves a lock the
    next cannot clear. The CLI passes `--directory` for the same reason.

    Conda deployment is on by default and nothing in the pipeline turns it off:
    the package ships the pipeline only, and each rule's tools come from the env
    files `prepare()` wrote. `conda_prefix` decides where those environments
    land, and it wants to be the same directory on every run — see
    `cli.default_conda_prefix()`.

    `deploy=False` exists for one caller: the test that drives a hand-written
    Snakefile with no `conda:` directives, to check the runner's event field
    names against real Snakemake. Enabling deployment makes Snakemake require
    conda even when no rule asks for an environment, and CI has no conda. No
    production path passes it, and a generated Snakefile has a directive on
    every rule.

    `stop`, once set, ends the stream within `STOP_POLL` whether or not an
    event arrives, and leaves Snakemake's thread to whoever set it. The TUI
    sets it on quitting: the stream is read on a thread `App.run()` joins on
    the way out, and waiting for the next event there was waiting for the
    running job to finish — 58 s on thylakoid on 2026-10-12.
    """
    events: queue.Queue[Event | None] = queue.Queue()
    handler = _Capture(events.put)
    logger = logging.getLogger("snakemake")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    def work() -> None:
        try:
            if profile is not None:
                from snakemake.cli import args_to_api, parse_args

                argv = _profile_argv(snakefile, profile, cores, workdir,
                                     dry_run, keep_going, rerun_incomplete,
                                     conda_prefix, deploy)
                parser, parsed = parse_args(argv)
                queue_retries(parsed)
                if args_to_api(parsed, parser):
                    events.put(Event("done"))
                else:
                    # `args_to_api` has already printed the exception through
                    # Snakemake's own printer; it returns a bool and nothing
                    # more, so there is no text to pass on here.
                    events.put(Event("error", message=(
                        "Snakemake reported a failed run — see the log under "
                        f"{workdir}/.snakemake/log/ for what the queue said")))
                return

            from snakemake.api import SnakemakeApi
            from snakemake.settings.types import (
                DAGSettings,
                DeploymentSettings,
                ExecutionSettings,
                OutputSettings,
                ResourceSettings,
            )
            # Not a snakemake.settings export: the enum lives in the executor
            # plugin interface package, which settings/types.py imports it from.
            from snakemake_interface_executor_plugins.settings import DeploymentMethod

            deployment = DeploymentSettings(
                deployment_method={DeploymentMethod.CONDA},
                conda_prefix=conda_prefix,
            ) if deploy else None

            with SnakemakeApi(OutputSettings(printshellcmds=False, quiet={"all"})) as api:
                workflow = api.workflow(
                    resource_settings=ResourceSettings(cores=cores),
                    deployment_settings=deployment,
                    snakefile=snakefile,
                    workdir=workdir,
                )
                dag = workflow.dag(
                    dag_settings=DAGSettings(force_incomplete=rerun_incomplete))
                # A dry run is an executor plugin, not a setting: passing
                # `dryrun=True` to ExecutionSettings raises TypeError, which
                # the handler below would have reported as a failed run.
                dag.execute_workflow(
                    executor="dryrun" if dry_run else "local",
                    execution_settings=ExecutionSettings(keep_going=keep_going),
                )
            events.put(Event("done"))
        except SystemExit as exc:
            # Snakemake's own CLI refuses a command line by printing a message
            # and calling `exit()`. That raises SystemExit, which is not an
            # Exception, so the handler below never saw it and `threading`
            # discarded it without a word — leaving a stream with no `done` and
            # no `error` at all.
            #
            # Which the TUI then read as *success*. In a fresh output directory
            # it said "Nothing ran"; in one holding earlier results it said
            # every selected tool was already up to date and wrote a report,
            # over a run Snakemake had refused to start. Measured with
            # `--profile <a directory holding no config.yaml>`: zero events.
            #
            # Only the profile branch can get here, and `--profile` is the only
            # part of that argv not built above, so the message says so rather
            # than guessing. Snakemake's own message names every directory it
            # searched, but it goes to stderr, which is not where a TUI user is
            # looking.
            blame = (f"the only part of it that came from you is "
                     f"--profile {profile!r}" if profile is not None else
                     "its own message went to stderr")
            events.put(Event("error", message=(
                f"Snakemake rejected the command line (exit {exc.code}) — "
                f"{blame}. Running the same command without --tui shows what "
                "it said.")))
        except Exception as exc:  # surfaced to the user, not swallowed
            events.put(Event("error", message=f"{type(exc).__name__}: {exc}"))
        finally:
            events.put(None)

    thread = threading.Thread(target=work, daemon=True, name=THREAD)
    thread.start()
    try:
        while stop is None or not stop.is_set():
            try:
                event = events.get(timeout=STOP_POLL)
            except queue.Empty:
                continue
            if event is None:
                break
            yield event
    finally:
        logger.removeHandler(handler)
        if stop is None or not stop.is_set():
            thread.join(timeout=5)
