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


@dataclass(frozen=True)
class Event:
    """One thing that happened during a run."""

    kind: str  # started | job_started | job_finished | job_error | progress | done | error
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


def run(snakefile: Path, cores: int, workdir: Path | None = None,
        dry_run: bool = False, keep_going: bool = False,
        rerun_incomplete: bool = True,
        conda_prefix: Path | None = None,
        deploy: bool = True) -> Iterator[Event]:
    """Execute the workflow, yielding events as they happen.

    Snakemake runs on a worker thread so the caller — a TUI, usually — keeps
    its own loop responsive.

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
    """
    events: queue.Queue[Event | None] = queue.Queue()
    handler = _Capture(events.put)
    logger = logging.getLogger("snakemake")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    def work() -> None:
        try:
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
        except Exception as exc:  # surfaced to the user, not swallowed
            events.put(Event("error", message=f"{type(exc).__name__}: {exc}"))
        finally:
            events.put(None)

    thread = threading.Thread(target=work, daemon=True, name="snakemake")
    thread.start()
    try:
        while True:
            event = events.get()
            if event is None:
                break
            yield event
    finally:
        logger.removeHandler(handler)
        thread.join(timeout=5)
