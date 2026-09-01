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
    """Turns Snakemake log records into `Event`s on a queue."""

    def __init__(self, sink: Callable[[Event], None]) -> None:
        super().__init__()
        self.sink = sink

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
        rule = data.get("rule_name")
        jobid = data.get("jobid")

        if name == "workflow_started":
            self.sink(Event("started"))
        elif name == "job_info":
            self.sink(Event("job_started", rule=rule, jobid=jobid))
        elif name == "job_finished":
            self.sink(Event("job_finished", rule=rule, jobid=jobid))
        elif name in ("job_error", "error"):
            self.sink(Event("job_error", rule=rule, jobid=jobid,
                            message=str(record.getMessage())[:400]))
        elif name == "progress":
            self.sink(Event("progress", done=data.get("done"), total=data.get("total")))


def run(snakefile: Path, cores: int, workdir: Path | None = None,
        dry_run: bool = False) -> Iterator[Event]:
    """Execute the workflow, yielding events as they happen.

    Snakemake runs on a worker thread so the caller — a TUI, usually — keeps
    its own loop responsive.
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
                ExecutionSettings,
                OutputSettings,
                ResourceSettings,
            )

            with SnakemakeApi(OutputSettings(printshellcmds=False, quiet={"all"})) as api:
                workflow = api.workflow(
                    resource_settings=ResourceSettings(cores=cores),
                    snakefile=snakefile,
                    workdir=workdir,
                )
                dag = workflow.dag(dag_settings=DAGSettings())
                if dry_run:
                    dag.execute_workflow(
                        execution_settings=ExecutionSettings(dryrun=True))
                else:
                    dag.execute_workflow()
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
