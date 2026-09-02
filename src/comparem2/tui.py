"""Keyboard-driven interface.

The point is that a non-bioinformatician should be able to see what will run,
what it will cost to install, and how far along it is, without reading a manual
or a scrolling wall of Snakemake output.

Tool state comes from `catalogue.py` and progress from `runner.Event`, so this
module knows nothing about Snakemake and nothing about which tools exist.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widgets import DataTable, Footer, Header, ProgressBar, RichLog, Static

from .catalogue import CATALOGUE
from .report import render_report
from .runner import Event, run
from .snakefile import prepare

PENDING, RUNNING, DONE, FAILED, SKIPPED, NOT_RUN = "·", "▸", "✓", "✗", "–", "○"
LABEL = {PENDING: "pending", RUNNING: "running", DONE: "done",
         FAILED: "failed", SKIPPED: "not selected", NOT_RUN: "not run"}

# The selection column. These were `[x]`, `[+]` and `[ ]` and rendered as
# nothing at all: a DataTable cell given a `str` is parsed as Rich markup, and
# `[x]` is a tag, not text. The selection UI had no visible selection.
MARK_ON, MARK_DEP, MARK_OFF = "▣", "▨", "▢"


class ComparemTUI(App):
    """Pick tools, watch them run, open the report."""

    CSS = """
    Screen { layout: vertical; }
    #cost { padding: 0 1; color: $text-muted; }
    #panes { height: 1fr; }
    DataTable { width: 46%; border: round $primary; }
    RichLog { width: 1fr; border: round $primary; padding: 0 1; }
    ProgressBar { padding: 0 1; }
    """

    BINDINGS = [
        ("space", "toggle", "Select/deselect"),
        ("a", "all", "Select all"),
        ("n", "none", "Select none"),
        ("r", "start", "Run"),
        ("q", "quit", "Quit"),
    ]

    def __init__(self, inputs: list[Path], workdir: Path, databases: Path,
                 samples: tuple[str, ...], cores: int,
                 selected: list[str] | None = None,
                 overrides: dict[str, tuple[tuple[str, str], ...]] | None = None,
                 launcher: Sequence[str] | None = None,
                 keep_going: bool = False) -> None:
        super().__init__()
        self.inputs = inputs
        self.workdir = workdir
        self.databases = databases
        self.samples = samples
        self.cores = cores
        self.overrides = overrides
        self.launcher = launcher
        self.keep_going = keep_going
        # Seeded from `--until` when given. Selecting everything by default
        # puts gtdbtk's 141.4 GB one keypress away, so a user who named the
        # tools they want on the command line gets exactly those.
        self.selected: set[str] = set(selected) if selected else {t.name for t in CATALOGUE}
        self.state: dict[str, str] = {t.name: PENDING for t in CATALOGUE}
        self.running = False
        self.cost_text = ""

    # --- layout ----------------------------------------------------

    def compose(self) -> ComposeResult:
        yield Header(show_clock=False)
        yield Static(id="cost")
        with Horizontal(id="panes"):
            yield DataTable(cursor_type="row", zebra_stripes=True)
            yield RichLog(highlight=False, markup=True, wrap=True)
        yield ProgressBar(total=100, show_eta=False)
        yield Footer()

    def on_mount(self) -> None:
        self.title = f"CompareM2 v3 — {len(self.samples)} assemblies"
        table = self.query_one(DataTable)
        # Explicit column keys: update_cell() matches on the key, not the label.
        table.add_column(" ", key="sel", width=3)
        table.add_column("Tool", key="tool", width=12)
        table.add_column("Status", key="status", width=13)
        table.add_column("What it does", key="summary")
        for tool in CATALOGUE:
            table.add_row(MARK_ON, tool.name, LABEL[PENDING], tool.summary, key=tool.name)
        # Not refresh_cost(): the selection can differ from "everything" before
        # a key has been pressed, because `--until` seeds it, and the rows have
        # to say so.
        self.sync_table()
        self.query_one(RichLog).write(
            "[dim]space[/] select · [dim]r[/] run · [dim]q[/] quit")

    # --- selection -------------------------------------------------

    def refresh_cost(self) -> None:
        chosen = sorted(self.selected)
        closure = CATALOGUE.closure(chosen) if chosen else []
        known = CATALOGUE.install_size(chosen) if chosen else 0
        unknown = CATALOGUE.unmeasured(chosen) if chosen else []
        cost = f"{known / 1e9:.1f} GB" if known else "no databases"
        if unknown:
            cost += f" + {len(unknown)} of unknown size ({', '.join(d.name for d in unknown)})"
        pulled = len(closure) - len(chosen)
        extra = f", {pulled} pulled in as dependencies" if pulled > 0 else ""
        self.cost_text = (
            f"{len(closure)} tools selected{extra} — databases to download: {cost}")
        self.query_one("#cost", Static).update(self.cost_text)

    def _row_key(self) -> str | None:
        table = self.query_one(DataTable)
        if table.cursor_row < 0:
            return None
        return str(table.get_row_at(table.cursor_row)[1])

    def action_toggle(self) -> None:
        name = self._row_key()
        if name is None or self.running:
            return
        self.selected.symmetric_difference_update({name})
        self.sync_table()

    def action_all(self) -> None:
        if not self.running:
            self.selected = {t.name for t in CATALOGUE}
            self.sync_table()

    def action_none(self) -> None:
        if not self.running:
            self.selected.clear()
            self.sync_table()

    def sync_table(self) -> None:
        table = self.query_one(DataTable)
        closure = {t.name for t in CATALOGUE.closure(sorted(self.selected))} if self.selected else set()
        for tool in CATALOGUE:
            mark = MARK_ON if tool.name in self.selected else (
                MARK_DEP if tool.name in closure else MARK_OFF)
            table.update_cell(tool.name, "sel", mark)
            if not self.running:
                state = self.state[tool.name] if tool.name in closure else SKIPPED
                table.update_cell(tool.name, "status", LABEL[state])
        self.refresh_cost()

    # --- running ---------------------------------------------------

    def action_start(self) -> None:
        if self.running or not self.selected:
            return
        self.running = True
        chosen = sorted(self.selected)
        for tool in CATALOGUE.closure(chosen):
            self.state[tool.name] = PENDING
        self.run_worker(self.execute(chosen), thread=True, exclusive=True)

    async def execute(self, chosen: list[str]) -> None:
        log = self.query_one(RichLog)
        snakefile = prepare(CATALOGUE, chosen, self.workdir, self.databases,
                            self.samples, overrides=self.overrides,
                            launcher=self.launcher)
        names = [t.name for t in CATALOGUE.closure(chosen)]

        self.call_from_thread(log.write, f"[bold]Running {len(names)} tools[/]")
        for event in run(snakefile, self.cores, workdir=self.workdir,
                         keep_going=self.keep_going):
            self.call_from_thread(self.apply_event, event)

        # Safe to read self.state here: call_from_thread blocks until the UI
        # thread has applied the update, so every event above has landed.
        self.call_from_thread(self.settle, names)
        done = [n for n in names if self.state[n] == DONE]
        failed = [n for n in names if self.state[n] == FAILED]

        if not done:
            reason = f" Failed: {', '.join(failed)}." if failed else ""
            self.call_from_thread(
                log.write,
                f"[bold red]Nothing ran.[/]{reason} No report written.")
        else:
            if failed:
                self.call_from_thread(
                    log.write,
                    f"[yellow]{len(failed)} of {len(names)} failed:[/] {', '.join(failed)}")
            report = render_report(CATALOGUE, chosen, self.workdir,
                                   self.databases, self.samples)
            self.call_from_thread(log.write, f"[bold green]Report:[/] {report}")
        self.running = False

    def settle(self, names: list[str]) -> None:
        """Once the run is over, a tool that never started reads `not run`.

        Leaving twelve rows at `pending` after the workflow had aborted was how
        a total failure came to look like a run still in progress.

        A row still marked `running` is left alone: that job did start, and the
        stream ended before saying how it went. Calling that `not run` would be
        a false statement rather than an unknown one.
        """
        table = self.query_one(DataTable)
        for name in names:
            if self.state[name] == PENDING:
                self.state[name] = NOT_RUN
                table.update_cell(name, "status", LABEL[NOT_RUN])

    # NB: not `on_event` — Textual reserves that for its own event dispatch,
    # and overriding it swallows every framework message.
    def apply_event(self, event: Event) -> None:
        log = self.query_one(RichLog)
        table = self.query_one(DataTable)

        if event.kind == "job_started" and event.rule:
            self.mark(table, event.rule, RUNNING)
            log.write(f"[cyan]▸[/] {event.rule}")
        elif event.kind == "job_finished" and event.rule:
            self.mark(table, event.rule, DONE)
            log.write(f"[green]✓[/] {event.rule}")
        elif event.kind == "job_error":
            if event.rule:
                self.mark(table, event.rule, FAILED)
            log.write(f"[red]✗ {event.rule or 'error'}[/] {event.message}")
        elif event.kind == "progress" and event.total:
            self.query_one(ProgressBar).update(
                total=event.total, progress=event.done or 0)
        elif event.kind == "error":
            log.write(f"[red]{event.message}[/]")
        elif event.kind == "done":
            log.write("[bold green]Finished[/]")

    def mark(self, table: DataTable, rule: str, state: str) -> None:
        # Rule names replace '-' with '_' for Snakemake; map back.
        name = rule if rule in self.state else rule.replace("_", "-")
        if name in self.state:
            self.state[name] = state
            table.update_cell(name, "status", LABEL[state])


def launch(inputs: list[Path], workdir: Path, databases: Path,
           samples: tuple[str, ...], cores: int,
           selected: list[str] | None = None,
           overrides: dict[str, tuple[tuple[str, str], ...]] | None = None,
           launcher: Sequence[str] | None = None,
           keep_going: bool = False) -> None:
    ComparemTUI(inputs, workdir, databases, samples, cores, selected,
                overrides, launcher, keep_going).run()
