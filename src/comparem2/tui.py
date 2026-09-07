"""Keyboard-driven interface.

The point is that a non-bioinformatician should be able to see what will run,
what it will cost to install, and how far along it is, without reading a manual
or a scrolling wall of Snakemake output.

Tool state comes from `catalogue.py` and progress from `runner.Event`, so this
module knows nothing about Snakemake and nothing about which tools exist.

Opening state comes from the output directory. A tool's declared outputs are
the same thing Snakemake's resumability is decided on, so what the table says
on startup is what a re-run would actually skip.
"""

from __future__ import annotations

from pathlib import Path

from textual.app import App, ComposeResult
from textual.containers import Horizontal
# textual, not rich: rich is textual's dependency rather than this package's.
from textual.markup import escape
from textual.widgets import DataTable, Footer, Header, ProgressBar, RichLog, Static

from .catalogue import CATALOGUE
from .cli import any_outputs_exist, run_settings
from .report import render_report
from .runner import Event, run
from .snakefile import prepare
from .tools import completion

PENDING, RUNNING, DONE, FAILED, SKIPPED, NOT_RUN = "·", "▸", "✓", "✗", "–", "○"
# Results found on disk when the interface opened, as against results this
# session produced. Distinct from DONE on purpose: "done" is something the user
# watched happen, and claiming it for a file left by a run last week would put
# this session's name on someone else's output.
EXISTS, PARTIAL = "◆", "◐"
LABEL = {PENDING: "pending", RUNNING: "running", DONE: "done",
         FAILED: "failed", SKIPPED: "not selected", NOT_RUN: "not run",
         EXISTS: "already run", PARTIAL: "part-finished"}

# The selection column. These were `[x]`, `[+]` and `[ ]` and rendered as
# nothing at all: a DataTable cell given a `str` is parsed as Rich markup, and
# `[x]` is a tag, not text. The selection UI had no visible selection.
MARK_ON, MARK_DEP, MARK_OFF = "▣", "▨", "▢"


class ComparemTUI(App):
    """Pick tools, watch them run, open the report."""

    CSS = """
    Screen { layout: vertical; }
    #cost { padding: 0 1; color: $text-muted; }
    #where { padding: 0 1; }
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
                 samples: tuple[str, ...], cores: int | None,
                 selected: list[str] | None = None,
                 overrides: dict[str, tuple[tuple[str, str], ...]] | None = None,
                 keep_going: bool = False,
                 conda_prefix: Path | None = None,
                 command: str | None = None,
                 profile: str | None = None) -> None:
        super().__init__()
        self.inputs = inputs
        self.workdir = workdir
        self.databases = databases
        self.samples = samples
        self.cores = cores
        self.overrides = overrides
        self.keep_going = keep_going
        self.conda_prefix = conda_prefix
        # A Snakemake profile directory or name. When set, jobs go to a queue
        # instead of this machine, and the run this interface is watching is a
        # frontend process waiting on sbatch — which is the case the progress
        # display is most useful in.
        self.profile = profile
        # Passed in rather than read from sys.argv here: the CLI already
        # renders it, and the TUI's job is to display what it was given.
        self.command = command
        # Seeded from `--until` when given. Selecting everything by default
        # puts gtdbtk's 60.8 GB one keypress away, so a user who named the
        # tools they want on the command line gets exactly those.
        self.selected: set[str] = set(selected) if selected else {t.name for t in CATALOGUE}
        # What is already in the output directory, and the live state of this
        # session. Kept apart: the second is overwritten by every event, and the
        # first is the answer to "what did I run here last time" — which the
        # table has to keep giving for a tool the user has since deselected.
        self.disk: dict[str, str] = self.scan()
        self.state: dict[str, str] = dict(self.disk)
        self.running = False
        self.cost_text = ""

    def scan(self) -> dict[str, str]:
        """Read the output directory: which tools already have results there.

        Opening the interface on a directory that had been run in showed
        fourteen rows of `pending`, so the only way to find out what was already
        done was to run it again and watch Snakemake skip things.
        """
        found = {}
        for tool in CATALOGUE:
            done = completion(tool, self.workdir, self.databases, self.samples)
            found[tool.name] = EXISTS if done.done else (
                PARTIAL if done.partial else PENDING)
        return found

    # --- layout ----------------------------------------------------

    def compose(self) -> ComposeResult:
        yield Header(show_clock=False)
        yield Static(self.where_text(), id="where")
        yield Static(id="cost")
        with Horizontal(id="panes"):
            yield DataTable(cursor_type="row", zebra_stripes=True)
            yield RichLog(highlight=False, markup=True, wrap=True)
        yield ProgressBar(total=100, show_eta=False)
        yield Footer()

    def where_text(self) -> str:
        """The four locations this run depends on, and where each came from.

        Databases and tool environments are settable by environment variable,
        which makes them the two settings most likely to be wrong without
        anyone noticing: both are exported once in a shell profile and never
        looked at again, and either one pointing somewhere unexpected costs a
        re-download or a re-solve rather than an error.
        """
        rows = run_settings(self.workdir, self.databases, self.conda_prefix,
                            self.profile)
        width = max(len(what) for what, _, _ in rows)
        # Escaped: a path may contain '[', which Rich reads as a markup tag.
        return "\n".join(
            f"[dim]{what.ljust(width)}[/] {escape(where)}  [dim]{escape(origin)}[/]"
            for what, where, origin in rows)

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
        log = self.query_one(RichLog)
        log.write("[dim]space[/] select · [dim]r[/] run · [dim]q[/] quit")
        # Said once, in words, because the table's own answer is spread over
        # fourteen rows: a directory that has been run in before is the case
        # where "what still needs doing" is the first question.
        done = [n for n, s in self.disk.items() if s == EXISTS]
        part = [n for n, s in self.disk.items() if s == PARTIAL]
        if done:
            log.write(f"[bold]{len(done)} of {len(CATALOGUE)} tools[/] already ran "
                      "in this directory — [dim]a run re-uses their output[/]")
        if part:
            # Missing one declared output is what makes Snakemake re-run a rule,
            # so this is a statement about what pressing `r` will do.
            log.write(f"[yellow]part-finished, and will be redone:[/] "
                      f"{', '.join(sorted(part))}")

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
                # An unselected tool still reports results it has on disk.
                # Whether it is selected is what the mark column is for, and
                # "not selected" over a finished analysis reads as "missing".
                state = (self.state[tool.name] if tool.name in closure
                         else (self.disk[tool.name] if self.disk[tool.name] != PENDING
                               else SKIPPED))
                table.update_cell(tool.name, "status", LABEL[state])
        self.refresh_cost()

    # --- running ---------------------------------------------------

    def action_start(self) -> None:
        if self.running or not self.selected:
            return
        self.running = True
        chosen = sorted(self.selected)
        # Reset to what is on disk, not to `pending`: Snakemake will skip a tool
        # whose outputs are current, so no event will ever arrive for it, and
        # `pending` would be a claim that its results are not there.
        self.disk = self.scan()
        for tool in CATALOGUE.closure(chosen):
            self.state[tool.name] = self.disk[tool.name]
        self.run_worker(self.execute(chosen), thread=True, exclusive=True)

    async def execute(self, chosen: list[str]) -> None:
        log = self.query_one(RichLog)
        snakefile = prepare(CATALOGUE, chosen, self.workdir, self.databases,
                            self.samples, overrides=self.overrides)
        names = [t.name for t in CATALOGUE.closure(chosen)]

        self.call_from_thread(log.write, f"[bold]Running {len(names)} tools[/]")
        # Snakemake's own "Creating conda environment" lines are quietened in
        # the API path, so a first run would otherwise look hung while the
        # environments solve.
        self.call_from_thread(
            log.write, f"[dim]deploying tool environments in {self.conda_prefix}"
                       " — first run only[/]")
        for event in run(snakefile, self.cores, workdir=self.workdir,
                         keep_going=self.keep_going,
                         conda_prefix=self.conda_prefix,
                         profile=self.profile):
            self.call_from_thread(self.apply_event, event)

        # Safe to read self.state here: call_from_thread blocks until the UI
        # thread has applied the update, so every event above has landed.
        self.call_from_thread(self.settle, names)
        done = [n for n in names if self.state[n] == DONE]
        failed = [n for n in names if self.state[n] == FAILED]
        # Is there anything to report? Asked of the files, and of the same
        # function the CLI asks, so the two paths cannot disagree about it.
        # Snakemake emits no job events for a rule it skips, so a directory that
        # was already current produced none at all, every row settled to
        # `not run`, and the report was withheld over a complete set of outputs.
        have = any_outputs_exist(chosen, self.workdir, self.databases, self.samples)

        if not have:
            reason = f" Failed: {', '.join(failed)}." if failed else ""
            self.call_from_thread(
                log.write,
                f"[bold red]Nothing ran.[/]{reason} No report written.")
        else:
            if not done:
                self.call_from_thread(
                    log.write,
                    "[dim]Nothing to do — every selected tool's output was "
                    "already up to date.[/]")
            if failed:
                self.call_from_thread(
                    log.write,
                    f"[yellow]{len(failed)} of {len(names)} failed:[/] {', '.join(failed)}")
            report = render_report(CATALOGUE, chosen, self.workdir,
                                   self.databases, self.samples,
                                   command=self.command)
            self.call_from_thread(log.write, f"[bold green]Report:[/] {report}")
        self.running = False

    def settle(self, names: list[str]) -> None:
        """Once the run is over, a tool that never started reads `not run`.

        Leaving twelve rows at `pending` after the workflow had aborted was how
        a total failure came to look like a run still in progress.

        A row still marked `running` is left alone: that job did start, and the
        stream ended before saying how it went. Calling that `not run` would be
        a false statement rather than an unknown one.

        The disk is re-read first, so a tool that finished without its finish
        event arriving — and one whose outputs a failed run part-wrote — reads
        as what it left behind rather than as `not run`.
        """
        table = self.query_one(DataTable)
        self.disk = self.scan()
        for name in names:
            if self.state[name] == PENDING:
                self.state[name] = (self.disk[name] if self.disk[name] != PENDING
                                    else NOT_RUN)
                table.update_cell(name, "status", LABEL[self.state[name]])

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
           samples: tuple[str, ...], cores: int | None,
           selected: list[str] | None = None,
           overrides: dict[str, tuple[tuple[str, str], ...]] | None = None,
           keep_going: bool = False,
           conda_prefix: Path | None = None, command: str | None = None,
           profile: str | None = None) -> None:
    ComparemTUI(inputs, workdir, databases, samples, cores, selected,
                overrides, keep_going, conda_prefix, command, profile).run()
