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

import shlex
from pathlib import Path
from time import monotonic

from textual.app import App, ComposeResult
from textual.containers import Horizontal
# textual, not rich: rich is textual's dependency rather than this package's.
from textual.markup import escape
from textual.screen import ModalScreen
from textual.timer import Timer
from textual.widgets import DataTable, Footer, Header, ProgressBar, RichLog, Static

from .catalogue import CATALOGUE
from .cli import any_outputs_exist, lock_files, run_settings, unlock
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

# Liveness. A run's first minutes are silent — the DAG, then conda solving six
# environments — and a still screen in that state is indistinguishable from a
# hung one. The frames are driven by a Textual timer on the UI thread, which is
# the property that makes them worth having: if the interface is genuinely
# blocked, the spinner stops with it rather than reassuring the user it hasn't.
SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
SPINNER_INTERVAL = 0.1  # 10 fps reads as motion; the clock only needs 1 Hz

# What the activity line says when no job is running. Three distinct cases, and
# they are worth distinguishing: the first is where minutes go on a first run,
# and the last is where they go at the end of a long one.
STARTING = "starting up — the DAG, and tool environments on a first run"
WAITING = "no job running — waiting on Snakemake"
REPORTING = "collecting outputs and writing the report"


def elapsed_text(seconds: float) -> str:
    """`3s`, `4m 12s`. Whole seconds, so a short run does not read `0m 3s`."""
    total = int(seconds)
    return f"{total}s" if total < 60 else f"{total // 60}m {total % 60:02d}s"


class ConfirmUnlock(ModalScreen[bool]):
    """Ask before clearing a lock, and say what cannot be known.

    Snakemake refuses to start on a locked directory, and the way out is
    `--unlock` — which means quitting the interface, finding the right
    `--output`, and typing a second command. So the interface offers it.

    It asks first, and this is the one dialog in the application, because the
    question is not "are you sure" but "is something else running": a lock file
    carries a list of paths and no PID, so neither the user nor this code can
    tell a dead run's lock from a live one's. Clearing a live one puts two
    Snakemake processes on the same outputs.
    """

    BINDINGS = [
        ("y", "yes", "Unlock"),
        ("n", "no", "Cancel"),
        ("escape", "no", "Cancel"),
    ]

    def __init__(self, workdir: Path, locks: int) -> None:
        super().__init__()
        self.workdir = workdir
        self.locks = locks

    def compose(self) -> ComposeResult:
        files = "1 lock file" if self.locks == 1 else f"{self.locks} lock files"
        yield Static(
            f"[bold]Unlock the output directory?[/]\n\n"
            f"{escape(str(self.workdir))}\n"
            f"[dim]{files} in .snakemake/locks/[/]\n\n"
            "Snakemake locks a directory while it runs and leaves the lock "
            "behind if it was killed. A lock file lists paths and no process "
            "id, so [bold]this cannot tell you whether that run is still "
            "alive[/] — check that nothing else is writing here before you "
            "say yes. Two Snakemake runs over one directory corrupt each "
            "other's outputs.\n\n"
            "[dim]y unlock · n cancel[/]",
            id="confirm")

    def action_yes(self) -> None:
        self.dismiss(True)

    def action_no(self) -> None:
        self.dismiss(False)


class ConfirmQuit(ModalScreen[str]):
    """Ask before leaving, and say what happens to the jobs.

    Quitting mid-run is silent about the one thing that matters, and what it is
    silent about was measured rather than assumed (Textual 8.2.8, a probe with
    this module's shape): `App.run()` returns in 1.52 s, the worker's next
    `call_from_thread` raises `App is not running` about a second later — so
    Snakemake stops being driven almost at once — and the child process
    **survived its parent** and finished its work twenty seconds on. Nothing
    downstream starts, no report is written, and Snakemake's lock stays on the
    output directory.

    So the dialog says which of the two situations the user is in, because they
    need different things done about them: with a profile the jobs are in a
    queue that does not care that this process is gone, and without one they
    are child processes of it. `s` acts on either — see `cancel.py` for why
    each half has to be done by hand.

    Returns `quit`, `stop` or `stay`.
    """

    BINDINGS = [
        ("y", "quit", "Quit"),
        ("s", "stop", "Quit and stop the jobs"),
        ("n", "stay", "Stay"),
        ("escape", "stay", "Stay"),
    ]

    def __init__(self, in_flight: bool, running: int,
                 profile: str | None) -> None:
        super().__init__()
        # Two different questions: whether a run is going at all, and how many
        # of its jobs have started. A first run spends its first minutes
        # solving conda environments with nothing started, and "0 jobs running"
        # would read as "nothing to lose" at exactly the wrong moment.
        self.in_flight = in_flight
        self.running = running
        self.profile = profile

    def compose(self) -> ComposeResult:
        yield Static(self.text(), id="confirm")

    def text(self) -> str:
        if not self.in_flight:
            # Thin on purpose. There is no consequence to state, and "are you
            # sure" is the dialog this application does not otherwise have.
            return ("[bold]Quit?[/]\n\nNothing is running.\n\n"
                    "[dim]y quit · n stay[/]")
        # `started` under a profile, because a job this interface has not been
        # told about may still be sitting in the queue: what is counted here is
        # what Snakemake said began, which is a floor and not a total.
        state = "started" if self.profile else "running"
        jobs = ("no job has started yet" if not self.running else
                f"1 job {state}" if self.running == 1 else
                f"{self.running} jobs {state}")
        where = (f"[dim]{jobs} · profile {escape(self.profile)}[/]"
                 if self.profile else f"[dim]{jobs} · on this machine[/]")
        what = (
            "Jobs already submitted [bold]keep running[/] — the queue does not "
            "care that this process is gone. Nothing further is submitted."
            if self.profile else
            "Jobs already started are child processes of this one and are "
            "[bold]not killed on the way out[/] — they keep running "
            "unattended. Nothing further starts."
        )
        stop = ("[dim]s[/] quit and cancel this run's queued jobs"
                if self.profile else
                "[dim]s[/] quit and stop them")
        return (
            f"[bold]Quit while a run is going?[/]\n\n{where}\n\n"
            f"{what} No report is written, and the output directory stays "
            f"locked — [dim]u releases it[/].\n\n"
            f"[dim]y[/] quit, leave the run going\n{stop}\n[dim]n[/] stay")

    def action_quit(self) -> None:
        self.dismiss("quit")

    def action_stop(self) -> None:
        # Nothing to stop when nothing is running, and dismissing as `stop`
        # would send the caller looking for processes that were never there.
        self.dismiss("stop" if self.in_flight else "quit")

    def action_stay(self) -> None:
        self.dismiss("stay")


class ComparemTUI(App):
    """Pick tools, watch them run, open the report."""

    CSS = """
    Screen { layout: vertical; }
    #cost { padding: 0 1; color: $text-muted; }
    #where { padding: 0 1; }
    #panes { height: 1fr; }
    #activity { padding: 0 1; height: 1; }
    DataTable { width: 46%; border: round $primary; }
    RichLog { width: 1fr; border: round $primary; padding: 0 1; }
    ProgressBar { padding: 0 1; }
    /* The indeterminate bar is $error by default — red, for a run that is
       merely still going. */
    Bar > .bar--indeterminate { color: $accent; }

    /* The cursor has to be visible on a terminal that has eight colours and no
       more, which is the normal case rather than a corner one: tmux ships
       `default-terminal screen`, an eight-colour TERM, so every run inside tmux
       over SSH is one. Textual's blurred cursor is $primary at 30% alpha, and
       an alpha blend is what defeats the downgrade — #0178D44C over #1E1E1E is
       #153854, which lands on ANSI 8 against a surface on ANSI 0, and
       `block-cursor-blurred-text-style` is `none`, so there is nothing else to
       tell them apart. The command palette is where it showed: its list is
       `can_focus=False` and therefore always drawn blurred, and the selected
       command was invisible. `reverse` is an SGR attribute rather than a
       colour, so it survives every colour system, and every theme — including
       the two whose surface is `ansi_default`, where no colour choice could. */
    DataTable > .datatable--cursor,
    OptionList > .option-list--option-highlighted {
        background: $surface; color: $primary; text-style: reverse;
    }
    /* Which widget has focus is carried by an attribute too, for the same
       reason: a second colour would have the same problem as the first. */
    DataTable:focus > .datatable--cursor,
    OptionList:focus > .option-list--option-highlighted {
        text-style: bold reverse;
    }

    ConfirmUnlock { align: center middle; }
    ConfirmQuit { align: center middle; }
    #confirm {
        width: 64; padding: 1 2;
        border: round $warning; background: $surface;
    }
    """

    BINDINGS = [
        ("space", "toggle", "Select/deselect"),
        ("a", "all", "Select all"),
        ("n", "none", "Select none"),
        ("r", "start", "Run"),
        ("u", "unlock", "Unlock directory"),
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
        # Seeded from `--until` when given, and otherwise empty: the user
        # chooses. Selecting all fourteen by default put gtdbtk's 60.8 GB
        # download one keypress from a user who had not read the table yet, and
        # it made the interface a confirmation step rather than a choice. `a`
        # is still one key away for anyone who does want everything.
        self.selected: set[str] = set(selected) if selected else set()
        # What is already in the output directory, and the live state of this
        # session. Kept apart: the second is overwritten by every event, and the
        # first is the answer to "what did I run here last time" — which the
        # table has to keep giving for a tool the user has since deselected.
        self.disk: dict[str, str] = self.scan()
        self.state: dict[str, str] = dict(self.disk)
        self.running = False
        self.cost_text = ""
        # The animation. `phase` is what the line says when no rule is running,
        # and it has to be maintained rather than assumed: after the last event
        # the worker is still rendering the report, and "starting up" would then
        # be a false statement at the moment a user is most likely reading it.
        self.frame = 0
        self.phase = STARTING
        self.started_at: float | None = None
        self.activity_timer: Timer | None = None
        # Whether Snakemake has said how many jobs there are. Until it has, the
        # bar has nothing to show, and a bar parked at 0% was the other half of
        # this same complaint.
        self.progress_seen = False
        # The name the SLURM executor plugin submits every job of this run
        # under, learned from its own log line. Without it a queue cannot be
        # cancelled — see `cancel.stop_slurm()`.
        self.slurm_run_id: str | None = None
        # What the user chose on the way out, read by `launch()` after the
        # interface has given the terminal back. Recorded at the moment of the
        # decision rather than looked up afterwards: `self.running` is cleared
        # by the worker thread, which is still finishing as the app exits, so
        # reading it after `run()` returns is a race.
        self.left_mid_run = False
        self.stop_requested = False

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
        yield Static(id="activity")
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
            table.add_row(MARK_OFF, tool.name, LABEL[PENDING], tool.summary, key=tool.name)
        # Not refresh_cost(): nothing is selected unless `--until` seeded it,
        # and the marks and statuses have to say so before anything is drawn.
        self.sync_table()
        log = self.query_one(RichLog)
        log.write("[dim]space[/] select · [dim]a[/] all · [dim]r[/] run · [dim]q[/] quit")
        if not self.selected:
            # Said in words, because an empty selection is the one state where
            # the table looks the same whether the interface is waiting for the
            # user or has decided there is nothing to do.
            log.write("[bold]Nothing is selected[/] — pick the analyses you "
                      "want, then press [dim]r[/]")
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
        # Said on opening rather than only on `r`: a lock is the one condition
        # that makes everything else this table promises impossible, and the
        # user is reading the log at that moment anyway.
        if lock_files(self.workdir):
            log.write("[bold yellow]The output directory is locked.[/] A run is "
                      "either still going or was killed before it could clean "
                      "up — [dim]press u to release the lock[/]")

    # --- leaving ----------------------------------------------------

    # Overriding Textual's own `quit` action rather than binding `q` to a new
    # one, so that everything which quits goes through the question: the `q`
    # key, the command palette's Quit entry, and ctrl+c — which in Textual 8
    # does not quit but points at whichever key runs the `quit` action.
    async def action_quit(self) -> None:
        running = sum(1 for tool in CATALOGUE
                      if self.state.get(tool.name) == RUNNING)
        self.push_screen(ConfirmQuit(self.running, running, self.profile),
                         self.departing)

    def departing(self, decision: str | None) -> None:
        if decision == "stay":
            return
        self.left_mid_run = self.running
        self.stop_requested = decision == "stop"
        self.exit()

    # --- the lock ---------------------------------------------------

    def action_unlock(self) -> None:
        log = self.query_one(RichLog)
        if self.running:
            # This run holds the lock. Clearing it would be clearing our own.
            log.write("[dim]not while a run is going — that lock is this run's[/]")
            return
        locks = lock_files(self.workdir)
        if not locks:
            log.write(f"[dim]not locked: {escape(str(self.workdir))}[/]")
            return
        self.push_screen(ConfirmUnlock(self.workdir, len(locks)), self.unlocked)

    def unlocked(self, confirmed: bool | None) -> None:
        """What the dialog decided. Runs on the UI thread; the work does not."""
        if confirmed:
            # A subprocess, so it goes to a thread: it is quick, but a blocked
            # UI thread is the state the spinner exists to make impossible.
            self.run_worker(self.do_unlock, thread=True)
        else:
            self.query_one(RichLog).write("[dim]left locked[/]")

    def do_unlock(self) -> None:
        problem = unlock(self.workdir)
        log = self.query_one(RichLog)
        if problem:
            self.call_from_thread(log.write, f"[red]{escape(problem)}[/]")
            return
        # Re-read rather than assume: `snakemake --unlock` exiting 0 is not the
        # same statement as "there is no lock now", and the difference decides
        # whether pressing `r` is about to fail.
        left = len(lock_files(self.workdir))
        if left:
            self.call_from_thread(
                log.write, f"[yellow]still locked[/] — {left} lock files remain "
                           f"in {escape(str(self.workdir))}/.snakemake/locks/")
        else:
            self.call_from_thread(log.write, "[green]Unlocked.[/] [dim]r runs[/]")

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

    # --- liveness --------------------------------------------------

    def activity_text(self) -> str:
        """One line: a frame, what is running, and how long it has been.

        The elapsed clock does the same work as the spinner and survives a
        screenshot, which is the form the question usually arrives in.
        """
        frame = SPINNER[self.frame % len(SPINNER)]
        active = [t.name for t in CATALOGUE if self.state.get(t.name) == RUNNING]
        shown = ", ".join(active[:3]) if active else self.phase
        if len(active) > 3:
            shown += f" +{len(active) - 3} more"
        since = elapsed_text(monotonic() - self.started_at) if self.started_at else ""
        return f"[cyan]{frame}[/] {shown} [dim]· {since}[/]"

    def refresh_activity(self) -> None:
        self.query_one("#activity", Static).update(self.activity_text())

    def advance_activity(self) -> None:
        self.frame += 1
        self.refresh_activity()

    def stop_activity(self) -> None:
        """Stop the animation, and say why it stopped.

        A still spinner and a hung interface look identical, so the line is
        replaced rather than frozen mid-frame.
        """
        if self.activity_timer is not None:
            self.activity_timer.stop()
            self.activity_timer = None
        took = elapsed_text(monotonic() - self.started_at) if self.started_at else ""
        self.query_one("#activity", Static).update(
            f"[dim]not running — the last run took {took}[/]" if took else "")
        if not self.progress_seen:
            # An indeterminate bar animates for as long as it exists. Left
            # going after the run it would be exactly the false signal the
            # spinner is here to avoid.
            self.query_one(ProgressBar).update(total=100, progress=0)

    # --- running ---------------------------------------------------

    def action_start(self) -> None:
        if self.running:
            return
        if not self.selected:
            # `r` on an empty selection used to be a silent no-op, which is the
            # worst answer available now that empty is the opening state: the
            # key that runs things appears to do nothing at all.
            self.query_one(RichLog).write(
                "[yellow]Nothing selected[/] — [dim]space picks the tool under "
                "the cursor, a selects all[/]")
            return
        if lock_files(self.workdir):
            # Refused here rather than left to Snakemake: from inside the run
            # this arrives as a WorkflowError several seconds in, by which time
            # the interface has already claimed to be starting up.
            self.query_one(RichLog).write(
                "[bold yellow]Not started — the output directory is locked.[/] "
                "[dim]press u to release the lock[/]")
            return
        self.running = True
        chosen = sorted(self.selected)
        # Reset to what is on disk, not to `pending`: Snakemake will skip a tool
        # whose outputs are current, so no event will ever arrive for it, and
        # `pending` would be a claim that its results are not there.
        self.disk = self.scan()
        for tool in CATALOGUE.closure(chosen):
            self.state[tool.name] = self.disk[tool.name]
        self.frame = 0
        self.phase = STARTING
        self.progress_seen = False
        self.started_at = monotonic()
        # Indeterminate until Snakemake says how many jobs there are: a bar
        # that cannot know its total should not draw itself at 0%.
        self.query_one(ProgressBar).update(total=None)
        self.refresh_activity()
        self.activity_timer = self.set_interval(SPINNER_INTERVAL,
                                                self.advance_activity)
        self.run_worker(self.execute(chosen), thread=True, exclusive=True)

    async def execute(self, chosen: list[str]) -> None:
        log = self.query_one(RichLog)
        # try/finally around the whole run: the animation is a claim that work
        # is happening, so it has to come down even on a path that raises.
        try:
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

            # Rendering the report reads every output and takes seconds on a
            # real run, with no events left to arrive. Naming that phase is
            # what keeps the spinner's line true to the end.
            self.phase = REPORTING
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
                # `and not failed`: `have` is satisfied by *any* output in the
                # closure, including one an earlier run left, so a run in which
                # every job failed can reach here with nothing done — and the
                # line below would then contradict the failure line under it.
                if not done and not failed:
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
        finally:
            self.running = False
            self.call_from_thread(self.stop_activity)

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
            # Whatever the wait is from here on, it is not the first solve.
            self.phase = WAITING
        elif event.kind == "job_finished" and event.rule:
            self.mark(table, event.rule, DONE)
            log.write(f"[green]✓[/] {event.rule}")
        elif event.kind == "job_error":
            if event.rule:
                self.mark(table, event.rule, FAILED)
            log.write(f"[red]✗ {event.rule or 'error'}[/] {event.message}")
        elif event.kind == "progress" and event.total:
            self.progress_seen = True
            self.query_one(ProgressBar).update(
                total=event.total, progress=event.done or 0)
        elif event.kind == "slurm_run_id":
            self.slurm_run_id = event.message
            # Shown, not just kept: it is the handle on the queue from any
            # other terminal too — `squeue --name <id>`, `scancel --name <id>`
            # — which matters most in the case this interface cannot help
            # with, a frontend that died with the run still queued.
            log.write(f"[dim]SLURM run id {escape(event.message)}[/]")
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


def departure(mid_run: bool, stop: bool, profile: str | None,
              slurm_run_id: str | None, workdir: Path) -> list[str]:
    """What to do and to say once the interface has closed.

    Deliberately not done inside the app. Cancelling means `scancel`, which the
    plugin's own code allows a minute for, and a process tree that gets three
    seconds to take a SIGTERM — both of them on a screen the user has just
    asked to leave. Doing it here means the terminal is already back, the
    outcome is ordinary text they can scroll to, and none of it can leave a
    half-torn-down Textual display behind.

    Split from `launch()` so it can be tested: it is the branch that decides
    whether a queue gets cancelled, and it should not need a terminal to check.
    """
    if not mid_run:
        return []
    if not stop:
        note = ["Left the run going."]
        if profile:
            note.append("  Jobs already submitted stay in the queue. Nothing "
                        "further will be submitted and no report was written.")
            if slurm_run_id:
                note.append(f"  squeue --name {slurm_run_id}   "
                            f"scancel --name {slurm_run_id}")
        else:
            note.append("  Jobs already started keep running unattended. "
                        "Nothing further will start and no report was written.")
        # Quoted, because this line is meant to be pasted: an output directory
        # with a space in it otherwise prints a command that reads the tail as
        # a positional. `_invocation()` and `snakefile._rule` already quote for
        # the same reason, and `tests/E._faecium/116_2 duplicate.fna` exists to
        # keep that honest.
        note.append(f"  The output directory is still locked: "
                    f"comparem2 --unlock --output {shlex.quote(str(workdir))}")
        return note

    from .cancel import stop_local, stop_slurm

    lines = []
    if profile:
        if slurm_run_id:
            lines.append(stop_slurm(slurm_run_id))
        else:
            # No id means the plugin never announced one: nothing had been
            # submitted yet, or the executor is not the SLURM plugin at all —
            # cluster-generic for PBS and SGE goes through this same branch and
            # has no equivalent handle. Either way, guessing at job ids would
            # be worse than saying so.
            lines.append("No SLURM run id had been announced, so the queue was "
                         "left alone — either nothing was submitted yet, or "
                         "this profile's executor is not the SLURM plugin. "
                         "Check with: squeue -u $USER")
    # Run in both cases. Under a profile the analyses are in the queue, but the
    # four database downloads are `localrules` and run here on the login node —
    # a 60.8 GB GTDB fetch is a child of this process, not a job.
    lines.append(stop_local())
    lines.append(f"The output directory is still locked: "
                 f"comparem2 --unlock --output {shlex.quote(str(workdir))}")
    return lines


def launch(inputs: list[Path], workdir: Path, databases: Path,
           samples: tuple[str, ...], cores: int | None,
           selected: list[str] | None = None,
           overrides: dict[str, tuple[tuple[str, str], ...]] | None = None,
           keep_going: bool = False,
           conda_prefix: Path | None = None, command: str | None = None,
           profile: str | None = None) -> None:
    app = ComparemTUI(inputs, workdir, databases, samples, cores, selected,
                      overrides, keep_going, conda_prefix, command, profile)
    app.run()
    for line in departure(app.left_mid_run, app.stop_requested, profile,
                          app.slurm_run_id, workdir):
        print(line)
