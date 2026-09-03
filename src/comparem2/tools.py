"""The tool contract.

Every analysis in CompareM2 v3 is one `Tool` value. Adding or removing an
analysis means adding or removing a spec — the runner, CLI, TUI and report all
read from this contract and none of them need to know which tools exist.

A tool declares its outputs rather than writing wherever it likes. That single
choice buys resumability (a job whose outputs are present and current is
skipped), progress reporting, and the report's knowledge of what to render.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class Scope(Enum):
    """How often a tool runs."""

    GENOME = "genome"  # once per input assembly
    SET = "set"  # once over the whole set of assemblies


@dataclass(frozen=True)
class Database:
    """A reference database a tool needs before it can run.

    `size` is measured, never estimated: it is the `content-length` of `url`.
    The total across selected tools is shown to the user before anything is
    downloaded, because install weight is the constraint v3 exists to respect.

    `fetch` is what makes that promise real. Until 2026-09-02 this class carried
    a `url` that no code read: the total was announced and then nothing
    downloaded anything, so a run failed minutes later inside a tool with a
    tool-specific error. Declaring a size without a way to satisfy it is worse
    than declaring neither.
    """

    name: str
    url: str  # or the command that fetches it, when there is no static URL
    size: int | None = None  # bytes; None means not yet measured — never guess
    sha256: str | None = None
    # Whichever digest the source publishes, recorded rather than computed. The
    # fetch does not check it: hashing 141.4 GB costs minutes, and this exists
    # to tell a truncated download from a corrupt one after something failed.
    md5: str | None = None
    # How to fetch it, given the database root. Returns steps, each an argv
    # list — same discipline as `Tool.command`, and for the same reason.
    fetch: Callable[[Path], Sequence[Sequence[str]]] | None = None
    # What the fetch itself needs installed. A download step is a rule like any
    # other, so it gets an environment like any other — which it must, because
    # two of these fetches run a tool binary rather than curl (`bakta_db
    # download`, `amrfinder -u`) and the other two need curl and tar. Omitting
    # it means the download rule runs in the pipeline's own environment, where
    # none of those exist.
    conda: tuple[str, ...] = ()
    # Which named environment those packages constitute. Same contract as
    # `Tool.environment`: the name selects the file, the packages are its
    # content, and a name used twice must mean the same packages both times.
    environment: str = "main"
    # Path relative to the database root whose presence means "ready". This is
    # the download rule's declared output, so Snakemake skips a database that
    # is already complete and re-runs one that is half-finished. Prefer a real
    # file the tool itself needs over a stamp: a stamp can outlive the data.
    ready: str = ""
    # Set when the fetch cannot write under the database root because the tool
    # forbids it. amrfinder is the only one — `amrfinder -u` rejects `-d` with
    # "only operates on the default database directory", verified 2026-09-02 —
    # so its data lands in the conda prefix and `ready` is only a stamp. Which
    # means it does not survive the environment being rebuilt.
    out_of_tree: bool = False

    def marker_root(self, databases: Path, workdir: Path) -> Path:
        """Where this database's readiness marker lives.

        Under `--databases` for everything that stores its data there. For an
        out-of-tree database it is the **run's output directory**, and that is
        not cosmetic: AMRFinder's data lives in the tool's conda prefix, so a
        marker under `--databases` outlives the thing it describes. Rebuild the
        environment and the marker still says ready while the data is gone.

        Verified the hard way 2026-09-02: a thirteen-tool run failed with
        `No valid AMRFinder database is found` against a four-hour-old stamp,
        because the environment had been rebuilt in between. Per-run readiness
        means at most one run can be wrong about it, and a fresh output
        directory is always right.
        """
        return workdir if self.out_of_tree else databases

    def ready_path(self, databases: Path, workdir: Path) -> Path:
        """The file that proves this database is usable."""
        root = self.marker_root(databases, workdir)
        return root / (self.ready or f"{self.name}/.fetched")

    @property
    def human_size(self) -> str:
        if self.size is None:
            return "unmeasured"
        n = float(self.size)
        for unit in ("B", "kB", "MB", "GB", "TB"):
            if n < 1000 or unit == "TB":
                return f"{n:.1f} {unit}" if unit != "B" else f"{n:.0f} B"
            n /= 1000
        raise AssertionError("unreachable")


@dataclass(frozen=True)
class Context:
    """Everything a tool needs to build its command line.

    Inputs are canonicalised: before anything runs, each input assembly is
    linked to `<workdir>/samples/<sample>/<sample>.fna`, and every tool reads
    from there. That is what lets the same code render both a real path and a
    Snakemake wildcard path — build a Context with `sample="{sample}"` and
    every path comes out templated.
    """

    workdir: Path  # output directory for this run
    databases: Path  # root of the database directory
    # `int` when running for real. During rule generation this is the literal
    # string "{threads}", so Snakemake substitutes the count it actually
    # granted — the same substitution trick as `sample="{sample}"`.
    threads: int | str
    samples: tuple[str, ...]  # every sample name, in input order
    sample: str | None = None  # set for Scope.GENOME, None for Scope.SET
    # Tool defaults merged with any user overrides. Carried over from v2's
    # `set_<tool>--<flag>: <value>` passthrough, which let users reach any
    # tool argument without the pipeline having to know about it.
    params: tuple[tuple[str, str], ...] = ()

    def args(self) -> list[str]:
        """Parameters as a flat argument list, in declaration order.

        A flag with an empty value is emitted bare, so `--verbose: ""` works.
        """
        out: list[str] = []
        for flag, value in self.params:
            out.append(flag)
            if value != "":
                out.append(value)
        return out

    @property
    def assembly(self) -> Path:
        """This sample's canonical assembly. Only valid for Scope.GENOME."""
        if self.sample is None:
            raise ValueError("assembly is only defined for Scope.GENOME tools")
        return self.sample_out(self.sample, f"{self.sample}.fna")

    @property
    def assemblies(self) -> list[Path]:
        """Every sample's canonical assembly, in input order."""
        return [self.sample_out(s, f"{s}.fna") for s in self.samples]

    def out(self, *parts: str) -> Path:
        """A path inside this tool's output directory."""
        if self.sample is None:
            return self.workdir.joinpath(*parts)
        return self.sample_out(self.sample, *parts)

    def sample_out(self, sample: str, *parts: str) -> Path:
        """A path inside another sample's output directory.

        Set-scope tools need this to collect per-genome results from the tools
        they depend on — the pangenome reading every genome's annotation, say.
        """
        return self.workdir.joinpath("samples", sample, *parts)


@dataclass(frozen=True)
class Tool:
    """One analysis step."""

    name: str
    summary: str  # one line; shown in the TUI and above the report section
    scope: Scope
    # The full package list of the environment this tool runs in — not just its
    # own package. Several tools share one list, which is the point: Snakemake
    # deploys by content, so rules whose env files are identical are one
    # environment on disk.
    conda: tuple[str, ...]
    command: Callable[[Context], Sequence[str]]
    outputs: Callable[[Context], Sequence[Path]]
    needs: tuple[str, ...] = ()  # names of tools that must finish first
    database: Database | None = None
    threads: int = 1
    # Default arguments, overridable per run. Values carried from v2's
    # config.yaml so v3 reproduces v2's behaviour unless told otherwise.
    params: tuple[tuple[str, str], ...] = ()
    # Which named environment `conda` constitutes. This is what the rule's
    # `conda:` directive points at, so it is also the count the user pays for:
    # two names means two solves and two copies on disk. A name used twice must
    # carry the same packages both times, and a test enforces that — otherwise
    # one file would be written twice with different content and whichever rule
    # rendered last would silently decide what the other one ran in.
    environment: str = "main"
    # Some tools write their result to stdout rather than taking an -o flag.
    # Commands stay as argument lists — never hand-built shell strings — so the
    # redirect is declared here and added by whatever runs the command.
    stdout_to_output: bool = False
    # Input files this tool needs that are not another tool's output: a mapping
    # of path to content, rendered from the Context. `prepare()` writes them and
    # the generated rule declares them as inputs, so they are part of the DAG
    # rather than something a shell step has to remember to create.
    #
    # GTDB-Tk is the only user: once inputs are canonicalised each genome sits
    # in its own directory, so `--genome_dir` cannot be used and it takes a
    # two-column `--batchfile` instead. Writing that from the Context is exact —
    # the sample list is already known before anything runs.
    files: Callable[[Context], Mapping[Path, str]] | None = None
    # Steps to run after the command, each an argument list, same discipline as
    # `Database.fetch`. For turning what a tool actually writes into what the
    # spec declared — GTDB-Tk emits `bac120` and `ar53` summaries separately and
    # the pipeline declares one table.
    #
    # This is not a licence for arbitrary shell work: a tool needing several of
    # these is a tool whose spec is lying about what it does.
    post: Callable[[Context], Sequence[Sequence[str]]] | None = None
    # Environment variables the tool needs, given its Context. Some tools take
    # their database location only this way: GTDB-Tk reads GTDBTK_DATA_PATH and
    # has no equivalent flag, so without this its `--databases` value would be
    # silently ignored. Also the route to FastTreeMP, which takes its thread
    # count from OMP_NUM_THREADS rather than an argument.
    env: Callable[[Context], Sequence[tuple[str, str]]] | None = None

    def __post_init__(self) -> None:
        if not self.summary:
            raise ValueError(f"{self.name}: summary is required — it is what the report shows")



class Registry:
    """The set of tools available to a run."""

    def __init__(self, tools: Sequence[Tool] = ()) -> None:
        self._tools: dict[str, Tool] = {}
        for tool in tools:
            self.add(tool)

    def add(self, tool: Tool) -> Tool:
        if tool.name in self._tools:
            raise ValueError(f"duplicate tool name: {tool.name}")
        self._tools[tool.name] = tool
        return tool

    def __contains__(self, name: object) -> bool:
        return name in self._tools

    def __getitem__(self, name: str) -> Tool:
        return self._tools[name]

    def __iter__(self):
        return iter(self._tools.values())

    def __len__(self) -> int:
        return len(self._tools)

    def databases(self, selected: Sequence[str] | None = None) -> list[Database]:
        """Distinct databases needed by `selected` and its dependencies."""
        seen: dict[str, Database] = {}
        for tool in self.closure(selected):
            if tool.database is not None:
                seen.setdefault(tool.database.name, tool.database)
        return list(seen.values())

    def install_size(self, selected: Sequence[str] | None = None) -> int:
        """Total known bytes to download before `selected` can run.

        Databases whose size has not been measured contribute nothing here, so
        this is a lower bound. Call `unmeasured()` alongside it and say so —
        a total that silently omits an unknown is worse than no total.
        """
        return sum(db.size for db in self.databases(selected) if db.size is not None)

    def unmeasured(self, selected: Sequence[str] | None = None) -> list[Database]:
        """Databases in `selected` whose download size is not yet known."""
        return [db for db in self.databases(selected) if db.size is None]

    def closure(self, selected: Sequence[str] | None = None) -> list[Tool]:
        """`selected` plus everything it depends on, in declaration order.

        Ordering the work is Snakemake's job, not ours. This exists so the CLI
        and TUI can answer "what will actually run, and how much do I have to
        download first" before anything is downloaded or submitted.
        """
        wanted = list(self._tools) if selected is None else list(selected)
        needed: set[str] = set()
        stack = list(wanted)
        while stack:
            name = stack.pop()
            if name in needed:
                continue
            if name not in self._tools:
                raise KeyError(f"unknown tool: {name}")
            needed.add(name)
            stack.extend(self._tools[name].needs)
        return [t for t in self if t.name in needed]
