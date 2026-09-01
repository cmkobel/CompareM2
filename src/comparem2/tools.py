"""The tool contract.

Every analysis in CompareM2 v3 is one `Tool` value. Adding or removing an
analysis means adding or removing a spec — the runner, CLI, TUI and report all
read from this contract and none of them need to know which tools exist.

A tool declares its outputs rather than writing wherever it likes. That single
choice buys resumability (a job whose outputs are present and current is
skipped), progress reporting, and the report's knowledge of what to render.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
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
    """

    name: str
    url: str
    size: int  # bytes, as served
    sha256: str | None = None

    @property
    def human_size(self) -> str:
        n = float(self.size)
        for unit in ("B", "kB", "MB", "GB", "TB"):
            if n < 1000 or unit == "TB":
                return f"{n:.1f} {unit}" if unit != "B" else f"{n:.0f} B"
            n /= 1000
        raise AssertionError("unreachable")


@dataclass(frozen=True)
class Context:
    """Everything a tool needs to build its command line."""

    workdir: Path  # output directory for this run
    databases: Path  # root of the database directory
    threads: int
    assemblies: tuple[Path, ...]  # every input assembly
    sample: str | None = None  # set for Scope.GENOME, None for Scope.SET

    @property
    def assembly(self) -> Path:
        """The single input assembly. Only valid for Scope.GENOME."""
        if self.sample is None:
            raise ValueError("assembly is only defined for Scope.GENOME tools")
        return next(a for a in self.assemblies if a.stem == self.sample)

    def out(self, *parts: str) -> Path:
        """A path inside this tool's output directory."""
        if self.sample is None:
            return self.workdir.joinpath(*parts)
        return self.workdir.joinpath("samples", self.sample, *parts)


@dataclass(frozen=True)
class Tool:
    """One analysis step."""

    name: str
    summary: str  # one line; shown in the TUI and above the report section
    scope: Scope
    conda: tuple[str, ...]  # conda package specs, e.g. ("bioconda::checkm2=1.1.0",)
    command: Callable[[Context], Sequence[str]]
    outputs: Callable[[Context], Sequence[Path]]
    needs: tuple[str, ...] = ()  # names of tools that must finish first
    database: Database | None = None
    threads: int = 1

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
        """Total bytes to download before `selected` can run."""
        return sum(db.size for db in self.databases(selected))

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
