"""Command-line entry point.

Canonicalises inputs, generates the workflow, hands off to Snakemake, and
renders the report. `--tui` hands the same arguments to the TUI, which shares
`snakefile.prepare()` with this path rather than building the workflow itself.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from .catalogue import CATALOGUE
from .report import render_report
from .snakefile import prepare


def default_databases() -> Path:
    """Where databases go when `-d` is not given.

    Home-relative, so it is the same location for every run, and outside any
    checkout so deleting a checkout does not cost a re-download. The default
    used to be `./databases`, which meant the same command run from two
    directories silently fetched a second copy of up to 143 GB.

    `$COMPAREM2_DATABASES` exists because home is the wrong place for 143 GB on
    a cluster with a home quota, and it is the only way to move the location
    once for every run — `-d` has to be retyped on each invocation.
    """
    override = os.environ.get("COMPAREM2_DATABASES")
    if override:
        return Path(override).expanduser()
    return Path.home() / ".comparem2" / "databases"


def slug(stem: str) -> str:
    """A sample name safe to use as a Snakemake wildcard and a path component.

    Spaces and shell metacharacters in filenames are normal in real data — v2's
    own test set contains `116_2 duplicate.fna` — and a wildcard containing a
    space silently produces a broken rule rather than an error.
    """
    cleaned = "".join(c if (c.isalnum() or c in "-_.") else "_" for c in stem)
    return cleaned.strip("_.") or "sample"


def parse_overrides(settings: list[str]) -> dict[str, tuple[tuple[str, str], ...]]:
    """Turn `--set tool--flag=value` into per-tool parameter tuples.

    Overrides are merged onto the tool's defaults: naming one flag replaces
    only that flag, which is what v2's config-file passthrough did. A flag with
    an empty value is passed bare.

    The flag keeps whatever dashes it was written with, so both `skani-c=125`
    and `treecluster--threshold=0.1` work. An earlier version split on the
    first `--` and prepended `--` to whatever followed, which could not express
    a single-dash flag at all: `skani-c=125` was rejected outright, and
    `skani--c=125` silently produced `-c 70 --c 125` — failing to replace the
    default *and* inventing a flag skani does not have.

    Tool names may themselves contain a dash (`snp-dists`), so the tool is
    matched against the catalogue longest-first rather than found by splitting.
    """
    known = sorted((t.name for t in CATALOGUE), key=len, reverse=True)
    merged: dict[str, dict[str, str]] = {}
    for setting in settings:
        tool = next((n for n in known if setting.startswith(n + "-")), None)
        if tool is None:
            raise SystemExit(
                f"--set needs TOOL-FLAG=VALUE with a known tool, got {setting!r}"
            )
        flag, sep, value = setting[len(tool):].partition("=")
        merged.setdefault(tool, dict(CATALOGUE[tool].params))[flag] = value if sep else ""

    return {
        tool: tuple(flags.items())
        for tool, flags in merged.items()
    }


def canonicalise(inputs: list[Path], workdir: Path) -> tuple[str, ...]:
    """Link each input to <workdir>/samples/<sample>/<sample>.fna.

    Everything downstream reads from this layout, so tools never see the user's
    directory structure and sample names come from exactly one place.
    """
    samples: list[str] = []
    for path in inputs:
        sample = slug(path.stem)
        if sample != path.stem:
            print(f"note: {path.name!r} -> sample {sample!r}", file=sys.stderr)
        if sample in samples:
            raise SystemExit(
                f"duplicate sample name {sample!r} from {path} — "
                "rename the input, sample names must be unique"
            )
        target = workdir / "samples" / sample / f"{sample}.fna"
        target.parent.mkdir(parents=True, exist_ok=True)
        # A stale link is not the same as a missing one: `exists()` follows the
        # link and reports False when it dangles, so a plain `symlink_to` would
        # then raise FileExistsError. This happens whenever a results directory
        # outlives a move of the inputs.
        if target.is_symlink() and target.resolve() != path.resolve():
            target.unlink()
        if not target.exists():
            if target.is_symlink():
                target.unlink()
            target.symlink_to(path.resolve())
        samples.append(sample)
    return tuple(samples)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="comparem2", description="A wide view of a set of assemblies.")
    p.add_argument("inputs", nargs="+", type=Path, help="assembly FASTA files")
    p.add_argument("-o", "--output", type=Path, default=Path("results_comparem2"))
    p.add_argument("-d", "--databases", type=Path, default=None,
                   help="where databases live; shared across runs (default: "
                        f"{default_databases()}, overridden for every run by "
                        "$COMPAREM2_DATABASES)")
    p.add_argument("-t", "--cores", type=int, default=4)
    p.add_argument("--until", nargs="*", default=None, metavar="TOOL",
                   help="run only these tools and their dependencies")
    p.add_argument("--set", action="append", default=[], metavar="TOOL--FLAG=VALUE",
                   help="override a tool argument, e.g. --set treecluster--threshold=0.1 "
                        "(v2 spelled this set_treecluster--threshold in config.yaml)")
    p.add_argument("--tui", action="store_true", help="interactive keyboard interface")
    p.add_argument("--isolated-launcher", default=None, metavar="CMD",
                   help="how to enter an isolated tool's environment; {tool} is "
                        "substituted, e.g. 'pixi run -e {tool}'")
    p.add_argument("--keep-going", action="store_true",
                   help="keep running independent tools after one fails")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--report-only", action="store_true",
                   help="re-render the report from existing outputs")
    args = p.parse_args(argv)

    missing = [str(i) for i in args.inputs if not i.exists()]
    if missing:
        raise SystemExit(f"no such file: {', '.join(missing)}")

    # Absolute, because Snakemake is given this as its working directory (see
    # below) and every generated path has to survive that.
    workdir: Path = args.output.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    databases: Path = (args.databases or default_databases()).expanduser().resolve()
    samples = canonicalise(args.inputs, workdir)

    tools = CATALOGUE.closure(args.until)
    print(f"{len(samples)} assemblies, {len(tools)} tools", file=sys.stderr)

    # Only report what is actually going to be fetched. Announcing a total that
    # includes databases already on disk is how "databases: 143.2 GB" came to
    # be printed before a run that downloaded nothing at all.
    pending = [db for db in CATALOGUE.databases(args.until)
               if not db.ready_path(databases).exists()]
    if pending:
        known = sum(db.size for db in pending if db.size is not None)
        unmeasured = [db for db in pending if db.size is None]
        size = f"{known / 1e9:.1f} GB" if known else ""
        if unmeasured:
            size += (" + " if size else "") + f"{len(unmeasured)} of unknown size"
        names = ", ".join(db.name for db in pending)
        # Say where, because the default is no longer in sight of the cwd and
        # this is the line that precedes a download of up to 143 GB.
        print(f"to download: {names} ({size}) -> {databases}", file=sys.stderr)

    overrides = parse_overrides(args.set)
    launcher = args.isolated_launcher.split() if args.isolated_launcher else None

    if args.tui:
        if args.dry_run:
            raise SystemExit(
                "--dry-run and --tui do not combine: the tool list is already "
                "the dry run, and it shows the download size too")
        from .tui import launch

        launch(args.inputs, workdir, databases, samples, args.cores,
               selected=args.until, overrides=overrides, launcher=launcher,
               keep_going=args.keep_going)
        return 0

    snakefile = prepare(CATALOGUE, args.until, workdir, databases, samples,
                        overrides=overrides, launcher=launcher)

    if not args.report_only:
        cmd = [
            "snakemake", "--snakefile", str(snakefile), "--cores", str(args.cores),
            "--rerun-incomplete",
            # Snakemake locks its working directory, not its output paths, so
            # without this every run in one checkout shares `./.snakemake` —
            # two runs collided even with different --output, and a killed run
            # left a lock the next could not clear. Pointing it at the output
            # directory makes the lock mean what it should: runs writing to
            # different places are independent, runs writing to the same place
            # correctly refuse to overlap. This is why workdir is resolved to
            # an absolute path above.
            "--directory", str(workdir),
        ]
        if args.keep_going:
            cmd.append("--keep-going")
        if args.dry_run:
            cmd.append("--dry-run")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            return result.returncode
        if args.dry_run:
            return 0

    report = render_report(CATALOGUE, args.until, workdir, databases, samples)
    print(f"report: {report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
