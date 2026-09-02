"""Command-line entry point.

Canonicalises inputs, generates the workflow, hands off to Snakemake, and
renders the report. The TUI will drive the same functions.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from .catalogue import CATALOGUE
from .report import render_report
from .snakefile import render, render_envs


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
    """
    merged: dict[str, dict[str, str]] = {}
    for setting in settings:
        if "--" not in setting:
            raise SystemExit(f"--set needs TOOL--FLAG=VALUE, got {setting!r}")
        tool, _, rest = setting.partition("--")
        flag, sep, value = rest.partition("=")
        if tool not in CATALOGUE:
            raise SystemExit(f"--set names unknown tool {tool!r}")
        merged.setdefault(tool, dict(CATALOGUE[tool].params))["--" + flag] = value if sep else ""

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
    p.add_argument("-d", "--databases", type=Path, default=Path("databases"))
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

    workdir: Path = args.output
    workdir.mkdir(parents=True, exist_ok=True)
    samples = canonicalise(args.inputs, workdir)

    tools = CATALOGUE.closure(args.until)
    unmeasured = CATALOGUE.unmeasured(args.until)
    known = CATALOGUE.install_size(args.until)
    size = f"{known / 1e9:.1f} GB" if known else "0 GB"
    if unmeasured:
        size += f" + {len(unmeasured)} database(s) of unknown size"
    print(f"{len(samples)} assemblies, {len(tools)} tools, databases: {size}", file=sys.stderr)

    if args.tui:
        from .tui import launch

        launch(args.inputs, workdir, args.databases, samples, args.cores)
        return 0

    build = workdir / ".comparem2"
    (build / "envs").mkdir(parents=True, exist_ok=True)
    snakefile = build / "Snakefile"
    overrides = parse_overrides(args.set)
    snakefile.write_text(
        render(CATALOGUE, args.until, workdir, args.databases, samples,
               overrides=overrides,
               launcher=args.isolated_launcher.split() if args.isolated_launcher else None)
    )
    for name, text in render_envs(CATALOGUE, args.until).items():
        (build / "envs" / name).write_text(text)

    if not args.report_only:
        cmd = [
            "snakemake", "--snakefile", str(snakefile), "--cores", str(args.cores),
            "--rerun-incomplete",
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

    report = render_report(CATALOGUE, args.until, workdir, args.databases, samples)
    print(f"report: {report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
