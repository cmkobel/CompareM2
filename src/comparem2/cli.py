"""Command-line entry point.

Canonicalises inputs, generates the workflow, hands off to Snakemake, and
renders the report. `--tui` hands the same arguments to the TUI, which shares
`snakefile.prepare()` with this path rather than building the workflow itself.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from . import __version__
from . import demo
from .catalogue import CATALOGUE
from .report import render_report
from .snakefile import prepare, render_envs
from .tools import Context, Scope


def _invocation() -> str:
    """The command that produced this run, for the report's provenance line.

    Quoted, so a path with a space round-trips into something a reader can
    paste back. `sys.argv[0]` is an absolute interpreter path under pixi and
    tells the reader nothing, so it is reduced to the program name.
    """
    return shlex.join([Path(sys.argv[0]).name, *sys.argv[1:]])


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


def default_conda_prefix() -> Path:
    """Where Snakemake deploys the tools' environments.

    Shared and home-relative, for the same reason the database default is:
    solved environments outlive any one run, and a per-run default would
    re-solve and re-download them for every set of genomes.

    It is more load-bearing than it looks. Snakemake addresses a deployed
    environment by md5(realpath(envs_dir) + env file content), so moving this
    directory invalidates every environment — and with them AMRFinder's
    database, which lives inside its environment rather than under
    `--databases`.
    """
    override = os.environ.get("COMPAREM2_CONDA_PREFIX")
    if override:
        return Path(override).expanduser()
    return Path.home() / ".comparem2" / "envs"


def missing_conda() -> str | None:
    """The one preflight left, and the one that actually fires.

    Snakemake deploys every tool, so nothing checks for tools on PATH any more
    — but that hands the whole job to `conda`, which Snakemake looks up by name
    and reports as `Error running conda info` from inside DAG construction,
    naming neither PATH nor what the caller should do. Cost a failed run on
    2026-09-03, where conda existed but was a pixi global outside the
    environment's PATH.

    Strictly `conda`, not mamba or micromamba: Snakemake shells out to
    `conda info --json` regardless of which one solves, so accepting a
    micromamba-only machine here would pass the check and fail the run.
    """
    return None if shutil.which("conda") else "conda"


def any_outputs_exist(selected: list[str] | None, workdir: Path,
                      databases: Path, samples: tuple[str, ...]) -> bool:
    """Did any selected tool leave a complete set of declared outputs?

    The question behind "is there anything to report". Asked of the declared
    outputs rather than of the exit code, because a run can fail one tool and
    finish twelve.
    """
    for tool in CATALOGUE.closure(selected):
        contexts = ([Context(workdir, databases, tool.threads, samples, s)
                     for s in samples]
                    if tool.scope is Scope.GENOME else
                    [Context(workdir, databases, tool.threads, samples, None)])
        for ctx in contexts:
            outputs = list(tool.outputs(ctx))
            if outputs and all(Path(o).exists() for o in outputs):
                return True
    return False


def invocation_dir() -> Path:
    """The directory the command was typed in, which is not always the cwd.

    `pixi run cm2 ...` executes the task from the workspace manifest root, not
    from the shell's directory, so `pixi run cm2 *.fna` in a subdirectory
    arrives as relative paths that resolve against the wrong place: every input
    is reported missing, including ones plainly listed by `ls` in that same
    shell. Pixi sets `$INIT_CWD` to the directory the task was launched from
    (verified against pixi 0.78.0), and honouring it makes a relative path mean
    what it looked like it meant.

    Outside pixi the variable is unset and the cwd is already right.
    """
    init = os.environ.get("INIT_CWD")
    if init and Path(init).is_dir():
        return Path(init)
    return Path.cwd()


def resolve(path: Path, base: Path) -> Path:
    """A user-supplied path made absolute, relative to `base`.

    `base / path` leaves an already-absolute path untouched, so this changes
    nothing except the relative case it exists for.
    """
    return (base / path.expanduser()).resolve()


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


def setup_environments(selected: list[str] | None, databases: Path,
                       conda_prefix: Path) -> int:
    """Deploy the tool environments now, so the first real run does not.

    Snakemake builds a missing environment during DAG construction, before the
    first job — so a first run sits silent for as long as the solves take, and
    on a cluster that is the point at which someone kills it. This is the same
    work, asked for on purpose.

    `--conda-create-envs-only` is Snakemake's own mode, so this creates exactly
    what a run would and nothing else. Two things make it usable as a *setup*
    step, both measured 2026-09-03:

    - **It needs no real data.** A DAG has to be constructible, which means the
      per-sample FASTA at the bottom of it must exist — so a 26-byte stub is
      written and thrown away. Nothing reads it, because nothing runs.
    - **It needs no databases.** Every database path in the DAG is some rule's
      output, so the graph closes with `--databases` pointing at a directory
      that does not exist. Setup therefore works before any of the 62.5 GB
      has been fetched.

    The workdir is temporary and that is safe, because Snakemake hashes a
    deployed environment from the **conda prefix** and the env file's content,
    not from where the env file sits: verified by pointing two runs with
    different `--output` at one prefix and watching the second build nothing.
    `--conda-prefix` is the one argument here that has to match what later runs
    use, since its realpath *is* in the hash.
    """
    envs = render_envs(CATALOGUE, selected)
    print(f"deploying {len(envs)} tool environments in {conda_prefix}",
          file=sys.stderr)

    scratch = Path(tempfile.mkdtemp(prefix="comparem2-setup-"))
    try:
        stub = scratch / "stub.fna"
        stub.write_text(">c1\nACGTACGTACGTACGTACGTAAGCTT\n")
        workdir = scratch / "workdir"
        workdir.mkdir()
        samples = canonicalise([stub], workdir)
        snakefile = prepare(CATALOGUE, selected, workdir, databases, samples)
        return subprocess.run([
            sys.executable, "-m", "snakemake",
            "--snakefile", str(snakefile),
            "--directory", str(workdir),
            "--cores", "1",
            "--software-deployment-method", "conda",
            "--conda-prefix", str(conda_prefix),
            "--conda-create-envs-only",
        ]).returncode
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="comparem2", description="A wide view of a set of assemblies.")
    # `*` rather than `+` because --setup takes none. A bare invocation is
    # checked below and still says what is missing.
    p.add_argument("inputs", nargs="*", type=Path, help="assembly FASTA files")
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
    p.add_argument("--version", action="version", version=f"CompareM2 {__version__}")
    p.add_argument("--conda-prefix", type=Path, default=None,
                   help="where the tools' environments are deployed; shared "
                        f"across runs (default: {default_conda_prefix()}, "
                        "overridden for every run by $COMPAREM2_CONDA_PREFIX)")
    p.add_argument("--setup", action="store_true",
                   help="deploy the tool environments and exit, so the first "
                        "real run does not pay for them; takes no assemblies, "
                        "and needs no databases")
    p.add_argument("--demo", action="store_true",
                   help="run the four database-free analyses over six bundled "
                        "Enterococcus faecium plasmids; takes no assemblies, "
                        "and downloads nothing")
    p.add_argument("--keep-going", action="store_true",
                   help="keep running independent tools after one fails")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--unlock", action="store_true",
                   help="release a stale lock left by a killed run, then exit")
    p.add_argument("--report-only", action="store_true",
                   help="re-render the report from existing outputs")
    args = p.parse_args(argv)

    # Every path the user typed is relative to where they typed it, which under
    # `pixi run` is not the cwd — see invocation_dir().
    base = invocation_dir()

    if args.setup:
        # Rejected rather than ignored: a command naming genomes that are never
        # read should not look like it analysed them.
        if args.inputs:
            raise SystemExit("--setup takes no assemblies: it deploys the tool "
                             "environments and exits")
        if args.dry_run or args.tui or args.report_only or args.unlock:
            raise SystemExit("--setup does not combine with --dry-run, --tui, "
                             "--report-only or --unlock")
        if missing_conda():
            raise SystemExit(
                "not on PATH: conda\n"
                "--setup deploys the tool environments with conda, so conda "
                "has to be available. Install it, or add its directory to "
                "PATH — a pixi global lives in ~/.pixi/bin.")
        return setup_environments(
            args.until,
            resolve(args.databases or default_databases(), base),
            resolve(args.conda_prefix or default_conda_prefix(), base))

    if args.demo:
        # Rejected rather than merged, for --setup's reason: a command naming
        # genomes that are never analysed should not look like it analysed them.
        if args.inputs:
            raise SystemExit("--demo takes no assemblies: it runs the bundled "
                             "ones. Drop --demo to analyse your own.")
        # The assemblies land under the output directory, so `--demo` leaves
        # everything it made in one place a reader can inspect and delete.
        args.inputs = demo.extract(resolve(args.output, base) / "demo_assemblies")
        # Fixed rather than defaulted: the bundled inputs are plasmids, and
        # CheckM2 would report a correct and useless completeness near zero on
        # them. An explicit --until is still honoured, because someone asking
        # for a specific tool on known inputs has a reason.
        if args.until is None:
            args.until = list(demo.TOOLS)
        print(f"demo: {len(args.inputs)} bundled plasmids -> "
              f"{resolve(args.output, base)}\n"
              f"      '{demo.DUPLICATE_AS}' is '{demo.DUPLICATE_OF}' again, so "
              "the two should come out identical", file=sys.stderr)

    if not args.inputs:
        raise SystemExit("no assemblies given — pass one or more FASTA files, "
                         "--demo to run the bundled ones, or --setup to deploy "
                         "the tool environments")

    inputs = [resolve(i, base) for i in args.inputs]

    # Report the resolved path, not what was typed: "no such file: 116_2.fna"
    # is unhelpful precisely when the file is sitting there in the directory
    # the user is looking at, and the absolute form says where we did look.
    missing = [str(i) for i in inputs if not i.exists()]
    if missing:
        raise SystemExit(f"no such file: {', '.join(missing)}")

    # Absolute, because Snakemake is given this as its working directory (see
    # below) and every generated path has to survive that.
    workdir: Path = resolve(args.output, base)
    workdir.mkdir(parents=True, exist_ok=True)
    databases: Path = resolve(args.databases or default_databases(), base)
    conda_prefix: Path = resolve(args.conda_prefix or default_conda_prefix(), base)
    samples = canonicalise(inputs, workdir)

    tools = CATALOGUE.closure(args.until)
    print(f"{len(samples)} assemblies, {len(tools)} tools", file=sys.stderr)

    # Say it now rather than from inside DAG construction. Skipped where nothing
    # will be deployed: --report-only reads output that already exists and
    # --unlock only clears a lock. A --dry-run *is* checked, because Snakemake
    # queries conda while building the DAG.
    if not (args.report_only or args.unlock) and missing_conda():
        raise SystemExit(
            "not on PATH: conda\n"
            "CompareM2 runs each tool in an environment Snakemake deploys, so "
            "conda has to be available. Install it, or add its directory to "
            "PATH — a pixi global lives in ~/.pixi/bin.")

    # Only report what is actually going to be fetched. Announcing a total that
    # includes databases already on disk is how "databases: 143.2 GB" came to
    # be printed before a run that downloaded nothing at all.
    pending = [db for db in CATALOGUE.databases(args.until)
               if not db.ready_path(databases, workdir).exists()]
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

    # One file per environment, so this is the number Snakemake will build: two
    # for the full catalogue. Worth saying because the first run pays for both
    # solves and every later run pays for neither.
    envs = len(render_envs(CATALOGUE, args.until))
    print(f"tool environments: {envs} in {conda_prefix}"
          + ("" if conda_prefix.exists() else " (none built yet)"),
          file=sys.stderr)

    overrides = parse_overrides(args.set)

    if args.tui:
        if args.dry_run:
            raise SystemExit(
                "--dry-run and --tui do not combine: the tool list is already "
                "the dry run, and it shows the download size too")
        from .tui import launch

        launch(inputs, workdir, databases, samples, args.cores,
               selected=args.until, overrides=overrides,
               keep_going=args.keep_going,
               conda_prefix=conda_prefix, command=_invocation())
        return 0

    snakefile = prepare(CATALOGUE, args.until, workdir, databases, samples,
                        overrides=overrides)

    if args.unlock:
        # Snakemake locks the output directory, and a run that died without
        # releasing it — SIGKILL, a lost node, a power cut — leaves the next
        # one refusing to start. Snakemake's own message names `--unlock`, but
        # it was not a flag CompareM2 had, so the only way out was to know that
        # a Snakefile sits in `<output>/.comparem2/` and to invoke Snakemake
        # against it by hand. Found by killing a 60.8 GB download.
        return subprocess.run([sys.executable, "-m", "snakemake",
                               "--snakefile", str(snakefile),
                               "--directory", str(workdir), "--unlock"]).returncode

    if not args.report_only:
        cmd = [
            # `sys.executable -m`, not a bare `snakemake`: Snakemake is this
            # package's own dependency, so the right one is the one installed
            # beside the interpreter that is running. Calling it by name found
            # nothing at all when a conda-installed `comparem2` was invoked by
            # absolute path without its environment activated — a
            # FileNotFoundError traceback out of subprocess, three lines after
            # announcing what it was about to do.
            sys.executable, "-m", "snakemake",
            "--snakefile", str(snakefile), "--cores", str(args.cores),
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
            # `--software-deployment-method` is the Snakemake 8+ spelling;
            # `--use-conda` survives as a deprecated alias and is not used here.
            # Unconditional: this is the only way a tool arrives.
            "--software-deployment-method", "conda",
            "--conda-prefix", str(conda_prefix),
        ]
        if args.keep_going:
            cmd.append("--keep-going")
        if args.dry_run:
            cmd.append("--dry-run")
        result = subprocess.run(cmd)
        if args.dry_run:
            return result.returncode
        failure = result.returncode
    else:
        failure = 0

    # A partial run still has a product. `--keep-going` exists to finish the
    # tools that can finish, and returning here withheld the report from twelve
    # tools because a thirteenth failed — measured on a real thirteen-tool run
    # where amrfinder died and no report was written at all. The TUI already
    # behaved this way; the CLI did not.
    #
    # The one case that stays silent is a run that produced nothing: a report
    # describing no results is worse than no report, which is the rule the TUI
    # settled on.
    if failure and not any_outputs_exist(args.until, workdir, databases, samples):
        print("nothing ran; no report written", file=sys.stderr)
        return failure

    report = render_report(CATALOGUE, args.until, workdir, databases, samples,
                           command=_invocation())
    print(f"report: {report}", file=sys.stderr)
    if failure:
        print("some tools failed — the report covers the ones that did not",
              file=sys.stderr)
        return failure
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
