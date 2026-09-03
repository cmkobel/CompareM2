"""Generate a Snakefile from `Tool` specs.

This is the load-bearing idea of v3: tools are declared once, and the workflow
is derived. Snakemake still owns everything it is good at — DAG resolution,
resumability, retries, SLURM submission — but no rule is written by hand, so
adding a tool cannot mean editing a workflow.

The trick that makes it work is that `Context` renders paths from a sample
name, so building one with `sample="{sample}"` yields Snakemake wildcards from
exactly the same code that yields real paths at report time.
"""

from __future__ import annotations

import shlex
from collections.abc import Sequence
from pathlib import Path

from .tools import Context, Registry, Scope, Tool

WILDCARD = "{sample}"
THREADS = "{threads}"


def _ctx(workdir: Path, databases: Path, samples: tuple[str, ...], tool: Tool,
         overrides: dict[str, tuple[tuple[str, str], ...]] | None = None) -> Context:
    """A Context that renders wildcard paths for per-genome tools."""
    overrides = overrides or {}
    return Context(
        workdir=workdir,
        databases=databases,
        threads=THREADS,
        samples=samples,
        sample=WILDCARD if tool.scope is Scope.GENOME else None,
        params=overrides.get(tool.name, tool.params),
    )


def _inputs(ctx: Context, tool: Tool, registry: Registry) -> list[str]:
    """Files this tool needs: its assembly, plus its dependencies' outputs.

    A per-genome tool depends on its own sample's copy of an upstream output.
    A set-scope tool that depends on a per-genome tool needs *every* sample's
    copy, so the wildcard is expanded here rather than left to Snakemake.

    A tool's database is an input too, which is what makes downloads part of
    the DAG rather than a separate phase: Snakemake fetches it once, in
    parallel with unrelated work, and skips it when it is already there.
    """
    files = [str(ctx.assembly)] if tool.scope is Scope.GENOME else [str(a) for a in ctx.assemblies]

    if tool.database is not None:
        files.append(str(tool.database.ready_path(ctx.databases, ctx.workdir)))

    # Files the tool declares and `prepare()` writes. Declaring them as inputs
    # is what stops a rule from depending on something invisible: GTDB-Tk's
    # batchfile was named on its command line by a rule that never created it.
    if tool.files is not None:
        files += [str(p) for p in tool.files(ctx)]

    for dep_name in tool.needs:
        dep = registry[dep_name]
        if dep.scope is Scope.GENOME and tool.scope is Scope.SET:
            for sample in ctx.samples:
                dep_ctx = Context(ctx.workdir, ctx.databases, dep.threads, ctx.samples, sample)
                files += [str(o) for o in dep.outputs(dep_ctx)]
        else:
            dep_ctx = Context(ctx.workdir, ctx.databases, dep.threads, ctx.samples, ctx.sample)
            files += [str(o) for o in dep.outputs(dep_ctx)]
    return files


def _rule(tool: Tool, workdir: Path, databases: Path, samples: tuple[str, ...],
          registry: Registry,
          overrides: dict[str, tuple[tuple[str, str], ...]] | None = None) -> str:
    ctx = _ctx(workdir, databases, samples, tool, overrides)
    outputs = [str(o) for o in tool.outputs(ctx)]

    command = " ".join(shlex.quote(a) for a in tool.command(ctx))
    if tool.stdout_to_output:
        if len(outputs) != 1:
            raise ValueError(
                f"{tool.name}: stdout_to_output needs exactly one output, got {len(outputs)}"
            )
        command += f" > {shlex.quote(outputs[0])}"

    # Steps that run after the tool, to turn what it writes into what the spec
    # declared. Quoted the same way, so a glob reaches the step as a literal
    # rather than being expanded by the rule's shell against the pre-run
    # directory — which is empty of the files it is meant to match.
    post = [" ".join(shlex.quote(a) for a in step)
            for step in (tool.post(ctx) if tool.post else ())]

    # shlex.quote escapes the braces of {sample}; unescape so Snakemake sees it.
    for token in (WILDCARD, THREADS):
        command = command.replace(f"'{token}'", token)
        post = [p.replace(f"'{token}'", token) for p in post]

    log = str(ctx.out("logs", f"{tool.name}.log"))

    # Inside a shell block Snakemake exposes wildcards as `wildcards.sample`,
    # not bare `sample` — bare braces are only resolved in input/output paths.
    def shellify(text: str) -> str:
        return text.replace(WILDCARD, "{wildcards.sample}")

    lines = [
        f"rule {tool.name.replace('-', '_')}:",
        "    input:",
        *[f"        {_q(f)}," for f in _inputs(ctx, tool, registry)],
        "    output:",
        *[f"        {_q(o)}," for o in outputs],
        f"    log: {_q(log)}",
        f"    threads: {tool.threads}",
        # Unconditional: Snakemake deploying the tool is the only way a tool
        # arrives. There is no mode in which the pipeline expects one on PATH.
        # Named by *environment*, not by tool — thirteen rules point at the same
        # file, which is what makes them one environment on disk.
        f"    conda: {_q('envs/' + tool.environment + '.yaml')}",
        "    shell:",
        f'        """',
        f"        mkdir -p $(dirname {{log}})",
        f"        exec > {{log}} 2>&1",
        *[f"        export {name}={shlex.quote(value)}"
          for name, value in (tool.env(ctx) if tool.env else ())],
        f"        mkdir -p {shellify(' '.join(_dirnames(outputs)))}",
        f"        {shellify(command)}",
        *[f"        {shellify(step)}" for step in post],
        f'        """',
        "",
    ]
    return "\n".join(lines)


def _download_rule(db, databases: Path, workdir: Path) -> str:
    """One rule per database, whose output is the file that proves it arrived.

    Downloads are rules rather than a separate `--download` phase so that they
    inherit what Snakemake already does well: a database that is present is
    skipped, a half-finished one is redone, and fetching runs alongside work
    that does not depend on it.

    Being a rule has a second consequence: a download needs an environment of
    its own, because two of these fetches run a tool binary (`bakta_db`,
    `amrfinder -u`) that the pipeline's own environment does not contain, and
    the other two need curl and tar, which nothing guarantees either.
    """
    ready = str(db.ready_path(databases, workdir))
    # The fetch is handed the same root its marker lives under, so an
    # out-of-tree database writes its stamp beside this run's results rather
    # than into a shared directory it does not actually populate.
    root = db.marker_root(databases, workdir)
    steps = db.fetch(root) if db.fetch else []
    if not steps:
        raise ValueError(f"database {db.name} declares no fetch")
    if not db.conda:
        raise ValueError(
            f"database {db.name} declares no conda packages, so its download "
            "rule cannot be given an environment")
    return "\n".join([
        f"rule download_{db.name.replace('-', '_')}:",
        "    output:",
        f"        {_q(ready)},",
        f"    log: {_q(str(db.marker_root(databases, workdir) / 'logs' / ('download_' + db.name + '.log')))}",
        "    threads: 1",
        f"    conda: {_q('envs/' + db.environment + '.yaml')}",
        "    shell:",
        '        """',
        "        mkdir -p $(dirname {log})",
        "        exec > {log} 2>&1",
        *[f"        {' '.join(shlex.quote(a) for a in step)}" for step in steps],
        '        """',
        "",
    ])


def _q(path: str) -> str:
    """Quote a path for a Snakefile, leaving {sample} as a live wildcard."""
    return '"' + path + '"'


def _dirnames(outputs: list[str]) -> list[str]:
    seen: dict[str, None] = {}
    for o in outputs:
        seen.setdefault(str(Path(o).parent), None)
    return [shlex.quote(d).replace(f"'{WILDCARD}'", WILDCARD) for d in seen]


def render(registry: Registry, selected: list[str] | None, workdir: Path,
           databases: Path, samples: tuple[str, ...],
           overrides: dict[str, tuple[tuple[str, str], ...]] | None = None) -> str:
    """The whole Snakefile, as text."""
    tools = registry.closure(selected)

    targets: list[str] = []
    for tool in tools:
        ctx = _ctx(workdir, databases, samples, tool, overrides)
        if tool.scope is Scope.GENOME:
            for sample in samples:
                s_ctx = Context(workdir, databases, tool.threads, samples, sample)
                targets += [str(o) for o in tool.outputs(s_ctx)]
        else:
            targets += [str(o) for o in tool.outputs(ctx)]

    header = [
        "# Generated by comparem2 — do not edit.",
        "# Every rule below is derived from src/comparem2/catalogue.py.",
        "",
        "rule all:",
        "    input:",
        *[f"        {_q(t)}," for t in targets],
        "",
    ]
    body = [_download_rule(db, databases, workdir)
            for db in registry.databases(selected)]
    body += [_rule(t, workdir, databases, samples, registry, overrides)
             for t in tools]
    return "\n".join(header) + "\n" + "\n".join(body)


def declared_files(registry: Registry, selected: list[str] | None, workdir: Path,
                   databases: Path, samples: tuple[str, ...]) -> dict[Path, str]:
    """Every `Tool.files` entry, with real paths rather than wildcards.

    The rule's `input:` gets these as `{sample}` templates from the same
    callable; this is the concrete side of that, one entry per sample for a
    per-genome tool.
    """
    out: dict[Path, str] = {}
    for tool in registry.closure(selected):
        if tool.files is None:
            continue
        contexts = ([Context(workdir, databases, tool.threads, samples, s)
                     for s in samples]
                    if tool.scope is Scope.GENOME else
                    [Context(workdir, databases, tool.threads, samples, None)])
        for ctx in contexts:
            out.update(tool.files(ctx))
    return out


def _env_file(packages: Sequence[str]) -> str:
    """A conda environment file, as text.

    Byte-for-byte stability matters beyond tidiness: Snakemake addresses a
    deployed environment by md5(realpath(envs_dir) + this content), so two rules
    whose env files are identical share one environment on disk — which is how
    amrfinder's download rule and its analysis rules end up looking at the same
    $CONDA_PREFIX, and how the same environment is reused across runs.
    """
    return ("channels:\n  - conda-forge\n  - bioconda\ndependencies:\n"
            + "\n".join(f"  - {p}" for p in packages) + "\n")


def render_envs(registry: Registry, selected: list[str] | None) -> dict[str, str]:
    """One conda environment file per *environment*, not per rule.

    So `len()` of this is the number of environments the user pays for — two for
    the full catalogue — rather than the number of rules that need one, which is
    eighteen. Content addressing would deduplicate identical files anyway; naming
    them by environment means the count is visible in the directory rather than
    inferred from md5s.

    A name carrying two different package lists is a bug, not a merge: one file
    would be written twice and whichever rule rendered last would silently
    decide what the other one ran in. Raised here rather than tested for,
    because this is the only place that could notice.
    """
    envs: dict[str, str] = {}

    def add(name: str, packages: Sequence[str], owner: str) -> None:
        text = _env_file(packages)
        if name in envs and envs[name] != text:
            raise ValueError(
                f"environment {name!r} is declared with two different package "
                f"lists; {owner} disagrees with an earlier rule")
        envs[name] = text

    for tool in registry.closure(selected):
        add(f"{tool.environment}.yaml", tool.conda, tool.name)
    for db in registry.databases(selected):
        if db.conda:
            add(f"{db.environment}.yaml", db.conda, f"database {db.name}")
    return envs


def prepare(registry: Registry, selected: list[str] | None, workdir: Path,
            databases: Path, samples: tuple[str, ...],
            overrides: dict[str, tuple[tuple[str, str], ...]] | None = None) -> Path:
    """Write the Snakefile and its env files; return the Snakefile path.

    Every entry point goes through here. They did not once: the TUI wrote the
    Snakefile itself and skipped the env files, so the `conda:` directive on
    checkm2's rule pointed at a file that did not exist. Snakemake ran the job,
    then failed *recording metadata* for it with a WorkflowError — which aborts
    the workflow rather than the job, so twelve tools with no relation to
    checkm2 never started. A missing 200-byte file killed a 13-tool run.
    """
    for path, text in declared_files(registry, selected, workdir, databases,
                                     samples).items():
        path.parent.mkdir(parents=True, exist_ok=True)
        # Only when it changed. These are rule *inputs*, so rewriting an
        # identical file would move its mtime and re-run the tool that reads it
        # on every invocation — for GTDB-Tk, a re-run measured in hours because
        # a four-line file was touched.
        if not path.exists() or path.read_text() != text:
            path.write_text(text)

    build = workdir / ".comparem2"
    (build / "envs").mkdir(parents=True, exist_ok=True)
    snakefile = build / "Snakefile"
    snakefile.write_text(render(registry, selected, workdir, databases, samples,
                                overrides=overrides))
    for name, text in render_envs(registry, selected).items():
        (build / "envs" / name).write_text(text)
    return snakefile
