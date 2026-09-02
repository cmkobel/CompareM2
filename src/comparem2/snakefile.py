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
        files.append(str(tool.database.ready_path(ctx.databases)))

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
          registry: Registry, per_rule_conda: bool,
          overrides: dict[str, tuple[tuple[str, str], ...]] | None = None,
          launcher: Sequence[str] | None = None) -> str:
    ctx = _ctx(workdir, databases, samples, tool, overrides)
    outputs = [str(o) for o in tool.outputs(ctx)]

    argv = list(tool.command(ctx))
    if tool.isolated and launcher:
        # An isolated tool lives in its own environment, so its command needs a
        # prefix that enters that environment. Templated on {tool} so the same
        # mechanism serves pixi (`pixi run -e {tool}`), conda (`conda run -n
        # {tool}`) or a container runtime.
        argv = [part.replace("{tool}", tool.name) for part in launcher] + argv
    command = " ".join(shlex.quote(a) for a in argv)
    if tool.stdout_to_output:
        if len(outputs) != 1:
            raise ValueError(
                f"{tool.name}: stdout_to_output needs exactly one output, got {len(outputs)}"
            )
        command += f" > {shlex.quote(outputs[0])}"

    # shlex.quote escapes the braces of {sample}; unescape so Snakemake sees it.
    for token in (WILDCARD, THREADS):
        command = command.replace(f"'{token}'", token)

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
        # v3 targets one solved environment with every tool on PATH, so no
        # conda directive by default. Opt in only to isolate a tool that
        # cannot co-solve — the situation antiSMASH would have created.
        *([f"    conda: {_q('envs/' + tool.name + '.yaml')}"]
          if (per_rule_conda or tool.isolated) else []),
        "    shell:",
        f'        """',
        f"        mkdir -p $(dirname {{log}})",
        f"        exec > {{log}} 2>&1",
        *[f"        export {name}={shlex.quote(value)}"
          for name, value in (tool.env(ctx) if tool.env else ())],
        f"        mkdir -p {shellify(' '.join(_dirnames(outputs)))}",
        f"        {shellify(command)}",
        f'        """',
        "",
    ]
    return "\n".join(lines)


def _download_rule(db, databases: Path) -> str:
    """One rule per database, whose output is the file that proves it arrived.

    Downloads are rules rather than a separate `--download` phase so that they
    inherit what Snakemake already does well: a database that is present is
    skipped, a half-finished one is redone, and fetching runs alongside work
    that does not depend on it.
    """
    ready = str(db.ready_path(databases))
    steps = db.fetch(databases) if db.fetch else []
    if not steps:
        raise ValueError(f"database {db.name} declares no fetch")
    return "\n".join([
        f"rule download_{db.name.replace('-', '_')}:",
        "    output:",
        f"        {_q(ready)},",
        f"    log: {_q(str(databases / 'logs' / ('download_' + db.name + '.log')))}",
        "    threads: 1",
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
           per_rule_conda: bool = False,
           overrides: dict[str, tuple[tuple[str, str], ...]] | None = None,
           launcher: Sequence[str] | None = None) -> str:
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
    body = [_download_rule(db, databases) for db in registry.databases(selected)]
    body += [_rule(t, workdir, databases, samples, registry, per_rule_conda,
                   overrides, launcher)
             for t in tools]
    return "\n".join(header) + "\n" + "\n".join(body)


def render_envs(registry: Registry, selected: list[str] | None) -> dict[str, str]:
    """One conda environment file per tool.

    v3 targets a single solved environment, but Snakemake wants a file per rule
    and writing them per tool keeps the option of splitting one out later —
    which is exactly what antiSMASH would have needed.
    """
    envs = {}
    for tool in registry.closure(selected):
        packages = "\n".join(f"  - {p}" for p in tool.conda)
        envs[f"{tool.name}.yaml"] = (
            "channels:\n  - conda-forge\n  - bioconda\ndependencies:\n" + packages + "\n"
        )
    return envs


def prepare(registry: Registry, selected: list[str] | None, workdir: Path,
            databases: Path, samples: tuple[str, ...],
            overrides: dict[str, tuple[tuple[str, str], ...]] | None = None,
            launcher: Sequence[str] | None = None) -> Path:
    """Write the Snakefile and its env files; return the Snakefile path.

    Every entry point goes through here. They did not once: the TUI wrote the
    Snakefile itself and skipped the env files, so the `conda:` directive on
    checkm2's rule pointed at a file that did not exist. Snakemake ran the job,
    then failed *recording metadata* for it with a WorkflowError — which aborts
    the workflow rather than the job, so twelve tools with no relation to
    checkm2 never started. A missing 200-byte file killed a 13-tool run.
    """
    build = workdir / ".comparem2"
    (build / "envs").mkdir(parents=True, exist_ok=True)
    snakefile = build / "Snakefile"
    snakefile.write_text(render(registry, selected, workdir, databases, samples,
                                overrides=overrides, launcher=launcher))
    for name, text in render_envs(registry, selected).items():
        (build / "envs" / name).write_text(text)
    return snakefile
