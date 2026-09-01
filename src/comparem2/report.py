"""The report.

Once the tools are commodity, this is the product — so it has two rules.

First, **every tool gets a section**. A generic table renderer means a tool
cannot silently produce nothing, which is how v2 ended up running antiSMASH,
InterProScan, IQ-TREE, FastTree and TreeCluster without displaying any of them.
A specific renderer is an improvement on the fallback, never a prerequisite.

Second, **sections appear only when their outputs exist**, so a partial run
still produces a readable document.

Output is one self-contained HTML file with no external assets, so it survives
being emailed, copied off a cluster, or opened offline.
"""

from __future__ import annotations

import csv
import html
from pathlib import Path

from .tools import Context, Registry, Scope, Tool

CSS = """
:root { --fg:#1a1a1a; --mut:#666; --line:#e3e3e3; --accent:#2b6cb0; --bg:#fff; }
@media (prefers-color-scheme: dark) {
  :root { --fg:#e8e8e8; --mut:#9a9a9a; --line:#333; --accent:#7fb3e8; --bg:#161616; }
}
* { box-sizing: border-box; }
body { font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
       color: var(--fg); background: var(--bg); margin: 0 auto; padding: 3rem 1.5rem;
       max-width: 62rem; }
h1 { font-size: 1.75rem; margin: 0 0 .25rem; letter-spacing: -.02em; }
h2 { font-size: 1.15rem; margin: 2.75rem 0 .2rem; letter-spacing: -.01em; }
.sub { color: var(--mut); margin: 0 0 2rem; }
.summary { color: var(--mut); margin: 0 0 .9rem; font-size: .93rem; }
table { border-collapse: collapse; width: 100%; font-size: .87rem; font-variant-numeric: tabular-nums; }
th, td { text-align: left; padding: .4rem .6rem; border-bottom: 1px solid var(--line); }
th { font-weight: 600; color: var(--mut); font-size: .8rem; text-transform: uppercase;
     letter-spacing: .04em; }
td.n { text-align: right; }
pre { background: color-mix(in srgb, var(--fg) 5%, transparent); padding: .8rem;
      border-radius: 6px; overflow-x: auto; font-size: .8rem; }
.missing { color: var(--mut); font-style: italic; }
footer { margin-top: 4rem; color: var(--mut); font-size: .8rem;
         border-top: 1px solid var(--line); padding-top: 1rem; }
"""


def _read_tsv(path: Path, limit: int | None = None) -> list[list[str]]:
    with path.open() as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    return rows[:limit] if limit else rows


def _table(rows: list[list[str]], header: list[str] | None = None) -> str:
    if not rows:
        return '<p class="missing">No rows.</p>'
    head = header or rows[0]
    body = rows if header else rows[1:]
    cells = "".join(f"<th>{html.escape(c)}</th>" for c in head)
    out = [f"<table><thead><tr>{cells}</tr></thead><tbody>"]
    for row in body:
        tds = "".join(
            f'<td class="n">{html.escape(c)}</td>' if _numeric(c) else f"<td>{html.escape(c)}</td>"
            for c in row
        )
        out.append(f"<tr>{tds}</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def _numeric(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


# --- Tool-specific sections ----------------------------------------
# Anything without an entry here falls back to a plain table of its output.


def _section_seqkit(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Assembly statistics, derived from the per-contig table."""
    rows = [["Genome", "Contigs", "Total length", "Largest", "N50", "GC %"]]
    for sample in ctx.samples:
        path = ctx.sample_out(sample, "seqkit", "contigs.tsv")
        if not path.exists():
            continue
        lengths, gc = [], []
        for rec in _read_tsv(path):
            if len(rec) >= 3 and _numeric(rec[1]):
                lengths.append(int(float(rec[1])))
                gc.append(float(rec[2]))
        if not lengths:
            continue
        lengths.sort(reverse=True)
        total = sum(lengths)
        run, n50 = 0, lengths[-1]
        for length in lengths:
            run += length
            if run >= total / 2:
                n50 = length
                break
        mean_gc = sum(g * n for g, n in zip(gc, lengths)) / total if total else 0
        rows.append([sample, f"{len(lengths)}", f"{total:,}", f"{lengths[0]:,}",
                     f"{n50:,}", f"{mean_gc:.1f}"])
    return _table(rows[1:], header=rows[0]) if len(rows) > 1 else '<p class="missing">No data.</p>'


def _section_mashtree(tool: Tool, ctx: Context, workdir: Path) -> str:
    path = ctx.out("mashtree", "mashtree.newick")
    text = path.read_text().strip() if path.exists() else ""
    return (
        "<p class=\"summary\">Newick, rendered as text until the tree viewer exists.</p>"
        f"<pre>{html.escape(text)}</pre>"
    )


SECTIONS = {
    "seqkit": _section_seqkit,
    "mashtree": _section_mashtree,
}


def _fallback(tool: Tool, ctx: Context, workdir: Path) -> str:
    """A plain table of the tool's first output. Every tool gets at least this."""
    outputs = list(tool.outputs(ctx))
    for path in outputs:
        if path.exists() and path.suffix in {".tsv", ".txt", ".Rtab"}:
            rows = _read_tsv(path, limit=51)
            note = ""
            if len(rows) == 51:
                rows = rows[:50]
                note = '<p class="summary">First 50 rows.</p>'
            return note + _table(rows)
    shown = ", ".join(html.escape(p.name) for p in outputs)
    return f'<p class="missing">Produced {shown}; no table renderer yet.</p>'


def render_report(registry: Registry, selected: list[str] | None, workdir: Path,
                  databases: Path, samples: tuple[str, ...],
                  title: str | None = None) -> Path:
    """Write the report and return its path."""
    parts = [
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
        f"<title>{html.escape(title or workdir.name)}</title><style>{CSS}</style></head><body>",
        f"<h1>{html.escape(title or workdir.name)}</h1>",
        f'<p class="sub">{len(samples)} assemblies</p>',
    ]

    shown = 0
    for tool in registry.closure(selected):
        ctx = Context(workdir, databases, tool.threads, samples,
                      sample=samples[0] if tool.scope is Scope.GENOME else None)
        produced = [p for p in tool.outputs(ctx) if p.exists()]
        if tool.scope is Scope.GENOME:
            produced = [
                p for s in samples
                for p in tool.outputs(Context(workdir, databases, tool.threads, samples, s))
                if p.exists()
            ]
        if not produced:
            continue  # partial runs stay readable
        shown += 1
        renderer = SECTIONS.get(tool.name, _fallback)
        parts += [
            f"<h2>{html.escape(tool.name)}</h2>",
            f'<p class="summary">{html.escape(tool.summary)}</p>',
            renderer(tool, ctx, workdir),
        ]

    if not shown:
        parts.append('<p class="missing">No results yet.</p>')
    parts.append(
        f"<footer>CompareM2 v3 — {shown} of {len(registry.closure(selected))} "
        "tools produced output.</footer></body></html>"
    )

    path = workdir / "report.html"
    path.write_text("".join(parts))
    return path
