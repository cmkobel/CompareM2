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


class _Node:
    __slots__ = ("name", "length", "children", "x", "y")

    def __init__(self, name: str = "", length: float = 0.0) -> None:
        self.name = name
        self.length = length
        self.children: list[_Node] = []
        self.x = 0.0
        self.y = 0.0


def parse_newick(text: str) -> _Node:
    """Parse a newick string into a tree.

    Deliberately dependency-free: pulling in a phylogenetics library to draw a
    tree would undo the install-weight discipline the whole redesign is about.
    """
    pos = 0

    def parse_node() -> _Node:
        nonlocal pos
        node = _Node()
        if text[pos] == "(":
            pos += 1
            while True:
                node.children.append(parse_node())
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                break
        start = pos
        # Stop at structural characters only. `:` must NOT be here — it
        # separates the label from its branch length and belongs to the label.
        while pos < len(text) and text[pos] not in ",);":
            pos += 1
        label = text[start:pos]
        if ":" in label:
            name, _, length = label.partition(":")
            node.name = name.strip("'\"")
            try:
                node.length = float(length)
            except ValueError:
                node.length = 0.0
        else:
            node.name = label.strip("'\"")
        return node

    return parse_node()


def _leaves(node: _Node) -> list[_Node]:
    return [node] if not node.children else [x for c in node.children for x in _leaves(c)]


PALETTE = ["#2b6cb0", "#b7791f", "#2f855a", "#9b2c2c", "#6b46c1", "#0987a0", "#b83280"]


def draw_tree(root: _Node, clusters: dict[str, str] | None = None,
              width: int = 720, row: int = 26) -> str:
    """A rectangular phylogram as inline SVG."""
    for i, leaf in enumerate(_leaves(root)):
        leaf.y = i * row + row / 2

    def place(node: _Node, x: float) -> float:
        node.x = x + node.length
        if node.children:
            for child in node.children:
                place(child, node.x)
            node.y = (node.children[0].y + node.children[-1].y) / 2
        return node.x

    place(root, 0.0)
    depth = max(n.x for n in _walk(root)) or 1.0
    label_room = 260
    scale = (width - label_room) / depth

    colours: dict[str, str] = {}
    if clusters:
        for name in sorted(set(clusters.values())):
            colours[name] = PALETTE[len(colours) % len(PALETTE)]

    parts: list[str] = []
    for node in _walk(root):
        px = (node.x - node.length) * scale
        parts.append(
            f'<path d="M{px:.1f},{node.y:.1f} H{node.x * scale:.1f}" '
            'fill="none" stroke="currentColor" stroke-width="1.4"/>'
        )
        if node.children:
            top, bottom = node.children[0].y, node.children[-1].y
            parts.append(
                f'<path d="M{node.x * scale:.1f},{top:.1f} V{bottom:.1f}" '
                'fill="none" stroke="currentColor" stroke-width="1.4"/>'
            )
        elif node.name:
            cluster = (clusters or {}).get(node.name)
            fill = colours.get(cluster, "currentColor") if cluster else "currentColor"
            tag = f" · {html.escape(cluster)}" if cluster else ""
            parts.append(
                f'<text x="{node.x * scale + 8:.1f}" y="{node.y + 4:.1f}" '
                f'font-size="13" fill="{fill}">{html.escape(node.name)}{tag}</text>'
            )

    height = len(_leaves(root)) * row + row
    return (
        f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" '
        f'role="img" aria-label="phylogenetic tree">{"".join(parts)}</svg>'
    )


def _walk(node: _Node):
    yield node
    for child in node.children:
        yield from _walk(child)


def _clusters(workdir: Path) -> dict[str, str]:
    """TreeCluster assignments, when they exist, keyed by sample name."""
    path = workdir / "treecluster" / "treecluster.tsv"
    if not path.exists():
        return {}
    rows = _read_tsv(path)
    return {r[0]: r[1] for r in rows[1:] if len(r) >= 2}


def _section_tree(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Render whichever newick this tool declares — not a hardcoded path.

    Both mashtree and fasttree produce a newick, so the renderer takes the path
    from the tool's own `outputs`, which is what the contract is for.
    """
    newicks = [p for p in tool.outputs(ctx) if p.suffix == ".newick"]
    path = newicks[0] if newicks else None
    if path is None or not path.exists():
        return '<p class="missing">No tree.</p>'
    text = path.read_text().strip()
    try:
        svg = draw_tree(parse_newick(text), _clusters(workdir))
    except (IndexError, ValueError):
        # A malformed tree should degrade to the newick, not lose the section.
        return f"<pre>{html.escape(text)}</pre>"
    note = ""
    if _clusters(workdir):
        note = '<p class="summary">Leaves coloured by TreeCluster assignment.</p>'
    return note + svg


def _sample_of(value: str, samples: tuple[str, ...]) -> str:
    """Map a tool's file-path column back to a sample name.

    Tools that take file paths echo them back verbatim, which makes for an
    unreadable report. Inputs are canonicalised to samples/<name>/<name>.fna,
    so the stem is the sample name.
    """
    stem = Path(value).stem
    return stem if stem in samples else value


def _section_mlst(tool: Tool, ctx: Context, workdir: Path) -> str:
    path = ctx.out("mlst", "mlst.tsv")
    if not path.exists():
        return '<p class="missing">No results.</p>'
    rows = []
    for rec in _read_tsv(path):
        if not rec:
            continue
        scheme = rec[1] if len(rec) > 1 else ""
        st = rec[2] if len(rec) > 2 else ""
        alleles = ", ".join(rec[3:]) if len(rec) > 3 else ""
        rows.append([_sample_of(rec[0], ctx.samples), scheme, st, alleles])
    return _table(rows, header=["Genome", "Scheme", "Sequence type", "Alleles"])


def _section_skani(tool: Tool, ctx: Context, workdir: Path) -> str:
    """ANI as a shaded matrix — the densest wide-view summary in the report."""
    path = ctx.out("skani", "ani.tsv")
    if not path.exists():
        return '<p class="missing">No results.</p>'
    rows = _read_tsv(path)
    if len(rows) < 2:
        return '<p class="missing">Not enough genomes to compare.</p>'

    names, matrix = [], []
    for rec in rows[1:]:  # first line is the genome count
        if len(rec) < 2:
            continue
        names.append(_sample_of(rec[0], ctx.samples))
        matrix.append([float(v) if _numeric(v) else None for v in rec[1:]])

    if not names:
        return '<p class="missing">No comparisons.</p>'

    values = [v for row in matrix for v in row if v is not None and v < 100]
    low = min(values) if values else 95.0

    head = "".join(f"<th>{html.escape(n)}</th>" for n in names)
    out = [f"<table><thead><tr><th></th>{head}</tr></thead><tbody>"]
    for name, row in zip(names, matrix):
        cells = []
        for value in row:
            if value is None:
                cells.append("<td></td>")
                continue
            # Shade across the observed range so differences stay visible even
            # when every genome is the same species.
            frac = 0.0 if low >= 100 else max(0.0, min(1.0, (value - low) / (100 - low)))
            alpha = 0.08 + 0.5 * frac
            cells.append(
                f'<td class="n" style="background:rgba(43,108,176,{alpha:.2f})">'
                f"{value:.2f}</td>"
            )
        out.append(f"<tr><th>{html.escape(name)}</th>{''.join(cells)}</tr>")
    out.append("</tbody></table>")
    return f'<p class="summary">Shaded from {low:.2f}% to 100%.</p>' + "".join(out)


def _bar(value: float, colour: str, width: int = 90) -> str:
    """A bar drawn with a div, so it survives being emailed as one file."""
    pct = max(0.0, min(100.0, value))
    return (
        f'<span style="display:inline-block;width:{width}px;height:8px;'
        f'background:color-mix(in srgb, currentColor 12%, transparent);'
        f'border-radius:4px;vertical-align:middle;margin-right:.5rem">'
        f'<span style="display:block;width:{pct:.0f}%;height:100%;'
        f'background:{colour};border-radius:4px"></span></span>'
    )


def _section_checkm2(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Completeness and contamination — the 'is this genome usable' answer."""
    path = ctx.out("checkm2", "quality_report.tsv")
    if not path.exists():
        return '<p class="missing">No results.</p>'
    rows = _read_tsv(path)
    if len(rows) < 2:
        return '<p class="missing">No genomes assessed.</p>'
    header = rows[0]
    try:
        i_name = header.index("Name")
        i_comp = header.index("Completeness")
        i_cont = header.index("Contamination")
    except ValueError:
        return _table(rows)

    out = ["<table><thead><tr><th>Genome</th><th>Completeness</th>"
           "<th>Contamination</th></tr></thead><tbody>"]
    for rec in rows[1:]:
        if len(rec) <= max(i_name, i_comp, i_cont):
            continue
        comp = float(rec[i_comp]) if _numeric(rec[i_comp]) else 0.0
        cont = float(rec[i_cont]) if _numeric(rec[i_cont]) else 0.0
        # MIMAG-style reading: high quality is >90% complete, <5% contaminated.
        comp_colour = "#2f855a" if comp >= 90 else ("#b7791f" if comp >= 50 else "#9b2c2c")
        cont_colour = "#2f855a" if cont < 5 else ("#b7791f" if cont < 10 else "#9b2c2c")
        out.append(
            f"<tr><td>{html.escape(rec[i_name])}</td>"
            f'<td class="n">{_bar(comp, comp_colour)}{comp:.1f}%</td>'
            f'<td class="n">{_bar(min(cont * 5, 100), cont_colour)}{cont:.2f}%</td></tr>'
        )
    out.append("</tbody></table>")
    return ('<p class="summary">Green is high quality by MIMAG conventions: '
            "&gt;90% complete, &lt;5% contaminated. Contamination bars are "
            'scaled &times;5 so small values stay visible.</p>') + "".join(out)


def _section_bakta(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Feature counts per genome, read straight from the GFF3."""
    rows = []
    for sample in ctx.samples:
        gff = ctx.sample_out(sample, "bakta", f"{sample}.gff3")
        if not gff.exists():
            continue
        counts: dict[str, int] = {}
        with gff.open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) > 2:
                    counts[fields[2]] = counts.get(fields[2], 0) + 1
        rows.append([
            sample,
            f"{counts.get('CDS', 0):,}",
            f"{counts.get('tRNA', 0):,}",
            f"{counts.get('rRNA', 0):,}",
            f"{sum(v for k, v in counts.items() if k not in ('CDS', 'tRNA', 'rRNA', 'region')):,}",
        ])
    if not rows:
        return '<p class="missing">No annotations.</p>'
    return _table(rows, header=["Genome", "CDS", "tRNA", "rRNA", "Other features"])


def _section_amrfinder(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Resistance classes per genome, as a presence matrix.

    Per-genome tables would mean one section per sample; across a wide set the
    useful question is which classes appear where.
    """
    per_genome: dict[str, dict[str, int]] = {}
    classes: dict[str, None] = {}
    for sample in ctx.samples:
        path = ctx.sample_out(sample, "amrfinder", "amrfinder.tsv")
        if not path.exists():
            continue
        rows = _read_tsv(path)
        if not rows:
            continue
        header = rows[0]
        try:
            i_class = header.index("Class")
        except ValueError:
            continue
        counts: dict[str, int] = {}
        for rec in rows[1:]:
            if len(rec) > i_class and rec[i_class]:
                name = rec[i_class]
                counts[name] = counts.get(name, 0) + 1
                classes.setdefault(name, None)
        per_genome[sample] = counts

    if not per_genome:
        return '<p class="missing">No resistance genes detected.</p>'

    order = sorted(classes)
    head = "".join(f"<th>{html.escape(c.title())}</th>" for c in order)
    out = [f"<table><thead><tr><th>Genome</th><th>Total</th>{head}</tr></thead><tbody>"]
    for sample, counts in per_genome.items():
        cells = "".join(
            f'<td class="n">{counts[c]}</td>' if c in counts else '<td class="n">·</td>'
            for c in order
        )
        out.append(
            f"<tr><td>{html.escape(sample)}</td>"
            f'<td class="n">{sum(counts.values())}</td>{cells}</tr>'
        )
    out.append("</tbody></table>")
    return "".join(out)


def draw_pangenome(names: list[str], patterns: list[tuple[tuple[bool, ...], int]],
                   width: int = 720, row: int = 22, label_room: int = 150) -> str:
    """The pangenome presence matrix, as inline SVG.

    Genes sharing a presence pattern are drawn as one block whose width is
    proportional to how many genes share it. That is lossless — with N genomes
    there are at most 2**N - 1 patterns — and it keeps the figure small for a
    set where a per-gene column would mean tens of thousands of rects.

    Blocks run from most to least common, so the core sits on the left and the
    genome-specific genes on the right.
    """
    total = sum(count for _, count in patterns) or 1
    span = width - label_room
    height = len(names) * row + 26

    parts: list[str] = []
    for i, name in enumerate(names):
        y = i * row
        parts.append(
            f'<text x="{label_room - 8}" y="{y + row / 2 + 4:.1f}" text-anchor="end" '
            f'font-size="12" fill="currentColor">{html.escape(name)}</text>'
        )

    x = float(label_room)
    for pattern, count in patterns:
        w = max(span * count / total, 0.6)  # never let a rare pattern vanish
        for i, present in enumerate(pattern):
            if not present:
                continue
            y = i * row + 2
            parts.append(
                f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" height="{row - 4}" '
                f'fill="#2b6cb0" fill-opacity="0.85"/>'
            )
        x += w

    parts.append(
        f'<line x1="{label_room}" y1="{len(names) * row + 2}" x2="{width}" '
        f'y2="{len(names) * row + 2}" stroke="currentColor" stroke-opacity="0.25"/>'
    )
    parts.append(
        f'<text x="{label_room}" y="{len(names) * row + 18}" font-size="11" '
        f'fill="currentColor" fill-opacity="0.6">core</text>'
    )
    parts.append(
        f'<text x="{width}" y="{len(names) * row + 18}" text-anchor="end" '
        f'font-size="11" fill="currentColor" fill-opacity="0.6">genome-specific</text>'
    )
    return (f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" '
            f'role="img" aria-label="pangenome presence matrix">{"".join(parts)}</svg>')


def _section_panaroo(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Pangenome structure: the presence matrix, the overlaps, and the split."""
    path = ctx.out("panaroo", "gene_presence_absence.Rtab")
    if not path.exists():
        return '<p class="missing">No pangenome.</p>'
    rows = _read_tsv(path)
    if len(rows) < 2:
        return '<p class="missing">No genes.</p>'

    names = rows[0][1:]
    n = len(names)
    if n == 0:
        return '<p class="missing">No genomes in the matrix.</p>'

    tally: dict[tuple[bool, ...], int] = {}
    per_genome = [0] * n
    core = soft = shell_ = cloud = 0
    for rec in rows[1:]:
        pattern = tuple(v.strip() == "1" for v in rec[1:1 + n])
        if len(pattern) != n:
            continue
        tally[pattern] = tally.get(pattern, 0) + 1
        for i, present in enumerate(pattern):
            if present:
                per_genome[i] += 1
        frac = sum(pattern) / n
        if frac >= 0.99:
            core += 1
        elif frac >= 0.95:
            soft += 1
        elif frac >= 0.15:
            shell_ += 1
        else:
            cloud += 1

    total = sum(tally.values())
    # Most-shared first, so the matrix reads core on the left.
    ordered = sorted(tally.items(), key=lambda kv: (-sum(kv[0]), -kv[1]))

    matrix = draw_pangenome(names, ordered)

    # The overlaps, as an explicit table. Ordered by how many genes share the
    # pattern, which is the question "what do these genomes have in common".
    by_size = sorted(tally.items(), key=lambda kv: -kv[1])
    shown, rest = by_size[:12], by_size[12:]
    overlap_rows = []
    for pattern, count in shown:
        members = [names[i] for i, p in enumerate(pattern) if p]
        label = "all genomes" if len(members) == n else ", ".join(members)
        overlap_rows.append([
            label, f"{len(members)}", f"{count:,}", f"{100 * count / total:.1f}%",
        ])
    if rest:
        overlap_rows.append([
            f"other patterns ({len(rest)})", "–", f"{sum(c for _, c in rest):,}",
            f"{100 * sum(c for _, c in rest) / total:.1f}%",
        ])

    overlaps = _table(
        overlap_rows, header=["Present in", "Genomes", "Gene clusters", "Share"])

    summary = _table(
        [["Core (≥99%)", f"{core:,}"], ["Soft core (95–99%)", f"{soft:,}"],
         ["Shell (15–95%)", f"{shell_:,}"], ["Cloud (<15%)", f"{cloud:,}"],
         ["Total", f"{total:,}"]],
        header=["Partition", "Genes"],
    )
    counts = _table(
        [[name, f"{per_genome[i]:,}"] for i, name in enumerate(names)],
        header=["Genome", "Genes in pangenome"],
    )
    return (
        f'<p class="summary">{total:,} gene clusters across {n} genomes, in '
        f"{len(tally)} distinct presence patterns. Each block is one pattern, "
        "its width proportional to the number of genes sharing it.</p>"
        + matrix
        + '<p class="summary">Which genomes share which genes.</p>' + overlaps
        + summary + counts
    )


def _section_carveme(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Model size per genome, counted straight out of the SBML.

    Genome-scale metabolic models are the one capability no comparable pipeline
    offers, so the report should say something about them rather than noting a
    file exists. Counted by streaming the XML — a 4 MB model per genome is not
    worth a parser dependency.
    """
    rows = []
    for sample in ctx.samples:
        path = ctx.sample_out(sample, "carveme", f"{sample}.xml")
        if not path.exists():
            continue
        reactions = species = genes = 0
        with path.open() as fh:
            for line in fh:
                reactions += line.count("<reaction ")
                species += line.count("<species ")
                genes += line.count("<fbc:geneProduct ")
        rows.append([sample, f"{reactions:,}", f"{species:,}", f"{genes:,}"])
    if not rows:
        return '<p class="missing">No models.</p>'
    return (_table(rows, header=["Genome", "Reactions", "Metabolites", "Genes"])
            + '<p class="summary">Draft models in SBML, ready for flux balance '
              "analysis. Not gap-filled for a specific medium.</p>")


SECTIONS = {
    "carveme": _section_carveme,
    "seqkit": _section_seqkit,
    "checkm2": _section_checkm2,
    "bakta": _section_bakta,
    "amrfinder": _section_amrfinder,
    "mashtree": _section_tree,
    "mlst": _section_mlst,
    "skani": _section_skani,
    "panaroo": _section_panaroo,
    "fasttree": _section_tree,  # same phylogram renderer, its own newick
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
