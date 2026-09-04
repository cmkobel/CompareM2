"""The report.

Once the tools are commodity, this is the product — so it has four rules.

First, **every tool gets a section**. A generic table renderer means a tool
cannot silently produce nothing, which is how v2 ended up running antiSMASH,
InterProScan, IQ-TREE, FastTree and TreeCluster without displaying any of them.
A specific renderer is an improvement on the fallback, never a prerequisite.

Second, **sections appear only when their outputs exist**, so a partial run
still produces a readable document. The same rule holds one level down: a
per-genome table lists every genome and marks the ones a tool produced nothing
for, because a row that is simply absent reads as a genome that was never
submitted.

Third, **genomes are listed in one order everywhere** — the order they were
given on the command line. Tools emit their rows in whatever order suits them;
echoing that back means no two sections can be read across.

Fourth, **every view has a form that survives 100 genomes**. A numeric matrix
is unreadable past about a dozen columns and a per-gene column is unrenderable
past a few thousand, so the matrix renderers switch to a heatmap and the
pangenome figure bins to sub-pixel columns. Nothing here may scale as O(genomes
x genes) in elements emitted.

Output is one self-contained HTML file with no external assets, so it survives
being emailed, copied off a cluster, or opened offline.
"""

from __future__ import annotations

import csv
import html
import math
import re
from datetime import datetime, timezone
from itertools import count
from pathlib import Path

from . import __version__
# The panel is data the wrapper and the report have to agree on, so it lives in
# one place. Importable from here because `biosynthesis.py` keeps its `reframed`
# imports inside the functions that solve — this process has no solver.
from .biosynthesis import ABSENT, DE_NOVO, LB, NO_ROUTE, PANEL, UPSTREAM
from .guidance import GUIDANCE, citations
from .tools import Context, Registry, Scope, Tool

# The content column, in CSS pixels: 62rem of body minus 1.5rem of padding on
# each side. Figures are drawn in these units and scaled with `max-width`, so
# one SVG unit is one pixel at full width and a `font-size="12"` in a figure is
# the same size as 12px of body text. They used to be drawn 720 wide and scaled
# up by 1.31, which is why their labels came out larger than the prose.
WIDTH = 944

CSS = """
:root { --fg:#1a1a1a; --mut:#666; --line:#e3e3e3; --accent:#2b6cb0; --bg:#fff;
        --measure:40em; }
@media (prefers-color-scheme: dark) {
  :root { --fg:#e8e8e8; --mut:#9a9a9a; --line:#333; --accent:#7fb3e8; --bg:#161616; }
}
* { box-sizing: border-box; }
body { font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
       color: var(--fg); background: var(--bg); margin: 0 auto; padding: 3rem 1.5rem;
       max-width: 62rem; }

/* Prose is capped at a readable measure while tables and figures keep the full
   column. The body used to be a single 62rem block, which set every sentence at
   a median of 105 characters — around 1.5x the width at which a line is
   comfortable to track back from. Measured in em rather than rem so the cap
   resolves against each element's own size: at a fixed pixel width the .85rem
   citations ran 20 characters longer per line than the body text. */
p, dl, ul, ol { max-width: var(--measure); }

h1 { font-size: 1.75rem; margin: 0 0 .25rem; letter-spacing: -.02em; }
h2 { font-size: 1.15rem; margin: 3rem 0 .2rem; letter-spacing: -.01em;
     scroll-margin-top: 1rem; }
h3 { font-size: .82rem; margin: 1.9rem 0 .5rem; color: var(--mut);
     text-transform: uppercase; letter-spacing: .05em; font-weight: 600; }
.dots { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        letter-spacing: .18em; white-space: nowrap; }
code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .9em; }
.sub { color: var(--mut); margin: 0 0 .15rem; }
.meta { color: var(--mut); margin: 0 0 1.4rem; font-size: .82rem; }
.summary { color: var(--mut); margin: 0 0 .9rem; font-size: .93rem; }
.note { color: var(--mut); font-size: .82rem; margin: .5rem 0 0; }

/* One line per section, so a 13-section document that runs several screens
   deep can be entered anywhere and linked to by section. */
nav.toc { margin: 0 0 1rem; font-size: .84rem; line-height: 1.9; max-width: none; }
nav.toc a { color: var(--accent); text-decoration: none; }
nav.toc a:hover { text-decoration: underline; }
nav.toc .sep { color: var(--line); margin: 0 .45rem; }

/* Tables size to their contents. At width:100% a four-column table spread its
   values across the full column, which put a genome name and its value 700px
   apart with nothing in between. */
.scroll { overflow-x: auto; max-width: 100%; margin: .1rem 0; }
table { border-collapse: collapse; width: auto; max-width: 100%;
        font-size: .87rem; font-variant-numeric: tabular-nums; }
th, td { text-align: left; padding: .4rem .8rem; border-bottom: 1px solid var(--line);
         white-space: nowrap; }
th:first-child, td:first-child { padding-left: 0; }
th:last-child, td:last-child { padding-right: 0; }
/* Except where the cell is shaded: trimming the padding off a background makes
   the last column's numbers look clipped by their own colour. */
table.matrix td:last-child { padding-right: .8rem; }
/* No text-transform: these headers are sometimes the tool's own column names,
   and matrix row labels are genome names. Uppercasing data turned
   `116_2_duplicate` into `116_2_DUPLICATE`, a label that no longer matched the
   name the same genome carried three sections earlier. */
thead th { font-weight: 600; color: var(--mut); font-size: .8rem; letter-spacing: .03em; }
tbody th { font-weight: 500; }
th.n, td.n { text-align: right; }
/* Matrix shading. The ramp lives here rather than baked into each cell so the
   two colour schemes can have different ones — an alpha that reads as a faint
   tint over white is nothing over #161616. The floor cell used to sit at a
   contrast ratio of 1.067 against the dark page, where 1.0 is literally
   indistinguishable; it is 1.115 now, matching light mode's 1.114. Both ramps
   span the same range (top ~2.3 against the background) and every shade stays
   above 4.5:1 against the text printed on it — the binding constraint is amber
   at full shade in dark, which is why the span is .58 and not more. */
:root { --shade-floor:.08; --shade-span:.50; }
@media (prefers-color-scheme: dark) { :root { --shade-floor:.13; --shade-span:.58; } }
.s0{--shade:0} .s1{--shade:.125} .s2{--shade:.25} .s3{--shade:.375} .s4{--shade:.5}
.s5{--shade:.625} .s6{--shade:.75} .s7{--shade:.875} .s8{--shade:1}
td.shade { background: rgba(var(--rgb),
           calc(var(--shade-floor) + var(--shade-span) * var(--shade))); }
path.shade { fill-opacity: calc(var(--shade-floor) + var(--shade-span) * var(--shade)); }
stop.shade { stop-opacity: calc(var(--shade-floor) + var(--shade-span) * var(--shade)); }

/* Cells in a table the report did not design: one long value must not set the
   column width for the whole table. Hovering gives the full text. */
table.clip td, table.clip th { max-width: 28ch; overflow: hidden;
                               text-overflow: ellipsis; }
/* An inline-block inside a nowrap cell gets its own wrapping context, which is
   how the one prose-shaped column in the report (mlst's allele list) wraps
   without letting every numeric column wrap too. */
.wrapcell { display: inline-block; white-space: normal; max-width: 26rem;
            vertical-align: top; }
pre { background: color-mix(in srgb, var(--fg) 5%, transparent); padding: .8rem;
      border-radius: 6px; overflow-x: auto; font-size: .8rem; max-width: 100%; }
.missing { color: var(--mut); font-style: italic; }
svg { max-width: 100%; height: auto; display: block; margin: .2rem 0 .4rem; }
/* Figures are laid out against measured glyph widths for this stack, so they
   have to be drawn in it — an SVG <text> otherwise falls back to the browser
   default sans and every computed gutter is wrong. */
svg text { font-family: inherit; }
footer { margin-top: 4rem; color: var(--mut); font-size: .8rem;
         border-top: 1px solid var(--line); padding-top: 1rem; }

/* The explanatory block. Collapsed by default: someone who already knows the
   tools should not have to scroll past a page of prose to reach their data,
   and someone who does not needs it one click away rather than in a manual. */
details.about { margin: .3rem 0 1.2rem; font-size: .88rem; }
details.about > summary { cursor: pointer; color: var(--accent); font-size: .82rem;
    text-transform: uppercase; letter-spacing: .04em; font-weight: 600;
    list-style: none; padding: .2rem 0; }
details.about > summary::-webkit-details-marker { display: none; }
/* Literal glyphs, not CSS escapes: in `content: "\\25B8 "` the trailing space
   terminates the hex escape instead of being part of the string, which ran the
   marker straight into the word after it. */
details.about > summary::before { content: "▸ "; }
details.about[open] > summary::before { content: "▾ "; }
details.about > div { border-left: 2px solid var(--line); padding: .1rem 0 .1rem 1rem;
    margin-top: .5rem; }
details.about p { margin: .55rem 0; }
details.about .method { color: var(--mut); }
/* h4, not h3: these label a panel that is itself subordinate to the section,
   and h3 is the section's own subheading level. */
details.about h4 { font-size: .78rem; margin: 1.5rem 0 .4rem; color: var(--mut);
    text-transform: uppercase; letter-spacing: .05em; font-weight: 600; }
details.about dl { margin: .55rem 0; }
details.about dt { font-weight: 600; margin-top: .7rem; font-size: .84rem; }
details.about dd { margin: .15rem 0 0; color: var(--mut); max-width: var(--measure); }
details.about ul { margin: .35rem 0; padding-left: 1.1rem; color: var(--mut); }
details.about li { margin: .3rem 0; }
details.about .cite { font-size: .82rem; color: var(--mut); margin-top: .9rem; }
a { color: var(--accent); }
.refs { font-size: .85rem; }
.refs li { margin: .5rem 0; }
.refs .aside { color: var(--mut); font-style: italic; }
"""


def _read_tsv(path: Path, limit: int | None = None) -> list[list[str]]:
    with path.open() as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    return rows[:limit] if limit else rows


def _scroll(table: str) -> str:
    """Wrap a table so an over-wide one scrolls instead of overflowing the page.

    The amrfinder matrix used to run 1,713px wide inside a 992px column because
    one resistance class is named `Lincosamide/Macrolide/Streptogramin`, pushing
    a horizontal scrollbar onto the whole document.
    """
    return f'<div class="scroll">{table}</div>'


# How wide a cell may get in a table the report did not design. A tool without
# a renderer can emit anything: GTDB-Tk's summary is twenty columns, four of
# them a seven-rank taxonomy string, and one long cell sets its column's width
# — 6,161px of table inside a 944px scroll strip. Clipping keeps every column
# reachable and the scroll finite, and the full value stays in the title.
_UNDESIGNED_CELL_CHARS = 28


def _table(rows: list[list[str]], header: list[str] | None = None,
           clip: bool = False) -> str:
    if not rows:
        return '<p class="missing">No rows.</p>'
    head = header or rows[0]
    body = rows if header else rows[1:]
    # A header sits over its column, so it takes the column's alignment. These
    # used to be left-aligned over right-aligned numbers, which put `CONTIGS`
    # 193px to the left of every value beneath it.
    #
    # Column 0 is exempt: in every table here it is an identifier — a genome, a
    # partition, a cluster number — and an identifier that happens to parse as a
    # number is still an identifier, not a quantity to line up the decimal of.
    numeric = {
        i for i in range(1, len(head))
        if body and all(_numeric(r[i]) for r in body if len(r) > i and r[i] not in _ABSENT)
        and any(_numeric(r[i]) for r in body if len(r) > i)
    }
    def full(text: str) -> str:
        """A title carrying the whole value, for a cell the CSS will clip."""
        return f' title="{html.escape(text)}"' if clip and len(text) > _UNDESIGNED_CELL_CHARS else ""

    cells = "".join(
        f'<th class="n"{full(c)}>{html.escape(c)}</th>' if i in numeric
        else f"<th{full(c)}>{html.escape(c)}</th>"
        for i, c in enumerate(head)
    )
    opening = '<table class="clip">' if clip else "<table>"
    out = [f"{opening}<thead><tr>{cells}</tr></thead><tbody>"]
    for row in body:
        tds = "".join(
            f'<td class="n"{full(c)}>{html.escape(c)}</td>'
            if i in numeric or (i and _numeric(c))
            else f"<td{full(c)}>{html.escape(c)}</td>"
            for i, c in enumerate(row)
        )
        out.append(f"<tr>{tds}</tr>")
    out.append("</tbody></table>")
    return _scroll("".join(out))


# What counts as a number *for display purposes* — i.e. what should be
# right-aligned in a table. `float()` is the wrong test in both directions:
# it rejects "2,091" (so formatted numbers ended up left-aligned) and accepts
# "116_2" (Python allows underscores in numeric literals), which right-aligned
# a sample name as though it were 1162.
_NUMBER = re.compile(r"^[-+]?(\d{1,3}(,\d{3})+|\d+)(\.\d+)?%?$")

# What a cell holds when the tool produced nothing for that genome. Kept out of
# the numeric-column test so one absent genome does not left-align a column of
# numbers.
_ABSENT = {"—", "·"}


def _numeric(value: str) -> bool:
    return bool(_NUMBER.match(value.strip()))


def _raw_table(rows: list[list[str]], header: list[str],
               numeric_columns: set[int] | None = None) -> str:
    """A table whose cells already contain markup.

    `_table` escapes everything, which is right for values read out of a tool's
    output. This variant is for cells the report itself builds — bars, dots —
    so callers must escape any tool-derived text they embed.
    """
    numeric_columns = numeric_columns or set()
    cells = "".join(
        f'<th class="n">{html.escape(c)}</th>' if i in numeric_columns
        else f"<th>{html.escape(c)}</th>"
        for i, c in enumerate(header)
    )
    out = [f"<table><thead><tr>{cells}</tr></thead><tbody>"]
    for row in rows:
        tds = "".join(
            f'<td class="n">{c}</td>' if i in numeric_columns else f"<td>{c}</td>"
            for i, c in enumerate(row)
        )
        out.append(f"<tr>{tds}</tr>")
    out.append("</tbody></table>")
    return _scroll("".join(out))


# --- ordering ------------------------------------------------------
# Genomes appear in input order in every section. Tools emit rows in whatever
# order suits them — mlst by scheme, snp-dists by its own matrix, TreeCluster
# by cluster — and echoing that back meant five different row orders across
# eleven tables in a four-genome report.


def _in_sample_order(found: dict[str, list[str]], samples: tuple[str, ...],
                     columns: int, absent_cell: str = "—"
                     ) -> tuple[list[list[str]], int]:
    """Rows for every sample, in input order, with the absent ones marked.

    Returns the rows and how many samples the tool produced nothing for. A
    genome that a tool skipped used to vanish from its table, so a reader
    comparing four genomes saw three rows and no indication why. `absent_cell`
    is the filler: plain for `_table`, which escapes, and markup for
    `_raw_table`, which does not.
    """
    rows: list[list[str]] = []
    absent = 0
    for sample in samples:
        row = found.get(sample)
        if row is None:
            rows.append([sample, *[absent_cell] * (columns - 1)])
            absent += 1
        else:
            rows.append(row)
    for name, row in found.items():  # anything the tool named that we did not submit
        if name not in samples:
            rows.append(row)
    return rows, absent


def _absent_note(count: int, total: int) -> str:
    if not count:
        return ""
    # The noun agrees with the total, not the count: "1 of 3 genomes".
    word = "genome" if total == 1 else "genomes"
    return (f'<p class="note">{count} of {total} {word} produced no output for '
            "this tool; those rows are marked —.</p>")


# --- tool-specific sections ----------------------------------------
# Anything without an entry here falls back to a plain table of its output.


# Viridis, sampled at nine stops and interpolated. Hardcoded because importing
# matplotlib to obtain nine colours would cost more than every analysis tool in
# the pipeline combined.
VIRIDIS = ["#440154", "#482777", "#3f4a8a", "#31678e", "#26838f",
           "#1f9d8a", "#6cce5a", "#b6de2b", "#fee825"]


def viridis(fraction: float) -> str:
    """Colour for a value in [0, 1], linearly interpolated between stops."""
    f = max(0.0, min(1.0, fraction)) * (len(VIRIDIS) - 1)
    i = min(int(f), len(VIRIDIS) - 2)
    t = f - i
    a, b = VIRIDIS[i], VIRIDIS[i + 1]
    channels = (
        round(int(a[1 + 2 * k:3 + 2 * k], 16) * (1 - t)
              + int(b[1 + 2 * k:3 + 2 * k], 16) * t)
        for k in range(3)
    )
    return "#" + "".join(f"{c:02x}" for c in channels)


_ids = count()


def _gid(prefix: str) -> str:
    """A document-unique id, for SVG gradients and their `url(#...)` refs."""
    return f"{prefix}{next(_ids)}"


# Glyph widths in em for the body font stack, measured once with
# `CanvasRenderingContext2D.measureText` at 100px and grouped to the nearest
# 0.02 em — a third of a pixel at the sizes figures use. Figures set
# `font-family: inherit` so they are drawn in the font these describe.
#
# They are metrics for one resolution of that stack, `-apple-system` on macOS;
# a reader on Linux gets DejaVu Sans and slightly different numbers. That is
# why `_text_width` adds slack and why every gutter computed from it is a
# minimum rather than an exact fit.
_EM: dict[str, float] = {}
for _width, _chars in (
    (0.21, " ,.:;ijl|I"),
    (0.27, "!'"),
    (0.28, "/\\"),
    (0.31, "frt"),
    (0.32, "()[]{}"),
    (0.40, '"*'),
    (0.44, "1-"),
    (0.47, "sxz"),
    (0.49, "?Jkvy"),
    (0.50, "`ac"),
    (0.51, "e"),
    (0.53, "L_nuF"),
    (0.55, "7bdghpqoE"),
    (0.59, "2357PSRT"),
    (0.60, "0489BK#$+<=>^~"),
    (0.62, "69AVXYZ"),
    (0.68, "CD&"),
    (0.70, "GHNU"),
    (0.73, "OQw"),
    (0.81, "%m"),
    (0.85, "M@"),
    (0.92, "W"),
):
    for _c in _chars:
        _EM[_c] = _width
_DEFAULT_EM = 0.60  # anything outside printable ASCII


def _text_width(text: str, size: float) -> float:
    """Rendered width of `text` at `size`, in user units.

    The per-character term is not slack: at the sizes figures use, glyph
    advances are quantised and come out consistently wider than the linear
    scale from the 100px metrics — by 0.3 to 0.8px per character across the
    strings this was checked against. Summing the em table alone under-measured
    a genome label by 6%, which is what pushed the leftmost name in the contig
    figure two pixels outside its own viewBox. The 3% on top is the slack.
    """
    ems = sum(_EM.get(ch, _DEFAULT_EM) for ch in text)
    return (ems * size + 0.55 * len(text)) * 1.03


def _clamp(value: float, low: float, high: float) -> float:
    return max(low, min(high, value))


def _svg(body: str, width: int, height: float, label: str) -> str:
    return (f'<svg viewBox="0 0 {width} {height:.0f}" width="{width}" '
            f'height="{height:.0f}" role="img" aria-label="{html.escape(label)}">'
            f"{body}</svg>")


def _bp(value: float) -> str:
    if value >= 1e6:
        return f"{value / 1e6:.3g} Mb"
    if value >= 1e3:
        return f"{value / 1e3:.3g} kb"
    return f"{value:.0f}"


def _ticks(maximum: float, target: int = 4) -> list[float]:
    """Round tick positions at or below `maximum`.

    Quartering the largest genome gives labels like `819.114 kb`. Steps are
    picked from 1/2/2.5/5 x 10^k instead, so the axis reads 0, 1 Mb, 2 Mb, 3 Mb.
    """
    if maximum <= 0:
        return [0.0]

    rough = maximum / max(target, 1)
    magnitude = 10 ** math.floor(math.log10(rough))
    for multiple in (1, 2, 2.5, 5, 10):
        step = multiple * magnitude
        if rough <= step:
            break
    out, value = [], 0.0
    while value <= maximum + step * 1e-9:
        out.append(value)
        value += step
    return out


# One rect per fasta record is right for four genomes and is 20,000 elements
# for a hundred draft assemblies. Past this budget the smallest records in each
# genome are merged into one trailing segment, which keeps their total length
# and mean GC on the figure and their exact count in the table beneath it.
_CONTIG_RECT_BUDGET = 3000


def _merge_tail(contigs: list[tuple[int, float]], keep: int
                ) -> tuple[list[tuple[int, float]], int]:
    """The `keep - 1` longest records in file order, plus the rest as one segment.

    The kept records stay where they sit in the file, which is the whole point
    of the figure; the merged segment is a documented departure appended at the
    end, carrying the combined length and length-weighted mean GC of everything
    too small to draw.
    """
    if len(contigs) <= keep:
        return contigs, 0
    cutoff = sorted((length for length, _ in contigs), reverse=True)[keep - 2]
    kept: list[tuple[int, float]] = []
    rest: list[tuple[int, float]] = []
    for record in contigs:
        if record[0] >= cutoff and len(kept) < keep - 1:
            kept.append(record)
        else:
            rest.append(record)
    total = sum(length for length, _ in rest)
    mean_gc = (sum(gc * length for length, gc in rest) / total) if total else 0.0
    return [*kept, (total, mean_gc)], len(rest)


def draw_contigs(rows: list[tuple[str, float, list[tuple[int, float]]]],
                 width: int = WIDTH) -> str:
    """Contig sizes per genome, coloured by GC content.

    One row per genome, one rect per fasta record, width proportional to the
    record's length and fill given by its GC. Ported from v2's
    `sequence_lengths` figure, which is the fastest way to see a fragmented
    assembly, an outsized genome, or a contig whose GC does not match its
    neighbours.
    """
    if not rows:
        return '<p class="missing">No contigs.</p>'

    n = len(rows)
    row_h = _clamp(round(1600 / n), 8, 18)
    font = _clamp(row_h - 4, 7, 11)
    label_room = min(
        max(_text_width(f"{s} ({gc:.1f}%)", font) for s, gc, _ in rows) + 10,
        width * 0.42,
    )
    longest = max(sum(l for l, _ in contigs) for _, _, contigs in rows) or 1
    gcs = [gc for _, _, contigs in rows for _, gc in contigs]
    lo, hi = (min(gcs), max(gcs)) if gcs else (0.0, 1.0)
    if hi - lo < 1e-9:
        lo, hi = lo - 1, hi + 1

    # Gutter for the colour key: the gap before the bar, the bar, and whichever
    # is wider of its tick labels and the "GC%" caption beneath it. Guessing a
    # constant here clipped the top tick of a two-digit scale.
    legend_gap, legend_bar, tick_gap = 20.0, 11.0, 5.0
    ticks = [f"{value:.0f}" for value in (lo, (lo + hi) / 2, hi)]
    legend_room = max(
        legend_gap + legend_bar + tick_gap + max(_text_width(t, 10) for t in ticks),
        legend_gap + _text_width("GC%", 10),
    ) + 4
    span = width - label_room - legend_room

    budget = max(20, _CONTIG_RECT_BUDGET // n)
    axis_h = 30
    height = n * row_h + axis_h + 8
    parts: list[str] = []
    merged = 0

    for i, (sample, overall_gc, contigs) in enumerate(rows):
        y = i * row_h
        parts.append(
            f'<text x="{label_room - 8:.1f}" y="{y + row_h / 2 + font / 3:.1f}" '
            f'text-anchor="end" font-size="{font:.0f}" fill="currentColor">'
            f"{html.escape(sample)} ({overall_gc:.1f}%)</text>"
        )
        drawn, folded = _merge_tail(contigs, budget)
        merged += folded
        x = float(label_room)
        for length, gc in drawn:
            w = span * length / longest
            parts.append(
                f'<rect x="{x:.2f}" y="{y + 2:.1f}" width="{max(w, 0.4):.2f}" '
                f'height="{row_h - 4:.1f}" fill="{viridis((gc - lo) / (hi - lo))}">'
                f"<title>{html.escape(sample)}: {length:,} bp, {gc:.1f}% GC</title>"
                "</rect>"
            )
            x += w

    # x axis
    base = n * row_h + 4
    parts.append(
        f'<line x1="{label_room:.1f}" y1="{base:.1f}" x2="{label_room + span:.1f}" '
        f'y2="{base:.1f}" stroke="currentColor" stroke-opacity="0.35"/>'
    )
    for value in _ticks(longest):
        x = label_room + span * value / longest
        parts.append(
            f'<line x1="{x:.1f}" y1="{base:.1f}" x2="{x:.1f}" y2="{base + 4:.1f}" '
            'stroke="currentColor" stroke-opacity="0.35"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{base + 16:.1f}" text-anchor="middle" font-size="10" '
            f'fill="currentColor" fill-opacity="0.7">{_bp(value)}</text>'
        )

    # GC legend
    gid = _gid("gc")
    lx = label_room + span + legend_gap
    ly = 4
    lh = min(n * row_h - 8, 150)
    parts.append(
        f'<defs><linearGradient id="{gid}" x1="0" y1="1" x2="0" y2="0">'
        + "".join(
            f'<stop offset="{j / (len(VIRIDIS) - 1):.3f}" stop-color="{c}"/>'
            for j, c in enumerate(VIRIDIS)
        )
        + "</linearGradient></defs>"
    )
    parts.append(
        f'<rect x="{lx:.1f}" y="{ly}" width="{legend_bar:.0f}" height="{lh:.1f}" '
        f'fill="url(#{gid})"/>'
    )
    for frac, text in zip((0.0, 0.5, 1.0), ticks):
        y = ly + lh - lh * frac
        parts.append(
            f'<text x="{lx + legend_bar + tick_gap:.1f}" y="{y + 3.5:.1f}" '
            f'font-size="10" fill="currentColor" fill-opacity="0.7">{text}</text>'
        )
    parts.append(
        f'<text x="{lx:.1f}" y="{ly + lh + 14:.1f}" font-size="10" fill="currentColor" '
        'fill-opacity="0.7">GC%</text>'
    )

    svg = _svg("".join(parts), width, height,
               "contig sizes coloured by GC content")
    if merged:
        svg += (f'<p class="note">The smallest records are folded into one '
                f"trailing segment per genome — {merged:,} in total — carrying "
                "their combined length and mean GC. Exact counts are in the "
                "table below.</p>")
    return svg


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


PALETTE = ["#2b6cb0", "#b7791f", "#2f855a", "#9b2c2c", "#6b46c1", "#0987a0", "#b83280",
           "#3182ce", "#975a16", "#276749", "#742a2a", "#553c9a", "#086f83", "#97266d"]


def draw_tree(root: _Node, clusters: dict[str, str] | None = None,
              width: int = WIDTH) -> str:
    """A rectangular phylogram as inline SVG."""
    leaves = _leaves(root)
    n = len(leaves)
    row = _clamp(round(1700 / max(n, 1)), 11, 26)
    font = _clamp(row - 5, 8, 13)

    for i, leaf in enumerate(leaves):
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

    colours: dict[str, str] = {}
    if clusters:
        for name in sorted(set(clusters.values())):
            colours[name] = PALETTE[len(colours) % len(PALETTE)]

    def label_of(leaf: _Node) -> str:
        cluster = (clusters or {}).get(leaf.name)
        return f"{leaf.name} · {cluster}" if cluster else leaf.name

    # Measured, not reserved: a constant 260-unit gutter left 30% of the figure
    # blank for labels like `116_2` and would have clipped a long strain name.
    label_room = min(
        max((_text_width(label_of(leaf), font) for leaf in leaves), default=0) + 12,
        width * 0.5,
    )
    scale = (width - label_room) / depth

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
                f'<text x="{node.x * scale + 8:.1f}" y="{node.y + font / 3:.1f}" '
                f'font-size="{font:.0f}" fill="{fill}">'
                f"{html.escape(node.name)}{tag}</text>"
            )

    return _svg("".join(parts), width, n * row + row / 2, "phylogenetic tree")


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


def _section_seqkit(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Contig sizes coloured by GC, plus the statistics derived from them."""
    found: dict[str, list[str]] = {}
    figure_rows: list[tuple[str, float, list[tuple[int, float]]]] = []

    for sample in ctx.samples:
        path = ctx.sample_out(sample, "seqkit", "contigs.tsv")
        if not path.exists():
            continue
        # Kept in file order, so the figure shows each fasta record where it
        # actually sits rather than an ordering the report invented.
        contigs = [
            (int(float(rec[1])), float(rec[2]))
            for rec in _read_tsv(path)
            if len(rec) >= 3 and _numeric(rec[1]) and _numeric(rec[2])
        ]
        if not contigs:
            continue

        total = sum(length for length, _ in contigs)
        mean_gc = (sum(gc * length for length, gc in contigs) / total) if total else 0.0
        figure_rows.append((sample, mean_gc, contigs))

        lengths = sorted((length for length, _ in contigs), reverse=True)
        run, n50 = 0, lengths[-1]
        for length in lengths:
            run += length
            if run >= total / 2:
                n50 = length
                break
        found[sample] = [sample, f"{len(lengths)}", f"{total:,}", f"{lengths[0]:,}",
                         f"{n50:,}", f"{mean_gc:.1f}"]

    if not found:
        return '<p class="missing">No data.</p>'

    rows, absent = _in_sample_order(found, ctx.samples, 6)
    return (
        '<p class="summary">Each bar is one genome, each segment one fasta '
        "record, sized by length and coloured by its GC content. Hover a "
        "segment for its length.</p>"
        + draw_contigs(figure_rows)
        + "<h3>Assembly statistics</h3>"
        + _table(rows, header=["Genome", "Contigs", "Total length", "Largest",
                               "N50", "GC %"])
        + _absent_note(absent, len(ctx.samples))
    )


def _section_mlst(tool: Tool, ctx: Context, workdir: Path) -> str:
    path = ctx.out("mlst", "mlst.tsv")
    if not path.exists():
        return '<p class="missing">No results.</p>'
    found: dict[str, list[str]] = {}
    for rec in _read_tsv(path):
        if not rec:
            continue
        name = _sample_of(rec[0], ctx.samples)
        scheme = rec[1] if len(rec) > 1 else ""
        st = rec[2] if len(rec) > 2 else ""
        alleles = ", ".join(rec[3:]) if len(rec) > 3 else ""
        found[name] = [name, scheme, st, alleles]
    rows, absent = _in_sample_order(found, ctx.samples, 4)
    # The allele list is the one cell in the report that should wrap: it is
    # prose-shaped, and nowrap would widen the table by several hundred pixels.
    body = _raw_table(
        [[html.escape(r[0]), html.escape(r[1]), html.escape(r[2]),
          f'<span class="wrapcell">{html.escape(r[3])}</span>']
         for r in rows],
        header=["Genome", "Scheme", "Sequence type", "Alleles"],
        numeric_columns={2},
    )
    return body + _absent_note(absent, len(ctx.samples))


# --- matrices ------------------------------------------------------
# A numeric cell needs about 70px, so a matrix stops fitting the 944px column
# at roughly a dozen genomes. Past that the same data is drawn as a heatmap,
# which is legible to a hundred and beyond.
_NUMERIC_MATRIX_MAXIMUM = 12

# Shades in a heatmap ramp. Nine is what a sequential ColorBrewer scale uses
# and about what the eye can order; it also decides how often neighbouring
# cells merge into one subpath, so raising it costs bytes for resolution
# nobody can read.
_SHADES = 8


def _read_skani_matrix(path: Path, samples: tuple[str, ...]
                       ) -> tuple[list[str], list[list[float | None]]]:
    """A skani triangle matrix: a genome count, then names down the side."""
    names: list[str] = []
    matrix: list[list[float | None]] = []
    for rec in _read_tsv(path)[1:]:  # first line is the genome count
        if len(rec) < 2:
            continue
        names.append(_sample_of(rec[0], samples))
        matrix.append([float(v) if _numeric(v) else None for v in rec[1:]])
    return names, matrix


def _read_labelled_matrix(path: Path, samples: tuple[str, ...]
                          ) -> tuple[list[str], list[list[float | None]]]:
    """A matrix whose first row is a header: snp-dists writes its own version
    into the corner cell, which the generic fallback rendered as a column
    heading reading `snp-dists 1.2.0`."""
    rows = _read_tsv(path)
    if len(rows) < 2:
        return [], []
    names: list[str] = []
    matrix: list[list[float | None]] = []
    for rec in rows[1:]:
        if len(rec) < 2:
            continue
        names.append(_sample_of(rec[0], samples))
        matrix.append([float(v) if _numeric(v) else None for v in rec[1:]])
    return names, matrix


def _reorder(names: list[str], matrix: list[list[float | None]],
             samples: tuple[str, ...]) -> tuple[list[str], list[list[float | None]]]:
    """Permute a square matrix into input order, rows and columns alike."""
    if len(names) != len(matrix) or any(len(r) != len(names) for r in matrix):
        return names, matrix  # not square; leave it as the tool wrote it
    order = [names.index(s) for s in samples if s in names]
    order += [i for i in range(len(names)) if i not in order]
    return ([names[i] for i in order],
            [[matrix[i][j] for j in order] for i in order])


def _numeric_matrix(names: list[str], matrix: list[list[float | None]],
                    rgb: str, fmt, frac) -> str:
    """The matrix as a shaded table of numbers. Small sets only."""
    head = "".join(f'<th class="n">{html.escape(n)}</th>' for n in names)
    # The hue is the only part of the shading the markup carries; how dark each
    # step is belongs to the stylesheet, which is what lets dark mode differ.
    out = [f'<table class="matrix" style="--rgb:{rgb}">'
           f"<thead><tr><th></th>{head}</tr></thead><tbody>"]
    for name, row in zip(names, matrix):
        cells = []
        for value in row:
            if value is None:
                cells.append("<td></td>")
                continue
            shade = round(_clamp(frac(value), 0.0, 1.0) * _SHADES)
            cells.append(f'<td class="n shade s{shade}">{fmt(value)}</td>')
        out.append(f"<tr><th>{html.escape(name)}</th>{''.join(cells)}</tr>")
    out.append("</tbody></table>")
    return _scroll("".join(out))


def _grid(row_labels: list[str], col_labels: list[str],
          fracs: list[list[float | None]], rgb: str, *,
          row_notes: list[str] | None = None,
          cell_values: list[list[str]] | None = None,
          width: int = WIDTH, label: str = "matrix") -> str:
    """A heatmap as inline SVG, emitted as one path per shade.

    A hundred genomes against a hundred genomes is ten thousand cells. As
    `<rect>` elements that is close to a megabyte of markup in a document whose
    whole point is being small enough to email, so shades are quantised to 17
    steps and runs of equal shade along a row are merged into one subpath.
    """
    if not row_labels or not col_labels:
        return '<p class="missing">Nothing to show.</p>'

    font = 10.0
    left = min(max(_text_width(s, font) for s in row_labels) + 10, width * 0.32)
    right = 0.0
    if row_notes:
        right = max(_text_width(s, font) for s in row_notes) + 12
    cell = _clamp((width - left - right) / len(col_labels), 5, 46)
    grid_w = cell * len(col_labels)
    # Cells are wider than they are tall, like the table rows they stand in for.
    row_h = _clamp(cell * 0.62, 9, 26)
    # Column labels are rotated, so a long resistance-class name costs vertical
    # space (cheap) instead of horizontal (which is what overflowed the page).
    show_cols = cell >= 7
    top = (min(max(_text_width(s, font) for s in col_labels) + 10, 260.0)
           if show_cols else 6.0)

    parts: list[str] = []
    for j, name in enumerate(col_labels):
        if not show_cols:
            break
        x = left + cell * j + cell / 2
        parts.append(
            f'<text x="{x:.1f}" y="{top - 6:.1f}" font-size="{min(font, cell - 1):.1f}" '
            f'fill="currentColor" fill-opacity="0.75" text-anchor="start" '
            f'transform="rotate(-90 {x:.1f} {top - 6:.1f})">{html.escape(name)}</text>'
        )

    # Integer cell edges, so a subpath is `M112,300h8v6h-8z` rather than
    # `M112.4,300.2h8.3v6.2h-8.3z`. Snapped rather than rounded independently,
    # so neighbouring cells still tile exactly instead of leaving hairlines.
    def edge(j: int) -> int:
        return int(left) + round(cell * j)

    def row_edge(i: int) -> int:
        return int(top) + round(row_h * i)

    buckets: dict[int, list[str]] = {}
    for i, row in enumerate(fracs):
        y0, y1 = row_edge(i), row_edge(i + 1)
        j = 0
        while j < len(col_labels):
            value = row[j] if j < len(row) else None
            b = None if value is None else int(round(_clamp(value, 0.0, 1.0) * _SHADES))
            k = j + 1
            while k < len(col_labels):
                nxt = row[k] if k < len(row) else None
                nb = None if nxt is None else int(round(_clamp(nxt, 0.0, 1.0) * _SHADES))
                if nb != b:
                    break
                k += 1
            if b is not None:
                x0, x1, h = edge(j), edge(k), y1 - y0
                buckets.setdefault(b, []).append(
                    f"M{x0},{y0}h{x1 - x0}v{h}h{x0 - x1}z")
            j = k

    for b in sorted(buckets):
        parts.append(
            f'<path class="shade s{b}" d="{"".join(buckets[b])}" fill="rgb({rgb})"/>'
        )

    if cell_values and cell >= 24:
        for i, row in enumerate(cell_values):
            for j, text in enumerate(row):
                if not text:
                    continue
                parts.append(
                    f'<text x="{left + cell * j + cell / 2:.1f}" '
                    f'y="{top + i * row_h + row_h / 2 + 3.4:.1f}" text-anchor="middle" '
                    f'font-size="{min(11, row_h - 6):.0f}" fill="currentColor">'
                    f"{html.escape(text)}</text>"
                )

    for i, name in enumerate(row_labels):
        y = top + i * row_h + row_h / 2 + font / 3
        parts.append(
            f'<text x="{left - 6:.1f}" y="{y:.1f}" text-anchor="end" '
            f'font-size="{min(font, row_h - 1):.1f}" fill="currentColor">'
            f"{html.escape(name)}</text>"
        )
        if row_notes and i < len(row_notes):
            parts.append(
                f'<text x="{left + grid_w + 8:.1f}" y="{y:.1f}" '
                f'font-size="{min(font, row_h - 1):.1f}" fill="currentColor" '
                f'fill-opacity="0.75">{html.escape(row_notes[i])}</text>'
            )

    height = top + len(row_labels) * row_h + 4
    return _svg("".join(parts), width, height, label)


def _ramp(rgb: str, low_label: str, high_label: str) -> str:
    """A horizontal shade key, so a stretched scale states its own bounds."""
    gid = _gid("ramp")
    # Same classes as the cells they key, so the legend cannot drift from the
    # matrix when the ramp changes.
    stops = "".join(
        f'<stop class="shade s{i}" offset="{i / _SHADES:.3f}" '
        f'stop-color="rgb({rgb})"/>' for i in range(_SHADES + 1)
    )
    body = (f'<defs><linearGradient id="{gid}" x1="0" y1="0" x2="1" y2="0">{stops}'
            "</linearGradient></defs>"
            f'<rect x="0" y="4" width="150" height="9" fill="url(#{gid})"/>'
            f'<text x="0" y="26" font-size="10" fill="currentColor" '
            f'fill-opacity="0.75">{html.escape(low_label)}</text>'
            f'<text x="150" y="26" font-size="10" text-anchor="end" '
            f'fill="currentColor" fill-opacity="0.75">{html.escape(high_label)}</text>')
    return (f'<svg viewBox="0 0 200 30" width="200" height="30" role="img" '
            f'aria-label="shading key">{body}</svg>')


def _matrix_block(names: list[str], matrix: list[list[float | None]], rgb: str, *,
                  fmt, frac, low_label: str, high_label: str, label: str) -> str:
    """A numeric table for a small set, the same data as a heatmap for a large."""
    if len(names) <= _NUMERIC_MATRIX_MAXIMUM:
        return _numeric_matrix(names, matrix, rgb, fmt, frac)
    fracs = [[None if v is None else frac(v) for v in row] for row in matrix]
    return (_grid(names, names, fracs, rgb, label=label)
            + _ramp(rgb, low_label, high_label))


def _section_skani(tool: Tool, ctx: Context, workdir: Path) -> str:
    """ANI and aligned fraction, as two shaded matrices.

    Both, not just ANI: skani reports an identity once alignment covers as
    little as ~15% of a genome, so a high ANI between genomes of very different
    size or completeness can reflect a shared plasmid or conserved core rather
    than whole-genome relatedness. `triangle --full-matrix` writes the aligned
    fraction to `<output>.af` and v3 ignored it until 2026-09-02.
    """
    path = ctx.out("skani", "ani.tsv")
    if not path.exists():
        return '<p class="missing">No results.</p>'
    names, matrix = _reorder(*_read_skani_matrix(path, ctx.samples), ctx.samples)
    if not names:
        return '<p class="missing">Not enough genomes to compare.</p>'

    def floor_of(m: list[list[float | None]]) -> float:
        values = [v for row in m for v in row if v is not None and v < 100]
        return min(values) if values else 95.0

    low = floor_of(matrix)
    parts = [
        f'<p class="summary">Average nucleotide identity, shaded from '
        f"{low:.2f}% to 100%.</p>",
        _matrix_block(
            names, matrix, "43,108,176",
            fmt=lambda v: f"{v:.2f}",
            frac=lambda v: 0.0 if low >= 100 else (v - low) / (100 - low),
            low_label=f"{low:.2f}%", high_label="100%",
            label="average nucleotide identity"),
    ]

    af_path = ctx.out("skani", "ani.tsv.af")
    if af_path.exists():
        af_names, af_matrix = _reorder(
            *_read_skani_matrix(af_path, ctx.samples), ctx.samples)
        if af_names:
            af_low = floor_of(af_matrix)
            parts += [
                "<h3>Aligned fraction</h3>",
                '<p class="summary">How much of the genome in each <em>row</em> '
                "aligned to the genome in each column — so unlike the matrix "
                "above, this one is not symmetric. Read it together with the "
                "ANI: a high identity over a small aligned fraction is not "
                f"whole-genome relatedness. Shaded from {af_low:.2f}% to 100%.</p>",
                _matrix_block(
                    af_names, af_matrix, "183,121,31",
                    fmt=lambda v: f"{v:.2f}",
                    frac=lambda v, lo=af_low: 0.0 if lo >= 100 else (v - lo) / (100 - lo),
                    low_label=f"{af_low:.2f}%", high_label="100%",
                    label="aligned fraction"),
            ]

    return "".join(parts)


def _section_snp_dists(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Pairwise SNP counts, shaded like the ANI matrix rather than dumped raw.

    Two distance matrices in one document rendered two different ways is a
    reader's problem, not an implementation detail: this one used to arrive
    through the generic fallback, unshaded and headed `snp-dists 1.2.0`.
    """
    path = ctx.out("snp-dists", "snp-dists.tsv")
    if not path.exists():
        return '<p class="missing">No results.</p>'
    names, matrix = _reorder(
        *_read_labelled_matrix(path, ctx.samples), ctx.samples)
    if not names:
        return '<p class="missing">Not enough genomes to compare.</p>'

    values = [v for row in matrix for v in row if v is not None]
    hi = max(values) if values else 1.0
    return (
        f'<p class="summary">SNP differences across the core-genome alignment, '
        f"shaded darkest at 0 and lightest at {hi:,.0f}. The scale is stretched "
        "to this set, so read the numbers rather than the shades when comparing "
        "with another run.</p>"
        + _matrix_block(
            names, matrix, "43,108,176",
            fmt=lambda v: f"{v:,.0f}",
            frac=lambda v: 1.0 - (v / hi if hi else 0.0),
            low_label=f"{hi:,.0f} SNPs", high_label="0 SNPs",
            label="pairwise SNP distances")
    )


def _section_treecluster(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Cluster sizes and membership.

    The raw two-column file arrived through the fallback as a table headed
    `SequenceName` / `ClusterNumber` that restated, in words, what the tree
    directly above it already showed in colour. What the tree cannot show at a
    hundred genomes is how many clusters there are and how big they get.
    """
    path = ctx.out("treecluster", "treecluster.tsv")
    if not path.exists():
        return '<p class="missing">No results.</p>'
    assignments = _clusters(workdir)
    if not assignments:
        return '<p class="missing">No assignments.</p>'

    members: dict[str, list[str]] = {}
    for sample in ctx.samples:  # input order inside each cluster, too
        cluster = assignments.get(sample)
        if cluster is not None:
            members.setdefault(cluster, []).append(sample)
    for name, cluster in assignments.items():
        if name not in ctx.samples:
            members.setdefault(cluster, []).append(name)

    # TreeCluster writes -1 for a genome it placed in no cluster. That is a
    # result, not a cluster called "-1", and it sorts last.
    def key(item: tuple[str, list[str]]) -> tuple[int, int, str]:
        return (1 if item[0] in ("-1", "") else 0, -len(item[1]), item[0])

    rows = []
    for cluster, names in sorted(members.items(), key=key):
        shown = ", ".join(names[:8])
        if len(names) > 8:
            shown += f", and {len(names) - 8} more"
        label = "unclustered" if cluster in ("-1", "") else cluster
        rows.append([label, f"{len(names)}", shown])

    real = [c for c in members if c not in ("-1", "")]
    singletons = sum(1 for c in real if len(members[c]) == 1)
    summary = (f"{len(real)} cluster{'s' if len(real) != 1 else ''} across "
               f"{sum(len(v) for v in members.values())} genomes")
    if singletons:
        summary += f", {singletons} of them a single genome"
    if len(real) == 1 and not singletons:
        summary = ("Every genome fell into one cluster at this threshold, so the "
                   "cut separates nothing in this set")
    return (f'<p class="summary">{summary}. Leaves in the trees above are '
            "coloured by these assignments.</p>"
            + _table(rows, header=["Cluster", "Genomes", "Members"]))


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

    found: dict[str, list[str]] = {}
    for rec in rows[1:]:
        if len(rec) <= max(i_name, i_comp, i_cont):
            continue
        comp = float(rec[i_comp]) if _numeric(rec[i_comp]) else 0.0
        cont = float(rec[i_cont]) if _numeric(rec[i_cont]) else 0.0
        # MIMAG-style reading: high quality is >90% complete, <5% contaminated.
        comp_colour = "#2f855a" if comp >= 90 else ("#b7791f" if comp >= 50 else "#9b2c2c")
        cont_colour = "#2f855a" if cont < 5 else ("#b7791f" if cont < 10 else "#9b2c2c")
        name = _sample_of(rec[i_name], ctx.samples)
        found[name] = [
            name,
            f"{_bar(comp, comp_colour)}{comp:.1f}%",
            f"{_bar(min(cont * 5, 100), cont_colour)}{cont:.2f}%",
        ]

    rows_out, absent = _in_sample_order(found, ctx.samples, 3)
    body = _raw_table(
        [[html.escape(name), *cells] for name, *cells in rows_out],
        header=["Genome", "Completeness", "Contamination"], numeric_columns={1, 2})
    return ('<p class="summary">Green is high quality by MIMAG conventions: '
            "&gt;90% complete, &lt;5% contaminated. Contamination bars are "
            'scaled &times;5 so small values stay visible.</p>'
            + body + _absent_note(absent, len(ctx.samples)))


# GTDB's seven ranks, in order, keyed by the one-letter prefix its
# classification strings use: `d__Bacteria;p__Bacillota;…;s__Enterococcus faecium`.
_GTDB_RANKS = (("d", "Domain"), ("p", "Phylum"), ("c", "Class"), ("o", "Order"),
               ("f", "Family"), ("g", "Genus"), ("s", "Species"))


def _gtdb_ranks(classification: str) -> list[str] | None:
    """Split a classification into seven values, blanks included.

    `None` for anything carrying no `x__` prefix at all: GTDB-Tk writes a bare
    `Unclassified Bacteria` in some failure modes, and that is a value to
    display, not a lineage to parse.

    A rank that is present but *empty* stays empty on purpose. A bare `s__`
    means GTDB-Tk would not commit to a species, which guidance.py explains is
    a result rather than a failure — so it has to survive parsing to be shown.
    """
    tags = [tag for tag, _ in _GTDB_RANKS]
    values = [""] * len(_GTDB_RANKS)
    found = False
    for part in classification.split(";"):
        tag, sep, value = part.strip().partition("__")
        if not sep or tag not in tags:
            continue
        values[tags.index(tag)] = value.strip()
        found = True
    return values if found else None


def _gtdb_shared_depth(lineages: list[list[str]]) -> int:
    """How many leading ranks every genome agrees on, ignoring blanks.

    A set of one species shares all seven; a set spanning two phyla shares one.
    The shared part is stated once above the table so the columns can show only
    where the genomes actually differ — the same reason the pangenome figure
    collapses identical presence patterns.
    """
    depth = 0
    for rank in range(len(_GTDB_RANKS)):
        values = {lineage[rank] for lineage in lineages}
        if len(values) != 1 or not lineages[0][rank]:
            break
        depth += 1
    return depth


def _gtdb_ani(ani: str, radius: str) -> str:
    """ANI to the closest reference, judged against *that reference's* radius.

    Which is the point guidance.py makes: the species cut-off is per-reference,
    not a global 95%, so the comparison worth drawing is against the radius the
    tool reported next to it.
    """
    if not _numeric(ani):
        return "—"
    value = float(ani)
    if not _numeric(radius):
        return f"{value:.2f}%"
    limit = float(radius)
    colour = "#2f855a" if value >= limit else "#b7791f"
    return (f'<span style="color:{colour}" title="species radius {limit:.2f}%">'
            f"{value:.2f}%</span>")


def _section_gtdbtk(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Where each genome sits in GTDB, and how far down the name is trustworthy.

    Reads the merged summary that `Tool.post` writes from GTDB-Tk's separate
    `bac120` and `ar53` files, so archaea and bacteria arrive in one table.
    """
    path = ctx.out("gtdbtk", "gtdbtk.summary.tsv")
    if not path.exists():
        return '<p class="missing">No results.</p>'
    rows = _read_tsv(path)
    if len(rows) < 2:
        return '<p class="missing">No genomes classified.</p>'

    at = {name.strip(): i for i, name in enumerate(rows[0])}
    if "user_genome" not in at or "classification" not in at:
        # Not the summary this was written against. Show the file rather than
        # guess at it: the columns are documented, but this one has never been
        # seen here and GTDB-Tk has changed its output between major versions.
        return _table(rows, clip=True)

    def value(rec: list[str], column: str) -> str:
        i = at.get(column)
        return rec[i].strip() if i is not None and len(rec) > i else ""

    found: dict[str, dict[str, str]] = {}
    lineages: dict[str, list[str] | None] = {}
    for rec in rows[1:]:
        name = _sample_of(value(rec, "user_genome"), ctx.samples)
        lineages[name] = _gtdb_ranks(value(rec, "classification"))
        found[name] = {
            "raw": value(rec, "classification"),
            "ani": _gtdb_ani(value(rec, "closest_genome_ani"),
                             value(rec, "closest_genome_reference_radius")),
            "af": value(rec, "closest_genome_af") or "—",
            "method": value(rec, "classification_method") or "—",
        }

    parsed = [lineage for lineage in lineages.values() if lineage is not None]
    shared = _gtdb_shared_depth(parsed) if parsed else 0
    # Never show no ranks at all: a set of one species agrees on all seven, and
    # the species is the column a reader came for.
    start = min(shared, len(_GTDB_RANKS) - 1)
    shown = list(range(start, len(_GTDB_RANKS)))

    cells: dict[str, list[str]] = {}
    for name, lineage in lineages.items():
        if lineage is None:
            # Unparseable, so nothing is dropped: the raw value goes in the
            # first rank column and the rest are marked absent. Marked as
            # missing rather than set like the others, because that column has
            # a rank for its heading and GTDB-Tk's `Unclassified Bacteria` is a
            # domain-level failure — under a `Species` heading, in the same
            # type as `Enterococcus faecium` above it, it reads as a species.
            raw = html.escape(found[name]["raw"]) or "—"
            ranks = [f'<span class="missing" title="GTDB-Tk returned no '
                     f'parseable lineage for this genome">{raw}</span>']
            ranks += ["—"] * (len(shown) - 1)
        else:
            ranks = [html.escape(lineage[i]) or "—" for i in shown]
        cells[name] = [html.escape(name), *ranks, found[name]["ani"],
                       html.escape(found[name]["af"]), html.escape(found[name]["method"])]

    columns = 1 + len(shown) + 3
    rows_out, absent = _in_sample_order(cells, ctx.samples, columns)
    body = _raw_table(
        rows_out,
        header=["Genome", *[label for _, label in
                            (_GTDB_RANKS[i] for i in shown)], "ANI", "AF", "Method"],
        numeric_columns={columns - 3, columns - 2},
    )

    notes = []
    if start:
        # `start`, not `shared`: when every genome is one species all seven
        # ranks are shared, and naming the species here as well as in the
        # column below it says the same thing twice.
        lineage = " › ".join(parsed[0][:start])
        word = "genome" if len(parsed) == 1 else "genomes"
        notes.append(
            f'<p class="summary">All {len(parsed)} classified {word} share '
            f"<strong>{html.escape(lineage)}</strong>. The columns below start "
            "where they diverge.</p>")
    notes.append(
        '<p class="summary">ANI is identity to the closest reference, green '
        "when it clears <em>that reference's own</em> species radius rather "
        "than a global 95% — hover for the radius. Read it with AF: a high "
        "identity over a small aligned fraction is not evidence of one "
        "species.</p>")

    # Two different outcomes, counted separately. A blank species is a genome
    # GTDB-Tk placed but would not name; an unparseable classification is one
    # it did not place. Folding them into one count made the note account for
    # half the empty-looking cells in the column.
    tail = []
    unnamed = sum(1 for lineage in parsed if not lineage[-1])
    if unnamed:
        word = "genome carries" if unnamed == 1 else "genomes carry"
        tail.append(f'<p class="note">{unnamed} {word} no species name. GTDB-Tk '
                    "found no reference within the species radius, which says the "
                    "genome is novel relative to this release — not that the "
                    "assignment failed.</p>")
    unplaced = sum(1 for lineage in lineages.values() if lineage is None)
    if unplaced:
        word = "genome" if unplaced == 1 else "genomes"
        tail.append(f'<p class="note">{unplaced} {word} came back with no '
                    "parseable lineage at all, shown in italic where a rank "
                    "would be. That is not the same as a missing species: "
                    "GTDB-Tk placed no lineage for it.</p>")
    return ("".join(notes) + body + "".join(tail)
            + _absent_note(absent, len(ctx.samples)))


def _section_bakta(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Feature counts per genome, read straight from the GFF3."""
    found: dict[str, list[str]] = {}
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
        found[sample] = [
            sample,
            f"{counts.get('CDS', 0):,}",
            f"{counts.get('tRNA', 0):,}",
            f"{counts.get('rRNA', 0):,}",
            f"{sum(v for k, v in counts.items() if k not in ('CDS', 'tRNA', 'rRNA', 'region')):,}",
        ]
    if not found:
        return '<p class="missing">No annotations.</p>'
    rows, absent = _in_sample_order(found, ctx.samples, 5)
    return (_table(rows, header=["Genome", "CDS", "tRNA", "rRNA", "Other features"])
            + _absent_note(absent, len(ctx.samples)))


def _section_amrfinder(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Resistance classes per genome, as a presence grid.

    Per-genome tables would mean one section per sample; across a wide set the
    useful question is which classes appear where. Drawn rather than tabulated
    because the class names are long — `Lincosamide/Macrolide/Streptogramin` is
    one of them — and a column per class overflowed the page at four genomes.
    Rotated, a long name costs vertical space instead.

    A row with no class is not a class. `--plus` also reports virulence genes,
    and AMRFinderPlus writes them out with `Class=NA` — a literal string, not an
    empty field, so accepting any non-empty `Class` gave them a column of their
    own. It sorted between Mercury and Quaternary Ammonium and drew as `Na`, in
    the styling of a real class: in the seven-genome *S. aureus* run of
    2026-09-02 it took 121 of 220 rows, and because shading is scaled to the
    peak, one genome's 22 virulence hits pushed every genuine class into the
    bottom third of the ramp. It also reversed the row totals — TW20 with 29
    real hits printed a lower total than N315 with 27. They are counted and
    named in the summary instead, because `--plus` asks for them deliberately
    and the per-genome TSV is where they are legible.
    """
    per_genome: dict[str, dict[str, int]] = {}
    classes: dict[str, None] = {}
    unclassed = 0
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
                if name.strip().upper() == "NA":
                    unclassed += 1
                    continue
                counts[name] = counts.get(name, 0) + 1
                classes.setdefault(name, None)
        per_genome[sample] = counts

    # Named rather than merely subtracted: a reader who counts the rows of the
    # TSV and the cells of this grid should be told where the difference went.
    note = (f" A further {unclassed} hit{'s' if unclassed != 1 else ''} "
            "carried no resistance class and are not shown — `--plus` also "
            "reports virulence genes, which have none. They are in each "
            "genome's own table." if unclassed else "")

    if not per_genome:
        return '<p class="missing">No resistance genes detected.</p>'
    if not classes:
        return ('<p class="summary">Every genome was searched and none carried a '
                f"gene with a resistance class assigned.{note}</p>")

    order = sorted(classes)
    absent = len(ctx.samples) - len(per_genome)
    peak = max((c for counts in per_genome.values() for c in counts.values()), default=1)

    # Every genome gets a row, including one the tool produced nothing for: an
    # empty row and a — for its total says "not run", where leaving it out says
    # "carries no resistance genes".
    fracs, values, totals = [], [], []
    for sample in ctx.samples:
        counts = per_genome.get(sample)
        if counts is None:
            fracs.append([None] * len(order))
            values.append([""] * len(order))
            totals.append("—")
            continue
        fracs.append([(counts[c] / peak) if c in counts else None for c in order])
        values.append([str(counts.get(c, "")) for c in order])
        totals.append(str(sum(counts.values())))

    return (
        f'<p class="summary">{len(order)} resistance class'
        f"{'es' if len(order) != 1 else ''} across {len(per_genome)} genome"
        f"{'s' if len(per_genome) != 1 else ''}. Shading is the number of hits in "
        f"that class, darkest at {peak}; the figure after each row is the "
        f"genome's total.{note}</p>"
        + _grid(list(ctx.samples), [c.title() for c in order],
                fracs, "43,108,176", row_notes=totals, cell_values=values,
                label="resistance classes per genome")
        + _absent_note(absent, len(ctx.samples))
    )


# A hard ceiling on the elements the pangenome figure may emit, across all
# rows. Sub-pixel gap merging alone is not a bound — it depends on the data
# interleaving favourably, and a set whose patterns alternate on the pixel grid
# produced a megabyte of markup for sixty genomes. Beyond this budget the
# coarsest distinctions survive and the finest are painted over, which is the
# right way round: they were narrower than a pixel to begin with.
_PANGENOME_SUBPATH_BUDGET = 12_000


def _coarsen(runs: list[list[float]], cap: int) -> list[list[float]]:
    """Merge the narrowest gaps until at most `cap` runs remain."""
    if len(runs) <= cap:
        return runs
    gaps = sorted(runs[i + 1][0] - runs[i][1] for i in range(len(runs) - 1))
    threshold = gaps[len(runs) - cap - 1]
    merged = [list(runs[0])]
    for run in runs[1:]:
        if run[0] - merged[-1][1] <= threshold:
            merged[-1][1] = run[1]
        else:
            merged.append(list(run))
    return merged


def draw_pangenome(names: list[str], patterns: list[tuple[tuple[bool, ...], int]],
                   width: int = WIDTH) -> str:
    """The pangenome presence matrix, as inline SVG.

    Genes sharing a presence pattern are drawn as one block whose width is
    proportional to how many genes share it. Blocks run from most to least
    common, so the core sits on the left and the genome-specific genes on the
    right.

    Two things make that survive a hundred genomes. Block widths are exact
    fractions of the span — the old code gave every block a 0.6px minimum, so a
    set with thousands of patterns ran the figure several times past its own
    viewBox and everything after the first screenful was drawn outside it and
    never seen. And a row is one `<path>` whose subpaths are the maximal runs
    where that genome is present, rather than one `<rect>` per block, which at
    100 genomes x 10,000 patterns would be a million elements.
    """
    if not names or not patterns:
        return '<p class="missing">No pangenome.</p>'

    n = len(names)
    row = _clamp(round(1500 / n), 6, 22)
    font = _clamp(row - 4, 7, 12)
    label_room = min(max(_text_width(s, font) for s in names) + 10, width * 0.32)
    span = width - label_room
    total = sum(count for _, count in patterns) or 1
    height = n * row + 26
    cap = max(24, _PANGENOME_SUBPATH_BUDGET // n)

    edges: list[tuple[float, float]] = []
    x = float(label_room)
    for _, cnt in patterns:
        w = span * cnt / total
        edges.append((x, x + w))
        x += w

    parts: list[str] = []
    coarsened = 0
    for i, name in enumerate(names):
        y = i * row
        parts.append(
            f'<text x="{label_room - 8:.1f}" y="{y + row / 2 + font / 3:.1f}" '
            f'text-anchor="end" font-size="{font:.0f}" fill="currentColor">'
            f"{html.escape(name)}</text>"
        )
        runs: list[list[float]] = []
        for (x0, x1), (pattern, _) in zip(edges, patterns):
            if i >= len(pattern) or not pattern[i]:
                continue
            # Merge across sub-pixel gaps. The alternative is dropping blocks
            # too narrow to see, which would hide exactly the genome-specific
            # genes the right-hand end of the figure exists to show.
            if runs and x0 - runs[-1][1] < 0.35:
                runs[-1][1] = x1
            else:
                runs.append([x0, x1])
        if len(runs) > cap:
            runs = _coarsen(runs, cap)
            coarsened += 1
        if not runs:
            continue
        d = []
        for x0, x1 in runs:
            w = max(x1 - x0, 0.5)
            x0 = min(x0, label_room + span - w)
            d.append(f"M{x0:.2f},{y + 2:.1f}h{w:.2f}v{row - 4:.1f}h{-w:.2f}z")
        parts.append(
            f'<path d="{"".join(d)}" fill="#2b6cb0" fill-opacity="0.85"/>')

    parts.append(
        f'<line x1="{label_room:.1f}" y1="{n * row + 2:.1f}" x2="{width}" '
        f'y2="{n * row + 2:.1f}" stroke="currentColor" stroke-opacity="0.25"/>'
    )
    parts.append(
        f'<text x="{label_room:.1f}" y="{n * row + 18:.1f}" font-size="11" '
        f'fill="currentColor" fill-opacity="0.6">core</text>'
    )
    parts.append(
        f'<text x="{width}" y="{n * row + 18:.1f}" text-anchor="end" '
        f'font-size="11" fill="currentColor" fill-opacity="0.6">genome-specific</text>'
    )
    svg = _svg("".join(parts), width, height, "pangenome presence matrix")
    if coarsened:
        svg += ('<p class="note">Blocks narrower than a pixel are painted over '
                f"in {coarsened} of {n} rows, so a gap that thin reads as "
                "present. The counts below are exact.</p>")
    return svg


# Below this many genomes, the conventional pangenome bins are arithmetically
# unreachable and the table lies by omission. Soft core (95–99%) needs
# (N−1)/N >= 0.95, i.e. N >= 20; Cloud (<15%) needs 1/N < 0.15, i.e. N >= 7.
# At 20 all four bins can hold something, so that is the switch-over.
_PARTITION_VOCABULARY_MINIMUM = 20

# Above this many genomes the dot column is wider than the page and its legend
# is a list as long as the genome set, so the overlap table is dropped: the
# partition table answers the same question by frequency instead of by name.
_OVERLAP_VOCABULARY_MAXIMUM = 12


def _partitions(shared: dict[int, int], n: int, total: int) -> str:
    """The core/accessory split — exact counts on small sets, bins on large ones.

    The Core/Soft core/Shell/Cloud vocabulary is what the literature uses, but
    its boundaries are fractions of N, so on a handful of genomes two of the
    four bins cannot contain anything and every accessory cluster piles into
    Shell. Printing a structurally-empty row as though it were a measurement is
    worse than not printing it, so below `_PARTITION_VOCABULARY_MINIMUM` the
    table shows what was actually counted: clusters against the number of
    genomes sharing them.
    """
    if n >= _PARTITION_VOCABULARY_MINIMUM:
        bins = [("Core (≥99%)", 0.99, 1.01), ("Soft core (95–99%)", 0.95, 0.99),
                ("Shell (15–95%)", 0.15, 0.95), ("Cloud (<15%)", 0.0, 0.15)]
        rows = []
        for label, low, high in bins:
            count = sum(c for k, c in shared.items() if low <= k / n < high)
            rows.append([label, f"{count:,}",
                         f"{100 * count / total:.1f}%" if total else "—"])
        rows.append(["Total", f"{total:,}", "100.0%"])
        return _table(rows, header=["Partition", "Gene clusters", "Share"])

    rows = []
    for k in range(n, 0, -1):
        count = shared.get(k, 0)
        if k == n:
            label = f"All {n} genomes (core)"
        elif k == 1:
            label = "1 genome only"
        else:
            label = f"{k} of {n} genomes"
        rows.append([label, f"{count:,}", f"{100 * count / total:.1f}%" if total else "—"])
    rows.append(["Total", f"{total:,}", "100.0%"])
    return _table(rows, header=["Shared by", "Gene clusters", "Share"])


def _section_panaroo(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Pangenome structure: the presence matrix, the overlaps, and the split."""
    path = ctx.out("panaroo", "gene_presence_absence.Rtab")
    if not path.exists():
        return '<p class="missing">No pangenome.</p>'
    rows = _read_tsv(path)
    if len(rows) < 2:
        return '<p class="missing">No genes.</p>'

    file_names = rows[0][1:]
    n = len(file_names)
    if n == 0:
        return '<p class="missing">No genomes in the matrix.</p>'
    # Panaroo emits its columns in its own order; put them back in input order
    # so the figure's rows line up with every other section's.
    order = [file_names.index(s) for s in ctx.samples if s in file_names]
    order += [i for i in range(n) if i not in order]
    names = [file_names[i] for i in order]

    tally: dict[tuple[bool, ...], int] = {}
    per_genome = [0] * n
    private = [0] * n
    # Clusters counted by how many genomes carry them. This is the raw quantity;
    # both partition views below are derived from it, so neither can disagree
    # with the other or with the matrix.
    shared: dict[int, int] = {}
    for rec in rows[1:]:
        raw = rec[1:1 + n]
        if len(raw) != n:
            continue
        pattern = tuple(raw[i].strip() == "1" for i in order)
        tally[pattern] = tally.get(pattern, 0) + 1
        k = sum(pattern)
        for i, present in enumerate(pattern):
            if present:
                per_genome[i] += 1
                if k == 1:
                    private[i] += 1
        shared[k] = shared.get(k, 0) + 1

    total = sum(tally.values())
    if not total:
        return '<p class="missing">No genes.</p>'
    # Most-shared first, so the matrix reads core on the left.
    ordered = sorted(tally.items(), key=lambda kv: (-sum(kv[0]), -kv[1]))

    parts = [
        f'<p class="summary">{total:,} gene clusters across {n} genomes, in '
        f"{len(tally):,} distinct presence patterns. Each block below is one "
        "pattern, its width proportional to the number of genes sharing it.</p>",
        draw_pangenome(names, ordered),
    ]

    if n <= _OVERLAP_VOCABULARY_MAXIMUM:
        # Overlaps. The presence pattern is shown as dots in a fixed genome
        # order rather than a comma-separated name list: the lists grow
        # unreadable past a handful of genomes, and dots line up column-wise so
        # patterns are scannable and match the row order of the figure above.
        by_size = sorted(tally.items(), key=lambda kv: -kv[1])
        top, rest = by_size[:12], by_size[12:]
        overlap_rows = []
        for pattern, cnt in top:
            dots = "".join("●" if p else "·" for p in pattern)
            members = [names[i] for i, p in enumerate(pattern) if p]
            if len(members) == n:
                label = "all genomes"
            elif len(members) == 1:
                label = f"only {html.escape(members[0])}"
            else:
                label = ", ".join(html.escape(m) for m in members)
            share = 100 * cnt / total
            overlap_rows.append([
                f'<span class="dots">{dots}</span>', label, f"{cnt:,}",
                _bar(share, "#2b6cb0", width=60) + f"{share:.1f}%",
            ])
        if rest:
            other = sum(c for _, c in rest)
            overlap_rows.append([
                "", f"{len(rest):,} further patterns", f"{other:,}",
                _bar(100 * other / total, "#2b6cb0", width=60)
                + f"{100 * other / total:.1f}%",
            ])
        legend = " ".join(
            f'<span class="dots">{"·" * i}●{"·" * (n - i - 1)}</span> {html.escape(name)}'
            for i, name in enumerate(names)
        )
        parts += [
            "<h3>Shared gene content</h3>",
            f'<p class="summary">Reading order: {legend}</p>',
            _raw_table(overlap_rows,
                       header=["Pattern", "Present in", "Gene clusters", "Share"],
                       numeric_columns={2, 3}),
        ]

    parts.append("<h3>Pangenome partitions</h3>")
    if n < _PARTITION_VOCABULARY_MINIMUM:
        parts.append(
            f'<p class="summary">With {n} genomes the usual Core/Soft core/Shell/Cloud '
            "bins are fractions that cannot all be reached, so this shows the exact "
            "count instead. See the notes above.</p>")
    parts.append(_partitions(shared, n, total))

    parts += [
        "<h3>Genes per genome</h3>",
        '<p class="summary">Private clusters are those found in this genome and '
        "no other.</p>",
        _table(
            [[name, f"{per_genome[i]:,}", f"{100 * per_genome[i] / total:.1f}%",
              f"{private[i]:,}"]
             for i, name in enumerate(names)],
            header=["Genome", "Gene clusters", "Of pangenome", "Private"],
        ),
    ]
    return "".join(parts)


def _section_carveme(tool: Tool, ctx: Context, workdir: Path) -> str:
    """Model size per genome, counted straight out of the SBML.

    Genome-scale metabolic models are the one capability no comparable pipeline
    offers, so the report should say something about them rather than noting a
    file exists. Counted by streaming the XML — a 4 MB model per genome is not
    worth a parser dependency.
    """
    found: dict[str, list[str]] = {}
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
        found[sample] = [sample, f"{reactions:,}", f"{species:,}", f"{genes:,}"]
    if not found:
        return '<p class="missing">No models.</p>'
    rows, absent = _in_sample_order(found, ctx.samples, 4)
    return ('<p class="summary">Draft models in SBML, ready for flux balance '
            "analysis. Not gap-filled for a specific medium.</p>"
            + _table(rows, header=["Genome", "Reactions", "Metabolites", "Genes"])
            + _absent_note(absent, len(ctx.samples)))


# Dark is a dependency, not an ability: the cells worth looking at are the ones
# a genome has to be given. Reversing this put nine tenths of a prototrophic
# genome's row at full shade and made the two interesting cells invisible.
#
# `de_novo` is one step up from the floor rather than on it, because a cell at
# the floor is nearly indistinguishable from the unshaded blank that means "not
# in this model" — and "makes it itself" against "was never in the network" is
# the last pair of meanings that should be allowed to look alike.
_SHADE_BY_VERDICT = {DE_NOVO: 0.12, UPSTREAM: 0.55, NO_ROUTE: 1.0}
# A glyph as well as a shade, so the grid survives being printed, and a reader
# who cannot separate two ambers still gets the four states apart. `·` for
# absent is the same mark the tables use for a cell a tool produced nothing for.
_GLYPH_BY_VERDICT = {DE_NOVO: "", UPSTREAM: "~", NO_ROUTE: "×", ABSENT: "·"}

_MEDIA_COLUMNS = (("M9", "M9"), ("M9[-O2]", "M9 anaerobic"),
                  ("LB", "LB"), ("LB[-O2]", "LB anaerobic"),
                  ("complete", "Complete"))

# How many compound names are worth spelling out in a note. Past this the
# sentence is longer than the grid it describes.
_BIOSYNTHESIS_NAME_BUDGET = 12


def _biosynthesis_naming(varying: list) -> str:
    """Name whichever of the two sets is short enough to be useful.

    Which columns to look at is the point of the sentence. At one species most
    of the panel varies and the short list is the *uniform* one; across phyla it
    is the other way round; at a hundred genomes neither is worth printing.
    """
    if not varying:
        return "."
    if len(varying) <= _BIOSYNTHESIS_NAME_BUDGET:
        return ": " + ", ".join(c.name for c in varying) + "."
    same = [c for c in PANEL if c not in varying]
    if len(same) <= _BIOSYNTHESIS_NAME_BUDGET:
        return ". The same in every genome: " + ", ".join(c.name for c in same) + "."
    return "."


def _section_biosynthesis(tool: Tool, ctx: Context, workdir: Path) -> str:
    """What each genome can build, and what it has to be given.

    Three states per compound rather than two, and a grid rather than a table:
    32 compounds is more columns than a table holds, and the comparison a
    reader wants is down a column — which genome differs from the others.
    """
    verdicts: dict[str, dict[str, str]] = {}
    media: dict[str, dict[str, list[str]]] = {}
    for sample in ctx.samples:
        panel = ctx.sample_out(sample, "biosynthesis", f"{sample}.tsv")
        if panel.exists():
            verdicts[sample] = {r[0]: r[3] for r in _read_tsv(panel)[1:] if len(r) > 3}
        table = ctx.sample_out(sample, "biosynthesis", f"{sample}.media.tsv")
        if table.exists():
            media[sample] = {r[0]: r for r in _read_tsv(table)[1:] if len(r) > 3}
    if not verdicts:
        return '<p class="missing">No capability tables.</p>'

    counts: dict[str, list[str]] = {}
    for sample, found in verdicts.items():
        tally = {k: 0 for k in (DE_NOVO, UPSTREAM, NO_ROUTE, ABSENT)}
        for verdict in found.values():
            if verdict in tally:
                tally[verdict] += 1
        counts[sample] = [sample, str(tally[DE_NOVO]), str(tally[UPSTREAM]),
                          str(tally[NO_ROUTE]), str(tally[ABSENT])]
    rows, absent = _in_sample_order(counts, ctx.samples, 5)

    # Only the genomes that produced a table get a grid row, and the summary
    # above already says which did not — a blank grid row would read as a
    # genome that makes everything.
    named = [s for s in ctx.samples if s in verdicts]
    fracs, cells, notes = [], [], []
    for sample in named:
        found = verdicts[sample]
        fracs.append([_SHADE_BY_VERDICT.get(found.get(c.bigg, ABSENT)) for c in PANEL])
        cells.append([_GLYPH_BY_VERDICT.get(found.get(c.bigg, ABSENT), "") for c in PANEL])
        notes.append(str(sum(1 for c in PANEL if found.get(c.bigg) == NO_ROUTE)))

    varying = [c for c in PANEL
               if len({verdicts[s].get(c.bigg, ABSENT) for s in named}) > 1]

    parts = [
        f'<p class="summary">{len(PANEL)} building blocks per genome. '
        "<em>De novo</em> is a complete route from a minimal medium; "
        "<em>upstream</em> means the route is there but another compound on the "
        "list is what is missing; <em>no route</em> is what the genome has to be "
        "given. These describe the draft model, not the organism — see the notes "
        "above.</p>",
        _table(rows, header=["Genome", "De novo", "Upstream", "No route",
                             "Not in model"]),
        _absent_note(absent, len(ctx.samples)),
        "<h3>Building blocks per genome</h3>",
        '<p class="summary">Darkest is <em>no route</em> (×), mid is '
        "<em>upstream</em> (~), lightest is <em>de novo</em>, and an unshaded "
        "cell (·) is a compound this model does not contain. The glyphs appear "
        "where the grid is wide enough for them. The figure after each row is "
        "the genome's <em>no route</em> count.</p>",
        _grid(named, [c.name for c in PANEL], fracs, "183,121,31",
              row_notes=notes, cell_values=cells,
              label="building blocks each genome can make"),
    ]
    if len(named) > 1:
        parts.append(
            f'<p class="note">{len(varying)} of {len(PANEL)} compounds differ '
            f"between these {len(named)} genomes"
            + _biosynthesis_naming(varying)
            + " A column that is the same everywhere carries no comparative "
            "information, and can be a real lineage character as easily as a gap "
            "in the reaction database.</p>")

    if media:
        # Parsed defensively: the report is meant to survive a partial run, and
        # these files can be read while a rule is still writing them.
        def cell(sample: str, medium: str, column: int) -> float | None:
            row = media.get(sample, {}).get(medium)
            if row is None or len(row) <= column:
                return None
            try:
                return float(row[column])
            except ValueError:
                return None

        grew = [s for s in named
                if any((cell(s, m, 3) or 0.0) > 0 for m, _ in _MEDIA_COLUMNS[:4])]
        present = [int(v) for v in (cell(s, "LB", 2) for s in media) if v is not None]
        mrows = {}
        for sample in media:
            mrows[sample] = [sample] + [media[sample].get(key, ["", "", "", "—"])[3]
                                        for key, _ in _MEDIA_COLUMNS]
        ordered, missing = _in_sample_order(mrows, ctx.samples, 6)
        parts += [
            "<h3>Growth on the reference media</h3>",
            '<p class="summary">The feasibility check to make before running flux '
            "balance analysis on these models, in h⁻¹, with every compound capped "
            "at 10 mmol/gDW/h. Not a growth rate to quote — the cap applies to "
            f"oxygen too. {len(grew)} of {len(named)} genomes grow on any defined "
            "medium here; the complete medium is the control that the model is "
            "feasible at all.</p>",
            _table(ordered, header=["Genome", *[label for _, label in _MEDIA_COLUMNS]]),
            _absent_note(missing, len(ctx.samples)),
        ]
        if present:
            lo, hi = min(present), max(present)
            span = f"{lo}" if lo == hi else f"{lo}–{hi}"
            parts.append(
                '<p class="note">A zero on a rich medium is usually missing '
                "transport, not missing metabolism: these models carry exchange "
                f"reactions for {span} of LB's {len(LB)} compounds.</p>")

    return "".join(parts)


SECTIONS = {
    "carveme": _section_carveme,
    "biosynthesis": _section_biosynthesis,
    "seqkit": _section_seqkit,
    "checkm2": _section_checkm2,
    "gtdbtk": _section_gtdbtk,
    "bakta": _section_bakta,
    "amrfinder": _section_amrfinder,
    "mashtree": _section_tree,
    "treecluster": _section_treecluster,
    "mlst": _section_mlst,
    "skani": _section_skani,
    "snp-dists": _section_snp_dists,
    "panaroo": _section_panaroo,
    "fasttree": _section_tree,  # same phylogram renderer, its own newick
}


def _about(tool: Tool) -> str:
    """The explanatory block for one tool, or nothing if it has no guidance.

    Every string here is tool-independent prose from `guidance.py`, but it is
    escaped anyway: the alternative is a module where some fields may contain
    markup and some may not, which is the kind of distinction that stops being
    true the first time someone adds an entry.
    """
    entry = GUIDANCE.get(tool.name)
    if entry is None:
        return ""

    parts = [
        f"<p>{html.escape(entry.blurb)}</p>",
        f'<p class="method">{html.escape(entry.method)}</p>',
    ]

    if entry.reading:
        items = "".join(
            f"<dt>{html.escape(label)}</dt><dd>{html.escape(text)}</dd>"
            for label, text in entry.reading
        )
        parts.append(f"<h4>Reading this section</h4><dl>{items}</dl>")

    if entry.caveats:
        items = "".join(f"<li>{html.escape(c)}</li>" for c in entry.caveats)
        parts.append(f"<h4>What this cannot tell you</h4><ul>{items}</ul>")

    papers = [entry.citation, *entry.also]
    cites = "; ".join(
        f'{html.escape(c.text)} <a href="{html.escape(c.url)}">doi:{html.escape(c.doi)}</a>'
        for c in papers
    )
    parts.append(f'<p class="cite">{cites}</p>')

    return ('<details class="about"><summary>What this is, and how to read it'
            f'</summary><div>{"".join(parts)}</div></details>')


def _methods(names: list[str]) -> str:
    """Every paper behind the tools that actually ran.

    At the end rather than the top, and flat rather than per-tool, because its
    job is to be pasted into a manuscript's methods section — which is also why
    a tool that produced no output must not appear in it.
    """
    papers = citations(names)
    if not papers:
        return ""
    items = []
    for c in papers:
        note = f' <span class="aside">{html.escape(c.note)}</span>' if c.note else ""
        items.append(
            f'<li>{html.escape(c.text)}. '
            f'<a href="{html.escape(c.url)}">doi:{html.escape(c.doi)}</a>{note}</li>'
        )
    return (
        '<h2 id="methods">Methods and citations</h2>'
        '<p class="summary">Every tool that produced output above, and the '
        "underlying methods it runs. Please cite these alongside CompareM2.</p>"
        f'<ol class="refs">{"".join(items)}</ol>'
    )


def _fallback(tool: Tool, ctx: Context, workdir: Path, samples: int = 0) -> str:
    """A plain table of the tool's first output. Every tool gets at least this."""
    outputs = list(tool.outputs(ctx))
    # One row per genome plus a header is the floor: a 50-row cap silently
    # truncated any per-genome fallback on a set bigger than that.
    cap = max(51, samples + 1)
    for path in outputs:
        if path.exists() and path.suffix in {".tsv", ".txt", ".Rtab"}:
            rows = _read_tsv(path, limit=cap + 1)
            notes = []
            if len(rows) > cap:
                rows = rows[:cap]
                notes.append(f"First {cap - 1:,} rows.")
            if any(len(c) > _UNDESIGNED_CELL_CHARS for r in rows for c in r):
                notes.append("Long values are clipped — hover a cell for the "
                             f"whole text, or read {html.escape(path.name)} "
                             "for the table itself.")
            note = f'<p class="note">{" ".join(notes)}</p>' if notes else ""
            return _table(rows, clip=True) + note
    shown = ", ".join(html.escape(p.name) for p in outputs)
    return f'<p class="missing">Produced {shown}; no table renderer yet.</p>'


def render_report(registry: Registry, selected: list[str] | None, workdir: Path,
                  databases: Path, samples: tuple[str, ...],
                  title: str | None = None, command: str | None = None) -> Path:
    """Write the report and return its path."""
    sections: list[str] = []
    toc: list[str] = []
    shown = 0
    ran: list[str] = []

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
        ran.append(tool.name)
        renderer = SECTIONS.get(tool.name)
        body = (renderer(tool, ctx, workdir) if renderer
                else _fallback(tool, ctx, workdir, len(samples)))
        name = html.escape(tool.name)
        toc.append(f'<a href="#{name}">{name}</a>')
        sections += [
            f'<h2 id="{name}">{name}</h2>',
            f'<p class="summary">{html.escape(tool.summary)}</p>',
            _about(tool),
            body,
        ]

    if not shown:
        sections.append('<p class="missing">No results yet.</p>')
    methods = _methods(ran)
    if methods:
        toc.append('<a href="#methods">methods</a>')

    heading = html.escape(title or workdir.name)
    total = len(registry.closure(selected))
    # Provenance at the top, not the bottom: this file is meant to be copied off
    # a cluster and read somewhere else, where when it was made and by what are
    # the first two questions.
    meta = [f"CompareM2 {__version__}",
            datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")]
    if command:
        meta.append(f"<code>{html.escape(command)}</code>")

    parts = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{heading}</title><style>{CSS}</style></head><body>",
        f"<h1>{heading}</h1>",
        f'<p class="sub">{len(samples)} assemblies &middot; '
        f"{shown} of {total} tools produced output</p>",
        f'<p class="meta">{" &middot; ".join(meta)}</p>',
    ]
    if toc:
        parts.append(
            '<nav class="toc">'
            + '<span class="sep">/</span>'.join(toc)
            + "</nav>"
        )
    parts += sections
    parts.append(methods)
    parts.append(f"<footer>CompareM2 {__version__} — {heading}</footer></body></html>")

    path = workdir / "report.html"
    path.write_text("".join(parts))
    return path
