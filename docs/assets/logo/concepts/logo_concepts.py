#!/usr/bin/env python3
"""CompareM2 logo concepts — "compared contigs".

Every concept is a function returning SVG elements inside a 64x64 box, drawn
in two colours: `c1` for the thing the comparison finds (shared / homologous /
the operator) and `c2` for everything else. Colours match report.py:
--accent #2b6cb0 (light) / #7fb3e8 (dark).

Usage: python logo_concepts.py OUTDIR
Writes per-concept SVG+PNG (light, dark, badge, lockup) and a preview sheet.
"""
import math
import os
import sys
from pathlib import Path

import cairosvg
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.ttLib import TTFont

# Inter (OFL), same location make_logo.py expects; override with INTER_FONT=...
FONT = Path(os.environ.get("INTER_FONT", "/tmp/inter/extras/otf/Inter-SemiBold.otf"))
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "out")
OUT.mkdir(parents=True, exist_ok=True)

ACCENT_L, TINT_L = "#2b6cb0", "#7fb3e8"      # light theme: accent + tint
ACCENT_D, TINT_D = "#7fb3e8", "#3f80c4"      # dark theme
FG_L, FG_D = "#1a1a1a", "#e8e8e8"
BG_D = "#161616"
BADGE_BG, BADGE_C1, BADGE_C2 = "#2b6cb0", "#ffffff", "#bcd8f5"

H, RX = 10, 5  # contig pill height / radius


def pill(x, w, y, fill, h=H, rx=RX):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}"/>'


# ---------------------------------------------------------------- concepts --
def concept_existing(c1, c2, uid=""):
    """The Sept-7 draft: three rows of contigs, alternating tints."""
    rows = [(11, [(8, 27), (40, 16)]), (27, [(8, 13), (26, 30)]), (43, [(8, 36)])]
    out = []
    for i, (y, segs) in enumerate(rows):
        col = c1 if i % 2 == 0 else c2
        out += [pill(x, w, y, col) for x, w in segs]
    return "\n".join(out)


def concept_column(c1, c2, uid="col"):
    """A. Conserved column: three genomes broken at different places; one
    region (x 24-40) is present in all of them and lights up in c1."""
    rows = [(10, [(8, 36), (48, 8)]), (27, [(8, 10), (22, 34)]), (44, [(8, 48)])]
    bx, bw = 24, 16
    out = ["<defs>"]
    for i, (y, segs) in enumerate(rows):
        for j, (x, w) in enumerate(segs):
            out.append(f'<clipPath id="{uid}{i}{j}"><rect x="{x}" y="{y}" width="{w}" height="{H}" rx="{RX}"/></clipPath>')
    out.append("</defs>")
    for i, (y, segs) in enumerate(rows):
        for j, (x, w) in enumerate(segs):
            out.append(pill(x, w, y, c2))
            if x <= bx and x + w >= bx + bw:  # the shared block lives in this contig
                out.append(f'<rect x="{bx}" y="{y}" width="{bw}" height="{H}" fill="{c1}" clip-path="url(#{uid}{i}{j})"/>')
    return "\n".join(out)


def concept_ribbons(c1, c2, uid="rib"):
    """B. Alignment ribbons: two genomes, two homologous blocks that have
    swapped places; translucent bands connect each block to its match."""
    y1, y2 = 8, 46
    top = [(8, 24, c1), (36, 20, c2)]       # P, Q
    bot = [(8, 20, c2), (32, 24, c1)]       # Q', P'
    out = []
    # ribbons first (under the pills)
    def ribbon(ax, aw, bx, bw, col, inset=4):
        ax, aw, bx, bw = ax + inset, aw - 2 * inset, bx + inset, bw - 2 * inset
        ya, yb = y1 + H - 1, y2 + 1
        m = (ya + yb) / 2
        return (f'<path fill="{col}" fill-opacity="0.30" d="M{ax},{ya} L{ax+aw},{ya} '
                f'C{ax+aw},{m} {bx+bw},{m} {bx+bw},{yb} L{bx},{yb} C{bx},{m} {ax},{m} {ax},{ya} Z"/>')
    out.append(ribbon(8, 24, 32, 24, c1))
    out.append(ribbon(36, 20, 8, 20, c2))
    out += [pill(x, w, y1, col) for x, w, col in top]
    out += [pill(x, w, y2, col) for x, w, col in bot]
    return "\n".join(out)


def concept_rings(c1, c2, uid="ring"):
    """C. Concentric rings (BRIG-style): two circular genomes drawn as rings
    of contig arcs, outer in c1, inner in c2."""
    cx, cy = 32, 32
    sw = 8
    def arc(r, a0, a1, col):
        a0r, a1r = math.radians(a0), math.radians(a1)
        x0, y0 = cx + r * math.cos(a0r), cy + r * math.sin(a0r)
        x1, y1 = cx + r * math.cos(a1r), cy + r * math.sin(a1r)
        large = 1 if (a1 - a0) > 180 else 0
        return (f'<path d="M{x0:.2f},{y0:.2f} A{r},{r} 0 {large} 1 {x1:.2f},{y1:.2f}" '
                f'fill="none" stroke="{col}" stroke-width="{sw}" stroke-linecap="round"/>')
    out = []
    # outer ring r=24: three contigs. gap ≈ 27° leaves ~4px daylight after round caps
    ro, gi = 24, 27
    outer = [(-90, 20), (20 + gi, 150), (150 + gi, 270 - gi)]
    for a0, a1 in outer:
        out.append(arc(ro, a0, a1, c1))
    # inner ring r=12: two contigs, gap ≈ 55°
    ri, gj = 12, 55
    inner = [(-30, 110), (110 + gj, 330 - gj)]
    for a0, a1 in inner:
        out.append(arc(ri, a0, a1, c2))
    return "\n".join(out)


def concept_equals(c1, c2, uid="eq"):
    """D. The comparison operator: two genomes stacked so the rows read as
    '=', broken into contigs at different places."""
    h = 12
    out = [pill(8, 28, 16, c1, h=h, rx=6), pill(41, 15, 16, c1, h=h, rx=6),
           pill(8, 15, 36, c2, h=h, rx=6), pill(28, 28, 36, c2, h=h, rx=6)]
    return "\n".join(out)


CONCEPTS = [
    ("0 existing draft", "existing", concept_existing),
    ("A conserved column", "column", concept_column),
    ("B alignment ribbons", "ribbons", concept_ribbons),
    ("C concentric rings", "rings", concept_rings),
    ("D equals mark", "equals", concept_equals),
]


# ---------------------------------------------------------------- wordmark --
font = TTFont(str(FONT))
glyphset = font.getGlyphSet()
cmap = font.getBestCmap()
UPM = font["head"].unitsPerEm
hmtx = font["hmtx"]
TRACK = -0.01 * UPM


def text_paths(text, size):
    """[(path_d, x_offset_user_units)], total width in user units."""
    s = size / UPM
    x = 0.0
    parts = []
    for ch in text:
        g = cmap[ord(ch)]
        pen = SVGPathPen(glyphset)
        glyphset[g].draw(pen)
        parts.append((pen.getCommands(), x * s))
        x += hmtx[g][0] + TRACK
    return parts, x * s, s


def wordmark(fg, accent, size, x0, baseline):
    svg = []
    p1, w1, s = text_paths("Compare", size)
    p2, w2, _ = text_paths("M2", size)
    for d, dx in p1:
        svg.append(f'<path fill="{fg}" transform="translate({x0+dx:.2f},{baseline}) scale({s:.5f},{-s:.5f})" d="{d}"/>')
    for d, dx in p2:
        svg.append(f'<path fill="{accent}" transform="translate({x0+w1+dx:.2f},{baseline}) scale({s:.5f},{-s:.5f})" d="{d}"/>')
    return "\n".join(svg), w1 + w2


# ---------------------------------------------------------------- wrappers --
def wrap(inner, w, h, vb=None):
    vb = vb or f"0 0 {w} {h}"
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb}" width="{w}" height="{h}">\n{inner}\n</svg>\n'


def icon_svg(fn, c1, c2, bg=None, uid=""):
    bg_el = f'<rect width="64" height="64" rx="14" fill="{bg}"/>\n' if bg else ""
    return wrap(bg_el + fn(c1, c2, uid=uid), 64, 64)


def lockup_svg(fn, fg, accent, c1, c2, bg=None, uid=""):
    size, icon, pad, gap = 44, 40, 12, 12
    wm, wm_w = wordmark(fg, accent, size, pad + icon + gap, pad + 33)
    W, Hh = pad + icon + gap + wm_w + pad, icon + 2 * pad
    bg_el = f'<rect width="{W:.0f}" height="{Hh}" fill="{bg}"/>\n' if bg else ""
    body = (f'{bg_el}<g transform="translate({pad},{pad+2}) scale({icon/64:.4f})">\n{fn(c1, c2, uid=uid)}\n</g>\n{wm}')
    return wrap(body, f"{W:.0f}", Hh)


def write(name, svg, scale=4):
    (OUT / f"{name}.svg").write_text(svg)
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(OUT / f"{name}.png"), scale=scale)


# ---------------------------------------------------------------- sheet -----
def sheet():
    """One SVG: a row per concept — light icon at 64/32/16, badge at 64/32/16,
    lockup; then the same on a dark panel."""
    ROW, LABEL_W = 132, 190
    LIGHT_W, DARK_W = 900, 720
    W = LABEL_W + LIGHT_W + DARK_W
    Hh = 56 + ROW * len(CONCEPTS)
    el = [f'<rect width="{W}" height="{Hh}" fill="#ffffff"/>',
          f'<rect x="{LABEL_W+LIGHT_W}" width="{DARK_W}" height="{Hh}" fill="{BG_D}"/>']
    font_css = 'font-family="DejaVu Sans, sans-serif"'
    def label(x, y, s, col="#666", size=11, anchor="start"):
        el.append(f'<text x="{x}" y="{y}" {font_css} font-size="{size}" fill="{col}" text-anchor="{anchor}">{s}</text>')
    # column headers
    hx = LABEL_W
    for x, s in [(0, "icon 64 / 32 / 16"), (230, "badge 64 / 32 / 16"), (460, "lockup (light)")]:
        label(hx + x, 34, s)
    for x, s in [(0, "icon on dark 64 / 32 / 16"), (230, "badge"), (330, "lockup (dark)")]:
        label(LABEL_W + LIGHT_W + 24 + x, 34, s, col="#9a9a9a")

    def place(inner, x, y, size, uid):
        el.append(f'<g transform="translate({x},{y}) scale({size/64:.4f})">{inner}</g>')

    for i, (title, slug, fn) in enumerate(CONCEPTS):
        y0 = 56 + i * ROW
        label(20, y0 + 40, title, col="#1a1a1a", size=15)
        if fn.__doc__:
            words, lines, cur = fn.__doc__.split(), [], ""
            for w_ in words:
                if len(cur) + len(w_) > 30:
                    lines.append(cur); cur = w_
                else:
                    cur = (cur + " " + w_).strip()
            lines.append(cur)
            for k, ln in enumerate(lines[:4]):
                label(20, y0 + 60 + 13 * k, ln, size=9.5)
        # light panel
        x = LABEL_W
        place(fn(ACCENT_L, TINT_L, uid=f"l{slug}"), x, y0 + 20, 64, slug)
        place(fn(ACCENT_L, TINT_L, uid=f"l2{slug}"), x + 80, y0 + 52, 32, slug)
        place(fn(ACCENT_L, TINT_L, uid=f"l3{slug}"), x + 124, y0 + 68, 16, slug)
        for size, dx, dy in [(64, 230, 20), (32, 310, 52), (16, 354, 68)]:
            place(f'<rect width="64" height="64" rx="14" fill="{BADGE_BG}"/>' + fn(BADGE_C1, BADGE_C2, uid=f"b{size}{slug}"), x + dx, y0 + dy, size, slug)
        # lockup light (scaled to height 44)
        lk = lockup_svg(fn, FG_L, ACCENT_L, ACCENT_L, TINT_L, uid=f"lk{slug}")
        inner = lk.split(">", 1)[1].rsplit("</svg>", 1)[0]
        el.append(f'<g transform="translate({x+460},{y0+30}) scale(0.72)">{inner}</g>')
        # dark panel
        xd = LABEL_W + LIGHT_W + 24
        place(fn(ACCENT_D, TINT_D, uid=f"d{slug}"), xd, y0 + 20, 64, slug)
        place(fn(ACCENT_D, TINT_D, uid=f"d2{slug}"), xd + 80, y0 + 52, 32, slug)
        place(fn(ACCENT_D, TINT_D, uid=f"d3{slug}"), xd + 124, y0 + 68, 16, slug)
        place(f'<rect width="64" height="64" rx="14" fill="{BADGE_BG}"/>' + fn(BADGE_C1, BADGE_C2, uid=f"db{slug}"), xd + 230, y0 + 20, 64, slug)
        lkd = lockup_svg(fn, FG_D, ACCENT_D, ACCENT_D, TINT_D, uid=f"lkd{slug}")
        inner = lkd.split(">", 1)[1].rsplit("</svg>", 1)[0]
        el.append(f'<g transform="translate({xd+330},{y0+30}) scale(0.72)">{inner}</g>')
        if i < len(CONCEPTS) - 1:
            el.append(f'<line x1="16" x2="{W-16}" y1="{y0+ROW-8}" y2="{y0+ROW-8}" stroke="#e3e3e3"/>')
    return wrap("\n".join(el), W, Hh)


if __name__ == "__main__":
    for title, slug, fn in CONCEPTS:
        write(f"{slug}-icon", icon_svg(fn, ACCENT_L, TINT_L, uid="a"))
        write(f"{slug}-icon-dark", icon_svg(fn, ACCENT_D, TINT_D, uid="b"))
        write(f"{slug}-icon-badge", icon_svg(fn, BADGE_C1, BADGE_C2, bg=BADGE_BG, uid="c"), scale=8)
        write(f"{slug}-logo", lockup_svg(fn, FG_L, ACCENT_L, ACCENT_L, TINT_L, uid="d"))
        write(f"{slug}-logo-dark", lockup_svg(fn, FG_D, ACCENT_D, ACCENT_D, TINT_D, uid="e"))
    write("concept-sheet", sheet(), scale=1.6)
    print("wrote", len(list(OUT.iterdir())), "files to", OUT)
