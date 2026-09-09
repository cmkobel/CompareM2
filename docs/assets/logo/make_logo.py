#!/usr/bin/env python3
"""Generate CompareM2 logo concepts as SVG (text outlined to paths) and PNG.

Usage: python make_logo.py [outdir] [concept ...]
Concepts: synteny, dotplot, pangenome (default: all).
Requires: fonttools, cairosvg, and Inter (OFL) at FONT.
"""
import sys
from pathlib import Path

import cairosvg
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.ttLib import TTFont

FONT = "/tmp/inter/extras/otf/Inter-SemiBold.otf"
args = sys.argv[1:]
OUT = Path(args[0] if args else "out")
OUT.mkdir(parents=True, exist_ok=True)
WANTED = args[1:] or ["contigs"]

# Colours shared with report.py (--accent light/dark, --fg light/dark)
THEMES = {
    "light": dict(accent="#2b6cb0", tint="#7fb3e8", muted="#c9d3de", fg="#1a1a1a",
                  ribbon=("#a9c9ec", "#d3e4f6")),
    "dark":  dict(accent="#7fb3e8", tint="#3f80c4", muted="#454b52", fg="#e8e8e8",
                  ribbon=("#35608d", "#28466a")),
    "badge": dict(accent="#ffffff", tint="#bcd8f5", muted="#6a9bd3", fg="#ffffff", bg="#2b6cb0",
                  ribbon=("#6ea0d8", "#4f88c8")),
}


def rrect(x, y, w, h, fill, r=3):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}"/>'


# --------------------------------------------------------------------------
# Concept A — synteny: two genomes, homologous blocks share a colour and are
# joined by ribbons; one pair is rearranged (the crossing), one block is unique.
def synteny(c):
    top_y, bot_y, h = 8, 46, 10
    top = [(6, 18, c["accent"]), (28, 14, c["tint"]), (46, 12, c["muted"])]
    bot = [(6, 14, c["tint"]), (24, 18, c["accent"]), (46, 12, c["muted"])]
    els = []
    # ribbons first (solid, theme-specific light tones), lighter one underneath
    y0, y1 = top_y + h, bot_y
    m = (y0 + y1) / 2
    for (tx, tw, _), (bx, bw, _), f in [(top[1], bot[0], c["ribbon"][1]), (top[0], bot[1], c["ribbon"][0])]:
        els.append(
            f'<path fill="{f}" d="M{tx},{y0} C{tx},{m} {bx},{m} {bx},{y1} '
            f'L{bx + bw},{y1} C{bx + bw},{m} {tx + tw},{m} {tx + tw},{y0} Z"/>')
    els += [rrect(x, top_y, w, h, f) for x, w, f in top]
    els += [rrect(x, bot_y, w, h, f) for x, w, f in bot]
    return "\n  ".join(els)


# Concept B — dot plot: whole-genome alignment of one genome against another;
# forward matches on the diagonal, one inverted segment on the anti-diagonal.
def dotplot(c):
    sw = 7
    segs = [((9, 55), (23, 41), c["accent"]),
            ((27, 25), (39, 37), c["tint"]),      # inversion
            ((43, 21), (55, 9), c["accent"])]
    els = [f'<path d="M6,6 V58 H58" fill="none" stroke="{c["muted"]}" stroke-width="3" '
           f'stroke-linecap="round" stroke-linejoin="round"/>']
    for (x0, y0), (x1, y1), f in segs:
        els.append(f'<line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y1}" stroke="{f}" '
                   f'stroke-width="{sw}" stroke-linecap="round"/>')
    return "\n  ".join(els)


# Concept C — pangenome: three genomes as rows, gene clusters as columns that
# line up. Columns present in every genome (core) are solid accent; accessory
# clusters are the tint; a missing cluster leaves the gap.
def pangenome(c):
    cols = [(6, 12), (22, 14), (40, 8), (52, 8)]           # x, width
    present = [(1, 1, 1, 1), (1, 1, 0, 1), (1, 0, 1, 1)]   # rows x cols
    core = [all(r[j] for r in present) for j in range(len(cols))]
    els = []
    for i, row in enumerate(present):
        y = 9 + i * 17
        for j, (x, w) in enumerate(cols):
            if row[j]:
                els.append(rrect(x, y, w, 12, c["accent"] if core[j] else c["tint"], r=4))
            else:
                els.append(rrect(x + w / 2 - 2, y + 4, 4, 4, c["muted"], r=2))
    return "\n  ".join(els)


# Concept D — contigs (chosen): three genomes as rows of contigs, lengths and
# fragmentation differ; rows alternate accent/tint so neighbours read apart.
ROWS = [(11, [(8, 27), (40, 16)]),
        (27, [(8, 13), (26, 30)]),
        (43, [(8, 36)])]
BAR_H = 10


def contigs(c):
    els = []
    for i, (y, segs) in enumerate(ROWS):
        col = c["accent"] if i % 2 == 0 else c["tint"]
        els += [rrect(x, y, w, BAR_H, col, r=BAR_H / 2) for x, w in segs]
    return "\n  ".join(els)


CONCEPTS = {"contigs": contigs, "synteny": synteny, "dotplot": dotplot, "pangenome": pangenome}
# ink bounds (left, top, right, bottom) in the 64-unit box, used to align the
# mark to the wordmark: bottom of ink sits on the baseline, top on cap height
INK = {"contigs": (8, ROWS[0][0], 56, ROWS[-1][0] + BAR_H),
       "synteny": (6, 8, 58, 56), "dotplot": (4.5, 4.5, 59.5, 59.5), "pangenome": (6, 9, 60, 55)}


def icon_svg(draw, c, size=64):
    bg = f'<rect width="64" height="64" rx="14" fill="{c["bg"]}"/>\n  ' if c.get("bg") else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="{size}" height="{size}">\n'
            f'  {bg}{draw(c)}\n</svg>\n')


# --- wordmark: Inter SemiBold outlines ------------------------------------
font = TTFont(FONT)
glyphset = font.getGlyphSet()
cmap = font.getBestCmap()
UPM = font["head"].unitsPerEm
hmtx = font["hmtx"]
TRACK = -0.01 * UPM


def text_path(text, size):
    s = size / UPM
    x, parts = 0.0, []
    for ch in text:
        g = cmap[ord(ch)]
        pen = SVGPathPen(glyphset)
        glyphset[g].draw(pen)
        parts.append((pen.getCommands(), x))
        x += hmtx[g][0] + TRACK
    return parts, x * s, s


def wordmark(fg, accent, size, x0, baseline):
    svg = []
    parts, w1, s = text_path("Compare", size)
    parts2, w2, _ = text_path("M2", size)
    for d, xf in parts:
        svg.append(f'<path fill="{fg}" transform="translate({x0 + xf * s:.2f},{baseline}) scale({s:.5f},{-s:.5f})" d="{d}"/>')
    for d, xf in parts2:
        svg.append(f'<path fill="{accent}" transform="translate({x0 + w1 + xf * s:.2f},{baseline}) scale({s:.5f},{-s:.5f})" d="{d}"/>')
    return "\n  ".join(svg), w1 + w2


CAP = font["OS/2"].sCapHeight / UPM   # Inter: 0.727 em


def lockup_svg(name, c):
    draw = CONCEPTS[name]
    size, pad, gap = 44, 12, 12
    cap_px = CAP * size
    left, top, right, bottom = INK[name]
    sc = cap_px / (bottom - top)              # ink height == cap height
    baseline = pad + cap_px
    tx, ty = pad - left * sc, baseline - bottom * sc
    mark_w = (right - left) * sc
    x0 = pad + mark_w + gap
    wm, wm_w = wordmark(c["fg"], c["accent"], size, x0, baseline)
    W, H = x0 + wm_w + pad, baseline + pad          # descender (0.2 em) fits inside pad
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.1f} {H:.1f}" width="{W:.1f}" height="{H:.1f}">\n'
            f'  <g transform="translate({tx:.3f},{ty:.3f}) scale({sc:.5f})">\n  {draw(c)}\n  </g>\n'
            f'  {wm}\n</svg>\n')


def write(name, svg, png_scale=4):
    (OUT / f"{name}.svg").write_text(svg)
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(OUT / f"{name}.png"), scale=png_scale)


for name in WANTED:
    draw = CONCEPTS[name]
    p = "comparem2" if name == "contigs" else name   # the chosen mark gets the project name
    write(f"{p}-icon", icon_svg(draw, THEMES["light"]))
    write(f"{p}-icon-dark", icon_svg(draw, THEMES["dark"]))
    write(f"{p}-icon-badge", icon_svg(draw, THEMES["badge"]), png_scale=8)
    write(f"{p}-logo", lockup_svg(name, THEMES["light"]))
    write(f"{p}-logo-dark", lockup_svg(name, THEMES["dark"]))
    print("wrote", p)

# Preview sheet: one row per concept, light + dark + small sizes
rows = ""
for name in WANTED:
    name = "comparem2" if name == "contigs" else name
    rows += f"""
<div style="display:grid;grid-template-columns:1fr 1fr;border-bottom:1px solid #ccc">
  <div style="background:#fff;padding:28px 36px;display:flex;gap:32px;align-items:center">
    <span style="width:90px;color:#666;font:14px sans-serif">{name}</span>
    <img src="{name}-logo.svg" height="60"> <img src="{name}-icon.svg" width="64">
    <img src="{name}-icon-badge.svg" width="48"> <img src="{name}-icon-badge.svg" width="24"> <img src="{name}-icon-badge.svg" width="16">
  </div>
  <div style="background:#161616;padding:28px 36px;display:flex;gap:32px;align-items:center">
    <img src="{name}-logo-dark.svg" height="60"> <img src="{name}-icon-dark.svg" width="64"> <img src="{name}-icon-badge.svg" width="48">
  </div>
</div>"""
(OUT / "preview.html").write_text(f'<html><body style="margin:0">{rows}</body></html>')
