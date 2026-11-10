#!/usr/bin/env bash
# The arm probe: CompareM2 run whole on this machine, from this checkout, and
# its numbers set beside thylakoid's. For a machine CI cannot stand in for —
# the laptop under macOS, or Asahi Linux at 16 KiB pages — not for routine use.
#
#     tests/environments/arm_probe.sh [workdir]      # default ~/cm2-arm-probe
#
# Needs pixi on PATH, and about 10 GB in workdir and ~/.comparem2 together:
# the five tool environments, and the bakta, CheckM2 and AMRFinder databases
# (about 3.4 GB of downloads). GTDB-Tk is left out: no gtdbtk >=2.7 installs on
# arm until bioconda-recipes#71426 is merged.
#
# It runs, in order: `--demo` twice and `check.py demo`, as the CI job does;
# then the thirteen other tools over tests/E._faecium. On macOS, `comparem2`
# refuses to start, so the probe calls `main()` with that check replaced. That
# is the one difference from what a user runs, and the reason this is a probe.
#
# It never stops on a failed step: a probe that dies on the first failure
# tells you about one tool. Everything lands in workdir, logs included, and
# the summary at the end says which expectations held.

set -u
here=$(cd "$(dirname "$0")/../.." && pwd)
work=${1:-$HOME/cm2-arm-probe}
mkdir -p "$work"

case "$(uname -s)-$(uname -m)" in
  Darwin-arm64)  platform=osx-arm64 ;;
  Linux-aarch64) platform=linux-aarch64 ;;
  Linux-x86_64)  platform=linux-64 ;;
  *) echo "no conda platform for $(uname -sm)"; exit 1 ;;
esac
echo "$(uname -sm), page size $(getconf PAGESIZE), conda platform $platform"
echo "checkout $(git -C "$here" rev-parse --short HEAD), workdir $work"

# The pipeline's own environment, for this platform. pixi.toml lists linux-64
# only, so it is copied with the platform swapped and its tasks dropped.
mkdir -p "$work/pixi"
sed -e "s/^platforms = .*/platforms = [\"$platform\"]/" \
    -e '/^\[tasks\]/,/^\[activation.env\]/{/^\[activation.env\]/!d;}' \
    "$here/pixi.toml" > "$work/pixi/pixi.toml"
pixi install --manifest-path "$work/pixi/pixi.toml" || exit 1

# `main()` with the platform refusal switched off, and nothing else changed.
cm2() {
  (cd "$here" && PYTHONPATH="$here/src" pixi run --manifest-path "$work/pixi/pixi.toml" \
    python -c 'import sys, comparem2.cli as c
c.unsupported_platform = lambda: None
sys.exit(c.main(sys.argv[1:]))' "$@")
}

timed() {  # timed <log> <command...>: run, log, report exit and seconds
  local log=$1; shift
  local start; start=$(date +%s)
  "$@" > "$log" 2>&1
  local rc=$?
  echo "  exit $rc in $(( $(date +%s) - start )) s — $log"
  return $rc
}

echo "== demo, twice"
timed "$work/demo1.log" cm2 --demo --output "$work/demo"
timed "$work/demo2.log" cm2 --demo --output "$work/demo"
(cd "$here" && pixi run --manifest-path "$work/pixi/pixi.toml" \
  python tests/environments/check.py demo "$work/demo" --rerun-log "$work/demo2.log")

echo "== thirteen tools over tests/E._faecium"
(cd "$here" && unzip -o -q tests/E._faecium/fna.zip -d tests/E._faecium)
timed "$work/full.log" cm2 "$here"/tests/E._faecium/*.fna --output "$work/full" \
  --until seqkit checkm2 bakta amrfinder mlst mashtree treecluster skani \
          panaroo snp-dists fasttree carveme biosynthesis
grep -E "error in rule|Error in rule" "$work/full.log" | sort -u

echo "== against thylakoid"
pixi run --manifest-path "$work/pixi/pixi.toml" python - "$work/full" <<'EOF'
"""Each line: what this run gave, what thylakoid gave, and whether they agree.

Thylakoid's numbers are the ones STATUS.md records for these four genomes
(linux-64; bakta 1.12.1, checkm2 1.1.0, panaroo with the -aS fix, which 1.8.1
carries). CarveMe's reaction counts are printed but not judged: the MILP has
many optima, and which one a solver build returns is not a defect.
"""
import re
import sys
from pathlib import Path

out = Path(sys.argv[1])
S = ["116_2", "116_2_duplicate", "E8202", "SRR24"]
held = failed = 0


def expect(what, got, want):
    global held, failed
    ok = got == want
    held, failed = held + ok, failed + (not ok)
    print(f"  {'ok  ' if ok else 'DIFF'} {what}: {got}" + ("" if ok else f" (thylakoid: {want})"))


def read(path):
    p = out / path
    return p.read_text() if p.exists() else ""


def rows(path):
    return [l.split("\t") for l in read(path).splitlines()]


for s in S:
    gff = read(f"samples/{s}/bakta/{s}.gff3")
    cds = sum(1 for l in gff.splitlines() if l.split("\t")[2:3] == ["CDS"])
    expect(f"bakta CDS {s}", cds, {"E8202": 3180, "SRR24": 2815}.get(s, 2588))
faa = [read(f"samples/{s}/bakta/{s}.faa") for s in S[:2]]
expect("bakta .faa duplicate pair identical", faa[0] == faa[1] and bool(faa[0]), True)

for s, n in zip(S, (7, 7, 11, 8)):
    expect(f"amrfinder genes {s}", max(len(rows(f"samples/{s}/amrfinder/amrfinder.tsv")) - 1, 0), n)

st = {Path(r[0]).stem: r[2] for r in rows("mlst/mlst.tsv") if len(r) > 2}
expect("mlst", [st.get(s) for s in S], ["32", "32", "117", "78"])

q = {r[0]: (r[1], r[2]) for r in rows("checkm2/quality_report.tsv")[1:]}
expect("checkm2", [q.get(s) for s in S],
       [("100.0", "0.47"), ("100.0", "0.47"), ("100.0", "3.42"), ("100.0", "1.54")])

ani = rows("skani/ani.tsv")[1:]
names = [Path(r[0]).stem for r in ani]
if all(s in names for s in S):
    m = {(a, b): ani[names.index(a)][1 + names.index(b)] for a in S for b in S}
    expect("skani pair / E8202 / SRR24", [m[("116_2", "116_2_duplicate")], m[("116_2", "E8202")],
                                          m[("116_2", "SRR24")]], ["100.00", "99.14", "98.91"])
else:
    expect("skani", "missing", "present")

stats = {l.split("\t")[0]: l.split("\t")[-1] for l in read("panaroo/summary_statistics.txt").splitlines()}
expect("panaroo clusters / core", (stats.get("Total genes"), stats.get("Core genes")), ("3778", "2092"))

d = rows("snp-dists/snp-dists.tsv")
if d:
    head = d[0][1:]
    m = {(r[0], head[i]): r[1 + i] for r in d[1:] for i in range(len(head))}
    # The three strains' distances sorted: which pair is which is in STATUS.md,
    # and the alignment's row order is not stable from run to run.
    expect("snp-dists duplicate pair", m.get(("116_2", "116_2_duplicate")), "0")
    expect("snp-dists between strains, sorted",
           sorted(int(m.get(p, -1)) for p in (("116_2", "E8202"), ("116_2", "SRR24"), ("E8202", "SRR24"))),
           [1699, 1893, 2645])
else:
    expect("snp-dists", "missing", "present")

a, b = (re.escape(x) + r":0(\.0+)?" for x in S[:2])
expect("fasttree duplicate pair joined at 0", bool(re.search(rf"\(({a},{b}|{b},{a})\)", read("fasttree/fasttree.newick"))), True)

want = {"116_2": (11, 6, 11, 2), "116_2_duplicate": (11, 6, 11, 2),
        "E8202": (10, 9, 9, 2), "SRR24": (10, 7, 11, 2)}
for s in S:
    r = rows(f"samples/{s}/biosynthesis/{s}.tsv")
    col = r[0].index("verdict") if r else 0
    v = [x[col] for x in r[1:]]
    got = tuple(v.count(k) for k in ("de_novo", "upstream", "none", "absent"))
    expect(f"biosynthesis {s} de novo/upstream/none/absent", got, want[s])
    xml = read(f"samples/{s}/carveme/{s}.xml")
    print(f"       carveme {s}: {xml.count('<reaction ')} reactions")

print(f"{held} held, {failed} differ")
EOF
