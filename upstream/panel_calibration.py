#!/usr/bin/env python3
"""What the biosynthesis panel answers for genomes whose answer is known.

The calibration is a table of expected answers, one per genome CarveMe
bundles, and what this reports is where a draft departs from it. It used to be
a count — "29 of 30, for the curated `iML1515` and for drafts of five
prototrophs" — and the count hid two things. Three of those five were CarveMe
1.5.0 drafts that ship with the package, not drafts carved the way this
pipeline carves. And 29 is not every prototroph's right answer: an organism
with no menaquinone pathway should not make menaquinone-8.

**What each genome should not make, and why:**
- `adocbl` for every genome: not in the carving universe at all
  (`biosynthesis.UNIVERSE_CEILING`).
- `mqn8` for *P. aeruginosa* and *R. solanacearum*: neither bundled proteome
  annotates a menaquinone-pathway gene, only the `UbiE` methyltransferase that
  ubiquinone synthesis shares. *K. pneumoniae* and *S. aureus* annotate the
  whole Men pathway, checked the same way. Ubiquinone-9 is *P. aeruginosa*'s
  major quinone (Vo et al. 2020, J Biol Chem, PMC7335794), and Q-8 the
  predominant one in the *Ralstonia* species described since 2023.
- `asn__L` for *R. solanacearum*: its proteome annotates no asparagine
  synthetase, and does annotate the GatCAB amidotransferase, which makes
  Asn-tRNA from Asp-tRNA without free asparagine. A model has no
  representation of that route, so "no route to free asparagine" is the
  faithful reading of the genome. It is listed as expected so that it is not
  counted as a carving error, and the reason says why.
- *M. genitalium* is the negative control: a genome-reduced parasite, which
  must not come out making anything from a minimal medium.

Anything else not *de novo* is a **false dependency**, and anything listed
that is *de novo* is a **false prototrophy**.

Measured 2026-11-11 on thylakoid, CarveMe 1.6.6 through this repo's
`carve_scip.py`, ReFramed 1.6.0:

    Ecoli_K12_MG1655       29 de novo   as expected
    Bsubtilis_168          29 de novo   as expected
    Paeruginosa_PAO1       28 de novo   as expected
    Soneidensis_MR1        28 de novo   false dependency: pro__L
    Rsolanacearum_GMI1000  27 de novo   as expected
    M_genitalium_G37        0 de novo   as expected (negative control)

*S. oneidensis*'s proline is not traced. Its draft carries `GLU5K`, `G5SD`
and `P5CR`, the same as the *E. coli* draft, which does make proline.

**The inputs ship with CarveMe.** `carveme/data/benchmark/fasta/` holds the
proteomes of the five organisms the CarveMe paper benchmarks against plus
*M. genitalium*, so nothing has to be downloaded and nothing annotated.
Carving took 24–48 s a genome on thylakoid, and 336 s for *M. genitalium*.

Self-contained, like `carve_scip.py` and `biosynthesis.py` — runs under the
tool environment's own python. It imports `biosynthesis` by path rather than
from `comparem2`, because under `--use-conda` the package is not installed in
the environment that has ReFramed.

**Run it with that environment activated, not just with its interpreter**:
`carve` shells out to DIAMOND, so the environment's `bin` has to be on `PATH`
or the carve step fails with "Unable to run diamond". Models already in the
workdir are reused, so scoring alone needs no DIAMOND at all.

    conda activate <the carveme env>
    python3 panel_calibration.py --workdir /tmp/cal
    python3 panel_calibration.py --workdir /tmp/cal --model iML1515.xml

A `--model` has no expected answer and is scored without one.
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIOSYNTHESIS = HERE.parent / "src" / "comparem2" / "biosynthesis.py"
CARVE_SCIP = HERE.parent / "src" / "comparem2" / "carve_scip.py"

CEILING = "not in the carving universe"
NO_MENAQUINONE = "no menaquinone pathway annotated; respires with ubiquinone"
INDIRECT_ASN = "no asparagine synthetase annotated; Asn-tRNA by GatCAB"

# Every bundled genome, named as the files are, and what it should not make
# beyond the universe ceiling. See the module docstring for each reason.
EXPECTED = {
    "Ecoli_K12_MG1655": {},
    "Bsubtilis_168": {},
    "Paeruginosa_PAO1": {"mqn8": NO_MENAQUINONE},
    "Soneidensis_MR1": {},
    "Rsolanacearum_GMI1000": {"mqn8": NO_MENAQUINONE, "asn__L": INDIRECT_ASN},
}
NEGATIVE_CONTROL = "M_genitalium_G37"
GENOMES = (*EXPECTED, NEGATIVE_CONTROL)


def load_biosynthesis():
    spec = importlib.util.spec_from_file_location("biosynthesis", BIOSYNTHESIS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bundled_fasta() -> Path:
    """Where CarveMe keeps the benchmark proteomes in the installed package."""
    import carveme

    return Path(carveme.__file__).parent / "data" / "benchmark" / "fasta"


def carve(faa: Path, output: Path) -> None:
    """Through `carve_scip.py`, so this measures the pipeline's own path."""
    if output.exists():
        print(f"  {output.name} exists, reusing", file=sys.stderr)
        return
    subprocess.run([sys.executable, str(CARVE_SCIP),
                    "--faa", str(faa), "--output", str(output)], check=True)


def expected_not(bio, name: str) -> dict[str, str] | None:
    """What `name` should not make, with reasons; None for no expectation."""
    if name not in EXPECTED:
        return None
    return {**{b: CEILING for b in bio.UNIVERSE_CEILING}, **EXPECTED[name]}


def departures(bio, found: dict[str, str], name: str) -> list[str]:
    """Where a genome's verdicts depart from what it should answer."""
    if name == NEGATIVE_CONTROL:
        return [f"false prototrophy: {b}" for b, v in found.items()
                if v == bio.DE_NOVO]
    expected = expected_not(bio, name)
    if expected is None:
        return []
    out = [f"false dependency: {b} ({v})" for b, v in found.items()
           if v != bio.DE_NOVO and b not in expected]
    out += [f"false prototrophy: {b}" for b in expected
            if found.get(b) == bio.DE_NOVO]
    return out


def score(bio, path: Path) -> dict[str, str]:
    from reframed import load_cbmodel

    return bio.verdicts(bio._Probe(load_cbmodel(str(path), flavor="bigg")))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="panel_calibration", description=__doc__)
    p.add_argument("--workdir", type=Path, required=True,
                   help="where to carve into; existing models are reused")
    p.add_argument("--model", type=Path, action="append", default=[],
                   help="an extra model to score, e.g. a curated one. Repeatable")
    p.add_argument("--genome", action="append", default=[],
                   help=f"override the bundled set (default: {', '.join(GENOMES)})")
    args = p.parse_args(argv)

    bio = load_biosynthesis()
    args.workdir.mkdir(parents=True, exist_ok=True)

    models = list(args.model)
    fasta = bundled_fasta()
    for name in (args.genome or GENOMES):
        faa = fasta / f"{name}.faa"
        if not faa.exists():
            print(f"no bundled proteome {faa}", file=sys.stderr)
            return 1
        # Into the workdir, and not next to the proteome: `carve_scip.py` links
        # its input into the output's directory.
        output = args.workdir / f"{name}.xml"
        carve(faa, output)
        models.append(output)

    print(f"{'model':<22} {'de novo':>8}  departures from the expected answer")
    total = 0
    for path in models:
        name = Path(path).stem
        found = score(bio, Path(path))
        off = departures(bio, found, name)
        total += len(off)
        de_novo = sum(v == bio.DE_NOVO for v in found.values())
        known = (name == NEGATIVE_CONTROL) or name in EXPECTED
        print(f"{name:<22} {de_novo:>8}  "
              + ("; ".join(off) if off else
                 "as expected" if known else "no expectation"))
    print(f"\n{total} departure(s) across {len(models)} model(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
