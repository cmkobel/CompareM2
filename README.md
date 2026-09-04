# CompareM2 v3

[![unit tests](https://github.com/cmkobel/CompareM2/actions/workflows/unit.yaml/badge.svg)](https://github.com/cmkobel/CompareM2/actions/workflows/unit.yaml) [![https://doi.org/10.1093/bioinformatics/btaf517](https://img.shields.io/badge/doi%20%28OUP%29-10.1093%2Fbioinformatics%2Fbtaf517-blue.svg)](https://doi.org/10.1093/bioinformatics/btaf517)

```bash
conda install -c conda-forge -c bioconda comparem2=3
comparem2 *.fna
```

> On bioconda since 2026-09-04. `=3` is a fuzzy match on the v3 series, so 3.0.1
> and later still satisfy it; it is there so that a resolver conflict fails
> loudly instead of quietly giving you v2 — see
> [installation](https://comparem2.readthedocs.io/en/latest/10%20installation/).
> The paper describes v2, whose documentation remains at
> [comparem2.readthedocs.io/en/stable](https://comparem2.readthedocs.io/en/stable/).

CompareM2 takes microbial genome assemblies — isolates or MAGs — and produces a
single portable HTML report comparing them. v3 keeps only that philosophy from
v2: **easy to install, easy to run, easy to interpret.**

## What is different in v3

**A wide view, not a deep one.** v3 is triage across many assemblies rather than
a deep dive into each. Deep functional and metabolic analysis is where DRAM2 and
nf-core/funcscan are already strong; a fast wide sweep is not well served, and
it is what v2's benchmark result — near-linear scaling with input count —
actually supports.

That choice is what makes the rest possible:

| | v2 | v3 |
| --- | --- | --- |
| Analyses | 30+ | **14** |
| Conda environments | 25 | **2** — thirteen tools co-solve, CheckM2 cannot |
| Software on disk | — | **7.7 GB** measured, deployed on first run |
| Databases | 7 (incl. 84 GB Bakta full) | **4** (Bakta light) |
| Report | R, RMarkdown, pandoc | **Python**, one HTML file |
| Runtime | R + Python | **Python only** |
| Unit tests | none | **200** |

**Every tool is one declarative spec.** There is no hand-written Snakefile —
`src/comparem2/catalogue.py` holds the 14 specs and the workflow is generated
from them, so the CLI, the TUI and the report all read the same source of truth.

**The report explains itself.** Each section carries what the tool does, how to
read the specific columns on screen, and what the result cannot tell you — every
number quoted from the tool's own paper and checked against it — plus a
citations list of every paper behind the tools that ran.

## The 14 tools

| | |
| --- | --- |
| **Quality** | SeqKit (contig lengths, GC, N50), CheckM2 (completeness, contamination) |
| **Taxonomy** | GTDB-Tk |
| **Annotation** | Bakta |
| **Screening** | AMRFinderPlus, MLST |
| **Relatedness** | Mashtree, TreeCluster, skani (all-against-all ANI) |
| **Pangenome** | Panaroo, snp-dists, FastTree |
| **Metabolism** | CarveMe (genome-scale metabolic models), biosynthesis (which building blocks each genome can make) |

## Status

Linux-only; the tools are `linux-64`. Unit tests run anywhere.

**14 of 14** tool command lines have been executed end to end on real genomes,
all of them under the conda deployment that is now the only way a tool arrives.
GTDB-Tk was the last, and it took six defects and a database release change to
get there — its rule had never been run, and two of its steps were described in
comments that were not true.

[`STATUS.md`](STATUS.md) has the per-tool table. It tracks *execution*, never
installation — a clean `pixi install` says nothing about whether a tool runs,
and two tools have resolved to builds that installed cleanly and crashed.

For development, from a checkout rather than the package:

```bash
python -m pytest tests/unit -q   # 200 tests, no pixi required

pixi install                     # linux only
pixi run test-fast               # 4 genomes, no databases needed
pixi run cm2 --help
```

## Citation

The paper describes v2:

Kobel C.M., Aho V.T.E., Øyås O., Nørskov-Lauritsen N., Woodcroft B.J., Pope P.B.
CompareM2 is a genomes-to-report pipeline for comparing microbial genomes.
*Bioinformatics* 41(9), btaf517 (2025).
[doi:10.1093/bioinformatics/btaf517](https://doi.org/10.1093/bioinformatics/btaf517)

Please also cite the tools you used. The report's own "Methods and citations"
section lists them for the run you did.

## Links

- **Design**: [`DESIGN.md`](DESIGN.md) — what v3 is and why
- **Decision log**: [`DECISIONS.md`](DECISIONS.md) — how it got here, including what was reversed and what went wrong
- **Status**: [`STATUS.md`](STATUS.md) — what has actually been run
- **v3 documentation**: [comparem2.readthedocs.io/en/latest](https://comparem2.readthedocs.io/en/latest/)
- **v2 documentation**: [comparem2.readthedocs.io/en/stable](https://comparem2.readthedocs.io/en/stable/) — what the paper describes
- **Issues**: [github.com/cmkobel/CompareM2/issues](https://github.com/cmkobel/CompareM2/issues)
