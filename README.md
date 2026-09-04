# CompareM2

[![unit tests](https://github.com/cmkobel/CompareM2/actions/workflows/unit.yaml/badge.svg)](https://github.com/cmkobel/CompareM2/actions/workflows/unit.yaml) [![https://doi.org/10.1093/bioinformatics/btaf517](https://img.shields.io/badge/doi%20%28OUP%29-10.1093%2Fbioinformatics%2Fbtaf517-blue.svg)](https://doi.org/10.1093/bioinformatics/btaf517)

```bash
conda install -c conda-forge -c bioconda comparem2
comparem2 *.fna
```

No genomes to hand? `comparem2 --demo` runs six bundled *Enterococcus faecium*
plasmids and needs no databases.

CompareM2 takes microbial genome assemblies — isolates or MAGs, from any
sequencing technology — and produces a single portable HTML report comparing
them: **easy to install, easy to run, easy to interpret.**

## Highlights

**Built for many genomes at once.** Fourteen analyses across a whole set of
assemblies, in one report, at a cost that scales near-linearly with how many you
give it — the benchmark result the design is built on. A hundred genomes is an
ordinary input, not a stress test.

**Every tool is one declarative spec.** There is no hand-written Snakefile —
`src/comparem2/catalogue.py` holds the 14 specs and the workflow is generated
from them, so the CLI, the TUI and the report all read the same source of truth.

**An interpretative report.** Each section carries what the tool does, how to
read the specific columns on screen, and what the result cannot tell you — every
number quoted from the tool's own paper and checked against it — plus a
citations list of every paper behind the tools that ran.

At a glance:

| | |
| --- | --- |
| Analyses | **14** |
| Conda environments | **2** — thirteen tools co-solve, CheckM2 cannot |
| Software on disk | **7.7 GB** measured, deployed on first run |
| Databases | **4**, fetched automatically (Bakta light) |
| Report | **one** self-contained HTML file |
| Runtime | **Python** only |
| Unit tests | **212** |

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
all of them under the conda deployment that is the only way a tool arrives.
GTDB-Tk was the last, and it took six defects and a database release change to
get there — its rule had never been run, and two of its steps were described in
comments that were not true.

[`STATUS.md`](STATUS.md) has the per-tool table. It tracks *execution*, never
installation — a clean `pixi install` says nothing about whether a tool runs,
and two tools have resolved to builds that installed cleanly and crashed.

For development, from a checkout rather than the package:

```bash
python -m pytest tests/unit -q   # 212 tests, no pixi required

pixi install                     # linux only
pixi run test-fast               # 4 genomes, no databases needed
pixi run cm2 --help
```

## Citation

Kobel C.M., Aho V.T.E., Øyås O., Nørskov-Lauritsen N., Woodcroft B.J., Pope P.B.
CompareM2 is a genomes-to-report pipeline for comparing microbial genomes.
*Bioinformatics* 41(9), btaf517 (2025).
[doi:10.1093/bioinformatics/btaf517](https://doi.org/10.1093/bioinformatics/btaf517)

Please also cite the tools you used. The report's own "Methods and citations"
section lists them for the run you did.

## Links

- **Documentation**: [comparem2.readthedocs.io](https://comparem2.readthedocs.io)
- **Design**: [`DESIGN.md`](DESIGN.md) — what CompareM2 is and why it is shaped this way
- **Decision log**: [`DECISIONS.md`](DECISIONS.md) — how it got here, including what was reversed and what went wrong
- **Status**: [`STATUS.md`](STATUS.md) — what has actually been run
- **Issues**: [github.com/cmkobel/CompareM2/issues](https://github.com/cmkobel/CompareM2/issues)
