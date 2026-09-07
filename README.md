# CompareM2

[![unit tests](https://github.com/cmkobel/CompareM2/actions/workflows/unit.yaml/badge.svg)](https://github.com/cmkobel/CompareM2/actions/workflows/unit.yaml) [![https://doi.org/10.1093/bioinformatics/btaf517](https://img.shields.io/badge/doi%20%28OUP%29-10.1093%2Fbioinformatics%2Fbtaf517-blue.svg)](https://doi.org/10.1093/bioinformatics/btaf517)

```bash
pixi global install --channel conda-forge --channel bioconda comparem2
# or: conda install -c conda-forge -c bioconda comparem2

comparem2 *.fna
```

CompareM2 takes microbial genome assemblies — isolates or MAGs, from any
sequencing technology — and produces a single portable HTML report comparing
them: **easy to install, easy to run, easy to interpret.**

No genomes to hand? `comparem2 --demo` runs six bundled *Enterococcus faecium*
plasmids and needs no databases.

## The 14 analyses

| | |
| --- | --- |
| **Quality** | SeqKit (contig lengths, GC, N50), CheckM2 (completeness, contamination) |
| **Taxonomy** | GTDB-Tk |
| **Annotation** | Bakta |
| **Screening** | AMRFinderPlus, MLST |
| **Relatedness** | Mashtree, TreeCluster, skani (all-against-all ANI) |
| **Pangenome** | Panaroo, snp-dists, FastTree |
| **Metabolism** | CarveMe (genome-scale metabolic models), biosynthesis (which building blocks each genome can make) |

All fourteen run across the whole set and land in one report, at a cost that
scales near-linearly with how many assemblies you give it. A hundred genomes is
an ordinary input, not a stress test.

Each section of the report says what the tool does, how to read the specific
columns on screen, and what the result *cannot* tell you — every number quoted
from the tool's own paper and checked against it — and the report ends with a
citation list covering exactly the tools that ran.

## Installing

**Linux only**, because the analysis tools are `linux-64`. The package is the
pipeline alone: Snakemake deploys the tools into two conda environments
(**7.7 GB**) the first time they are needed, and fetches four databases
(**62.5 GB** measured, 60.8 GB of it GTDB-Tk) as the workflow reaches them.
Both defaults are shared across runs and both are movable, which matters on a
cluster with a home quota.

[The documentation](https://comparem2.readthedocs.io) covers all of that, plus
running a subset to skip the 60.8 GB, HPC, and how to read each analysis.

## Status

All **14 of 14** tool command lines have been executed end to end on real
genomes, under the conda deployment that is the only way a tool arrives.
[`STATUS.md`](STATUS.md) has the per-tool table; it tracks *execution*, never
installation, because two tools have resolved to builds that installed cleanly
and crashed on first use.

## Development

There is no hand-written Snakefile. `src/comparem2/catalogue.py` holds the 14
tool specs and the workflow is generated from them, so the CLI, the TUI and the
report all read one source of truth — which is why the unit tests, not an
end-to-end run, are the primary instrument.

```bash
pip install pytest pytest-asyncio textual   # what CI installs; no pixi needed
python -m pytest tests/unit -q              # 221 tests, ~2.5 s

pixi install                                # linux only
pixi run test-fast                          # 4 genomes, no databases needed
pixi run comparem2 --help
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
