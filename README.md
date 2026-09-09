<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/logo/comparem2-logo-dark.png">
    <img src="docs/assets/logo/comparem2-logo.png" alt="CompareM2" width="420">
  </picture>
</h1>

<p align="center">
  <a href="https://github.com/cmkobel/CompareM2/actions/workflows/unit.yaml"><img src="https://github.com/cmkobel/CompareM2/actions/workflows/unit.yaml/badge.svg" alt="unit tests"></a>
  <a href="https://doi.org/10.1093/bioinformatics/btaf517"><img src="https://img.shields.io/badge/doi%20%28OUP%29-10.1093%2Fbioinformatics%2Fbtaf517-blue.svg" alt="doi:10.1093/bioinformatics/btaf517"></a>
</p>

```bash
pixi global install --channel conda-forge --channel bioconda comparem2
# or: conda install -c conda-forge -c bioconda comparem2

comparem2 *.fna
```

Give CompareM2 a set of microbial genome assemblies and it hands back a single
portable HTML report comparing them. Isolates or MAGs, from any sequencing
technology. It aims to be easy to install, easy to run, and easy to interpret.

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

All fourteen run over the whole set and land in the same report. Cost scales
close to linearly with the number of assemblies, so a hundred genomes is an
ordinary input.

Every section explains what its tool does and how to read the columns on
screen, then says what the result cannot tell you. The numbers in those
explanations are quoted from each tool's own paper and checked against it. At
the end comes a citation list covering exactly the tools that ran.

**[See a real one](https://comparem2.readthedocs.io/en/latest/07%20an%20example%20report/)**:
all fourteen tools over eight *Streptococcus mitis* group genomes, including
what the report gets wrong on that set and how it says so.

## Installing

Linux only, since the analysis tools are `linux-64`. What you install is the
pipeline, not the tools. Snakemake deploys those into six conda environments
(1.4 GB) the first time a rule needs one, and fetches four databases (62.5 GB
measured, of which GTDB-Tk alone is 60.8 GB) as the workflow reaches them. Both
locations default to somewhere under `~/.comparem2`, are shared between runs,
and can be moved. You will want to move them on a cluster where home is under
quota.

[The documentation](https://comparem2.readthedocs.io) covers all of that, plus
how to run a subset and skip the 60.8 GB, how to submit to a queue, and how to
read each analysis.

## Status

Every one of the 14 tool command lines has been executed end to end on real
genomes, under the same conda deployment a user gets. [`STATUS.md`](STATUS.md)
has the per-tool table. It records execution and never installation, because
two tools have resolved to builds that installed cleanly and then crashed the
first time they ran.

## Development

There is no hand-written Snakefile. The 14 tool specs live in
`src/comparem2/catalogue.py` and the workflow is generated from them, so the
CLI, the TUI and the report all read one source of truth. That is why the unit
tests are the primary instrument here: a wrong spec yields a Snakefile that
parses cleanly and builds the wrong DAG, which an end-to-end run catches slowly
if at all.

```bash
pip install pytest pytest-asyncio textual   # what CI installs; no pixi needed
python -m pytest tests/unit -q              # 274 tests, ~6 s

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
- **Design**: [`DESIGN.md`](DESIGN.md). What CompareM2 is and why it is shaped this way.
- **Decision log**: [`DECISIONS.md`](DECISIONS.md). How it got here, including what was reversed and what went wrong.
- **Status**: [`STATUS.md`](STATUS.md). What has actually been run.
- **Issues**: [github.com/cmkobel/CompareM2/issues](https://github.com/cmkobel/CompareM2/issues)
