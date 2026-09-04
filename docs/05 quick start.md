# Quick start

## 1) Get it

Linux only — the analysis tools are `linux-64`.

```bash
conda install -c conda-forge -c bioconda comparem2
```

To work on CompareM2 rather than just use it, install from git with
[pixi](https://pixi.prefix.dev/latest/#installation) instead — see
[installation](10 installation.md). Every `cm2` below becomes `pixi run cm2`.

## 2) Check it works

A real run over four *Enterococcus faecium* genomes, which needs no databases
because none of the four tools it runs has one:

```bash
cm2 tests/E._faecium/*.fna --until seqkit mashtree treecluster skani
```

That writes `results_comparem2/report.html`. Open it. From a git checkout,
`pixi run test-fast` does the same thing with the genomes already unpacked, and
`pixi run pytest` runs the 200 unit tests — pure Python, no tools or databases.

!!! note "The first run builds the tool environments"
    The fourteen tools do not come with the package — Snakemake deploys them
    into two conda environments, 7.7 GB, the first time they are needed. That
    happens before any job starts, so the first run is quiet for about a minute
    while it works. `cm2 --setup` does it deliberately and takes no assemblies;
    see [installation](10 installation.md).

## 3) Run it on your own genomes

```bash
cm2 path/to/*.fna
```

CompareM2 prints what it is about to do, including how much database it needs to
download, before downloading anything:

```
4 assemblies, 14 tools
to download: checkm2, gtdb, bakta-light, amrfinder (62.5 GB + 2 of unknown size) -> ~/.comparem2/databases
```

**60.8 GB of that 62.5 GB is GTDB-Tk.** If you do not need taxonomic
assignment, skip it and the download collapses to under 2 GB:

```bash
cm2 *.fna --until seqkit checkm2 bakta amrfinder mlst \
                  mashtree treecluster skani panaroo snp-dists fasttree
```

## 4) Read the report

`results_comparem2/report.html` — one self-contained file with no external
assets, so it survives being emailed or copied off a cluster.

Each section has a **"What this is, and how to read it"** block, collapsed by
default. Open it before drawing a conclusion from the numbers: it says what the
columns mean, what the tool's own error is, and what the result cannot show.

## Useful next steps

```bash
# What would run, without running it
cm2 *.fna --dry-run

# Only one analysis and its dependencies
cm2 *.fna --until fasttree

# Change a tool's argument
cm2 *.fna --set treecluster--threshold=0.1

# Interactive interface
cm2 *.fna --tui

# More cores
cm2 *.fna --cores 24
```

See [usage](20 usage.md) for the rest.
