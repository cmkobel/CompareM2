# Quick start

!!! warning "v3 is in development"
    There is no released package yet. Install from the `v3` branch.

## 1) Get it

Linux only — the analysis tools are `linux-64`. With
[pixi](https://pixi.prefix.dev/latest/#installation):

```bash
git clone -b v3 https://github.com/cmkobel/comparem2.git
cd comparem2
pixi install
```

## 2) Check it works

No databases and no tools needed for this — the unit tests are pure Python:

```bash
pixi run pytest
```

Then a real run over four *Enterococcus faecium* genomes shipped with the
repository. This needs no databases either, because none of the four tools it
runs has one:

```bash
pixi run test-fast
```

That produces `cm2_test-fast/report.html`. Open it.

## 3) Run it on your own genomes

```bash
pixi run cm2 path/to/*.fna
```

CompareM2 prints what it is about to do, including how much database it needs to
download, before downloading anything:

```
4 assemblies, 13 tools, databases: 143.2 GB + 2 database(s) of unknown size
```

**143 GB is not a typo, and 141 GB of it is GTDB-Tk.** If you do not need
taxonomic assignment, skip it and the download collapses to under 2 GB:

```bash
pixi run cm2 *.fna --until seqkit checkm2 bakta amrfinder mlst \
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
pixi run cm2 *.fna --dry-run

# Only one analysis and its dependencies
pixi run cm2 *.fna --until fasttree

# Change a tool's argument
pixi run cm2 *.fna --set treecluster--threshold=0.1

# Interactive interface
pixi run cm2 *.fna --tui

# More cores
pixi run cm2 *.fna --cores 24
```

See [usage](20 usage.md) for the rest.
