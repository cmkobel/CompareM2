# Quick start

Linux only — the analysis tools are `linux-64`.

```bash
conda install -c conda-forge -c bioconda comparem2
```

Working on CompareM2 rather than using it? Install from git with pixi instead —
see [installation](10 installation.md), and read every `cm2` below as
`pixi run cm2`.

## 1) Run it

On your own assemblies:

```bash
cm2 *.fna
```

Or, to see a report before you commit any genomes of your own, on six
*Enterococcus faecium* plasmids that ship with the package:

```bash
cm2 --demo
```

That downloads nothing and runs the four analyses that need no database. One of
the seven inputs is another one again under a different name, so the tree and
the ANI matrix have a pair in them that must come out identical.

CompareM2 says what it is about to do, and how much database it needs, before
fetching anything:

```
4 assemblies, 14 tools
to download: checkm2, gtdb, bakta-light, amrfinder (62.5 GB + 2 of unknown size) -> ~/.comparem2/databases
```

**60.8 GB of that 62.5 GB is GTDB-Tk alone.** If you do not need taxonomic
assignment, name the analyses you do want instead and the download falls under
2 GB — see [running a subset](20 usage.md#running-a-subset).

The first run is then quiet for about a minute while Snakemake builds the two
tool environments. That is expected, and `cm2 --setup` does it in advance.

## 2) Read the report

`results_comparem2/report.html` — one self-contained file with no external
assets, so it survives being emailed or copied off a cluster.

Each section carries a **"What this is, and how to read it"** block, collapsed
by default. Open it before drawing a conclusion from the numbers: it says what
the columns mean, what the tool's own error is, and what the result cannot show.

---

[Installation](10 installation.md) covers databases, the tool environments and
HPC. [Usage](20 usage.md) covers the full CLI, the TUI and passthrough
parameters.
