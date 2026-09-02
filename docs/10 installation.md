# Installation

!!! warning "v3 is in development"
    No bioconda package and no container image yet. Both existed for v2 and
    will return; for now, install from the branch.

## Requirements

  - **Linux.** The analysis tools are `linux-64` only. On macOS you can run the
    unit tests and render reports, but not the pipeline.
  - **[pixi](https://pixi.prefix.dev/latest/#installation)** for the
    environment.
  - **Disk.** 1.58 GB of software, plus databases — see below. GTDB-Tk alone is
    141.4 GB.
  - **RAM.** GTDB-Tk's classify step is the peak; its own paper reports under
    55 GB for the v2 divide-and-conquer placement. Without GTDB-Tk, far less.

## Install

```bash
git clone -b v3 https://github.com/cmkobel/comparem2.git
cd comparem2
pixi install
pixi run pytest        # 75 unit tests, no databases needed
```

!!! note "Moving the directory invalidates the environment"
    Conda bakes the absolute prefix into shebangs and RPATHs, so after moving
    the checkout you need `rm -rf .pixi && pixi install`.

## Environments

v3 solves **one** environment holding twelve of the thirteen tools. v2 shipped
25.

The exception is **CheckM2**, which pins DIAMOND 2.1.x while current Bakta needs
2.2.x — they cannot co-solve. CheckM2 gets its own environment, and because
Snakemake's shell does not inherit an interactive PATH, it needs an absolute
launcher:

```bash
pixi run cm2 *.fna \
  --isolated-launcher "$HOME/.pixi/bin/pixi run -e {tool}"
```

`{tool}` is substituted with the tool name. Without this, the CheckM2 step fails
with `command not found`.

## Databases

Downloaded on first use into `databases/` (change with `-d`). Sizes are measured
from `content-length`, not estimated:

| Database | Size | Needed by |
|---|---|---|
| GTDB r226 | **141.4 GB** | `gtdbtk` |
| CheckM2 | 1.7 GB | `checkm2` |
| Bakta light | ~1.3 GB (documented, unmeasured) | `bakta` |
| AMRFinder | unmeasured | `amrfinder` |

**GTDB-Tk is 91% of the total.** It stays in the default path because it is the
authoritative answer to "what is this genome", but if you do not need taxonomy,
leaving it out is the single biggest saving available:

```bash
pixi run cm2 *.fna --until seqkit checkm2 bakta amrfinder mlst \
                           mashtree treecluster skani panaroo snp-dists fasttree
```

CompareM2 prints the total it is about to download before it starts, and says
explicitly when a database's size is unknown rather than silently omitting it
from the sum.

Bakta uses the **light** database (~1.3 GB / 3.9 GB on disk) rather than v2's
full one (30 GB / 84 GB). That saves 29 GB for less specific functional
annotation, which a wide view can absorb — but note the Bakta paper's
annotation-quality figures are measured on the full database and do not
transfer.

## HPC

Snakemake's SLURM executor plugin is installed. Point CompareM2 at more cores
and let Snakemake submit:

```bash
pixi run cm2 *.fna --cores 64
```

Job submission is Snakemake's business, not CompareM2's — the generated
Snakefile lives at `<output>/.comparem2/Snakefile` and takes any Snakemake
profile you already use.

## Verification status

A clean `pixi install` says nothing about whether the pipeline runs. Both Bakta
and Panaroo once resolved to years-old builds that installed cleanly and crashed
on first use, which is why every tool now carries a minimum version.

`DESIGN.md` in the repository tracks which command lines have actually been
**executed** on real genomes. At the time of writing, 12 of 13 have; GTDB-Tk is
outstanding, and two commands changed after their verification runs and need
re-checking.
