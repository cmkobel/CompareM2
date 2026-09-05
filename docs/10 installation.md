# Installation

CompareM2 is on Bioconda, so **pixi** or **conda** installs it. However it
arrives, the package is the pipeline alone — the fourteen analysis tools are
not in it, and Snakemake deploys them into two conda environments the first
time they are needed. There is no flag for whether to do that; it always
happens, into `--conda-prefix`.

## Requirements

  - **Linux.** The analysis tools are `linux-64` only. On macOS you can run the
    unit tests and render reports, but not the pipeline.
  - **[pixi](https://pixi.prefix.dev/latest/#installation)** or **conda**, to
    install the package. You do not need a conda of your own: `conda` is a
    dependency *of* the package, because Snakemake shells out to it to build
    the tool environments, so it arrives with CompareM2 either way.
  - **Disk.** **7.7 GB** of tool environments, plus databases — see below.
    GTDB-Tk alone is 60.8 GB to download and 94 GB unpacked.
  - **RAM.** GTDB-Tk's classify step is the peak; its own paper reports under
    55 GB for GTDB-Tk 2's divide-and-conquer placement. Without GTDB-Tk, far less.

## With pixi

To install globally, so that `comparem2` is on `PATH` everywhere:

```bash
pixi global install --channel conda-forge --channel bioconda comparem2
comparem2 *.fna
```

Or to add it to a workspace, alongside whatever else that project needs:

```bash
pixi workspace channel add conda-forge
pixi workspace channel add bioconda
pixi add comparem2
pixi run comparem2 *.fna
```

## With conda

Into the environment you have active:

```bash
conda install -c conda-forge -c bioconda comparem2
comparem2 *.fna
```

Or into a new one:

```bash
conda create -n comparem2 -c conda-forge -c bioconda comparem2
conda activate comparem2
```

!!! note "The channels are given explicitly on purpose"
    Bioconda's own instructions omit `-c`/`--channel` because they assume you
    have already
    [set the channels up](https://bioconda.github.io/#usage). Naming them makes
    the command work on a machine that has not.

## What you get, and the first run

The package is `noarch: python` and contains **the pipeline and none of the
fourteen tools**. Its only dependencies are Python, Snakemake and its two
executor plugins, Textual, and `conda` itself.

That last one is why a pixi install works without a conda of your own, and it
was worth checking rather than assuming: on a machine with **no conda on
`PATH` at all**, both pixi routes ran a dry run to completion (measured
2026-09-04, macOS). A `pixi global install` exposes only `comparem2` and `cm2`,
but its trampoline puts the environment's own `bin` on `PATH` for the process,
and the bundled `conda` is in there.

The first run builds the two tool environments before any job starts, which
takes about a minute on a warm package cache and longer on a machine that has
to download them. `comparem2 --setup` does it deliberately, with no assemblies
and no databases needed, so the first real run does not pay for it.

## From git, for development

```bash
git clone https://github.com/cmkobel/CompareM2.git
cd CompareM2
pixi install
pixi run pytest        # 212 unit tests, no databases and no tools needed
pixi run cm2 --help
```

`pixi install` builds the pipeline's environment only — the fourteen tools are
deliberately not in it, so that there is exactly one way a tool can arrive and
it is the way a conda install uses too.

!!! note "Moving the directory invalidates the environment"
    Conda bakes the absolute prefix into shebangs and RPATHs, so after moving
    the checkout you need `rm -rf .pixi && pixi install`.

## The two tool environments

| Environment | Packages | DIAMOND | Contents | Deployed |
| --- | ---: | --- | --- | ---: |
| `main` | 568 | 2.2.5 | thirteen tools, plus curl and tar | 6.0 GB |
| `checkm2` | 127 | 2.1.11 | CheckM2 alone | 1.8 GB |

Together they are **7.7 GB** measured; the two figures are rounded, so they add
to 7.8. Eighteen rules point at those two. **CheckM2 is separate because it pins
DIAMOND 2.1.x while current Bakta needs 2.2.x** — they cannot co-solve, which
is also the answer to why the tools are not simply dependencies of the conda
package: no single environment can hold all fourteen.

They go to `~/.comparem2/envs`, shared across runs, moved with `--conda-prefix`
or `$COMPAREM2_CONDA_PREFIX`.

!!! warning "Keep that path stable"
    Snakemake identifies an environment by a hash that includes the *realpath*
    of the prefix, so moving it rebuilds both — and re-fetches AMRFinder's
    database, which lives inside one of them.

    A **relative** `--conda-prefix` is the sharp edge: it resolves against the
    directory you typed the command in, so the same relative path from two
    directories is two prefixes and two 7.7 GB builds. The default is
    home-relative and safe.

### Building them up front

Snakemake builds a missing environment during DAG construction, *before the
first job* — so a first run is silent for as long as that takes, which on a
cluster is the point at which someone kills it. `--setup` asks for the same
work deliberately:

```bash
comparem2 --setup
```

It takes **no assemblies and needs no databases**, and measured 61.7 s into a
fresh prefix. Run against a prefix that already has them, it builds nothing and
returns in about two seconds. The only thing that has to match a later run is
`--conda-prefix`, however you set it.

Two caveats. `--setup` is one-time per *catalogue*, not per machine: changing
any tool's version pin changes the environment file, which changes the hash,
which rebuilds. And `--setup --until <subset>` is not a cheaper setup — the
environment is the whole thirteen-tool `main` whatever subset you name.

## Databases

**Downloaded automatically.** Each database is a step in the workflow, so it is
fetched once, skipped if already present, and re-fetched if a previous attempt
was interrupted. Nothing needs to be placed by hand.

| Database | Download | On disk | Needed by |
|---|---:|---:|---|
| GTDB r232 | **60.8 GB** | 94 GB | `gtdbtk` |
| CheckM2 | 1.7 GB | 2.9 GB | `checkm2` |
| Bakta light | 1.3 GB *(documented)* | 4.0 GB | `bakta` |
| AMRFinder | unmeasured | — | `amrfinder` |

Downloads are measured from `content-length`; on-disk figures are measured
after extraction. The **62.5 GB** CompareM2 prints is GTDB plus CheckM2 and
nothing else — Bakta's 1.3 GB is Bakta's own documented figure and AMRFinder
publishes none, so both are counted as "of unknown size" rather than folded into
a total that would then look measured.

**Plan volumes around 101 GB, not 62.5** — extraction inflates GTDB from
60.8 GB to 94 GB.

They go to `~/.comparem2/databases` unless you pass `-d` or set
`$COMPAREM2_DATABASES`. Before running, CompareM2 prints what is actually
missing — not a total that includes what you already have:

```
2 assemblies, 5 tools
to download: amrfinder (1 of unknown size) -> /home/you/.comparem2/databases
```

!!! warning "AMRFinder ignores `-d`"
    `amrfinder -u` refuses the `--database` option — *"only operates on the
    default database directory"* — so its data lands inside the deployed
    environment rather than under `--databases`. CompareM2 records that the
    update ran, and that record lives with the run, so a rebuilt environment no
    longer leaves a marker claiming data that is gone. Rebuilding the
    environment does mean fetching it again: 241 MB, timed at 26 and 27 s in
    two measurements.

**GTDB-Tk is 97% of the measured download total**, and the release matters: 2.7
accepts only r232 and refuses r226, which is also why the figure is 60.8 GB
rather than the 141.4 GB it was until r232 replaced FastANI's reference genomes
with skani sketches. It stays in the default path because it is the
authoritative answer to "what is this genome", but if you do not need taxonomy,
leaving it out is the single biggest saving available:

```bash
comparem2 *.fna --until seqkit checkm2 bakta amrfinder mlst mashtree \
                        treecluster skani panaroo snp-dists fasttree carveme \
                        biosynthesis
```

Bakta uses the **light** database (1.3 GB / 4.0 GB on disk) rather than the
full one (30 GB / 84 GB on disk). That saves 29 GB for less specific functional
annotation, which a wide view can absorb — but note the Bakta paper's
annotation-quality figures are measured on the full database, a 53 GB version of
it at the time, and do not transfer.

## HPC

Snakemake's SLURM executor plugin ships with the package. Point CompareM2 at
more cores and let Snakemake submit:

```bash
comparem2 *.fna --cores 64
```

Job submission is Snakemake's business, not CompareM2's — the generated
Snakefile lives at `<output>/.comparem2/Snakefile` and takes any Snakemake
profile you already use.

Two things worth setting once on a cluster, because home directories have
quotas and neither default belongs there:

```bash
export COMPAREM2_DATABASES=/scratch/you/cm2-databases
export COMPAREM2_CONDA_PREFIX=/scratch/you/cm2-envs
comparem2 --setup            # build the environments on the login node
```

## Verification status

A clean install says nothing about whether the pipeline runs. Both Bakta and
Panaroo once resolved to years-old builds that installed cleanly and crashed on
first use, which is why every tool carries a minimum version.

`STATUS.md` in the repository tracks which command lines have actually been
**executed** on real genomes. All fourteen have, under conda deployment, in one
31-step run. GTDB-Tk was the most instructive: its rule carried six defects that
only running it could show, including a database release the installed tool
refuses.
