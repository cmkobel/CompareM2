# Installation

CompareM2 is on Bioconda, so **pixi** or **conda** installs it. However it
arrives, the package is the pipeline alone — the fourteen analysis tools are
not in it, and Snakemake deploys them into six conda environments the first
time they are needed. There is no flag for whether to do that; it always
happens, into `--conda-prefix`.

## Requirements

  - **Linux.** The analysis tools are `linux-64` only. On macOS you can run the
    unit tests and render reports, but not the pipeline.
  - **[pixi](https://pixi.prefix.dev/latest/#installation)** or **conda**, to
    install the package. You do not need a conda of your own: `conda` is a
    dependency *of* the package, because Snakemake shells out to it to build
    the tool environments, so it arrives with CompareM2 either way. Checked on
    a machine with no conda on `PATH` at all — both pixi routes ran a dry run
    to completion (2026-09-04).
  - **Disk.** **1.4 GB** of tool environments, plus databases — see below.
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

## The six tool environments

The package is `noarch: python` and depends only on Python, Snakemake and its
two executor plugins, Textual, and `conda`. Everything the analyses need is
deployed on first use into six environments, grouped by dependency ecosystem:

| Environment | Tools | On disk |
| --- | --- | ---: |
| `basic` | seqkit, skani, snp-dists, fasttree, treecluster, curl, tar | 55 MB |
| `perl` | mashtree, mlst, panaroo | 543 MB |
| `annotation` | bakta, amrfinder | 83 MB |
| `gtdbtk` | gtdbtk | 62 MB |
| `carveme` | carveme, biosynthesis | 102 MB |
| `checkm2` | checkm2 | 645 MB |

**1.4 GB for all six**, measured on GenomeDK 2026-09-07. Conda hardlinks
packages that several environments share, so the total is well under the sum of
independent installs — and correspondingly, the per-environment figures do not
add up to it.

**Why six and not one.** CheckM2 pins DIAMOND 2.1.x while current Bakta needs
2.2.x, so no single environment ever held all fourteen. The rest of the split
is newer: in September 2026 the thirteen-tool environment stopped solving
because the Perl stack broke upstream and took every other tool with it. Tools
that share an environment share its fate, so they are grouped by what they
actually depend on. An environment per *tool* would be the opposite mistake.

They go to `~/.comparem2/envs`, shared across runs, moved with `--conda-prefix`
or `$COMPAREM2_CONDA_PREFIX`.

!!! warning "Keep that path stable"
    Snakemake identifies an environment by a hash that includes the *realpath*
    of the prefix, so moving it rebuilds all six — and re-fetches AMRFinder's
    database, which lives inside `annotation`.

    A **relative** `--conda-prefix` is the sharp edge: it resolves against the
    directory you typed the command in, so the same relative path from two
    directories is two prefixes and two full builds. The default is
    home-relative and safe.

### Building them up front

Snakemake builds a missing environment during DAG construction, *before the
first job* — so a first run is silent for a minute or more, which on a cluster
is the point at which someone kills it. `--setup` asks for the same work
deliberately:

```bash
comparem2 --setup
```

It takes **no assemblies and needs no databases**, and measured 61.7 s into a
fresh prefix; against a prefix that already has them it builds nothing and
returns in about two seconds. The only thing that has to match a later run is
`--conda-prefix`, however you set it.

Two caveats. `--setup` is one-time per *catalogue*, not per machine: changing
any tool's version pin changes the environment file, which changes the hash,
which rebuilds. And `--setup --until <subset>` builds only the environments
that subset needs — `--until seqkit skani` builds `basic` alone.

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

The **62.5 GB** CompareM2 prints is GTDB plus CheckM2 and nothing else:
downloads are measured from `content-length`, Bakta's 1.3 GB is Bakta's own
documented figure and AMRFinder publishes none, so both are counted as "of
unknown size" rather than folded into a total that would then look measured.
**Plan volumes around 101 GB, not 62.5** — extraction inflates GTDB from
60.8 GB to 94 GB.

Databases go to `~/.comparem2/databases` unless you pass `-d` or set
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
    environment does mean fetching it again: 241 MB, timed at 26 and 27 s.

### The two choices behind those numbers

**GTDB is 97% of the download**, and its release matters as much as its size:
GTDB-Tk 2.7 accepts only r232 and refuses r226 — which is also why the figure
is 60.8 GB rather than the 141.4 GB it was until r232 replaced FastANI's
reference genomes with skani sketches. It stays in the default path because it
is the authoritative answer to "what is this genome", but if you do not need
taxonomy, leaving it out is the single biggest saving available — see
[running a subset](20 usage.md#running-a-subset).

**Bakta uses the light database** (1.3 GB / 4.0 GB on disk) rather than the
full one (30 GB / 84 GB on disk). That saves 29 GB for less specific functional
annotation, which a wide view can absorb — but note the Bakta paper's
annotation-quality figures are measured on the full database, a 53 GB version of
it at the time, and do not transfer.

## HPC

Three things, in this order.

**1. Put the databases and the environments somewhere with room.** Home
directories have quotas and neither default belongs there — and both paths must
be visible from the compute nodes, not just the login node.

```bash
export COMPAREM2_DATABASES=/scratch/you/comparem2-databases
export COMPAREM2_CONDA_PREFIX=/scratch/you/comparem2-envs
comparem2 --setup            # build the environments on the login node
```

Run `--setup` on the login node before submitting anything, so the first job to
start does not build all six environments inside its own allocation. It also has
to be the login node on a cluster whose compute nodes have no outbound network,
which is common — CompareM2 has no way to fetch a conda package from a node
that cannot reach `conda.anaconda.org`.

Database downloads are handled for you: the four `download_*` rules are
declared `localrules`, so they run wherever Snakemake is — the login node —
rather than being submitted. That is deliberate, and it is what stops the
60.8 GB GTDB fetch from being sent to a node with no route to the internet.

**2. Write a profile.** Job submission is Snakemake's, not CompareM2's:
a profile directory holding a `config.yaml` names the executor and carries the
account, the partition and the resource defaults. CompareM2 ships none — your
cluster's account name and partitions are not ours to guess.

```yaml
# ~/.config/snakemake/slurm/config.yaml
executor: slurm
jobs: 200                      # max jobs in the queue at once
local-cores: 4                 # for rules that stay on the login node
latency-wait: 30               # shared filesystems are not instantaneous

default-resources:
  slurm_account: YOUR_ACCOUNT
  slurm_partition: normal
  mem_mb: 8000
  runtime: "4h"                # quote it — see the warning below

# CompareM2's rules declare threads but not memory or walltime, so the
# defaults above apply to all of them. Override the greedy one by name:
set-resources:
  gtdbtk:
    mem_mb: 64000
    runtime: "12h"
```

!!! danger "Quote `runtime`, and give it a unit"
    A bare number here is read as **seconds**, and `runtime` is in minutes — so
    `runtime: 240` asks for four minutes, not four hours, and every job dies at
    the walltime. Verified against Snakemake 9.16.3 on GenomeDK: `runtime: 60`
    parses to `Resource("runtime", 1)`, while `runtime: "12h"` parses to 720.

**Those numbers are starting points, not measurements.** GTDB-Tk's follows the
GTDB-Tk 2 paper's under-55 GB figure for divide-and-conquer placement, with
headroom; nothing else here has been measured on a cluster at all. Run
`sacct -o JobName,MaxRSS,Elapsed -j <jobid>` afterwards and correct them —
`MaxRSS` is the number that decides whether the next run is killed.

**3. Point CompareM2 at it.**

```bash
comparem2 *.fna --profile ~/.config/snakemake/slurm
```

A bare name works too — `--profile slurm` searches `~/.config/snakemake`, the
way Snakemake's own `--profile` does. That directory is the default place for a
profile because it is Snakemake's default place, not ours:
`appdirs.AppDirs("snakemake", "snakemake")`, so `$XDG_CONFIG_HOME/snakemake` on
Linux, with `/etc/xdg/snakemake` searched too for a profile your sysadmin
installed for everyone. The flag is a passthrough, so anything Snakemake
supports works: `snakemake-executor-plugin-slurm` for SLURM and
`snakemake-executor-plugin-cluster-generic` for PBS, SGE and LSF are both
already installed.

**Or set it once and forget the flag.** `$SNAKEMAKE_PROFILE` is Snakemake's own
variable, and CompareM2 follows it, so this belongs next to the two exports in
step 1:

```bash
export SNAKEMAKE_PROFILE=slurm       # a name under ~/.config/snakemake, or a path
comparem2 *.fna                      # submits, and says so
```

Every run from that shell then goes to the queue. Because that is easy to
forget, a run that a variable submitted prints where it is going —
`execution: /home/you/.config/snakemake/slurm ($SNAKEMAKE_PROFILE)` — and the
TUI's header says the same in its `execution` row.

`--profile` beats the variable, and **`--profile none` runs locally in spite of
it**, which is what you want for a four-genome check on the login node.

`--tui` works with a profile: the interface runs on the login node and shows
jobs starting and finishing as the queue runs them.

!!! note "`--setup` and `--unlock` ignore the profile on purpose"
    Both are login-node bookkeeping — solving environments, clearing a lock —
    so they run locally whether or not `$SNAKEMAKE_PROFILE` is exported.

!!! warning "`--cores` is not a submission setting"
    `comparem2 *.fna --cores 64` starts 64 processes **on the machine you typed
    it on**. It submits nothing. On a login node that is a way to get an email
    from your sysadmin. Queue submission needs a profile.

    Under a profile, `-t/--cores` is left to the profile unless you pass it,
    because a number on the command line overrides the profile's own.

## From git, for development

```bash
git clone https://github.com/cmkobel/CompareM2.git
cd CompareM2
pixi install
pixi run pytest        # 269 unit tests, no databases and no tools needed
pixi run comparem2 --help
```

`pixi install` builds the pipeline's environment only — the fourteen tools are
deliberately not in it, so that there is exactly one way a tool can arrive and
it is the way a conda install uses too.

!!! note "Moving the directory invalidates the environment"
    Conda bakes the absolute prefix into shebangs and RPATHs, so after moving
    the checkout you need `rm -rf .pixi && pixi install`.
