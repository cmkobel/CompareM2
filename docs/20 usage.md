# Usage

```bash
comparem2 <assemblies>... [options]            # installed with conda
pixi run comparem2 <assemblies>... [options]   # from a git checkout
```

The examples below use the first form; from a git checkout, read each one as
`pixi run comparem2 …`. A conda install also puts `cm2` on your `PATH` as a
shorter alias, interchangeable with `comparem2` everywhere.

Assemblies are passed as paths, and the shell expands the glob:

```bash
comparem2 genomes/*.fna
```

Relative paths mean what they look like they mean, from any directory —
including under `pixi run`, which would otherwise resolve them against the
workspace root rather than your shell's directory. Inputs, `--output` and
`--databases` are all resolved against where you typed the command, so results
land next to the genomes.

## Options

| Option | Default | What it does |
|---|---|---|
| `-o`, `--output` | `results_comparem2` | output directory |
| `-d`, `--databases` | `~/.comparem2/databases` | where databases live |
| `-t`, `--cores` | `4` | cores for Snakemake; with `--profile`, left to the profile unless given |
| `--profile DIR` | — | Snakemake profile, for submitting to a cluster queue |
| `--until TOOL...` | *(all)* | run only these tools and their dependencies |
| `--set TOOL-FLAG=VALUE` | — | override a tool argument; repeatable |
| `--tui` | off | interactive keyboard interface |
| `--conda-prefix DIR` | `~/.comparem2/envs` | where the tool environments go |
| `--setup` | off | build those environments and exit; takes no assemblies |
| `--demo` | off | run the bundled plasmids; takes no assemblies |
| `--keep-going` | off | keep running independent tools after a failure |
| `--dry-run` | off | show what would run |
| `--report-only` | off | re-render the report from existing outputs |
| `--unlock` | off | release a stale lock on `--output` and exit; takes no assemblies |
| `--version` | | print the version and exit |

There is no flag for *whether* to deploy the tools. Snakemake always does, into
`--conda-prefix` — see [Installation](10 installation.md).

`--profile` is a passthrough to Snakemake and is what makes jobs go to a queue
rather than to this machine; `--cores` never submits anything. See
[HPC](10 installation.md#hpc) for a worked profile.

## Running a subset

`--until` takes tool names and pulls in whatever they need:

```bash
comparem2 *.fna --until fasttree     # runs bakta, panaroo, fasttree
comparem2 *.fna --until seqkit skani # runs just those two
```

There are no fixed presets to memorise: name what you want and the
prerequisites follow. Two combinations worth knowing:

```bash
# Fast, and no databases at all
--until seqkit mashtree treecluster skani

# Everything except the 60.8 GB GTDB download
--until seqkit checkm2 bakta amrfinder mlst mashtree treecluster skani \
        panaroo snp-dists fasttree carveme biosynthesis
```

The second is the single biggest saving available — GTDB-Tk is 60.8 GB of the
download, and everything else together is roughly 3 GB.

## Passthrough parameters

Any argument can be forwarded to any tool:

```bash
comparem2 *.fna \
  --set treecluster--threshold=0.1 \
  --set skani-c=125 \
  --set bakta--gram=+
```

Write the flag exactly as the tool spells it, dashes and all — that is why
`treecluster--threshold` has two and skani's `-c` has one. Naming one flag
replaces only that flag; the tool's other defaults stay. A flag with no value is
passed bare: `--set bakta--force=`.

Every tool's defaults are listed on
[what analyses does it do](30 what analyses does it do.md), generated from the
specs so they cannot drift from what actually runs.

!!! note "Two things worth overriding"
    `--set skani-c=125` if all your genomes are complete isolates — the default
    of 70 is the more accurate setting for fragmented MAGs but costs runtime.

    `--set treecluster--threshold=…` if the clusters look wrong. The threshold
    dominates the answer: TreeCluster's own paper moved from 181,574 clusters to
    10,112 by sweeping it, so re-run at a couple of nearby values before
    reporting a grouping.

## The TUI

```bash
comparem2 *.fna --tui
```

A keyboard interface over the same run: per-tool progress, the download size
before anything is fetched, and failures as they happen. It drives Snakemake
through its logger plugin system rather than scraping stdout, so the events are
structured.

`space` selects and deselects the tool under the cursor, `a` and `n` select all
and none, `r` runs, `u` releases a lock, `q` quits. `▣` is chosen, `▨` is pulled
in as a dependency of something chosen, `▢` is off.

**Nothing is selected when it opens.** Choosing the analyses is what the
interface is for, and selecting all fourteen by default put GTDB-Tk's 60.8 GB
download one keypress from someone who had not read the table yet. `a` is still
one key away if you do want everything. `--until` seeds the selection, so
`--tui --until mashtree treecluster` opens with exactly those two chosen.

### What it says before you press anything

**Which analyses already ran here.** Open it on a directory you have run in
before and the status column is filled in from what is on disk: `already run`
where every declared output is present, `part-finished` where some are missing.
Those are read from the same files Snakemake decides resumability on, so the
table is what a re-run would actually skip — a `part-finished` tool is one that
will be redone. A tool you have not selected still reports what it has; whether
it is selected is what the mark column says.

**Where the run's four locations come from.** Above the table:

```
output    /faststorage/project/x/run/results_comparem2  given
databases /faststorage/project/x/comparem2_databases    $COMPAREM2_DATABASES
tool envs /home/carl/.comparem2/envs                    default
execution local                                         default
```

The right-hand column is the *origin*, not the value: `default`, `given` for
something you typed, or the name of the environment variable the value came
from. `given, overriding $COMPAREM2_DATABASES` means the variable is set and a
`-d` beat it — which is worth seeing, because a variable exported in `.bashrc`
months ago and silently overridden looks exactly like no variable at all. Both
of these decide whether existing work gets re-used: a databases directory that
is not the one holding your 62.5 GB re-downloads it, and a moved
`$COMPAREM2_CONDA_PREFIX` re-solves every tool environment.

**Whether the output directory is locked**, and `u` clears it. A run that was
killed leaves a lock Snakemake refuses to start on — see [After a run is
killed](#after-a-run-is-killed) for what that is. `r` on a locked directory
refuses to start rather than letting the run fail several seconds in, and `u`
does the same thing `--unlock` does without leaving the interface.

It asks first, and the question is not "are you sure" but "is anything else
running": a lock file lists paths and carries no process id, so neither you nor
CompareM2 can tell a dead run's lock from a live one's, and clearing a live
one's puts two Snakemake processes on the same outputs. Check before saying
yes.

### What it says while it runs

Above the progress bar, one line answers "is this still going, or has it
frozen":

```
⠹ bakta, gtdbtk · 4m 12s
```

A spinner, what is running now, and how long since you pressed `r`. When no job
is running the line says why — `starting up — the DAG, and tool environments on
a first run` is the long one, because a first run solves six conda environments
before anything else happens and Snakemake's own output for that is quietened
under the interface. The progress bar has no total to draw until Snakemake
reports one, so it pulses rather than sitting at 0%.

The animation is driven by the interface's own event loop, which is the point:
if it is genuinely wedged, the spinner stops with it rather than reassuring you
it hasn't. When the run ends the line is replaced by `not running — the last run
took 4m 12s`, so a still screen never has to be interpreted.

Every other flag works the same way with `--tui` as without it — `--until` seeds
the selection, and `--set`, `--keep-going` and `-d` are all honoured:

```bash
comparem2 *.fna --tui --until mashtree treecluster
```

`--dry-run` is refused with `--tui`, because the tool list is already the dry
run and it shows the download size too.

## The bundled demo

```bash
comparem2 --demo
```

Six *Enterococcus faecium* plasmids ship inside the package — 461 KB, the only
non-Python file in it — so this needs no genomes of your own, no databases and
no network beyond the tool environments themselves. They are extracted to
`<output>/demo_assemblies/`, where you can look at them and delete them.

It runs `seqkit`, `mashtree`, `treecluster` and `skani`: the four analyses that
need no database. That list is fixed rather than defaulted, because the inputs
are **plasmids** — CheckM2 would report a completeness near zero, correctly and
uselessly, since it is looking for a chromosome's marker genes. Naming
`--until` yourself still overrides it, on the assumption that you have a reason.

A seventh input is the sixth one again as `116_2 duplicate.fna`. It costs
nothing to ship and it gives the report something to check itself against: the
pair must come out at 0.00000 mash distance and 100.00% ANI, and the space in
the filename exercises sample-name canonicalisation on the way.

## Sample names

Every input is linked to `<output>/samples/<name>/<name>.fna`, and that is what
every tool reads. The name comes from the filename stem with anything outside
`[A-Za-z0-9._-]` replaced by `_`, because a space in a filename otherwise
produces a silently broken workflow rule. CompareM2 tells you when it renames:

```
note: '116_2 duplicate.fna' -> sample '116_2_duplicate'
```

Two inputs that reduce to the same name is an error, not a silent overwrite.

## Where databases go

Databases are shared across runs, not stored per-run, so deleting a checkout
does not cost a re-download. Precedence is `-d`, then `$COMPAREM2_DATABASES`,
then `~/.comparem2/databases` — and a home directory is the wrong place for
101 GB on a cluster with a quota:

```bash
export COMPAREM2_DATABASES=/evo/postdoc/cm2-databases
```

Whichever location wins is printed before anything is fetched, listing only
what is actually missing:

```
to download: checkm2, gtdb, bakta-light, amrfinder (62.5 GB + 2 of unknown size) -> /evo/postdoc/cm2-databases
```

Two databases are not under this root, and cannot be:

- **AMRFinder** rejects `-d` on update (`amrfinder -u -d <dir>` exits with *"only
  operates on the default database directory"*), so its data lands in
  `$CONDA_PREFIX` and only a marker file is recorded here.
- **GTDB-Tk** has no flag for its database at all; it is passed
  `GTDBTK_DATA_PATH=<root>/gtdb` instead.

Sizes are in [Installation](10 installation.md#databases).

## Re-rendering the report

The report is regenerated on every run, but you can rebuild it alone — useful
after a partial run, or when only the report code changed:

```bash
comparem2 *.fna --report-only
```

Sections appear only when their outputs exist, so a partial run still gives a
readable document.

## After a run is killed

Snakemake locks the output directory, so a run that died without releasing it —
SIGKILL, a lost node, a power cut — leaves the next one refusing to start:

```
LockException: Directory cannot be locked.
```

Nothing is wrong with the results. Release the lock and carry on:

```bash
comparem2 --unlock         # add -o if the run wrote somewhere else
comparem2 *.fna            # picks up where it stopped
```

The lock belongs to the output directory, not to the assemblies, so `--unlock`
needs only `-o`. Naming the assemblies as well is accepted — adding the flag to
the command that just died is the obvious move — and they are not read.

Under `--tui` this is `u`, which asks before clearing. Either way, make sure no
other run is writing to that directory first: the lock is what stops two
Snakemake processes from corrupting each other's outputs, and nothing in it
says whether the process that made it is still alive.

Downloads resume rather than restart: a killed GTDB fetch continues its partial
tarball instead of fetching 60.8 GB again.

## Output layout

```
results_comparem2/
├── report.html                     the product; self-contained
├── <tool>/…                        whole-set results
├── logs/<tool>.log                 one per whole-set step
├── samples/<name>/<name>.fna       canonical link to your input
├── samples/<name>/<tool>/…         per-genome results
├── samples/<name>/logs/<tool>.log  one per per-genome step
└── .comparem2/
    ├── Snakefile                   generated from the tool specs
    ├── envs/                       the six generated conda env files
    └── gtdbtk_batchfile.tsv        generated input: genome path, genome id
```

A log sits beside the results it describes, so a per-genome tool leaves one log
per genome rather than one for the run. Database downloads are the exception:
they log to `<databases>/logs/`, next to the data instead of next to the run —
all but AMRFinder's, which has no directory of its own under `--databases` and
writes to `logs/download_amrfinder.log` here.

The generated `Snakefile` is a normal Snakemake workflow. If something fails,
that file plus the matching log is where to look — and you can run Snakemake
against it directly with any profile you already use.
