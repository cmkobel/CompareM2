# Usage

```bash
comparem2 <assemblies>... [options]      # installed with conda
pixi run cm2 <assemblies>... [options]   # from a git checkout
```

`cm2` is a second entry point for the same program, so `comparem2` and `cm2`
are interchangeable. The examples below use the plain form; from a git checkout
prefix them with `pixi run`.

Assemblies are passed as paths, and the shell expands the glob:

```bash
cm2 genomes/*.fna
```

Relative paths mean what they look like they mean, from any directory. A pixi
task runs from the workspace root rather than from your shell's directory, so
`pixi run cm2 *.fna` in a subdirectory would otherwise look for the files
somewhere else entirely; CompareM2 resolves inputs, `--output` and
`--databases` against `$INIT_CWD`, which pixi sets to where the command was
typed. Results land next to the genomes.

## Options

| Option | Default | What it does |
|---|---|---|
| `-o`, `--output` | `results_comparem2` | output directory |
| `-d`, `--databases` | `~/.comparem2/databases` | where databases live |
| `-t`, `--cores` | `4` | cores for Snakemake |
| `--until TOOL...` | *(all)* | run only these tools and their dependencies |
| `--set TOOL-FLAG=VALUE` | — | override a tool argument; repeatable |
| `--tui` | off | interactive keyboard interface |
| `--conda-prefix DIR` | `~/.comparem2/envs` | where the tool environments go |
| `--setup` | off | build those environments and exit; takes no assemblies |
| `--demo` | off | run the bundled plasmids; takes no assemblies |
| `--keep-going` | off | keep running independent tools after a failure |
| `--dry-run` | off | show what would run |
| `--report-only` | off | re-render the report from existing outputs |
| `--unlock` | off | release a stale lock left by a killed run |
| `--version` | | print the version and exit |

There is no flag for *whether* to deploy the tools. Snakemake always does, into
`--conda-prefix` — see [Installation](10 installation.md).

## The bundled demo

```bash
cm2 --demo
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

## Where databases go

Databases are shared across runs, not stored per-run: the default is
`~/.comparem2/databases`, the same location whatever directory you invoke from
and outside any checkout, so deleting a checkout does not cost a re-download.
The full set is 62.5 GB of measured downloads plus two of unmeasured size,
and 60.8 GB of that is GTDB alone.

A home directory is the wrong place for that on a cluster with a home quota, so
set the location once:

```bash
export COMPAREM2_DATABASES=/evo/postdoc/cm2-databases
```

Precedence is `-d`, then `$COMPAREM2_DATABASES`, then `~/.comparem2/databases`.
Whichever wins is printed before anything is fetched:

```
to download: checkm2, gtdb, bakta-light, amrfinder (62.5 GB + 2 of unknown size) -> /evo/postdoc/cm2-databases
```

Only what is actually missing is listed, so this shrinks to nothing once the
databases are in place.

Two databases are not under this root, and cannot be:

- **AMRFinder** rejects `-d` on update (`amrfinder -u -d <dir>` exits with *"only
  operates on the default database directory"*), so its data lands in
  `$CONDA_PREFIX` and only a marker file is recorded here.
- **GTDB-Tk** has no flag for its database at all; it is passed
  `GTDBTK_DATA_PATH=<root>/gtdb` instead.

## Sample names

Every input is linked to `<output>/samples/<name>/<name>.fna`, and that is what
every tool reads. The name comes from the filename stem with anything outside
`[A-Za-z0-9._-]` replaced by `_`, because a space in a filename otherwise
produces a silently broken workflow rule. CompareM2 tells you when it renames:

```
note: '116_2 duplicate.fna' -> sample '116_2_duplicate'
```

Two inputs that reduce to the same name is an error, not a silent overwrite.

## Running a subset

`--until` takes tool names and pulls in whatever they need:

```bash
cm2 *.fna --until fasttree     # runs bakta, panaroo, fasttree
cm2 *.fna --until seqkit skani # runs just those two
```

There are no fixed presets to memorise: name what you want and the
prerequisites follow. Some useful combinations:

```bash
# Fast, no databases at all
--until seqkit mashtree treecluster skani

# Everything except the 60.8 GB GTDB download
--until seqkit checkm2 bakta amrfinder mlst mashtree treecluster skani \
        panaroo snp-dists fasttree carveme biosynthesis
```

## Passthrough parameters

Any argument can be forwarded to any tool:

```bash
cm2 *.fna \
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
cm2 *.fna --tui
```

A keyboard interface over the same run: per-tool progress, the download size
before anything is fetched, and failures as they happen. It drives Snakemake
through its logger plugin system rather than scraping stdout, so the events are
structured.

`space` selects and deselects the tool under the cursor, `a` and `n` select all
and none, `r` runs, `q` quits. `▣` is chosen, `▨` is pulled in as a dependency
of something chosen, `▢` is off.

Every other flag works the same way with `--tui` as without it — `--until` seeds
the selection, and `--set`, `--keep-going` and `-d` are all honoured:

```bash
cm2 *.fna --tui --until mashtree treecluster
```

Without `--until` everything is selected, which includes GTDB-Tk and its 60.8
GB. Read the size on the second line before pressing `r`.

`--dry-run` is refused with `--tui`, because the tool list is already the dry
run and it shows the download size too.

## Re-rendering the report

The report is regenerated on every run, but you can rebuild it alone — useful
after a partial run, or when only the report code changed:

```bash
cm2 *.fna --report-only
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
cm2 *.fna --unlock
cm2 *.fna            # picks up where it stopped
```

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
    ├── envs/                       the two generated conda env files
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
