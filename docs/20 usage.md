# Usage

```bash
pixi run cm2 <assemblies>... [options]
```

Assemblies are passed as paths, not as a glob string. In v2 this was
`--config input_genomes="*.fna"`; now the shell expands it:

```bash
pixi run cm2 genomes/*.fna
```

## Options

| Option | Default | What it does |
|---|---|---|
| `-o`, `--output` | `results_comparem2` | output directory |
| `-d`, `--databases` | `databases` | where databases live |
| `-t`, `--cores` | `4` | cores for Snakemake |
| `--until TOOL...` | *(all)* | run only these tools and their dependencies |
| `--set TOOL--FLAG=VALUE` | — | override a tool argument; repeatable |
| `--tui` | off | interactive keyboard interface |
| `--isolated-launcher CMD` | — | how to enter an isolated tool's environment |
| `--keep-going` | off | keep running independent tools after a failure |
| `--dry-run` | off | show what would run |
| `--report-only` | off | re-render the report from existing outputs |

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
pixi run cm2 *.fna --until fasttree     # runs bakta, panaroo, fasttree
pixi run cm2 *.fna --until seqkit skani # runs just those two
```

There are no fixed pseudo-targets like v2's `fast` and `isolate`. The
dependency closure replaces them: name what you want and the prerequisites
follow. Some useful combinations:

```bash
# Fast, no databases at all
--until seqkit mashtree treecluster skani

# Everything except the 141 GB GTDB download
--until seqkit checkm2 bakta amrfinder mlst mashtree treecluster skani \
        panaroo snp-dists fasttree carveme
```

## Passthrough parameters

Any argument can be forwarded to any tool:

```bash
pixi run cm2 *.fna \
  --set treecluster--threshold=0.1 \
  --set skani-c=125 \
  --set bakta--gram=+
```

Naming one flag replaces only that flag — the tool's other defaults stay. A flag
with no value is passed bare: `--set bakta--force=`.

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
pixi run cm2 *.fna --tui
```

A keyboard interface over the same run: per-tool progress, the exact command
each step is running, and failures as they happen. It drives Snakemake through
its logger plugin system rather than scraping stdout, so the events are
structured.

## Re-rendering the report

The report is regenerated on every run, but you can rebuild it alone — useful
after a partial run, or when only the report code changed:

```bash
pixi run cm2 *.fna --report-only
```

Sections appear only when their outputs exist, so a partial run still gives a
readable document.

## Output layout

```
results_comparem2/
├── report.html                    the product; self-contained
├── samples/<name>/<name>.fna      canonical link to your input
├── samples/<name>/<tool>/…        per-genome results
├── <tool>/…                       whole-set results
└── .comparem2/
    ├── Snakefile                  generated from the tool specs
    ├── envs/                      generated conda envs for isolated tools
    └── log/                       one log per step
```

The generated `Snakefile` is a normal Snakemake workflow. If something fails,
that file plus the matching log in `.comparem2/log/` is where to look — and you
can run Snakemake against it directly with any profile you already use.
