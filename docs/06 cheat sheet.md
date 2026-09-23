# CompareM2 cheat sheet

Linux only. `cm2` is an alias for `comparem2` everywhere. In a git checkout or
pixi workspace, read every command below as `pixi run comparem2 …`.

## 1. Point at the genomes

```bash
cd /path/to/genomes        # relative paths resolve from here, even under pixi run
ls *.fna                   # the shell expands the glob; this is exactly what gets passed
```

Name no assemblies at all and CompareM2 uses the `*.fna`, `*.fa`, `*.fasta`
and `*.fas` files in this directory, and prints which directory that was.

## 2. Run

```bash
comparem2                          # the *.fna, *.fa, *.fasta and *.fas files in this directory
comparem2 *.fna                    # all 14 tools -> ./results_comparem2/report.html
comparem2 *.fna -o run1            # output directory
comparem2 --demo                   # 6 bundled plasmids, no genomes or databases needed
comparem2 *.fna --dry-run          # what would run
comparem2 *.fna --tui              # keyboard interface; nothing selected until you press space/a then r
```

Before fetching anything it prints the cost:

```
4 assemblies, 14 tools
to download: checkm2, gtdb, bakta-light, amrfinder (62.5 GB + 2 of unknown size) -> ~/.comparem2/databases
tool environments: 6 in ~/.comparem2/envs (none built yet)
```

The first run is quiet for about a minute while those six conda environments
solve. That is normal.

## 3. Run a subset: the big saving

`--until` names tools and pulls in their prerequisites.

```bash
comparem2 *.fna --until fasttree               # runs bakta, panaroo, fasttree
comparem2 *.fna --until seqkit mashtree treecluster skani   # fast, zero databases

# everything except the 60.8 GB GTDB download (~3 GB total instead of 62.5 GB)
comparem2 *.fna --until seqkit checkm2 bakta amrfinder mlst mashtree treecluster \
                       skani panaroo snp-dists fasttree carveme biosynthesis
```

Tool names: `seqkit checkm2 gtdbtk bakta amrfinder mlst mashtree treecluster
skani panaroo snp-dists fasttree carveme biosynthesis`.

Dependencies: `amrfinder`, `panaroo`, `carveme` ← `bakta`; `snp-dists`,
`fasttree` ← `panaroo`; `treecluster` ← `mashtree`; `biosynthesis` ← `carveme`.

Dependencies are selected automatically.


## 4. Output

```
results_comparem2/
├── report.html                    self-contained; survives being emailed off a cluster
├── <tool>/…                       whole-set results
├── samples/<name>/<tool>/…        per-genome results
├── logs/… and samples/<name>/logs/…
└── .comparem2/Snakefile           the generated workflow, if you need to debug
```

Sample names come from the filename stem with anything outside `[A-Za-z0-9._-]`
replaced by `_`, and the rename is printed.

```bash
comparem2 *.fna --report-only      # re-render the report from what is already on disk
comparem2 *.fna --keep-going       # don't stop the other tools when one fails
```


## 5. `--set`: forwarding arguments to a tool

`--set <tool><flag>=<value>`, repeatable. Spell the flag exactly as the tool
does, dashes and all: TreeCluster's long options take two, skani's `-c` takes
one. Naming one flag replaces only that flag; the tool's other defaults stay.

```bash
# the two worth knowing
comparem2 *.fna --set skani-c=125                  # complete isolates; default 70 suits fragmented MAGs
comparem2 *.fna --set treecluster--threshold=0.1   # default 0.05; the threshold dominates the clustering

# several at once
comparem2 *.fna \
  --set treecluster--method=avg_clade \
  --set treecluster--threshold=0.1 \
  --set skani-c=125 \
  --set mashtree--genomesize=3000000 \
  --set bakta--gram=+

# a flag with no value is passed bare
comparem2 *.fna --set bakta--force=
```

Defaults you are overriding: `bakta --force`; `mashtree --genomesize 5000000
--mindepth 5 --kmerlength 21 --sketch-size 10000`; `treecluster --method
max_clade --threshold 0.05`; `skani -c 70`. Every tool's full command line is in
[what analyses does it do](30 what analyses does it do.md), generated from the
specs so it cannot drift from what runs.

The full CLI, the TUI, cluster profiles and what to do after a run is killed are
in [usage](20 usage.md).
