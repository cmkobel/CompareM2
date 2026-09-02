# CompareM2 v3 — status

What is currently true of a real run. This file changes whenever something is
re-run, which is why it is not in [DESIGN.md](DESIGN.md) — decisions should not
need editing because a tool was verified again.

Last updated **2026-09-02**, from runs on thylakoid.

## Tool verification: 12 of 13

**Verified means executed** on the four *E. faecium* test genomes and producing
correct-looking output. It never means "installed" — see
[DECISIONS.md](DECISIONS.md#a-solve-is-not-a-working-environment).

| Tool | Verified | Evidence |
| ---- | -------- | -------- |
| seqkit | 09-02 | per-contig lengths and GC; identical stats for the duplicate pair |
| checkm2 | 09-02 | 100% complete all four, via `--isolated-launcher`; `--database_path` wants the `.dmnd` **file**, not its directory |
| bakta | 09-02 | 4/4 against db v6.0 light; db `software-min` 1.11 against bakta 1.12.1 |
| amrfinder | 09-02 | 7 and 11 genes; database fetched by the pipeline itself |
| mlst | 09-02 | ST32 / ST32 / ST117 / ST78 — duplicates agree |
| mashtree | 09-02 | duplicate pair at distance 0.00000 |
| treecluster | 09-02 | needed `--threshold`; v2's defaults adopted |
| skani | 09-02 | `-c 70` accepted; duplicates at 100.00% ANI; also writes an `.af` matrix |
| panaroo | 09-02 | 3,780 clusters, 2,091 core; core alignment 1,934,948 bp |
| snp-dists | 09-02 | 0 SNPs between the duplicate pair, identical gap counts |
| fasttree | 09-02 | at `threads=1`; duplicates at 0.0 branch length |
| carveme | 09-02 | 4/4 SBML models, ~4 MB each; ~9 min per genome |
| **gtdbtk** | **never** | needs the 141.4 GB download. Its command would also have failed with the database present, until `GTDBTK_DATA_PATH` was wired up |

### The standing cross-check
The test set contains `116_2.fna` and `116_2 duplicate.fna` — the same genome
twice, and a filename with a space. **Any tool that treats them differently is
wrong.** Every verified tool agrees: 0.00000 mash distance, 100.00% ANI, 0 SNPs,
identical contig statistics, identical CheckM2 values, model sizes within 30
bytes (the sample name is embedded in the SBML).

### Divergence sanity
A report whose numbers cannot be checked is not a result. 6,190 SNPs between
116_2 and E8202 over a 1,934,948 bp core alignment is 0.320% — the right order
for two strains of one species.

## Also verified

- `pixi install` — both environments, from a manifest that had never been
  installed
- **clean clone** — cloned fresh, both environments built, tests pass,
  `test-fast` 8/8
- `--isolated-launcher` — CheckM2 from its own environment
- Automatic database download — `download_amrfinder` fetched database
  `2026-08-07.1` in 26 s as a workflow step
- Snakemake lock now lives in the output directory
- **`--tui`, end to end** — `--until mashtree treecluster --tui` under a real
  terminal: two tools selected from the command line, both `done`, report
  written, clean exit. And headlessly through `run_test()` with `seqkit` and
  `checkm2`, so the isolated launcher goes through the TUI path too. Five
  defects had to be fixed first, see [DECISIONS.md](DECISIONS.md).
- 100 unit tests, ~1.7 s
- CI green on GitHub, 18–22 s per run
- `mkdocs build --strict`
- `docs/generate.py --check` — generated pages current

## Environments

Re-solved after sylph was removed. Measured, `linux-64`:

| Environment | Packages | DIAMOND | Contents |
| ----------- | -------: | ------- | -------- |
| `default` | 568 | 2.2.5 | bakta 1.12.1 and eleven other tools |
| `checkm2` | 127 | 2.1.11 | checkm2 1.1.0 alone |

The two DIAMOND versions are the conflict that forced the split, so this is the
isolation working rather than an accident. Both minimum-version pins held —
bakta 1.12.1 and panaroo 1.8.0, not the 1.8.1 and 1.1.2 builds that install
cleanly and crash.

Earlier figure, kept for comparison: 1,142 packages / 1.58 GB, measured with
sylph still in the set.

## Databases

| Database | Download | On disk | Note |
| -------- | -------: | ------: | ---- |
| GTDB r226 | 141.4 GB | — | measured; 91% of the total, never downloaded |
| CheckM2 | 1.7 GB | 2.9 GB | measured |
| Bakta light | 1.3 GB | 4.0 GB | download figure is Bakta's documented one; on-disk measured |
| AMRFinder | unmeasured | — | version `2026-08-07.1`, 26 s to fetch. Lands in `$CONDA_PREFIX`, **not** under `--databases` |

Software is 1.58 GB against 143 GB of data, and GTDB-Tk is essentially all of it.

## Where things live on thylakoid

24 cores, 125 GB RAM, 931 GB free on `/evo`.

| | |
| --- | --- |
| checkout | `/evo/postdoc/comparem2` — a real clone of branch `v3` |
| databases | `/evo/postdoc/cm2-databases` — 6.9 GB, deliberately **outside** any checkout so deleting a checkout does not cost a re-download |
| pixi | `/home/thylakoid/.pixi/bin/pixi` |

```bash
cd /evo/postdoc/comparem2
export PATH=$HOME/.pixi/bin:$PATH

pixi run pytest          # 92 tests, no tools or databases needed
pixi run test-fast       # 4 genomes, no databases needed

pixi run cm2 my/*.fna \
  --output results_myrun \
  -d /evo/postdoc/cm2-databases \
  --cores 24 \
  --isolated-launcher "$HOME/.pixi/bin/pixi run -e {tool}"
```

**Two flags are not optional there.** Without `-d`, the 6.9 GB is downloaded
again. Without an **absolute** launcher path, CheckM2 fails with
`command not found`, because Snakemake's shell does not inherit an interactive
PATH.

To skip the 141 GB:

```bash
--until seqkit checkm2 bakta amrfinder mlst mashtree treecluster skani \
        panaroo snp-dists fasttree carveme
```

## Known broken or unfinished

- **gtdbtk has never run.** Needs the download, and `gtdbtk.summary.tsv` needs a
  concatenation step that does not exist — GTDB-Tk writes `bac120` and `ar53`
  summaries separately.
- **The AMRFinder stamp does not survive `pixi install`.** Its data lives in the
  conda prefix; rebuild the environment and the database is fetched again.
- **`--tui` has not been run against a failing workflow interactively.** The
  "Nothing ran / no report" path is covered by unit tests and was reached once
  by accident, but not driven by hand since.
- **`/evo/postdoc/cm2v3`** is the old rsync scratch directory, now redundant,
  holding an 8.5 GB pixi environment that can be deleted.
- **The git remote points at `cmkobel/comparem2`** while GitHub has renamed the
  repository to `cmkobel/CompareM2`. It redirects, so nothing is broken, but
  every push prints a notice.

## Deliberately left alone

- **The 69 MB of test zips in git.** Real input data, in HEAD, referenced by
  `tests/README.md`. If they ever have to go it is a git-lfs or
  external-hosting decision, not a cleanup one.
- **The `v2` branch.** It is what the paper describes; leave it.
- **The 35.4 MB of PDFs in history.** See
  [DECISIONS.md](DECISIONS.md#lives-on-branch-v3-in-cmkobelcomparem2) — purging
  them would mean rewriting `master`.
