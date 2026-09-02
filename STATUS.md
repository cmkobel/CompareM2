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
| **gtdbtk** | **in progress** | GTDB r232 downloading since 16:21 (60.8 GB, ~5.6 h at the measured 3.0 MB/s), then classify runs unattended. **Three defects fixed first** — see below |

### GTDB-Tk: three defects, none of which a solve would show
Found 2026-09-02 while making the tool runnable. Each was invisible because the
rule had never been executed, and the first two were described in comments as
though they already existed:

1. **The batchfile was never written.** `--batchfile` named a path no rule
   created, so the command would have died on its first line. Now `Tool.files`.
2. **The summaries were never merged.** GTDB-Tk writes `bac120` and `ar53`
   separately, so the declared `gtdbtk.summary.tsv` could not appear and the
   job would have failed on a missing output even with a working command. Now
   `Tool.post`.
3. **`--skip_ani_screen` does not exist in 2.7.2.** It belonged to the versions
   that screened with Mash; 2.7 screens with skani and dropped the flag.
   `unrecognized arguments` was the only way to find this, short of reading the
   installed tool's `--help` — which is how it was found.

And a fourth, upstream of all of them: **the database was the wrong release.**
The catalogue pointed at r226; `gtdbtk/config/common.py` in 2.7.2 reads
`COMPATIBLE_REF_DATA_VERSIONS = ['r232']`, so the tool would have rejected the
data *after* a 141.4 GB download. r232 is 60.8 GB, less than half, because it
replaced FastANI's reference genomes with skani sketches.

### The gtdbtk report section is written, against documented columns only
`_section_gtdbtk` exists and renders: the lineage every genome shares stated
once, then a column per rank where they diverge, ANI coloured against *that
reference's* species radius rather than a global 95%, and a note when a genome
carries no species name. Twelve tests, including a real render.

**Every input to those tests is synthetic.** The columns come from GTDB-Tk's
own documentation, not from a file produced here, so the renderer checks its
headers and falls back to the generic table if they are not the ones it
expects. Four things to check against the real `gtdbtk.summary.tsv`, in order
of how quietly they fail:

1. **Does r232 emit `closest_genome_reference_radius` under that exact name?**
   If not, the ANI colouring degrades to plain text — the correct failure, and
   invisible unless looked for.
2. **How long are the real `classification_method` strings?** That column drives
   the table width; the synthetic worst case was 51 characters.
3. **Does `user_genome` hold the batchfile's genome id** (the sample name) or a
   path? `_sample_of` handles both, but only the first joins cleanly to every
   other section.
4. **Is the merged file the two domains concatenated once**, with one header —
   `merge-tsv`'s own contract, but it has never run on real GTDB-Tk output.

### The report has never covered more than four tools
Checked 2026-09-02, because "12 of 13 verified" invites the wrong inference.
Those twelve were verified across *separate* runs, and the largest report on
thylakoid has four sections — seqkit, checkm2, mlst, skani (31 KB, 12:38);
another has two, bakta and carveme. No run has produced a report covering the
catalogue.

Nothing suggests it would fail. But the report is the product, and a
genomes-to-report pipeline published without one whole report having been
produced is the gap least worth defending.

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
- **The package builds and installs.** `pip install --no-deps .` into a bare
  venv (Python 3.14, macOS): both entry points present, `comparem2 --version`
  and `cm2 --version` print `CompareM2 3.0.0.dev0`, and the metadata carries
  `License-Expression: GPL-3.0-or-later`.
- **The PATH preflight fires.** The same installed package, no tools on PATH:
  `not on PATH: seqkit (seqkit), mashtree (mashtree)` and exit 1, before any
  Snakemake call.
- **`--use-conda` reaches Snakemake with the right flags** — Snakemake 9.26.1
  echoed `--software-deployment-method conda --conda-prefix …` and began
  building the DAG. See below for where that run stops.
- **`--tui`, end to end** — `--until mashtree treecluster --tui` under a real
  terminal: two tools selected from the command line, both `done`, report
  written, clean exit. And headlessly through `run_test()` with `seqkit` and
  `checkm2`, so the isolated launcher goes through the TUI path too. Five
  defects had to be fixed first, see [DECISIONS.md](DECISIONS.md).
- 147 unit tests, ~1.6 s — 13 added for the conda-deployment path, 8 for the
  steps around GTDB-Tk's command, 16 for the report rewrite
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
| GTDB r232 | 60.8 GB | — | measured; 97% of the measured total. **r226 was wrong**: gtdbtk 2.7.2 accepts only r232 |
| CheckM2 | 1.7 GB | 2.9 GB | measured |
| Bakta light | 1.3 GB | 4.0 GB | download figure is Bakta's documented one; on-disk measured |
| AMRFinder | unmeasured | — | version `2026-08-07.1`, 26 s to fetch. Lands in `$CONDA_PREFIX`, **not** under `--databases` |

Software is 1.58 GB against 62.5 GB of measured data, and GTDB-Tk is 60.8 GB
of it. It was 143 GB until the release changed.

## Where things live on thylakoid

24 cores, 125 GB RAM, 914 GB free on `/evo` (measured 2026-09-02, before the
60.8 GB download).

| | |
| --- | --- |
| checkout | `/evo/postdoc/CompareM2` — a real clone, on `master`. **The path followed the repository rename**; the old lowercase path does not exist |
| databases | `/evo/postdoc/cm2-databases` — 6.9 GB, deliberately **outside** any checkout so deleting a checkout does not cost a re-download |
| pixi | `/home/thylakoid/.pixi/bin/pixi` |

```bash
cd /evo/postdoc/CompareM2
export PATH=$HOME/.pixi/bin:$PATH
export COMPAREM2_DATABASES=/evo/postdoc/cm2-databases

pixi run pytest          # 147 tests, no tools or databases needed
pixi run test-fast       # 4 genomes, no databases needed

pixi run cm2 my/*.fna \
  --output results_myrun \
  --cores 24 \
  --isolated-launcher "$HOME/.pixi/bin/pixi run -e {tool}"
```

**The export belongs in `.bashrc`, and the launcher path must be absolute.**
The default database directory is `~/.comparem2/databases`, which on thylakoid
is under `/home` and not the `/evo` volume the existing copy sits on, so
without `$COMPAREM2_DATABASES` (or `-d`, which has to be retyped every run) the
6.9 GB is fetched a second time, into home. Without an **absolute**
launcher path, CheckM2 fails with `command not found`, because Snakemake's
shell does not inherit an interactive PATH.

To skip the 60.8 GB:

```bash
--until seqkit checkm2 bakta amrfinder mlst mashtree treecluster skani \
        panaroo snp-dists fasttree carveme
```

## Packaging

The code side of the bioconda package is done; the release is not. What exists:
`pyproject.toml`, the `comparem2`/`cm2` entry points, `--use-conda` and
`--conda-prefix` through both execution paths, and a draft recipe in
`recipe/`. The model is *pipeline only* — see
[DESIGN.md](DESIGN.md#two-deployment-models).

| | |
| --- | --- |
| environments a full conda run builds | **14** (17 rules needing one, deduplicated by content) |
| published recipe today | `comparem2` **2.16.2**, `noarch: generic`, maintainer `cmkobel` |
| version here | `3.0.0.dev0` — the release needs a real one, in three files |

**A conda deployment has never been executed.** The dry run above stops at
`conda: command not found`, because this laptop has no conda and the tools are
`linux-64` regardless. That run is thylakoid work, and it is the one thing
between here and a recipe that can be trusted:

```bash
pixi run pip install -e .          # or a plain venv
comparem2 tests/E._faecium/*.fna --until seqkit mashtree \
  --use-conda --conda-prefix /evo/postdoc/cm2-envs
```

Two things to watch when it runs: whether AMRFinder's database survives into
its analysis rules (the shared-environment reasoning in DESIGN.md, untested),
and how long fourteen solves actually take.

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
- **No bioconda package**, but the blocker has moved: what remains is a tag, a
  sha256, the PR, and one verified `--use-conda` run. See *Packaging* above and
  `recipe/README.md`. **A hand-built container image is no longer planned** —
  decided 2026-09-02, see [DECISIONS.md](DECISIONS.md).
- **One unit test is failing in the working tree**, and it is not the packaging
  work: `test_pangenome_matrix_compresses_identical_patterns` counts `<rect>`
  elements, while the in-progress `report.py` figure rewrite draws one `<path>`
  per genome. `HEAD` itself is green at 110 tests.
- **The thylakoid checkout carries uncommitted `src/` changes**, rsynced from
  the laptop for the GTDB-Tk run rather than pulled. Once the work is committed,
  `git checkout . && git pull` there to get back to a clean tracked state.

## Deliberately left alone

- **The 69 MB of test zips in git.** Real input data, in HEAD, referenced by
  `tests/README.md`. If they ever have to go it is a git-lfs or
  external-hosting decision, not a cleanup one.
- **v2 itself.** Not preserved on a branch, deliberately — a decision taken on
  2026-09-02 when v3 was merged to `master`. `origin/v2` is a 2019-era branch,
  2,008 commits behind the pre-merge `master`, and is *not* what the paper
  describes; the paper-era state is the pre-merge `master` tip, in history and
  on the `ai-1` branch. It is also **tagged**: `origin` carries 39 tags up to
  `v2.16.2` (`ad352f2`), which is what bioconda's published recipe builds from.
  An earlier version of this line said the `v2.*` tags stopped at `v2.9.1` —
  that was this checkout's stale tag list, not the remote's; `git fetch --tags`
  brings the rest. Read the docs at
  [readthedocs `/en/stable/`](https://comparem2.readthedocs.io/en/stable/)
  rather than trying to reconstruct a v2 checkout.
- **The 35.4 MB of PDFs in history.** See
  [DECISIONS.md](DECISIONS.md#lives-on-branch-v3-in-cmkobelcomparem2) — purging
  them would mean rewriting `master`.
