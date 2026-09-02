# Repo cleanup plan

Written 2026-09-02. Delete this file when the last phase is done and record the
outcome in DESIGN.md.

Goal: this branch tracks v3 and nothing else.

Two phases: land the content, then verify it on linux. No history rewrite —
see below.

---

## Decided: history is not rewritten

A full `git filter-repo` repack was chosen first and **reversed the same day**
once the cost was measured, which is the right call. What the measurement showed:

| Target | Size | Refs affected |
| ------ | ---: | ------------- |
| `papers/` — accidental, only ever on `v3` | 37.4 MB, 17 blobs | 1 commit on 1 branch |
| raw genome FASTA — predates the branch split | 157.2 MB, 79 blobs | `master`, `v3`, `new`, `ai-1` |

The FASTAs sit in commits *below* the branch point, so purging them would
rewrite `master` — every published SHA of CompareM2 v2, across all 10 remote
branches, breaking every fork and clone. And scoping to `--refs v3` recovers
almost nothing, because those blobs stay reachable from `master`. It was
genuinely all-or-nothing, and the "all" side spends the published history of a
paper to save disk.

**So: no rewrite, no force-push.** `papers/` is untracked going forward by an
ordinary commit and the 35.4 MB stays in history. If it ever has to go, the
moment to do it is alongside some other deliberate rewrite, not on its own.

Everything below is Phase 1 and 2 only — ordinary commits, all reversible.

---

## Phase 1 — content — DONE 2026-09-02

Landed as three commits on `v3`: `1d84072` (untrack papers), `cd0c941`
(guidance + defect fixes + drop sylph), `c7e8d91` (remove v2, rewrite entry
points and docs). Not pushed — held until Phase 2 passes on linux.

Two things came out differently from the plan below, both recorded in the
commit messages:

- **`docs/30 what analyses does it do.md` is generated too**, not just the
  citation page. It renders each tool's real command line from the spec, so it
  cannot describe something that does not run. `docs/generate.py --check` fails
  in CI if either page is stale.
- **`--set` had a defect the skani fix exposed.** `parse_overrides` split on
  the first `--` and prepended `--` to the rest, so a single-dash flag was
  inexpressible: `skani-c=125` was rejected and `skani--c=125` produced
  `-c 70 --c 125`. Fixed, with tests.

<details><summary>The plan as written</summary>

### 1.1 Land what is already staged
Currently in the index: `src/comparem2/guidance.py` and `tests/unit/test_v3.py`
added, the sylph removal and three defect fixes, the narrowed `.gitignore`, the
DESIGN.md entries, and 15 staged deletions untracking `papers/`.

Split into two commits so the defect fixes are not buried in the guidance work:

```bash
git commit -m "papers: untrack the accidentally committed PDFs"   # the 15 deletions + .gitignore
git commit -m "report: explain every tool from its own paper"      # guidance.py, report.py, tests, DESIGN.md
```

This does not reclaim the 35.4 MB — the blobs stay reachable from `952a7d0`.
It stops the directory growing and makes the branch correct, which is all that
was ever on offer without a rewrite.

### 1.2 Delete v2
~104 files, 5.5 MB. Nothing in v3 reads any of it: v3 generates its own
`envs/*.yaml` into a build directory, so `workflow/envs/` is not a dependency.
Verified by grep before writing this.

```bash
git rm -r workflow dynamic_report profile config resources
git rm comparem2 Dockerfile environment.yaml changelog.txt \
       comparem2.Rproj readme-development.md
```

| Path | Files | Why it goes |
| ---- | ----: | ----------- |
| `workflow/` | 40 | Snakefile + 11 `.smk` + 25 env yamls, all generated from specs in v3 |
| `dynamic_report/` | 22 | the R/RMarkdown report; v3's report is Python |
| `profile/` | 15 | conda/apptainer × local/HPC profiles; v3 uses Snakemake's SLURM executor |
| `resources/` | 8 | includes 5.2 MB `ko00001.json` for clusterProfiler, which is dropped with R |
| `config/` | 1 | `config.yaml`; v3 takes `--set tool--flag=value` |
| `comparem2` | 1 | the 30 KB v2 launcher, replaced by `src/comparem2/cli.py` |
| `Dockerfile` | 1 | built from the 25-env conda tree |
| `environment.yaml` | 1 | superseded by `pixi.toml` |
| `changelog.txt` | 1 | v2 releases; start fresh or move to GitHub releases |
| `comparem2.Rproj` | 1 | RStudio project, no R left |
| `readme-development.md` | 1 | v2 dev instructions |

**Keep:** `src/`, `tests/` (8 input zips, README, unit tests), `DESIGN.md`,
`LICENSE`, `citation.cff`, `CONTRIBUTORS.md`, `pixi.*`, the CompareM2 paper PDF.

### 1.3 Rewrite `pixi.toml`
Its only task is `test-fast`, which shells out to `./comparem2` — broken the
moment 1.2 lands. DESIGN.md already claims `pixi run cm2` is the entry point and
no such task exists.

- add `cm2 = "python -m comparem2.cli"` (or a `[project]` script entry)
- add `pytest = "pytest tests/unit -q"`
- drop what looks v2-only: `mamba`, `conda`, `pandas`
- add `textual` (imported by `tui.py`) and `pytest`, `pytest-asyncio`
- keep `snakemake = "<10"` and the cluster-generic executor plugin

Then confirm on linux: `rm -rf .pixi && pixi install && pixi run pytest`.

### 1.4 Replace CI
Delete the 6 v2 workflows; add one job running the unit suite on
`ubuntu-latest`. The 75 tests are pure Python and finish in 0.7 s — no
bioinformatics dependencies, no databases, no test genomes needed. Without this,
tracking the tests in Phase 1.1 buys nothing.

Keep a manually-triggered (`workflow_dispatch`) end-to-end job for later, once
there is a linux runner with the databases.

### 1.5 Rewrite `CLAUDE.md` and `README.md`
`CLAUDE.md` is the urgent one: it currently documents v2's three-layer design,
`workflow/rules/*.smk`, 25 conda environments and "bump the version in three
places". Every one of those is false for v3, and it is loaded into context each
session, so it will actively mislead the next one.

### 1.6 Rewrite `docs/`
Keep `docs/`, `mkdocs.yml` and `.readthedocs.yaml` on this branch and rewrite
the content for v3, so docs and code stay in step where v3 lives.

| File | Change |
| ---- | ------ |
| `30 what analyses does it do.md` | 30+ tools → 13; drop the N-dependent table (v3 has a dependency closure); new pseudo-targets |
| `20 usage.md` | `--config set_x=y` → `--set tool--flag=value`; new CLI and TUI |
| `10 installation.md` | pixi; one environment plus isolated CheckM2, not 25; new database totals |
| `99 citation.md` | regenerate from `guidance.py` — it is already the single source of truth for citations |
| `05 quick start.md`, `index.md` | rewrite around the 13-tool set |

`99 citation.md` should be **generated**, not hand-maintained: `GUIDANCE`
already holds every citation with its DOI and the `citations()` helper
deduplicates them. Hand-editing it is how v2's list drifted.

Until this lands, readthedocs keeps building from `master`, so users see v2 docs
for v2 — which is correct, not a regression.

### 1.7 Push
`git push origin v3` — an ordinary fast-forward, no force. Held until Phase 2
passes on linux, since the branch is published.

</details>

---

## Phase 2 — verify on linux — DONE 2026-09-02

Run on thylakoid (`/evo/postdoc/cm2v3`, 24 cores) on the four *E. faecium*
genomes. `cm2v3` is gitignored scratch inside the `postdoc` repo, so the source
was rsynced rather than pushed.

- [x] `pixi install` — both environments, from a manifest never installed before
- [x] all four `pixi run` tasks resolve (`cm2`, `pytest`, `unpack`, `test-fast`)
- [x] 81 unit tests, 0.89 s
- [x] `docs/generate.py --check` — pages current
- [x] `pixi run cm2 --help`
- [x] `pixi run test-fast` — 8/8 steps, report rendered, 4 guidance blocks
- [x] **`skani -c 70` runs** — the one flag taken from a paper and never
      executed. Accepted, logged, duplicates at 100.00% ANI
- [x] `fasttree` at `threads=1` — duplicates at 0.0 branch length
- [x] panaroo exact-count table on real data — 2,091 core of 3,780
- [x] `snp-dists` — 0 between the duplicate pair, over a 1,934,948 bp alignment
- [x] `bakta` — 4/4 against db v6.0 light
- [x] `mlst` — ST32/ST32/ST117/ST78
- [x] **`--isolated-launcher`** — CheckM2 from its own environment, never
      exercised end to end before
- [x] `carveme` — 4 SBML models, ~4 MB each; ~9 min per genome
- [x] **`amrfinder`** — after the download mechanism landed: database
      `2026-08-07.1` fetched in 26 s as a workflow step, then 7 and 11 genes
- [x] `mkdocs build --strict`
- [x] CI green — 20 s and 22 s on the two pushes
- [x] `git push origin v3` — `b6a5f3b`
- [x] **clean clone works** — cloned to `/evo/postdoc/comparem2`, both
      environments installed, 92 tests, `test-fast` 8/8, and the runbook
      command with `-d` and `--isolated-launcher` ran 8/8 including CheckM2
- [ ] `gtdbtk` — never run, needs the 141.4 GB download

Three defects the run found, all fixed and recorded in DESIGN.md:

1. **No database downloads existed.** `Database.url` was read by nothing.
   Downloads are now Snakemake rules.
2. **GTDB-Tk's database was unreachable** — it reads `GTDBTK_DATA_PATH` and
   nothing set it, so `--databases` was ignored for 91% of the install weight.
3. **Snakemake locked the working directory, not `--output`.** Fixed with
   `--directory`.

## Where things live on thylakoid

- checkout: `/evo/postdoc/comparem2` (branch `v3`, a real clone)
- databases: `/evo/postdoc/cm2-databases` — 6.9 GB, deliberately outside any
  checkout so deleting a checkout does not cost a re-download. Pass it with
  `-d`.
- `/evo/postdoc/cm2v3` is the old rsync scratch directory and is now redundant;
  it holds an 8.5 GB pixi environment that can go.

Two flags are not optional there:

```bash
pixi run cm2 <assemblies>... \
  -d /evo/postdoc/cm2-databases \
  --isolated-launcher "$HOME/.pixi/bin/pixi run -e {tool}"
```

Without `-d`, the 6.9 GB is downloaded again. Without an **absolute** launcher
path, CheckM2 fails with `command not found`, because Snakemake's shell does
not inherit an interactive PATH.

## Not in scope

- **The 69 MB of test zips.** Real input data, in HEAD, referenced by the test
  README. If they ever need to go it is a `git-lfs` or external-hosting
  decision, not a cleanup one.
- **Deleting the `v2` branch.** It is the tag on what the paper describes;
  leave it.
