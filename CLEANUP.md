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

## Phase 1 — content (ordinary commits, no rewrite, fully reversible)

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

---

## Phase 2 — verify on linux

Not done until these pass on thylakoid
(`/evo/postdoc/cm2v3`, 24 cores):

- [ ] `rm -rf .pixi && pixi install` — the environment still solves without sylph
- [ ] `pixi run pytest` — 75 tests
- [ ] `pixi run cm2 --dry-run` on the four *E. faecium* genomes
- [ ] **`skani -c 70` actually runs** — added on macOS from the paper, never
      executed; the flag spelling is unverified
- [ ] `fasttree` still produces a tree at `threads=1`
- [ ] panaroo section renders the small-N exact-count table on a real run
- [ ] CI green on GitHub

---

## Not in scope

- **The 69 MB of test zips.** Real input data, in HEAD, referenced by the test
  README. If they ever need to go it is a `git-lfs` or external-hosting
  decision, not a cleanup one.
- **Deleting the `v2` branch.** It is the tag on what the paper describes;
  leave it.
