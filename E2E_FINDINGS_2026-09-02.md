# End-to-end test findings — 2026-09-02

Scratch notes from one run, kept in git but **not maintained** and not one of
the three canonical docs. Every count and line reference below is as of
2026-09-02 and several have moved since — the pipeline is 14 tools rather than
13, the unit suite 212 tests rather than 171, and deployment is two shared
environments rather than one per tool. Read it as a dated record, and check
[STATUS.md](STATUS.md) before acting on anything in it.

None of findings 2–5 has been folded into `STATUS.md` or `DECISIONS.md`, and
finding 3 was still true of `catalogue.py` on 2026-09-05: snp-dists and fasttree
are handed `core_gene_alignment.aln`, not the filtered file.

The run: a fresh `git clone` of `496c8b3` on thylakoid, `pixi install`, seven
*Staphylococcus aureus* genomes downloaded from NCBI RefSeq that morning, all
13 tools, `--cores 24`.

Ground truth (PubMLST typing, published AMR genotypes, the Holden et al. 2010
ST239/TW20 recombination) and the report were checked by independent agents
that re-derived numbers from the raw output files rather than trusting each
other's write-ups. Two workflow runs, both resumable if useful later:
- ground truth: `wf_3d31841d-06c`
- verification + adversarial refute panel: `wf_7cdeacb7-5d0`

Scripts are at
`/Users/physarum/.claude/projects/-Users-physarum-comparem2-new/a75a872a-f608-47ac-98d0-a11d6b3e59eb/workflows/scripts/`.
Raw output lives at `thylakoid:/evo/postdoc/cm2-e2e/results_full/`, mirrored
locally (partial, key files only) at `/tmp/cm2e2e/results_full/`.

## What worked, unconditionally

- Fresh clone → `pixi install` → `pixi run pytest` (171 passed, 2.54 s) →
  `test-fast` (8/8) → full 7-genome, 13-tool run: 39/39 Snakemake steps,
  15 min 25 s, 156,517-byte self-contained report (zero `<script>` blocks).
- `amrfinder` under `--use-conda` — the one path STATUS.md listed as
  untested — works: exit 0, identical results to the pixi run, 2 shared
  content-addressed environments for 5 tool rules.
- All 7 MLST sequence types match PubMLST exactly (8/250/36/5/254/239/707),
  independently re-typed from the raw PubMLST allele tables.
- The published ST239 (TW20) chimeric-genome signature reproduces
  quantitatively from the pipeline's own SNP matrix: two independent
  calculations give a 598 kb and 647 kb CC30-derived block against Holden et
  al. 2010's published 635 kb (within 2%).
- One checker recomputed ~300 printed report values against their source
  TSVs/GFF3s/Rtab/newicks/SBMLs: zero arithmetic errors, zero fallback
  renderers, zero missing samples.

## Blocker

**1. CarveMe: 3 of 7 metabolic models are solver-timeout artifacts, and
nothing in the report or logs says so.**
MRSA252, N315, and the draft genome (GCF_026920295.1) hit CarveMe's
hardcoded 600 s SCIP limit (`carving.py:181`, `limits/time=600`) — wall time
613–619 s vs 22–30 s for the other four. Reaction counts: 1,236–1,263 vs
1,609–1,636 (a 22–24% gap) despite comparable or *higher* DIAMOND evidence
going in (MRSA252 has more annotated CDS and more BiGG-gene hits than
NCTC8325, yet 346 fewer reactions). Raising the limit to 3,600 s does not
recover it — re-solved N315 still returned `status=timelimit`, 7 reactions
*fewer* than the 600 s run. `carve_scip.py`'s presolver-off fix (validated
in STATUS.md on one *E. faecium* genome) is necessary but not sufficient;
this run shows it fails to converge on 43% of a fresh 7-genome batch. Every
`carveme.log` (converged and truncated alike) contains only one identical
line (`carve_scip: presolving/milp/maxrounds=0`) — no SCIP status, gap, or
wall-clock reaches any artifact a user could inspect.
*Adversarial panel: 1/2 lenses upheld (survives — 2-lens threshold is
majority-of-2). Dissent argued a careful human reader scanning the 7-row
table could notice the gap by eye; conceded a script or downstream FBA
pipeline reading the raw SBML has no signal at all.*
*One sub-claim did not reproduce under adversarial check: "54 of 114
recoverable reactions match the USA300 BiGG model specifically" — a reviewer
recomputed 25, not 54, under every definition tried. Does not touch the core
finding (the 600 s ceiling, the ~23% reaction-count gap, the failure of a
6× time extension to help).*

**Possible next step:** try the PyPI SCIP wheel (no PaPILO) that
`carve_scip.py`'s own docstring already points at as the alternative fix.

## Important — proven, adversarial panel upheld (2/2 unless noted)

**2. AMRFinder's heatmap bins 121 virulence-gene hits into a fake resistance
class literally labelled "Na".**
`report.py`'s `_section_amrfinder` (~line 1465–1531) admits any non-empty
`Class` string with no `Type` filter. AMRFinderPlus writes `Class=NA` for
every `Type=VIRULENCE` row (121 of 220 rows total, 55%, zero exceptions
either direction). `NA` sorts alphabetically between "Mercury" and
"Quaternary Ammonium" and renders as a column titled "Na", identical styling
to real classes. The shading peak (22, N315's virulence count) is more than
3× the largest genuine resistance/stress class value anywhere (6), so real
classes are compressed into the bottom third of the color ramp. One lens
found the consequence worse than claimed: TW20 (29 real resistance/stress
hits) prints a lower row-end "total" (48) than N315 (49, only 27 real hits)
purely from virulence-gene noise — the totals column reverses a genuine
biological comparison.
*Where:* `src/comparem2/report.py` (`_section_amrfinder`).

**3. `catalogue.py` feeds fasttree and snp-dists the unfiltered core
alignment**, though Panaroo already writes (as a side effect of the `-a
core` flag already used) and documents `core_gene_alignment_filtered.aln`
as "recommended for building core genome phylogenies". Re-running the
pinned FastTree 2.2.0 and snp-dists 0.8.2/1.2.0 on the filtered file:
branch lengths down 26–71%, pairwise SNP counts down 20–60% across the full
7×7 matrix. Topology unaffected. `catalogue.py:454-455,467,478-481`. None of
DESIGN.md/DECISIONS.md/STATUS.md mention this choice.
*One lens downgraded to "less-severe": guidance.py already carries explicit
caveats against reading snp-dists/fasttree numbers as calibrated distances,
and no pair in this dataset sits near a realistic same-strain/outbreak
threshold, so no reader decision in this specific run is shown to flip.*

**4. `mlst.tsv` and skani's `ani.tsv`/`.af` key rows by full absolute file
path, not sample name** — the only 2 of 13 tools that break this
convention (contrast: checkm2, gtdbtk, mashtree, treecluster, snp-dists,
panaroo all key by bare sample name). `join -1 1 -2 1 <(sort mlst.tsv)
<(sort checkm2/quality_report.tsv)` returns 0 rows, exit 0, no warning.
Masked in the HTML report only by `report.py`'s `_sample_of()`
(`Path(value).stem`) — the raw declared/resumable artifact itself is not
portable if the output directory is copied or renamed.

**5. `--until <rule>` silently truncates `report.html` to just the narrowed
tool set, even against an already-complete run with nothing stale.**
`args.until` feeds both the job DAG (`cli.py:286`) and the report renderer
(`cli.py:413`, `render_report(CATALOGUE, args.until, ...)`). Reproduced
directly: running `--until seqkit` against the complete `results_full`
(156,517 bytes, 13 sections) rewrote it to 57,337 bytes, 1 section — while
Snakemake reported "Nothing to be done" for all 39 jobs (no re-execution, no
staleness). No warning, exit 0. A plain full re-run restores it byte-for-byte.
The narrowed report's own summary line reads "1 of 1 tools produced
output" — asserting completeness rather than disclosing omission.
*One lens downgraded to "less-severe": this matches the project's own
long-standing, deliberate design (STATUS.md repeatedly reports narrowed
`--until` runs producing narrowed reports as the expected outcome); `cli.py`
does print "7 assemblies, 1 tools" before running, which an attentive user
would see.*

**6. STATUS.md's "`pixi install` — both environments" is false as
written.** A bare `pixi install` on a fresh clone builds only `default`
(pixi's own output: "✔ The default environment has been installed.",
singular). Directory ctimes on the actual fresh-clone checkout: `.pixi/envs/
default` at 19:36:44, `.pixi/envs/checkm2` at 19:38:58 — 2m14s later, built
only once something ran `pixi run -e checkm2` separately. Not a functional
break (the on-demand build is what makes the two-DIAMOND-version isolation
work at all, and it does — 42 s including the build), but the sentence
overclaims.
*1/2 lenses upheld (survives). Dissent: STATUS.md is an execution log, not
install instructions; every actual how-to document (README, docs/10, pixi.toml's
own inline comment) already gets this right, and pixi's own terminal output
self-corrects the misunderstanding immediately if anyone acts on it — argues
this is a documentation nit, not an "important" defect.*

## Minor / documentation — checker-reported, not sent through the adversarial panel

7. STATUS.md's CarveMe fix claim ("Fixed... 50 s and 1,743 reactions") is
   scoped to one *E. faecium* genome and doesn't hold 43% of the time in
   this run (STATUS.md lines 29, 119–183). Same underlying issue as
   Blocker #1 — file this as the documentation half of that fix.
8. Draft genome (GCF_026920295.1) is missing `fosB`, present in all six
   other genomes — plausibly a fragmentation artifact of the 266-contig
   assembly rather than a real lineage loss; not distinguishable from the
   report.
9. Draft genome shows two `tet(38)` rows — possibly a real duplication or an
   assembly-graph redundancy artifact; not investigated further (would need
   contig-overlap/read-mapping analysis).
10. mashtree's guidance (`guidance.py`) lists assembly-quality causes for an
    "odd placement" but omits recombination as a cause, even though the
    identical caveat already exists in fasttree's guidance ~600 lines later.
11. Draft's 266-contig fragmentation drops it out of otherwise-universal
    gene clusters at 3–15× the rate of the other six genomes, deflating the
    reported core-gene count. Report's own guidance states the general
    caveat but doesn't connect it to this specific genome.
12. Several outputs a real user would want are undeclared, so
    Snakemake's resumability contract doesn't cover them and the report
    can't render them: Bakta's genome diagram (PNG/SVG), GenBank flat file,
    richer per-genome JSON; Panaroo's pangenome graph (`final_graph.gml`).
13. Three database-download rules lack Snakemake provenance metadata even
    on a clean dry-run (not fully characterised — flagged for follow-up).
14. Documented test count (167, in CLAUDE.md and STATUS.md) is stale;
    `pixi run pytest` actually collects 171 — and commit 496c8b3's own
    commit message already says 171. The prose was never updated to match.
15. No documentation discloses that the "3.8 s `pixi install`" and "clean
    clone" timings were measured against a warm pixi/rattler package cache
    (~22 GB pre-existing on thylakoid from other checkouts). A separate,
    informal cold-cache probe against the same commit gives ~26 s and a
    genuine ~6.9 GB network download for `default` alone — the cold path
    reads as untested in STATUS.md when it has, informally, been checked
    and is ~7× slower.

## Retired suspicions (checked, not defects)

- The draft genome's NCBI strain label "RN6432" looked inconsistent with a
  laboratory-derivative-of-NCTC-8325 assumption in the task brief — resolved:
  the assumption was wrong, not the pipeline. RN6432 is Novick's designation
  for the unrelated "Smith diffuse" strain (NARSA NRS148, ST707); the
  pipeline's MLST/ANI/SNP placement is correct.
- TW20's basal placement outside the CC8 clade in both mashtree and fasttree
  is the documented, expected consequence of an unmasked mosaic genome, not
  a computational error in either tool.
- Panaroo's "Soft core (95–99%) = 0" is arithmetically forced at N=7 (no
  integer k satisfies 6.65 ≤ k < 6.93) — and `report.py` already special-
  cases N<20 to show exact counts instead, so this never reaches the reader
  as a false finding; it only appears in Panaroo's own raw
  `summary_statistics.txt`.
- An apparent 3.7% gap between one checker's own supplementary SNP count for
  TW20–MRSA252 and the pipeline's `snp-dists.tsv` figure was fully explained
  by 27 core genes the checker had excluded as paralogous; re-running the
  same scan on the exact file snp-dists consumed matched exactly (33,263).

## Not verified this pass

- Whether a commercial MILP solver (CPLEX/Gurobi) solves the 3 slow CarveMe
  instances quickly — no license available on thylakoid.
- Whether this specific SCIP timeout pathology (as opposed to the
  PaPILO-presolve pathology `carve_scip.py` already documents) is a known
  upstream CarveMe/SCIP issue — worth a literature/issue-tracker check.
- Full verification of every guidance.py prose caveat against its cited
  paper — only caveats intersecting the numeric sections actually audited
  were checked this pass.
