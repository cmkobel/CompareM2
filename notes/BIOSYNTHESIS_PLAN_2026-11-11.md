# The biosynthesis review's nine findings — plan, 2026-11-11

## Outcome, 2026-11-11

Executed on branch `biosynthesis-review`, seven commits, `8583617`..`f5b8906`.
Every commit was checked on all 55 review models on thylakoid before it went
in. Across those models, the branch's tables differ from master's in two
places: `iSB619`'s `carbon` cell empties, and `iAF692` completes instead of
crashing.

**Where it departed from this plan:**
- **Commit 3 was rejected on its measurement.** Closing every non-exchange
  boundary moved no verdict in 55 models, and it stopped `iJN1463` and
  `iYS854` growing on any medium: their sinks stand in for a peptide and a
  protein the network cannot make. What shipped instead is the negative
  controls as `upstream/biosynthesis_controls.py`, plus a scope sentence in
  `medium_constraints`.
- **`iYS854`'s leak is not its sinks**, as the plan says below. With every
  sink shut it still makes 22 compounds without nitrogen, through `SRK`, a
  reaction that does not balance.
- **The solver rebuild is a precaution, not a repair.** On `iAF692` the
  outputs are identical with and without it.
- **The starved note names only the infeasibility** for a model with no
  solution on M9. The source check fails there too, so it is not a second
  reason.
- **The calibration prose moved into commit 6**, with the calibration
  itself. Commit 7 kept the rest of the prose.
- **The `iML1515` citation note** still said "the 29-of-30 figure". The
  example report's re-render is what caught it.

The plan as written follows, unedited.

Written against `d0dee44` (master), after a review that ran the module on 55
models on thylakoid: 13 curated BiGG models, 17 pipeline drafts, the 6 CarveMe
1.5.0 benchmark drafts, and 19 drafts carved that day through `carve_scip.py`.
Data and the throwaway instruments are in `/evo/postdoc/cm2-biosynth-review-1111`
(`bin/controls.py`, `bin/cobra_check.py`, `bin/summarise.py`; results in `out/`).

**What held, and must keep holding through every step below:**
- cobrapy, written from the docstring's definitions and solved with HiGHS (53
  models) and GLPK (34), agreed on all 1,590 verdicts. M9 and LB growth agreed
  to 4 decimals.
- None of the 42 CarveMe drafts made a panel compound on an empty medium, or on
  M9 without C, N, S, P, Fe or Co.
- No draft has a supply-capable boundary reaction other than `R_EX_*_e`, a
  forced flux, a mis-oriented exchange, or a non-optimal solve.
- Master reproduces the stored rc-3.5.1 *E. faecium* tables byte for byte.

So nothing here changes a verdict on a CarveMe draft. That is the regression
gate for the refactor: **all 42 draft panel TSVs byte-identical before and
after.**

## Shape of the change

Most findings are symptoms of three structural habits, and fixing those makes
the individual fixes small.

1. **Three drain mechanisms where one would do.** `add_demands`,
   `add_source_demands` and `add_precursor_demands` are the same function three
   times, with two prefixes and a three-way `_drain_ub` assembly. Amino acids are
   both panel compounds and biomass precursors, so they get two drains each.
2. **The probe's question is spread over two calls.**
   `probe.maximum(drain, probe.medium(compounds, drain)) > min_flux` appears five
   times, and `min_flux` is threaded through five signatures. Every caller wants
   one boolean: can this model make metabolite X on medium Y?
3. **The report reads the media TSV by column index.** `raw(sample, medium, 6)`,
   `len(row) > 5` and `cell(s, m, 3)` are scattered through
   `_section_biosynthesis`. "Did this model grow" is computed once for the
   summary sentence and not reused by `_closest_note`. Finding 1 is exactly that
   gap.

## Commits, in order

### 1. Refactor `biosynthesis.py`'s probe — no behaviour change

- **One drain per metabolite:** `_Probe.drains: dict[met_id, rxn_id]`, built by a
  single `_drain(model, met)` that keeps the collision check. Its ids are
  `R_CM2_DM_<met minus M_>`, e.g. `R_CM2_DM_trp__L_c`. Panel compounds, source
  forms and biomass precursors all register through it, which deduplicates the
  amino acids. Only one drain is ever open, so answers cannot move.
- **`probe.reaches(met, compounds) -> bool`** is the only question asked.
  `min_flux` and `max_uptake` live on the probe, from the CLI. `maximum` and
  `medium` become private.
- **`verdicts(probe) -> dict[bigg, verdict]`** in panel order, and
  `rescuers(probe, verdicts)`. `main()` builds the TSV rows from `PANEL` once.
- **The media become data:** `REFERENCE_MEDIA = (("M9", M9, True),
  ("M9[-O2]", …, True), ("LB", LB, False), …)`, where the flag says whether
  `SOURCES` is checked. That replaces `name.startswith("M9")`.
- Update `upstream/panel_calibration.py.report()` and the two test stubs
  (`_StubProbe`, `_RescueProbe`), which get simpler because each implements
  `reaches`.

**Gate:** the 42 draft panel and media TSVs come out byte-identical, and the
unit tests pass.

### 2. A solver that fails is a status, not a crash (finding 6)

- `_solve` catches the exception, tallies `ERROR = "Error"`, returns `(ERROR,
  None)` and **rebuilds the solver**. Rebuilding is not optional: reframed's
  `Solver.solve` has no `try/finally`, so the failed call's temporary bounds
  stay applied, and SCIP is left mid-solve. `flux_from` already maps it to 0,
  and `growth_cell` gives `error`. `main()`'s existing "N of M solves were not
  optimal" line reports it, with the first message.
- In the report, the unsolved-complete-medium note's "usual cause is a
  maintenance reaction" applies to `infeasible` only. Say something true for
  `error`.
- Test: a probe whose `_fba` raises once. It is tallied, the verdict is zero,
  the solver is reset, and the next solve succeeds.
- **Evidence:** `iAF692` dies today on the `atp_c` precursor drain on M9 with
  SCIP "unresolved numerical troubles in LP", and must complete afterwards.

### 3. The medium is total: nothing else supplies (finding 7)

`medium_constraints` promises "exactly `compounds` available, and nothing
else", which is false for a model with sinks. In `_Probe.__init__`, classify
every one-sided reaction that is neither `R_EX_*_e` nor ours. Pass the bounds
that shut its supply direction into `medium_constraints` as a fourth mapping:
`lb = max(lb, 0)` for `A -->`, `ub = min(ub, 0)` for `--> A`. It stays a pure
function, and gets a test.

- **Drafts:** no change. The audit found no such reaction in 42.
- **Why at all:** `panel_calibration.py --model` exists to score curated models,
  and curated models have sinks. `iYS854` makes 22 compounds without nitrogen
  through `SK_apoACP_c` and `SK_2fe1s_c`, and `iIT341` makes cysteine from
  nothing through `SK_ahcys_c`.
- Re-run the negative controls on the curated set. Leaks gone, and any verdict
  that moves gets recorded.
- **Not done:** forced lower bounds (`ATPM`) stay the model's. They are
  physiology, not supply. Commit 4 makes the report say when they break M9.

### 4. Carbon arrives as glucose *or* glucose 6-phosphate (finding 5)

- `SOURCES` becomes `Source(element, compound, forms)`:
  `("carbon", "glc__D", ("glc__D", "g6p"))` and
  `("nitrogen", "nh4", ("nh4",))`. An element is reached when any form's drain
  carries flux.
- g6p is the first cytoplasmic form on the PTS route as well as after
  hexokinase. `iSB619` reaches g6p and pyruvate but never free glucose, so today
  it reads "no carbon" beside 25 *de novo*.
- **Rejected: probing pyruvate instead.** It is downstream, so a model could read
  "carbon unreachable" with compounds still *de novo*. That would falsify the
  note's "every compound falling to upstream or no route at once". The entry
  forms keep that sentence true.
- **Not covered, and the comment says so:** a draft whose glucose enters only by
  periplasmic oxidation to gluconate. No measured model does.
- `_starved_note` widens from "cannot reach an element" to "M9 is unusable for
  this model". The reasons are an unreachable element, or an M9 growth cell
  holding a status: `iNF517` is infeasible on M9 and reads *no route* for
  everything today.
- Reword "a property of the reconstruction and not of the organism" to what is
  known: the column says nothing about that genome. Either the organism does
  not live on glucose and ammonium, or the draft lost the route. Both occur:
  the Methanoflorens drafts make nothing from CO2 + H2, acetate, formate or
  methanol either, and the pneumococci lost `NH4tpp`.
- Tests: replace `test_the_source_elements_are_m9s_own_and_carry_their_element`
  with the new shape. Add a stub case where only g6p is reached and carbon is
  not reported.

### 5. The report reads its tables by name (finding 1)

- One tolerant reader per table, keyed by `PANEL_HEADER` and `MEDIA_HEADER`
  (already importable), returning a small `NamedTuple` per medium: growth
  float-or-None, raw cell, missing, unreachable, reached, blocked. A short or
  half-written row simply has empty fields. That is the defensive reading
  today's code does by index, in one place.
- `grew` is one predicate, used by the summary sentence and by
  `_closest_note`, which now covers only models that grew on no defined medium.
  The two sentences can no longer contradict each other.
- Test: a model growing on M9 that fails M9[-O2] gets no closest-medium line.
  The existing all-fail and all-grow tests stay.
- **Evidence:** rendered today, the note called the 1.5.0 *E. coli* draft (LB
  4.83) a model that "did not grow". Five of six benchmark drafts trigger it,
  and so would any strict aerobe.

### 6. Calibration is a table of expected answers, not a count (finding 2)

- `upstream/panel_calibration.py` scores **all six** bundled proteomes by
  default, against a table of what each should *not* make, with a reason per
  entry:
  - `adocbl`, ceiling, for every genome;
  - `mqn8` for *P. aeruginosa* and *R. solanacearum*, which respire with
    ubiquinone only. **Check a source before committing either.**
  - `asn__L` for *R. solanacearum*, which has no asparagine synthetase and does
    have GatCAB;
  - *M. genitalium* as the negative control.

  It prints deviations both ways: false dependencies and false prototrophy.
  Today it would print one, *S. oneidensis* `pro__L`, untraced. Its draft has
  `GLU5K`, `G5SD` and `P5CR`, the same as the *E. coli* draft that makes
  proline.
- `tests/environments/check.py` compares *E. coli*'s non-de-novo **set** to
  `{"adocbl"}` instead of counting to 29. A draft that loses one and gains
  another passes a count. Update CLAUDE.md's "must come out 29 of 30" line in
  the same commit.
- **Commit the two review instruments**, so 0-of-1,590 and 0-of-42 can be
  re-measured:
  - `upstream/biosynthesis_controls.py`: negative-control media and the boundary
    audit, on the new `reaches` API;
  - `upstream/biosynthesis_crosscheck.py`: the cobrapy re-implementation. It
    needs a venv with `cobra highspy osqp`, and the docstring says so.

  Not the RefSeq proteomes: they need NCBI, so their numbers go in STATUS.

### 7. The prose says what was measured (findings 2, 3, 4, 8, 9)

- **Module docstring, rewritten rather than amended.** It describes the module
  as it now is. Dated narrative already in DECISIONS goes. The calibration
  paragraph and the asparagine paragraph are replaced. The superseded "BiGG
  universe reaches all 32" paragraph becomes the `universe_bacteria` result
  (29 of 30, `adocbl` absent). `none` is defined to include salvage-only
  cofactors.
- **`guidance.py`:**
  - *De novo*: today's numbers, from commit 6's table, and why 29 is not every
    prototroph's target.
  - *No route*: ThDP and NAD in *S. aureus* are reachable from thiamine and
    nicotinate, which the panel never offers. The verdict says "must be given
    the vitamin", not "no path at all".
  - *Asparagine*: no annotated synthetase and GatCAB present in both
    *S. aureus*, *A. baumannii* and *R. solanacearum*. The free-asparagine demand
    is the biomass reaction's, and curated `iYS854` says the same. **Before
    writing what the organism needs, check AsnRS and the literature**;
    otherwise say only what was measured.
  - Measurements labelled with run and date, per the sourcing rule. No
    literature comparisons (C. difficile, H. pylori) unless their papers are
    checked first.
- **`report.py`:** the summary's "no route means there is no route at all"
  (`:2297`), and the stale "32 compounds" (`:2222`).
- **Bounds:** pass `use_infinity=True` explicitly in `main()`, with the reason:
  it is CarveMe's own convention, and verdicts were measured invariant. The
  complete-medium figure is therefore higher than COBRA tools give:
  `iML1515` 30.76 against 20.73 at ±1000. **Rejected: `False`.** It would move
  every complete-medium number in the record and change no verdict.
- **Fold in `PAPERS_PLAN_2026-11-11.md` item 1.** It rewrites the same growth
  caveat and needs the same thylakoid session.
- `python docs/generate.py`. **Re-render `docs/assets/example-report.html`** on
  GenomeDK with `--report-only` from `results_showcase`, if a warmed
  ControlMaster makes that possible. It has been patched twice; commits 4 and 5
  change what it says, not just how. If GenomeDK is out of reach, patch it and
  say so in STATUS, as before.

## Records

- **STATUS.md** gets a new section, "biosynthesis review, executed 2026-11-11":
  the cross-check, the negative controls, the fresh 1.6.6 table (bundled six
  plus the RefSeq seven), curated against draft for the same organisms, and the
  `iAF692` / `iSB619` / `iYS854` cases.
- **DECISIONS.md** gets one dated entry, which records the rejections so they
  are not retried:
  - a fifth `salvage` verdict, which would bring back the vitamin forms the
    one-member-per-family rule keeps out;
  - dropping `asn__L` like `btn`/`q8`: it varies (de novo in *E. coli*,
    *K. pneumoniae*, *L. monocytogenes*), so it fails that test;
  - pyruvate as the carbon probe;
  - `use_infinity=False`;
  - an alternative-carbon fallback for non-glucose organisms, measured useless
    on the six archaea.
- **CLAUDE.md:** the test count, after each commit that adds tests, and the
  check.py line.

## Verification

- **Unit tests:** pixi will not install on the laptop, so use a venv with CI's
  own line (`pytest pytest-asyncio textual "snakemake>=9,<10"`, then `-e .`),
  plus `docs/generate.py --check`.
- **After each of commits 1–4,** re-run the controls, the summary and the HiGHS
  cross-check on all 55 models on thylakoid. Expected:
  - draft TSVs identical throughout;
  - `iAF692` completes (2);
  - `iYS854` and `iIT341` pass the controls (3);
  - `iSB619`'s `unreachable` empties (4).
- Render the report from today's tables after commit 5. `ecol`, `paer` and
  `bsub` should drop out of the closest-medium note.
- After pushing, `gh workflow run environments.yaml` for check.py's carve. That
  is Carl's call, since it is a long run.
