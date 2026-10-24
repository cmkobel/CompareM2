# Evidence: pplacer alpha20's crash, GTDB-Tk's alpha19 pin, and GTDB-Tk on Apple silicon

Behind [gtdbtk-pplacer-pin-pr.md](gtdbtk-pplacer-pin-pr.md) and
[gtdbtk-parallel-spawn-pr.md](gtdbtk-parallel-spawn-pr.md). Measured
2026-10-05 and 2026-10-12, on GenomeDK (linux-64) and the laptop (macOS,
Apple M5, 32 GiB); every claim re-checked 2026-10-12 before either draft was
considered ready (last section).

The question: why GTDB-Tk pins `pplacer =1.1.alpha19` (bioconda-recipes#61248,
GTDB-Tk 2.6.0), and whether the reason still holds for alpha22 and alpha23,
which — unlike alpha19 — have `osx-arm64` builds.

**Where the raw material is.** Neither location is in git.

- Laptop, `~/postdoc/cm2-macos/gtdbtk/`: `alpha-crash/pplacer/` (a clone of
  matsen/pplacer), `alpha-crash/yiran/` (pplacer#391's data and the scripts run
  on it), `alpha-crash/gtdb/` (the stress test), `alpha-crash/solve/` (the
  osx-arm64 environment, the r232 mock database, the Python 3.14 test),
  `alpha-crash/trees/` (published GTDB trees), `pplacer-compare/` (genome
  selection, `compare.py`, comparison outputs), `time.log` / `run.log` (the
  2026-10-05 ANI-screen run), `pr/` (working copies of both diffs).
  `alpha-crash/FINDINGS.md` was the working version of this note.
- GenomeDK, `/faststorage/project/PM_group/carl/comparem2/pplacer-compare/`:
  `env_alpha{19,20,22,23,23s}/` (pixi environments), `run_alpha{19,19b,20,22,23s}/`
  (the 234-genome runs), `shared/` (identify + align, done once), `yiran/`,
  `gtdb/`, `mock/py314/`. Job logs in `../slurm_logs/<jobid>.out`.

## Why the pin exists

GTDB-Tk 2.4.1 and 2.5.2 allowed `pplacer >=1.1.alpha17` and so got alpha20,
and users reported failures: Ecogenomics/GTDBTk#668, #670, #674. GTDB-Tk
2.6.0's release notes: "(#670, #674, #668) Fixed an issue where GTDB-Tk would
crash when using pplacer v1.1.alpha20 … by fixing pplacer to v1.1.alpha19",
and "a fixed version for skani (v0.3.1) and pplacer (v1.1.alpha19) to i) ensure
reproducibility of results". pchaumeil on #668, 2025-11-27: "I will enforce the
`pplacer` version to `v1.1.alpha19` for the time being while investigating what
is happening with `v1.1.alpha20`."

**The three are not the same failure.**

- **#668**, where the pin was decided: one metabat bin, `--skip_ani_screen`,
  pplacer (`dev`, i.e. alpha20) ran 3 h on class-level tree 5 and wrote a JSON
  with no placements; `guppy tog` then failed (`… has no placements!`).
- **#670**: the reporter shows only GTDB-Tk's `PplacerException`. The
  `Child_error` there is a comment by jianshu93 (2025-10-19): died at genome
  149 of 371 in class-level tree 4.
- **#674**: `Dune__exe__Multiprocessing.Child_error(_)`, in the reporter's
  comment of 2025-11-17 (518 genomes; always at genome ~53–56 of 70 in tree 3).

The users were on r220 / r226, not r232. #670's `conda list` shows
`pplacer 1.1.alpha20 hd563303_0` with `gsl 2.8`, the only linux-64 alpha20 build.

## What the builds are

| build | from | OCaml | GSL | linking |
| --- | --- | --- | --- | --- |
| alpha19 `h9ee0642_2` | the 2016 upstream release binary | — | its own | static |
| alpha20 `hd563303_0` | bioconda, from source | 4.14.2 | 2.8 | dynamic |
| alpha22 `hd563303_0` | bioconda, from source | 4.14.2 | 2.8 | dynamic |
| alpha23 `h01ede22_0` | bioconda, from source, 2026-10-10 | 4.14.2 | 2.7 (`libgsl.so.25`) | dynamic |

alpha19 is linux-64 only on bioconda; alpha20, 22 and 23 also have
linux-aarch64 and osx-arm64. The Linux builds print `--version` as `dev`; the
osx-arm64 alpha23 is upstream's release binary and prints `v1.1.alpha23`.

**Placement logic did not change between alpha19 and alpha20.** The diff over
the placement-relevant files (382 changed lines) is a port to current OCaml,
gsl-ocaml and Batteries: `Pervasives` → `Stdlib`, `string` → `bytes`,
`Gsl_min` → `Gsl.Min`, bigarray macro renames in `linear_c.c`, unused functions
deleted. The assertion that fails, `assert(prox_bl >= 0.)` in `three_tax.ml`
`set_dist_bl`, is in alpha19's source too (807f6f3, the alpha19 tag). Both
versions install GSL's error handler (`ppatteries.ml:784`).

The OCaml 5 multiprocessing bug (pplacer#392 / #393, Marshal header 20 → 16
bytes) is not GTDB-Tk's crash on Linux — bioconda's Linux builds are OCaml
4.14.2 — but would hang the osx-arm64 alpha20 binary, which is upstream's OCaml
5 build. Fixed in alpha21.

## Result 1: the `Child_error` is an assertion, fixed in alpha22

pplacer#391's own data (Dropbox link in the issue, live 2026-10-12; files
renamed by the reporter, md5s match once `CONTENTS.json` points at them) and
command, `alpha-crash/yiran/run_yiran.sh`, GenomeDK, 8 cores:

| build | `-j 8` | `-j 0` (single process) |
| --- | --- | --- |
| alpha19 | exit 0, 103 s | exit 0, 625 s |
| alpha20 | **exit 2** after 3,717 of 4,886: `Child_error(_)` | **exit 2** after 3,711: `File "pplacer_src/three_tax.ml", line 43, characters 4-10: Assertion failed` |
| alpha22 | exit 0, 80 s | — |
| alpha23 | exit 0, 81 s | exit 0, 463 s |

`Child_error` is raised in the parent when a worker sends back an `Exception`
(`multiprocessing.ml`, `obj_received`), so it hides what the worker raised.

**One read is enough** (`one_read.sh`): the 932 reference rows plus
`SRR7629282.41436329`, `-j 0`. alpha20 exits 2 on the same assertion; alpha19,
22 and 23 exit 0 and all place it on edge 245, logL −51014.536820, distal
5.029e-3, pendant 1.999993.

Over all 7,825 reads (4,886 unique), alpha22 and alpha23 against alpha19: the
same reads placed, **0 top-edge differences**, max |ΔlogL| on the same edge
5e-6, 2,744 rows bit-identical. alpha22 and alpha23 are identical to each other.

The fix is pplacer#394 (40509ec8, first in alpha22, 2026-01-28): clamp
`dist_bl` to `[0, cut_bl]`, because GSL's Brent minimiser can return a value
just outside its bounds. alpha21 (5f4b5186, #391's first fix) also catches
`Minimization.FindStartFailure` and `InvalidStartValues` — in the first
optimisation pass only (`core.ml:287–303` at v1.1.alpha23).

**Why alpha19 survives the same read is not established** — its 2016 build
might have assertions compiled out, or its older static GSL might never step
outside the bounds here. The #391 reference tree has 46 zero-length edges,
where that interval is `[0, 0]`.

## Result 2: #668's "no placements" is a dead worker, with alpha19 as well

pplacer's event loop (`multiprocessing.ml`, `event_loop`, the same in alpha19,
alpha20 and alpha23) treats a closed pipe as a finished worker
(`End_of_file -> FDM.remove`). A worker that dies by signal — the OOM killer,
a segfault — is dropped and its query lost, and pplacer still exits 0.

`yiran/kill_test.sh`: the #391 data at `-j 2`, SIGKILL to one worker 40 s in.

| build | result |
| --- | --- |
| alpha19 | exit 0, 4,885 of 4,886 placement records (7,824 of 7,825 names) |
| alpha23 | exit 0, 4,885 of 4,886 (7,824 of 7,825) |
| alpha20 | confounded: dies on the line-43 assertion anyway |

In a one-genome run that is a JSON with no placements — #668's symptom. The pin
does not guard against it. **Why #668's worker died is unknown.**

## Result 3: GTDB-Tk's own data

**r232's trees.** None of GTDB-Tk's r232 reference trees has a zero-length
edge. Most class-level trees bottom out at FastTree's 5e-4; the exceptions are
class-level tree 1 (29 edges below 1e-5, minimum 1e-18), tree 4 (5) and ar53
(22, at 1e-6). Backbone minimum 2.16e-3.

**234 genomes through `gtdbtk classify --place_species`** (2026-10-05 and
-12): 232 r232 non-representatives — `gtdb_representative = f` for all 232 in
the r232 metadata; 192 bacterial, 40 archaeal; 116 isolates 95.89–100 %
complete, 116 MAGs 52.94–80 %; contamination ≤ 4.99 %; drawn per class-level
tree with a fixed seed by `pplacer-compare/select_genomes.py` — plus the
*E. faecium* pair `116_2` / `116_2_duplicate`. identify and align were run once
(`prep.sh`); only `classify` differs between builds.

| against alpha19 | classifications differing | top edges differing (of 428) | max \|ΔlogL\| | classify wall |
| --- | --- | --- | --- | --- |
| alpha19, second run | 0 | 0 | — | 6,047 s (first run 6,080 s) |
| alpha20 | 0 | 0 | 0.0013 | 3,708 s |
| alpha22 | 0 | 0 | 0.0013 | 4,307 s |
| alpha23 (`env_alpha23s`) | 0 | 0 | 0.0013 | 4,377 s |

alpha23 against alpha22: no non-zero difference in any placement field. Two
genomes differ in `other_related_references` in every comparison, including
alpha19 against its own second run, so that column is run-to-run variation and
not pplacer. alpha20 did not crash. Outputs:
`pplacer-compare/compare-2026-10-12.txt`,
`pplacer-compare/compare-2026-10-12-with-alpha23.txt`.

*Environments:* env_alpha20 and 22 differ from env_alpha19 in pplacer, GSL
(2.7 → 2.8) and fastani / hmmer build strings; hmmer does not enter `classify`.
env_alpha23 solved to skani 0.4.0 against 0.3.2, which can change
classifications, so the 234-genome alpha23 run used `env_alpha23s` (skani 0.3.2;
otherwise Python 3.13.15 → 3.13.16 and numpy 2.5.3 → 2.5.4). A first alpha23
run with skani 0.4.0 was cancelled 15 min in.

**472 queries built to push the optimiser to its bounds**, `gtdb/make_queries.py`
(fixed seed; byte-identical query files across builds) through GTDB-Tk's own
command (`pplacer -m wag -j N -c <refpkg> -o <json> <msa>`), `gtdb/run_gtdb.sh`:
copies of the leaves on edges < 1e-5 as they are and masked to 50 / 80 / 90 %,
20 random leaves as they are and 80 % masked, 20 with 20 % of residues
substituted, and 20 half-and-half chimeras.

| refpkg | queries (unique) | alpha20 | same top edge as alpha19 | max \|ΔlogL\| |
| --- | --- | --- | --- | --- |
| ar53 | 176 (176) | exit 0, 159 s | 176 of 176 | 0.0114 |
| class-level 4 | 100 (97) | exit 0, 485 s | 100 of 100 | 0.0002 |
| class-level 1 | 196 (179) | exit 0, 474 s | 194 of 196 | 0.0005 |

alpha19, 22 and 23 also exit 0 on all three. The two class-level-1 queries are
ties, not disagreements: equal logL and LWR 1/7 on each of seven kept edges; for
one the seven are the same set in another order, for the other the tie spans at
least nine edges and a different seven were kept. alpha20, 22 and 23 are
identical to each other on all three trees. (The SLURM jobs read FAILED, exit 1:
the script's last `grep` matched nothing under `pipefail`. Every pplacer call
exited 0.)

**Inconclusive — the users' release.** In the *published* trees
(`alpha-crash/trees/`): bac120 r220 has 57 edges < 1e-5, r226 none (minimum
5e-4), r232 40 zero-length; ar53 r220 6 zero-length, r226 4, r232 none. GTDB-Tk's
class-level trees are re-estimated (r232's package trees differ from r232's
published tree), so this says nothing certain about the r226 package trees the
users ran.

## Result 4: GTDB-Tk on Apple silicon, and the same bug on Python 3.14

gtdbtk 2.7.2's bioconda run list with `pplacer =1.1.alpha23` solves with
`pixi lock` on linux-64 (52 packages), linux-aarch64 (52) and osx-arm64 (49)
(`alpha-crash/solve/`).

GTDB's r232 mock database (25.1 MB; `classify_wf -x gz --place_species -f` per
its HOWTO; `solve/mock/`):

| run | result |
| --- | --- |
| osx-arm64, Python 3.13, default `spawn` | hangs at "Running Prodigal", 0 % CPU; faulthandler shows `biolib_lite/parallel.py:186` waiting, and `AttributeError: 'Parallel' object has no attribute '__process_manager'` in a child |
| same, forced `fork` (`gtdbtk_fork.py`) | exit 0, 6.0 s, peak RSS 233 MB, pplacer step included |
| same, `spawn`, the two methods renamed | exit 0, pplacer step included; summary byte-identical to the `fork` run |
| linux-64 (GenomeDK), Python 3.14.8, default `forkserver`, as released | `AttributeError`, then hangs; killed at 600 s |
| same, renamed | exit 0, 24 s, same classification |

The cause, reproduced in eight lines outside GTDB-Tk: a bound method pickles as
`(getattr, (obj, '__name'))` with the *unmangled* name, so a double-underscore
method passed as `mp.Process(target=…)` cannot be found in a `spawn` or
`forkserver` child — exit 1 under both, 0 under `fork`; a single-underscore
method exits 0 under all three. It is the error in Ecogenomics/GTDBTk#619
(macOS Ventura) and #683 (Python 3.14), which bioconda-recipes#61837
(2026-01-13) worked around by capping Python below 3.14. The three
`target=self.__…` lines in `parallel.py` are the only private methods passed to
`mp.Process` anywhere in GTDB-Tk; `master`'s `parallel.py` equals 2.7.2's.

The ANI-screen-only path never reaches this: on 2026-10-05, two *E. faecium*
genomes on the laptop went through `classify_wf` in 20.27 s ("Identify and
Align steps will be skipped"), maximum RSS 9,155,149,824 B — macOS `time -l`
reports the largest single process.

The mock database's classification is meaningless by design (its HOWTO says
so). This shows the pipeline completes on arm. **No real genome has been placed
by pplacer on arm.**

## Re-check on GTDB-Tk 2.8.0, 2026-10-24

GTDB-Tk 2.8.0 was released 2026-10-13 and replaced 2.7.2 in the bioconda
recipe the same day (bioconda-recipes#70223, build 0). Same r232 reference
package; the changes are to the ANI screen and species assignment in
`classify` / `classify_wf`. It touches neither `pplacer.py` nor anything under
`biolib_lite/`, and its `parallel.py` is byte-identical to 2.7.2's. The recipe
still pins `pplacer =1.1.alpha19` and `python <3.14`, and its run requirements
are 2.7.2's line for line.

- **Solve.** The 2.8.0 run list with `pplacer =1.1.alpha23`
  (`alpha-crash/solve-2.8.0/`): linux-64 52 packages, linux-aarch64 52,
  osx-arm64 49 — the same counts as 2.7.2's. Python 3.13.17, skani 0.4.0,
  GSL 2.7 on all three.
- **Mock database, osx-arm64** (`solve-2.8.0/mock/run_mock_mac.sh`; gtdbtk
  2.8.0 pip-installed `--no-deps`, Python 3.13.17, pplacer `v1.1.alpha23`;
  each run killed at 180 s):

  | run | result |
  | --- | --- |
  | as released, `spawn` | hangs after "Running Prodigal"; `AttributeError: 'Parallel' object has no attribute '__process_manager'`; killed at 180 s |
  | as released, forced `fork` | exit 0, 5.9 s, max RSS 228 MB |
  | `spawn`, the two methods renamed | exit 0, 5.2 s, "Running pplacer v1.1.alpha23"; summary byte-identical to the `fork` run |

- **Python 3.14 on Linux: not re-run on 2.8.0.** Prepared in
  `solve-2.8.0/mock/py314-280/` (Python 3.14.9 in the lock); needs a GenomeDK
  login. The 2.7.2 result above stands.
- **Not re-run, deliberately:** the 234-genome and 472-query comparisons. They
  compare pplacer builds on the same alignments under one GTDB-Tk version, and
  2.8.0 does not change how pplacer is called (`pplacer.py` untouched). It does
  rewrite top-hit selection in `files/marker/tophit.py` and `copy_number.py`
  (`sorted(…, reverse=True)[0]` → `max(…)`), which feeds the alignment; that
  should pick the same hit, and is not checked.

## Not shown

- That #670 and #674 crashed on the line-43 assertion rather than on one of the
  two exceptions alpha21 catches. No GTDB-Tk input tried (234 genomes, 472
  queries) crashed alpha20; it would take one of their inputs.
- Why #668's worker died.
- Why alpha19 survives the read that kills alpha20.
- From the source, never observed: alpha23's final branch-length pass
  (`core.ml:373–381`) still catches only `Gsl.Error.Gsl_exn`, so
  `FindStartFailure` / `InvalidStartValues` raised there would still surface as
  `Child_error`; and a query whose every location is skipped ends the run at
  `failwith "empty results for %s!"` (`core.ml:333–335`). A plausible trigger
  for the first: #394's clamp (`three_tax.ml:41–48`) can leave the distal
  length exactly on a bound, and `brent` raises `InvalidStartValues` for a
  start not strictly inside its interval (`minimization.ml:12–13`); whether the
  state is ever on a bound when the first pass ends is not established. Both
  are in [pplacer-final-pass-issue.md](pplacer-final-pass-issue.md).
  *Corrected 2026-10-24:* this bullet said such a query would reach
  `assert(sorted_ml_results <> [])` (`core.ml:353`). The `failwith` at 333
  comes first, so that assert is unreachable on this path; the consequence is
  the same.

## Verification pass, 2026-10-12

Every sentence of both drafts re-checked against the runs, the source at the
tags, and the threads, before either was considered ready. What it changed:

- #668 had been described as a `Child_error` crash; it is "no placements"
  (Result 2 came out of this). #670's `Child_error` had been attributed to
  #674's reporter; it is jianshu93's.
- The 234 genomes had been called 194 + 40 GTDB non-representatives; they are
  232 plus the *E. faecium* pair.
- 455 stress queries had mixed in a deduplicated count; there were 472.
- alpha21's catch had been described as covering every pass; it covers the
  first.
- The alpha23 environment had skani 0.4.0; the 234-genome run was redone with
  0.3.2 before "alpha22 and alpha23 identical throughout" went into the draft.
- The 9.2 GB figure is the largest process, not the total.
- "The placement code is the same" had been checked over four files; the full
  diff was then read.

Confirmed as written: the #61248 diff (`>=1.1.alpha17` → `=1.1.alpha19`); the
2.6.0 release-note wording; line 43 at v1.1.alpha20 and columns 4–10; which tag
first carries each fix; the final pass calling `Three_tax.optimize` →
`Minimization.brent`, which raises both exceptions; #61837 as the source of
`python <3.14`; no other private `mp.Process` target; `master` = 2.7.2 for
`parallel.py`; byte-identical mock summaries.
