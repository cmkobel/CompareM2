# Draft: pplacer's final branch-length pass does not catch what the first pass does

**Unsent.** Destination: an issue on
[matsen/pplacer](https://github.com/matsen/pplacer), addressed to matsen, who
wrote both #391's fix (5f4b5186, alpha21) and #394's (40509ec8, alpha22). Send
it **after** [gtdbtk-pplacer-pin-pr.md](gtdbtk-pplacer-pin-pr.md): that text
says "I'll take it up with pplacer", and this is that follow-up. It can go the
same day or later; nothing in the pin PR waits on it.

An issue, not a PR: the change below has not been compiled (pplacer needs OCaml
4.14 and its opam stack, which this laptop does not have), and the trigger is
read from the source, never observed. An issue with the sketch asks matsen to
judge it; a PR would ask him to trust it.

Drafted 2026-10-24 against `master` = `v1.1.alpha23` (f60eef0f), from the
local clone at `~/postdoc/cm2-macos/gtdbtk/alpha-crash/pplacer/`.

## The issue as drafted

**Title:** `Final branch-length pass can still raise FindStartFailure / InvalidStartValues as Child_error`

**Body:**

````markdown
Thanks for the quick fixes in #391 and #394. I've been checking whether GTDB-Tk could move from alpha19 to alpha23 (bioconda-recipes PR NNN), and while reading `core.ml` I noticed one path that #391's fix doesn't reach. I haven't seen it happen, so this is from the source only.

alpha21 catches `Minimization.FindStartFailure` and `InvalidStartValues` in the first optimisation pass and skips the location ([core.ml#L285-L310](https://github.com/matsen/pplacer/blob/f60eef0fe2e536dba0264060c056a41993ab6ae2/pplacer_src/core.ml#L285-L310)). The final branch-length refinement calls the same `safe_ml_optimize_location`, but catches only `Gsl.Error.Gsl_exn` ([core.ml#L369-L382](https://github.com/matsen/pplacer/blob/f60eef0fe2e536dba0264060c056a41993ab6ae2/pplacer_src/core.ml#L369-L382)). Either exception raised there would leave the worker and come back as `Child_error`, ending the run.

One way it could happen: the final pass starts Brent from the first pass's branch lengths. Since #394, `set_dist_bl` clamps the distal length to `[0, cut_bl]` ([three_tax.ml#L41-L48](https://github.com/matsen/pplacer/blob/f60eef0fe2e536dba0264060c056a41993ab6ae2/pplacer_src/three_tax.ml#L41-L48)), so it can sit exactly on a bound, and `start_finder_prep` raises `InvalidStartValues` when the start is not strictly inside the interval ([minimization.ml#L11-L13](https://github.com/matsen/pplacer/blob/f60eef0fe2e536dba0264060c056a41993ab6ae2/pplacer_src/minimization.ml#L11-L13)). I haven't checked whether Brent ever leaves the state on the bound at the end of the first pass.

Something like this would treat them as the `Gsl_exn` case already is, keeping the first-pass result (not compiled):

```diff
             with
               | Gsl.Error.Gsl_exn(_,warn_str) ->
                 dprintf
                   "Warning: GSL problem with final branch length optimization for location %d. %s\n"
                   loc
                   warn_str;
-                initial)
+                initial
+              | Minimization.FindStartFailure
+              | Minimization.InvalidStartValues _ ->
+                dprintf
+                  "Warning: couldn't refine branch lengths for location %d; kept the first-pass values.\n"
+                  loc;
+                initial)
```

A related question, rather than a bug report: if every location for a query is skipped, `failwith "empty results for %s!"` ([core.ml#L333-L335](https://github.com/matsen/pplacer/blob/f60eef0fe2e536dba0264060c056a41993ab6ae2/pplacer_src/core.ml#L333-L335)) ends the whole run, not just that query. In a GTDB-Tk run that's every genome in the batch. Is that intended?

For context, across 234 GTDB genomes, 472 queries built to stress the optimiser on GTDB's short-edged trees, and the 7,825 reads from #391, alpha23 never raised either exception and agreed with alpha19 on every top edge apart from two exact ties. So this is a narrow path, not something I've hit.

Happy to open a PR if you'd prefer.

Best, Carl
````

## Why the claims hold

All at `v1.1.alpha23` (f60eef0f, which is also `master` on 2026-10-24).

- **First pass catches both**: `core.ml:285–310`, the `try … with` around
  `safe_ml_optimize_location (initial_tolerance prefs) loc` (line 288), arms
  `Gsl_exn`, `FindStartFailure` (297), `InvalidStartValues` (303), each
  returning `None` so the location is skipped.
- **Final pass catches one**: `core.ml:373–381`, `try (loc,
  safe_ml_optimize_location final_tolerance loc) with Gsl_exn -> … initial`.
  `safe_ml_optimize_location` (272) itself retries only on `Gsl_exn`, and
  `ml_optimize_location` (255) catches only `ExceededMaxIter`.
- **Both exceptions are reachable from it**: `Three_tax.optimize` →
  `optimize_pend_bl` / `optimize_dist_bl` → `optimize_something` →
  `Minimization.brent` with the default `bisection_start_finder`, which calls
  `start_finder_prep` (raises `InvalidStartValues` when `left >= raw_start ||
  raw_start >= right`, `minimization.ml:12–13`) and raises `FindStartFailure`
  after `max_iters` (line 36). `brent` catches only `FoundMin`.
- **The clamp can produce a start on the bound**: `three_tax.ml:41–48` sets
  `dist` to `max 0. (min dist_bl cut_bl)`; `optimize_dist_bl` (68–71) then
  starts from `get_dist_bl tt` with bounds `0.` and `cut_bl`. Whether the
  stored value is ever exactly on a bound when the first pass ends is **not
  established** — the state is whatever `set_fun` was last called with, and
  `start_finder_prep` itself evaluates `f left` and `f right`.
- **An exception in a worker becomes `Child_error`**: `multiprocessing.ml:244`,
  `305`, `397` (`Fatal_exception exn -> raise (Child_error exn)`).
- **Empty results end the run**: `core.ml:333–335`. This corrects
  [gtdbtk-pplacer-evidence.md](gtdbtk-pplacer-evidence.md) as written on
  2026-10-12, which named `assert(sorted_ml_results <> [])` (line 353) as what
  such a query would reach; the `failwith` comes first, so the assert is
  unreachable on that path. The consequence — the run ends — is the same.
- **"Never raised either exception"**: the 234-genome, 472-query and #391 runs
  in the evidence note all exited 0 with alpha23, and pplacer prints a warning
  when the first pass skips a location. *Not checked:* whether any of those
  logs contain the first-pass "Skipped" warnings. Re-check before sending, or
  soften the sentence to "never failed".

## What is left

- **Before sending:** the `NNN` is the bioconda pin PR's number. Re-check that
  `master` is still f60eef0f (otherwise re-read the three sites and update the
  permalinks), and grep the alpha23 run logs on GenomeDK and the laptop for
  `couldn't find valid optimization start` / `invalid optimization bounds`
  (see the last bullet above).
- **What it commits Carl to:** "Happy to open a PR" — a compiled, tested
  version of the sketch if matsen asks. That needs an OCaml 4.14 build
  environment, which bioconda's recipe (`recipes/pplacer/build.sh`) is the
  quickest route to.
- **After sending:** the issue number goes in [README.md](README.md) and
  [../STATUS.md](../STATUS.md).
