# Draft: move gtdbtk's pplacer pin from alpha19 to alpha23

**Unsent.** Destination: a PR against
[bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes),
`recipes/gtdbtk/meta.yaml`, from a fork branch, addressed to the recipe
maintainers (pchaumeil, aaronmussig). Send it **after**
[gtdbtk-parallel-spawn-pr.md](gtdbtk-parallel-spawn-pr.md) is open, and replace
`#NNN` below with that PR's number. The two go out the same day; do not wait
for the first to be merged or released.

Drafted 2026-10-12 and checked claim by claim the same day against the runs,
the pplacer source and the GitHub threads. Every number below is in
[gtdbtk-pplacer-evidence.md](gtdbtk-pplacer-evidence.md).

A PR rather than an email, decided 2026-10-12: the pin is pchaumeil's own line
(#61248) in a recipe he and aaronmussig maintain, so a PR goes straight to the
people who can approve it, and they only have to click. It is worded as a
request, with an offer to close — what to avoid is a PR merged *around* them.

## The PR as drafted

**Title:** `Update gtdbtk: pin pplacer 1.1.alpha23 instead of alpha19`

**Body:**

````markdown
This moves gtdbtk's pplacer pin from `=1.1.alpha19` to `=1.1.alpha23` and bumps the build number.

@pchaumeil @aaronmussig, you pinned alpha19 in #61248 because of the pplacer failures under alpha20 (Ecogenomics/GTDBTk#668, #670, #674), and Pierre wrote on #668 that it was for the time being, while you looked into alpha20. I've looked into it, and I think the crash is found and fixed in pplacer, so I'd like to ask whether you'd be happy to move the pin. If you'd rather keep alpha19, I'll close this.

**The crash.** #674, and a comment on #670, show `Dune__exe__Multiprocessing.Child_error(_)`, which is how pplacer reports any exception in one of its worker processes, so it hides the real error. The data linked from matsen/pplacer#391, which reported the same message, reproduces it with bioconda's alpha20 build. With `-j 0` (no workers) the real error shows:

```
File "pplacer_src/three_tax.ml", line 43, characters 4-10: Assertion failed
```

That's `assert(prox_bl >= 0.)` in `set_dist_bl`. matsen/pplacer#394 (alpha22) clamps the distal length to the edge instead, because GSL's Brent minimiser can return a value just outside its bounds. One read from that dataset is enough: alpha20 exits 2, while alpha19, alpha22 and alpha23 all exit 0 and place it on the same edge with the same likelihood. alpha21 also catches the minimiser's `FindStartFailure` and `InvalidStartValues` in the first optimisation pass, which would otherwise surface as the same `Child_error`.

**#668 is a different failure.** There pplacer finished but wrote no placements. pplacer drops a worker that dies without raising an exception, one killed by the OOM killer for instance, and carries on without its query: with one worker killed partway through the #391 data, alpha19 and alpha23 both exit 0 with 4,885 of 4,886 placements. That is what a dead worker looks like in a one-genome run, with alpha19 as well, but I can't tell why the worker died in #668.

**Compared with alpha19**, on linux-64 with bioconda's builds, in environments that differ in pplacer (and the GSL it links) and otherwise only by patch or build versions:

- 234 genomes through GTDB-Tk 2.7.2's `gtdbtk classify --place_species`: 232 r232 non-representatives (192 bacterial and 40 archaeal; half isolates, half MAGs at 53–80 % completeness) and two *E. faecium* genomes of ours. No classification differs, and all 428 placements land on the same top edge.
- 472 queries built to push the optimiser to its bounds, run with GTDB-Tk's own pplacer command against ar53 and class-level trees 1 and 4, the ones with edges shorter than 1e-5: copies of the leaves on those edges, masked to 50–90 %, plus mutated and chimeric alignment rows. All land on alpha19's top edge, except two whose weight is spread evenly over tied edges, where the builds break the tie differently.
- All 7,825 reads in the #391 data land on alpha19's top edge.

alpha22 and alpha23 gave identical placements throughout, and an exact pin keeps the reproducibility the alpha19 pin was there for.

**What I couldn't show.** None of the GTDB-Tk inputs crashed alpha20, so I can't say whether #670 and #674 hit this assertion or one of the two exceptions alpha21 catches. One path is still open in alpha23: the final branch-length pass (`core.ml:373`) catches only GSL errors, so those two exceptions could still escape from there. I haven't seen that happen, and I'll take it up with pplacer.

**Why I'm asking.** alpha19 only exists for linux-64, while alpha23 has bioconda builds for linux-64, linux-aarch64 and osx-arm64, and this recipe is `noarch: python`, so with this change gtdbtk's run requirements solve on all three. On Apple silicon it also needs one code change, which I've sent as Ecogenomics/GTDBTk#NNN: `biolib_lite/parallel.py` hands the private methods `__producer` and `__process_manager` to `mp.Process`, and under name mangling they can't be pickled. That's the `'Parallel' object has no attribute '__process_manager'` in Ecogenomics/GTDBTk#619 and #683. macOS starts processes with `spawn` and Python 3.14 on Linux with `forkserver`, so the same bug is behind the `python <3.14` cap from #61837. With that change, GTDB-Tk 2.8.0 and alpha23 run your r232 mock database end to end on an Apple silicon Mac, pplacer step included. Full placements still need the memory they need on Linux, so a 32 GB laptop will only get through genomes the ANI screen classifies. That's still very useful to be able to do: on my 32 GB laptop, two *E. faecium* genomes went through the ANI screen in 20 s, with the largest process at 9.2 GB.

I'm happy to share the scripts and outputs.

Best, Carl
````

## The diff

Against `recipes/gtdbtk/meta.yaml` as of 2026-10-24 (gtdbtk 2.8.0, build 0,
from bioconda-recipes#70223 on 2026-10-13; blob `6e27d0d`). The working copy
is `~/postdoc/cm2-macos/gtdbtk/pr/meta-2.8.0.yaml`, diff `recipe-2.8.0.diff`.
The 2.7.2 versions (`meta.yaml`, `recipe.diff`) are kept beside them. 2.8.0's
run requirements are 2.7.2's line for line, so only the version, hash and
build number moved.

```diff
@@ -10,7 +10,7 @@
   sha256: 42e9b04db4490ff12ce5288d07d76b0148649a6b019cd32c7d4728b5d61779e2

 build:
-  number: 0
+  number: 1
   noarch: python
   entry_points:
     - gtdbtk = gtdbtk.__main__:main
@@ -28,7 +28,7 @@
     - dendropy >=4.1.0
     - prodigal >=2.6.2
     - hmmer 3.*
-    - pplacer =1.1.alpha19
+    - pplacer =1.1.alpha23
     - fastani >=1.32
     - fasttree >=2.1.9
     - numpy >=1.9.0
```

## What the text commits Carl to

- **`#NNN` has to exist** when this goes up: the GTDB-Tk PR is sent first.
- **"I'll take it up with pplacer"**: an issue to
  [matsen/pplacer](https://github.com/matsen/pplacer) about the final
  branch-length pass (`core.ml:373` at v1.1.alpha23) catching only `Gsl_exn`,
  and the related `failwith "empty results"` (`core.ml:333–335`) that ends the
  run if every location of a query is skipped. **Drafted 2026-10-24:**
  [pplacer-final-pass-issue.md](pplacer-final-pass-issue.md). (This line named
  `assert(sorted_ml_results <> [])` at `core.ml:353` until 2026-10-24; the
  `failwith` comes first, so the assert is unreachable on that path.)
- **"I'm happy to share the scripts and outputs"**: nothing is public. They
  are in `~/postdoc/cm2-macos/gtdbtk/` and on GenomeDK under
  `pplacer-compare/` (paths in the evidence note).

## What is left

- **Before sending, re-check three things that can move:** the recipe
  (`gh api repos/bioconda/bioconda-recipes/contents/recipes/gtdbtk/meta.yaml`
  — still 2.8.0, build 0, blob `6e27d0d`?), whether pplacer has released past
  alpha23 (then decide whether to pin that, and re-run the comparisons), and
  whether either open path in the "What I couldn't show" paragraph has been
  fixed upstream.
- **Re-checked 2026-10-24.** GTDB-Tk 2.8.0 (2026-10-13) had replaced 2.7.2 in
  the recipe, which still pins `pplacer =1.1.alpha19` and `python <3.14`;
  pplacer is still at alpha23 (`master` = the tag, f60eef0f), and neither open
  path is fixed. The 2.8.0 run list with alpha23 solves on linux-64 (52
  packages), linux-aarch64 (52) and osx-arm64 (49), as 2.7.2's did, and the
  mock-database run on the Mac was repeated on 2.8.0 (evidence note, *Re-check
  on GTDB-Tk 2.8.0*). The 234-genome comparison was not re-run: it is between
  pplacer builds on the same alignments under one GTDB-Tk, which the text now
  names, and 2.8.0 does not touch `pplacer.py`.
- **After sending:** the PR number goes in the table in [README.md](README.md)
  and in [../STATUS.md](../STATUS.md), and the `#NNN` in this file.
- **Unverified, and the text does not claim it:** no real genome has been placed
  by pplacer on arm. The mock database's classification is meaningless by
  design; it shows the pipeline completes.
