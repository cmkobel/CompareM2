# Draft: move gtdbtk's pplacer pin from alpha19 to alpha23

**Sent 2026-11-09 as [bioconda-recipes#71426](https://github.com/bioconda/bioconda-recipes/pull/71426)**, with the body
below and the title `Update gtdbtk: pin pplacer 1.1.alpha23 instead of alpha19`.
Branch `cmkobel/bioconda-recipes:gtdbtk-pplacer-alpha23`, one commit `4efc2ad`
on `d962af0`, made through the API: that base was already in the fork via
#71400's branch, so the new ref adds no workflow change and needs no
`workflow` scope. The recipe blob was still `6e27d0d` on upstream `master`
(`7da5106`). Destination: a PR against
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

@pchaumeil @aaronmussig, you pinned alpha19 in #61248 after the failures under alpha20 in Ecogenomics/GTDBTk#668, Ecogenomics/GTDBTk#670 and Ecogenomics/GTDBTk#674, and Pierre wrote on [#668](https://github.com/Ecogenomics/GTDBTk/issues/668) that it was for the time being, while you investigated alpha20. I've looked into it, and I'd like to ask whether you'd be happy to move the pin. If you'd rather keep alpha19, I'll close this.

- **The `Child_error` in [#670](https://github.com/Ecogenomics/GTDBTk/issues/670) and [#674](https://github.com/Ecogenomics/GTDBTk/issues/674)** is how pplacer reports any exception in a worker. The one I could reproduce, with the data from matsen/pplacer#391, is a failed assertion in alpha20, fixed in alpha22 (matsen/pplacer#394). No GTDB-Tk input I tried crashed alpha20, so I can't say theirs was the same one.
- **The genome attached to [#668](https://github.com/Ecogenomics/GTDBTk/issues/668)** goes through alpha19, alpha20 and alpha23 alike, with the reporter's pplacer settings (`--pplacer_cpus 5`, every genome placed): the same edge and likelihood in class-level tree 5, and the same classification. That's on r232; I haven't tried r226.
- **alpha23 places like alpha19** on everything I've tried: 234 genomes through `gtdbtk classify`, on 2.7.2 and again on 2.8.0 (no classification differs, all 428 placements on the same edge), 472 queries built to push the optimiser to its bounds, and the 7,825 reads from matsen/pplacer#391. An exact pin keeps pplacer at one version, as before.
- **Why I'm asking:** alpha19 exists only for linux-64, while alpha23 has bioconda builds for linux-64, linux-aarch64 and osx-arm64, so with this change gtdbtk installs on all three. On linux-aarch64 that was enough for unmodified GTDB-Tk 2.8.0 with alpha23 to run your r232 mock database end to end, pplacer included. macOS also needs Ecogenomics/GTDBTk#723, which you merged this week (thank you). Once a release carries it, I can send the PR that lifts the `python <3.14` cap from #61837 as well.

On Linux, bioconda's alpha23 reports `pplacer --version` as `dev`, as alpha20 does; #71400 fixes that (CI green, awaiting review), and this pin picks up the rebuild without a further change.

The scripts and outputs are in [this gist](https://gist.github.com/cmkobel/527f5b02932f5eac7d6604fb52515f8d).

<details>
<summary>Details</summary>

**The crash.** With `-j 0` (no workers), alpha20 on the matsen/pplacer#391 data shows the real error:

```
File "pplacer_src/three_tax.ml", line 43, characters 4-10: Assertion failed
```

That's `assert(prox_bl >= 0.)`. matsen/pplacer#394 replaced it with a clamp, attributing the failure to GSL's Brent minimiser returning values just outside its bounds. One read from that dataset is enough: alpha20 exits 2, while alpha19, alpha22 and alpha23 exit 0 and place it identically. alpha21 also catches the minimiser's `FindStartFailure` and `InvalidStartValues` in the first optimisation pass.

**[#668](https://github.com/Ecogenomics/GTDBTk/issues/668)** is a different failure: pplacer finished but wrote no placements, and `tog` failed on the empty file. pplacer silently drops a worker that dies without raising (one killed by the OOM killer, say), and its query is lost. With one worker killed partway through the [#391](https://github.com/matsen/pplacer/issues/391) data, alpha19 and alpha23 both exit 0 with 4,885 of 4,886 placements, so the pin doesn't guard against this symptom.

**The comparisons**, on linux-64 with bioconda's builds, in environments that differ in pplacer and the GSL it links, and otherwise only in patch or build versions (alpha20 included wherever it didn't crash):

- 234 genomes through `classify --place_species`, on GTDB-Tk 2.7.2 (alpha19, alpha20, alpha22, alpha23) and again on 2.8.0 (alpha19, alpha23): 232 r232 non-representatives (192 bacterial, 40 archaeal; half isolates, half MAGs at 53–80 % completeness) and one *E. faecium* genome of ours, twice.
- 472 queries, run with GTDB-Tk's own pplacer command on ar53 and class-level trees 1 and 4, the trees with edges shorter than 1e-5: copies of the leaves on those edges, masked to 50–90 %, plus random, mutated and chimeric rows. All land on alpha19's top edge except two exact ties.
- alpha22 and alpha23 placed identically throughout.

**What I couldn't show.** alpha23's final branch-length pass ([`core.ml:373`](https://github.com/matsen/pplacer/blob/v1.1.alpha23/pplacer_src/core.ml#L373)) still catches only GSL errors, as both of alpha19's passes do, so `FindStartFailure` and `InvalidStartValues` could escape there. I haven't seen it happen, and I'll raise it with pplacer.

**Apple silicon.** With [#723](https://github.com/Ecogenomics/GTDBTk/pull/723) applied, GTDB-Tk 2.8.0 and alpha23 run your r232 mock database end to end on a Mac, pplacer included. Full placements still need the memory they need on Linux, but genomes the ANI screen classifies need far less: two took 9.2 GB at most on a 32 GB laptop, with GTDB-Tk 2.7.2.

</details>

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

- **Done 2026-10-27:** the GTDB-Tk PR is open as #723, and the body above cites it.
  **Merged 2026-11-09**; the body now says so, and that it is in the next release.
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
- **Re-checked 2026-11-09**, after #723 merged: the recipe is still blob
  `6e27d0d` (2.8.0, build 0, `pplacer =1.1.alpha19`, `python <3.14`) with no
  open gtdbtk PR; pplacer is still alpha23 (`f60eef0f`), with no issue since
  #396; neither open path is fixed. #619 and #683 are closed.
- **Shortened and re-verified 2026-11-09**, 757 → 581 words. Every claim
  checked again against its primary source: #61248's diff and author, the
  2.6.0 release notes, pchaumeil's #668 comment of 2025-11-27, the
  `Child_error` in jianshu93's #670 comment and minhasbushra's #674 comment,
  #668's "has no placements!", #391 and #394, `three_tax.ml:43` at
  v1.1.alpha20, the `FindStartFailure` / `InvalidStartValues` catch absent at
  alpha20 and present at alpha21 (`core.ml:297–303`), `core.ml:373` catching
  only `Gsl_exn` at alpha23, `multiprocessing.ml`'s `Child_error` wrapper and
  `End_of_file -> FDM.remove`, bioconda's pplacer subdirs (alpha19 linux-64
  only; alpha23 linux-64, linux-aarch64, osx-arm64), #61837 as the source of
  `python <3.14`, #723 merged into `staging`, and the error text in #619 and
  #683. The numbers match the evidence note. A fresh `pixi lock` of the
  recipe's run list with `pplacer ==1.1.alpha23` solves on all three
  platforms (51 / 51 / 48 packages without gtdbtk itself). One sentence was
  dropped as imprecise rather than fixed: the mangled methods *are* pickled,
  under their unmangled name, and it is the lookup in the child that fails.
- **Devil's-advocate pass, 2026-11-09, which made it longer again (707
  words) and found two things that block sending.** Fixed in the text: bare
  `#668`, `#670`, `#674` and `#391` in a bioconda-recipes PR autolink to
  bioconda-recipes PRs from 2015–16 (pathtools, cage, pybedtools, seqbuster)
  and would post backlinks on them — every cross-repo reference is now
  qualified or linked; "the crash is found" became "a crash … that gives the
  same error", since #670/#674 are not shown to be it; the Brent overshoot is
  attributed to #394, which we did not establish; "two *E. faecium* genomes"
  is one genome twice; "exact pin keeps the reproducibility" became "keeps
  pplacer at one version", because osx-arm64's alpha23 is upstream's OCaml 5
  binary, not bioconda's source build; the 9.2 GB figure is labelled 2.7.2.
  **Blocking, open:** (1) #668's genome is attached to the issue
  (`bin.1.txt`, 5.0 MB, 6 contigs, a metabat bin, r226, class-level tree 5,
  `--skip_ani_screen`) and has never been run here. The reporter says the
  pip + alpha19 environment (also Python 3.10, hmmer 3.3.2) put it through,
  and pchaumeil could not reproduce. The #668 paragraph now says only that.
  (2) bioconda's Linux alpha20 *and* alpha23 print `pplacer --version` as
  `dev` (run here with `pixi exec`, on linux-aarch64 — this laptop is
  aarch64; linux-64 from the evidence note and the #668 jobs), the string pchaumeil used on #668 to spot
  the bad build; pplacer#397 is *in* alpha23 (the tag is one commit past its
  merge) and reads `PPLACER_VERSION` at build time; bioconda's
  `recipes/pplacer/build.sh` never sets it. (This note first said #397 was
  unreleased, copying STATUS.md's 2026-10-05 line, written two days before
  alpha23.) The fix is one line in that recipe, no pplacer change. Also
  found: alpha19 catches only `Gsl_exn` in both passes (`core.ml:274, 289,
  362` at v1.1.alpha19), so alpha23's final-pass gap is one alpha19 has
  everywhere; the body now says so. Disclosed in a new
  paragraph; fixing it in `recipes/pplacer` first would be better.
- **#668's genome run, 2026-11-09**, GenomeDK jobs 7741880–83, in
  `pplacer-compare/issue668/` (`run.sh`; the issue's `bin.1.txt`, md5
  `47a07fe9…`, as `genomes/bin.1.fa`). Identify + align once in env_alpha19,
  then `gtdbtk classify --place_species --cpus 10 --pplacer_cpus 5` in
  env_alpha19 / 20 / 23s against r232. All exit 0; each writes one placement
  to the backbone and one to **class-level tree 5** (the reporter's r226 tree
  number too); edge 3917 / 28954, logL −19182386.0294 / −16607059.6889, LWR 1
  in all three; backbone distal and pendant differ in the sixth decimal
  (alpha19 0.245169 / 0.0624994, alpha20 and 23 0.245168 / 0.0624987), tree
  5 identical; all three *s__Mycobacterium abscessus*, ANI 97.28 to
  GCF_000069185.1. Classify wall 1,011 s (alpha19), 655 s, 651 s; MaxRSS
  ~130 GB each. linux-64 `pplacer --version`: alpha19
  `v1.1.alpha19-0-g807f6f3`, alpha20 and alpha23 `dev`. **Not tried: r226.**
  The #668 paragraph now says this instead of "I can't rule out".
- **Restructured 2026-11-09** after the #668 run: the ask and four bullets
  up top (the version a maintainer can approve from), everything else in a
  collapsed `<details>`. *One thing this makes worse* is gone on the
  assumption that bioconda-recipes#71400 merges first; if it has not, put it
  back. `LINK` is a placeholder for the published scripts and outputs. The
  linux-aarch64 sentence is evidence *Result 6*.
- **2.8.0 re-run and the gist, 2026-11-09.** The 234-genome comparison again
  on GTDB-Tk 2.8.0 (GenomeDK jobs 7752342–44, `pplacer-compare/v280/`, the
  October locks of env_alpha19 and env_alpha23s with gtdbtk 2.8.0 pip-installed
  `--no-deps`; identify and align redone on 2.8.0): 0 of 234 classifications,
  0 of 428 top edges differ, max |ΔlogL| 0.0013 — the 2.7.2 figures exactly;
  classify 7,108 s (alpha19) and 4,852 s (alpha23). Evidence published as a
  public gist, https://gist.github.com/cmkobel/527f5b02932f5eac7d6604fb52515f8d
  (31 files, README maps each claim to its files; #391's data linked, not
  republished; no SLURM wrappers, which carry an email address). Every
  comparison table in it was regenerated from the run directories on
  GenomeDK and matches the evidence note's numbers.
- **Sent-now variant, 2026-11-09:** one sentence on the `dev` version string
  and #71400 added after the bullets, so the PR can go before #71400 merges
  (the two are independent: `=1.1.alpha23` fixes no build number, and conda
  takes the highest build). Also: "the reporter's settings" narrowed to
  "pplacer settings" (GTDB-Tk 2.7.2 and r232 here, 2.4.1 and r226 there), and a
  repeated sentence in the #668 details paragraph dropped. **Once #71400
  merges:** confirm `pplacer --version` prints `v1.1.alpha23` from bioconda on
  linux-64 and linux-aarch64, and comment on the pin PR if it is still open.
- **After sending:** the PR number goes in the table in [README.md](README.md)
  and in [../STATUS.md](../STATUS.md), and the `#NNN` in this file.
- **Unverified, and the text does not claim it:** no real genome has been placed
  by pplacer on arm. The mock database's classification is meaningless by
  design; it shows the pipeline completes.
