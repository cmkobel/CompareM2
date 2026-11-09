# Draft: make GTDB-Tk's `Parallel` work under `spawn` and `forkserver`

**Sent 2026-10-27 as [Ecogenomics/GTDBTk#723](https://github.com/Ecogenomics/GTDBTk/pull/723)**, to `staging`, with the text below. **Merged 2026-11-09 at 08:08Z** by pchaumeil, as `5f13497` on `staging`, with: "the fix will be available in the next GTDB-Tk release. This should also resolve #619 and #683". Not released as of that day; the latest is 2.8.0. Reframed before sending to lead with Linux on Python 3.14, which breaks on x86 today, so it stands even if macOS is not their concern. Destination: a PR against
[Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk), **`staging`**, from
[cmkobel/GTDBTk:parallel-spawn](https://github.com/cmkobel/GTDBTk/tree/parallel-spawn)
(commit `70497b7`, pushed 2026-10-27, one commit on `staging`'s `e552e6e`).
`staging` and not `master`: their CONTRIBUTING.md asks for branches based on
`staging` and PRs to it, and `parallel.py` is the same blob on both (`cbf0882`).
**Send this one first**:
it is the easy yes (it fixes Python 3.14 for their Linux users and asks nothing
of them), and [gtdbtk-pplacer-pin-pr.md](gtdbtk-pplacer-pin-pr.md) cites its
number.

Found 2026-10-12 while checking the pin PR's macOS claim: with the pin lifted,
gtdbtk installs on Apple silicon, then hangs at "Running Prodigal" on any genome
the ANI screen does not classify. The ANI-screen-only path is not affected
(identify is skipped). Evidence:
[gtdbtk-pplacer-evidence.md](gtdbtk-pplacer-evidence.md), *Result 4*.

## The PR as drafted

**Title:** `Fix Parallel hanging on Python 3.14 and macOS`

**Body:**

````markdown
Hi,

GTDB-Tk hangs on Python 3.14. At "Running Prodigal" a child process dies with

```
AttributeError: 'Parallel' object has no attribute '__process_manager'
```

and the main process waits for it forever. This is the error in #683, and I believe it's why the bioconda recipe caps Python below 3.14 (bioconda/bioconda-recipes#61837).

The cause is in `biolib_lite.Parallel`, which passes its two private methods, `__producer` and `__process_manager`, to `mp.Process`. That works with `fork`, but with `spawn` and `forkserver` the target has to be pickled, and Python pickles a bound method by its unmangled name, so the child can't find it. Python 3.14 changed the default on Linux from `fork` to `forkserver`, which is why this shows up now. macOS has defaulted to `spawn` since Python 3.8, which is #619 and how I ran into it.

The fix is just renaming the two methods to `_producer` and `_process_manager`. Nothing else in the package refers to them. I've also added a small test that runs `Parallel` under each start method. It times out under `spawn` and `forkserver` without the fix.

I've tested it with 2.8.0 on the mock database (`classify_wf -x gz --place_species -f`):

- Linux (x86-64), Python 3.14: hangs as released (I killed it after 10 minutes), finishes in 30 s with the fix.
- macOS (Apple silicon), Python 3.13, with pplacer 1.1.alpha23 since alpha19 has no arm build: hangs as released, finishes with the fix, and gives the same summary as forcing `fork`.

Let me know if you'd like anything done differently.

Kind regards, Carl
````

## The diff

`master`'s `parallel.py` is byte-identical to 2.7.2's (checked 2026-10-12) and
to 2.8.0's (checked 2026-10-24: 2.8.0 changes nothing under `biolib_lite/`),
and the diff applies to it unchanged. Working copies: `~/postdoc/cm2-macos/gtdbtk/pr/`
(`parallel.py.orig`, `parallel.py`, `gtdbtk-parallel.diff`).

```diff
--- a/gtdbtk/biolib_lite/parallel.py
+++ b/gtdbtk/biolib_lite/parallel.py
@@ -87,7 +87,7 @@
 
         self.cpus = cpus
 
-    def __producer(self, producer_callback, producer_queue, consumer_queue):
+    def _producer(self, producer_callback, producer_queue, consumer_queue):
         """Process data items with producer processes.
 
         Parameters
@@ -108,7 +108,7 @@
 
             consumer_queue.put(rtn)
 
-    def __process_manager(self, producer_callback, producer_queue, consumer_queue):
+    def _process_manager(self, producer_callback, producer_queue, consumer_queue):
         """Setup producer processes and manage completion of processes.
 
         Parameters
@@ -122,7 +122,7 @@
         """
 
         try:
-            producer_proc = [mp.Process(target=self.__producer, args=(
+            producer_proc = [mp.Process(target=self._producer, args=(
                 producer_callback, producer_queue, consumer_queue)) for _ in range(self.cpus)]
 
             for p in producer_proc:
@@ -169,7 +169,7 @@
 
         try:
             consumer_queue = mp.Queue()
-            manager_proc = mp.Process(target=self.__process_manager, args=(
+            manager_proc = mp.Process(target=self._process_manager, args=(
                 producer, producer_queue, consumer_queue))
 
             manager_proc.start()
@@ -244,7 +244,7 @@
         data_items = sum(1 for _ in read_seq(seq_file))
         try:
             consumer_queue = mp.Queue()
-            manager_proc = mp.Process(target=self.__process_manager, args=(
+            manager_proc = mp.Process(target=self._process_manager, args=(
                 producer, producer_queue, consumer_queue))
 
             manager_proc.start()
```

## Why the claims hold

- **The mechanism**, reproduced in isolation: a bound private method reduces to
  `(getattr, (obj, '__work'))`, the unmangled name, so unpickling fails. As a
  Process target it exits 1 under `spawn` and `forkserver` and 0 under `fork`;
  a single-underscore method exits 0 under all three.
- **"The only private methods passed to `mp.Process`"**: every other `target=`
  in the package (`pplacer.py`, `prodigal.py`, `fastani.py`, `skani.py`,
  `hmm_aligner.py`, `tigrfam_search.py`, `pfam_search.py`) is a
  single-underscore method, and the two `mp.Pool` calls in `pipeline/align.py`
  map module-level functions.
- **#619** is macOS Ventura (Intel); **#683** is Python 3.14 with
  `forkserver.py` in the traceback, closed by bioconda-recipes#61837 capping
  Python below 3.14.

## What is left

- **Re-checked 2026-10-24, against GTDB-Tk 2.8.0** (released 2026-10-13, the
  day after this was drafted). The macOS run was repeated on 2.8.0 with Python
  3.13.17: as released, the same `AttributeError` and a hang at "Running
  Prodigal", killed at 180 s; renamed, exit 0 in 5.2 s, pplacer alpha23 step
  included, summary byte-identical to the `fork` run (5.9 s). The Linux run
  was repeated on 2.8.0 the same day (GenomeDK job 5869055, Python 3.14.9,
  `solve-2.8.0/mock/py314-280/`, remote copy under
  `pplacer-compare/mock/py314-280/`): as released, the `AttributeError` and a
  hang, killed at 600 s; renamed, exit 0 in 30 s (2.7.2: 24 s), same
  classification. Both bullets in the text are now 2.8.0's.
- **Scope.** Tested on `classify_wf` only, which covers identify, align and
  classify with pplacer. Other subcommands (`de_novo_wf`, `infer`, …) were not
  run under `spawn`; the text claims only `classify_wf`.
- **The unit test, added 2026-10-27** because CONTRIBUTING.md asks for tests
  "where applicable". It sums `abs` over `[-1, -2, 3]` with two workers, one
  subprocess per start method, 60 s timeout. Measured: with the rename, OK on
  the laptop (macOS, Python 3.14.3, 0.3–0.7 s) and on thylakoid (Linux x86_64,
  Python 3.11.2 and 3.11.16, 0.2 s); against `staging` as is, `spawn` and
  `forkserver` time out and `fork` passes, on both machines (120 s). No orphaned
  processes after the timeouts. Written to run on 3.6, which their
  `__python_requires__` still admits (no `capture_output`/`text`); not run on
  anything older than 3.11. Their CI runs no tests — only the docs build.
- **CONTRIBUTING.md also asks for an issue first.** None opened: #619 and #683
  already describe the bug, and the PR cites both. Both are closed.
- **Not done by this PR:** lifting the recipe's `python <3.14` cap. That is
  theirs, once a release carries the fix.
- **After sending:** record the number in [README.md](README.md), in
  [../STATUS.md](../STATUS.md), and in place of `#NNN` in the pin draft.
