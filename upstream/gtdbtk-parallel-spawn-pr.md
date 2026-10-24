# Draft: make GTDB-Tk's `Parallel` work under `spawn` and `forkserver`

**Unsent.** Destination: a PR against
[Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk), `master`, file
`gtdbtk/biolib_lite/parallel.py`, from a fork branch. **Send this one first**:
it is the easy yes (it fixes Python 3.14 for their Linux users and asks nothing
of them), and [gtdbtk-pplacer-pin-pr.md](gtdbtk-pplacer-pin-pr.md) cites its
number.

Found 2026-10-12 while checking the pin PR's macOS claim: with the pin lifted,
gtdbtk installs on Apple silicon, then hangs at "Running Prodigal" on any genome
the ANI screen does not classify. The ANI-screen-only path is not affected
(identify is skipped). Evidence:
[gtdbtk-pplacer-evidence.md](gtdbtk-pplacer-evidence.md), *Result 4*.

## The PR as drafted

**Title:** `Make Parallel work under spawn and forkserver (macOS, Python 3.14)`

**Body:**

````markdown
`biolib_lite.Parallel` passes two private methods, `self.__producer` and `self.__process_manager`, as `mp.Process` targets. Under `fork` that works, but under `spawn` and `forkserver` the target has to be pickled, and a bound method is pickled by its unmangled name. The child can't find it and dies with

```
AttributeError: 'Parallel' object has no attribute '__process_manager'
```

while the parent waits for it forever. macOS defaults to `spawn`, and Python 3.14 changed the Linux default to `forkserver`, so this is the error in #619 and #683, and the reason the bioconda recipe caps Python below 3.14 (bioconda/bioconda-recipes#61837).

This renames the two methods to `_producer` and `_process_manager`. Nothing else in the package refers to them, and they are the only private methods passed to `mp.Process`.

Tested on the r232 mock database (`classify_wf -x gz --place_species -f`):

- Linux, GTDB-Tk 2.7.2, Python 3.14.8 (`forkserver`): as released, the `AttributeError` and then a hang, killed after 600 s. With this change, exit 0 in 24 s.
- macOS on Apple silicon, GTDB-Tk 2.8.0, Python 3.13 (`spawn`), with pplacer 1.1.alpha23: as released, it hangs at "Running Prodigal". With this change, exit 0, with a summary byte-identical to the same run forced to `fork`.

Best, Carl
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
  included, summary byte-identical to the `fork` run (5.9 s). **The Linux
  Python 3.14 run is still 2.7.2's** — the 2.8.0 job is prepared
  (`solve-2.8.0/mock/py314-280/`) and needs a GenomeDK login. If it runs,
  change the first bullet to 2.8.0 and its numbers; if not, the text is
  accurate as it stands.
- **Scope.** Tested on `classify_wf` only, which covers identify, align and
  classify with pplacer. Other subcommands (`de_novo_wf`, `infer`, …) were not
  run under `spawn`; the text claims only `classify_wf`.
- **Not done by this PR:** lifting the recipe's `python <3.14` cap. That is
  theirs, once a release carries the fix.
- **After sending:** record the number in [README.md](README.md), in
  [../STATUS.md](../STATUS.md), and in place of `#NNN` in the pin draft.
