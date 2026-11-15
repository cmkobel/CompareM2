# Draft: carry CheckM2#164 as a bioconda recipe patch

**Unsent.** Destination: a PR against
[bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
`recipes/checkm2/`, from a branch on `cmkobel/bioconda-recipes`, pushed over
SSH as for [pplacer-recipe-version-pr.md](pplacer-recipe-version-pr.md). This is
the first option in DECISIONS.md, *2026-11-10 — CheckM2 on macOS: what to do is
paused*. Drafted 2026-11-11; sending it settles part of that pause, so it is
Carl's call.

## The sequence

1. **Comment on [#164](https://github.com/chklovski/CheckM2/pull/164) first**
   (text below). Bioconda's patching guideline asks for evidence of an effort
   upstream. This gives #164's author and chklovski our measurements, and the
   bioconda PR can link to it.
2. **Open the bioconda PR**: `meta.yaml` changed as below, plus two new files,
   `mp-spawn.patch` (#164 verbatim) and `run_test.py`.
3. **When CI is green**, comment `@BiocondaBot please add label`, as for pplacer.
4. **After merge**: check the published build on all three platforms, then
   comment on [#97](https://github.com/chklovski/CheckM2/issues/97), where users
   are still downgrading Python to get round this, to say that bioconda's build 3
   fixes it. Then the CompareM2 follow-ups at the end of this file.

## Review

chklovski is the recipe's only listed maintainer and is not a member of the
bioconda org, so a bioconda member reviews and merges. The PR body below
cc's chklovski as the recipe maintainer and wwood as a CheckM2 author.
**m-nikolaidis** wrote #164; credit them in the body.

## What is true, and how it was checked (2026-11-11, the laptop under Asahi, linux-aarch64)

- **The recipe** (`meta.yaml` blob `5b9f986`, last changed by #69888 on
  2026-10-04) builds `noarch: python` from the 1.1.0 GitHub tarball at build 2,
  with no patches. linux-64, linux-aarch64 and osx-arm64 all install that one
  build, so one patch reaches all three.
- **#164**: head `3f45368` on `m-nikolaidis/CheckM2:macos_arm64_issue`, one
  file, +15 −13. The four `mp.Process` targets lose one leading underscore, and
  five `== None` become `is None`. `gh pr diff 164` has sha256
  `57ac9762…dcdb73`. `main` is two commits past 1.1.0, both in `main.py`, so
  `predictQuality.py` is the same in both. Those four are the only name-mangled
  `mp.Process` targets in CheckM2.
- **conda-build 26.9.1 builds the patched recipe and passes its tests**
  (`--override-channels -c conda-forge -c bioconda`, 3 min 26 s, 2.1 GB peak).
  conda-build applied the patch at level 1 with no offsets and no fuzz, and
  reported it reversible and minimal (`RA-MD1LOVE`). `checkm2 --help` and
  `run_test.py` pass, and the test prints `3 protein files, metadata for 3`.
  **Without the `patches:` entry the same recipe fails its test**, with #97's
  `AttributeError` for `__set_up_prodigal_thread` and `__reportProgress`. So the
  test checks the fix rather than the install. Work directory:
  `~/.cache/cm2-checkm2-build/`, deletable.
- **The real entry point, under `spawn`.** `bin/checkm2`, unmodified, with a
  `sitecustomize.py` on `PYTHONPATH` that sets `spawn` in every interpreter,
  including the workers. This reproduces the macOS failure on Linux exactly.

  | | as released (build 2) | #164 |
  | --- | --- | --- |
  | `checkm2 testrun` | exit 1, 25 s: `No protein files were generated`, #97 word for word | exit 0, 115 s: `Test run successful!` |
  | `checkm2 predict`, four genomes | exit 1, 26 s | exit 0, 116 s, report byte-identical to the 11-10 run |

  #104's `module '__main__' has no attribute '__spec__'` does not occur. That
  report came from an egg install run through `pkg_resources.run_script`, which
  the conda package's entry point does not use.
- **Results do not move.** `predict` over the four `tests/E._faecium` genomes,
  run through a wrapper that sets the start method: the patched copy's
  `quality_report.tsv` is byte-identical under `spawn` (82 s), `forkserver`
  (77 s) and `fork` (75 s) to the pipeline's unpatched run under `fork` on
  2026-11-10. Completeness 100.0, contamination 0.47 / 0.47 / 3.42 / 1.54.
- **Linux is not exposed today, and why.** checkm2 1.1.0's own
  `tensorflow >=2.17,<2.20` has no build for Python 3.13 or 3.14 on any of the
  three platforms. `pixi lock` finds no solution for `python 3.13.*` or
  `3.14.*`, and solves `3.12.*`. Linux forks below 3.14, so the bug reaches
  Linux the day the TensorFlow range admits 3.14.
- **macOS arm64**, 2026-11-10 (STATUS.md, *macOS arm64, end to end*): the same
  four renames, though not #164 byte for byte, gave the same numbers.
  *Not run:* #164 or the built package on macOS. See *Before sending*.
- **Demand**: #97 (2024-02-27) has five more users reporting the same error,
  the latest 2026-03-10. The workarounds posted there are Python 3.7 with older
  numpy, pandas and scipy, and a script that rewrites the installed file.

## Comment for #164 (step 1)

````markdown
Thanks for this. I've tested it as a bioconda recipe patch on 1.1.0, and plan to carry it there until it's in a release, so that conda users on macOS get it now.

Forcing the `spawn` start method on Linux reproduces #97 exactly, which makes it testable without a Mac. With this PR's diff applied to 1.1.0:

| | 1.1.0 | with #164 |
| --- | --- | --- |
| `checkm2 testrun` under `spawn` | `No protein files were generated` | `Test run successful!` |
| `predict` on 4 genomes under `spawn` / `forkserver` / `fork` | fails / fails / passes | `quality_report.tsv` byte-identical to 1.1.0 under `fork`, all three |

This matters on Linux too, eventually: Python 3.14 makes `forkserver` the default there.
````

## The PR as drafted

**Title:** `Update checkm2: patch multiprocessing so predict runs on macOS`

**Body:**

````markdown
CheckM2 1.1.0's `predict` and `testrun` fail on macOS. `predictQuality.py` passes four name-mangled methods (`__set_up_prodigal_thread` and three others) to `mp.Process`. Under `spawn`, the macOS default, the worker processes can't look them up by name, every worker dies with `AttributeError: 'Predictor' object has no attribute '__set_up_prodigal_thread'`, and the run ends at `No protein files were generated`. This is chklovski/CheckM2#97, open since 2024-02. Five more users have reported it there; the workarounds are downgrading to Python 3.7 or editing the installed file.

The fix is chklovski/CheckM2#164 by @m-nikolaidis, not yet reviewed upstream. This carries its diff unchanged as `mp-spawn.patch` and bumps the build number. The patch can come out in the first release that includes #164. The package is noarch, so this one build fixes osx-arm64. Linux isn't affected today because it forks, but Python 3.14 makes `forkserver` the Linux default, which fails the same way.

`run_test.py` runs the two multiprocessing stages, gene calling and metadata, under `spawn` on the bundled test genomes, without the DIAMOND database. Built locally with conda-build, it passes with the patch and fails with the error above without it.

Results are unchanged. On four *Enterococcus faecium* genomes, the patched package gives a `quality_report.tsv` byte-identical to 1.1.0's under `spawn`, `forkserver` and `fork`, and `checkm2 testrun` passes under `spawn`.

@chklovski as the recipe maintainer, and @wwood as a CheckM2 author.
````

## The diff

Against `recipes/checkm2/` as of 2026-11-11 (`meta.yaml` blob `5b9f986`), plus
two new files. `mp-spawn.patch` is `gh pr diff 164 --repo chklovski/CheckM2` at
`3f45368`, verbatim. `run_test.py` is below, exactly as built.

```diff
--- a/recipes/checkm2/meta.yaml
+++ b/recipes/checkm2/meta.yaml
@@ -9,9 +9,13 @@
 source:
   url: https://github.com/chklovski/{{ name }}/archive/{{ version }}.tar.gz
   sha256: {{ sha256 }}
+  patches:
+    # chklovski/CheckM2#164, for #97: name-mangled mp.Process targets fail
+    # under spawn (macOS) and forkserver (Linux, Python >=3.14).
+    - mp-spawn.patch
 
 build:
-  number: 2
+  number: 3
   noarch: python
   script:
     - sed -i.bak 's|find_packages|find_namespace_packages|g' setup.py
```

### `run_test.py`

conda-build runs it on its own. It calls two private methods, so an upstream
refactor can break it, and it comes out with the patch. In the build log,
conda-build's own `running run_test.py` line repeats once per worker, because
`spawn` re-imports the file in each worker. That is expected, not a fault.

```python
"""Gene calling and metadata under spawn: the two stages chklovski/CheckM2#97
broke, both of which run before DIAMOND opens the database."""
import multiprocessing as mp
import os
import sys
import tempfile

if __name__ == "__main__":
    mp.set_start_method("spawn")
    from checkm2.defaultValues import DefaultValues
    from checkm2.predictQuality import Predictor

    genomes = DefaultValues.TESTRUN_GENOMES
    stand_in = os.path.join(genomes, sorted(os.listdir(genomes))[0])  # existence is all it checks
    p = Predictor(genomes, tempfile.mkdtemp(), ".tst", 2, False, stand_in)
    p._Predictor__run_prodigal(None)
    faa = sorted(os.path.join(p.prodigal_folder, f) for f in os.listdir(p.prodigal_folder))
    meta = p._Predictor__calculate_metadata(faa)
    print(f"{len(faa)} protein files, metadata for {len(meta)}")
    sys.exit(0 if len(faa) == len(meta) == 3 else 1)
```

## Before sending

- #164 still at `3f45368` and the recipe still at blob `5b9f986`. If #164 has
  moved, re-take the diff and rebuild.
- **On macOS, build the recipe natively**: the same `conda-build` command, run
  there, applies the patch and runs `run_test.py` under the platform's own
  `spawn`. This is the one platform the PR is for that it has not yet run on, and
  a reviewer may ask. The recipe directory to rebuild from is three files: the
  patched `meta.yaml`, `mp-spawn.patch` and `run_test.py`. Optionally also
  `checkm2 testrun` from that build.
- The comment on #164 is posted, and the PR body links to it.

## After merge

- `pixi exec --spec 'checkm2==1.1.0=*_3' -- checkm2 --help` resolves on
  linux-64, linux-aarch64 and osx-arm64, and `run_test.py` passes in it, on the
  laptop under Asahi and under macOS.
- Comment on #97 with the fixed build.

## What it commits CompareM2 to, once merged

- **A floor on the patched build** in `CHECKM2_ENV`, such as
  `bioconda::checkm2[version='>=1.1.0',build_number='>=3']`. Changing the spec
  is also what makes an existing install rebuild instead of keeping build 2.
  *Unverified:* that `pixi lock` and `render_envs` accept the bracket form. The
  floor test looks only for `>=`.
- **The same probe in `check.py`**, beside the zero-vector check, which never
  reaches gene calling. It fails on build 2 on every platform, so it goes in with
  the floor, not before.
- The macOS refusal is pipeline-wide (`cli.py`, `unsupported_platform`) and is
  not settled by this alone.
