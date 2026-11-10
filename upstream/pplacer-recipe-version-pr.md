# Draft: bioconda's Linux pplacer builds report `--version` as `dev`

**Sent 2026-11-09 as [bioconda-recipes#71400](https://github.com/bioconda/bioconda-recipes/pull/71400)**, with the text below and the
title `Update pplacer: set the version string on Linux` (bioconda's template asks
for "Update" first). Branch `cmkobel/bioconda-recipes:pplacer-version-string`,
one commit `421daae` on upstream `d962af0`, pushed over SSH: the gh token has
no `workflow` scope, and the fork's `master` has diverged from upstream (sync
conflicts), so a branch made through the API from upstream's tip was refused.
Once CI is green, comment `@BiocondaBot please add label`. Destination: a PR against
[bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes),
`recipes/pplacer/build.sh` and `meta.yaml`, from a branch on
`cmkobel/bioconda-recipes`. Send it **before**
[gtdbtk-pplacer-pin-pr.md](gtdbtk-pplacer-pin-pr.md): once it is merged, the
pin PR's paragraph *One thing this makes worse* comes out, because gtdbtk's
`pplacer =1.1.alpha23` does not fix a build number and so picks up the rebuild.

Found 2026-11-09 in the devil's-advocate pass over the pin PR. On #668
pchaumeil read `pplacer version: dev` in the log as the sign of a non-standard
build, and alpha20 and alpha23 both print it.

## What is true, and how it was checked (2026-11-09)

- **pplacer #397 is in alpha23.** GitHub's compare API: `v1.1.alpha23` is one
  commit ahead of #397's merge commit `29f6885`, and `master` is identical to
  the tag. Since #397, `common_src/dune` writes `version.ml` at build time from
  `$PPLACER_VERSION`, else `git -C %{workspace_root} describe`, else `dev`.
  STATUS.md's 2026-10-05 line "merged 2026-09-02 and not yet in a release" was
  true then and stopped being true with alpha23 on 2026-10-07.
- **The recipe never sets it.** `build.sh` (blob `e3ed111`) runs `dune build`
  on the GitHub tag tarball, which has no `.git`. `conda-dynamic-linking.patch`
  touches only the top-level `dune`, not `common_src/dune`.
- **Observed `dev`**: `pixi exec` of bioconda's `pplacer ==1.1.alpha20` and
  `==1.1.alpha23` on this laptop, which is **linux-aarch64**
  (`h9b329fd_0` for alpha23). linux-64 is recorded in the evidence note
  ("the Linux builds print `--version` as `dev`") and is printed again by the
  #668 jobs' `run.sh classify` on GenomeDK.
- **The mechanism works**: a minimal dune project holding alpha23's
  `common_src/dune` verbatim plus a one-line executable, outside any git
  repository, conda-forge `ocaml-dune` and `ocaml 4.14`: unset prints `dev`,
  `PPLACER_VERSION=v1.1.alpha23` prints `v1.1.alpha23`
  (scratchpad `dunetest2/`).
- **osx-arm64 already passes the new test**: the package
  (`pplacer-1.1.alpha23-h4f1019d_0.conda`, upstream's release binary) has the
  string `v1.1.alpha23` in `bin/pplacer` and no `dev` line (`strings`).
- **Not done:** a full bioconda build of the recipe. The added test is the
  check, and bioconda's CI runs it on linux-64, linux-aarch64 and osx-arm64.

## The PR as drafted

**Title:** `Update pplacer: set the version string on Linux`

**Body:**

````markdown
On Linux, bioconda's pplacer 1.1.alpha23 reports its version as `dev`.

Since matsen/pplacer#397, which is in alpha23, `common_src/dune` writes the version at build time from `PPLACER_VERSION`, falling back to `git describe` and then to `dev`. The source tarball has no `.git` and the recipe doesn't set the variable, so the Linux builds get `dev`. The osx-arm64 package is upstream's release binary and already reports `v1.1.alpha23`.

This exports `PPLACER_VERSION="v${PKG_VERSION}"` before `dune build`, bumps the build number, and adds a test that `pplacer --version` prints `v1.1.alpha23`, which runs on all three platforms.

It matters downstream: GTDB-Tk logs `pplacer --version`, and its maintainers read that line to tell pplacer builds apart (on Ecogenomics/GTDBTk#668, `dev` was how alpha20 showed up). gtdbtk pins alpha19 at the moment, and I'm about to ask its maintainers to move to alpha23.

@colinbrislawn, since you did the alpha23 update in #70019.
````

## The diff

Against `recipes/pplacer/` as of 2026-11-09 (`build.sh` blob `e3ed111`,
`meta.yaml` blob `ffd011c`, last changed by #70019 on 2026-10-10).

```diff
--- a/recipes/pplacer/build.sh
+++ b/recipes/pplacer/build.sh
@@ -47,7 +47,9 @@
         && echo "All MCL libraries built successfully!"
 
     cd ..
-    # Build pplacer
+    # Build pplacer. The tarball has no .git, so without this --version
+    # prints "dev" (pplacer's common_src/dune reads it at build time).
+    export PPLACER_VERSION="v${PKG_VERSION}"
     eval $(opam env) \
         && dune build
 
--- a/recipes/pplacer/meta.yaml
+++ b/recipes/pplacer/meta.yaml
@@ -20,7 +20,7 @@
     sha256: d2502e40a95ff06b6fd4c380d682c9727b8413c157a6bfb93a87470f9d6fe9f6  # [osx and arm64]
     
 build:
-  number: 0
+  number: 1
   skip: true  # [osx and x86_64]
   run_exports:
     - {{ pin_subpackage(name , max_pin=None) }}
@@ -45,6 +45,7 @@
 test:
   commands:
     - pplacer --help > /dev/null
+    - pplacer --version | grep -qx "v{{ version }}"
     - guppy to_csv --help > /dev/null
     - rppr check --help > /dev/null
 
```

## CI, 2026-11-09

All checks green on the first push: Lint, Linux Tests, OSX-64 Tests, CircleCI
`build and test (ARM)`. The Linux log shows `export PPLACER_VERSION=v1.1.alpha23`
in the build and the new test `pplacer --version | grep -qx v1.1.alpha23`
passing on `pplacer-1.1.alpha23-h01ede22_1` — the fix verified in a real
bioconda build, not only in the minimal dune project. `@BiocondaBot please add
label` posted the same day.

## Merged, 2026-11-10

Approved and merged by bgruening at 06:49Z ("Danke!"). bioconda published
`pplacer-1.1.alpha23` build 1 for linux-64 (`h01ede22_1`), linux-aarch64
(`h9b329fd_1`) and osx-arm64 (`h4f1019d_1`) within four minutes, and
`pixi exec … pplacer==1.1.alpha23=*_1 -- pplacer --version` prints
`v1.1.alpha23` on GenomeDK (linux-64) and on the laptop (linux-aarch64).

Before approving he commented on the `export` line: "That is a nice fake ;)
The problem is that the installed version then always reports `dev` isn't it?
Can you create an upstream issue that fixes that, please." Done the same
morning as [matsen/pplacer#399](https://github.com/matsen/pplacer/issues/399),
"`--version` reports `dev` when built from a release tarball": the fallback
chain, the bioconda workaround, and two fixes — a `VERSION` file read before
`dev`, or `export-subst` with `$Format:%(describe:tags)$`. The second was
checked on setuptools-scm's GitHub tarball (tag `vcs-versioning-v2.6.0`): the
first `.git_archival.txt` came out expanded, a second copy in the same archive
did not, because git expands one `%(describe)` per archive. Replied in his
thread with the answer (before: always `dev`; now `v1.1.alpha23`) and the
issue link.

## What is left

- **Before sending:** re-check both blobs, and that no other pplacer PR is open.
- **After sending:** the number goes in [README.md](README.md) and here.
- **After it merges:** confirm `pixi exec … pplacer ==1.1.alpha23` prints
  `v1.1.alpha23` on linux-64 and linux-aarch64, then delete *One thing this
  makes worse* from the pin PR.
