# Draft: drop `prokka` from panaroo's run dependencies

**Sent 2026-10-05 as [bioconda-recipes#69901](https://github.com/bioconda/bioconda-recipes/pull/69901)**, after the author approved by email 2026-10-04
(below). The text was reworded before sending — it leads with Apple silicon and
drops the email paragraph — so the PR, not the draft below, is the record of
what was said. Destination: a PR against
[bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes),
`recipes/panaroo/meta.yaml`, from a fork branch.

Measured 2026-09-03 on macOS 26.2 / Apple silicon, panaroo 1.8.0; re-run
2026-09-28 on macOS 26.7 against the published intbitset. This was the second of
two changes panaroo needs on arm. The first,
[intbitset-feedstock-pr.md](intbitset-feedstock-pr.md), **merged 2026-09-25**,
so this line is now the only blocker.

## The PR as drafted (the sent text differs; see #69901)

Everything below this section is the working note it was written from.

**Title:** `Update panaroo: drop prokka from run requirements`

**Body:**

> panaroo's recipe lists `prokka` as a run requirement, but panaroo does not
> need it:
>
> - upstream's `setup.py` does not declare it;
> - `panaroo/prokka.py` is a GFF parser and runs nothing;
> - the only call to the prokka binary is in the separate `run_prokka` helper
>   (`panaroo/run_prokka.py`).
>
> That one line is what keeps panaroo off Apple silicon. prokka requires
> `tbl2asn-forever`, which repackages NCBI's retired x86 tbl2asn and has no
> `osx-arm64` build. intbitset, the other blocker, has had `osx-arm64` builds on
> conda-forge since conda-forge/intbitset-feedstock#21.
>
> **Cost:** `conda install panaroo` no longer brings prokka on any platform, so
> `run_prokka` users need `conda install prokka` alongside. All eight test
> commands, `run_prokka --help` included, pass without prokka on PATH.
>
> @gtonkinhill — as agreed by email. You added prokka here in #21341; this
> removes it rather than making it optional. If you'd rather keep a version
> floor for `run_prokka` users, `run_constrained: - prokka >=1.14` does that
> without installing it; happy to switch.
>
> **Tested** on osx-arm64 (macOS 26, Apple silicon), this recipe built locally
> and everything else from conda-forge and bioconda: panaroo 1.8.0 on three
> *E. faecium* genomes, one a byte-identical copy of another, exits 0 with
> 3,569 gene clusters; the copy pair is identical in every cluster and 0 SNPs
> apart in the core alignment.

Every claim in it was checked 2026-09-28. The three facts about panaroo are
from v1.8.0's source, 2026-09-03. The eight test commands plus
`import panaroo` exit 0 with neither `prokka` nor `tbl2asn` on PATH. The run is
the re-run in `../STATUS.md`: 141 packages, 191.66 s, and
`gene_presence_absence.Rtab` byte-identical to 09-03. That the copy pair differs
in 0 of 3,569 clusters is from 09-03; the identical Rtab carries it forward.

After CI passes, bioconda wants `@BiocondaBot please add label` as a comment.
CI builds `noarch` on linux-64 only, so the arm result above is ours, not CI's.

## The diff

```diff
 build:
-  number: 0
+  number: 1
   noarch: python
@@
     - mash
-    - prokka
     - cd-hit
```

## Why

`prokka` requires `tbl2asn-forever >=25.7`, which bioconda ships for
`linux-64`, `linux-aarch64` and `osx-64`. It repackages NCBI's **prebuilt**
tbl2asn binary, so there is nothing to rebuild for Apple silicon; bioconda's
`table2asn` has no arm build either. (NCBI's own `table2asn` is arm64 on mac,
and prokka runs on arm behind a shim — but packaging that is two further PRs;
see [prokka-arm-table2asn.md](prokka-arm-table2asn.md).) That one dependency makes panaroo uninstallable on
`osx-arm64`, and because panaroo's core gene alignment is the input to both
`snp-dists` and `fasttree`, it costs three tools rather than one.

## Why removing it is correct, not a workaround

Four independent reasons, in descending order of how much they should count:

1. **Upstream does not declare it.** panaroo 1.8.0's `setup.py`
   `install_requires` is `networkx, gffutils, BioPython, joblib, tqdm, edlib,
   scipy, numpy, matplotlib, scikit-learn, plotly, dendropy, intbitset,
   biocode`. No prokka. The recipe's `prokka` line is a bioconda addition, not
   an upstream requirement.
2. **`panaroo` never invokes the binary.** `panaroo/prokka.py` is a GFF
   *parser* — `process_prokka_input`, no `subprocess` — imported by
   `__main__.py:12` and `integrate.py:14`. The only shell-out in the package is
   `panaroo/run_prokka.py:136`, which belongs to the separate `run_prokka`
   console script: a convenience wrapper for *producing* the GFFs panaroo then
   reads.
3. **The recipe's own tests still pass.** All eight `test: commands:`,
   `run_prokka --help` among them, exit 0 in an environment with no prokka on
   PATH. Verified 2026-09-03 in a 141-package `osx-arm64` environment built
   without it.
4. **prokka is end-of-life.** Release 1.15.6 (2025-12-14) is upstream's
   declared last release, and its notes recommend bakta. A live tool is being
   held off a platform by an archived one.

## What it costs, stated plainly

The recipe is `noarch: python`, so this applies on every platform, not just
arm: a `conda install panaroo` on Linux stops pulling prokka too. Anyone using
the `run_prokka` helper then has to `conda install prokka` alongside — which
still works on Linux and `osx-64`, and which is where that dependency belongs,
since it is a dependency of one optional entry point rather than of panaroo.

**The line is Gerry Tonkin-Hill's own.** He added `prokka`, with `mash` and
`intbitset`, in bioconda-recipes#21341 (2020-04-08, panaroo 1.2.0), and the
recipe lists no `recipe-maintainers` — so a bioconda member could merge this
without him, and should not. He merged intbitset-feedstock#21 himself on
2026-09-25, the day he was asked.

**He agreed, 2026-10-04.** Carl asked by email on 2026-09-28 ("Panaroo on
Macos?"); the reply: *"If you're willing to submit a pull request to remove the
prokka dependency (or make it optional) that would be great!"* So either form
is approved; the diff below is the removal.

If he wants to keep a floor for `run_prokka` users, conda's form of an optional
dependency is `run_constrained: - prokka >=1.14` — it bounds prokka when
installed and never pulls it in. Not tested whether bioconda's linter accepts it
on a `noarch: python` recipe.

`additional-platforms` is not an alternative here. Both panaroo (`noarch:
python`) and prokka (`noarch: generic`) are architecture-independent already;
what is missing is not a build but a dependency closure. Nor can the dependency
be made conditional — selectors do not apply to `noarch` run requirements.

## Evidence that the result works

panaroo 1.8.0 installed into an `osx-arm64` conda environment holding its
conda dependencies minus `intbitset` and `prokka`, with `intbitset` 4.1.2 from
its PyPI arm wheel — python 3.13.15, `Mach-O 64-bit bundle arm64`.

Input: prodigal 2.6.3 GFF3 with the FASTA appended, over CompareM2's
*E. faecium* test set — **2,587 / 2,587 / 3,192** CDS, the third genome being a
byte-identical duplicate of the first.

`panaroo --clean-mode strict -a core -t 4 --remove-invalid-genes` exited 0 in
**180.96 s**:

| | |
| --- | --- |
| gene clusters | 3,569 |
| core (99–100%) | 2,146 |
| shell (15–95%) | 1,423 |
| core alignment | 1,952,790 columns |
| **duplicate pair, clusters differing** | **0 of 3,569** |
| `snp-dists` 1.2.0, duplicate pair | **0** SNPs (4,111 to E8202) |
| `FastTree` 2.2.0 `-nt -gtr` | duplicate pair at branch length 0.0 |

Every number matches the run recorded in `../STATUS.md` from a hand-built
environment the day before, which used a locally compiled intbitset instead of
the wheel.

`--remove-invalid-genes` is needed because the input is prodigal rather than
bakta, whose partial genes panaroo rejects; it drops E8202 from 3,192 to 3,128.
So these gene counts are **not** comparable to CompareM2's Linux runs, which
annotate with bakta. What this shows is that panaroo runs correctly on arm, not
that the pipeline's numbers reproduce there.

## What is left

- Sent 2026-10-05 from `cmkobel/bioconda-recipes`, branch `panaroo-drop-prokka`,
  commit `ac02471` on bioconda master `7074fb0`. Next: once CI passes, comment
  `@BiocondaBot please add label`. The PR does not @-mention Gerry and the
  recipe lists no maintainers, so GitHub will not notify him.
- Resolved 2026-10-04: the author's approval, by email.
- Not a blocker: a bakta-annotated arm run, so that gene counts are comparable
  to the Linux reference. bakta runs on arm (STATUS.md, 2026-09-04), but no
  bakta GFF has been fed to arm panaroo yet.
- Resolved 2026-09-25: the intbitset PR, merged.
- Optional, separate: a bioconda issue about recipes still pulling EOL prokka.
  Several recipes will have the same problem, and a one-line fix per recipe is
  not the general answer.
