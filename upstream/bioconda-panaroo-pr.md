# Draft: drop `prokka` from panaroo's run dependencies

**Unsent.** Destination: a PR against
[bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes),
`recipes/panaroo/meta.yaml`.

Measured 2026-09-03 on macOS 26.2 / Apple silicon, panaroo 1.8.0. This is the
second of the two changes CompareM2 needs for panaroo on macOS; the other is
[intbitset-feedstock-pr.md](intbitset-feedstock-pr.md). **Neither is sufficient
alone** — with an arm `intbitset` in a local channel and this line still
present, `panaroo>=1.5` still does not solve on `osx-arm64`.

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
tbl2asn binary, so there is nothing to rebuild for Apple silicon; `table2asn`
has no arm build either. That one dependency makes panaroo uninstallable on
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

## What this draft still needs before it is sent

- A bakta-annotated arm run, so the counts *are* comparable to the Linux
  reference. That needs bakta's 1.3 GB light database on the laptop and is the
  natural next step; it is not a blocker for this PR.
- A decision on whether to open this before or alongside the intbitset PR.
  They are independent — neither blocks review of the other — but this one is
  pointless in isolation, so the PR body should link the other and say so.
- Whether to raise the EOL-prokka point as a separate bioconda issue. Several
  recipes will have the same problem, and a one-line fix per recipe is not the
  general answer.
