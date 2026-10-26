# Draft: `aS=AS` typo in `panaroo/cdhit.py` makes arm64 skip gene-family collapsing

**Sent 2026-10-20 as [gtonkinhill/panaroo#382](https://github.com/gtonkinhill/panaroo/pull/382)**, verbatim from [panaroo-cdhit-aS-pr.md](panaroo-cdhit-aS-pr.md), from `cmkobel:fix-aS-typo` (commit `963bac6`). Destination: a PR against
[gtonkinhill/panaroo](https://github.com/gtonkinhill/panaroo), `panaroo/cdhit.py`,
from a fork branch. Measured 2026-10-19 and 10-20; the full record is in
`STATUS.md`, *panaroo installs on arm, and its cleaning stage does nothing
there*.

Why it is urgent rather than cosmetic: bioconda-recipes#69901, which the
author approved, made panaroo installable on Apple silicon on 2026-10-19, and
every arm install since then gives a different pangenome from x86 without any
warning.

## The change

Three one-character edits, all in `panaroo/cdhit.py` on master (checked
2026-10-20):

```diff
@@ def align_dna_cdhit  (line 333)
-                  aS=AS,
+                  aS=aS,
@@ def iterative_cdhit  (line 417, dna branch)
-                          aS=AS,
+                          aS=aS,
@@ def iterative_cdhit  (line 431, protein branch)
-                      aS=AS,
+                      aS=aS,
```

## It was reported once, and half-fixed

The typo was in all five call sites from panaroo's first visible commit
(`9ee331b`, 2019-07-02). `iterative_cdhit` (`bba6191`, 2019-07-21) inherited it.
Every call passes `aL=aL, AL=AL, aS=AS` and no `AS=` at all, which reads like
`aS=aS, AS=AS` collapsed into one argument. That is a guess; the history
before 2019-07-02 is not in the repository.

[#303](https://github.com/gtonkinhill/panaroo/issues/303) (2024-08-02) reported
it in `cluster_nodes_cdhit` only. The author confirmed the typo on 2024-08-16
but judged it harmless, *"as Panaroo only uses global alignment
(use_local=False), which overrides this parameter"*. `1662aed` (2024-09-20)
fixed the two reported lines and left the other three.

**The "global overrides it" reasoning does not hold for CD-HIT 4.8.1.**
`ControlShortCoverage` enables coverage control on `short_coverage > 0.0`
alone (`cdhit-common.c++:2698`). The rejection at 3375–3377 is gated only on
that flag. `global_identity` changes nothing but `len_eff` (2704). Every
measurement below ran in global mode (`use_local=False`, so no `-G 0`), and
`-aS` changed the output on both platforms. On x86 the change is small enough
to look harmless, which may be why the reasoning seemed right. The PR has to
answer that reasoning, or it reads as a re-report of a closed issue.

Lines 193 and 206 (`cluster_nodes_cdhit`) already read `aS=aS`. Neither caller
of `iterative_cdhit` (`clean_network.py:114, 126`) passes `aS` or `AS`, so after
the fix CD-HIT gets the default `-aS 0.0`. `align_dna_cdhit` is imported by
`find_missing.py` and never called, so line 333 changes no output today.

## The PR text

Kept separately, ready to paste, in [panaroo-cdhit-aS-pr.md](panaroo-cdhit-aS-pr.md):
commit message, title and body. This file is the working note behind it.

## What follows once it is merged

1. **A panaroo release.** bioconda's autobump bot picks it up; the recipe no
   longer lists `prokka`, so arm installs keep working.
2. **If the release is slow:** a `patches:` entry on bioconda's 1.8.0 recipe
   (build 2) would fix arm installs sooner. Only with the author's agreement,
   and only after the upstream PR exists.
3. **CompareM2:** no change is needed for correctness, since it runs on
   linux-64 only. But its results move by the x86 delta above, so raise
   `PERL_ENV`'s floor to the fixed version and record the change in
   `STATUS.md`. Then a report never mixes the two behaviours.

Not proposed: a guard in CD-HIT rejecting `-aS > 1`. That would have made this
fail loudly, but CD-HIT's last release was 4.8.1 in 2019, and the panaroo
change alone is sufficient.
