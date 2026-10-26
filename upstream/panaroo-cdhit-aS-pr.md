**Sent 2026-10-20 as [gtonkinhill/panaroo#382](https://github.com/gtonkinhill/panaroo/pull/382)**, exactly as below.

**Branch:** `fix-aS-typo` on cmkobel/panaroo, against `master`

**Commit message:**

```
Pass aS instead of AS to cd-hit in iterative_cdhit and align_dna_cdhit

Same typo as the one fixed in cluster_nodes_cdhit for #303. cd-hit was
getting -aS 99999999, which overflows int when cd-hit turns it into a
minimum alignment length. On arm64 that makes cd-hit reject every
alignment, so collapse_families never merges anything there.
```

**Title:** Fix remaining aS=AS typos in cdhit.py (different results on arm64)

**Body:**

Hi Gerry,

This fixes the last three `aS=AS` typos in `cdhit.py`. It's the same one as #303: 1662aed fixed it in `cluster_nodes_cdhit`, but `iterative_cdhit` (both branches) and `align_dna_cdhit` still pass `AS` as `aS`, so the cd-hit calls in the cleaning step run with `-aS 99999999`.

In #303 you mentioned that global alignment overrides this parameter. As far as I can tell that isn't the case in cd-hit 4.8.1. Coverage control is switched on whenever `-aS` is above 0, regardless of `-G` (`cdhit-common.c++`, line 2698), and the check is applied to global alignments too.

On x86 this is mostly harmless, which is probably why nobody has noticed. cd-hit computes `(int)(len * 99999999)`, which overflows. On x86 the result ends up negative, so the coverage filter lets everything through. On arm64 the conversion saturates to INT_MAX instead, so every alignment fails the coverage check, and the iterative clustering in `collapse_families` doesn't merge anything.

I ran into it after the bioconda change (bioconda-recipes#69901) made panaroo installable on Apple silicon. Here are the same four *E. faecium* Bakta GFF3s, panaroo 1.8.0, `--clean-mode strict -a core`:

| | x86-64 (Linux) | arm64 (macOS) |
|---|---|---|
| 1.8.0 | 3,780 total, 2,091 core | 4,060 total, 1,974 core |
| this PR | 3,778 total, 2,092 core | 3,778 total, 2,092 core |

With the fix, both platforms give the same clusters and an identical `gene_presence_absence.Rtab`. I see the same thing on linux-aarch64, so it's the architecture, not macOS.

Running cd-hit directly on the combined proteins at `-c 0.9` shows it more clearly. `-aS 0.0` gives 3,578 clusters on both platforms. `-aS 99999999` gives 3,582 on x86 and 11,174 on arm64, which is close to no clustering at all.

Note that the fix also changes x86 output slightly (2 fewer genes in the example above). Any positive `-aS` turns on cd-hit's coverage branch, and that moved a few short sequences in this data.

`align_dna_cdhit` doesn't seem to be called anywhere at the moment, so that change is just for consistency.

Happy to share the test data or run anything else if it helps.

Carl
