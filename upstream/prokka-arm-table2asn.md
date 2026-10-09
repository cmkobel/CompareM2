# prokka on Apple silicon: tbl2asn → table2asn

**Not sent; nothing drafted against a recipe yet.** Measured 2026-09-28 on the
laptop, macOS 26.7 / Apple silicon. Working material, outside git:
`~/postdoc/cm2-macos/prokka/` — the shim, both lockfiles, every run output.

This is the fallback for panaroo on arm, not the plan. panaroo never invokes
prokka ([bioconda-panaroo-pr.md](bioconda-panaroo-pr.md)), so deleting one line
from panaroo's recipe gets the same install that this route needs two upstream
changes to reach. It matters only if that line is refused — and panaroo's
author approved removing it on 2026-10-04; the PR is [bioconda-recipes#69901](https://github.com/bioconda/bioconda-recipes/pull/69901), sent 2026-10-05.

## What blocks prokka, and what does not

prokka 1.15.6 is `noarch: generic` perl. Its bioconda `run:` deps, **minus
`tbl2asn-forever`**, install natively on `osx-arm64`: 222 packages, perl
5.32.1, blast 2.17.0, hmmer 3.4, infernal 1.1.5, aragorn 1.2.41, barrnap 1.10.6,
prodigal 2.6.3, minced 0.4.2 (which brings openjdk 25). `prokka --setupdb` and
`--depends` pass.

`tbl2asn-forever` cannot get an arm build: it repackages NCBI's last prebuilt
tbl2asn 25.7 (x86_64 on mac) and date-fakes it with libfaketime via
`DYLD_INSERT_LIBRARIES`. NCBI no longer distributes tbl2asn at all — its FTP
readme says so and points to `table2asn`. Aside: the **linux-aarch64**
`tbl2asn-forever` build ships an x86-64 ELF (`file bin/real-tbl2asn`), so it
installs on aarch64 but cannot run there natively. Not executed; see below for
why its test passes anyway.

`table2asn` from NCBI has been **arm64 on mac since 2025-03-03**: the 2024-10-23
`mac.table2asn` is x86_64, 2025-03-03 and 2025-09-16 are arm64. The current one
is 1.29.324 and links only `/usr/lib` system dylibs (libiconv, libbz2, libz,
libxml2, libxslt, libexslt, libresolv, libc++, libSystem). bioconda's
`table2asn` recipe is 1.28.1179 (2024-10-23), linux-64 and osx-64 only; the one
bump to 2025-03-03 (bioconda-recipes#56474) was closed unmerged with no comment.
That the arm binary is why is an inference.

prokka upstream will not fix this: 1.15.6 is the declared last release, and
issues #666, #669, #678, #702 about tbl2asn/table2asn are open with no PR.

## Flag translation

prokka calls tbl2asn at `bin/prokka:1429` and probes its version with
`tbl2asn - | grep '^tbl2asn'` against `MINVER 24.3` (line 171).

| tbl2asn (prokka passes) | table2asn |
| --- | --- |
| `-a r10k` — runs of 10+ N are gaps, known length | `-gaps-min 10` — no `rNk` form exists |
| `-M n` | `-M n` |
| `-M b` — prokka's choice above 10,000 contigs | **none**; the shim drops it. Untested |
| `-Z <file>` | `-Z` is a switch; the report goes to `<prefix>.dr` |
| `-V b`, `-l paired-ends`, `-N`, `-y`, `-i` | same |
| exit 0 with validator errors | **exit 2** with validator errors, outputs still written |

The exit code is the one that bites: prokka's `runcmd` dies on non-zero, so an
unshimmed table2asn fails every genome at the last step. NCBI's readme does not
document the code; "2 means validator errors" is inferred from it tracking the
`.val` file across `-Z` / `-V b` variants. The shim maps 2 to 0 and relies on
prokka's next step, which reads the `.gbf`, to fail a run that really failed.

## Executed

Shim at `shim/tbl2asn` (sha256 of the NCBI binary beside it `94f276a1…`).
`tests/E._faecium/116_2.fna` and `116_2 duplicate.fna`, `--cpus 8`:

- **exit 0 both, 42.8 s and 43.4 s wall.** 4 contigs, 2,585 CDS, 67 tRNA,
  18 rRNA, 1 tmRNA each.
- Duplicate pair **identical** in `.faa`, `.ffn`, `.gff`, `.tbl` and `.gbk`
  after normalising the prefix and date lines. Locus-tag prefix `JEENLPHM` in
  both.

Against the real thing: tbl2asn-forever 25.7.2f from `osx-64`, **under
Rosetta**, on the same `.fsa`/`.tbl` pair:

- Same **2,675 features at identical locations**; **2,608** with identical
  qualifiers.
- The other **67 are every tRNA**: tbl2asn turns `product tRNA-Leu(gag)` into
  `/product="tRNA-Leu"` plus `/note="tRNA-Leu(gag)"`; table2asn keeps the
  product and **drops the anticodon**. It survives in `.tbl`, `.gff`, `.tsv`.
  *Why*, from the `.sqn`: both tools store only `ext tRNA { aa ncbieaa 76 }`
  (Leu) — neither turns `(gag)` into structured data. tbl2asn also keeps the
  raw string as the feature's `comment`, which the flatfile prints as `/note`.
  table2asn does not: `x_TrnaToAaString` in the C++ toolkit's
  `src/objtools/readers/readfeat.cpp` strips `tRNA-`, cuts at the first of
  `-,;:()=\'_~`, and the tRNA `product` case then consumes the qualifier. It
  keeps the original only for `fMet`, `iMet` and `Ile2`. Read on `master`
  2026-09-28, not at the 1.29.324 tag. NCBI's model has a place for this —
  an `anticodon` qualifier carrying a *location* — but prokka builds the
  product as `$x[1].$x[4]` (`bin/prokka:527`) and discards Aragorn's anticodon
  offset `$x[3]`, so emitting one would need a prokka patch. The shim could
  restore tbl2asn's `/note` by adding a `note` line to the `.tbl` instead.
  Untested.
- table2asn **adds a `gene` feature per CDS**, 2,585 of them; tbl2asn emits none.
- `.val`: tbl2asn 4,141 `GeneXrefWithoutGene` errors, table2asn 86. Advisory
  either way.

Nothing downstream of the `.gbk` was run, and `.sqn` was not compared.

## What it would take

Step one is the same whatever comes after: **bioconda's `table2asn` needs an
`osx-arm64` build.** Bump the recipe to 1.29.324 (2025-09-16), add
`extra: additional-platforms: [osx-arm64]`, and `skip: True  # [osx and not
arm64]`. NCBI ships no x86 mac binary after 2024-10-23, so osx-64 keeps the
1.28.1179 builds already on the channel. The recipe's `install_name_tool`
libbz2 line should be a no-op on arm, where the binary already names
`/usr/lib/libbz2.1.0.dylib`. Not built. The recipe is maintained in practice by
mencian, who opened the one bump past 2024-10-23 (#56474) and closed it 32
minutes later. That the arm binary was the reason is an inference.

**table2asn does not expire the way tbl2asn did**, which is what makes this
possible at all. The 2024-10-23 binary, 23 months old today, run under Rosetta
on the same input: it prints "more than 1 year old", exits 2 like 1.29.324, and
writes a `.gbf` with the same 2,585 CDS. No faketime needed.

Step two can go in one of three places:

| | change | prokka recipe | Linux / osx-64 output | persuade |
| --- | --- | --- | --- | --- |
| **A** | `tbl2asn-forever` gets an `osx-arm64` build that installs the shim and `run: - table2asn  # [osx and arm64]` | **untouched** — 1.15.6's `tbl2asn-forever >=25.7` is satisfied | **unchanged** | bioconda reviewers only; the recipe is already per-platform, so selectors work |
| B | prokka drops `noarch`, and uses table2asn plus a `bin/prokka` patch on arm only | rebuilt per platform, `skip-lints: should_be_noarch_generic` | unchanged | reviewers of a heavily used recipe |
| C | prokka swaps to `table2asn` everywhere plus the patch | `noarch` kept | **changes**: tRNA anticodon note lost, `gene` features added, linux-aarch64 no longer installable | same, with a behaviour change to defend |

**A is the smallest**, and it also makes the panaroo recipe change
unnecessary: stock prokka, and therefore stock panaroo, would install on arm.
Its cost is honesty in naming. On one platform a package called
`tbl2asn-forever` 25.7.2f would really be `table2asn` 1.29.324 behind a
translator, and a reviewer may reasonably balk at that. The shim also has to
map `--help` to `-help`: `table2asn --help` exits 1. Its recipe test
`tbl2asn-test` only checks that the output does *not* contain "more than a year
old" — which, reading the script, is also why the linux-aarch64 build passes
with an x86 binary that cannot run there. On arm the test should run the shim
against a real input instead.

Nothing in A, B or C has been built as a recipe. What is verified is the shim
and the binaries, above. Proving A the way intbitset was proved means building
both recipes locally, then installing **unmodified** bioconda prokka and panaroo
on `osx-arm64` against that channel, and running both.

**Rosetta is not a route.** An all-`osx-64` prokka environment works today —
tbl2asn-forever ran above — but Apple's developer notice makes macOS 27 the
last release that runs Intel apps generally
(https://developer.apple.com/news/?id=w5ngl9k2).
