# Fixing the known bugs before 3.5.0 — plan, 2026-09-24

> **Status, 2026-10-09:** items 0–4 implemented, and the six findings of the
> 2026-09-25 review below are addressed (allowlist narrowed to `_ . / -`;
> `--setup` skips the output check; all three refusal cases assert the
> message; guidance says "from 3.5.0"; the stale sentences fixed). Version
> bumped to 3.5.0. SLURM check on GenomeDK passed: `--demo` 11 of 11, the
> showcase 51 of 51 after one resume, and its report replaces the docs copy
> (STATUS.md, *The 3.5.0 release check, on GenomeDK*). Issue #156's panaroo
> fixes and the re-run defect went in too. Committed as one release commit
> and tagged v3.5.0 on 2026-10-09, on Carl's instruction. Item 5 deferred, as
> Carl decided on 2026-09-25. The draft tag message is at the end.

Written against `7a6f5b9` (master, 18 commits past v3.4.0), after the release
checks recorded in `STATUS.md` under *The 3.5.0 release check*. Each item says
what the fix is, why it is the small one, and what proves it.

The principle throughout: most of these are two pieces of code disagreeing, or
a line that guards against something Snakemake already handles. The fix is to
delete the second opinion, not add a third.

## 0. CarveMe cannot find a solver in a fresh environment — fix first

**Found today, and the published 3.4.0 has it too.** A fresh solve of
`CARVEME_ENV` picks `scip 10.1.0`. conda-forge's `pyscipopt 6.2.1` accepts
that version but links `libscip.so.10.0`, so `import pyscipopt` fails, and
reframed reports `No solver available.` The failure only appears when carveme
runs: 4 of 4 genomes failed. Environments built before the drift carry scip
10.0.3 and still work.

**Fix: one pin.** `CARVEME_ENV = ("bioconda::carveme>=1.6.6",
"conda-forge::scip>=10.0.3,<10.1")`. It has a floor, as the spec test
requires, and it puts the solver build explicitly into the pinned surface,
which `carve_scip.py` already says it is. The 9.8 s / 601 s measurements in
that docstring were taken on 10.0.x, so an upper bound is honest anyway: a
new SCIP minor version should be re-measured before it is let in. **Verified
today in a scratch clone:** the pinned environment solved to 10.0.3, reframed
found `scip`, and carveme + biosynthesis finished 9 of 9 steps with verdicts
equal to 2026-09-11's.

Also: `DECISIONS.md` (why there is an upper bound, and when to lift it), the
environment list in `docs/` if it prints specs, and the fact that this makes
**3.5.0 urgent rather than optional.** Once the fix is tagged, an issue or
release note telling 3.4.0 users to upgrade is worth the minute.

Worth reporting upstream: pyscipopt's conda-forge feedstock should pin `scip`
to the minor version it was built against. That is outward-facing, so it is
Carl's to send.

## 1. A space in `--output` breaks every rule — fix before tagging

**Reproduced today on thylakoid:** `comparem2 --demo --output "~/rc/with space"`
exits 1 with `nothing ran; no report written`. It also creates a stray empty
`~/rc/with/` **outside** the output directory, and `space/logs/` inside it. GNU
`dirname` accepts several operands, so
`mkdir -p $(dirname /…/with space/logs/x.log)` makes two directories.

**Fix — delete, then quote.** Both shell blocks (`snakefile.py:122-123`, tool
rules; `:169-170`, `download_*` rules) open with

```
mkdir -p $(dirname {log})
exec > {log} 2>&1
```

- **Delete the `mkdir` line.** Snakemake creates a log's parent directory before
  the job runs. Checked on 9.26.1 today with a toy rule whose log sat at
  `log dir/sub/a.log`, in a directory whose own path had a space: the log was
  written with no `mkdir` in the rule.
- **`exec > {log:q} 2>&1`.** `:q` is Snakemake's shell quoting, checked in the
  same toy rule. Everything else in the block already goes through
  `shlex.quote`; `{log}` was the one unquoted path.

**Proof.** Add a test that renders the full Snakefile with a workdir containing
a space, gets the formatted shell out of Snakemake (item 2's instrument with
`--printshellcmds`), and asserts that every `exec >` line `shlex.split`s to
exactly three tokens. Then re-run the demo from the fixed tree with
`--output "…/with space"` on thylakoid.

**Not fixed by this, unverified:** a `{`, `}` or `"` in `--output`. `_q()`
wraps paths in bare double quotes, and Snakemake reads braces as format fields,
so both *probably* break. Nobody has tried it. If they do break, one check in
`cli.py` refusing those three characters, with a sentence saying why, is
cheaper than escaping at every site.

## 2. No test builds a DAG from a generated Snakefile — land with item 1

It is item 1's instrument, and it closes the failure CLAUDE.md names. Measured
on 2026-09-22: `prepare()`, then `snakemake --dry-run --cores 1` over two
one-contig FASTAs, builds the full 24-job DAG in 0.75 s. It must be
`prepare()`, not `render()`: `prepare()` writes the GTDB-Tk batchfile the DAG
declares as an input.

- One test, `importorskip("snakemake")`, asserting job counts per rule (genome
  scope ×2, set scope ×1, four downloads). Parametrise the workdir over a plain
  path and one with a space.
- CI: add `snakemake-minimal` to the `pip install` line in
  `.github/workflows/unit.yaml`. It costs about 30 s, and it also un-skips the
  two `runner.run()` tests that CI currently skips.

## 3. snp-dists and FastTree read the unfiltered core alignment — Carl's call

This is a result-changing decision, not a code defect. Panaroo's own
documentation recommends `core_gene_alignment_filtered.aln` for core-genome
phylogenies. On seven *S. aureus* (2026-09-02), switching to it cut **branch
lengths by 26–71% and pairwise SNPs by 20–60%, with the topology unchanged.**

**Recommendation: switch both,** in one commit. Declare the filtered file as a
Panaroo output, point both commands at it, and say in `guidance.py` which
alignment the numbers come from. Two tools reading different alignments would
be worse than either choice. It touches `catalogue.py`, `guidance.py`,
`docs/30 what analyses does it do.md`, `STATUS.md` (both tools go back to
unverified until run) and `DECISIONS.md`. It is not needed for 3.5.0. If you
would rather keep the unfiltered file, the fix is one sentence in
`guidance.py` and a `DECISIONS.md` entry saying it was a choice.

## 4. `_profile_argv` ties conda deployment to a prefix being set — small

`runner.py:203` passes `conda_prefix if deploy else None`, so "no prefix" and
"don't deploy" are the same `None`. The profile branch then drops
`--software-deployment-method conda`, where the API branch keeps it. No user
path reaches this, because `cli.main` always resolves a prefix.

**Fix:** pass `deploy` to `_profile_argv` as its own argument. Then emit
`--software-deployment-method conda` when `deploy` is true, and `--conda-prefix`
when a prefix is given. That matches the API branch, and the docstring sentence
becomes true of both. Add one test of the `deploy=True, conda_prefix=None`
argv.

## 5. The end of a run walks the output tree three times — defer

`settle()`'s `scan()`, `any_outputs_exist()` and `render_report()` each call
`completion()` per tool, and `scan()` runs on the UI thread. **Unmeasured:**
nobody has timed it at scale. The fix (one `{tool: Completion}` map computed in
a worker, shared by all three) is a refactor across `tui.py`, `cli.py` and
`report.py`. Measure first, on a few hundred genomes on a network filesystem;
do not do it for 3.5.0.

## Left as they are

- **AMRFinder's database in the conda prefix.** An accepted cost (27 s
  refetch), not a defect. Today's fresh-prefix dry run showed it as the one
  database "to download".
- **`--tui` against a failing workflow.** This needs a terminal and a person.
  Run `cm2 --demo --tui --output "…/with space"` on the *unfixed* tree: it is a
  ready-made failing workflow.

## Order

0. Item 0, alone. It is the one a user is hitting now.
1. Items 1 + 2, one commit, test first so it fails on `7a6f5b9`.
2. Item 4.
3. Bump to 3.5.0 in `__init__.py`, `pixi.toml`, `citation.cff`,
   `recipe/meta.yaml`.
4. The SLURM release check on GenomeDK. It needs Carl to warm the SSH socket
   (`hpc_login.sh`); non-interactive login is refused.
5. Tag, and do not tag again until the bioconda PR has merged (three days
   after it opens). A newer tag before then overwrites the PR, which is how
   3.1.0 and 3.2.0 never got published.

Item 3 whenever it is decided; item 5 after it is measured.

## Review, 2026-09-25 — before the bump

A review of the uncommitted tree, with three checks run on thylakoid (as
`ghrunner`, scratch under `~/rc/`).

1. **The allowlist admits a comma, and Panaroo refuses one.** A full run with
   `--output 'kørsel+a=b,c:d@e%f'` got 27 of 31 steps. Panaroo's
   `get_gene_sequences` raises on `,` in a GFF path (`if ',' in
   gff_file_name`), because it writes paths into `gene_data.csv`. The 11 tools
   before it accepted every character, but `+ = : @ %` never reached Panaroo,
   snp-dists or FastTree. **The narrower set `[^\w./-]` is verified:**
   `--output 'kørsel-2026.09_x'` ran 31 of 31 in 255 s, with the same SNPs.
   The fix is to narrow the regex, and the test that asserts
   `/data/a+b=c,d@e:f%` is allowed, the `--help` line, `docs/20 usage.md`,
   DESIGN.md and DECISIONS.md.
2. **Resuming an existing output directory re-runs all 31 jobs.** Snakemake
   reports "Code has changed since last execution" for every rule, because
   the shell template changed. It is correct, it is not avoidable without
   dropping the `code` rerun trigger (which also catches real command
   changes), and it costs a full GTDB-Tk and Bakta pass on upgrade. It belongs
   in the release notes. snp-dists and FastTree do re-run, so a resumed run
   never keeps an unfiltered matrix.
3. **`--report-only` over a run made before this change mislabels the SNP
   matrix and tree.** The new guidance says "filtered", but the data are not.
   The checked-in `docs/assets/example-report.html` is self-consistent (old
   guidance baked in), but it and `docs/07 an example report.md` (65 SNPs,
   1,183,231 columns) are unfiltered numbers, and a re-render would mislabel
   them. Phrase the guidance as "from 3.5.0", and re-run the showcase on
   GenomeDK alongside the SLURM check.
4. **`comparem2 --setup` from a directory with a space is refused** over the
   default `--output`, which `--setup` never uses (it builds in a temp
   directory). Skip the output check for `--setup`.
5. **Two of the refusal test's three cases assert only `SystemExit`.** The
   `--setup` case passes for the wrong reason when conda is missing, and with
   conda present and the check broken it would deploy environments. Assert
   the message in all three.
6. **Text:** STATUS.md's "uncommitted at the time of writing". The
   `run_hook` docstring and DESIGN.md's sentence about surviving a spaced
   output directory describe a case the CLI now refuses. The snp-dists blurb
   still says the alignment is "the genes shared by nearly all genomes".

Checked and fine: the scip pin; `{log:q}` with the `mkdir` gone;
`_profile_argv`; the DAG test; and the guidance percentages, which reproduce
exactly when both alignments are taken from one Panaroo run (8,193 → 2,645,
6,190 → 1,893, 8,907 → 1,699). Still unverified: CI with Snakemake
installed, where the two `runner.run()` tests run for the first time.

## Draft tag message for v3.5.0, 2026-10-09

The release notes are the annotated tag message, as for 3.4.0. Edit the
numbers if the showcase re-run changes any of them.

```
Version 3.5.0: panaroo survives selenoproteins, CarveMe has a solver again

Panaroo no longer aborts on Bakta's selenoprotein genes (issue #156). Bakta
writes a selenocysteine CDS read through its TGA codon, panaroo took that
for a premature stop, and the whole pangenome branch (panaroo, snp-dists,
FastTree) failed for any set containing, for instance, E. coli. Panaroo now
runs with --remove-invalid-genes, which only matters for sets that used to
fail. And --set now reaches every tool; for seqkit, checkm2, gtdbtk, mlst,
panaroo, snp-dists and FastTree it was accepted and silently ignored.

Upgrade from 3.4.0. A fresh 3.4.0 install resolves scip 10.1, which the
pyscipopt it ships with cannot load, so every carveme job fails with "No
solver available." (reconfirmed on a fresh --setup, 2026-10-08). 3.5.0 pins
scip >=10.0.3,<10.1.

Results change. snp-dists and FastTree now read Panaroo's filtered core
alignment, the one Panaroo recommends for phylogenies. On the eight-genome
Streptococcus showcase that is 11-21% fewer SNPs per pair from 8.9% less
alignment, and one of five tree splits moved; on the four E. faecium test
genomes, 68-81% fewer, and the closest pair changed. Counts and trees from
3.4.0 or earlier are not comparable with these.

--output and --databases may now hold only letters, digits and _ . / -.
Anything else is refused before anything is written, because GTDB-Tk,
CheckM2 and Panaroo pass paths to a shell of their own and Panaroo refuses a
comma. Until now, a space in --output made every rule fail.

Re-running panaroo into an existing output directory works again. Panaroo
1.8.0 leaves a resume manifest behind and refused to start while it was
there, so in 3.4.0 adding a genome to a run, or changing a panaroo setting,
failed panaroo, snp-dists and FastTree.

Resuming an output directory made by 3.4.0 re-runs every job, including
GTDB-Tk and Bakta, because the rule code changed.

Also: log paths are quoted; --profile and the API pass the same deployment
flags; CI now hands a generated Snakefile to Snakemake; and a weekly
workflow builds and exercises every tool environment.
```
