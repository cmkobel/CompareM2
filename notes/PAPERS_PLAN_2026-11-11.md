# The reference papers, and the one unsourced number — plan, 2026-11-11

Written against `cbae03b` (master). It started as a question: are the tools'
papers in the repo? The citations all are: every one of the 14 catalogue tools
has one in `guidance.py` and an entry in `papers/tools.bib` (28 entries). The
PDFs are not, and that is intended: `papers/*.pdf` has been gitignored since
2026-09-02. This checkout has only the CompareM2 paper, which is tracked
despite the ignore rule.

Two conclusions were reached in the same conversation. They are recorded
under *Decided* so the next session does not reopen them.

## 1. iML1515's "~0.88 h⁻¹" has no source — measure it, then reword

`guidance.py:1044`, in the biosynthesis caveat about the growth column:

> iML1515 comes out at 0.53 h⁻¹ on M9 glucose that way, against the ~0.88 h⁻¹
> it is usually simulated at with oxygen unconstrained.

The 0.88 appears nowhere in `SUMMARIES.md` or `DECISIONS.md`, and nothing
attributes it to Monk 2017 (whose PDF was never fetched). It is a number in the
report with no audit trail, which the sourcing rule in CLAUDE.md does not allow.
The 0.53 *is* ours: STATUS.md, 2026-09-10, `iML1515` M9 = 0.5315. But the
sentence does not say so.

**Fix: replace the claim about common practice with a second measurement on the
same model.** On thylakoid, in the carveme environment, load `iML1515.xml`
through `biosynthesis.py` exactly as `upstream/panel_calibration.py --model`
does. Solve on M9 twice: once as the module does it, with every compound
including O₂ capped at 10 mmol/gDW/h, which should reproduce 0.5315, and once
with `EX_o2_e` unconstrained. Then rewrite the sentence to cite both
numbers as this pipeline's measurement, for example *"…comes out at 0.53 h⁻¹ on
M9 glucose this way and at X h⁻¹ with oxygen unconstrained (iML1515, measured
on thylakoid, <date>)"*.

**Why not verify against Monk 2017 instead:** "usually simulated at" reads as
common practice, not as a figure the paper states. It may not be in the paper,
and then we would have to measure it anyway. A measurement on our own code
passes the rule outright and also labels the 0.53.

**Why not drop the number:** it is the one concrete figure that makes the
caveat land.

Ships with: the `guidance.py` sentence, a STATUS.md row next to the 2026-09-10
table, and a dated DECISIONS.md entry saying the 0.88 was unsourced and why
it is now a measurement. Grep `0.88` and `0\.53` afterwards.

**Unverified until it runs.** If the unconstrained value is not about 0.88,
the old sentence was wrong as well as unsourced, and the new one simply states
whatever comes out.

## 2. `papers/README.md` has four stale statements

- "The 13 catalogue tools" — there are 14, and biosynthesis is missing from
  the table. Its row cites CarveMe, with Orth 2010 (FBA) and Monk 2017
  (iML1515) as methods. Neither has a PDF or a licence check.
- "`tools.bib` holds 26 BibTeX entries" — 28 now (Orth and Monk were added).
- "a copy of the PDF already in the repo root" — that PDF moved to `papers/`
  in `7a6f5b9`.
- The sylph PDF "is still in this directory" — true of one machine, not of
  this checkout. Say that, or drop the paragraph.

## 3. `papers/fetch.sh` — the PDFs on demand

Turn *How these were fetched* into a script that reads the DOIs and fills the
gitignored directory. Publisher-direct where the publisher serves PDFs, and
Europe PMC's `?pdf=render` with a browser user-agent otherwise; the README
already lists which is which. Check each PDF with `pdftotext -f 1 -l 1`
against its title, as was done originally. Its job is that an agent that
needs to check a number against a paper runs one command instead of
re-learning Europe PMC's quirks. Not run yet, so not verified.

## 4. Untrack the CompareM2 PDF — Carl's call

`papers/Kobel et al. - 2024 - ….pdf` (1.2 MB) is the one PDF tracked against
the ignore rule. Untracking it does not shrink history, but it removes the
file from every release tarball the bioconda recipe downloads
(`recipe/meta.yaml:9` builds from the GitHub tag archive). The DOI is in
`citation.cff`. Optional.

## Decided, not to be redone

- **The PDFs stay out of git.** They are not useful agent context:
  `SUMMARIES.md` already holds every number with its verbatim quote, and a
  PDF costs tens of thousands of tokens through `pdftotext`. Tracking them
  would put about 50 MB into every source tarball and into history for good.
  Not all of them may be redistributed either: MAFFT is CC BY-NC, the
  published CheckM2 paper has no CC licence, and Orth and Monk were never
  checked. The only task that needs a paper is checking a number against it,
  and that needs one paper, fetched by DOI (item 3).
- **Biosynthesis gets no `SUMMARIES.md` section.** That file is the audit trail
  for claims taken from other people's papers; biosynthesis is ours, and its
  record is `biosynthesis.py`, `upstream/panel_calibration.py` and
  DECISIONS.md. The places where its guidance does borrow were checked: the
  two CarveMe quotes are in `SUMMARIES.md`'s carveme section, and the
  10 mmol/gDW/h protocol is at `DECISIONS.md:580`. The only gap is item 1.
  Orth and Monk are cited for the method and the reference model, not for
  numbers, so they need no summaries.

## Order

1 needs thylakoid; 2–4 do not. 2 and 3 go in one commit (the README describes
the script). 1 is its own commit with its DECISIONS and STATUS entries.
