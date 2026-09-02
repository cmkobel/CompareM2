# Reference papers for the v3 tool set

23 papers for the current tool set. Downloaded 2026-09-02. **The PDFs are not
tracked; this file, `SUMMARIES.md` and `tools.bib` are** (decided 2026-09-02).
The DOI tables below plus *How these were fetched* are what make the PDFs
recoverable, which is the reason they can stay out of git history.

A 24th PDF, `sylph - Shaw 2024 - sylph.pdf` (6.5 MB), is still in this directory
but is no longer part of the set: sylph left the catalogue 2026-09-02 because it
profiles metagenomic reads rather than assemblies (DECISIONS.md). Its summary
section and BibTeX entry are gone; the file itself is untracked and can be
deleted whenever.

`_comparem2 - Kobel 2025 - CompareM2.pdf` is a copy of the PDF already in the
repo root, which despite its `2024` filename is the version of record
(*Bioinformatics* 2025, 41(9), btaf517, 4 pp). Identical bytes, md5 `4618db50e19a`.

## Licences, if these are ever redistributed

Read from each PDF on disk — the text stamp or the embedded `creativecommons.org`
link, not inferred from the journal. 22 of the 23 could legally ship with
CompareM2; one could not. (Sylph's paper, checked at the same time, was also
CC BY 4.0.)

| Licence | n | Papers |
| ------- | - | ------ |
| CC BY 4.0 | 16 | comparem2, diamond, gtdb, mash, pyrodigal, snakemake, amrfinder-2021, bakta, carveme, gtdbtk, mashtree, mlst, panaroo, seqkit-2016, skani, treecluster |
| CC BY 3.0 / 2.0 | 2 | cd-hit (`by/3.0`), prodigal (`by/2.0`) |
| CC BY, version not printed | 2 | fasttree (PLOS ONE 2010), seqkit2 (iMeta/Wiley) |
| US Government work | 1 | amrfinder-2019 — "not subject to copyright protection in the United States. Foreign copyrights may apply." |
| CC BY-ND 4.0 | 1 | checkm2 preprint — verbatim redistribution *is* permitted; ND restricts derivatives, not distribution |
| CC BY-NC 3.0 | 1 | **mafft** (Katoh 2013) — the one blocker |

MAFFT is non-commercial only, and CompareM2 is GPL-3, which anyone may use
commercially. So MAFFT stays a citation and a link, never a bundled file.

Bundling any of the rest brings an attribution obligation with it: author,
title, licence, link, and a statement that no changes were made. The published
*Nature Methods* CheckM2 paper is also not redistributable (Springer Nature
text-and-data-mining licence, no CC) — only the preprint here is.

`tools.bib` holds 26 BibTeX entries, fetched from Crossref (or DataCite, for the
Zenodo record) rather than hand-typed. Every PDF was checked with
`pdftotext -f 1 -l 1` against its expected title, so none is a stray error page.

`SUMMARIES.md` is the long-form reading of all 23 papers (2026-09-02): per tool,
what it does, how to read each column, caveats, a table of every quantitative
claim with the verbatim quote it came from, and a *Not established* list of what
the papers do not answer. The condensed version is `src/comparem2/guidance.py`,
which is what the report renders — this file is the audit trail behind it.

## The 13 catalogue tools

Checked against the `CATALOGUE` registry in `src/comparem2/catalogue.py`, not
against DESIGN.md.

| Tool | Citation | DOI | PDF |
| ---- | -------- | --- | --- |
| seqkit | Shen et al. 2016, *PLOS ONE* 11:e0163962 | 10.1371/journal.pone.0163962 | yes |
| seqkit | Shen et al. 2024, *iMeta* 3:e191 (SeqKit2) | 10.1002/imt2.191 | yes |
| checkm2 | Chklovski et al. 2023, *Nat. Methods* 20:1203–1212 | 10.1038/s41592-023-01940-w | preprint (see below) |
| gtdbtk | Chaumeil et al. 2022, *Bioinformatics* 38:5315–5316 | 10.1093/bioinformatics/btac672 | yes |
| bakta | Schwengers et al. 2021, *Microb. Genom.* 7:000685 | 10.1099/mgen.0.000685 | yes |
| amrfinder | Feldgarden et al. 2021, *Sci. Rep.* 11:12728 | 10.1038/s41598-021-91456-0 | yes |
| amrfinder | Feldgarden et al. 2019, *AAC* 63:e00483-19 (validation) | 10.1128/AAC.00483-19 | yes |
| mlst | Jolley et al. 2018, *Wellcome Open Res.* 3:124 (PubMLST) | 10.12688/wellcomeopenres.14826.1 | yes |
| mashtree | Katz et al. 2019, *JOSS* 4:1762 | 10.21105/joss.01762 | yes |
| treecluster | Balaban et al. 2019, *PLOS ONE* 14:e0221068 | 10.1371/journal.pone.0221068 | yes |
| skani | Shaw & Yu 2023, *Nat. Methods* 20:1661–1665 | 10.1038/s41592-023-02018-3 | yes |
| panaroo | Tonkin-Hill et al. 2020, *Genome Biol.* 21:180 | 10.1186/s13059-020-02090-4 | yes |
| snp-dists | Seemann T., Zenodo (source archive) | 10.5281/zenodo.1411986 | **no paper exists** |
| fasttree | Price et al. 2010, *PLOS ONE* 5:e9490 | 10.1371/journal.pone.0009490 | yes |
| carveme | Machado et al. 2018, *Nucleic Acids Res.* 46:7542–7553 | 10.1093/nar/gky537 | yes |

Two tools have no paper of their own, so the row cites what should actually be
credited:

- **mlst** — `tseemann/mlst` is unpublished; it queries PubMLST, so Jolley et
  al. 2018 is the citation. Cite the GitHub repository alongside it.
- **snp-dists** — `tseemann/snp-dists` is unpublished and has no preprint. The
  Zenodo record is a source-code archive, not an article; there is nothing to
  download. Cite the DOI.

## Methods the tools run, which need citing in their own right

Prefixed `_` so they sort together. These are not catalogue entries, but the
pipeline executes or ships them and the report describes what they do.

| | Citation | DOI | Why it matters here |
| --- | -------- | --- | ------------------- |
| GTDB | Parks et al. 2022, *NAR* 50:D785–D794 | 10.1093/nar/gkab776 | The 141 GB database itself. GTDB-Tk's own docs ask for both citations. |
| Mash | Ondov et al. 2016, *Genome Biol.* 17:132 | 10.1186/s13059-016-0997-x | Mashtree is a wrapper; Mash is the distance measure. The JOSS note is 6 pp. |
| DIAMOND | Buchfink et al. 2021, *Nat. Methods* 18:366–368 | 10.1038/s41592-021-01101-x | Used by checkm2, bakta and carveme. Its version pin is why checkm2 is `isolated=True`. |
| Prodigal | Hyatt et al. 2010, *BMC Bioinformatics* 11:119 | 10.1186/1471-2105-11-119 | Bakta's gene caller. |
| Pyrodigal | Larralde 2022, *JOSS* 7:4296 | 10.21105/joss.04296 | The binding bakta actually calls — and the 3.x rename that broke bakta 1.8.1. |
| MAFFT | Katoh & Standley 2013, *MBE* 30:772–780 | 10.1093/molbev/mst010 | Builds panaroo's core alignment, which fasttree and snp-dists both consume. |
| CD-HIT | Fu et al. 2012, *Bioinformatics* 28:3150–3152 | 10.1093/bioinformatics/bts565 | Panaroo's clustering step. |
| Snakemake | Mölder et al. 2025, *F1000Research* 10:33 (**v3**) | 10.12688/f1000research.29032.3 | The executor. Supplied by Carl; v3 supersedes the v2 I first pulled. |
| CompareM2 | Kobel et al. 2025, *Bioinformatics* 41:btaf517 | 10.1093/bioinformatics/btaf517 | The pipeline itself. CC-BY 4.0. |

## Gaps

None outstanding. Three things are deliberately not the version of record:

1. **CheckM2 is the bioRxiv preprint, by decision (Carl, 2026-09-02).**
   10.1101/2022.07.11.499243 v1, CC-BY-ND, 24 pp, posted 2022-07-11; bioRxiv
   links it to the published DOI. The *Nature Methods* version is
   subscription-only — Crossref lists a Springer Nature text-and-data-mining
   licence and no Creative Commons, and `nature.com` serves HTML, not a PDF.
   **Cite the published paper (10.1038/s41592-023-01940-w); the PDF here is the
   preprint.** `tools.bib` carries both entries.
2. **Snakemake 2012** (Köster & Rahmann, *Bioinformatics* 28:2520–2522,
   10.1093/bioinformatics/bts480) — OUP returns 403 and it is not in PMC. Not a
   real gap: the F1000Research paper supersedes it and is here. BibTeX only.
3. **snp-dists and mlst** have no publications at all. Nothing to fetch.

## How these were fetched

Publisher-direct where the publisher serves PDFs (PLOS, Nature OA, Microbiology
Society, JOSS, Wellcome Open Research, F1000Research), and
`https://europepmc.org/articles/<PMCID>?pdf=render` with a browser user-agent
for the ones behind a CDN that blocks direct requests (Oxford, Wiley, BMC, ASM).
The Europe PMC REST `fullTextPDF` endpoint 404s and `pmc.ncbi.nlm.nih.gov`
serves a bot-check page; neither is usable.
