# An example report

**[Open the example report](assets/example-report.html)**: a real run, all
fourteen tools, eight genomes. Nothing in it is loaded from the network. It is
one self-contained HTML file, which is what every run produces.

If you would rather not download a database first, `comparem2 --demo` runs four
of the tools over six bundled *Enterococcus faecium* plasmids. This page is the
other end of the scale.

## The genomes

Eight complete genomes from the *Streptococcus mitis* group: six pneumococci
spanning serotypes and both pandemic multidrug-resistant lineages, plus the two
close relatives that make the species boundary interesting.

| Sample | Accession | Why it is here |
| --- | --- | --- |
| `Spn_TIGR4` | GCF_000006885.1 | Serotype 4, ST205, the first pneumococcal genome |
| `Spn_D39` | GCF_000014365.2 | Serotype 2, Avery's transformation strain |
| `Spn_R6` | GCF_000007045.1 | D39's unencapsulated laboratory derivative |
| `Spn_ATCC700669` | GCF_000026665.1 | PMEN1 / Spain²³F ST81; pandemic, penicillin-resistant |
| `Spn_Hungary19A` | GCF_000019265.1 | Serotype 19A ST268, a post-PCV7 replacement lineage |
| `Spn_P1031` | GCF_000019005.1 | Serotype 1 ST303, the meningitis-belt serotype |
| `Spseudo_IS7493` | GCF_000221985.1 | *S. pseudopneumoniae*, the "is it pneumococcus?" case |
| `Smitis_B6` | GCF_000027165.1 | *S. mitis*, penicillin-resistant commensal |

All eight are 2.04–2.25 Mb with 1,992–2,242 coding sequences, and CheckM2 puts
every one at 100.0% complete (*S. mitis* B6 at 99.97%) with 0.09–0.38%
contamination. Nothing below is an artefact of a poor assembly.

## Four things worth checking in it

### D39 and R6 are the same strain, and most of the report says so

R6 is a laboratory derivative of D39 that lost its capsule locus, and the
report independently recovers that: 100.00% ANI, 65 SNPs across a
1,183,231-column core alignment, the same ST595, branch lengths of 1.7e-5 and
4.2e-5 from their common ancestor, and one AMR gene each. A near-identical pair
is the cheapest sanity check you can put in a run, and it is worth doing with
your own data. See [what the report gets wrong
here](#what-the-report-gets-wrong-here).

### The species boundary is visible, and it is not at 95%

skani puts the six pneumococci at 98.31–98.99% ANI to each other,
*S. pseudopneumoniae* at 93.98–94.51%, under the conventional 95% line, and
*S. mitis* at 91.03–92.09%. GTDB-Tk resolves all three as separate species, and
colours each ANI against *that reference's own* radius instead of a global
threshold.

### A commensal carries the most resistance genes

AMRFinderPlus finds 9 resistance classes across the set. The pneumococci run 1,
1, 1, 2, 3 and 5 genes, with the MDR clones at the top. *S. mitis* B6, the
commensal, has 8, more than any pathogen here. The mitis group is a resistance
reservoir for pneumococcus, and this is the shape of it.

### MLST declines to answer for the wrong species

*S. pseudopneumoniae* is typed against the pneumococcal scheme, returns no
sequence type, and reports inexact alleles such as `ddl(656?)`. That is the
tool being right about its own limits: an MLST scheme is species-specific, and
a genome outside it gets no ST instead of a wrong one.

## What the report gets wrong here

The metabolism sections disagree with themselves on this set, and the D39/R6
pair is what shows it. The biosynthesis panel calls 18 of 32 compounds *de
novo* for D39 and 0 for R6, flipping 20 of the 32 in one direction, for two
genomes 65 SNPs apart. Their CarveMe models share 1,048 reactions, with 530
unique to R6 and 85 unique to D39.

The cause is in the model reconstruction and not in the panel, and it is
specific: four of the eight models cannot take up ammonium. M9's only nitrogen
source is ammonium, so a model missing a link in that uptake chain can build
nothing at all from a minimal medium. One hole, thirty-two zeros. The four
affected are exactly the four reporting 0 de novo:

| model | `EX_nh4_e` | `NH4tex` | `NH4tpp` | de novo |
| --- | :-: | :-: | :-: | --: |
| D39, P1031, *S. pseudopneumoniae*, *S. mitis* | ✓ | ✓ | ✓ | 17–18 |
| R6, TIGR4, Hungary19A | — | — | ✓ | **0** |
| ATCC700669 | ✓ | ✓ | — | **0** |

All eight carry glutamate dehydrogenase, aspartate transaminase and alanine
transaminase with identical bounds, so the enzymes were never the difference.
Only whether nitrogen could reach them. In the working models the route is
gene-associated, with glutamate and 2-oxoglutarate cycling catalytically
instead of being synthesised: neither is reachable from M9 in *any* of the
eight, which is the expected result for a lactic acid bacterium with no
oxidative TCA cycle.

So the zeros are a defect and the 17–18 are the honest answers here. The reason
the defect is not merely bad luck is that R6's model is the one that converged
to a certified optimum, in 15.4 s, with the most reactions of any model in the
set at 1,579. CarveMe's objective does not require the network to be able to
eat, and nothing downstream checks.

Giving the solver more time does not fix it, which is measured and not assumed:
re-solved with eighteen times the budget, D39 ran a full 3 h and came back with
one more reaction, 1,135 against 1,134. Convergence does not fix it either.
`Spn_P1031` reaches a certified optimum too, and answers 18.

Finding this changed the biosynthesis section, so the report now says it
itself. The section checks whether each model can reach the minimal medium's
carbon and nitrogen sources, and names those that cannot directly beneath the
*de novo* counts, so the qualification sits next to the number a reader would
quote. In this report that line reads *"4 of 8 models cannot reach a source
element from the minimal medium."*

The check is a flux probe and not a look at the model's exchange reactions,
which is why ATCC700669 is caught: it *has* an ammonium exchange and reads a
reassuring 18 of 20 medium compounds present. It simply cannot move ammonium
into the cytoplasm. R6 and D39 both read 17 of 20 while differing in which
three they lack, which is what let this hide behind a count.

So read the carveme reaction counts (1,037–1,579 here) and the biosynthesis
verdicts as properties of *the draft model*, not as measured phenotypes. The
other twelve sections are unaffected, since they read the assemblies and the
annotation and not the model.

We would rather ship a showcase that says this than one that hides it.

## How it was produced

On GenomeDK, submitting to SLURM through a Snakemake profile:

```bash
comparem2 showcase_genomes/*.fna -o results_showcase --profile profiles/showcase
```

That is the command in the report's own provenance header, with the cluster's
absolute path prefix removed so it reads as something you could type.

CompareM2 3.3.0, 53 jobs, six conda environments deployed by Snakemake.
Thirteen of the fourteen tools finished in 25 minutes of wall time
(12:44:14 → 13:09) across the queue. GTDB-Tk then waited on its database: a
one-time 60.8 GB download that took 2 h 15 m at a measured 10.8 MB/s, after
which classification itself took under a minute for all eight genomes. That
download is why [`comparem2 --setup`](10 installation.md) and the printed cost
estimate exist, and why you should consider
[running a subset](20 usage.md#running-a-subset) if you do not need taxonomy.

Every version and citation for the run is in the report's own *Methods and
citations* section, generated from what actually ran.
