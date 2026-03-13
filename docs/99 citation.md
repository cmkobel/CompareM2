

# Citing and alternatives

## Citing CompareM2

If you use CompareM2, please cite:

  - Kobel C.M., Aho V.T.E., Øyås O., Nørskov-Lauritsen N., Woodcroft B.J., Pope P.B. CompareM2 is a genomes-to-report pipeline for comparing microbial genomes. *Bioinformatics* 41(9), btaf517 (2025). [https://doi.org/10.1093/bioinformatics/btaf517](https://doi.org/10.1093/bioinformatics/btaf517)

CompareM2 is Open Access under the [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license.

## References for included tools

CompareM2 integrates many software packages. Please cite the relevant tools when publishing results. Contact carl.mathias.kobel near nmbu.no if something is missing.

  - Bakta: Schwengers O. et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics* 7(11). [doi:10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685)
  - TreeCluster: Balaban M. et al. (2019). TreeCluster: Clustering biological sequences using phylogenetic trees. *PLOS ONE* 14, e0221068.
  - antiSMASH: Blin K. et al. (2023). antiSMASH 7.0. *Nucleic Acids Res.* 51, W46–W50.
  - eggNOG-mapper: Cantalapiedra C.P. et al. (2021). eggNOG-mapper v2. *Mol. Biol. Evol.* 38, 5825–5829.
  - GTDB-Tk: Chaumeil P.-A. et al. (2022). GTDB-Tk v2. *Bioinformatics* 38, 5315–5316.
  - CheckM2: Chklovski A. et al. (2023). CheckM2. *Nat. Methods* 20, 1203–1212.
  - AMRFinderPlus: Feldgarden M. et al. (2019). *Antimicrob. Agents Chemother.* 63, 10.1128/aac.00483-19.
  - Mashtree: Katz L.S. et al. (2019). Mashtree. *J. Open Source Softw.* 4, 1762.
  - IQ-TREE: Minh B.Q. et al. (2020). IQ-TREE 2. *Mol. Biol. Evol.* 37, 1530–1534.
  - FastTree: Price M.N. et al. (2010). FastTree 2. *PLOS ONE* 5, e9490.
  - snp-dists: Seemann T. et al. (2025). tseemann/snp-dists, GitHub.
  - MLST: Seemann T. (2025). tseemann/mlst, GitHub.
  - Prokka: Seemann T. (2014). Prokka. *Bioinformatics* 30, 2068–2069.
  - SeqKit: Shen W. et al. SeqKit2: A Swiss army knife for sequence and alignment processing. *iMeta*, Wiley.
  - Panaroo: Tonkin-Hill G. et al. (2020). Panaroo. *Genome Biol.* 21, 180.
  - clusterProfiler: Wu T. et al. (2021). clusterProfiler 4.0. *The Innovation* 2, 100141.
  - dbCAN: Yin Y. et al. (2012). dbCAN. *Nucleic Acids Res.* 40, W445–W451.
  - InterProScan: Zdobnov E.M. & Apweiler R. (2001). InterProScan. *Bioinformatics* 17, 847–848.
  - gapseq: Zimmermann J. et al. (2021). gapseq. *Genome Biol.* 22, 81.


## Alternative tools

CompareM2 works strictly downstream of assembly and binning. Unlike some competing tools, it deliberately excludes read-level processing (mapping, assembly, binning), since these steps are highly dependent on sequencing technology. CompareM2 instead provides a portable platform for comparing genomes regardless of how they were generated.

Benchmarking in the CompareM2 paper showed it to be significantly faster than Tormes (which schedules jobs sequentially) and Bactopia (which is reads-based and must generate artificial reads for assembly-only input). CompareM2's running time scales approximately linearly with input count thanks to Snakemake's parallel job scheduling.

Pipelines with partial feature overlap (alphabetical):

  - [Anvi'o](https://anvio.org/)
  - [ASA3P](https://github.com/oschwengers/asap)
  - [Aviary](https://github.com/rhysnewell/aviary)
  - [Bactopia](https://github.com/bactopia/bactopia)
  - [DRAM](https://github.com/WrightonLabCSU/DRAM)
  - [Nullarbor](https://github.com/tseemann/nullarbor)
  - [Tormes](https://github.com/nmquijada/tormes)
  - [VEBA](https://github.com/jolespin/veba)


{!resources/footer.md!}
