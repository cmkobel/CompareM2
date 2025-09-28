

# Citing and alternatives

## Citing CompareM2

If you use CompareM2, you can support further funding by bumping the citation count on this one:

  - Kobel et al. CompareM2 is a genomes-to-report pipeline for comparing microbial genomes. Bioinformatics 41, btaf517 (2025) [https://doi.org/10.1093/bioinformatics/btaf517](https://doi.org/10.1093/bioinformatics/btaf517)

## References for the included tools

CompareM2 would not have existed, if it hadn't been for the integrated software packages and their databases. Please reach out to carl.mathias.kobel near nmbu.no if you think something is missing.

  - Andersen, T. O. et al. Metabolic influence of core ciliates within the rumen microbiome. ISME J. 17, 1128–1140 (2023).
  - Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification | Microbiology Society.
  - Balaban, M., Moshiri, N., Mai, U., Jia, X. & Mirarab, S. TreeCluster: Clustering biological sequences using phylogenetic trees. PLOS ONE 14, e0221068 (2019).
  - Blin, K. et al. antiSMASH 7.0: new and improved predictions for detection, regulation, chemical structures and visualisation. Nucleic Acids Res. 51, W46–W50 (2023).
  - Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P. & Huerta-Cepas, J. eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Mol. Biol. Evol. 38, 5825–5829 (2021).
  - Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P. & Parks, D. H. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. Bioinformatics 38, 5315–5316 (2022).
  - Chklovski, A., Parks, D. H., Woodcroft, B. J. & Tyson, G. W. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat. Methods 20, 1203–1212 (2023).
  - Feldgarden, M. et al. Validating the AMRFinder Tool and Resistance Gene Database by Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of Isolates. Antimicrob. Agents Chemother. 63, 10.1128/aac.00483-19 (2019).
  - Huang, W., Li, L., Myers, J. R. & Marth, G. T. ART: a next-generation sequencing read simulator. Bioinformatics 28, 593–594 (2012).
  - Katz, L. S. et al. Mashtree: a rapid comparison of whole genome sequence files. J. Open Source Softw. 4, 1762 (2019).
  - Lazear, M. R. Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale. J. Proteome Res. 22, 3652–3659 (2023).
  - Minh, B. Q. et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. Evol. 37, 1530–1534 (2020).
  - Price, M. N., Dehal, P. S. & Arkin, A. P. FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. PLOS ONE 5, e9490 (2010).
  - Schwengers O., Jelonek L., Dieckmann M. A., Beyvers S., Blom J., Goesmann A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11). https://doi.org/10.1099/mgen.0.000685
  - Seemann T, F Klötzl, AJ Page., tseemann/snp-dists Github (2025).
  - Seemann T., tseemann/mlst Github (2025).        
  - Seemann, T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 30, 2068–2069 (2014).
  - Shaffer, M. et al. DRAM for distilling microbial metabolism to automate the curation of microbiome function. Nucleic Acids Res. 48, 8883–8900 (2020).
  - Shen W, Sipos B, Zhao L. SeqKit2: A Swiss army knife for sequence and alignment processing, iMeta - Wiley Online Library.
  - Tonkin-Hill, G. et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol. 21, 180 (2020).
  - Wu, T. et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation 2, 100141 (2021).
  - Yin, Y. et al. dbCAN: a web resource for automated carbohydrate-active enzyme annotation. Nucleic Acids Res. 40, W445–W451 (2012).
  - Zdobnov, E. M. & Apweiler, R. InterProScan – an integration platform for the signature-recognition methods in InterPro. Bioinformatics 17, 847–848 (2001).
  - Zimmermann, J., Kaleta, C. & Waschina, S. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. Genome Biol. 22, 81 (2021).


## Alternative tools

What is unique about CompareM2 is that it works strictly downstream of assembling and binning. Many other tools also include all the steps necessary to turn raw reads into genome representatives, and then does varying degrees of biological characterization of these freshly created assemblies/bins/genomes. It is a conscious decision to exclude the raw read-dependent tools out of the equation for CompareM2. This is because read-mapping, assembling or even binning is highly dependent on the sequencing technology used and requires a highly specialized pipeline for each technological use case. CompareM2 takes a different approach which is to offer a portable and flexible platform where you can easily compare your genomes, no matter where they came from, regardless of the sequencing technology used to create them in the first place. Genome quality is only increasing and in the future we will not have to be worried when comparing pyrosequencing and single-molecule sequencing or hybrid approach based genomes in a single batch of CompareM2. 

Below we are listing some competing pipelines that partly overlap with the use cases of CompareM2. Sorted alphabetically.
  
  - [Anvi'o](https://anvio.org/)
  - [ASA³P](https://github.com/oschwengers/asap)
  - [Aviary](https://github.com/rhysnewell/aviary)
  - [Bactopia](https://github.com/bactopia/bactopia)
  - [DRAM](https://github.com/WrightonLabCSU/DRAM)
  - [Nullarbor](https://github.com/tseemann/nullarbor) 
  - [Tormes](https://github.com/nmquijada/tormes)
  - [VEBA](https://github.com/jolespin/veba)



{!resources/footer.md!}
