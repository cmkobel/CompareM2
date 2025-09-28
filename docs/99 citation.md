

# Citing and alternatives

## Citing CompareM2

If you use CompareM2, you can support further funding by bumping the citation count on this one:

  - Kobel et al. CompareM2 is a genomes-to-report pipeline for comparing microbial genomes. Bioinformatics 41, btaf517 (2025) [https://doi.org/10.1093/bioinformatics/btaf517](https://doi.org/10.1093/bioinformatics/btaf517)

## References for the included tools

CompareM2 would not have existed, if it hadn't been for the integrated software packages and their databases. Please reach out to carl.mathias.kobel near nmbu.no if you think something is missing.

  - Bacterial Antimicrobial Resistance Reference Gene ... (ID 313047) - BioProject - NCBI. https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047.
  - Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification | Microbiology Society. https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000685.
  - Balaban, M., Moshiri, N., Mai, U., Jia, X. & Mirarab, S. TreeCluster: Clustering biological sequences using phylogenetic trees. PLOS ONE 14, e0221068 (2019).
  - Baumer, B. & Udwin, D. R Markdown. WIREs Comput. Stat. 7, 167–177 (2015).
  - Blin, K. et al. antiSMASH 7.0: new and improved predictions for detection, regulation, chemical structures and visualisation. Nucleic Acids Res. 51, W46–W50 (2023).
  - Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P. & Huerta-Cepas, J. eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Mol. Biol. Evol. 38, 5825–5829 (2021).
  - Carattoli, A. & Hasman, H. PlasmidFinder and In Silico pMLST: Identification and Typing of Plasmid Replicons in Whole-Genome Sequencing (WGS). in Horizontal Gene Transfer: Methods and Protocols (ed. de la Cruz, F.) 285–294 (Springer US, New York, NY, 2020). doi:10.1007/978-1-4939-9877-7_20.
  - Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P. & Parks, D. H. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. Bioinformatics 38, 5315–5316 (2022).
  - Chklovski, A., Parks, D. H., Woodcroft, B. J. & Tyson, G. W. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat. Methods 20, 1203–1212 (2023).
  - Dykstra, D. Apptainer Without Setuid. EPJ Web Conf. 295, 07005 (2024).
  - Jette, M. A. & Wickberg, T. Architecture of the Slurm Workload Manager. in Job Scheduling Strategies for Parallel Processing (eds. Klusáček, D., Corbalán, J. & Rodrigo, G. P.) 3–23 (Springer Nature Switzerland, Cham, 2023). doi:10.1007/978-3-031-43943-8_1.
  - Jolley, K. A. & Maiden, M. C. BIGSdb: Scalable analysis of bacterial genome variation at the population level. BMC Bioinformatics 11, 595 (2010).
  - Katz, L. S. et al. Mashtree: a rapid comparison of whole genome sequence files. J. Open Source Softw. 4, 1762 (2019).
  - Liu, B., Zheng, D., Zhou, S., Chen, L. & Yang, J. VFDB 2022: a general classification scheme for bacterial virulence factors. Nucleic Acids Res. 50, D912–D917 (2022).
  - Miell, I. & Sayers, A. Docker in Practice, Second Edition. (Simon and Schuster, 2019).
  - Minh, B. Q. et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. Evol. 37, 1530–1534 (2020).
  - Mölder, F. et al. Sustainable data analysis with Snakemake. Preprint at https://doi.org/10.12688/f1000research.29032.2 (2021).
  - Page, A. J. et al. Roary: rapid large-scale prokaryote pan genome analysis. Bioinformatics 31, 3691–3693 (2015).
  - Price, M. N., Dehal, P. S. & Arkin, A. P. FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. PLOS ONE 5, e9490 (2010).
  - sanger-pathogens/assembly-stats. Pathogen Informatics, Wellcome Sanger Institute (2024).
  - Seemann, T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 30, 2068–2069 (2014).
  - Seemann, T. tseemann/abricate. (2024).
  - Seemann, T. tseemann/any2fasta. (2024).
  - Seemann, T. tseemann/mlst. (2024).
  - Seemann, T. tseemann/snp-dists. (2024).
  - SeqKit2: A Swiss army knife for sequence and alignment processing - Shen - iMeta - Wiley Online Library. https://onlinelibrary.wiley.com/doi/10.1002/imt2.191.
  - Smith, K. W. et al. A standardized nomenclature for resistance-modifying agents in the Comprehensive Antibiotic Resistance Database. Microbiol. Spectr. 11, e0274423 (2023).
  - Tonkin-Hill, G. et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol. 21, 180 (2020).
  - Wickham, H. et al. Welcome to the Tidyverse. J. Open Source Softw. 4, 1686 (2019).
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
