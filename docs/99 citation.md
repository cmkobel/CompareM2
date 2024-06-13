

# Citing and alternatives

## Citing CompareM2

If you use CompareM2, you can support further funding by bumping the citation count on this one:

  - Kobel CM. *CompareM2* **GitHub** https://github.com/cmkobel/comparem2

## References for the included tools

CompareM2 would not have existed, if it hadn't been for the integrated software packages and their databases. Please reach out to carl.mathias.kobel near nmbu.no if you think something is missing.

  - Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, Forster J, Lee S, Twardziok SO, Kanitz A, Wilm A, Holtgrewe M, Rahmann S, Nahnsen S, Köster J. Sustainable data analysis with Snakemake. F1000Res. 2021 Apr 19;10:33. doi: 10.12688/f1000research.29032.2. PMCID: PMC8114187.
  - Seemann T, Goncalves da Silva A, Bulach DM, Schultz MB, Kwong JC, Howden BP. Nullarbor Github https://github.com/tseemann/nullarbor
  - W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962. 
  - Richard Challis. (2017). rjchallis/assembly-stats 17.02 (17.02). Zenodo. https://doi.org/10.5281/zenodo.322347
  - Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323
  - Chklovski, A., Parks, D.H., Woodcroft, B.J. et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat Methods 20, 1203–1212 (2023). https://doi.org/10.1038/s41592-023-01940-w
  - Seemann T., Prokka: rapid prokaryotic genome annotation, Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063
  - Yin Y*, Mao X*, Yang JC, Chen X, Mao F and Xu Y, dbCAN: a web resource for automated carbohydrate-active enzyme annotation, Nucleic Acids Res. 2012 Jul;40(Web Server issue):W445-51 
  - P. Jones et al., Bioinformatics (2014), PMID: 24451626
  - Seemann T, Abricate, Github https://github.com/tseemann/abricate
  - Katz et al., (2019). Mashtree: a rapid comparison of whole genome sequence files. Journal of Open Source Software, 4(44), 1762, https://doi.org/10.21105/joss.01762
  - Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018 Sep 24;3:124. doi: 10.12688/wellcomeopenres.14826.1. PMID: 30345391; PMCID: PMC6192448.
  - Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics. 2019 Nov 15;36(6):1925–7. doi: 10.1093/bioinformatics/btz848. Epub ahead of print. PMID: 31730192; PMCID: PMC7703759.
  - Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
  - B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., 37:1530-1534. https://doi.org/10.1093/molbev/msaa015
  - mOTUpan: a robust Bayesian approach to leverage metagenome assembled genomes for core-genome estimation Moritz Buck, Maliheh Mehrshad, and Stefan Bertilsson bioRxiv 2021.06.25.449606; doi: https://doi.org/10.1101/2021.06.25.449606


## Alternative tools

What is unique about CompareM2 is that it works strictly downstream of assembling and binning. Many other tools also include all the steps necessary to turn raw reads into genome representatives, and then does varying degrees of biological characterization of these freshly created assemblies/bins/genomes. It is a conscious decision to exclude the raw read-dependent tools out of the equation for CompareM2. This is because read-mapping, assembling or even binning is highly dependent on the sequencing technology used and requires a highly specialized pipeline for each technological use case. CompareM2 takes a different approach which is to offer a portable and flexible platform where you can easily compare your genomes, no matter where they came from, regardless of the sequencing technology used to create them in the first place. Genome quality is only increasing and in the future we will not have to be worried when comparing pyrosequencing and single-molecule sequencing or hybrid approach based genomes in a single batch of CompareM2. 

Below we are listing some competing pipelines that partly overlap with the use cases of CompareM2.

  - [Aviary](https://github.com/rhysnewell/aviary)
  - [Bactopia](https://github.com/bactopia/bactopia)
  - [DRAM](https://github.com/WrightonLabCSU/DRAM)
  - [Nullarbor](https://github.com/tseemann/nullarbor) 
  - [Tormes](https://github.com/nmquijada/tormes)
  - [VEBA](https://github.com/jolespin/veba)



{!resources/footer.md!}
