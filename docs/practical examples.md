# Practical examples (tutorials)

!!! note
    This page is under construction. Will be finished before the end of august 2025.

CompareM2 contains different tools to cover varying aspects of microbial genome analysis. Here we showcase some different ways to use CompareM2 to achieve insights into the functional potentials of microbial genomes, and how they relate to one another.

CompareM2 has two main typical uses cases. 1) To characterize assembled isolates, for example, from a clinical setting. Here, anti-microbial resistance profiling and core genome phylogeny comes to mind. 2) To characterize complete binned microbial metagenomes. Here pan-genomes and mash-tree phylogenies come to mind.
These two use cases can of course be varied, but in these tutorials, we'll focus on these for now.

## Tutorial 1 - Publicly available *Enterococcus faecium*

Enterococcus faecium have spread on a hospital and the patients have been .. the 

The quickest way to get started with CompareM2 is to use publicly available genomes. Here, we'll investigate clusters of the *Enterococcus faecium* species, which are human commensals and pathogens, sometimes carrying vancomycin resistance.

First, we'll gain a quick phylogenetic overview over some publicly available genomes. A sufficient tree can be computed with mashtree.

After following the [installation instructions](https://comparem2.readthedocs.io/en/latest/10%20installation/), the following can be run on the command line

!!! hint
    You can generate a list of Refseq IDs of your own favorite microbial species by going to https://www.ncbi.nlm.nih.gov/datasets/genome/ and writing a species name in the search box. 

```bash
conda activate comparem2_launcher

mkdir E._faecium_public && cd E._faecium_public

comparem2 --config add_ncbi=GCF_900186865.1,GCF_028891345.1,GCF_017724035.1,GCF_014779555.2,GCF_040529065.1,GCF_900183975.1,GCF_013177435.1,GCF_009769165.1,GCF_001908725.1,GCF_040529055.1,GCF_014236795.1,GCF_028752115.1,GCF_016613535.2,GCF_041240605.1,GCF_041242975.1,GCF_041245805.1,GCF_033170385.1,GCF_031834635.1,GCF_040529045.1,GCF_030549615.1,GCF_027924565.1,GCF_016613475.2,GCF_006547045.1,GCF_002327105.1,GCF_018406605.1,GCF_018406645.1,GCF_014489475.1,GCF_004803715.2,GCF_006337105.1,GCF_007795095.1,GCF_015209725.1,GCF_015209745.1,GCF_907164845.1,GCF_024498955.1,GCF_003293695.1,GCF_000961215.1,GCF_008370835.2,GCF_013388375.1,GCF_905367715.1,GCF_017301775.1,GCF_002285515.1,GCF_008639345.1,GCF_002355975.1,GCF_014109845.1,GCF_034508895.1,GCF_041519315.1,GCF_000277895.2,GCF_009931595.1,GCF_004551575.1,GCF_014397415.1,GCF_009755585.1,GCF_041956365.1,GCF_043632995.1,GCF_000233915.3,GCF_001442805.1,GCF_025266575.1,GCF_009883735.1,GCF_030322905.1,GCF_021906995.1,GCF_003233695.1,GCF_001442745.1,GCF_040616535.1,GCF_030168915.1,GCF_014395425.1,GCF_903989455.1,GCF_023337365.3,GCF_046529185.1,GCF_031199375.1,GCF_040687735.1,GCF_032698475.1,GCF_032697525.1,GCF_040256885.1,GCF_022637495.1,GCF_000661955.1,GCF_007556775.1,GCF_031834595.1,GCF_031834615.1,GCF_021233355.1,GCF_019551435.1,GCF_041531565.1,GCF_007993735.1,GCF_000454545.1,GCF_900101175.1,GCF_013112235.1,GCF_013185915.1,GCF_010093235.1,GCF_004358165.1,GCF_003814555.1,GCF_014202435.1,GCF_015453285.1,GCF_900215535.1,GCF_004346265.1,GCF_001431525.1,GCF_010093305.1,GCF_016820515.1,GCF_009192945.1,GCF_000420545.1,GCF_010211745.1,GCF_000768345.1,GCF_000426365.1 --until mashtree
```

## Tutorial 2 - Metagenomic characterization

Assemblies of microbial genomes come in many different forms. A common way to investigate the diversity of ecological niches is to assemble and bin the microbes in a biological sample. By assembling a metagenome and binning it into discrete species representatives, it is possible to paint a picture of the functional potential they ?direct on their environment.



Here, we'll use genomes from the ?? study, which can be downloaded here.



## Tutorial 3 - Genus *Streptococcus* 

*Streptococcus* presents a diverse genus of Lactobacillales bacteria tolerating both aerobic and anaerobic growth. In addition to streptococcal pharyngitis (strep throat), certain Streptococcus species are responsible for many cases of pink eye, meningitis, bacterial pneumonia, endocarditis, erysipelas, and necrotizing fasciitis (the 'flesh-eating' bacterial infections). However, many streptococcal species are not pathogenic, and form part of the commensal human microbiota of the mouth, skin, intestine, and upper respiratory tract. Streptococci are also a necessary ingredient in producing Emmentaler ("Swiss") cheese.

Because of their relation to humans in both health and disease, their taxonomy has been subject to much scrutiny. Here, we will download and analyze all type species of genus *Streptococcus* with CompareM2, to analyze the species names and their taxonomic relations. Since CompareM2 has a built in feature to download all NCBI (GenBank and RefSeq) accessions, we just need a list of accessions to input on the command line. 

```bash


comparem2 --config add_ncbi=GCF... --until what

```



*Please reach out if you have comments or ideas about these tutorials*


{!resources/footer.md!}
