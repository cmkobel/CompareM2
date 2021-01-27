# .tabseq

You probably love fasta. But tell me - can you quickly concatenate all the genes from a fasta alignment of several samples? Can you quickly create a core genome from a list of fasta files containing the harbored genes? If not, then take a look at the _.tabseq_-format.

This is a very simple idea: Use the tab-separated-values format (.tsv or .tab) for sequence data.

Here I come up with some reasonable data-features and conversion scripts that should make it easy to get started ditching the old fasta-format.

I got tired of reformatting fasta and xmfa files again and again. Sometimes you want an alignment where the core genes are concatenated, sometimes you want each gene in order with its own header. It all depends on what software you want to run on your data. If you find that your .fasta or .xmfa data files have the wrong structure for the analysis you want to carry out, you will find yourself spending lots of time restructuring the file until it does.

The file format I propose here purposedly solves the aforementioned problems. 

This repository contains the documentation of the .tabseq-format as well as a series of conversion-scripts, that makes it easy to import and export between all imagineable structures and formats.

TODO: explain structures better.

## Formatting
As the name suggests, this file format is an heir of the tab-separated-value format. Each line contains one sequence trailed by a number of features.

Four features are generally required: `species`, `sample`, `gene` and `sequence`.

The `gene` feature specifies the name of the gene. If it is not needed for the data at hand, the value can be set to `NA` or just an empty string.

The header begins with a hash, to make it easy to filter out using command line tools:

_Note: in the following examples, the use of tabs is emphasized using tab-symbols_ `⇥`_._

**Example 1**: alignment of one gene from several samples 
```
#species⇥sample⇥part⇥sequence
Rhizobium leguminosarum⇥3789⇥rpoB⇥ATGTGCAGCCGATGATTCTACTAGTGC
Rhizobium leguminosarum⇥3790⇥rpoB⇥ATGGCAGCCGATGATTCTACTAGTGCT
Rhizobium leguminosarum⇥3792⇥rpoB⇥ATGCCAGCCGATGATTCTACTAGTGCT
```

**Example 2**: alignment of full genomes from several samples
```
#species⇥sample⇥part⇥sequence
Rhizobium leguminosarum⇥3789⇥NA⇥NNNNNNNNATGTGTGTTTATATAGATTAN
Rhizobium leguminosarum⇥3790⇥NA⇥NNNNNNNNATGTGNGTTTATATAGATTAN
Rhizobium leguminosarum⇥3792⇥NA⇥NNNNNNNNATGTCTGTTTATATAGATTAN
```

**Example 3.1**: All genes from one sample
```
#species⇥sample⇥part⇥sequence
Rhizobium leguminosarum⇥3789⇥gutA⇥ATGCGATGTGAGCACGCACAGCA
Rhizobium leguminosarum⇥3789⇥rpoB⇥ATGATATAGTGACTGACATGCAG
Rhizobium leguminosarum⇥3789⇥Glyt⇥ATGCTGATCTGCGCCACGTGAAA
```

Of course, if you feel like you need more features, just add them. These will be ignored by the bridge-scripts. 

**Example 3.2**: All genes from one sample with a custom feature added
```
#species⇥sample⇥part⇥position⇥sequence
Rhizobium leguminosarum⇥3789⇥gutA⇥10-33⇥ATGCGATGTGAGCACGCACAGCA
Rhizobium leguminosarum⇥3789⇥rpoB⇥44-67⇥ATGATATAGTGACTGACATGCAG
Rhizobium leguminosarum⇥3789⇥Glyt⇥98-75⇥ATGCTGATCTGCGCCACGTGAAA
```


If you just want to store a sequence, and nothing else, it is not a problem. Just be aware that an empty string is not equivalent to `NA`.


**Example 4**: just a sequence
```
#species⇥sample⇥part⇥sequence
⇥⇥NA⇥ATCTGCGTCACGACGTACGATAACGATCTCATTATGACATCTACTACGAT
```

And of course you can save aminoacids as well


**Example 5**: just a sequence of amino acids
```
#species⇥sample⇥part⇥sequence
⇥⇥NA⇥GSTTSAAVGSILSEEGVPINSGCO
```


## Conversion-scripts

In order to help converting between formats, I created (in progress) a series of scripts "conversion-scripts" that help with this.

At the time of writing, only one such scripts exists: `tabseq2fasta.py`
This file will convert the tabseq to fasta, by the arguments given.


# Installation of R-package (all platforms)
Install and load devtools, then install directly from github

```{R}
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)

install_github("cmkobel/tabseq")

# Then load the library in your script
library(tabseq)
```


# Installation of cli-tools (lunix)

Clone this repo, and expand your path to the `python/` directory.

Depends on python3.

```{sh}
cd ~
git clone git@github.com:cmkobel/tabseq.git

# Optionally, add the conversion scripts to your PATH-variable
echo "PATH=\$PATH:~/tabseq/python" >> ~/.bashrc
```

**Dependencies**:
 * git
 * python3
