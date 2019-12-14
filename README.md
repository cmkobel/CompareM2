```
┌─┐┌─┐┌─┐┌─┐┌┬┐┌┐ ┬ ┬ ┬  ┌─┐┌─┐┌┬┐┌─┐┌─┐┬─┐┌─┐┌┬┐┌─┐┬─┐KMA 
├─┤└─┐└─┐├┤ │││├┴┐│ └┬┘  │  │ ││││├─┘├─┤├┬┘├─┤ │ │ │├┬┘
┴ ┴└─┘└─┘└─┘┴ ┴└─┘┴─┘┴   └─┘└─┘┴ ┴┴  ┴ ┴┴└─┴ ┴ ┴ └─┘┴└─
```

The KMA Assembly Comparator is a simple pipeline consisting of independent tools for the analysis of microbial genomics. The idea is that the user provides a directory of assemblies of related bacterial strains, and the pipeline will run a series of generally relevant analyses.

The KMA Assembly Comparator supports all formats: (gb,fa,fq,gff,gfa,clw,sth) which can all be (.gz,bz2,zip)-compressed.

# Usage
In order to use the KMA Assebly Comparator, you have to install it such that you can call the main file (bash) in any directory on your system. For instance
```
export PATH=/path/to/assemblycomparator/:$PATH
```

Now, create a directory anywhere on your system containing the fasta-formatted assemblies you want to analyse, and run the main program `assemblycomparator`.

Then, the program will list a series of options, helping you to start the pipeline.


# Development

This project is still in development

## Phase 1
Get it working

A few subanalyses still need to be implemented, (see parenthesized programs in the below list)

## Phase 2
Get it working with a shared environment


### For each assembly
  * any2fasta (wide format support)
  * prokka (annotation)
  * kraken2 (species identification)
  * mlst (multi locus sequence typing)
  * (abricate) (virulence/resistance gene identification)
  * (snippy) (pairwise snp-calling)
  
  
### For each group
  * roary (pan and core genome)
  * panito (average nucleotide identity
  * FastTree (phylogenetic tree of core genome)
  * (IQ-tree) (phylogenetic tree of core genome)
  
  
  
  
