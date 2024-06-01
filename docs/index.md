# Assemblycomparator2

[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/cmkobel/assemblycomparator2/dry-run.yaml)](https://github.com/cmkobel/assemblycomparator2/actions/) [![Docker Pulls](https://img.shields.io/docker/pulls/cmkobel/assemblycomparator2?label=docker%20pulls)](#installation)  [![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/assemblycomparator2?label=Bioconda%20downloads&color=%2300CC00)](#installation) [![conda build](https://img.shields.io/conda/v/bioconda/assemblycomparator2)](#installation)

🧬 Assemblycomparator2 "asscom2" is a genomes-to-report pipeline. It accepts prokaryotic (bacterial and archaeal) genomic assemblies and compares them in many different ways. 

🦠 Being designed to analyze assemblies of both isolates and metagenomes (MAGs), it is useful for anyone working with microbial genomics.

💾 [Installing](#installation) Assemblycomparator2 on your system gives you access to [many](https://github.com/cmkobel/assemblycomparator2#what-analyses-does-it-do) powerful state-of-the-art tools for analysis of prokaryotic genomes which will accelerate your research. It is easy to use and can be used by non-bioinformaticians.

<img alt="asscom2 animation" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/623f6b42-2de6-457c-8f0d-3b3e5d646967">



👩‍🔬 Assemblycomparator2 integrates several analyses that yield scientific results about genomic assemblies on several levels: Quality control, annotation, function and species calling as well as comparative analyses like computation of core/pan genomes and phylogenetics. 

<img width="150" alt="snakemake logo" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/7188e748-9d37-43ae-a5d5-100e9560df1f">

🐍 Assemblycomparator2 works by calling a Snakemake workflow that can be easily modified to use different parameters for the  underlying tools.

<a href="https://github.com/cmkobel/assemblycomparator2/blob/master/readme-demos.md"><img height="192" alt="report document logo" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/e5f9b72c-2137-4850-8779-a5528d8ccbaf"></a>

📙 Central results are dynamically integrated in a compact portable report .html-document. It can be browsed in any web browser and can be easily shared as a single file. This report is generated even if some jobs in the pipeline fail. See [examples](readme-demos.md).

🧑‍💻 Assemblycomparator2 can be run either on a local workstation (recommended >= 64GiB RAM), or a HPC (high performance computing) cluster. Both  Apptainer/Singularity/Docker images and conda environment definitions are available for all dependent software to run.

Assemblycomparator2 will -depending on circumstances- be renamed to "Genomisk" or something completely different.

## Usage examples

Make a directory with the prokaryotic genomic assembly-files -or metagenomic bins- you want to investigate with Assemblycomparator2. 
Go into that directory in the terminal, and run the command `asscom2`. 
Assemblycomparator2 will then create a sub-directory, named "results_ac2/" containing a plethora of analysis results. 
  
  - Execute a "dry run". That is, to show what will be run without actually doing it.

    ```
    asscom2 --dry-run
    ```

  - Run Assemblycomparator2 on the genomes in the current directory:

    ```
    asscom2
    ```
    

##### A bit more advanced controls 

  - Run analyses that are relevant to metagenomic assemblies only (as opposed to isolates):

    ```
    asscom2 --until meta
    ```
    
  - Execute all jobs until one or more specific rules: (until implies including)
    
    ```
    asscom2 --until roary abricate
    ```
    
  - Select a specific MLST-scheme to use on all of the samples: (defaults to automatic)
    
    ```
    asscom2 --config mlst_scheme=hpylori
    ```
    
  - Select a specific roary blastp-identity: (default is 95)

    ```
    asscom2 --config roary_blastp_identity=90
    ```
      