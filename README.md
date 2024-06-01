# Assemblycomparator2
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/cmkobel/assemblycomparator2/dry-run.yaml)](https://github.com/cmkobel/assemblycomparator2/actions/) [![Docker Pulls](https://img.shields.io/docker/pulls/cmkobel/assemblycomparator2?label=docker%20pulls)](https://assemblycomparator2.readthedocs.io/en/latest/10%20installation/)  [![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/assemblycomparator2?label=Bioconda%20downloads&color=%2300CC00)](https://assemblycomparator2.readthedocs.io/en/latest/10%20installation/) [![conda build](https://img.shields.io/conda/v/bioconda/assemblycomparator2)](https://assemblycomparator2.readthedocs.io/en/latest/10%20installation/) [![Documentation Status](https://readthedocs.org/projects/assemblycomparator2/badge/?version=latest)](https://assemblycomparator2.readthedocs.io/en/latest/?badge=latest)

🧬 Assemblycomparator2 "asscom2" is a genomes-to-report pipeline. It accepts prokaryotic (bacterial and archaeal) genomic assemblies and compares them in many different ways. 

🦠 Being designed to analyze assemblies of both isolates and metagenomes (MAGs), it is useful for anyone working with microbial genomics.

💾 [Installing](https://assemblycomparator2.readthedocs.io/en/latest/10%20installation/) Assemblycomparator2 on your system gives you access to [many](https://assemblycomparator2.readthedocs.io/en/latest/what%20analyses%20does%20it%20do/) powerful state-of-the-art tools for analysis of prokaryotic genomes which will accelerate your research. It is easy to use and can be used by non-bioinformaticians.

<img alt="asscom2 animation" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/623f6b42-2de6-457c-8f0d-3b3e5d646967">


👩‍🔬 Assemblycomparator2 integrates several analyses that yield scientific results about genomic assemblies on several levels: Quality control, annotation, function and species calling as well as comparative analyses like computation of core/pan genomes and phylogenetics. 

<img width="150" alt="snakemake logo" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/7188e748-9d37-43ae-a5d5-100e9560df1f">

🐍 Assemblycomparator2 works by calling a Snakemake workflow that can be easily modified to use different parameters for the  underlying tools.

<a href="https://assemblycomparator2.readthedocs.io/en/latest/20%20usage/#demo-reports"><img height="192" alt="report document logo" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/e5f9b72c-2137-4850-8779-a5528d8ccbaf"></a>

📙 Central results are dynamically integrated in a compact portable report .html-document. It can be browsed in any web browser and can be easily shared as a single file. This report is generated even if some jobs in the pipeline fail. See [examples](https://assemblycomparator2.readthedocs.io/en/latest/20%20usage/#demo-reports).

🧑‍💻 Assemblycomparator2 can be run either on a local workstation (recommended >= 64GiB RAM), or a HPC (high performance computing) cluster. Both  Apptainer/Singularity/Docker images and conda environment definitions are available for all dependent software to run.


**The comprehensible documentation is available at [assemblycomparator2.readthedocs.org](https://assemblycomparator2.readthedocs.org).**

Please log any issues you have in the [issues tab of the repository](https://github.com/cmkobel/assemblycomparator2/issues).

Development is active and will continue. We're actively looking for collaborators to join and synergize this project. Pull requests are warmly encouraged.

---

Assemblycomparator2 genomes to report pipeline. Copyright (C) 2024 [Carl M. Kobel](https://github.com/cmkobel) [GNU GPL v3](https://github.com/cmkobel/assemblycomparator2/blob/master/LICENSE)

