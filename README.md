# CompareM2
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/cmkobel/comparem2/dry-run.yaml)](https://github.com/cmkobel/CompareM2/actions/workflows/dry-run.yaml) [![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/comparem2?label=Bioconda%20downloads&color=%2300CC00)](https://comparem2.readthedocs.io/en/latest/10%20installation/) [![Docker Pulls](https://img.shields.io/docker/pulls/cmkobel/comparem2?label=docker%20pulls)](https://comparem2.readthedocs.io/en/latest/10%20installation/) [![Documentation Status](https://readthedocs.org/projects/comparem2/badge/?version=latest)](https://comparem2.readthedocs.io/en/latest/?badge=latest) [![Conda Version](https://img.shields.io/conda/v/bioconda/comparem2)](https://anaconda.org/bioconda/comparem2) [![https://doi.org/10.1101/2024.07.12.603264](https://img.shields.io/badge/DOI-10.1101%2F2024.07.12.603264-blue.svg)](https://doi.org/10.1101/2024.07.12.603264) [![GitHub Repo stars](https://img.shields.io/github/stars/cmkobel/comparem2)](https://github.com/cmkobel/comparem2)

!!! note
    If you're looking for the original version of CompareM, a tool to calculate AAI and codon usage, please follow this link: [github.com/donovan-h-parks/CompareM](https://github.com/donovan-h-parks/CompareM)
    


🧬 CompareM2 is a genomes-to-report pipeline. It accepts prokaryotic (bacterial and archaeal) genomic assemblies and compares them in many different ways. 

🦠 Being designed to analyze assemblies of both isolates and metagenomes (MAGs), it is useful for anyone working with microbial genomics.

💾 [Installing](https://comparem2.readthedocs.io/en/latest/10%20installation/) CompareM2 on your system gives you access to many powerful state-of-the-art tools for analysis of prokaryotic genomes which will accelerate your research. It is easy to use and can be used by non-bioinformaticians.

👩‍🔬 CompareM2 integrates [several analyses](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/) that yield scientific results about genomic assemblies on several levels: Quality control, annotation, function and species calling as well as comparative analyses like computation of core/pan genomes and phylogenetics. 

🐍 CompareM2 works by calling a Snakemake workflow that can be easily modified to use [different parameters](https://comparem2.readthedocs.io/en/latest/20%20usage/#passthrough-arguments) for the  underlying tools.

<a href="https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#rendered-report"><img width="220" style="width: 220px" alt="report document logo" align="right" src="https://github.com/cmkobel/comparem2/assets/5913696/e5f9b72c-2137-4850-8779-a5528d8ccbaf"></a>

📄 Central results are dynamically integrated in a compact portable report .html-document. It can be browsed in any web browser and can be easily shared as a single file. This report is generated even if some jobs in the pipeline fail. [See examples](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#rendered-report).

🧑‍💻 CompareM2 can be run either on a local workstation (recommended >= 64GiB RAM), or a HPC (high performance computing) cluster. Both  Apptainer/Singularity/Docker images and conda environment definitions are available for all dependent software to run.

🙋 If you have any questions, issues or ideas about using CompareM2, please raise an issue [here](https://github.com/cmkobel/CompareM2/issues).

📙 **The comprehensive documentation is available at [CompareM2.readthedocs.io](https://comparem2.readthedocs.io).**

---

[CompareM2](https://github.com/cmkobel/comparem2) genomes-to-report pipeline. Copyright (C) 2024 [C. M. Kobel](https://github.com/cmkobel) GNU GPL v3


