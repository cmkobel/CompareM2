FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="6b78d3f87417ac666505c938d49127fa0b016be51e6e47917bd7d2c7bddecca6"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: conda_definitions/abricate.yaml
#   prefix: /conda-envs/a77979bb28302ba15eb49ca8f93197f7
#   name: abricate
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - abricate=1
RUN mkdir -p /conda-envs/a77979bb28302ba15eb49ca8f93197f7
COPY conda_definitions/abricate.yaml /conda-envs/a77979bb28302ba15eb49ca8f93197f7/environment.yaml

# Conda environment:
#   source: conda_definitions/any2fasta.yaml
#   prefix: /conda-envs/55f257b7a2d5325550dfe5652bd53c22
#   name: any2fasta
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - any2fasta=0
RUN mkdir -p /conda-envs/55f257b7a2d5325550dfe5652bd53c22
COPY conda_definitions/any2fasta.yaml /conda-envs/55f257b7a2d5325550dfe5652bd53c22/environment.yaml

# Conda environment:
#   source: conda_definitions/assembly-stats.yaml
#   prefix: /conda-envs/5c39a05942fff213bbc045e0d6dc61a3
#   name: assembly-stats
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - assembly-stats=1
RUN mkdir -p /conda-envs/5c39a05942fff213bbc045e0d6dc61a3
COPY conda_definitions/assembly-stats.yaml /conda-envs/5c39a05942fff213bbc045e0d6dc61a3/environment.yaml

# Conda environment:
#   source: conda_definitions/busco.yaml
#   prefix: /conda-envs/29ddfde03b648d05bcbe1b488c9d275e
#   name: busco
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - busco=5.7
RUN mkdir -p /conda-envs/29ddfde03b648d05bcbe1b488c9d275e
COPY conda_definitions/busco.yaml /conda-envs/29ddfde03b648d05bcbe1b488c9d275e/environment.yaml

# Conda environment:
#   source: conda_definitions/checkm2.yaml
#   prefix: /conda-envs/cc9543e16cbf71e94c901e40623c2b25
#   name: checkm2_conda
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - checkm2=1
RUN mkdir -p /conda-envs/cc9543e16cbf71e94c901e40623c2b25
COPY conda_definitions/checkm2.yaml /conda-envs/cc9543e16cbf71e94c901e40623c2b25/environment.yaml

# Conda environment:
#   source: conda_definitions/dbcan.yaml
#   prefix: /conda-envs/5029d7d5859d34f300af826f9c1e7465
#   name: dbcan
#   channels:
#     - anaconda
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.8
#     - dbcan=4
#     - anaconda::wget # Not sure if I should also add anaconda as a channel?
RUN mkdir -p /conda-envs/5029d7d5859d34f300af826f9c1e7465
COPY conda_definitions/dbcan.yaml /conda-envs/5029d7d5859d34f300af826f9c1e7465/environment.yaml

# Conda environment:
#   source: conda_definitions/diamond.yaml
#   prefix: /conda-envs/86344782fe4f6e7e8853a33d355eac38
#   name: diamond
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - diamond=2.1.8
RUN mkdir -p /conda-envs/86344782fe4f6e7e8853a33d355eac38
COPY conda_definitions/diamond.yaml /conda-envs/86344782fe4f6e7e8853a33d355eac38/environment.yaml

# Conda environment:
#   source: conda_definitions/fasttree.yaml
#   prefix: /conda-envs/785d03a6b41f3872cfe5d325544fd5a7
#   name: fasttree
#   channels:
#     - bioconda
#   dependencies:
#     - fasttree=2
RUN mkdir -p /conda-envs/785d03a6b41f3872cfe5d325544fd5a7
COPY conda_definitions/fasttree.yaml /conda-envs/785d03a6b41f3872cfe5d325544fd5a7/environment.yaml

# Conda environment:
#   source: conda_definitions/gtdbtk.yaml
#   prefix: /conda-envs/ea337c41a8bd47af414ed24d5d2cfced
#   name: gtdbtk
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - gtdbtk=2.3
RUN mkdir -p /conda-envs/ea337c41a8bd47af414ed24d5d2cfced
COPY conda_definitions/gtdbtk.yaml /conda-envs/ea337c41a8bd47af414ed24d5d2cfced/environment.yaml

# Conda environment:
#   source: conda_definitions/interproscan.yaml
#   prefix: /conda-envs/8614b56b87888c570d4ed7d9eda3b5be
#   name: interproscan
#   channels: 
#     - bioconda
#   dependencies:
#     - interproscan=5
#   
#   # conda install -c bioconda interproscan
RUN mkdir -p /conda-envs/8614b56b87888c570d4ed7d9eda3b5be
COPY conda_definitions/interproscan.yaml /conda-envs/8614b56b87888c570d4ed7d9eda3b5be/environment.yaml

# Conda environment:
#   source: conda_definitions/iqtree.yaml
#   prefix: /conda-envs/31753de779143be0587a030432205123
#   name: iqtree
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - iqtree=2
RUN mkdir -p /conda-envs/31753de779143be0587a030432205123
COPY conda_definitions/iqtree.yaml /conda-envs/31753de779143be0587a030432205123/environment.yaml

# Conda environment:
#   source: conda_definitions/mashtree.yaml
#   prefix: /conda-envs/19c3e22137382d27c7ff2a7f802d1b99
#   name: mashtree
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - mashtree=1
RUN mkdir -p /conda-envs/19c3e22137382d27c7ff2a7f802d1b99
COPY conda_definitions/mashtree.yaml /conda-envs/19c3e22137382d27c7ff2a7f802d1b99/environment.yaml

# Conda environment:
#   source: conda_definitions/mlst.yaml
#   prefix: /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b
#   name: mlst
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults    
#   dependencies:
#     - mlst=2
RUN mkdir -p /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b
COPY conda_definitions/mlst.yaml /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b/environment.yaml

# Conda environment:
#   source: conda_definitions/motulizer.yaml
#   prefix: /conda-envs/bd4163ffe28cb2a044a74e425b76c45d
#   name: motulizer
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - motulizer=0
RUN mkdir -p /conda-envs/bd4163ffe28cb2a044a74e425b76c45d
COPY conda_definitions/motulizer.yaml /conda-envs/bd4163ffe28cb2a044a74e425b76c45d/environment.yaml

# Conda environment:
#   source: conda_definitions/prokka.yaml
#   prefix: /conda-envs/b417c9814002b674416363d8b7b41887
#   name: prokka
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults # not sure about whether this one should be present
#   dependencies:
#     - prokka=1
#     - openjdk<=17.0.2
RUN mkdir -p /conda-envs/b417c9814002b674416363d8b7b41887
COPY conda_definitions/prokka.yaml /conda-envs/b417c9814002b674416363d8b7b41887/environment.yaml

# Conda environment:
#   source: conda_definitions/r-clusterProfiler.yaml
#   prefix: /conda-envs/c7ea168f1c2f891901d5ae5c900889dd
#   name: r-clusterProfiler
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bioconductor-clusterprofiler=4
#     - r-tidyverse
#     - r-jsonlite
RUN mkdir -p /conda-envs/c7ea168f1c2f891901d5ae5c900889dd
COPY conda_definitions/r-clusterProfiler.yaml /conda-envs/c7ea168f1c2f891901d5ae5c900889dd/environment.yaml

# Conda environment:
#   source: conda_definitions/roary_see-comments-in-this-file.yaml
#   prefix: /conda-envs/8de1358f988451556029ab664d67c6d0
#   name: roary_new
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#     - defaults
#   dependencies:
#     - roary=3
#   
#   
#   
#   # Warning: This only works if you set the channel priority to flexible. Not strict.
#   # conda config --set channel_priority false
#   # 
#   # Install roary
#   #
#   # Then set it back with
#   # conda config --set channel_priority strict
RUN mkdir -p /conda-envs/8de1358f988451556029ab664d67c6d0
COPY conda_definitions/roary_see-comments-in-this-file.yaml /conda-envs/8de1358f988451556029ab664d67c6d0/environment.yaml

# Conda environment:
#   source: conda_definitions/seqkit.yaml
#   prefix: /conda-envs/578271d635a58e711b2c022ba49a040e
#   name: seqkit
#   channels:
#     - bioconda
#   dependencies:
#     - seqkit=2
RUN mkdir -p /conda-envs/578271d635a58e711b2c022ba49a040e
COPY conda_definitions/seqkit.yaml /conda-envs/578271d635a58e711b2c022ba49a040e/environment.yaml

# Conda environment:
#   source: conda_definitions/snp-dists.yaml
#   prefix: /conda-envs/c003df5f674e91aedca01151ce198f86
#   name: snp-dists
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - snp-dists=0
RUN mkdir -p /conda-envs/c003df5f674e91aedca01151ce198f86
COPY conda_definitions/snp-dists.yaml /conda-envs/c003df5f674e91aedca01151ce198f86/environment.yaml

# Conda environment:
#   source: conda_definitions/wget.yaml
#   prefix: /conda-envs/3d07f847725ad902762b77365977b881
#   name: wget
#   channels:
#     - anaconda
#   dependencies:
#     - wget=1
RUN mkdir -p /conda-envs/3d07f847725ad902762b77365977b881
COPY conda_definitions/wget.yaml /conda-envs/3d07f847725ad902762b77365977b881/environment.yaml

# Conda environment:
#   source: report_subpipeline/conda_definitions/r-markdown.yaml
#   prefix: /conda-envs/2e8da1d04478713e88851b1dbb6cce04
#   name: r-markdown
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#   dependencies:
#     - r-base=4.2 #=4.2.2
#     - r-essentials
#     - r-tidyverse
#     - r-dt
#     - r-ape # For phylogenetic tree plotting. Maybe I should switch to ggtree?
#     
#     #- r-prettydoc
#     #- r-rmarkdown
#     #- r-phytools
#     #- r-gridextra
#   
#   # Removing most of these dependencies fixes the #74 issue. So the plan now is to purge as many of them as possible and then build a new master docker image and see if it fixes the problem. KISS, no bloat, is the philosophy.
RUN mkdir -p /conda-envs/2e8da1d04478713e88851b1dbb6cce04
COPY report_subpipeline/conda_definitions/r-markdown.yaml /conda-envs/2e8da1d04478713e88851b1dbb6cce04/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/a77979bb28302ba15eb49ca8f93197f7 --file /conda-envs/a77979bb28302ba15eb49ca8f93197f7/environment.yaml && \
    mamba env create --prefix /conda-envs/55f257b7a2d5325550dfe5652bd53c22 --file /conda-envs/55f257b7a2d5325550dfe5652bd53c22/environment.yaml && \
    mamba env create --prefix /conda-envs/5c39a05942fff213bbc045e0d6dc61a3 --file /conda-envs/5c39a05942fff213bbc045e0d6dc61a3/environment.yaml && \
    mamba env create --prefix /conda-envs/29ddfde03b648d05bcbe1b488c9d275e --file /conda-envs/29ddfde03b648d05bcbe1b488c9d275e/environment.yaml && \
    mamba env create --prefix /conda-envs/cc9543e16cbf71e94c901e40623c2b25 --file /conda-envs/cc9543e16cbf71e94c901e40623c2b25/environment.yaml && \
    mamba env create --prefix /conda-envs/5029d7d5859d34f300af826f9c1e7465 --file /conda-envs/5029d7d5859d34f300af826f9c1e7465/environment.yaml && \
    mamba env create --prefix /conda-envs/86344782fe4f6e7e8853a33d355eac38 --file /conda-envs/86344782fe4f6e7e8853a33d355eac38/environment.yaml && \
    mamba env create --prefix /conda-envs/785d03a6b41f3872cfe5d325544fd5a7 --file /conda-envs/785d03a6b41f3872cfe5d325544fd5a7/environment.yaml && \
    mamba env create --prefix /conda-envs/ea337c41a8bd47af414ed24d5d2cfced --file /conda-envs/ea337c41a8bd47af414ed24d5d2cfced/environment.yaml && \
    mamba env create --prefix /conda-envs/8614b56b87888c570d4ed7d9eda3b5be --file /conda-envs/8614b56b87888c570d4ed7d9eda3b5be/environment.yaml && \
    mamba env create --prefix /conda-envs/31753de779143be0587a030432205123 --file /conda-envs/31753de779143be0587a030432205123/environment.yaml && \
    mamba env create --prefix /conda-envs/19c3e22137382d27c7ff2a7f802d1b99 --file /conda-envs/19c3e22137382d27c7ff2a7f802d1b99/environment.yaml && \
    mamba env create --prefix /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b --file /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b/environment.yaml && \
    mamba env create --prefix /conda-envs/bd4163ffe28cb2a044a74e425b76c45d --file /conda-envs/bd4163ffe28cb2a044a74e425b76c45d/environment.yaml && \
    mamba env create --prefix /conda-envs/b417c9814002b674416363d8b7b41887 --file /conda-envs/b417c9814002b674416363d8b7b41887/environment.yaml && \
    mamba env create --prefix /conda-envs/c7ea168f1c2f891901d5ae5c900889dd --file /conda-envs/c7ea168f1c2f891901d5ae5c900889dd/environment.yaml && \
    mamba env create --prefix /conda-envs/8de1358f988451556029ab664d67c6d0 --file /conda-envs/8de1358f988451556029ab664d67c6d0/environment.yaml && \
    mamba env create --prefix /conda-envs/578271d635a58e711b2c022ba49a040e --file /conda-envs/578271d635a58e711b2c022ba49a040e/environment.yaml && \
    mamba env create --prefix /conda-envs/c003df5f674e91aedca01151ce198f86 --file /conda-envs/c003df5f674e91aedca01151ce198f86/environment.yaml && \
    mamba env create --prefix /conda-envs/3d07f847725ad902762b77365977b881 --file /conda-envs/3d07f847725ad902762b77365977b881/environment.yaml && \
    mamba env create --prefix /conda-envs/2e8da1d04478713e88851b1dbb6cce04 --file /conda-envs/2e8da1d04478713e88851b1dbb6cce04/environment.yaml && \
    mamba clean --all -y
