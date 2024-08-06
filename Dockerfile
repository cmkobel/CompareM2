FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="1012a21c168dba061542d315c54c19d5de65aa5c4817eedcaaf1733f1ad13ddc"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: dynamic_report/workflow/envs/r-markdown.yaml
#   prefix: /conda-envs/b06d3784eaabfc8abcf62a294f54bad9
#   name: r-markdown
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.1 #=4.2.2
#     - r-essentials
#     - conda-forge::r-tidyverse
#     - r-rmarkdown
#     - r-dt
#     - r-ape # For phylogenetic tree plotting. Maybe I should switch to ggtree?
#     
#     
#     #- r-prettydoc
#     #- r-rmarkdown
#     #- r-phytools
#     #- r-gridextra
#   
#   # Removing most of these dependencies fixes the #74 issue. So the plan now is to purge as many of them as possible and then build a new master docker image and see if it fixes the problem. KISS, no bloat, is the philosophy.
RUN mkdir -p /conda-envs/b06d3784eaabfc8abcf62a294f54bad9
COPY dynamic_report/workflow/envs/r-markdown.yaml /conda-envs/b06d3784eaabfc8abcf62a294f54bad9/environment.yaml

# Conda environment:
#   source: workflow/envs/abricate.yaml
#   prefix: /conda-envs/a77979bb28302ba15eb49ca8f93197f7
#   name: abricate
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - abricate=1
RUN mkdir -p /conda-envs/a77979bb28302ba15eb49ca8f93197f7
COPY workflow/envs/abricate.yaml /conda-envs/a77979bb28302ba15eb49ca8f93197f7/environment.yaml

# Conda environment:
#   source: workflow/envs/antismash.yaml
#   prefix: /conda-envs/2ae6d7c5468f4b0e569c09e067b92cfe
#   name: antismash
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - antismash=7
#   
#   # conda create -n antismash antismash
#   # conda activate antismash
#   # download-antismash-databases
#   # conda deactivate
RUN mkdir -p /conda-envs/2ae6d7c5468f4b0e569c09e067b92cfe
COPY workflow/envs/antismash.yaml /conda-envs/2ae6d7c5468f4b0e569c09e067b92cfe/environment.yaml

# Conda environment:
#   source: workflow/envs/any2fasta.yaml
#   prefix: /conda-envs/55f257b7a2d5325550dfe5652bd53c22
#   name: any2fasta
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - any2fasta=0
RUN mkdir -p /conda-envs/55f257b7a2d5325550dfe5652bd53c22
COPY workflow/envs/any2fasta.yaml /conda-envs/55f257b7a2d5325550dfe5652bd53c22/environment.yaml

# Conda environment:
#   source: workflow/envs/assembly-stats.yaml
#   prefix: /conda-envs/5c39a05942fff213bbc045e0d6dc61a3
#   name: assembly-stats
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - assembly-stats=1
RUN mkdir -p /conda-envs/5c39a05942fff213bbc045e0d6dc61a3
COPY workflow/envs/assembly-stats.yaml /conda-envs/5c39a05942fff213bbc045e0d6dc61a3/environment.yaml

# Conda environment:
#   source: workflow/envs/bakta.yaml
#   prefix: /conda-envs/349bbcba5bdfd7faab9097a8a3f7f6d2
#   name: bakta
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bakta=1
#   
#   # https://github.com/oschwengers/bakta?tab=readme-ov-file#bioconda
#   # conda install -c conda-forge -c bioconda bakta
RUN mkdir -p /conda-envs/349bbcba5bdfd7faab9097a8a3f7f6d2
COPY workflow/envs/bakta.yaml /conda-envs/349bbcba5bdfd7faab9097a8a3f7f6d2/environment.yaml

# Conda environment:
#   source: workflow/envs/busco.yaml
#   prefix: /conda-envs/29ddfde03b648d05bcbe1b488c9d275e
#   name: busco
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - busco=5.7
RUN mkdir -p /conda-envs/29ddfde03b648d05bcbe1b488c9d275e
COPY workflow/envs/busco.yaml /conda-envs/29ddfde03b648d05bcbe1b488c9d275e/environment.yaml

# Conda environment:
#   source: workflow/envs/checkm2.yaml
#   prefix: /conda-envs/cc9543e16cbf71e94c901e40623c2b25
#   name: checkm2_conda
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - checkm2=1
RUN mkdir -p /conda-envs/cc9543e16cbf71e94c901e40623c2b25
COPY workflow/envs/checkm2.yaml /conda-envs/cc9543e16cbf71e94c901e40623c2b25/environment.yaml

# Conda environment:
#   source: workflow/envs/dbcan.yaml
#   prefix: /conda-envs/8fd5cba1a415e69056b71dd889a203b8
#   name: dbcan
#   channels:
#     - anaconda
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.8
#     - dbcan=4 #.1.3 # I might want to fix this to 4.1.3 as I'm seeing problems with 4.1.4 ?
#     - anaconda::wget # Not sure if I should also add anaconda as a channel?
RUN mkdir -p /conda-envs/8fd5cba1a415e69056b71dd889a203b8
COPY workflow/envs/dbcan.yaml /conda-envs/8fd5cba1a415e69056b71dd889a203b8/environment.yaml

# Conda environment:
#   source: workflow/envs/eggnog.yaml
#   prefix: /conda-envs/5ab146b4f8a224935ba7154dde2d9065
#   name: eggnog
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bioconda::eggnog-mapper=2.1.12
#   
#   
#   # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Installation
#   # conda install -c bioconda -c conda-forge eggnog-mapper
RUN mkdir -p /conda-envs/5ab146b4f8a224935ba7154dde2d9065
COPY workflow/envs/eggnog.yaml /conda-envs/5ab146b4f8a224935ba7154dde2d9065/environment.yaml

# Conda environment:
#   source: workflow/envs/fasttree.yaml
#   prefix: /conda-envs/785d03a6b41f3872cfe5d325544fd5a7
#   name: fasttree
#   channels:
#     - bioconda
#   dependencies:
#     - fasttree=2
RUN mkdir -p /conda-envs/785d03a6b41f3872cfe5d325544fd5a7
COPY workflow/envs/fasttree.yaml /conda-envs/785d03a6b41f3872cfe5d325544fd5a7/environment.yaml

# Conda environment:
#   source: workflow/envs/gtdbtk.yaml
#   prefix: /conda-envs/df2a53a81d11717a831dda8590f53cd5
#   name: gtdbtk
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - gtdbtk=2.4
RUN mkdir -p /conda-envs/df2a53a81d11717a831dda8590f53cd5
COPY workflow/envs/gtdbtk.yaml /conda-envs/df2a53a81d11717a831dda8590f53cd5/environment.yaml

# Conda environment:
#   source: workflow/envs/interproscan.yaml
#   prefix: /conda-envs/8614b56b87888c570d4ed7d9eda3b5be
#   name: interproscan
#   channels: 
#     - bioconda
#   dependencies:
#     - interproscan=5
#   
#   # conda install -c bioconda interproscan
RUN mkdir -p /conda-envs/8614b56b87888c570d4ed7d9eda3b5be
COPY workflow/envs/interproscan.yaml /conda-envs/8614b56b87888c570d4ed7d9eda3b5be/environment.yaml

# Conda environment:
#   source: workflow/envs/iqtree.yaml
#   prefix: /conda-envs/31753de779143be0587a030432205123
#   name: iqtree
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - iqtree=2
RUN mkdir -p /conda-envs/31753de779143be0587a030432205123
COPY workflow/envs/iqtree.yaml /conda-envs/31753de779143be0587a030432205123/environment.yaml

# Conda environment:
#   source: workflow/envs/mashtree.yaml
#   prefix: /conda-envs/19c3e22137382d27c7ff2a7f802d1b99
#   name: mashtree
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - mashtree=1
RUN mkdir -p /conda-envs/19c3e22137382d27c7ff2a7f802d1b99
COPY workflow/envs/mashtree.yaml /conda-envs/19c3e22137382d27c7ff2a7f802d1b99/environment.yaml

# Conda environment:
#   source: workflow/envs/mlst.yaml
#   prefix: /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b
#   name: mlst
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults    
#   dependencies:
#     - mlst=2
RUN mkdir -p /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b
COPY workflow/envs/mlst.yaml /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b/environment.yaml

# Conda environment:
#   source: workflow/envs/panaroo.yaml
#   prefix: /conda-envs/75b0e966fbd3cc274cc2d44b43c772ea
#   name: panaroo
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.9
#     - panaroo>=1.3
#     
#   # https://gtonkinhill.github.io/panaroo/#/gettingstarted/installation
RUN mkdir -p /conda-envs/75b0e966fbd3cc274cc2d44b43c772ea
COPY workflow/envs/panaroo.yaml /conda-envs/75b0e966fbd3cc274cc2d44b43c772ea/environment.yaml

# Conda environment:
#   source: workflow/envs/prokka.yaml
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
COPY workflow/envs/prokka.yaml /conda-envs/b417c9814002b674416363d8b7b41887/environment.yaml

# Conda environment:
#   source: workflow/envs/r-clusterProfiler.yaml
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
COPY workflow/envs/r-clusterProfiler.yaml /conda-envs/c7ea168f1c2f891901d5ae5c900889dd/environment.yaml

# Conda environment:
#   source: workflow/envs/seqkit.yaml
#   prefix: /conda-envs/578271d635a58e711b2c022ba49a040e
#   name: seqkit
#   channels:
#     - bioconda
#   dependencies:
#     - seqkit=2
RUN mkdir -p /conda-envs/578271d635a58e711b2c022ba49a040e
COPY workflow/envs/seqkit.yaml /conda-envs/578271d635a58e711b2c022ba49a040e/environment.yaml

# Conda environment:
#   source: workflow/envs/snp-dists.yaml
#   prefix: /conda-envs/c003df5f674e91aedca01151ce198f86
#   name: snp-dists
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - snp-dists=0
RUN mkdir -p /conda-envs/c003df5f674e91aedca01151ce198f86
COPY workflow/envs/snp-dists.yaml /conda-envs/c003df5f674e91aedca01151ce198f86/environment.yaml

# Conda environment:
#   source: workflow/envs/treecluster.yaml
#   prefix: /conda-envs/460abd4b10197ad284c861f21f94f546
#   name: treecluster
#   channels: 
#     - defaults
#   dependencies:
#     - pip:
#       - treecluster
RUN mkdir -p /conda-envs/460abd4b10197ad284c861f21f94f546
COPY workflow/envs/treecluster.yaml /conda-envs/460abd4b10197ad284c861f21f94f546/environment.yaml

# Conda environment:
#   source: workflow/envs/wget.yaml
#   prefix: /conda-envs/3d07f847725ad902762b77365977b881
#   name: wget
#   channels:
#     - anaconda
#   dependencies:
#     - wget=1
RUN mkdir -p /conda-envs/3d07f847725ad902762b77365977b881
COPY workflow/envs/wget.yaml /conda-envs/3d07f847725ad902762b77365977b881/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/b06d3784eaabfc8abcf62a294f54bad9 --file /conda-envs/b06d3784eaabfc8abcf62a294f54bad9/environment.yaml && \
    mamba env create --prefix /conda-envs/a77979bb28302ba15eb49ca8f93197f7 --file /conda-envs/a77979bb28302ba15eb49ca8f93197f7/environment.yaml && \
    mamba env create --prefix /conda-envs/2ae6d7c5468f4b0e569c09e067b92cfe --file /conda-envs/2ae6d7c5468f4b0e569c09e067b92cfe/environment.yaml && \
    mamba env create --prefix /conda-envs/55f257b7a2d5325550dfe5652bd53c22 --file /conda-envs/55f257b7a2d5325550dfe5652bd53c22/environment.yaml && \
    mamba env create --prefix /conda-envs/5c39a05942fff213bbc045e0d6dc61a3 --file /conda-envs/5c39a05942fff213bbc045e0d6dc61a3/environment.yaml && \
    mamba env create --prefix /conda-envs/349bbcba5bdfd7faab9097a8a3f7f6d2 --file /conda-envs/349bbcba5bdfd7faab9097a8a3f7f6d2/environment.yaml && \
    mamba env create --prefix /conda-envs/29ddfde03b648d05bcbe1b488c9d275e --file /conda-envs/29ddfde03b648d05bcbe1b488c9d275e/environment.yaml && \
    mamba env create --prefix /conda-envs/cc9543e16cbf71e94c901e40623c2b25 --file /conda-envs/cc9543e16cbf71e94c901e40623c2b25/environment.yaml && \
    mamba env create --prefix /conda-envs/8fd5cba1a415e69056b71dd889a203b8 --file /conda-envs/8fd5cba1a415e69056b71dd889a203b8/environment.yaml && \
    mamba env create --prefix /conda-envs/5ab146b4f8a224935ba7154dde2d9065 --file /conda-envs/5ab146b4f8a224935ba7154dde2d9065/environment.yaml && \
    mamba env create --prefix /conda-envs/785d03a6b41f3872cfe5d325544fd5a7 --file /conda-envs/785d03a6b41f3872cfe5d325544fd5a7/environment.yaml && \
    mamba env create --prefix /conda-envs/df2a53a81d11717a831dda8590f53cd5 --file /conda-envs/df2a53a81d11717a831dda8590f53cd5/environment.yaml && \
    mamba env create --prefix /conda-envs/8614b56b87888c570d4ed7d9eda3b5be --file /conda-envs/8614b56b87888c570d4ed7d9eda3b5be/environment.yaml && \
    mamba env create --prefix /conda-envs/31753de779143be0587a030432205123 --file /conda-envs/31753de779143be0587a030432205123/environment.yaml && \
    mamba env create --prefix /conda-envs/19c3e22137382d27c7ff2a7f802d1b99 --file /conda-envs/19c3e22137382d27c7ff2a7f802d1b99/environment.yaml && \
    mamba env create --prefix /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b --file /conda-envs/9d60f0b222fb3b85a8e15fb48d16d15b/environment.yaml && \
    mamba env create --prefix /conda-envs/75b0e966fbd3cc274cc2d44b43c772ea --file /conda-envs/75b0e966fbd3cc274cc2d44b43c772ea/environment.yaml && \
    mamba env create --prefix /conda-envs/b417c9814002b674416363d8b7b41887 --file /conda-envs/b417c9814002b674416363d8b7b41887/environment.yaml && \
    mamba env create --prefix /conda-envs/c7ea168f1c2f891901d5ae5c900889dd --file /conda-envs/c7ea168f1c2f891901d5ae5c900889dd/environment.yaml && \
    mamba env create --prefix /conda-envs/578271d635a58e711b2c022ba49a040e --file /conda-envs/578271d635a58e711b2c022ba49a040e/environment.yaml && \
    mamba env create --prefix /conda-envs/c003df5f674e91aedca01151ce198f86 --file /conda-envs/c003df5f674e91aedca01151ce198f86/environment.yaml && \
    mamba env create --prefix /conda-envs/460abd4b10197ad284c861f21f94f546 --file /conda-envs/460abd4b10197ad284c861f21f94f546/environment.yaml && \
    mamba env create --prefix /conda-envs/3d07f847725ad902762b77365977b881 --file /conda-envs/3d07f847725ad902762b77365977b881/environment.yaml && \
    mamba clean --all -y
