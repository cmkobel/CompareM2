FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="22cb0e992026e0a51e7865cde3fe4f121ebee6884b03236720bb5eac7c59f150"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: dynamic_report/workflow/envs/r-markdown.yaml
#   prefix: /conda-envs/b5c6b988e4a325065292b3d680193d32
#   name: r-markdown
#   channels:
#     - conda-forge
#     - bioconda
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
RUN mkdir -p /conda-envs/b5c6b988e4a325065292b3d680193d32
COPY dynamic_report/workflow/envs/r-markdown.yaml /conda-envs/b5c6b988e4a325065292b3d680193d32/environment.yaml

# Conda environment:
#   source: workflow/envs/amrfinder.yaml
#   prefix: /conda-envs/98af290ab46d66d8dd29fcb7342c19b1
#   name: amrfinder
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - ncbi-amrfinderplus=4
RUN mkdir -p /conda-envs/98af290ab46d66d8dd29fcb7342c19b1
COPY workflow/envs/amrfinder.yaml /conda-envs/98af290ab46d66d8dd29fcb7342c19b1/environment.yaml

# Conda environment:
#   source: workflow/envs/antismash.yaml
#   prefix: /conda-envs/4ae94d2c4b9a52f46cfad5eb92750982
#   name: antismash
#   channels: 
#     - conda-forge
#     - bioconda
#     #- defaults
#   dependencies:
#     - antismash=7
#   
#   # conda create -n antismash antismash
#   # conda activate antismash
#   # download-antismash-databases
#   # conda deactivate
RUN mkdir -p /conda-envs/4ae94d2c4b9a52f46cfad5eb92750982
COPY workflow/envs/antismash.yaml /conda-envs/4ae94d2c4b9a52f46cfad5eb92750982/environment.yaml

# Conda environment:
#   source: workflow/envs/assembly-stats.yaml
#   prefix: /conda-envs/b116edf022cd78b810f7c18e1d4cb568
#   name: assembly-stats
#   channels:
#     - conda-forge
#     - bioconda
#     #- defaults
#   dependencies:
#     - assembly-stats=1
RUN mkdir -p /conda-envs/b116edf022cd78b810f7c18e1d4cb568
COPY workflow/envs/assembly-stats.yaml /conda-envs/b116edf022cd78b810f7c18e1d4cb568/environment.yaml

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
#   prefix: /conda-envs/06095a97a5b2a1ad8552af1df6312c10
#   name: dbcan
#   channels:
#     #- anaconda
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.8
#     - dbcan=4 #.1.3 # I might want to fix this to 4.1.3 as I'm seeing problems with 4.1.4 ?
#     - conda-forge::wget # Not sure if I should also add anaconda as a channel?
RUN mkdir -p /conda-envs/06095a97a5b2a1ad8552af1df6312c10
COPY workflow/envs/dbcan.yaml /conda-envs/06095a97a5b2a1ad8552af1df6312c10/environment.yaml

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
#   prefix: /conda-envs/5bf95c616c0ba0db0478590fc40d2505
#   name: fasttree
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - fasttree=2
RUN mkdir -p /conda-envs/5bf95c616c0ba0db0478590fc40d2505
COPY workflow/envs/fasttree.yaml /conda-envs/5bf95c616c0ba0db0478590fc40d2505/environment.yaml

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
#   prefix: /conda-envs/7d12a67dde56d160f55937194337f2c4
#   name: interproscan
#   channels: 
#     - conda-forge
#     - bioconda
#   dependencies:
#     - interproscan=5
#   
#   # conda install -c bioconda interproscan
RUN mkdir -p /conda-envs/7d12a67dde56d160f55937194337f2c4
COPY workflow/envs/interproscan.yaml /conda-envs/7d12a67dde56d160f55937194337f2c4/environment.yaml

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
#   prefix: /conda-envs/441fdf65234f80a1472c9e3dabc23e99
#   name: mlst
#   channels:
#     - conda-forge
#     - bioconda
#     #- defaults    
#   dependencies:
#     - mlst=2
RUN mkdir -p /conda-envs/441fdf65234f80a1472c9e3dabc23e99
COPY workflow/envs/mlst.yaml /conda-envs/441fdf65234f80a1472c9e3dabc23e99/environment.yaml

# Conda environment:
#   source: workflow/envs/panaroo.yaml
#   prefix: /conda-envs/ee827f9c5c29bae85e8f39537571ae58
#   name: panaroo
#   channels:
#     - conda-forge
#     - bioconda
#     #- defaults
#   dependencies:
#     - python=3.9
#     - panaroo>=1.3
#     
#   # https://gtonkinhill.github.io/panaroo/#/gettingstarted/installation
RUN mkdir -p /conda-envs/ee827f9c5c29bae85e8f39537571ae58
COPY workflow/envs/panaroo.yaml /conda-envs/ee827f9c5c29bae85e8f39537571ae58/environment.yaml

# Conda environment:
#   source: workflow/envs/prokka.yaml
#   prefix: /conda-envs/6c1c1cc7bda74d875fab714481417a27
#   name: prokka
#   channels:
#     - conda-forge
#     - bioconda
#     #- defaults # not sure about whether this one should be present
#   dependencies:
#     - prokka=1
#     - openjdk<=17.0.2
RUN mkdir -p /conda-envs/6c1c1cc7bda74d875fab714481417a27
COPY workflow/envs/prokka.yaml /conda-envs/6c1c1cc7bda74d875fab714481417a27/environment.yaml

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
#   prefix: /conda-envs/87a18131a84b118b7e4841b35afee564
#   name: seqkit
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - seqkit=2
RUN mkdir -p /conda-envs/87a18131a84b118b7e4841b35afee564
COPY workflow/envs/seqkit.yaml /conda-envs/87a18131a84b118b7e4841b35afee564/environment.yaml

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

RUN mamba env create --prefix /conda-envs/b5c6b988e4a325065292b3d680193d32 --file /conda-envs/b5c6b988e4a325065292b3d680193d32/environment.yaml && \
    mamba env create --prefix /conda-envs/98af290ab46d66d8dd29fcb7342c19b1 --file /conda-envs/98af290ab46d66d8dd29fcb7342c19b1/environment.yaml && \
    mamba env create --prefix /conda-envs/4ae94d2c4b9a52f46cfad5eb92750982 --file /conda-envs/4ae94d2c4b9a52f46cfad5eb92750982/environment.yaml && \
    mamba env create --prefix /conda-envs/b116edf022cd78b810f7c18e1d4cb568 --file /conda-envs/b116edf022cd78b810f7c18e1d4cb568/environment.yaml && \
    mamba env create --prefix /conda-envs/349bbcba5bdfd7faab9097a8a3f7f6d2 --file /conda-envs/349bbcba5bdfd7faab9097a8a3f7f6d2/environment.yaml && \
    mamba env create --prefix /conda-envs/cc9543e16cbf71e94c901e40623c2b25 --file /conda-envs/cc9543e16cbf71e94c901e40623c2b25/environment.yaml && \
    mamba env create --prefix /conda-envs/06095a97a5b2a1ad8552af1df6312c10 --file /conda-envs/06095a97a5b2a1ad8552af1df6312c10/environment.yaml && \
    mamba env create --prefix /conda-envs/5ab146b4f8a224935ba7154dde2d9065 --file /conda-envs/5ab146b4f8a224935ba7154dde2d9065/environment.yaml && \
    mamba env create --prefix /conda-envs/5bf95c616c0ba0db0478590fc40d2505 --file /conda-envs/5bf95c616c0ba0db0478590fc40d2505/environment.yaml && \
    mamba env create --prefix /conda-envs/df2a53a81d11717a831dda8590f53cd5 --file /conda-envs/df2a53a81d11717a831dda8590f53cd5/environment.yaml && \
    mamba env create --prefix /conda-envs/7d12a67dde56d160f55937194337f2c4 --file /conda-envs/7d12a67dde56d160f55937194337f2c4/environment.yaml && \
    mamba env create --prefix /conda-envs/31753de779143be0587a030432205123 --file /conda-envs/31753de779143be0587a030432205123/environment.yaml && \
    mamba env create --prefix /conda-envs/19c3e22137382d27c7ff2a7f802d1b99 --file /conda-envs/19c3e22137382d27c7ff2a7f802d1b99/environment.yaml && \
    mamba env create --prefix /conda-envs/441fdf65234f80a1472c9e3dabc23e99 --file /conda-envs/441fdf65234f80a1472c9e3dabc23e99/environment.yaml && \
    mamba env create --prefix /conda-envs/ee827f9c5c29bae85e8f39537571ae58 --file /conda-envs/ee827f9c5c29bae85e8f39537571ae58/environment.yaml && \
    mamba env create --prefix /conda-envs/6c1c1cc7bda74d875fab714481417a27 --file /conda-envs/6c1c1cc7bda74d875fab714481417a27/environment.yaml && \
    mamba env create --prefix /conda-envs/c7ea168f1c2f891901d5ae5c900889dd --file /conda-envs/c7ea168f1c2f891901d5ae5c900889dd/environment.yaml && \
    mamba env create --prefix /conda-envs/87a18131a84b118b7e4841b35afee564 --file /conda-envs/87a18131a84b118b7e4841b35afee564/environment.yaml && \
    mamba env create --prefix /conda-envs/c003df5f674e91aedca01151ce198f86 --file /conda-envs/c003df5f674e91aedca01151ce198f86/environment.yaml && \
    mamba env create --prefix /conda-envs/460abd4b10197ad284c861f21f94f546 --file /conda-envs/460abd4b10197ad284c861f21f94f546/environment.yaml && \
    mamba env create --prefix /conda-envs/3d07f847725ad902762b77365977b881 --file /conda-envs/3d07f847725ad902762b77365977b881/environment.yaml && \
    mamba clean --all -y
