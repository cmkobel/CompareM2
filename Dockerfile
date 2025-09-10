FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="0382cde30146d49c4e96ced9dd831418b2ac672f7368d69bfb0a5650f70417b7"

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
#   prefix: /conda-envs/6db8141981445cf76592e9714ec68471
#   name: amrfinder
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - ncbi-amrfinderplus =4
RUN mkdir -p /conda-envs/6db8141981445cf76592e9714ec68471
COPY workflow/envs/amrfinder.yaml /conda-envs/6db8141981445cf76592e9714ec68471/environment.yaml

# Conda environment:
#   source: workflow/envs/antismash.yaml
#   prefix: /conda-envs/3e17ccfe4ba16146395d22ee67a2b4bf
#   name: antismash
#   channels: 
#     - conda-forge
#     - bioconda
#     #- defaults
#   dependencies:
#     - antismash =7
#   
#   # conda create -n antismash antismash
#   # conda activate antismash
#   # download-antismash-databases
#   # conda deactivate
RUN mkdir -p /conda-envs/3e17ccfe4ba16146395d22ee67a2b4bf
COPY workflow/envs/antismash.yaml /conda-envs/3e17ccfe4ba16146395d22ee67a2b4bf/environment.yaml

# Conda environment:
#   source: workflow/envs/assembly-stats.yaml
#   prefix: /conda-envs/822f5402692dbf871bbe4c46418ebd90
#   name: assembly-stats
#   channels:
#     - conda-forge
#     - bioconda
#     #- defaults
#   dependencies:
#     - assembly-stats =1
RUN mkdir -p /conda-envs/822f5402692dbf871bbe4c46418ebd90
COPY workflow/envs/assembly-stats.yaml /conda-envs/822f5402692dbf871bbe4c46418ebd90/environment.yaml

# Conda environment:
#   source: workflow/envs/bakta.yaml
#   prefix: /conda-envs/54dddd6c96716032490ac50d79fbfca3
#   name: bakta
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bakta =1
#   
#   # https://github.com/oschwengers/bakta?tab=readme-ov-file#bioconda
#   # conda install -c conda-forge -c bioconda bakta
RUN mkdir -p /conda-envs/54dddd6c96716032490ac50d79fbfca3
COPY workflow/envs/bakta.yaml /conda-envs/54dddd6c96716032490ac50d79fbfca3/environment.yaml

# Conda environment:
#   source: workflow/envs/checkm2.yaml
#   prefix: /conda-envs/bdfdfd3801dae19174ec1853e1ee5c8d
#   name: checkm2_conda
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - checkm2 >=1.1 # Necessary for the current "14897628" database.
RUN mkdir -p /conda-envs/bdfdfd3801dae19174ec1853e1ee5c8d
COPY workflow/envs/checkm2.yaml /conda-envs/bdfdfd3801dae19174ec1853e1ee5c8d/environment.yaml

# Conda environment:
#   source: workflow/envs/dbcan.yaml
#   prefix: /conda-envs/0bdcada8f1eab9a3ce6fde9e0c178019
#   name: dbcan
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python =3.8
#     - dbcan =4 #.1.3 # I might want to fix this to 4.1.3 as I'm seeing problems with 4.1.4 ?
#     - conda-forge::wget # Not sure if I should also add anaconda as a channel?
RUN mkdir -p /conda-envs/0bdcada8f1eab9a3ce6fde9e0c178019
COPY workflow/envs/dbcan.yaml /conda-envs/0bdcada8f1eab9a3ce6fde9e0c178019/environment.yaml

# Conda environment:
#   source: workflow/envs/eggnog.yaml
#   prefix: /conda-envs/2493d9ca1cf103c1a5d590eabc55ee62
#   name: eggnog
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bioconda::eggnog-mapper =2.1.12
#   
#   # Eggnog-mapper (output files) seems to change a lot, so I'm fixating the version quite tightly.
#   
#   # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Installation
#   # conda install -c bioconda -c conda-forge eggnog-mapper
RUN mkdir -p /conda-envs/2493d9ca1cf103c1a5d590eabc55ee62
COPY workflow/envs/eggnog.yaml /conda-envs/2493d9ca1cf103c1a5d590eabc55ee62/environment.yaml

# Conda environment:
#   source: workflow/envs/fasttree.yaml
#   prefix: /conda-envs/16b1c29308bfbc29a57953ab2bb80d8f
#   name: fasttree
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - fasttree =2
RUN mkdir -p /conda-envs/16b1c29308bfbc29a57953ab2bb80d8f
COPY workflow/envs/fasttree.yaml /conda-envs/16b1c29308bfbc29a57953ab2bb80d8f/environment.yaml

# Conda environment:
#   source: workflow/envs/gapseq.yaml
#   prefix: /conda-envs/69a625824fd518a54832a144b821a75c
#   name: gapseq
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - gapseq
RUN mkdir -p /conda-envs/69a625824fd518a54832a144b821a75c
COPY workflow/envs/gapseq.yaml /conda-envs/69a625824fd518a54832a144b821a75c/environment.yaml

# Conda environment:
#   source: workflow/envs/gtdbtk.yaml
#   prefix: /conda-envs/7c44445e5cb959587ca3a58c9ef12b14
#   name: gtdbtk
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - gtdbtk =2.4
RUN mkdir -p /conda-envs/7c44445e5cb959587ca3a58c9ef12b14
COPY workflow/envs/gtdbtk.yaml /conda-envs/7c44445e5cb959587ca3a58c9ef12b14/environment.yaml

# Conda environment:
#   source: workflow/envs/interproscan.yaml
#   prefix: /conda-envs/f6fc38804e167a2fe16ed0c7d67cc1e7
#   name: interproscan
#   channels: 
#     - conda-forge
#     - bioconda
#   dependencies:
#     - interproscan =5
#   
#   # conda install -c bioconda interproscan
RUN mkdir -p /conda-envs/f6fc38804e167a2fe16ed0c7d67cc1e7
COPY workflow/envs/interproscan.yaml /conda-envs/f6fc38804e167a2fe16ed0c7d67cc1e7/environment.yaml

# Conda environment:
#   source: workflow/envs/iqtree.yaml
#   prefix: /conda-envs/ac4a767610baaaa9a655c4a9913ea512
#   name: iqtree
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - iqtree =2
RUN mkdir -p /conda-envs/ac4a767610baaaa9a655c4a9913ea512
COPY workflow/envs/iqtree.yaml /conda-envs/ac4a767610baaaa9a655c4a9913ea512/environment.yaml

# Conda environment:
#   source: workflow/envs/mashtree.yaml
#   prefix: /conda-envs/e3a79faa3f4cd94133c27c0a9454bf52
#   name: mashtree
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - mashtree =1
RUN mkdir -p /conda-envs/e3a79faa3f4cd94133c27c0a9454bf52
COPY workflow/envs/mashtree.yaml /conda-envs/e3a79faa3f4cd94133c27c0a9454bf52/environment.yaml

# Conda environment:
#   source: workflow/envs/mlst.yaml
#   prefix: /conda-envs/4063a2f2fb09255885f1c9b7f22a9ad7
#   name: mlst
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - mlst =2
RUN mkdir -p /conda-envs/4063a2f2fb09255885f1c9b7f22a9ad7
COPY workflow/envs/mlst.yaml /conda-envs/4063a2f2fb09255885f1c9b7f22a9ad7/environment.yaml

# Conda environment:
#   source: workflow/envs/ncbi_datasets.yaml
#   prefix: /conda-envs/caa70cfc9d0c61df7c2685fd3870bbae
#   name: ncbi_datasets
#   channels:
#     - conda-forge  
#   dependencies:
#     - conda-forge::ncbi-datasets-cli =18
#     - conda-forge::unzip
RUN mkdir -p /conda-envs/caa70cfc9d0c61df7c2685fd3870bbae
COPY workflow/envs/ncbi_datasets.yaml /conda-envs/caa70cfc9d0c61df7c2685fd3870bbae/environment.yaml

# Conda environment:
#   source: workflow/envs/panaroo.yaml
#   prefix: /conda-envs/6e9bfa118513c9018de18e1ac9dcfdcb
#   name: panaroo
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python =3.9
#     - panaroo >=1.3,<2
#     
#   # https://gtonkinhill.github.io/panaroo/#/gettingstarted/installation
RUN mkdir -p /conda-envs/6e9bfa118513c9018de18e1ac9dcfdcb
COPY workflow/envs/panaroo.yaml /conda-envs/6e9bfa118513c9018de18e1ac9dcfdcb/environment.yaml

# Conda environment:
#   source: workflow/envs/prokka.yaml
#   prefix: /conda-envs/5dbc3d7ce225bd1b7c2b1f47616e16fc
#   name: prokka
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - prokka =1
#     - openjdk <=17.0.2
RUN mkdir -p /conda-envs/5dbc3d7ce225bd1b7c2b1f47616e16fc
COPY workflow/envs/prokka.yaml /conda-envs/5dbc3d7ce225bd1b7c2b1f47616e16fc/environment.yaml

# Conda environment:
#   source: workflow/envs/r-clusterProfiler.yaml
#   prefix: /conda-envs/e89e20aa37ae0fd6a678586e1f1be456
#   name: r-clusterProfiler
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bioconductor-clusterprofiler =4
#     - r-tidyverse
#     - r-jsonlite
RUN mkdir -p /conda-envs/e89e20aa37ae0fd6a678586e1f1be456
COPY workflow/envs/r-clusterProfiler.yaml /conda-envs/e89e20aa37ae0fd6a678586e1f1be456/environment.yaml

# Conda environment:
#   source: workflow/envs/seqkit.yaml
#   prefix: /conda-envs/43923cc46cd18e9bdafa97781e938664
#   name: seqkit
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - seqkit =2
RUN mkdir -p /conda-envs/43923cc46cd18e9bdafa97781e938664
COPY workflow/envs/seqkit.yaml /conda-envs/43923cc46cd18e9bdafa97781e938664/environment.yaml

# Conda environment:
#   source: workflow/envs/snp-dists.yaml
#   prefix: /conda-envs/533a4cc9b04e69ccca21da322eead69b
#   name: snp-dists
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - snp-dists =0
RUN mkdir -p /conda-envs/533a4cc9b04e69ccca21da322eead69b
COPY workflow/envs/snp-dists.yaml /conda-envs/533a4cc9b04e69ccca21da322eead69b/environment.yaml

# Conda environment:
#   source: workflow/envs/treecluster.yaml
#   prefix: /conda-envs/20c1c58fb46486fe37cc7ea8d66b726f
#   name: treecluster
#   channels: 
#     - conda-forge
#     - bioconda
#   dependencies:
#     - pip:
#       - treecluster
RUN mkdir -p /conda-envs/20c1c58fb46486fe37cc7ea8d66b726f
COPY workflow/envs/treecluster.yaml /conda-envs/20c1c58fb46486fe37cc7ea8d66b726f/environment.yaml

# Conda environment:
#   source: workflow/envs/wget.yaml
#   prefix: /conda-envs/7b6ee336e45b07e84a6bf1104b851a24
#   name: wget
#   channels:
#     - conda-forge
#   dependencies:
#     - wget =1
RUN mkdir -p /conda-envs/7b6ee336e45b07e84a6bf1104b851a24
COPY workflow/envs/wget.yaml /conda-envs/7b6ee336e45b07e84a6bf1104b851a24/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/b5c6b988e4a325065292b3d680193d32 --file /conda-envs/b5c6b988e4a325065292b3d680193d32/environment.yaml && \
    mamba env create --prefix /conda-envs/6db8141981445cf76592e9714ec68471 --file /conda-envs/6db8141981445cf76592e9714ec68471/environment.yaml && \
    mamba env create --prefix /conda-envs/3e17ccfe4ba16146395d22ee67a2b4bf --file /conda-envs/3e17ccfe4ba16146395d22ee67a2b4bf/environment.yaml && \
    mamba env create --prefix /conda-envs/822f5402692dbf871bbe4c46418ebd90 --file /conda-envs/822f5402692dbf871bbe4c46418ebd90/environment.yaml && \
    mamba env create --prefix /conda-envs/54dddd6c96716032490ac50d79fbfca3 --file /conda-envs/54dddd6c96716032490ac50d79fbfca3/environment.yaml && \
    mamba env create --prefix /conda-envs/bdfdfd3801dae19174ec1853e1ee5c8d --file /conda-envs/bdfdfd3801dae19174ec1853e1ee5c8d/environment.yaml && \
    mamba env create --prefix /conda-envs/0bdcada8f1eab9a3ce6fde9e0c178019 --file /conda-envs/0bdcada8f1eab9a3ce6fde9e0c178019/environment.yaml && \
    mamba env create --prefix /conda-envs/2493d9ca1cf103c1a5d590eabc55ee62 --file /conda-envs/2493d9ca1cf103c1a5d590eabc55ee62/environment.yaml && \
    mamba env create --prefix /conda-envs/16b1c29308bfbc29a57953ab2bb80d8f --file /conda-envs/16b1c29308bfbc29a57953ab2bb80d8f/environment.yaml && \
    mamba env create --prefix /conda-envs/69a625824fd518a54832a144b821a75c --file /conda-envs/69a625824fd518a54832a144b821a75c/environment.yaml && \
    mamba env create --prefix /conda-envs/7c44445e5cb959587ca3a58c9ef12b14 --file /conda-envs/7c44445e5cb959587ca3a58c9ef12b14/environment.yaml && \
    mamba env create --prefix /conda-envs/f6fc38804e167a2fe16ed0c7d67cc1e7 --file /conda-envs/f6fc38804e167a2fe16ed0c7d67cc1e7/environment.yaml && \
    mamba env create --prefix /conda-envs/ac4a767610baaaa9a655c4a9913ea512 --file /conda-envs/ac4a767610baaaa9a655c4a9913ea512/environment.yaml && \
    mamba env create --prefix /conda-envs/e3a79faa3f4cd94133c27c0a9454bf52 --file /conda-envs/e3a79faa3f4cd94133c27c0a9454bf52/environment.yaml && \
    mamba env create --prefix /conda-envs/4063a2f2fb09255885f1c9b7f22a9ad7 --file /conda-envs/4063a2f2fb09255885f1c9b7f22a9ad7/environment.yaml && \
    mamba env create --prefix /conda-envs/caa70cfc9d0c61df7c2685fd3870bbae --file /conda-envs/caa70cfc9d0c61df7c2685fd3870bbae/environment.yaml && \
    mamba env create --prefix /conda-envs/6e9bfa118513c9018de18e1ac9dcfdcb --file /conda-envs/6e9bfa118513c9018de18e1ac9dcfdcb/environment.yaml && \
    mamba env create --prefix /conda-envs/5dbc3d7ce225bd1b7c2b1f47616e16fc --file /conda-envs/5dbc3d7ce225bd1b7c2b1f47616e16fc/environment.yaml && \
    mamba env create --prefix /conda-envs/e89e20aa37ae0fd6a678586e1f1be456 --file /conda-envs/e89e20aa37ae0fd6a678586e1f1be456/environment.yaml && \
    mamba env create --prefix /conda-envs/43923cc46cd18e9bdafa97781e938664 --file /conda-envs/43923cc46cd18e9bdafa97781e938664/environment.yaml && \
    mamba env create --prefix /conda-envs/533a4cc9b04e69ccca21da322eead69b --file /conda-envs/533a4cc9b04e69ccca21da322eead69b/environment.yaml && \
    mamba env create --prefix /conda-envs/20c1c58fb46486fe37cc7ea8d66b726f --file /conda-envs/20c1c58fb46486fe37cc7ea8d66b726f/environment.yaml && \
    mamba env create --prefix /conda-envs/7b6ee336e45b07e84a6bf1104b851a24 --file /conda-envs/7b6ee336e45b07e84a6bf1104b851a24/environment.yaml && \
    mamba clean --all -y
