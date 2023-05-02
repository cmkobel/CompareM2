# In string you define the names (without .yaml) of all of the environments. Then the script generates the dockerfiles.
string="abricate any2fasta assembly-stats bioawk busco checkm2_conda curl fasttree gtdbtk kraken2 mashtree mlst prokka r-clusterProfiler roary_see-comments-in-this-file roary r-tidyverse seqkit snp-dists"

for i in $string; do 
    echo $i

    echo """FROM mambaorg/micromamba:1.4.2
COPY --chown=\$MAMBA_USER:\$MAMBA_USER conda_definitions/$i.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \\
    micromamba clean --all --yes

## Inspired heavily from
# https://github.com/mamba-org/micromamba-docker#quick-start""" > Dockerfile_$i


    echo 
done    