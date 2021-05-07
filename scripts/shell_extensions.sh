#!/bin/bash



# And save it into your .bashrc
echo "export ASSCOM2_BASE=$ASSCOM2_BASE" >> ~/.bashrc 

# If you are going to use conda, it is a good idea to set the SNAKEMAKE_CONDA_PREFIX-variable, so the package installations can be reused between runs.
echo "export SNAKEMAKE_CONDA_PREFIX=${ASSCOM2_BASE}/conda_base" >> ~/.bashrc 

# For local setups (using Conda for jobs):
echo "alias assemblycomparator2='conda activate assemblycomparator2; snakemake --snakefile ${ASSCOM2_BASE}/snakefile --profile ${ASSCOM2_BASE}/configs/local/ --use-conda'" >> ~/.bashrc

# For HPC's with Slurm (using Singularity for jobs):
echo "alias assemblycomparator2='conda activate assemblycomparator2; snakemake --snakefile ${ASSCOM2_BASE}/snakefile --profile ${ASSCOM2_BASE}/configs/slurm/ --cluster-config ${ASSCOM2_BASE}/configs/cluster.yaml --use-singularity'" >> ~/.bashrc