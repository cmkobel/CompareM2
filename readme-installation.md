## Installation of AssemblyComparator2 (asscom2) on Linux

asscom2 can be installed in three steps: By downloading the code, setting up an alias in your user profile (~/.bashrc) that let's you launch the pipeline from any directory on your machine.

The only requisites for running asscom2 is that you have:
  - [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
  - git (can be installed with conda by typing `conda install -c anaconda git`)
  - [apptainer](https://apptainer.org/docs/user/main/quick_start.html#installation-request)

#### 0) First, check that you have the prerequisites available on your system:

```bash
which conda && conda --version
which git && git --version
which apptainer && apptainer --version
```

#### 1) Then download the asscom2 pipeline binary and set up an alias in your profile (~/.bashrc on most linux distributions)
```bash
git clone https://github.com/cmkobel/assemblycomparator2.git asscom2

cd asscom2
echo "export ASSCOM2_BASE=$(pwd)" >> ~/.bashrc

conda env create -n asscom2_launcher -f environment.yaml # Installs snakemake and mamba in a specific environment

```


#### 2) Finally, define the alias that will be used to launch asscom2 from any directory on your machine.
```bash
echo "alias asscom2='conda run --live-stream --name asscom2_launcher \
    snakemake \
        --snakefile \${ASSCOM2_BASE}/snakefile \
        --profile \${ASSCOM2_BASE}/apptainer/local \
        --configfile \${ASSCOM2_BASE}/config.yaml'" >> ~/.bashrc


```







