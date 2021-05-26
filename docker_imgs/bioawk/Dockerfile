FROM continuumio/miniconda3

ENV LC_ALL en_US.UTF-8



WORKDIR /app

# Make RUN commands use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]




RUN conda install -y -c conda-forge -c bioconda -c defaults bioawk


