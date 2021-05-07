#!/bin/bash

set -e

# Define (manually) where you want the kraken2 database to reside
#kraken2_database_dir=~/kraken2_database_dir
kraken2_database_dir=~/ClinicalMicrobio/faststorage/database/kraken2/k2_pluspf_20210127

# Preemptively create and enter this dir
echo "Creating and entering ${kraken2_database_dir}"
mkdir -p $kraken2_database_dir
cd $kraken2_database_dir

echo "Downloading database ..."
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20210127.tar.gz
tar -xvf k2_pluspf_20210127.tar.gz

#cd k2_pluspf_20210127

echo "Setting environment varible to current WD ..."
echo $(pwd)
echo "export ASSCOM_KRAKEN2_DB='$(pwd)'" >> ~/.bashrc

source ~/.bashrc


