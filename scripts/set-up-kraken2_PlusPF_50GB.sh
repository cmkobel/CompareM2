#!/bin/bash

#################################################################################################################
# This script downloads a kraken2 database, and sets a system variable so that assemblycomparator2 can find it. #
#################################################################################################################

# This flag makes sure the script stops if an error arises
set -e

# Select a title from this page: https://benlangmead.github.io/aws-indexes/k2
#database_title='k2_standard_8gb_20201202'   # Standard-8    Standard with DB capped at 8 GB
database_title='k2_pluspf_20210127'         # PlusPF        Standard plus protozoa & fungi (fixed from 12/2/20 version2) (50GB)
echo title $database_title

# Where to put the database
#base_dir=~/ClinicalMicrobio/faststorage/database/kraken2
base_dir=~/databases/kraken2
echo basedir $base_dir

# The web source to download
source=https://genome-idx.s3.amazonaws.com/kraken/${database_title}.tar.gz
echo source $source

# Where the database will end up
destination=${base_dir}/${database_title}
echo destination $destination
echo

# Preemptively create and enter the destination.
echo "Creating and entering ${destination}"
mkdir -p $destination

echo "Clearing the destination ..."
cd $destination && rm -f *.tar.gz
echo

echo "Downloading database ..."

wget https://genome-idx.s3.amazonaws.com/kraken/${database_title}.tar.gz
tar -xvf ${database_title}.tar.gz


echo "Setting environment variable to database destination: $(pwd)"
echo "export ASSCOM2_KRAKEN2_DB='$(pwd)'" >> ~/.bashrc




