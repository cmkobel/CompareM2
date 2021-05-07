
# Define (manually) where you want the kraken2 database to reside
kraken2_database_dir="~/kraken2_database_dir"

# Preemptively create and enter this dir
mkdir -p $kraken2_database_dir
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20210127.tar.gz
tar -xvf k2_pluspf_20210127.tar.gz

cd k2_pluspf_20210127

echo "export ASSCOM_KRAKEN2_DB='$(wd)'" >> ~/.bashrc
source ~/.bashrc