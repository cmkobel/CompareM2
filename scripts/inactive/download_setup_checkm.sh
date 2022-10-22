#if [ ! -f {params.directory}/checkm_OK.flag ]; then    
#
#        wget --directory-prefix={params.directory} https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz 
#        tar -xvf {params.directory}/checkm_data_2015_01_16.tar.gz -C {params.directory}
#        touch {params.directory}/checkm_OK.flag
#
#    fi
