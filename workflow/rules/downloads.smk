
# --- Downloads -----------------------------------------------------------------

# This rule runs once, downloading the busco dataset that is needed for rule busco.
# Make sure that this job is run on a node that has internet access.
rule busco_download:
    output:
        #touch("{base_variable}/databases/busco/busco_download_done.flag") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
        #database_representative = touch("{base_variable}/databases/busco/file_versions.tsv") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
        database_representative = touch(DATABASES + "/busco/comparem2_busco_database_representative.flag") # Should point to the directory where the following files reside: "file_versions.tsv  lineages/  placement_files/"
    conda: "../envs/busco.yaml"
    shell: """

        # If some previous batch of comparem2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output:q}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            # Busco is a bit stupid in the way that it requires an input file, but doesn't read it when you just download.
            touch dummy.fasta
            
            # https://busco.ezlab.org/busco_userguide.html#download-and-automated-update
            busco \
                --in dummy.fasta \
                --out dummy \
                --mode geno \
                --auto-lineage-prok \
                --force \
                --download_path $(dirname {output:q}) \
                --download prokaryota
            
            
            touch {output:q}

            # Clean up 
            rm dummy.fasta
        
        fi

    """



# Updated according to chklovski's idea in https://github.com/chklovski/CheckM2/issues/73#issuecomment-1744207103
rule checkm2_download:
    output:
        database_representative = DATABASES + "/checkm2/comparem2_checkm2_database_representative.flag",
    params:
        destination = DATABASES + "/checkm2"
    conda: "../envs/wget.yaml"
    shell: """


        # If some previous batch of comparem2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output:q}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            url="https://zenodo.org/records/5571251/files/checkm2_database.tar.gz"

            wget --no-check-certificate -O "{params.destination}/checkm2_database.tar.gz" "$url"

            tar \
                -xvf "{params.destination}/checkm2_database.tar.gz" \
                --directory "{params.destination}"

            rm "{params.destination}/checkm2_database.tar.gz" || echo "failed to clean up"
            
            touch {output:q}
        
        fi

    """




rule dbcan_download:
    output:
        #database_representative = touch("{base_variable}/databases/dbcan/comparem2_dbcan_database_representative.flag"),
        database_representative = DATABASES + "/dbcan/comparem2_dbcan_database_representative.flag",
    conda: "../envs/dbcan.yaml"
    threads: 8
    shell: """
        
        # If some previous batch of comparem2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output:q}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."            

            
            # Testing the new automatic downloader.
            dbcan_build \
                --cpus {threads} \
                --db-dir "$(dirname {output.database_representative})" \
                --clean

            # Comments on using the download code from https://github.com/linnabrown/run_dbcan (june 2023): I deleted the test ecoli files in the bottom, and added --continue, to make sure that not a .1 suffixed file is left over when retrying downloads.



            >&2 echo "dbcan setup completed"
            echo "Downloaded dbcan at $(date -Iseconds)" > $(dirname {output.database_representative})/info.txt

            mkdir -p $(dirname {output:q})
            touch {output:q}

        fi

    """




rule gtdb_download:
    output:
        database_representative = DATABASES + "/gtdb/comparem2_gtdb_database_representative.flag"
    conda: "../envs/wget.yaml"
    shell: """

        # https://ecogenomics.github.io/GTDBTk/installing/index.html

        # Pick a source file
        # Release 214
        # db_pick="https://ns9864k.web.sigma2.no/TheMEMOgroup/cmkobel/comparem2-assets/gtdb/release214/gtdbtk_data.tar.gz" # NMBU/MEMO mirror in norway
        #db_pick="https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz" # Official location
        #db_pick="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz" # Official alternative mirror
        
        # Release 220
        db_pick="https://ns9864k.web.sigma2.no/TheMEMOgroup/cmkobel/comparem2-assets/gtdb/release220/gtdbtk_r220_data.tar.gz" # NMBU/MEMO mirror in norway
        #db_pick="https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
        #db_pick="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"


        db_destination=$(dirname {output.database_representative})/gtdb_db.tar.gz

        # If some previous batch of comparem2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output:q}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            >&2 echo "Downlading $db_pick to $db_destination"
            mkdir -p $(dirname "$db_destination")

            wget --no-check-certificate -O "$db_destination" "$db_pick"

            >&2 echo "Decompressing ..."
            tar \
                -xf $db_destination \
                --directory $(dirname $db_destination)

            rm $db_destination || echo "Failed to clean up." 

            >&2 echo "gtdb DB setup completed"
            echo "Downloaded $db_pick at $(date -Iseconds)" > $(dirname $db_destination)/info.txt

            mkdir -p $(dirname {output:q})
            touch {output:q}

        fi

    """






rule bakta_download:
    output:
        database_representative = DATABASES + "/bakta/comparem2_bakta_database_representative.flag"
    conda: "../envs/bakta.yaml"
    shell: """

        # https://github.com/oschwengers/bakta?tab=readme-ov-file#database-download

        # If some previous batch of comparem2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output:q}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."
            bakta --version
            bakta_db download \
                --output $(dirname {output.database_representative}) \
                --type full # or light
                
            # Run Bakta using '--db /spaceface/shared_databases/comparem2_v2.5.5+/bakta/db' or set a BAKTA_DB environment variable: 'export BAKTA_DB=/spaceface/shared_databases/comparem2_v2.5.5+/bakta/db'

            
            >&2 echo "bakta DB setup completed"
            echo "Downloaded at $(date -Iseconds)" > $(dirname {output.database_representative})/info.txt

            touch {output:q}

        fi

    """




rule eggnog_download:
    output:
        database_representative = DATABASES + "/eggnog/comparem2_eggnog_database_representative.flag"
    conda: "../envs/eggnog.yaml"
    shell: """

        # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#setup

        # If some previous batch of comparem2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output:q}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."
            download_eggnog_data.py \
                -y \
                -f \
                --data_dir $(dirname {output.database_representative}) 
            
            >&2 echo "eggnog DB setup completed"
            echo "Downloaded at $(date -Iseconds)" > $(dirname {output.database_representative})/info.txt

            touch {output:q}

        fi

    """





rule antismash_download:
    output:
        database_representative = DATABASES + "/antismash/comparem2_antismash_database_representative.flag"
    conda: "../envs/antismash.yaml"
    shell: """

        # https://docs.antismash.secondarymetabolites.org/install/

        # If some previous batch of comparem2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output:q}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."
            download-antismash-databases \
                --database-dir $(dirname {output.database_representative:q})
            
            >&2 echo "antismash DB setup completed"
            echo "Downloaded at $(date -Iseconds)" > $(dirname {output.database_representative})/info.txt

            touch {output:q}

        fi

    """




