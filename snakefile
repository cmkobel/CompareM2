__author__ = 'Carl M. Kobel'

__version__ = "v2.5.12" 
# Places to bump besides here.
#  - changelog
#  - ./asscom2 binary
#  - report markdown bottom
#
#  And if the dockerfile is updated
#  - snakefile Dockerfile image pull
#  - report snakefile Dockerfile image pull so it reuses the same already downloaded image.


# May the data passing through this pipeline
# somehow help to bring just a little more peace 
# in this troubled world.

# Update Dockerfile with this command
# asscom2 --containerize > Dockerfile # And remove header-text.

# Update dag picture in documentation with this command (with anaconda/graphviz)
# asscom2 --forceall --rulegraph | dot -Tpdf > dag.pdf


import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
from shutil import copyfile

# For void_report
import subprocess
import datetime
        
# Newer version of asscom2 still uses the docker image v2.5.5
containerized: "docker://cmkobel/assemblycomparator2:v2.5.12" # Remember to copy the same version to the report_subpipeline/snakefile. I wonder if I can put this in the profile or config instead? 

# When executing, Snakemake will fail with a reasonable error message if the variables below are undefined.
envvars:
    "ASSCOM2_BASE",
    "ASSCOM2_PROFILE",
    "ASSCOM2_DATABASES"
    # What about parametrization of databases?

cwd = os.getcwd()
batch_title = cwd.split("/")[-1]
base_variable = os.environ['ASSCOM2_BASE'] # rename to ASSCOM2_BASE
DATABASES = os.environ['ASSCOM2_DATABASES'] # Defines where the databases are stored. One for all. when snakemake issue 262 is solved I'll make this more flexible for each rule.



print("/*") # for .dot exports used to generate dag visualizations.
print(f"                                                                   {__version__}")
print("         █████╗ ███████╗███████╗ ██████╗ ██████╗ ███╗   ███╗██████╗ ")
print("        ██╔══██╗██╔════╝██╔════╝██╔════╝██╔═══██╗████╗ ████║╚════██╗")
print("        ███████║███████╗███████╗██║     ██║   ██║██╔████╔██║ █████╔╝")
print("        ██╔══██║╚════██║╚════██║██║     ██║   ██║██║╚██╔╝██║██╔═══╝ ")
print("        ██║  ██║███████║███████║╚██████╗╚██████╔╝██║ ╚═╝ ██║███████╗")
print("        ╚═╝  ╚═╝╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝")
print("                       A.K.A. assemblycomparator2                   ")
print("                         Please log issues at:                      ")
print("              github.com/cmkobel/assemblycomparator2/issues         ")
print("                                                                    ")
print(f"    batch_title:           {batch_title}")
print(f"    base_variable:         {base_variable}                               ")
print(f"    databases:             {DATABASES}                             ")
print(f"    roary_blastp_identity: {config['roary_blastp_identity']} (default 95)")
print(f"    mlst_scheme:           {config['mlst_scheme']} (default automatic)   ")
#print(f"    kraken2 database:      {config['asscom2_kraken2_db']}                ")
#print(f"    gtdb:                  {config['gtdbtk_data_path']}                  ")



print()
print("    Available rules:")
print("    sequence_lengths prokka kraken2 dbcan interproscan kofam_scan")
print("    busco checkm2 diamond_kegg kegg_pathway roary snp_dists")
print("    assembly_stats gtdbtk abricate mlst mashtree fasttree")
print("    report downloads fast")

#results_directory = "output_asscom2"



results_directory = "results_ac2"



#reference = config["reference"]



# --- Read in relevant files in the current working directory ----------------

relative_wd = "."
#extension_whitelist = ["fna", "fa", "fas", "fasta", "seq"] # old 
extension_whitelist = ["fna", "fa", "fas", "fasta", "seq", "gb", "fq", "gff", "gfa", "clw", "sth", "gz", "bz2"] # From any2fasta

present_files = [f for f in listdir(relative_wd) if isfile(join(relative_wd,f))]


df = pd.DataFrame(data = {'input_file': present_files})

# Check that the directory is not empty.
if df.shape[0] == 0:
    raise Exception("Error: No fasta files in the current directory. Quitting ...")
    #raise Exception("Zero genomic files present.")
    sys.exit(1)



df = df[~df["input_file"].str.startswith(".", na = False)] # Remove hidden files
df['sample_raw'] = [".".join(i.split(".")[:-1]) for i in df['input_file'].tolist()] # Extract everything before the extension dot.
df['sample'] = df['sample_raw'].str.replace(' ','_').str.replace(',','_').str.replace('"','_').str.replace('\'','_') # Convert punctuation marks to underscores: Makes everything easier.
df['extension'] =  [i.split(".")[-1] for i in df['input_file'].tolist()] # Extract extension
df['input_file_fasta'] = results_directory + "/samples/" + df['sample'] + "/" + df['sample'] + ".fa" # This is where the input file is copied to in the first snakemake rule.

df = df[df['extension'].isin(extension_whitelist)] # Remove files with unsupported formats.
df['1-index'] = [i+1 for i in range(len(df))]

# Check that the directory is not empty, again.
if df.shape[0] == 0:
    raise Exception("Error: No fasta files in the current directory. Quitting ...(2)")
    #raise Exception("Zero genomic files present.")
    sys.exit(1)


#df_mini = df_mini.apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # strip whitespace and replace spaces with underscores.

  
# --- Displaying filtered dataframe ready for analysis ------------------------
print() # Padding
df = df.reset_index(drop = True)
#print(df[['input_file', 'sample', 'extension']])
#print(df[['input_file', 'extension', 'input_file_fasta']])
print(df[['1-index', 'sample', 'extension']].to_string(index = False))
print("//")
print()




# The DATABASES directory must exist, otherwise apptainer gets confused and throws the following:
# WARNING: skipping mount of /home/thylakoid/assemblycomparator2/adatabaseas: stat /home/thylakoid/assemblycomparator2/adatabaseas: no such file or directory
if not os.path.isdir(DATABASES):
    os.mkdir(DATABASES)


# --- Make sure the output directory exists. ----------------------------------
if not os.path.isdir(results_directory):
    os.mkdir(results_directory) # If running with local profile, the directory won't be created. This is necessary in the edge case that the user _only_ runs "--until report".



# The modification time of this file tells the report subpipeline whether it needs to run. Thus, void_report is called in the end of every successful rule.
#void_report = f"touch {results_directory}/.asscom2_void_report.flag"
void_report = f"date -Iseconds >> {results_directory}/.asscom2_void_report.flag"
#annotate_log_ = "2>&1 | while IFS= read -r line; do printf '[%s] %s\n' \"$(date '+%Y-%m-%d %H:%M:%S')\" \"$line\"; done"




localrules: metadata, checkm2_download, kraken2_download, dbcan_download, busco_download, gtdb_download, report, install_report_environment_aot

# --- Collect all targets. ------------------------------------------
rule all:
    input: expand([\
        "{results_directory}/metadata.tsv", \
        "{results_directory}/.install_report_environment_aot.flag", \
        "{results_directory}/assembly-stats/assembly-stats.tsv", \
        "{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", \
        "{results_directory}/samples/{sample}/busco/short_summary_extract.tsv", \
        "{results_directory}/checkm2/quality_report.tsv", \
        "{results_directory}/samples/{sample}/diamond_kegg/{sample}_diamond_kegg.tsv", \
        "{results_directory}/kegg_pathway/kegg_pathway_enrichment_analysis.tsv", \
        "{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv", \
        "{results_directory}/samples/{sample}/dbcan/overview.txt", \
        "{results_directory}/samples/{sample}/interproscan/{sample}_interproscan.tsv", \
        "{results_directory}/gtdbtk/gtdbtk.summary.tsv", \
        "{results_directory}/mlst/mlst.tsv", \
        "{results_directory}/abricate/card_detailed.tsv", \
        "{results_directory}/samples/{sample}/prokka/{sample}.gff", \
        "{results_directory}/roary/summary_statistics.txt", \
        "{results_directory}/fasttree/fasttree.newick", \
        "{results_directory}/iqtree/core_genome_iqtree.treefile", \
        "{results_directory}/snp-dists/snp-dists.tsv", \
        "{results_directory}/motulizer/motulizer_results.tsv", \
        "{results_directory}/motulizer/motupan_results.tsv", \
        "{results_directory}/mashtree/mashtree.newick"], \
        results_directory = results_directory, sample = df["sample"]) 

        
        # temporary disabled
        #"{results_directory}/roary/summary_statistics.txt", \
        #"{results_directory}/fasttree/fasttree.newick", \
        #"{results_directory}/snp-dists/snp-dists.tsv", \







# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
# Should I rename this to "rule any2fasta" just to make it more transparent?
rule copy:
    input: 
        genome = lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0],
    output: "{results_directory}/samples/{sample}/{sample}.fa"
    conda: "conda_definitions/any2fasta.yaml"
    threads: 1 # Weirdly, or bugly, there must be a thread n definition in the rule. Otherwise, the set-threads option (in the orion profile) will not be taken up. 
    resources:
        mem_mb = 256,
        runtime = "10m",
    shell: """

        any2fasta {input.genome:q} > {output:q}

    """  


# Write the df table to the directory for later reference.
# Why isn't this a run: instead of a shell: ?
rule metadata:
    input: expand("{results_directory}/samples/{sample}/{sample}.fa", results_directory = results_directory, sample = df["sample"]) # From rule copy
    output: "{results_directory}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t")
    resources:
        runtime = "1h"
    #run:
        #df.to_csv(str(output), index_label = "index", sep = "\t")
        #os.system(f"cp ${{ASSCOM2_BASE}}/scripts/{report_template_file_basename} {results_directory}")
    shell: """

        echo '''{params.dataframe}''' > {output}

        {void_report}
    """



# --- Downloads -----------------------------------------------------

# This rule runs once, downloading the busco dataset that is needed for rule busco_individual.
# Make sure that this job is run on a node that has internet access.
rule busco_download:
    output:
        #touch("{base_variable}/databases/busco/busco_download_done.flag") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
        #database_representative = touch("{base_variable}/databases/busco/file_versions.tsv") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
        database_representative = touch(DATABASES + "/busco/ac2_busco_database_representative.flag") # Should point to the directory where the following files reside: "file_versions.tsv  lineages/  placement_files/"
    conda: "conda_definitions/busco.yaml"
    shell: """

        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
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
                --download_path $(dirname {output}) \
                --download prokaryota
            
            
            touch {output}

            # Clean up 
            rm dummy.fasta
        
        fi

    """



# Updated according to chklovski's idea in https://github.com/chklovski/CheckM2/issues/73#issuecomment-1744207103
rule checkm2_download:
    output:
        database_representative = DATABASES + "/checkm2/ac2_checkm2_database_representative.flag",
    params:
        destination = DATABASES + "/checkm2"
    conda: "conda_definitions/curl.yaml"
    shell: """


        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            url="https://zenodo.org/records/5571251/files/checkm2_database.tar.gz"

            
            curl \
                $url \
                --output "{params.destination}/checkm2_database.tar.gz"
                

            tar \
                -xvf "{params.destination}/checkm2_database.tar.gz" \
                --directory "{params.destination}"

            rm "{params.destination}/checkm2_database.tar.gz" || echo "failed to clean up"
            
            touch {output}
        
        fi

    """








rule kraken2_download:
    output:
        #database_representative = touch("{base_variable}/databases/kraken2/ac2_kraken2_database_representative.flag"),
        database_representative = touch(DATABASES + "/kraken2/ac2_kraken2_database_representative.flag"),
    params:
        db_destination_smk = DATABASES + "/kraken2/kraken2_db.tar.gz"
    conda: "conda_definitions/curl.yaml"
    shell: """

        ## Pick a db from this list
        # https://benlangmead.github.io/aws-indexes/k2

        ## Shortcuts. Select no bigger than the size of your RAM
        db_pick="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz"      # Standard 49GB
        #db_pick="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz" # Standard  8GB
        #db_pick="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20230314.tar.gz" # Standard 16GB
                

        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            >&2 echo "Downlading $db_pick to {params.db_destination_smk}"
            #mkdir -p $(dirname "$db_destination")
            curl "$db_pick" \
                --output {params.db_destination_smk}

            >&2 echo "Decompressing ..."
            tar \
                -xvf {params.db_destination_smk} \
                --directory $(dirname {params.db_destination_smk})

            rm {params.db_destination_smk} || echo "Failed to clean up."

            >&2 echo "kraken2 DB setup completed"
            echo "Downloaded $db_pick at $(date -Iseconds)" > $(dirname {params.db_destination_smk})/info.txt

            mkdir -p $(dirname {output})
            touch {output}
            
        fi

    """





rule dbcan_download:
    output:
        #database_representative = touch("{base_variable}/databases/dbcan/ac2_dbcan_database_representative.flag"),
        database_representative = DATABASES + "/dbcan/ac2_dbcan_database_representative.flag",
    conda: "conda_definitions/dbcan.yaml"
    shell: """
        
        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."            

            cd $(dirname {output.database_representative}) \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08252022.tsv \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa && makeblastdb -in PUL.faa -dbtype prot \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.xlsx \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.txt \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz && tar xvf dbCAN-PUL.tar.gz \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm && hmmpress dbCAN_sub.hmm \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/V11/CAZyDB.08062022.fa && diamond makedb --in CAZyDB.08062022.fa -d CAZy \
                && wget --continue https://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt && mv dbCAN-HMMdb-V11.txt dbCAN.txt && hmmpress dbCAN.txt \
                && wget --continue https://bcb.unl.edu/dbCAN2/download/Databases/V11/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-1.hmm && hmmpress tf-1.hmm \
                && wget --continue http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-2.hmm && hmmpress tf-2.hmm \
                && wget --continue https://bcb.unl.edu/dbCAN2/download/Databases/V11/stp.hmm && hmmpress stp.hmm

            # Comments on using the download code from https://github.com/linnabrown/run_dbcan (june 2023): I deleted the test ecoli files in the bottom, and added --continue, to make sure that not a .1 suffixed file is left over when retrying downloads.



            >&2 echo "dbcan setup completed"
            echo "Downloaded dbcan at $(date -Iseconds)" > $(dirname {output.database_representative})/info.txt

            mkdir -p $(dirname {output})
            touch {output}

        fi

    """



# This rule is currently not in use since kegg_diamond is both better and faster.
rule kofam_download:
    output:
        database_representative = DATABASES + "/kofam/ac2_kofam_database_representative.flag",
    conda: "conda_definitions/curl.yaml" # Now it also has unzip which necessary here.
    shell: """


        # TODO: Make a way to check for internet access instead of just crashing. Same for the others


        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output.database_representative}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output.database_representative}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            # Enter the directory so we don't have to pass a lot of the same paths in this script.
            cd $(dirname {output.database_representative})


            # 1) Clone git repository
            wget \
                --continue \
                https://github.com/takaram/kofam_scan/archive/refs/heads/master.zip 

            unzip master.zip \
                -d github-takaram-kofam_scan

            echo "Downloaded takaram/kofam_scan from github at $(date)" > github-takaram-kofam_scan/ac2_info.txt


    
            # 2) Download FTP stuff

            wget --continue ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
            gunzip ko_list.gz

            wget --continue ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
            tar -xf profiles.tar.gz


            rm profiles.tar.gz || echo "Failed to clean up." 


            >&2 echo "kofam setup completed"
            echo "Downloaded kofam at $(date -Iseconds)" > $(dirname {output.database_representative})/info.txt

            mkdir -p $(dirname {output.database_representative})
            touch {output.database_representative}

        fi

    """



rule gtdb_download:
    output:
        #database_representative = touch("{base_variable}/databases/gtdb/gtdb_download_done.flag") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
        #database_representative = touch("{base_variable}/databases/gtdb/ac2_gtdb_database_representative.flag")
        database_representative = DATABASES + "/gtdb/ac2_gtdb_database_representative.flag"
    conda: "conda_definitions/curl.yaml"
    shell: """

        # TODO: Make a way to check for internet access instead of just crashing. Same for busco and checkm2

        # https://ecogenomics.github.io/GTDBTk/installing/index.html

        # Pick a source file
        db_pick="https://ns9864k.web.sigma2.no/TheMEMOgroup/cmkobel/asscom2-assets/gtdb/release214/gtdbtk_data.tar.gz" # NMBU/MEMO mirror in norway
        #db_pick="https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz" # Official location
        #db_pick="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz" # Official alternative mirror


        db_destination=$(dirname {output.database_representative})/gtdb_db.tar.gz

        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            >&2 echo "Downlading $db_pick to $db_destination"
            mkdir -p $(dirname "$db_destination")
            curl "$db_pick" \
                --output "$db_destination"

            >&2 echo "Decompressing ..."
            tar \
                -xf $db_destination \
                --directory $(dirname $db_destination)

            rm $db_destination || echo "Failed to clean up." 

            >&2 echo "gtdb DB setup completed"
            echo "Downloaded $db_pick at $(date -Iseconds)" > $(dirname $db_destination)/info.txt

            mkdir -p $(dirname {output})
            touch {output}

        fi

    """


# ---- Other -------------------------------------------------------------------







# --- Targets for each sample below: --------------------------------

rule sequence_lengths:
    input:
        metadata = "{results_directory}/metadata.tsv",
        assembly = "{results_directory}/samples/{sample}/{sample}.fa", 
    output: "{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"
    threads: 1
    resources:
        runtime = "60m",
    conda: "conda_definitions/seqkit.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.sequence_lengths_individual.{sample}.tsv"
    shell: """

        seqkit fx2tab {input.assembly} -l -g -G -n -H \
        > {output}

    """






rule prokka:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        assembly = "{results_directory}/samples/{sample}/{sample}.fa"
    output:
        gff = "{results_directory}/samples/{sample}/prokka/{sample}.gff",
        faa = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # Used in dbcan, interproscan, diamond_kegg, motupan
        log = "{results_directory}/samples/{sample}/prokka/{sample}.log",
        tsv = "{results_directory}/samples/{sample}/prokka/{sample}.tsv",
        gff_nofasta = "{results_directory}/samples/{sample}/prokka/{sample}.gff_nofasta",
    conda: "conda_definitions/prokka.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.prokka_individual.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
      
        # For debugging of minced
        minced --version 
        
        prokka \
            --cpus {threads} \
            --force \
            --rfam \
            --compliant \
            --outdir {wildcards.results_directory}/samples/{wildcards.sample}/prokka \
            --prefix {wildcards.sample} {input.assembly} \
        | tee {output.log} 

        # I don't remember what I'm actually using this output for?
        # Remove fasta from gff and add sample label
        gff_fasta_start=$(grep --line-number --extended-regexp "^##FASTA" {output.gff} | cut -f1 -d:)
        head --lines $(($gff_fasta_start-1)) {output.gff} \
        > {output.gff_nofasta}

        {void_report}

    """



# I see that much time is spent just reading the database from disk. It is much more efficient to just read the database once, and then iterate through each sample. This can be accomplished be catting all files together, keeping track of the sample names.
# rule kraken2:
#     input: 
#         metadata = "{results_directory}/metadata.tsv",
#         fasta = df["input_file_fasta"].tolist(),
#         database = expand("{base_variable}/databases/kraken2/ac2_kraken2_database_representative.flag", base_variable = base_variable),
#     output: 
#         report = "{results_directory}/kraken2/kraken2_report.tsv",
#         full = "{results_directory}/kraken2/kraken2_full.tsv",
#     params: 
#         #asscom2_kraken2_db = config["asscom2_kraken2_db"],
#         temporary_concatenation = "{results_directory}/kraken2/temporary_concatenation.fa",
#         base_variable = base_variable,
#     conda: "conda_definitions/kraken2.yaml"
#     container: "docker://cmkobel/kraken2" # was disabled already
#     benchmark: "{results_directory}/benchmarks/benchmark.kraken2_all.tsv"
#     threads: 8
#     resources:
#         mem_mb = 64000,
#     shell: """

#         db_path="{params.base_variable}/databases/kraken2"
#         echo using kraken2 database $db_path

#         # Purge old concatenation file
#         test -f {params.temporary_concatenation} && rm {params.temporary_concatenation}

#         # Concatenate assemblies with an &&&& separator in the header
#         for fasta in {input.fasta:q}; do
#             bn=$(basename $fasta)
#             echo "Catting ${{bn}} ..."

#             sed "s/>.*/&\&\&\&\&\"$bn\"/" $fasta \
#             >> {params.temporary_concatenation}

#         done


#         # Run kraken2
#         # https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
#         kraken2 \
#             --threads {threads} \
#             --db $db_path \
#             --report {output.report} \
#             --use-mpa-style \
#             --use-names \
#             {params.temporary_concatenation} \
#             > {output.full}

#         # Argument on confidence parameter https://www.biostars.org/p/402619/ --confidence 0.1 \


#         # Clean up
#         rm {params.temporary_concatenation}

#         {void_report}
#     """

rule kraken2:
    input: 
        metadata = expand("{results_directory}/metadata.tsv", results_directory = results_directory),
        assembly = "{results_directory}/samples/{sample}/{sample}.fa",
        #database = expand("{base_variable}/databases/kraken2/hash.k2d", base_variable = base_variable),
        #database = expand("{base_variable}/databases/kraken2/ac2_kraken2_database_representative.flag", base_variable = base_variable),
        database_representative = DATABASES + "/kraken2/ac2_kraken2_database_representative.flag"
    output: 
        report = "{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv",
        full = "{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_full.tsv",
    conda: "conda_definitions/kraken2.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.kraken2_individual.{sample}.tsv"
    threads: 2
    resources:
        mem_mb = 75000,
    shell: """

        echo using kraken2 database $(dirname {input.database_representative})

        # Run kraken2
        # https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
        kraken2 \
            --threads {threads} \
            --db $(dirname {input.database_representative}) \
            --confidence 0.1 \
            --report {output.report} \
            --report-minimizer-data \
            {input.assembly} \
            > {output.full}

        # Argument on confidence parameter https://www.biostars.org/p/402619/

    """

rule dbcan: # I can't decide whether this rule should really be called "run_dbcan", since that is the name of the software.
    input: 
        metadata = "{results_directory}/metadata.tsv",
        aminoacid = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # From prokka
        #database_representative = expand("{base_variable}/databases/dbcan/ac2_dbcan_database_representative.flag", base_variable = base_variable),
        database_representative = DATABASES + "/dbcan/ac2_dbcan_database_representative.flag"
    output: 
        overview_table = "{results_directory}/samples/{sample}/dbcan/overview.txt"
    params: 
        out_dir = "{results_directory}/samples/{sample}/dbcan"
    conda: "conda_definitions/dbcan.yaml" # Not sure if it should be called by a version number?
    # container: TODO # was disabled already
    benchmark: "{results_directory}/benchmarks/benchmark.dbcan.{sample}.tsv"
    threads: 8
    resources: 
        mem_mb = 8000
    shell: """

        # It seems to be necessary to set all the cpu thread counts manually.

        export HMMER_NCPU={threads}

        run_dbcan \
            --dbcan_thread {threads} \
            --dia_cpu {threads} \
            --hmm_cpu {threads} \
            --tf_cpu {threads} \
            --stp_cpu {threads} \
            --db_dir $(dirname {input.database_representative}) \
            --out_dir {params.out_dir} \
            {input.aminoacid} \
            protein 

    """

# By default it runs TIGRFAM, Hamap, Pfam
rule interproscan:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        aminoacid = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # From prokka
        # database_representative # No external database is needed.
    output:
        tsv = "{results_directory}/samples/{sample}/interproscan/{sample}_interproscan.tsv",
    params:
        file_base = "{results_directory}/samples/{sample}/interproscan/{sample}_interproscan",
    conda: "conda_definitions/interproscan.yaml" # Not sure if it should be called by a version number?
    benchmark: "{results_directory}/benchmarks/benchmark.interproscan.{sample}.tsv"
    threads: 8
    resources: 
        mem_mb = 8000
    shell: """

        # https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html#command-line-options
        interproscan.sh \
            --applications TIGRFAM,Hamap,Pfam \
            --cpu {threads} \
            --output-file-base {params.file_base} \
            --disable-precalc \
            --formats TSV \
            --goterms \
            --iprlookup \
            --pathways \
            --seqtype p \
            --tempdir {resources.tmpdir} \
            --input {input.aminoacid}

    """



# Disabled because I think the diamond_kegg based results are much better, at least faster.
rule kofam_scan:
    input: 
        metadata = "{results_directory}/metadata.tsv", # For the report
        aminoacid = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # From prokka
        database_representative = base_variable + "/databases/kofam/ac2_kofam_database_representative.flag",
    output:
        full = "{results_directory}/samples/{sample}/kofam_scan/{sample}_kofam_scan_full.tsv",
        significant = "{results_directory}/samples/{sample}/kofam_scan/{sample}_kofam_scan_significant.tsv",
    params: 
        executable = base_variable + "/databases/kofam/github-takaram-kofam_scan/kofam_scan-master/exec_annotation",
        ko_list = base_variable + "/databases/kofam/ko_list",
        profile = base_variable + "/databases/kofam/profiles"
    conda: "conda_definitions/kofam_scan.yaml"
    # container: TODO # was disabled already
    benchmark: "{results_directory}/benchmarks/benchmark.kofam_scan.{sample}.tsv"
    threads: 4
    shell: """

        # Must use unique --tmp-dir because of https://github.com/takaram/kofam_scan/issues/8#issuecomment-619640654

        # https://github.com/takaram/kofam_scan/
        {params.executable} \
            --cpu {threads} \
            --ko-list {params.ko_list} \
            --profile {params.profile} \
            --format detail-tsv \
            --e-value 0.0001 \
            --tmp-dir $(dirname {output.full})/tmp \
            -o {output.full} \
            {input.aminoacid}


        # kofam_scan outputs way, way too many spurious results, so here I'm making a version that is filtered for "singificant" results only.
        echo "Making a filtered version with only the asterisk (*) marked alignments."
        grep \
            -E "^\*" \
            {output.full} \
            > {output.significant}
        # And then, one could consider deleting the "full" file altogether, and just stick with the significant results.


        # Clean up temporary directory
        rm -r $(dirname {output.full})/tmp



    """


rule busco:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        #busco_download = expand("{base_variable}/databases/busco/file_versions.tsv", base_variable = base_variable), # This is a bad idea, because it requires a complete reinstall if snakemake somehow removes the file, which is quite likely.
        database_representative = DATABASES + "/busco/ac2_busco_database_representative.flag", # Should point to the directory where the following files reside: "file_versions.tsv  lineages/  placement_files/"
        fasta = "{results_directory}/samples/{sample}/{sample}.fa",
    output: 
        flag = touch("{results_directory}/samples/{sample}/busco/busco_done.flag"),
        table_extract = "{results_directory}/samples/{sample}/busco/short_summary_extract.tsv"
    params:
        base_variable = base_variable,
        #results_directory = results_directory,
        database_path = DATABASES + "/busco", # Was {params.base_variable}/databases/busco
        out_dir = "{results_directory}/samples/{sample}/busco",
    conda: "conda_definitions/busco.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.busco_individual.{sample}.tsv"
    threads: 1 # Because run_sepp hangs for a long time, not doing anything, I'd rather have more processes started on my small CPU.
    resources:
        mem_mb = 8192,
        runtime = "6h",
    shell: """

        # Busco fails because of a problem with the sepp package. This doesn't really matter as we just want the completeness results.
        # But, this means that we need a hacky workaround to let this job exit gracefully (exit code 0) on the basis of whether any completeness results have been written to disk.
        # Hence, the actual exit code of busco, we will ignore.


        >&2 echo "Busco individual"
        # https://busco.ezlab.org/busco_userguide.html#offline
        timeout 3600 \
            busco \
                --cpu {threads} \
                --in {input.fasta} \
                --out {params.out_dir} \
                --mode geno \
                --auto-lineage-prok \
                --force \
                --tar \
                --download_path {params.database_path} \
                --offline || (>&2 echo "ac2: busco failed internally")


        >&2 echo "debug1"
        ### New: Using JSON


        # Cat all auto lineage results together or create empty file
        # The following cat command will fail if the glob doesn't resolve any files: This is the wanted behaviour.
        cat {wildcards.results_directory}/samples/{wildcards.sample}/busco/auto_lineage/*/short_summary.json \
        > {output.table_extract}_temp         
        
        >&2 echo "Results clearly must exist ... "
        
        # Extract relevant features
        cat {output.table_extract}_temp \
        | grep -oE "(\\"in\\"|\\"name\\"|\\"one_line_summary\\").+" \
        > {output.table_extract}

        >&2 echo "debug3"
        # Clean up
        rm {output.table_extract}_temp

        >&2 echo "debug4"

        {void_report}

    """


# --- Targets for the complete set below: ---------------------------

rule checkm2:
    input:
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/checkm2/ac2_checkm2_database_representative.flag",
        fasta = df["input_file_fasta"].tolist()
    output:
        table = touch("{results_directory}/checkm2/quality_report.tsv"),
    conda: "conda_definitions/checkm2.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.checkm2.tsv"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = "24h",
    params:
        rule_dir = results_directory + "/checkm2",
        base_variable = base_variable,
    shell: """

        checkm2 predict \
            --threads {threads} \
            --input {input.fasta} \
            --output-directory {params.rule_dir} \
            --extension .fa \
            --database_path $(dirname {input.database_representative})/CheckM2_database/uniref100.KO.1.dmnd \
            --force

        {void_report}

    """



# Use the same database as checkm2, but run on amino-acid files instead of dna.
# This will be used for a subsequent pathway enrichment analysis.
# Idea for speed up: concatenate all genomes together first, like they do in checkm2. Then we only need to load the database once.
rule diamond_kegg: # or uniref_ko?
    input: 
        metadata = "{results_directory}/metadata.tsv", # For the report
        aminoacid = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # From prokka
        #database_representative = base_variable + "/databases/checkm2/ac2_checkm2_database_representative.flag",
        database_representative = DATABASES + "/checkm2/ac2_checkm2_database_representative.flag",
        
    output:
        tsv = "{results_directory}/samples/{sample}/diamond_kegg/{sample}_diamond_kegg.tsv", # Or should it be named uniref-100?
    params: 
        query_cover = 85,
        subject_cover = 85, # 
        percent_id = 50, # 30 is probably fine for checkm2, but I feel like I'd rather have censored data than spurious results.
        evalue = "1e-05",
        blocksize = 2, # A value of 2 corresponds to running checkm2 in non-lowmem mode.
        database_path = DATABASES + "/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
    conda: "conda_definitions/diamond.yaml" 
    resources:
        mem_mb = 20000, # Seems to use around 18G at max.
        runtime = "1h",
    # container: TODO # was disabled already
    benchmark: "{results_directory}/benchmarks/benchmark.diamond_kegg.{sample}.tsv"
    threads: 8
    shell: """

        # Inspired from https://github.com/chklovski/CheckM2/blob/319dae65f1c7f2fc1c0bb160d90ac3ba64ed9457/checkm2/diamond.py#L79
    
        # blastp: Align amino acid query sequences against a protein reference database

        diamond blastp \
            --outfmt 6 \
            --max-target-seqs 1 \
            --query {input.aminoacid}  \
            --out {output.tsv}  \
            --threads {threads}  \
            --db {params.database_path} \
            --query-cover {params.query_cover}  \
            --subject-cover {params.subject_cover}  \
            --id {params.percent_id}  \
            --evalue {params.evalue} \
            --block-size {params.blocksize}

    """



# Runs per batch, should maybe be moved in this document.
rule kegg_pathway:
    input: 
        kegg_asset = base_variable + "/assets/ko00001.json", # Downloaded from kegg.jp
        #diamond = rules.diamond_kegg.output,
        kegg_diamond = expand("{results_directory}/samples/{sample}/diamond_kegg/{sample}_diamond_kegg.tsv", sample = df["sample"], results_directory = results_directory)
    output: 
        #"{results_directory}/samples/{sample}/kegg_pathway/{sample}_kegg_pathway.tsv",
        diamond = "{results_directory}/kegg_pathway/kegg_pathway_enrichment_analysis.tsv"
    params:
        output_dir = "{results_directory}/kegg_pathway",
        script = base_variable + "/scripts/kegg_pathway_enrichment_analysis.R"
    conda: "conda_definitions/r-clusterProfiler.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.kegg_pathway.tsv"
    shell: """

        Rscript {params.script} \
            {input.kegg_asset} \
            {params.output_dir} \
            {input.kegg_diamond}  
            
        {void_report}

    """





def get_mem_roary(wildcards, attempt): 
    return [32000, 64000, 128000][attempt-1]


rule roary:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        gff = expand("{results_directory}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], results_directory = results_directory),
    output:
        analyses = ["{results_directory}/roary/summary_statistics.txt", "{results_directory}/roary/core_gene_alignment.aln", "{results_directory}/roary/gene_presence_absence.csv"]
    params:
        blastp_identity = int(config['roary_blastp_identity']), # = 95 # For clustering genes
        core_perc = 99,  # Definition of the core genome
    benchmark: "{results_directory}/benchmarks/benchmark.roary.tsv"
    threads: 16
    #retries: 2
    resources:
        #mem_mb = 32768,
        mem_mb = get_mem_roary,
        runtime = "24h",
    conda: "conda_definitions/roary_see-comments-in-this-file.yaml"
    shell: """
    
        # Since I reinstalled conda, I've had problems with "Can't locate Bio/Roary/CommandLine/Roary.pm in INC". Below is a hacky fix
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0
        

        # Silence parallel's citation pester:
        echo "will cite" | parallel --citation > /dev/null 2> /dev/null

        # Roary is confused by the way snakemake creates directories ahead of time.
        rm -r {wildcards.results_directory}/roary


        roary -a -r -e --mafft \
            -p {threads} \
            -i {params.blastp_identity} \
            -cd {params.core_perc} \
            -f {wildcards.results_directory}/roary \
            --group_limit 100000 \
            {input.gff}

        # Default group limit is 50000
            
        {void_report}

    """





rule snp_dists:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        aln = "{results_directory}/roary/core_gene_alignment.aln",
    output: "{results_directory}/snp-dists/snp-dists.tsv"
    conda: "conda_definitions/snp-dists.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.snp_dists.tsv"
    threads: 4
    shell: """

        snp-dists \
            -j {threads} \
            {input.aln} > {output}

        {void_report}
        
    """


rule motulizer:
    input:
        fnas = df["input_file_fasta"].tolist(),
    output: 
        tsv = "{results_directory}/motulizer/motulizer_results.tsv"
    params:
        version_file = "{results_directory}/motulizer/motulizer_ac2_version.txt"
    conda: "conda_definitions/motulizer.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.motulizer.tsv"
    threads: 1
    shell: """
        
        mOTUlize.py --version > {params.version_file} || echo "Catched exit code 1 when asking for the motulize version."

        mOTUlize.py --fnas {input.fnas} -o {output.tsv} 

    """

# motulizer and motupan could run together in the same rule, but I like how the first job can start right away and the other can run trailing prokka. Also (todo), in the future I might want to run motupan separately for each motulizer-mOTU.

rule motupan: 
    input: 
        faas =  expand("{results_directory}/samples/{sample}/prokka/{sample}.faa", results_directory = results_directory, sample = df["sample"]) # From prokka
    output: 
        tsv = "{results_directory}/motulizer/motupan_results.tsv"
    conda: "conda_definitions/motulizer.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.motupan.tsv"
    threads: 1
    shell: """
    
        # Same version as motulizer, so no need to save the version again.

        mOTUpan.py --faas {input.faas} -o {output.tsv} 

    """



rule assembly_stats:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{results_directory}/assembly-stats/assembly-stats.tsv"
    conda: "conda_definitions/assembly-stats.yaml"
    benchmark: "{results_directory}/benchmarks/assembly_stats.tsv"
    shell: """
        
        assembly-stats -t {input.fasta} > {output}

        {void_report}
    """


def get_mem_gtdbtk(wildcards, attempt): 
    return [150000, 300000, 400000, 500000][attempt-1]



# The rule is called by the software, whereas the results are called by the database. Is that confusing?
rule gtdbtk:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        #database_representative = expand("{base_variable}/databases/gtdb/ac2_gtdb_database_representative.flag", base_variable = base_variable),
        database_representative = DATABASES + "/gtdb/ac2_gtdb_database_representative.flag",
        
        fasta = df["input_file_fasta"].tolist(),
    output: "{results_directory}/gtdbtk/gtdbtk.summary.tsv"
    params:
        batchfile_content = df[['input_file_fasta', 'sample']].to_csv(header = False, index = False, sep = "\t"),
        out_dir = "{results_directory}/gtdbtk/",
        base_variable = base_variable,
        mash_db = f"{DATABASES}/gtdb_sketch/mash_db.msh",
        #gtdbtk_data_path = config["gtdbtk_data_path"],
    threads: 8
    #retries: 3
    resources:
        #mem_mb = 150000, # Last time I remember, it used 130000
        # mem_mb = get_mem_gtdbtk, # works great, but I want to use a limited amount to test a potential hardware upgrade
        mem_mb = 131072,
        runtime = "48h"
    conda: "conda_definitions/gtdbtk.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.gtdbtk.tsv"
    shell: """

        # TODO: Using skip-ani-screen is not optimal, as it possibly speeds up a lot.
        mkdir -p $(dirname {params.mash_db})

        # I need to find a neat way of setting these variables. Maybe the user has an older/newer version than what is hardcoded here. 
        export GTDBTK_DATA_PATH="$(dirname {input.database_representative})/release214/" # Should be defined from config file, and not be hardwired.
        

        # Create batchfile
        echo '''{params.batchfile_content}''' > {wildcards.results_directory}/gtdbtk/batchfile.tsv
        
        gtdbtk classify_wf \
            --mash_db {params.mash_db} \
            --batchfile {wildcards.results_directory}/gtdbtk/batchfile.tsv \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --keep_intermediates \
            --force

        # Homogenize database version number
        #cp {wildcards.results_directory}/gtdbtk/gtdbtk.bac120.summary.tsv {output}
        # New better version below that also incorporates archaea
        #cat {wildcards.results_directory}/gtdbtk/gtdbtk.*.summary.tsv {output}
        

        # Even better: Should be tested on originals
        echo -e "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings" \
        > {output}
        tail --quiet -n +2 {wildcards.results_directory}/gtdbtk/gtdbtk.*.summary.tsv \
        >> {output}
        


        {void_report}


    """ 



rule abricate:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output:
        ncbi_detailed = "{results_directory}/abricate/ncbi_detailed.tsv",
        ncbi_sum = "{results_directory}/abricate/ncbi_summarized.tsv",
        card_detailed = "{results_directory}/abricate/card_detailed.tsv",
        card_sum = "{results_directory}/abricate/card_summarized.tsv",
        plasmidfinder_detailed = "{results_directory}/abricate/plasmidfinder_detailed.tsv",
        plasmidfinder_sum = "{results_directory}/abricate/plasmidfinder_summarized.tsv",
        vfdb_detailed = "{results_directory}/abricate/vfdb_detailed.tsv",
        vfdb_sum = "{results_directory}/abricate/vfdb_summarized.tsv",
    conda: "conda_definitions/abricate.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.abricate.tsv"
    shell: """

        # TODO: update these databases

        abricate --db ncbi {input.fasta} > {output.ncbi_detailed}
        abricate --summary {output.ncbi_detailed} > {output.ncbi_sum}

        abricate --db card {input.fasta} > {output.card_detailed}
        abricate --summary {output.card_detailed} > {output.card_sum}
        
        abricate --db plasmidfinder {input.fasta} > {output.plasmidfinder_detailed}
        abricate --summary {output.plasmidfinder_detailed} > {output.plasmidfinder_sum}

        abricate --db vfdb {input.fasta} > {output.vfdb_detailed}
        abricate --summary {output.vfdb_detailed} > {output.vfdb_sum}

        {void_report}
    """



# Parse the mlst scheme for bash
if config["mlst_scheme"] == "automatic":
    mlst_scheme_interpreted = "",
else:
    mlst_scheme_interpreted = f"--scheme {config['mlst_scheme']}",
#print(f"Info: The mlst_scheme is set to <{mlst_scheme_interpreted}>") # Debug message.

rule mlst:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{results_directory}/mlst/mlst.tsv",
    params:
        mlst_scheme_interpreted = mlst_scheme_interpreted,
        list_ = "{results_directory}/mlst/mlst_schemes.txt", 
    threads: 4
    conda: "conda_definitions/mlst.yaml"
    benchmark: "{results_directory}/benchmarks/mlst.tsv"
    shell: """

        mlst \
            --threads {threads} {params.mlst_scheme_interpreted} \
            {input.fasta} \
            > {output}

        # Dump available mlst databases
        mlst --list > {params.list_}

        {void_report}
    """




rule mashtree:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: 
        tree = "{results_directory}/mashtree/mashtree.newick",
        dist = "{results_directory}/mashtree/mash_dist.tsv",
    conda: "conda_definitions/mashtree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.mashtree.tsv"
    threads: 16
    resources:
        mem_mb = 16000,
    shell: """

        mashtree \
            --numcpus {threads} \
            --outmatrix {output.dist} \
            {input.fasta} > {output.tree}

        {void_report}
    """ 

# TODO:
#rule mash_screen:

def get_mem_fasttree(wildcards, attempt): 
    return [16000, 32000, 64000, 0][attempt-1]


rule fasttree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = "{results_directory}/roary/core_gene_alignment.aln",
    output: "{results_directory}/fasttree/fasttree.newick"
    conda: "conda_definitions/fasttree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.fasttree.tsv"
    threads: 4
    retries: 2
    resources:
        mem_mb = get_mem_fasttree,
        runtime = "24h",
    shell: """

        OMP_NUM_THREADS={threads}

        FastTree \
            -nt \
            -gtr {input.fasta} \
        > {output} \
        2> {output}.log 

        {void_report}

    """



def get_mem_iqtree(wildcards, attempt): 
    return [16000, 32000, 64000, 128000][attempt-1]

rule iqtree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = "{results_directory}/roary/core_gene_alignment.aln",
    output: 
        newick = "{results_directory}/iqtree/core_genome_iqtree.treefile"
    params:
        version_file = "{results_directory}/iqtree/iqtree_ac2_version.txt"
    conda: "conda_definitions/iqtree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.iqtree.tsv"
    threads: 16
    retries: 3
    resources:
        mem_mb = get_mem_iqtree,
        runtime = "24h",
    shell: """

        iqtree --version > {params.version_file}

        iqtree \
            -s {input.fasta} \
            -m GTR \
            --boot 100 \
            --prefix $(dirname {output.newick})/core_genome_iqtree \
            -redo

        # {void_report} Not in the report yet.


    """


# rule fetch_report_template:
#     output: "{results_directory}/rmarkdown_template.rmd"
#     shell: """

#         cp $ASSCOM2_BASE/scripts/genomes_to_report_v2.Rmd {output}

#         {void_report}
#     """



# This rule might seem silly, but it makes sure that the report environment is ready to rock when the report subpipeline eventually is run: This has two pros:
#    1) The vastly faster mamba configuration in the asscom2 pipeline is used
#    2) The conda/mamba debugging is taken care of, without having to wait for jobs to finish on fresh installations.
# Since all snakemake conda environments are installed in $SNAKEMAKE_CONDA_PREFIX set to ${ASSCOM2_BASE}/conda_base, reuse is guaranteed.
rule install_report_environment_aot:
    output: touch(f"{results_directory}/.install_report_environment_aot.flag")
    conda: "report_subpipeline/conda_definitions/r-markdown.yaml"
    shell: """

        echo "Report conda environment OK ..."

    """

# Just a dummy rule if you wanna force the report
# assemblycomparator2 --until report
# It isn't enough to just touch the file. The report_subpipeline will not be triggered if the file is empty. Thus we add the date, and we have a nice debug log for seeing when the report was triggered.
# Will only but run if asked to. No need to use --forcerun, since snakemake states this in the output: "reason: Rules with neither input nor output files are always executed."
# Rule report does not depend on metadata, as the metadata is not interesting in itself.
rule report:
    shell: """
        
        {void_report}

    """



## Pseudo targets:

# Makes it easy to check that all databases are installed properly. Eventually for touching the database representatives in case of using prior installations.
rule downloads:
    input:
        DATABASES + "/checkm2/ac2_checkm2_database_representative.flag", 
        DATABASES + "/kraken2/ac2_kraken2_database_representative.flag",
        DATABASES + "/busco/ac2_busco_database_representative.flag", 
        DATABASES + "/dbcan/ac2_dbcan_database_representative.flag", 
        DATABASES + "/gtdb/ac2_gtdb_database_representative.flag"


# Blink-of-an-eye analysis
rule fast:
    input:
        expand(\
            ["{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", \
            "{results_directory}/assembly-stats/assembly-stats.tsv", \
            "{results_directory}/mashtree/mashtree.newick"], \
            results_directory = results_directory, \
            sample = df["sample"]) # TODO: define the expansion in each rule instead.


# Rules designed for bins of metagenomic origin
rule meta: 
    input:
        expand(\
            ["{results_directory}/metadata.tsv", \
            "{results_directory}/.install_report_environment_aot.flag", \
            "{results_directory}/assembly-stats/assembly-stats.tsv", \
            "{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", \
            "{results_directory}/samples/{sample}/busco/short_summary_extract.tsv", \
            "{results_directory}/checkm2/quality_report.tsv", \
            "{results_directory}/samples/{sample}/diamond_kegg/{sample}_diamond_kegg.tsv", \
            "{results_directory}/kegg_pathway/kegg_pathway_enrichment_analysis.tsv", \
            "{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv", \
            "{results_directory}/samples/{sample}/dbcan/overview.txt", \
            "{results_directory}/samples/{sample}/interproscan/{sample}_interproscan.tsv", \
            "{results_directory}/gtdbtk/gtdbtk.summary.tsv", \
            "{results_directory}/mlst/mlst.tsv", \
            "{results_directory}/samples/{sample}/prokka/{sample}.gff", \
            "{results_directory}/mashtree/mashtree.newick"], \
            results_directory = results_directory, \
            sample = df["sample"])


# Rules designed for cultured isolates
rule isolate:
    input: expand(\
        ["{results_directory}/metadata.tsv", \
        "{results_directory}/.install_report_environment_aot.flag", \
        "{results_directory}/assembly-stats/assembly-stats.tsv", \
        "{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", \
        "{results_directory}/samples/{sample}/diamond_kegg/{sample}_diamond_kegg.tsv", \
        "{results_directory}/kegg_pathway/kegg_pathway_enrichment_analysis.tsv", \
        "{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv", \
        "{results_directory}/gtdbtk/gtdbtk.summary.tsv", \
        "{results_directory}/mlst/mlst.tsv", \
        "{results_directory}/abricate/card_detailed.tsv", \
        "{results_directory}/samples/{sample}/prokka/{sample}.gff", \
        "{results_directory}/roary/summary_statistics.txt", \
        "{results_directory}/fasttree/fasttree.newick", \
        "{results_directory}/snp-dists/snp-dists.tsv", \
        "{results_directory}/mashtree/mashtree.newick"], \
        results_directory = results_directory, sample = df["sample"])




# A major todo is to find a way to make the report run as a conda or containerized job dependending on the use-conda/use-singularity setting in the config. The best way might be to have an environment variable that points to the wanted config. 

# Call the report subpipeline
# I wonder if adding the ampersands means that the report creation will be ^c cancellable? Not tested yet..
report_call = f"""
    
    # Had to chain these tests/if-statements in order to be able to cancel the running interactively.
    (
        test -f "{results_directory}/.asscom2_void_report.flag" \
        && test -f "{results_directory}/metadata.tsv" \
        && mkdir -p {results_directory}/logs \
        && snakemake \
            --snakefile $ASSCOM2_BASE/report_subpipeline/snakefile \
            --profile $ASSCOM2_PROFILE \
            --config results_directory=$(pwd)/{results_directory} base_variable={base_variable} batch_title={batch_title}
    ) \
    || echo "Info: Not calling the report_subpipeline as either the flag or metadata is missing. Run a rule that mandates a report section to generate the report."

    """

onsuccess:
    print("On success: calling report subpipeline ...")
    shell(report_call)

onerror:
    print("On error: calling report subpipeline ...")
    shell(report_call)





print("*/")




