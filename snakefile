

# snakemake --snakefile ~/assemblycomparator2/snakefile --profile ~/assemblycomparator2/configs/slurm/ --cluster-config ~/assemblycomparator2/configs/cluster.yaml 

__version__ = "v2.3.0"
__author__ = 'Oliver Kjærlund Hansen & Carl M. Kobel'

import os
from os import listdir
from os.path import isfile, join
#import yaml
import pandas as pd
import numpy as np
from shutil import copyfile
#import time
#import re
#from shutil import copyfile
#import re
cwd = os.getcwd()
batch_title = cwd.split("/")[-1]
print("/*")
print()
print("         █████╗ ███████╗███████╗ ██████╗ ██████╗ ███╗   ███╗██████╗ ")
print("        ██╔══██╗██╔════╝██╔════╝██╔════╝██╔═══██╗████╗ ████║╚════██╗")
print("        ███████║███████╗███████╗██║     ██║   ██║██╔████╔██║ █████╔╝")
print("        ██╔══██║╚════██║╚════██║██║     ██║   ██║██║╚██╔╝██║██╔═══╝ ")
print("        ██║  ██║███████║███████║╚██████╗╚██████╔╝██║ ╚═╝ ██║███████╗")
print("        ╚═╝  ╚═╝╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝")
print("                       A.K.A. assemblycomparator2                   ")
print("                         Please log issues at:                      ")
print("              github.com/cmkobel/assemblycomparator2/issues         ")
print()
print(f"            batch_title:            {batch_title}")
print(f"            roary_blastp_identity:  {config['roary_blastp_identity']} (default 95)")
print(f"            mlst_scheme:            {config['mlst_scheme']} (default automatic)")
print()



out_base_var = "output_asscom2"

base_variable = os.environ['ASSCOM2_BASE']


print('base_variable:', base_variable)


#reference = config["reference"]



# --- Read in relevant files in the current working directory ----------------

relative_wd = "."
#extension_whitelist = ["fna", "fa", "fas", "fasta", "seq"] # old 
extension_whitelist = ["fna", "fa", "fas", "fasta", "seq", "gb", "fq", "gff", "gfa", "clw", "sth", "gz", "bz2"]

present_files = [f for f in listdir(relative_wd) if isfile(join(relative_wd,f))]


df = pd.DataFrame(data = {'input_file': present_files})

# Check that the directory is not empty.
if df.shape[0] == 0:
    print("Error: No fasta files in the current directory. Quitting ...")
    raise Exception("Zero genomic files present.")



df = df[~df["input_file"].str.startswith(".", na = False)] # Remove hidden files
df['sample_raw'] = [".".join(i.split(".")[:-1]) for i in df['input_file'].tolist()] # Extract everything before the extension dot.
df['sample'] = df['sample_raw'].str.replace(' ','_')
df['extension'] =  [i.split(".")[-1] for i in df['input_file'].tolist()] # Extract extension
df['input_file_fasta'] = out_base_var + "/samples/" + df['sample'] + "/" + df['sample'] + ".fa" # This is where the input file is copied to in the first snakemake rule.

df = df[df['extension'].isin(extension_whitelist)] # Remove files with unsupported formats.

# Check that the directory is not empty, again.
if df.shape[0] == 0:
    print("Error: No fasta files in the current directory. Quitting ...(2)")
    raise Exception("Zero genomic files present.")


#df_mini = df_mini.apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # strip whitespace and replace spaces with underscores.

  
# --- Displaying filtered dataframe ready for analysis --------------

df = df.reset_index(drop = True)
#print(df[['input_file', 'sample', 'extension']])
print(df)
print("//")
print()




# --- Make sure the output directory exists. ---------------------------

try: 
    os.mkdir("output_asscom2")
except:
    pass



# The modification time of this file tells the report subpipeline whether it needs to run. Thus, void_report is called in the end of every successful rule.
#void_report = f"touch {out_base_var}/.asscom2_void_report.flag"
void_report = f"echo $(date --iso-8601=seconds) >> {out_base_var}/.asscom2_void_report.flag"






# --- Collect all targets. ------------------------------------------
rule all:
    input: expand([\
        "{out_base}/metadata.tsv", \
        "{out_base}/.install_report_environment_aot.flag", \
        "{out_base}/checkm2/checkm2.tsv", \
        "{out_base}/assembly-stats/assembly-stats.tsv", \
        "{out_base}/collected_results/sequence_lengths.tsv", \
        "{out_base}/collected_results/GC_summary.tsv", \
        "{out_base}/collected_results/prokka_summarized.txt", \
        "{out_base}/collected_results/kraken2_reports.tsv", \
        "{out_base}/collected_results/busco.tsv", \
        "{out_base}/roary/summary_statistics.txt", \
        "{out_base}/abricate/card_detailed.tsv", \
        "{out_base}/mashtree/mashtree.newick", \
        "{out_base}/mlst/mlst.tsv", \
        "{out_base}/fasttree/fasttree.newick", \
        "{out_base}/gtdbtk/gtdbtk.bac.summary.tsv", \
        "{out_base}/snp-dists/snp-dists.tsv"], \
        out_base = out_base_var) 




# Dummy test
rule test:
    output: ".test_done.flag"
    shell: """
        #sleep 5
        touch {output}

        {void_report}
    """


# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
rule copy:
    #input: "{sample}"
    input: 
        genome = lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0],
    output: "{out_base}/samples/{sample}/{sample}.fa"
    #log: "logs/{out_base}_{wildcards.sample}.out.log"
    container: "docker://pvstodghill/any2fasta"
    conda: "conda_definitions/any2fasta.yaml"
    resources:
        mem_mb = 512,
        runtime = "00:10:00"
    shell: """

        any2fasta "{input.genome}" > {output}

    """  


# Write the df table to the directory for later reference.
# Why isn't this a run: instead of a shell: ?
rule metadata:
    input: expand("{out_base}/samples/{sample}/{sample}.fa", out_base = out_base_var, sample = df["sample"]) # From rule copy
    output: "{out_base}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t")
    resources:
        runtime = "01:00:00"
    #run:
        #df.to_csv(str(output), index_label = "index", sep = "\t")
        #os.system(f"cp ${{ASSCOM2_BASE}}/scripts/{report_template_file_basename} {out_base_var}")
    shell: """

        echo '''{params.dataframe}''' > {output}


        {void_report}
    """


# Mac-only problem (possibly only on M1 which no one uses anyway..)
# Seems that checkm doesn't work on mac. pplacer does not exist, and I get errors:
#   AttributeError: 'MarkerGeneFinder' object has no attribute '__reportProgress'
#   AttributeError: 'MarkerGeneFinder' object has no attribute '__processBin'
#   [2022-09-27 10:35:55] INFO: Saving HMM info to file.
#   [2022-09-27 10:35:55] INFO: Calculating genome statistics for 3 bins with 1 threads:
#   ...
#   [2022-09-27 10:35:56] INFO: Extracting marker genes to align.
#   [2022-09-27 10:35:56] ERROR: Models must be parsed before identifying HMM hits.
#   I will have to consider if I will develop this on linux, or find an alternative.
# rules download_checkm and checkm have been disabled below
# rule checkm_download:
#     output:
#         flag = touch("{out_base}/.checkm_OK.flag")
#     conda: "conda_definitions/wget.yaml"
#     params:
#         directory = base_variable + "/databases/checkm/" # trailing slash??
#     shell: """
 
#     # Check if the database exists. 
#     # If it doesn't, download/untar the db
#     if [ ! -f {params.directory}/checkm_OK.flag ]; then    
 
#         mkdir -p {params.directory}
#         wget --directory-prefix={params.directory} https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz 
#         tar -xvf {params.directory}/checkm_data_2015_01_16.tar.gz -C {params.directory}
#         touch {params.directory}/checkm_OK.flag # This is just to make sure that the process went well.
 
#     fi

#     # If the flag exists already, then the .checkm_OK.flag will be touched immediately by snakemake
 
#     """
#
#rule checkm:
#    input: "{out_base}/.checkm_OK.flag"
#    output: touch("{out_base}/checkm/output")
#    conda: "conda_definitions/checkm.yaml"
#    params:
#            directory = base_variable + "/databases/checkm/"
#    shell: """
#
#        checkm data setRoot {params.directory}
#        
#
#
#        checkm lineage_wf /Users/kartoffel/assemblycomparator2/tests/E._faecium ./output
#
#    """



# --- CheckM2 --------------------------------------------------------

# checkm2 is a bit tricky since it requires a manual installation. If
# it turns out that it is neccessary to have the git download in the 
# same directory as where it is run, then it might not be feasible to 
# use it at all.
rule checkm2:
    input:
        metadata = "{out_base}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist()
    output:
        table = touch("{out_base}/checkm2/checkm2.tsv"),
        diamond = touch("{out_base}/checkm2/diamond_output/DIAMOND_RESULTS.tsv")
    conda: "conda_definitions/checkm2.yaml"
    threads: 1
    resources:
        mem_mb = 16000,
    params:
        rule_dir = out_base_var + "/checkm2",
        base_variable = base_variable
    shell: """

        {params.base_variable}/assets/checkm2 predict \
            --input {input.fasta} \
            --output-directory {params.rule_dir} \
            --extension .fa \
            --force

    """





# --- Targets for each sample below: --------------------------------

rule seqlen_individual:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output: "{out_base}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"
    container: "docker://cmkobel/bioawk"
    conda: "conda_definitions/bioawk.yaml"
    shell: """

        bioawk -v sam={wildcards.sample} -c fastx '{{ print sam, $name, length($seq) }}' < {input} \
        > {output}

    """


rule gc_summary_individual:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output: "{out_base}/samples/{sample}/statistics/{sample}_gc.tsv"
    container: "docker://rocker/tidyverse" # remember to add devtools
    conda: "conda_definitions/r-tidyverse.yaml" # like r-markdown, but much simpler.
    params: base_variable = base_variable
    shell: """


        Rscript $ASSCOM2_BASE/scripts/tabseq_gc.r $ASSCOM2_BASE/scripts/tabseq_tiny.r {input} \
        > {output} 2> {output}.fail || echo what



    """




rule prokka_individual:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output:
        gff = "{out_base}/samples/{sample}/prokka/{sample}.gff",
        log = "{out_base}/samples/{sample}/prokka/{sample}.log",
        tsv = "{out_base}/samples/{sample}/prokka/{sample}.tsv",
        summarized_txt = "{out_base}/samples/{sample}/prokka/{sample}_summary.txt",
        labelled_tsv = "{out_base}/samples/{sample}/prokka/{sample}_labelled.tsv",
        labelled_gff = "{out_base}/samples/{sample}/prokka/{sample}_labelled.gff"

    container: "docker://staphb/prokka"
    conda: "conda_definitions/prokka.yaml"
    benchmark: "{out_base}/benchmarks/benchmark.prokka_individual.{sample}.tsv"
    resources:
        mem_mb = 8192
    threads: 4
    shell: """
      
        # For debugging of minced
        minced --version 
        
        prokka \
            --cpus {threads} \
            --force \
            --outdir {wildcards.out_base}/samples/{wildcards.sample}/prokka \
            --prefix {wildcards.sample} {input} \
            | tee {output.log} 

        # Label summary file
        cat {output.log} \
            | grep "Found" \
            | grep -E "tRNAs|rRNAs|CRISPRs|CDS|unique" \
            | cut -d" " -f 3,4 \
            | awk -v sam={wildcards.sample} '{{ print sam " " $0 }}' \
            >> {output.summarized_txt} # jeg undrer mig over hvorfor den har to gt question mark

        # Label tsv file
        cat {output.tsv} \
            | awk -v sam={wildcards.sample} '{{ print $0 "\t" sam }}' \
            > {output.labelled_tsv}


        # Remove fasta from gff and add sample label
        gff_fasta_start=$(grep --line-number --extended-regexp "^##FASTA" {output.gff} | cut -f1 -d:)
        head --lines $((-1+$gff_fasta_start)) {output.gff} \
            | awk -v sam={wildcards.sample} '{{ print $0 "\t" sam }}' \
            > {output.labelled_gff}



    """



rule kraken2_individual:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output: "{out_base}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv"
    container: "docker://staphb/kraken2"
    conda: "conda_definitions/kraken2.yaml"
    threads: 2
    resources:
        mem_mb = 65536
    benchmark: "{out_base}/benchmarks/benchmark.kraken2_individual.{sample}.tsv"
    shell: """


        if [ ! -z $ASSCOM2_KRAKEN2_DB ]; then
            echo using kraken2 database $ASSCOM2_KRAKEN2_DB

            # Run kraken2
            kraken2 \
                --threads {threads} \
                --db $ASSCOM2_KRAKEN2_DB \
                --report {output}_tmp \
                {input} \
                > /dev/null

            # Put sample names in front
            cat {output}_tmp \
            | awk -v sam={wildcards.sample} '{{ print sam "\t" $0 }}' \
            > {output}

            # Remove temp file
            rm {output}_tmp

        else
            echo "The ASSCOM2_KRAKEN2_DB variable is not set, and thus the kraken2 rule and its jobs will not be run. Consider using the scripts/set_up_kraken2.sh script for downloading and linking the latest kraken2 database."
        fi

    """

# This rule runs once, downloading the busco dataset that is needed for rule busco_individual.
# Make sure that this job is run on a node that has internet access.
rule busco_download:
    output:
        touch("{base_variable}/databases/busco/busco_download_done.flag")
    conda: "conda_definitions/busco.yaml"
    shell: """

        
        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{wildcards.base_variable}/databases/busco/busco_download_done.flag" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            # Busco is a bit stupid in the way that it requires an input file, but doesn't read it when you just download.
            touch dummy.fasta
            
            # https://busco.ezlab.org/busco_userguide.html#download-and-automated-update
            busco \
                --in dummy.fasta \
                --out shouldnotbenecessarytosettheoutdirwhenjustdownloading \
                --mode geno \
                --auto-lineage-prok \
                --force \
                --download_path {wildcards.base_variable}/databases/busco \
                --download prokaryota

            # Info: You can also swap "--download prokaryota" with "--download virus" if you're feeling adventurous ...
            
            touch {output}

            # Clean up 
            rm dummy.fasta
        
        fi

    """

rule busco_individual:
    input: 
        busco_download = expand("{base_variable}/databases/busco/busco_download_done.flag", base_variable = base_variable),
        fasta = "{out_base}/samples/{sample}/{sample}.fa"
    output: 
        flag = touch("{out_base}/samples/{sample}/busco/busco_done.flag"),
        table_extract = "{out_base}/samples/{sample}/busco/short_summary_extract.tsv"
    params:
        base_variable = base_variable,
        #out_base = out_base_var,
        out_dir = "{out_base}/samples/{sample}/busco"
    conda: "conda_definitions/busco.yaml"
    threads: 1
    resources:
        mem_mb = 8192,
        runtime = "06:00:00",
    shell: """


        # https://busco.ezlab.org/busco_userguide.html#offline
        busco \
            --cpu {threads} \
            --in {input.fasta} \
            --out {params.out_dir} \
            --mode geno \
            --auto-lineage-prok \
            --force \
            --tar \
            --download_path {params.base_variable}/databases/busco \
            --offline 
            #--download prokaryota
            
        # Potentially remove the comment block below in the case that rule busco_download works out well.
        ## To set it up the first time, swap "--offline" with "--download prokaryota"
        ## This will not run the analysis, but just download to the path set.
        ## Info: Obviously make sure to only download with a single job, otherwise the downloads will overlap and corrupt..
        ## Info: Consider adding virus on top of downloading prokaryota



        # Grab the most important output with the following grob-grep 
        cat {wildcards.out_base}/samples/{wildcards.sample}/busco/short_summary.*.txt \
        | grep -E "# The lineage dataset is:|# Summarized benchmarking|C:" \
        > {output.table_extract}


    """


# --- Collect results among all samples -----------------------------

rule busco:
    input: 
        metadata = "{out_base}/metadata.tsv",
        #tables = expand("{out_base}/samples/{sample}/busco/run_bacteria_odb10/short_summary_extract.tsv", out_base = out_base_var, sample = df["sample"]),
        tables = expand("{out_base}/samples/{sample}/busco/short_summary_extract.tsv", out_base = out_base_var, sample = df["sample"]),
    output: "{out_base}/collected_results/busco.tsv"
    resources: 
        mem_mb = 128,
        runtime = "00:10:00"
    shell: """

        cat {input.tables} >> {output}

        {void_report} # TODO: Make a nice summary in the report.

    """


rule kraken2:
    input: expand("{out_base}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv", out_base = out_base_var, sample = df["sample"]),
    output: "{out_base}/collected_results/kraken2_reports.tsv",
    resources:
        runtime = "01:00:00"
    shell: """

        # kraken2
        echo -e "sample\tmatch_percent\tclade_mappings\tlevel_mappings\tlevel\ttaxonomic_id\tclade" \
        > {output}

        cat {input} >> {output} 

        {void_report}
    """


rule sequence_lengths:
    input: expand("{out_base}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", out_base = out_base_var, sample = df["sample"])
    output: "{out_base}/collected_results/sequence_lengths.tsv"
    resources:
        runtime = "01:00:00"
    shell: """

        # Sequence lengths
        echo -e "sample\trecord\tlength" \
        > {output}

        cat {input} >> {output} 

        {void_report}
    """

rule gc_summary:
    input: expand("{out_base}/samples/{sample}/statistics/{sample}_gc.tsv", out_base = out_base_var, sample = df["sample"])
    output: "{out_base}/collected_results/GC_summary.tsv"
    resources:
        runtime = "01:00:00"
    shell: """

        # Sequence lengths
        echo -e "sample\tpart\tlength\tGC" \
        > {output}

        cat {input} | grep -vE "^#" >> {output} # Append content without headers

        {void_report}
    """






rule prokka:
    input:
        metadata = "{out_base}/metadata.tsv",
        summarized_txt = expand("{out_base}/samples/{sample}/prokka/{sample}_summary.txt", out_base = out_base_var, sample = df["sample"]),
        labelled_tsv = expand("{out_base}/samples/{sample}/prokka/{sample}_labelled.tsv", out_base = out_base_var, sample = df["sample"]),
        labelled_gff = expand("{out_base}/samples/{sample}/prokka/{sample}_labelled.gff", out_base = out_base_var, sample = df["sample"]),

    output: 
        summarized_txt = "{out_base}/collected_results/prokka_summarized.txt",
        labelled_tsv = "{out_base}/collected_results/prokka_labelled.tsv",
        labelled_gff = "{out_base}/collected_results/prokka_labelled.gff",

    resources:
        runtime = "01:00:00"
    shell: """

        echo "sample value name" > {output.summarized_txt}
        cat {input.summarized_txt} >> {output.summarized_txt} 

        cat {input.labelled_tsv} > {output.labelled_tsv}

        cat {input.labelled_gff} > {output.labelled_gff}


        {void_report}
    """







rule sample_pathway_enrichment_analysis:
    input: "{out_base}/collected_results/prokka_labelled.tsv"
    output: "{out_base}/collected_results/sample_pathway_enrichment_analysis.tsv"
    conda: "conda_definitions/r-clusterProfiler.yaml"
    shell: """


        Rscript $ASSCOM2_BASE/scripts/sample_pathway_enrichment_analysis.R $ASSCOM2_BASE/assets/ko {input} \
            > {output}

        {void_report}
    """







# --- Targets for the complete set below: ---------------------------

def get_mem_roary(wildcards, attempt): 
    return [32000, 64000, 128000][attempt-1]

rule roary:
    input: 
        metadata = "{out_base}/metadata.tsv",
        gff = expand("{out_base}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], out_base = out_base_var)
    output: ["{out_base}/roary/summary_statistics.txt", "{out_base}/roary/core_gene_alignment.aln", "{out_base}/roary/gene_presence_absence.csv", "{out_base}/roary/roary_done.flag"]
    params:
        blastp_identity = int(config['roary_blastp_identity']), # = 95 # For clustering genes
        core_perc = 99  # Definition of the core genome
    #conda: "envs/roary.yml"
    threads: 16
    #retries: 2
    resources:
        #mem_mb = 32768,
        mem_mb = get_mem_roary,
        runtime = "23:59:59" # Well, fuck me if this doesn't work on PBS
    container: "docker://sangerpathogens/roary"
    conda: "conda_definitions/roary.yaml"
    shell: """
    
        
        # Since I reinstalled conda, I've had problems with "Can't locate Bio/Roary/CommandLine/Roary.pm in INC". Below is a hacky fix
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0
        

        # Silence parallel's citation pester:
        echo "will cite" | parallel --citation > /dev/null 2> /dev/null

        # Roary is confused by the way snakemake creates directories ahead of time.
        # So I will delete it manually here before calling roary.
        rm -r {wildcards.out_base}/roary

        roary -a -r -e --mafft \
            -p {threads} \
            -i {params.blastp_identity} \
            -cd {params.core_perc} \
            -f {wildcards.out_base}/roary \
            {input.gff} || echo roary failed

        touch {output}
                
        {void_report}
    """


rule snp_dists:
    input: 
        metadata = "{out_base}/metadata.tsv",
        aln = "{out_base}/roary/core_gene_alignment.aln"
    output: "{out_base}/snp-dists/snp-dists.tsv"
    conda: "conda_definitions/snp-dists.yaml"
    container: "docker://staphb/snp-dists"
    shell: """

        snp-dists {input.aln} > {output}

        {void_report}
    """




rule assembly_stats:
    input: 
        metadata = "{out_base}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist()
    output: "{out_base}/assembly-stats/assembly-stats.tsv"
    container: "docker://sangerpathogens/assembly-stats"
    conda: "conda_definitions/assembly-stats.yaml"
    shell: """
        
        assembly-stats -t {input.fasta} > {output}

        {void_report}
    """


def get_mem_gtdbtk(wildcards, attempt): 
    return [150000, 300000, 400000, 500000][attempt-1]

rule gtdbtk:
    input: 
        metadata = "{out_base}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist()
    output: "{out_base}/gtdbtk/gtdbtk.bac.summary.tsv"
    params:
        batchfile_content = df[['input_file_fasta', 'sample']].to_csv(header = False, index = False, sep = "\t"),
        out_dir = "{out_base}/gtdbtk/"
    threads: 8
    #retries: 3
    resources:
        #mem_mb = 150000 # Last time I remember, it used 130000
        mem_mb = get_mem_gtdbtk
    conda: "conda_definitions/gtdbtk.yaml"
    benchmark: "{out_base}/benchmarks/benchmark.gtdbtk.tsv"
    shell: """

        echo "GTDBTK_DATA_PATH is $GTDBTK_DATA_PATH"

        # Create batchfile
        echo '''{params.batchfile_content}''' > {wildcards.out_base}/gtdbtk/batchfile.tsv

        gtdbtk classify_wf -h
        
        gtdbtk classify_wf \
            --batchfile {wildcards.out_base}/gtdbtk/batchfile.tsv \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --keep_intermediates \
            --force

        # Homogenize database version number
        cp {wildcards.out_base}/gtdbtk/gtdbtk.bac120.summary.tsv {output}


        {void_report}


    """ 



rule abricate:
    input: 
        metadata = "{out_base}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist()
    output:
        card_detailed = "{out_base}/abricate/card_detailed.tsv",
        card_sum = "{out_base}/abricate/card_summarized.tsv",
        plasmidfinder_detailed = "{out_base}/abricate/plasmidfinder_detailed.tsv",
        plasmidfinder_sum = "{out_base}/abricate/plasmidfinder_summarized.tsv",
        ncbi_detailed = "{out_base}/abricate/ncbi_detailed.tsv",
        ncbi_sum = "{out_base}/abricate/ncbi_summarized.tsv",
        vfdb_detailed = "{out_base}/abricate/vfdb_detailed.tsv",
        vfdb_sum = "{out_base}/abricate/vfdb_summarized.tsv"

    container: "docker://staphb/abricate"
    conda: "conda_definitions/abricate.yaml"
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
    mlst_scheme_interpreted = ""
else:
    mlst_scheme_interpreted = f"--scheme {config['mlst_scheme']}"
#print(f"Info: The mlst_scheme is set to <{mlst_scheme_interpreted}>") # Debug message.

rule mlst:
    input: 
        metadata = "{out_base}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist()
    output: "{out_base}/mlst/mlst.tsv",
    params:
        mlst_scheme_interpreted = mlst_scheme_interpreted,
        list_ = "{out_base}/mlst/mlst_schemes.txt"
    container: "docker://staphb/mlst"
    conda: "conda_definitions/mlst.yaml"
    shell: """

        mlst {params.mlst_scheme_interpreted} {input.fasta} > {output}

        mlst --list > {params.list_}

        {void_report}
    """




rule mashtree:
    input: 
        metadata = "{out_base}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist()
    output: 
        tree = "{out_base}/mashtree/mashtree.newick",
        dist = "{out_base}/mashtree/mash_dist.tsv"
    container: "docker://staphb/mashtree"
    conda: "conda_definitions/mashtree.yaml"
    threads: 4
    resources:
        mem_mb = 16000
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
    return [8000, 64000][attempt-1]


rule fasttree:
    input:
        metadata = "{out_base}/metadata.tsv",
        fasta = "{out_base}/roary/core_gene_alignment.aln"
    output: "{out_base}/fasttree/fasttree.newick"
    container: "docker://staphb/fasttree"
    conda: "conda_definitions/fasttree.yaml"
    threads: 4
    #retries: 1
    resources:
        mem_mb = get_mem_fasttree,
        runntime = "23:59:59"
    shell: """

        OMP_NUM_THREADS={threads}

        FastTree -nt -gtr {input.fasta} > {output} 2> {output}.log || echo "fasttree failed"

        touch {output}

        {void_report}
    """


rule fetch_report_template:
    output: "{out_base}/rmarkdown_template.rmd"
    shell: """

        cp $ASSCOM2_BASE/scripts/genomes_to_report_v2.Rmd {output}

        {void_report}
    """


# rule report:
#     input:
#         roary = "{out_base}/roary/roary_done.flag", # fasttree depends on roary, so the roary dependency is not necessary.
#         fasttree = "{out_base}/fasttree/fasttree.newick", 
#         snp_dists = "{out_base}/snp-dists/snp-dists.tsv",
#         rmarkdown_template = "{out_base}/rmarkdown_template.rmd"
#     #output: "{out_base}/report.html"
#     output: "{out_base}/report_{batch_title}.html"
#     params:
#         #markdown_template_rmd = "rmarkdown_template.rmd", # "genomes_to_report_v2.Rmd"
#         markdown_template_html = "genomes_to_report_v2.html"
#     container: "docker://cmkobel/assemblycomparator2_report"
#     conda: "conda_definitions/r-markdown.yaml"
#     shell: """

#         cd {wildcards.out_base}

#         Rscript -e 'library(rmarkdown); rmarkdown::render("rmarkdown_template.rmd", "html_document")'

#         rm rmarkdown_template.rmd
#         mv rmarkdown_template.html ../{output}
        
#     """



# This rule might seem silly, but it makes sure that the report environment is ready to rock when the report subpipeline eventually is run: This has two pros:
#    1) The vastly faster mamba configuration in the asscom2 pipeline is used
#    2) The conda/mamba debugging is taken care of, without having to wait for jobs to finish on fresh installations.
# Since all snakemake conda environments are installed in $SNAKEMAKE_CONDA_PREFIX set to ${ASSCOM2_BASE}/conda_base, reuse is guaranteed.
rule install_report_environment_aot:
    output: touch("{out_base}/.install_report_environment_aot.flag")
    conda: "report_subpipeline/conda_definitions/r-markdown.yaml"
    shell: """

        echo OK

    """

# Just a dummy rule if you wanna force the report
# assemblycomparator2 --until report
# TODO: Should never run on the queue system
rule report:
    resources:
        runtime = "00:01:00"
    shell: """
        {void_report}
    """


# Call the report subpipeline
report_call = f"""
    mkdir -p output_asscom2/logs; \
    snakemake \
        --snakefile $ASSCOM2_BASE/report_subpipeline/snakefile \
        --cores 4 \
        --use-conda \
        --config out_base=$(pwd)/output_asscom2 base_variable={base_variable} batch_title={batch_title} 2> output_asscom2/logs/report.err.log 
    """

onsuccess:
    print("onsuccess: calling report subpipeline ...")
    shell(report_call)

onerror:
    print("onerror: calling report subpipeline ...")
    shell(report_call)





print("*/")







