# These rules make sure that all genomes, from ncbi or local, are ready and copied to a local location for further analysis.


# Used to generate the ncbi cache marker path for a single sample/accession. The assembly file has an unpredictable path so we use this .jsonl file instead.
ncbi_marker_path = lambda sample: f"{output_directory}/.ncbi_cache/accessions/{sample}/sequence_report.jsonl"

ncbi_cache_misses = [sample for sample in df[df["origin"] == "ncbi"]["sample"].tolist() if not os.path.isfile(ncbi_marker_path(sample))]   
#print("misses", ncbi_cache_misses) # debug

# rule get_ncbi downloads to a cache, from which copy copies from. The cache is flat and not arranged by download.
rule get_ncbi: # Per sample
    output: 
        marker = [ncbi_marker_path(sample) for sample in ncbi_cache_misses],
    conda: "../envs/ncbi_datasets.yaml"
    params:
        path_zip = f"{output_directory}/.ncbi_cache/download/ncbi_dataset.zip", # Deleted after decompression.
        path_decompressed = f"{output_directory}/.ncbi_cache/download/decompressed",
        accessions_tsv = f"{output_directory}/.ncbi_cache/download/accessions.tsv",
        accessions_joined_newline = "\n".join(ncbi_cache_misses), 
        n_accessions = len(ncbi_cache_misses)
    shell: """
    
        # Since this rule has no inputs or outputs when there is no missing accessions, we need to handle running without accessions.
        if [ {params.n_accessions} -gt 0 ]; then
            
            echo CM2 rule get_ncbi: Number of accessions to download: {params.n_accessions}
    
            mkdir -p {output_directory}/.ncbi_cache/download
            mkdir -p {output_directory}/.ncbi_cache/accessions/
            
            datasets --version > {output_directory}/.ncbi_cache/.software_version.txt

            echo '''{params.accessions_joined_newline}''' > {params.accessions_tsv}
            
            # Download 
            datasets \
                download genome \
                accession --inputfile {params.accessions_tsv} \
                --filename {params.path_zip} \
                --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report
                
            # Unzip archive
            echo "Unzipping ..."
            echo "A" | unzip -q {params.path_zip} -d {params.path_decompressed}
            
            # cp metadata into dirs
            for dir in {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/*/; do [ -d "$dir" ] && cp {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/assembly_data_report.jsonl "$dir" ; done
            for dir in {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/*/; do [ -d "$dir" ] && cp {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/dataset_catalog.json "$dir" ; done
            
            # Move to a flat hierarchy
            echo "Moving ..."
            
            
            # Instead, manually generate the file list and mention the ones that have missing sequence reports (like GCA_902373455.1 which is suppressed)
            ## declare an array variable


            # Sometimes an accession can not be downloaded from NCBI. This may be because it has been suppressed due to low biological quality, or because there may be a server error at NCBI. In any case, there is no need to waste the full download by failing this job, when there are errors. Therefore, the while loop below continues and creates empty sequence reports when genomes are non-downloadable. The downstream jobs will fail explicitly because there is no genome input, so the user will not have to worry about overlooking missing genomes.
            
            while IFS= read -r accession; do
            
                #echo "$accession"
                
                dir={output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/${{accession}}
                
                if [ -d $dir ]; then
                    #echo content existsis
                    # Copy content
                    # Create dir and copy all content 
                    mkdir -p {output_directory}/.ncbi_cache/accessions/${{accession}}
                    cp -f ${{dir}}/* {output_directory}/.ncbi_cache/accessions/${{accession}}
                    
                    # Add accession to list of successes (so the user can rerun the analysis without the ones that fail downloading.)
                    echo $accession >> {output_directory}/.ncbi_cache/successful_transfers.tsv
                else
                    # "Throw" warning and continue
                    # It is very important that we don't "waste" the download as it is time-consuming. Therefore it is better that the pipeline later fails because the genome is missing anyway. Therefore this rule succeeds and "continues" even when some genomes are not downloadable.
                    
                    warning_msg="Warning: rule get_ncbi: Accession $accession is not available for download. This may be because the genome has been suppressed due to low biological quality, or because there is an error at NCBI. Please inspect https://www.ncbi.nlm.nih.gov/datasets/genome/${{accession}}/ for further info. CompareM2 will proceed without this genome."
                    echo $warning_msg
                    echo "Please rerun CompareM2 without this accession. This can be easily done by using {output_directory}/.ncbi_cache/successful_transfers.tsv as input." 
                    echo $accession >> {output_directory}/.ncbi_cache/failed_transfers.tsv
                    
                    # Create empty file to continue analysis
                    mkdir -p results_comparem2/.ncbi_cache/accessions/${{accession}}
                    
                    warning_msg_jsonl="{{'assemblyAccession':'$accession','CM2_warning_message':'$warning_msg'}}"
                    echo $warning_msg_jsonl > results_comparem2/.ncbi_cache/accessions/${{accession}}/sequence_report.jsonl
                    
                fi
                
            done < {params.accessions_tsv}
            
            
            # Tidy up
            rm -r {output_directory}/.ncbi_cache/download
        
        else
            echo CM2 rule get_ncbi: Nothing to download
            exit 0
            
        fi
        
    """
    


def conditional_ncbi_dependency(wildcards):
    """ For rule copy
        Rule copy only needs the ncbi dependency if there are any ncbi accessions added to the add_ncbi config-parameter. 
        Because df["input_file"] for ncbi accessions does not exist yet (because of unpredictable path), we need to use a marker instead.
    """

    conditional_dependency_path = ncbi_marker_path(wildcards.sample)
    
    origin = df[df["sample"] == wildcards.sample]["origin"].tolist()[0]

    if origin == "ncbi":
        return conditional_dependency_path
    else:
        return list()


rule copy: # Per sample
    input:
        conditional_ncbi_dependency
        #df["input_file"].tolist()
        #lambda wildcards: ncbi_marker_path(sample) for sample in df[df["sample"] == wildcards.sample]["sample"] else list()
    output: 
        input_file_copy = ensure("{output_directory}/samples/{sample}/{sample}.fna", non_empty = True),
        log = "{output_directory}/samples/{sample}/{sample}.log",
    params: 
        origin = lambda wildcards: df[df["sample"] == wildcards.sample]["origin"].tolist(),
        input_file = lambda wildcards: df[df["sample"] == wildcards.sample]["input_file"].tolist(), 
    #resources: # Using resources actually worked quite well, but retries is just much simpler and quicker.
    #    downloads = 1
    log: "{output_directory}/samples/{sample}/{sample}.log"
    run: # Python code with conditionals to control gathering 
    
        if params.origin == ["ncbi"]: # Download the file over the internet.
            shell("""
            
            # First, make the file name predictable
            cp {output_directory}/.ncbi_cache/accessions/{wildcards.sample}/{wildcards.sample}*_genomic.fna {params.input_file}
            

                
            """)
            
        # Finally, for both:
        shell("""
                
            cp {params.input_file:q} {output.input_file_copy:q}
            md5sum {output.input_file_copy:q} > {output.log:q}
            echo "# Copied from $(realpath {params.input_file}) on $(date)" >> {output.log:q}
            
        """)
        