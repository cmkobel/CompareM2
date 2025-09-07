# These rules make sure that all genomes, from ncbi or local, are ready and copied to a local location for further analysis.


# Used to generate the ncbi cache marker path for a single sample/accession.
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
        accessions_joined_comma = ",".join(ncbi_cache_misses), 
        n_accessions = len(ncbi_cache_misses)
    shell: """
    
        # Since this rule has no inputs or outputs when there is no missing accessions, we need to handle running without accessions.
        if [ {params.n_accessions} -gt 0 ]; then
            
            echo CM2 rule get_ncbi: Number of accessions to download: {params.n_accessions}
    
            mkdir -p {output_directory}/.ncbi_cache/download
            mkdir -p {output_directory}/.ncbi_cache/accessions/
            
            datasets --version > {output_directory}/.ncbi_cache/.software_version.txt

            
            # Download 
            datasets \
                download genome \
                accession {params.accessions_joined_comma} \
                --filename {params.path_zip} \
                --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report
                
            # Unzip archive
            echo "A" | unzip -q {params.path_zip} -d {params.path_decompressed}
            
            
            # cp metadata into dirs
            for dir in {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/*/; do [ -d "$dir" ] && cp {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/assembly_data_report.jsonl "$dir" ; done
            for dir in {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/*/; do [ -d "$dir" ] && cp {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/dataset_catalog.json "$dir" ; done
            
            
            
            # Move to a flat hierarchy
            mv {output_directory}/.ncbi_cache/download/decompressed/ncbi_dataset/data/*/ {output_directory}/.ncbi_cache/accessions/
            
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
    """
    #conditional_dependency_path = "{output_directory}/samples/{sample}/ncbi/ncbi_dataset/data/dataset_catalog.json"
    #conditional_dependency_path = "{output_directory}/.ncbi_cache/accessions/{sample}/sequence_report.jsonl"
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
        