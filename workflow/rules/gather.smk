# These rules make sure that all genomes, from ncbi or local, are ready and copied to a local location for further analysis.

# Only runs if there is one or more ncbi accessions added.
# Finally I decided to run this in sample mode. So this needs to output sample specific files, and then it needs to use the resources load system to define maxiumum concurrent jobs. Using batch mode only makes sense if I implement my own cache, but that is way outside the scope of comparem2 and I'm just trying to make this work robustly, not prematurely quick..
# So the bottom line is that each genome needs to unpack in its own directory, and then the annotation rule needs to use that sample-specific path.
# This is getting better and better.
rule get_ncbi: # Per sample
    output: 
        marker = "{output_directory}/samples/{sample}/ncbi/ncbi_dataset/data/dataset_catalog.json"
        #"results_comparem2/samples/{wildcards.sample}/ncbi_dataset/sequence_report.jsonl"
        #accessions = "{output_directory}/samples/{sample}/ncbi/accessions.tsv"
    conda: "../envs/ncbi_datasets.yaml"
    params:
        path_zip = "{output_directory}/samples/{sample}/ncbi/ncbi_dataset.zip",
        path_decompressed = "{output_directory}/samples/{sample}/ncbi",
        #output_directory = output_directory,
        #accessions_joined_newline = "\n".join(df[df["origin"] == "ncbi"]["sample"].tolist()),
        #samples = df[df["origin"] == "ncbi"]["sample"].tolist(),
        #accession = df[df["origin"] == "ncbi"]["sample"].tolist()[0],
    resources: # Using resources actually worked quite well, but retries is just much simpler and quicker.
        downloads = 1 
    retries: 1 # Just in case of an intermittent server error.
    shell: """
        
        # Download 
        datasets \
            download genome \
            accession {wildcards.sample} \
            --filename {params.path_zip} \
            --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report
            
        # Unzip
        echo "A" | unzip -q {params.path_zip} -d {params.path_decompressed}
        
        
        # mv
        mv {wildcards.output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/{wildcards.sample}/* {wildcards.output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/
        
        # Tidy up
        rm -r {wildcards.output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/{wildcards.sample}
        rm {params.path_zip}
        

    """
    


def conditional_ncbi_dependency(wildcards):
    """ For rule copy
        Rule copy only needs the ncbi dependency if there are any ncbi accessions added to the add_ncbi config-parameter. 
    """
    conditional_dependency_path = "{output_directory}/samples/{sample}/ncbi/ncbi_dataset/data/dataset_catalog.json"

    origin = df[df["sample"] == wildcards.sample]["origin"].tolist()[0]

    if origin == "ncbi":
        return conditional_dependency_path
    else:
        return list()


rule copy: # Per sample
    input:
        conditional_ncbi_dependency
        #lambda wildcards: "{output_directory}/ncbi_dataset/ncbi_dataset.zip" if "ncbi" in df["origin"].tolist() else list()
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
            cp {wildcards.output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/{wildcards.sample}*_genomic.fna {params.input_file}

                
            """)
            
        # Finally, for both:
        shell("""
                
            cp {params.input_file:q} {output.input_file_copy:q}
            md5sum {output.input_file_copy:q} > {output.log:q}
            echo "# Copied from $(realpath {params.input_file}) on $(date)" >> {output.log:q}
            
        """)
        