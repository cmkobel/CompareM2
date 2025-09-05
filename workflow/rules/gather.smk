# These rules make sure that all genomes, from ncbi or local, are ready and copied to a local location for further analysis.

# Only runs if there is one or more ncbi accessions added.
rule ncbi_dataset: # Per batch
    output: 
        zip = "{output_directory}/ncbi_dataset/ncbi_dataset.zip", # Conditional dependency for rule gather
        accessions = "{output_directory}/ncbi_dataset/accessions.tsv" # And then this file will also be a dependency for copy, which will fix solve the original problem in #132
    conda: "../envs/ncbi_datasets.yaml"
    params:
        accessions_joined_newline = "\n".join(df[df["origin"] == "ncbi"]["sample"].tolist()),
        samples = df[df["origin"] == "ncbi"]["sample"].tolist(),
        output_directory = output_directory
    shell: """
            
        path_zip="{params.output_directory}/ncbi_dataset/ncbi_dataset.zip"
        path_decompressed="{params.output_directory}/ncbi_dataset/decompressed"
        
        mkdir -p $path_decompressed
        
        echo '''{params.accessions_joined_newline}''' > {output.accessions}
        
        # Download 
        datasets \
            download genome \
            accession --inputfile {output.accessions} \
            --filename $path_zip \
            --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report
            

            
        # Unzip
        echo "A" | unzip -q $path_zip -d $path_decompressed 

    """
    


rule copy: # Per sample
    input:
        lambda wildcards: "{output_directory}/ncbi_dataset/ncbi_dataset.zip" if "ncbi" in df["origin"].tolist() else list()
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
            cp results_comparem2/ncbi_dataset/decompressed/ncbi_dataset/data/{wildcards.sample}/{wildcards.sample}*_genomic.fna {params.input_file}
                
            """)
            
        # Finally, for both:
        shell("""
                
            cp {params.input_file:q} {output.input_file_copy:q}
            md5sum {output.input_file_copy:q} >> {output.log:q}
            #echo "# Copied from $(realpath {params.input_file}) on $(date)" >> {output.log:q}
            
        """)
        