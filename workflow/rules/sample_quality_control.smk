
# One new main copy that uses the origin column to decide how to gain the genome. Either locally from disk or via refseq ncbi-datasets-cli download.
# One problem I haven't solved is how to deal with inputs, as refseqqed genomes can't have an input
# I'm actually considering going back to the old system where the ncbi accessions had their own rule. Then it can also pool and run quicker, have its own environment and resource settings. #129

ncbi_input_files = df[df["origin"] == "refseq"]["input_file"].tolist() # Empty list in case of no ncbi accessions.

rule ncbi: # batch
    output: 
        ncbi_input_files
    #conda: "../envs/ncbi_datasets.yaml"
    params:
        accessions_joined_comma = ",".join(df[df["origin"] == "refseq"]["sample"].tolist()),
        samples = df[df["origin"] == "refseq"]["sample"].tolist(),
        output_directory = output_directory
    shell: """
            
        path_zip="{params.output_directory}/ncbi_dataset/ncbi_dataset.zip"
        path_decompressed="{params.output_directory}/ncbi_dataset/decompressed"
        
        mkdir -p $path_decompressed
        
        # Download 
        datasets \
            download genome \
            accession "{params.accessions_joined_comma}" \
            --filename $path_zip \
            --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report
            
        # Unzip
        echo "A" | unzip $path_zip -d $path_decompressed
        
        
        # Old school but works way of making the file names predictable.
        # This can also be moved into the gather rule, but then you need a flag for making the correct dependency.
        declare -a arr=({params.samples:q}) # Example: declare -a arr=("element1" "element2" "element3")
        for i in "${{arr[@]}}"; do
            source="results_comparem2/ncbi_dataset/decompressed/ncbi_dataset/data/${{i}}/${{i}}*_genomic.fna"
            destination="results_comparem2/refseq_download/decompressed/ncbi_dataset/data/${{i}}/${{i}}_assembly.fna"
            cp $source $destination || touch $destination # In case an accession is missing its assembly, we can let the others run.
        done
        
        # This output should mimick the one in the df["input_file"]

    """

rule gather: # Before this rule there needs to be something that downloads the ncbi package in a single rule. But how do I make the dependency? It needs to be a flag?
    input:
        df["input_file"].tolist()
    output: 
        assembly = ensure("{output_directory}/samples/{sample}/{sample}.fna", non_empty = True),
        log = "{output_directory}/samples/{sample}/{sample}.log",
    retries: 2
    params: 
        origin = lambda wildcards: df[df["sample"] == wildcards.sample]["origin"].tolist(),
        input_file = lambda wildcards: df[df["sample"] == wildcards.sample]["input_file"].tolist(), 
    #resources: # Using resources actually worked quite well, but retries is just much simpler and quicker.
    #    downloads = 1
    log: "{output_directory}/samples/{sample}/{sample}.log"
    run: # Python code with conditionals to control gathering 
    
        if params.origin == "refseq": # Download the file over the internet.
            shell("""
            
            
            # First, make the file name predictable
            cp {output_directory}/samples/{wildcards.sample}/refseq_download/ncbi_dataset/data/{wildcards.sample}/{wildcards.sample}_*_genomic.fna {params.input_file}

                
            """)
            
        # Finally, for both:
        shell("""
                
            cp {params.input_file:q} {output.assembly:q}
            #md5sum {output.assembly:q} > {output.log:q}
            #echo "# Copied from $(realpath {params.input_file}) on $(date)" >> {output.log:q}
            
        """)
        
    


rule sequence_lengths: # TODO: rename to seqkit?
    input:
        assembly = "{output_directory}/samples/{sample}/{sample}.fna", 
        #metadata = "{output_directory}/metadata.tsv", # Only batch jobs need metadata (for proper updating)
    output: "{output_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"
    threads: 1
    resources:
        runtime = "60m",
    conda: "../envs/seqkit.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.sequence_lengths_sample.{sample}.tsv"
    shell: """
    
        # Collect version number.
        echo "seqkit $(seqkit -h | grep 'Version:')" > "$(dirname {output})/.software_version.txt"
        
        seqkit fx2tab \
            {input.assembly:q} \
            -l \
            -g \
            -G \
            -n \
            -H \
        > {output:q}

    """


