
# One new main copy that uses the origin column to decide how to gain the genome. Either locally from disk or via refseq ncbi-datasets-cli download.
# One problem I haven't solved is how to deal with inputs, as refseqqed genomes can't have an input
rule gather:
    output: 
        assembly = ensure("{output_directory}/samples/{sample}/{sample}.fna", non_empty = True),
        log = "{output_directory}/samples/{sample}/{sample}.log",
    retries: 2
    params: 
        origin = lambda wildcards: df[df["sample"] == wildcards.sample]["origin"].tolist()[0],
        input_file = lambda wildcards: df[df["sample"] == wildcards.sample]["input_file"].tolist()[0], # Only used for genomes of local origin. Value should be NaN for refseq genomes.
    run: # Python code with conditionals to control gathering 
        if params.origin == "local": # Just copy the local file
            shell("""
                
                cp {params.input_file:q} {output.assembly:q}
                
            """)
        
        elif params.origin == "refseq": # Download the file over the internet.
            shell("""
            
            zip_path={output_directory}/samples/{wildcards.sample}/ncbi_{wildcards.sample}.zip
            
            datasets \
                download genome \
                accession {wildcards.sample} \
                --filename $zip_path \
                --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report
                
            echo "A" | unzip $zip_path -d {output_directory}/samples/{wildcards.sample}/refseq_download
            rm $zip_path
        
            # Link files to more reasonable paths
            # The asterisk is because the assembly name is unpredictable.
            cp {output_directory}/samples/{wildcards.sample}/refseq_download/ncbi_dataset/data/{wildcards.sample}/{wildcards.sample}_*_genomic.fna {output.assembly}
                
            """)
            
        # Finally, for both:
        shell("""
            
            md5sum {output.assembly:q} > {output.log:q}
            echo "# Copied from $(realpath {params.input_file}) on $(date)" >> {output.log:q}
        
        """)
   
    


rule sequence_lengths:
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


