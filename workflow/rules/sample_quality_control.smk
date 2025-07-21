# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
# Should I rename this to "rule any2fasta" just to make it more transparent? Or just remove it. I think it creates as many problems as it solves.




 
# Generates the input_file later used by rule copy.
rule get_refseq:
    output:
        fasta = touch("{output_directory}/.refseq_downloads/{sample}/{sample}.fna"),
    conda: "../envs/ncbi_datasets.yaml"
    shell: """

        # Reset directory in case of previous partial downloads.    
        rm -rf $(dirname {output.fasta})/* || echo "Continuing ..."

        datasets \
            download genome \
            accession {wildcards.sample} \
            --filename "$(dirname {output})/ncbi_{wildcards.sample}.zip" \
            --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report
            
        unzip "$(dirname {output.fasta})/ncbi_{wildcards.sample}.zip" -d "$(dirname {output.fasta})"
        
        # Make a predictable path (The assembly accession is unpredictable)
        cp $(dirname {output.fasta})/ncbi_dataset/data/{wildcards.sample}/{wildcards.sample}_*_genomic.fna $(dirname {output.fasta})/ncbi_dataset/data/{wildcards.sample}/{wildcards.sample}_assembly.fna
        
    """


# Makes an adjacent copy of the input file for later reference.
rule copy:
    input: 
        genome = lambda wildcards: df[df["sample"] == wildcards.sample]["input_file"].tolist()
    output: 
        fasta = "{output_directory}/samples/{sample}/{sample}.fna",
        md5sum = "{output_directory}/samples/{sample}/{sample}.md5.txt",
        log = "{output_directory}/samples/{sample}/{sample}.log",
    #conda: "../envs/assembly-stats.yaml"
    params:
        s = "e"
        #i = wildcards.input_file,
    threads: 1 # Weirdly, or bugly, there must be a thread n definition in the rule. Otherwise, the set-threads option (in the orion profile) will not be taken up. 
    resources:
        mem_mb = 256,
        runtime = "10m",
    shell: """
    
        echo {wildcards}
        cp {input.genome:q} {output.fasta:q}
        md5sum {output.fasta:q} > {output.md5sum:q}
        echo "Copied from $(realpath {input}) on $(date)" > {output.log}
        
    """  
    
    
   
    


rule sequence_lengths:
    input:
        metadata = "{output_directory}/metadata.tsv",
        assembly = "{output_directory}/samples/{sample}/{sample}.fna", 
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


