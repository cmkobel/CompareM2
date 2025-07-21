# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
# Should I rename this to "rule any2fasta" just to make it more transparent? Or just remove it. I think it creates as many problems as it solves.




 
# Generates the input_file later used by rule copy.
rule get_refseq:
    output:
        assembly = "{output_directory}/.refseq_downloads/{sample}/{sample}.fna",
        gff = "{output_directory}/.refseq_downloads/{sample}/ncbi_dataset/data/{sample}/renamed/{sample}.gff",
        faa = "{output_directory}/.refseq_downloads/{sample}/ncbi_dataset/data/{sample}/renamed/{sample}.faa",
        log = "{output_directory}/.refseq_downloads/{sample}/ncbi_dataset/data/{sample}/renamed/{sample}.log",
        ffn = "{output_directory}/.refseq_downloads/{sample}/ncbi_dataset/data/{sample}/renamed/{sample}.ffn", # better check this one 
        tsv = "{output_directory}/.refseq_downloads/{sample}/ncbi_dataset/data/{sample}/renamed/{sample}.tsv", 
    conda: "../envs/ncbi_datasets.yaml"
    shell: """


        datasets \
            download genome \
            accession {wildcards.sample} \
            --filename "$(dirname {output.assembly})/ncbi_{wildcards.sample}.zip" \
            --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report
            
        echo "A" | unzip "$(dirname {output.assembly})/ncbi_{wildcards.sample}.zip" -d "$(dirname {output.assembly})"
        
        # Make a predictable path (The assembly accession is unpredictable)
        cp $(dirname {output.assembly})/ncbi_dataset/data/{wildcards.sample}/{wildcards.sample}_*_genomic.fna {output.assembly}
        
        # Prepare annotation files for later.
        cp $(dirname {output.assembly})/ncbi_dataset/data/{wildcards.sample}/genomic.gff {output.gff}
        cp $(dirname {output.assembly})/ncbi_dataset/data/{wildcards.sample}/protein.faa {output.faa}
        cp $(dirname {output.assembly})/ncbi_dataset/data/{wildcards.sample}/sequence_report.jsonl {output.log}
        cp $(dirname {output.assembly})/ncbi_dataset/data/{wildcards.sample}/cds_from_genomic.fna {output.ffn}
        cp $(dirname {output.assembly})/ncbi_dataset/data/{wildcards.sample}/genomic.gtf {output.tsv}
        
        
        
        
    """


# Makes an adjacent copy of the input file for later reference.
rule copy:
    input: 
        genome = lambda wildcards: df[df["sample"] == wildcards.sample]["input_file"].tolist()
    output: 
        fasta = "{output_directory}/samples/{sample}/{sample}.fna",
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
        echo "Copied from $(realpath {input}) on $(date)" > {output.log:q}
        md5sum {output.fasta:q} >> {output.log:q}
        
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


