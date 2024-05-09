
# aka eggnog-mapper
rule eggnog:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/eggnog/ac2_eggnog_database_representative.flag",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{results_directory}/samples/{sample}/eggnog/{sample}.emapper.gff",
        ffn = "{results_directory}/samples/{sample}/eggnog/{sample}.emapper.fasta", # Used in dbcan, interproscan, diamond_kegg, motupan
        tsv = "{results_directory}/samples/{sample}/eggnog/{sample}.emapper.annotations",
    #params:
            
    conda: "../envs/eggnog.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.eggnog_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 8
    shell: """
        
        # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#basic-usage
        
        emapper.py \
            -m diamond \
            --data_dir $(dirname {input.database_representative}) \
            --itype genome \
            --override \
            --cpu {threads} \
            --output_dir "$(dirname {output.gff})/" \
            -o "{wildcards.sample}" \
            -i {input.assembly:q} 
            
        #touch {output} # DEBUG

        {void_report}

    """
    