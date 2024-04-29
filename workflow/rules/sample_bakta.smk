
rule bakta:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/bakta/ac2_bakta_database_representative.flag",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{results_directory}/samples/{sample}/bakta/{sample}.gff3",
        faa = "{results_directory}/samples/{sample}/bakta/{sample}.faa", # Used in dbcan, interproscan, diamond_kegg, motupan
        
        #gff_generic = "{results_directory}/samples/{sample}/annotation/{sample}.gff3",
    params:
        DATABASES = DATABASES
    conda: "../envs/bakta.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.bakta_individual.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
                
        bakta \
            --db {params.DATABASES}/bakta/db \
            --output $(dirname {output.gff}) \
            --threads {threads} \
            --force \
            {input.assembly}
        
        {void_report}

    """


