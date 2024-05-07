
rule bakta:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/bakta/ac2_bakta_database_representative.flag",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{results_directory}/samples/{sample}/bakta/{sample}.gff",
        faa = "{results_directory}/samples/{sample}/bakta/{sample}.faa",
        
        #gff_generic = "{results_directory}/samples/{sample}/annotation/{sample}.gff3",
    params:
        DATABASES = DATABASES
    conda: "../envs/bakta.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.bakta_sample.{sample}.tsv"
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
            
        # Optimize compatibility between prokka and bakta to make fluent use of either result easy.
        cp {output.gff}3 {output.gff}
        
        {void_report}

    """


