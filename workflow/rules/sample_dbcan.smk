
rule dbcan: # I can't decide whether this rule should really be called "run_dbcan", since that is the name of the software.
    input: 
        metadata = "{results_directory}/metadata.tsv",
        aminoacid = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # From prokka
        database_representative = DATABASES + "/dbcan/ac2_dbcan_database_representative.flag"
    output: 
        overview_table = "{results_directory}/samples/{sample}/dbcan/overview.txt",
        diamond_table = "{results_directory}/samples/{sample}/dbcan/diamond.out"
    params: 
        out_dir = "{results_directory}/samples/{sample}/dbcan"
    conda: "../envs/dbcan.yaml" # Not sure if it should be called by a version number?
    benchmark: "{results_directory}/benchmarks/benchmark.dbcan.{sample}.tsv"
    threads: 8
    resources: 
        mem_mb = 8000
    shell: """

        # It seems to be necessary to set all the cpu thread counts manually.

        export HMMER_NCPU={threads}
        
        run_dbcan -h > $(dirname {output.overview_table:q})/run_dbcan_version.txt

        run_dbcan \
            --dbcan_thread {threads} \
            --dia_cpu {threads} \
            --hmm_cpu {threads} \
            --tf_cpu {threads} \
            --stp_cpu {threads} \
            --db_dir $(dirname {input.database_representative:q}) \
            --out_dir {params.out_dir:q} \
            {input.aminoacid:q} \
            protein 

    """
