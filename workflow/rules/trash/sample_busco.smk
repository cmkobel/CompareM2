
rule busco:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        #busco_download = expand("{base_variable}/databases/busco/file_versions.tsv", base_variable = base_variable), # This is a bad idea, because it requires a complete reinstall if snakemake somehow removes the file, which is quite likely.
        database_representative = DATABASES + "/busco/ac2_busco_database_representative.flag", # Should point to the directory where the following files reside: "file_versions.tsv  lineages/  placement_files/"
        #fasta = "{results_directory}/samples/{sample}/{sample}.fna",
        faa = "{results_directory}/samples/{sample}/prokka/{sample}.faa",

        
    output: 
        flag = touch("{results_directory}/samples/{sample}/busco/busco_done.flag"),
        table_extract = "{results_directory}/samples/{sample}/busco/short_summary_extract.tsv"
    params:
        base_variable = base_variable,
        #results_directory = results_directory,
        database_path = DATABASES + "/busco", # Was {params.base_variable}/databases/busco
        out_dir = "{results_directory}/samples/{sample}/busco",
    conda: "../envs/busco.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.busco_sample.{sample}.tsv"
    threads: 1 # Because run_sepp hangs for a long time, not doing anything, I'd rather have more processes started on any CPU.
    resources:
        mem_mb = 8192,
        runtime = "1h",
    shell: """

        # Busco fails because of a problem with the sepp package. This doesn't really matter as we just want the completeness results.
        # But, this means that we need a hacky workaround to let this job exit gracefully (exit code 0) on the basis of whether any completeness results have been written to disk.
        # Hence, the actual exit code of busco, we will ignore.

        # https://busco.ezlab.org/busco_userguide.html#offline
        # Is the timeout bug fixed? Update: nope.
        timeout 1800 \
            busco \
                --cpu {threads} \
                --in {input.faa:q} \
                --out {params.out_dir} \
                --mode protein \
                --auto-lineage-prok \
                --force \
                --tar \
                --skip_bbtools \
                --download_path {params.database_path} \
                --offline || (>&2 echo "ac2: Busco failed internally.")

        # Cat all auto lineage results together or create empty file
        # The following cat command will fail if the glob doesn't resolve any files: This is the wanted behaviour.
        cat {wildcards.results_directory}/samples/{wildcards.sample}/busco/auto_lineage/*/short_summary.json \
        > "{output.table_extract}_temp"
        
        >&2 echo "ac2: Intermediate results have been collected, as they are still useful."
        
        # Extract relevant features
        cat "{output.table_extract}_temp" \
        | grep -oE "(\\"in\\"|\\"name\\"|\\"one_line_summary\\").+" \
        > {output.table_extract:q}

        # Clean up
        rm "{output.table_extract}_temp" || echo "no file"
        
        mv busco_*.log .snakemake/log || echo "no logs to move"

        {void_report}

    """
