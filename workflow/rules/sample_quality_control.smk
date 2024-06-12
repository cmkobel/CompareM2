# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
# Should I rename this to "rule any2fasta" just to make it more transparent? 
rule copy:
    input: 
        genome = lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0],
    output: 
        fasta = "{output_directory}/samples/{sample}/{sample}.fna",
        md5sum = "{output_directory}/samples/{sample}/{sample}.md5.txt",
    conda: "../envs/any2fasta.yaml"
    threads: 1 # Weirdly, or bugly, there must be a thread n definition in the rule. Otherwise, the set-threads option (in the orion profile) will not be taken up. 
    resources:
        mem_mb = 256,
        runtime = "10m",
    shell: """

        # Collect version number.
        any2fasta -v > "$(dirname {output.fasta})/.software_version.txt"
        
        any2fasta {input.genome:q} > {output.fasta:q}
        
        md5sum {output.fasta:q} > {output.md5sum:q}
        
    """  


rule assembly_stats:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{output_directory}/assembly-stats/assembly-stats.tsv"
    conda: "../envs/assembly-stats.yaml"
    benchmark: "{output_directory}/benchmarks/benchmarks.assembly_stats.tsv"
    shell: """
    
        # Collect version number.
        echo "assembly-stats $(assembly-stats -v)" > "$(dirname {output})/.software_version.txt"
        
        assembly-stats -t {input.fasta:q} > {output:q}

        {void_report}

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
        
        seqkit fx2tab {input.assembly:q} -l -g -G -n -H \
        > {output:q}

    """


rule busco:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        database_representative = DATABASES + "/busco/ac2_busco_database_representative.flag",
        faa = "{output_directory}/samples/{sample}/prokka/{sample}.faa",        
    output: 
        flag = touch("{output_directory}/samples/{sample}/busco/busco_done.flag"),
        table_extract = "{output_directory}/samples/{sample}/busco/short_summary_extract.tsv"
    params:
        base_variable = base_variable,
        database_path = DATABASES + "/busco",
        out_dir = "{output_directory}/samples/{sample}/busco",
    conda: "../envs/busco.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.busco_sample.{sample}.tsv"
    threads: 1 # Because run_sepp hangs for a long time, not doing anything, I'd rather have more processes started on any type CPU.
    resources:
        mem_mb = 8192,
        runtime = "1h",
    shell: """

        # Busco fails because of a problem with the sepp package. This doesn't really matter as we just want the completeness results.
        # But, this means that we need a hacky workaround to let this job exit gracefully (exit code 0) on the basis of whether any completeness results have been written to disk.
        # Hence, the actual exit code of busco, we will ignore.
        
        # Collect version number.
        busco -v > "$(dirname {output.table_extract})/.software_version.txt"
        
        # Collect database version.
        echo -e "$(date -Iseconds)\t{input.database_representative}" > "$(dirname {output.table_extract})/.database_version.txt"

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
        cat {wildcards.output_directory}/samples/{wildcards.sample}/busco/auto_lineage/*/short_summary.json \
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
