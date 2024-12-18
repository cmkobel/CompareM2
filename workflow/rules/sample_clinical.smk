
rule amrfinder:
    input:
        metadata = "{output_directory}/metadata.tsv",
        faa = "{output_directory}/samples/{sample}/.annotation/{sample}.faa",
        database_representative = DATABASES + "/amrfinder/comparem2_amrfinder_database_representative.flag"
    output: 
        table = "{output_directory}/samples/{sample}/amrfinder/{sample}_amrfinder.tsv",
    params:
        passthrough_parameters = passthrough_parameter_unpack("amrfinder")
    threads: 1
    resources:
        runtime = "60m",
    conda: "../envs/amrfinder.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.sequence_lengths_sample.{sample}.tsv"
    shell: """
    
        # Collect version number.
        echo "amrfinder $(amrfinder --version)" > "$(dirname {output.table})/.software_version.txt"
        
        amrfinder \
            -p {input.faa:q} \
            --database $(dirname {input.database_representative:q})/latest \
            {params.passthrough_parameters} \
            > {output.table}
        
        {void_report}
        
    """