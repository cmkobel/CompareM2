

rule gapseq:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna"
    output:
        pathways = "{results_directory}/samples/{sample}/gapseq/{sample}_pathways.tsv",
    conda: "../envs/gapseq.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.gapseq_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
        


        
        gapseq \
            doall \
            {input.assembly:q} \
            -K {threads}
        
        
        touch {output} # DEBUG

        {void_report}

    """
    