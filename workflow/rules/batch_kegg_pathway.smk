
rule kegg_pathway:
    input: 
        kegg_asset = base_variable + "/resources/ko00001.json", # Downloaded from kegg.jp
        kegg_diamond = expand("{results_directory}/samples/{sample}/diamond_kegg/{sample}_diamond_kegg.tsv", sample = df["sample"], results_directory = results_directory)
    output: 
        diamond = "{results_directory}/kegg_pathway/kegg_pathway_enrichment_analysis.tsv"
    params:
        output_dir = "{results_directory}/kegg_pathway",
        script = base_variable + "/workflow/scripts/kegg_pathway_enrichment_analysis.R"
    conda: "../envs/r-clusterProfiler.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.kegg_pathway.tsv"
    shell: """

        Rscript {params.script:q} \
            {input.kegg_asset:q} \
            {params.output_dir:q} \
            {input.kegg_diamond:q}  
            
        {void_report}

    """
