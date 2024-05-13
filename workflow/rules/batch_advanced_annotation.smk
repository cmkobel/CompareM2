
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
    
        # Collect version number.
        R -s -q -e "library(clusterProfiler); sessionInfo()"  | grep -P "R version|clusterProfiler" > "$(dirname {output.diamond}).software_version.txt"


        Rscript {params.script:q} \
            {input.kegg_asset:q} \
            {params.output_dir:q} \
            {input.kegg_diamond:q}  
            
        {void_report}

    """



def get_mem_gtdbtk(wildcards, attempt): 
    return [150000, 300000, 400000, 500000][attempt-1]



# The rule is called by the software, whereas the results are called by the database. Is that confusing?
rule gtdbtk:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/gtdb/ac2_gtdb_database_representative.flag",
        fasta = df["input_file_fasta"].tolist(),
    output: 
        tsv = "{results_directory}/gtdbtk/gtdbtk.summary.tsv"
    params:
        batchfile_content = df[['input_file_fasta', 'sample']].to_csv(header = False, index = False, sep = "\t"),
        out_dir = "{results_directory}/gtdbtk/",
        base_variable = base_variable, # not used?
        mash_db = f"{DATABASES}/gtdb_sketch_release220/mash_db.msh",
    threads: 8
    #retries: 3
    resources:
        # mem_mb = get_mem_gtdbtk, # works great, but I want to use a limited amount to test a potential hardware upgrade
        mem_mb = 131072,
        runtime = "48h"
    conda: "../envs/gtdbtk.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.gtdbtk.tsv"
    shell: """
    
        # Collect version number.
        gtdbtk -v > "$(dirname {output.tsv}).software_version.txt"
        
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.tsv})/.database_version.txt"

        # TODO: Using skip-ani-screen is not optimal, as it possibly speeds up a lot.
        mkdir -p $(dirname {params.mash_db:q})

        # I need to find a neat way of setting these variables. Maybe the user has an older/newer version than what is hardcoded here. 
        export GTDBTK_DATA_PATH="$(dirname {input.database_representative:q})/release220/" # Should be defined from config file, and not be hardwired.
        
        # Create batchfile
        echo '''{params.batchfile_content}''' > {wildcards.results_directory}/gtdbtk/batchfile.tsv
        
        gtdbtk classify_wf \
            --mash_db {params.mash_db:q} \
            --batchfile {wildcards.results_directory}/gtdbtk/batchfile.tsv \
            --out_dir {params.out_dir:q} \
            --cpus {threads} \
            --keep_intermediates \
            --force

        # Homogenize database version number
        #cp {wildcards.results_directory}/gtdbtk/gtdbtk.bac120.summary.tsv {output:q}
        # New better version below that also incorporates archaea
        #cat {wildcards.results_directory}/gtdbtk/gtdbtk.*.summary.tsv {output:q}
        

        # Even better: Should be tested on originals
        echo -e "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings" \
        > {output:q}
        tail --quiet -n +2 {wildcards.results_directory}/gtdbtk/gtdbtk.*.summary.tsv \
        >> {output:q}
        

        {void_report}

    """ 

