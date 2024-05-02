
def get_annotation_results(wildcards):
    """ 
        Reads the config annotator parameter and returns the necessary 
        annotation results for that annotator. This could have also been an 
        anonymous rule but I think this is more snakemakistic. 
        This could also have been an action in the end of the annotaton rules (prokka and bakta), but then it would become complicated when the user wanted to run both, to possibly compare the output.
    """
        
    raw = config["annotator"]
    parsed = str(raw).lower().strip()
    
    extra_message = "" # Extra message that can be added for debugging etc.
    
    if parsed == "bakta":
        annotator = "bakta"
    elif parsed == "prokka":
        annotator = "prokka"
    else: # base case, fall back to prokka.
        extra_message += " (default)"
        annotator = "prokka"
        
        
    print(f"Pipeline: Using the {annotator} annotator for sample \"{wildcards.sample}\"{extra_message}.")
    
    # Return should only contain one item, as I can't name them and I need to refer to a single one in rule annotate where I'm accessing its dirname().
    return f"{wildcards.results_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.gff"
    

# There is a non-critical bug. When this rule is run after changing the annotator, the inputs from the old annotator are deleted? How does that work? Does snakemake keep track of changes in the job dag and remove old inputs?
rule annotate:
    input: get_annotation_results
    output:
        gff = "{results_directory}/samples/{sample}/.annotation/{sample}.gff",
        faa = "{results_directory}/samples/{sample}/.annotation/{sample}.faa", # Used in dbcan, interproscan, diamond_kegg
        log = "{results_directory}/samples/{sample}/.annotation/{sample}.log",
        ffn = "{results_directory}/samples/{sample}/.annotation/{sample}.ffn",
        tsv = "{results_directory}/samples/{sample}/.annotation/{sample}.tsv",
        
    params:
        annotator = config['annotator']
    shell: """
        
        # Clean up / make ready.
        test -d $(dirname {output.gff}) && rm -rv $(dirname {output.gff}) # Remove potential snakemake created dir.
        test -h $(dirname {output.gff}) && rm -v $(dirname {output.gff}) # Remove potential old link.

        # Using a softlink means that the choice of annotator can be changed without loss of information. This is especially important when the locustags are "volatile".
        ln -sr $(dirname {input}) $(dirname {output.gff})
    
    """
    

rule prokka:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{results_directory}/samples/{sample}/prokka/{sample}.gff",
        faa = "{results_directory}/samples/{sample}/prokka/{sample}.faa",
        #ffn = "{results_directory}/samples/{sample}/prokka/{sample}.ffn",
        log = "{results_directory}/samples/{sample}/prokka/{sample}.log",
        tsv = "{results_directory}/samples/{sample}/prokka/{sample}.tsv",
        gff_nofasta = "{results_directory}/samples/{sample}/prokka/{sample}.gff_nofasta", # Might come in handy.
    conda: "../envs/prokka.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.prokka_individual.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
        
        prokka \
            --cpus {threads} \
            --force \
            --rfam \
            --compliant \
            --outdir {wildcards.results_directory}/samples/{wildcards.sample}/prokka \
            --prefix {wildcards.sample} {input.assembly:q} \
        | tee {output.log:q} 

        # Remove fasta from gff and add sample label
        gff_fasta_start=$(grep --line-number --extended-regexp "^##FASTA" {output.gff:q} | cut -f1 -d:)
        head --lines $(($gff_fasta_start-1)) {output.gff:q} \
        > {output.gff_nofasta:q}

        {void_report}

    """
    

rule bakta:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/bakta/ac2_bakta_database_representative.flag",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{results_directory}/samples/{sample}/bakta/{sample}.gff",
        faa = "{results_directory}/samples/{sample}/bakta/{sample}.faa",
        tsv = "{results_directory}/samples/{sample}/bakta/{sample}.tsv",
        
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
            
        # Optimize compatibility between prokka and bakta to make fluent use of either result easy.
        cp {output.gff}3 {output.gff}
        
        {void_report}

    """



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
    benchmark: "{results_directory}/benchmarks/benchmark.eggnog_individual.{sample}.tsv"
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
    