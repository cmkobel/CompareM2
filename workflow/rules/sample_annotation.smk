
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
    return [
        f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.gff",
        f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.faa",
        f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.log",
        f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.ffn",
        f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.tsv",
        
    ]
        
    

# There is a non-critical bug. When this rule is run after changing the annotator, the inputs from the old annotator are deleted? How does that work? Does snakemake keep track of changes in the job dag and remove old inputs?
# When the annotator choice is changed, the old files are deleted. they shouldn't though, so what I'm thinking is that these files should be linked as individual files instead of as a whole directory. Then the individual links can be deleted and new ones can be laid out. Problem solved, let's get cracking. Update: The problem with this solution is that it doesn't rerun the annotation step when the annotator is changed, but I guess that makes fine sense. Another update: Nu har jeg prøvet at skifte fra prokka til bakta, det ser ud til at den ikke laver nye links? Sandsynligvis fordi annotate slet ikke bliver kørt når man skriver --until bakta. Ja hvis man bare skriver --until annotate. After testing, I can confirm that it works well, but you will have to write --forcerun bakta if you want it to update. Conclusion: Now I've tried both approaches - both having a linked dir and individually linked files. I think a linked dir is cleaner, but the code is not because you need a lot of dirname commands, and you need to remove the (potential old) links in the end of bakta and prokka. The solution I have now, with individually linked files is simpler, so I think I'll stick with it.
rule annotate:
    input: get_annotation_results
    output: # These are mostly the outputs that are used downstream.
        dir = directory("{output_directory}/samples/{sample}/.annotation/"),
        gff = "{output_directory}/samples/{sample}/.annotation/{sample}.gff",
        faa = "{output_directory}/samples/{sample}/.annotation/{sample}.faa", # Used in dbcan, interproscan, diamond_kegg
        log = "{output_directory}/samples/{sample}/.annotation/{sample}.log",
        ffn = "{output_directory}/samples/{sample}/.annotation/{sample}.ffn",
        tsv = "{output_directory}/samples/{sample}/.annotation/{sample}.tsv",
    shell: """

        # Using a softlink means that the choice of annotator can be changed without loss of information. This is especially important when the locustags are "volatile".
        # Though using this solution means that the user can mix up results, so that if the user switches between annotators while running the downstream tools, a mix of annotators can be used. This is very powerful in the sense that the user can use different annotators for different downstream tools, but it is also dangerous in the sense that it can create confusion around which tool is being used. Maybe there should be a log somewhere, or I should go back to having the directory being linked so that the "old" annotation results are removed everytime the annotator is switched.
        
        ln -sr {input} {output.dir}
        
    """


rule prokka:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        assembly = "{output_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{output_directory}/samples/{sample}/prokka/{sample}.gff",
        faa = "{output_directory}/samples/{sample}/prokka/{sample}.faa",
        ffn = "{output_directory}/samples/{sample}/prokka/{sample}.ffn",
        log = "{output_directory}/samples/{sample}/prokka/{sample}.log",
        tsv = "{output_directory}/samples/{sample}/prokka/{sample}.tsv",
        gbk = "{output_directory}/samples/{sample}/prokka/{sample}.gbk",
        gff_nofasta = "{output_directory}/samples/{sample}/prokka/{sample}.gff_nofasta", # Might come in handy.
    params: 
        prokka_rfam = "--rfam" if interpret_true(config['prokka_rfam']) else "", # Set to true (default) or false in config.
        prokka_compliant = "--compliant" if interpret_true(config['prokka_compliant']) else "" # Set to true (default) or false in config.
    conda: "../envs/prokka.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.prokka_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
    
        # Collect version number.
        prokka --version > "$(dirname {output.gff})/.software_version.txt"
        
        prokka \
            --cpus {threads} \
            --force \
            {params.prokka_rfam} \
            {params.prokka_compliant} \
            --outdir {wildcards.output_directory}/samples/{wildcards.sample}/prokka \
            --prefix {wildcards.sample} {input.assembly:q} \
        | tee {output.log:q} 

        # Remove fasta from gff and add sample label
        sed '/^##FASTA/Q' {output.gff:q} > {output.gff_nofasta:q}

        {void_report}

    """
    

rule bakta:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        database_representative = DATABASES + "/bakta/ac2_bakta_database_representative.flag",
        assembly = "{output_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{output_directory}/samples/{sample}/bakta/{sample}.gff",
        faa = "{output_directory}/samples/{sample}/bakta/{sample}.faa",
        tsv = "{output_directory}/samples/{sample}/bakta/{sample}.tsv",
        log = "{output_directory}/samples/{sample}/bakta/{sample}.log",
        ffn = "{output_directory}/samples/{sample}/bakta/{sample}.ffn",
        #gff_generic = "{output_directory}/samples/{sample}/annotation/{sample}.gff3",
    params:
        DATABASES = DATABASES
    conda: "../envs/bakta.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.bakta_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
                        
        # Collect version number.
        bakta --version > "$(dirname {output.gff})/.software_version.txt"
        
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.gff})/.database_version.txt"
        
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
        metadata = "{output_directory}/metadata.tsv",
        database_representative = DATABASES + "/eggnog/ac2_eggnog_database_representative.flag",
        assembly = "{output_directory}/samples/{sample}/{sample}.fna" # Would it be possible to input the .faa from .annotation? That one also uses prodigal.
    output:
        ffn = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.genepred.fasta",
        gff = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.gff",
        hits = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.hits",
        orthologs = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.seed_orthologs",
        tsv = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.annotations",
    #params:
            
    conda: "../envs/eggnog.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.eggnog_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4 # Not sure if the underlying tools are capable of doing lots of parallel computation.
    shell: """
        
        # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#basic-usage
        
        # Collect version number.
        echo "eggnog $(emapper.py --data_dir $(dirname {input.database_representative}) --version) > "$(dirname {output.gff})/.software_version.txt"
            
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.gff})/.database_version.txt"
        
        emapper.py \
            -m diamond \
            --data_dir $(dirname {input.database_representative}) \
            --itype genome \
            --genepred prodigal \
            --override \
            --cpu {threads} \
            --output_dir "$(dirname {output.gff})/" \
            -o "{wildcards.sample}" \
            -i {input.assembly:q} 

        {void_report}

    """
    