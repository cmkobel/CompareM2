
    
def core_genome_if_exists(wildcards): 
    
    summary_file = checkpoints.roary.get(**wildcards).output["summary"]
    alignment_file = checkpoints.roary.get(**wildcards).output["alignment"]
    
    try:   
        with summary_file.open() as f:
        #with open(file) as f: 

            for line in f:
                line = line.strip()
                #print(f"line: {line}")
                
                # Find the line that states the size of the core genome.
                if "Core genes" in line:
                    #print(f"Pipeline: Parsing Roary core genome from line \"{line}\"")
                    splitted = line.split("\t")
                    #print(splitted)
                    core_genome_size = int(splitted[-1])
                    
                    #print(f"Pipeline: Parsed core genome size is: {repr(core_genome_size)}")
                    break
                
            # Return true when there is a core genome.
            if core_genome_size > 0:
                print(f"Pipeline: Core genome exists ({core_genome_size} genes).")
                return [alignment_file]
                
            else:
                print(f"Pipeline: Core genome is empty ({core_genome_size} genes).")
                return list() # Empty list
                
                
    # In case the file cannot be parsed, we're assuming that there is no core genome. Hence false.
    except Exception as e: 
        print(f"Pipeline Error: {e}") # Show errors in parsing the roary summary, but don't fail.
        print("Pipeline: Core genome size not known?")
        return list() # Empty list
def get_mem_roary(wildcards, attempt): 
    return [32000, 64000, 128000][attempt-1]


checkpoint roary: # Checkpoint, because some rules i.e. fasttree, iqtree, snp-dists should only run if the core genome is non-empty.
    input: 
        metadata = "{results_directory}/metadata.tsv",
        gff = expand("{results_directory}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], results_directory = results_directory),
    output:
        summary = "{results_directory}/roary/summary_statistics.txt",
        presence = "{results_directory}/roary/gene_presence_absence.csv",
        alignment = "{results_directory}/roary/core_gene_alignment.aln",
        #analyses = ["{results_directory}/roary/summary_statistics.txt", "{results_directory}/roary/core_gene_alignment.aln", "{results_directory}/roary/gene_presence_absence.csv"]
    params:
        blastp_identity = int(config['roary_blastp_identity']), # = 95 # For clustering genes
        core_perc = 99,  # Definition of the core genome
        group_limit = 100000, # Default group limit is 50000 # TODO put this and any other program argument into the config file like blastp_identity
    benchmark: "{results_directory}/benchmarks/benchmark.roary.tsv"
    threads: 16
    #retries: 2
    resources:
        #mem_mb = 32768,
        mem_mb = get_mem_roary,
        runtime = "24h",
    conda: "../envs/roary_see-comments-in-this-file.yaml" # deactivated 
    shell: """
    
        # Since I reinstalled conda, I've had problems with "Can't locate Bio/Roary/CommandLine/Roary.pm in INC". Below is a hacky fix
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0
        
        # Silence
        echo "will cite" | parallel --citation > /dev/null 2> /dev/null; echo "Please cite GNU Parallel, Ole Tange and Free Software Foundation."

        # Roary is confused by the way snakemake creates directories ahead of time.
        rm -r {wildcards.results_directory}/roary


        roary -a -r -e --mafft \
            -p {threads} \
            -i {params.blastp_identity} \
            -cd {params.core_perc} \
            -f {wildcards.results_directory}/roary \
            --group_limit {params.group_limit} \
            {input.gff:q}

            
        {void_report}

    """
