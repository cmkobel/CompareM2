rule mashtree:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: 
        tree = "{results_directory}/mashtree/mashtree.newick",
        dist = "{results_directory}/mashtree/mash_dist.tsv",
    conda: "../envs/mashtree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.mashtree.tsv"
    threads: 16
    resources:
        mem_mb = 16000,
    shell: """
    
        # Collect version number.
        mashtree -v > "$(dirname {output.tree})/.software_version.txt"
        
        mashtree \
            --numcpus {threads} \
            --outmatrix {output.dist:q} \
            {input.fasta:q} > {output.tree:q}

        {void_report}
        
    """ 

rule treecluster:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        newick = "{results_directory}/mashtree/mashtree.newick",
    output: 
        table_10 = "{results_directory}/treecluster/treecluster_threshold_0.10.tsv",
    params:
        treecluster_threshold = float(config['treecluster_threshold']) # Default is "0.045"
    conda: "../envs/treecluster.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.treecluster.tsv"
    threads: 8
    resources:
        mem_mb = 16000,
    shell: """
    
        # Collect version number.
        TreeCluster.py --version > "$(dirname {output.table_10})/.software_version.txt"
        
        # First compute with standard 10 (=1-0.90 ANI)
        TreeCluster.py \
          --input {input.newick} \
          --output {output.table_10} \
          --threshold 0.10
          
        # Secondly compute with custom (defaults to 0.045 (= 1-0.955 ANI))
        TreeCluster.py \
          --input {input.newick} \
          --output {wildcards.results_directory}/treecluster/treecluster_threshold_{params.treecluster_threshold}.tsv \
          --method max_clade \
          --threshold {params.treecluster_threshold}
              
        {void_report}
        
    """ 




def get_mem_fasttree(wildcards, attempt): 
    return [16000, 32000, 64000, 0][attempt-1]
rule fasttree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = core_genome_if_exists,
    output: 
        newick = "{results_directory}/fasttree/fasttree.newick"
    conda: "../envs/fasttree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.fasttree.tsv"
    threads: 4
    retries: 2
    resources:
        mem_mb = get_mem_fasttree,
        runtime = "24h",
    shell: """
    
        # Collect version number. (Is quite hard for fasttree)
        fasttree -expert 2&> .temp_fasttree_version.txt; grep 'Detailed usage' .temp_fasttree_version.txt > "$(dirname {output.newick})/.software_version.txt"; rm .temp_fasttree_version.txt || echo "error deleting .temp_fasttree_version.txt"

        OMP_NUM_THREADS={threads}

        fasttree \
            -nt \
            -gtr {input.fasta:q} \
        > {output.newick:q} \
        2> {output.newick:q}.log 

        {void_report}

    """



rule iqtree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = core_genome_if_exists,
    output: 
        newick = "{results_directory}/iqtree/core_genome_iqtree.treefile"
    params: 
        iqtree_bootstraps = int(config['iqtree_bootstraps'])
    conda: "../envs/iqtree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.iqtree.tsv"
    threads: 16
    retries: 3
    resources:
        mem_mb = 32000,
        runtime = "24h",
    shell: """

        # Collect version number.
        iqtree --version | grep version > "$(dirname {output.newick})/.software_version.txt"

        iqtree \
            -s {input.fasta:q} \
            -m GTR \
            --boot {params.iqtree_bootstraps} \
            --prefix $(dirname {output.newick:q})/core_genome_iqtree \
            -redo

        # {void_report} Not in the report yet.

    """
