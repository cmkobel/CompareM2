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
        mashtree -v > "$(dirname {output})/.software_version.txt"
        
        mashtree \
            --numcpus {threads} \
            --outmatrix {output.dist:q} \
            {input.fasta:q} > {output.tree:q}

        {void_report}
        
    """ 



def get_mem_fasttree(wildcards, attempt): 
    return [16000, 32000, 64000, 0][attempt-1]
rule fasttree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = core_genome_if_exists,
    output: "{results_directory}/fasttree/fasttree.newick"
    conda: "../envs/fasttree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.fasttree.tsv"
    threads: 4
    retries: 2
    resources:
        mem_mb = get_mem_fasttree,
        runtime = "24h",
    shell: """
    
        # Collect version number. (Is quite hard for fasttree)
        fasttree -expert 2&> .temp_fasttree_version.txt; grep 'Detailed usage' .temp_fasttree_version.txt > "$(dirname {output})/.software_version.txt"; rm .temp_fasttree_version.txt || echo "error deleting .temp_fasttree_version.txt"

        OMP_NUM_THREADS={threads}

        fasttree \
            -nt \
            -gtr {input.fasta:q} \
        > {output:q} \
        2> {output:q}.log 

        {void_report}

    """



rule iqtree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = core_genome_if_exists,
    output: 
        newick = "{results_directory}/iqtree/core_genome_iqtree.treefile"
    conda: "../envs/iqtree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.iqtree.tsv"
    threads: 16
    retries: 3
    resources:
        mem_mb = 32000,
        runtime = "24h",
    shell: """

        # Collect version number.
        iqtree --version > "$(dirname {output})/.software_version.txt"

        iqtree \
            -s {input.fasta:q} \
            -m GTR \
            --boot 100 \
            --prefix $(dirname {output.newick:q})/core_genome_iqtree \
            -redo

        # {void_report} Not in the report yet.


    """
