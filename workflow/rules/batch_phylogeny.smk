rule mashtree:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        fasta = df["input_file_copy"].tolist(),
    output: 
        tree = "{output_directory}/mashtree/mashtree.newick",
        dist = "{output_directory}/mashtree/mash_dist.tsv",
    params:
        passthrough_parameters = passthrough_parameter_unpack("mashtree")
    conda: "../envs/mashtree.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.mashtree.tsv"
    threads: 16
    resources:
        mem_mb = 16000,
    shell: """
    
        # Collect version number.
        mashtree -v > "$(dirname {output.tree})/.software_version.txt"
        
        mashtree \
            --numcpus {threads} \
            --outmatrix {output.dist:q} \
            {params.passthrough_parameters} \
            {input.fasta:q} \
        > {output.tree:q}
        
    """ 

# passthrough_parameter_unpack() misreads keys if this rule is named "mashtree_bootstrap". 
rule bootstrap_mashtree:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        fasta = df["input_file_copy"].tolist(),
    output: 
        tree_boot = "{output_directory}/mashtree/mashtree_bootstrap.newick",
    params:
        passthrough_parameters_bootstrap = passthrough_parameter_unpack("bootstrap_mashtree"),
        passthrough_parameters = passthrough_parameter_unpack("mashtree"),
    conda: "../envs/mashtree.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.mashtree.tsv"
    threads: 16 # sweet spot set on thylakoid. Essentially one thread per physical core.
    resources:
        mem_mb = 16000,
    shell: """
    
        # Collect version number.
        mashtree -v > "$(dirname {output.tree_boot})/.software_version.txt"
        
        mashtree_bootstrap.pl \
            {params.passthrough_parameters_bootstrap} \
            --numcpus {threads} \
            {input.fasta} \
            -- \
            {params.passthrough_parameters} \
        > {output.tree_boot}
        
        # Info:
        #  -- Is used to separate options for mashtree_bootstrap.pl and mashtree
    

        {void_report}
        
    """ 

rule treecluster:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        newick = "{output_directory}/mashtree/mashtree.newick",
    output: 
        table = "{output_directory}/treecluster/treecluster.tsv",
    params:
        passthrough_parameters = passthrough_parameter_unpack("treecluster")
    conda: "../envs/treecluster.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.treecluster.tsv"
    threads: 8
    resources:
        mem_mb = 16000,
    shell: """
    
        # Collect version number.
        TreeCluster.py --version > "$(dirname {output})/.software_version.txt"
        
        TreeCluster.py \
          --input {input.newick} \
          --output {output.table} \
          {params.passthrough_parameters}
          
        {void_report} # Implemented yet?
        
    """ 




rule fasttree:
    input:
        metadata = "{output_directory}/metadata.tsv",
        fasta = "{output_directory}/panaroo/core_gene_alignment_verified.aln",
    output: 
        newick = "{output_directory}/fasttree/fasttree.newick"
    params:
        passthrough_parameters = passthrough_parameter_unpack("fasttree"),
    conda: "../envs/fasttree.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.fasttree.tsv"
    threads: 4
    retries: 2
    resources:
        mem_mb = retry_mem(16000, 32000, 64000, 0),
        runtime = "24h",
    shell: """
    
        # Collect version number. (Is quite hard for fasttree)
        fasttree -expert 2&> .temp_fasttree_version.txt; grep 'Detailed usage' .temp_fasttree_version.txt > "$(dirname {output.newick})/.software_version.txt"; rm .temp_fasttree_version.txt || echo "error deleting .temp_fasttree_version.txt"

        OMP_NUM_THREADS={threads}

        fasttree \
            -nt \
            {params.passthrough_parameters} \
            {input.fasta:q} \
        > {output.newick:q} \
        2> {output.newick:q}.log 

        {void_report}

    """



rule iqtree:
    input:
        metadata = "{output_directory}/metadata.tsv",
        fasta = "{output_directory}/panaroo/core_gene_alignment_verified.aln",
    output: 
        newick = "{output_directory}/iqtree/core_genome_iqtree.treefile"
    params: 
        passthrough_parameters = passthrough_parameter_unpack("iqtree")
    conda: "../envs/iqtree.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.iqtree.tsv"
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
            {params.passthrough_parameters} \
            --prefix $(dirname {output.newick:q})/core_genome_iqtree \
            -redo

        # {void_report} Not in the report yet.

    """
