
from gwf import *



def stem(input):
    """ Routine that finds the stem of a file name (removes extension) """
    stem = '.'.join(str(input).split('.')[:-1])
    return stem

def initialize(title, source_dir, target_dir):
    """ Creates the output/{title} directory"""
    inputs = ""
    outputs = target_dir + "/output/" + title
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '0:02:00', 'queue': 'normal', 'account': 'clinicalmicrobio'}
    spec = f"""
mkdir -p {target_dir}/output/{title}


"""
    return inputs, outputs, options, spec


def copy(source, target_dir, target_file):
    """  Copies the contigs to folders in the output directory """
    inputs = source
    outputs = target_dir + '/' + target_file
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '0:05:00', 'queue': 'normal', 'account': 'clinicalmicrobio'}
    spec = f"""

mkdir -p {target_dir}
cp "{source}" {target_dir}/{target_file}

"""
    return inputs, outputs, options, spec
 

def kraken2(target_dir, title, name):
    inputs = target_dir + '/output/' + title + '/' + name + '/contigs.fa'
    outputs = target_dir + '/output/' + title + '/kraken2/' + name + '_report.txt'
    options = {'nodes': 1, 'cores': 1, 'memory': '8g', 'walltime': '1:00:00', 'queue': 'normal', 'account': 'clinicalmicrobio'}
    spec = f"""
cd {target_dir}/output/{title}

mkdir -p kraken2
cd kraken2



cp ../{name}/contigs.fa {name}.fa

kraken2 --db /project/ClinicalMicrobio/faststorage/database/minikraken2_v2_8GB_201904_UPDATE --report {name}_report.txt {name}.fa > /dev/null 2> /dev/null

rm {name}.fa


    """
    return inputs, outputs, options, spec




def prokka(target_dir, title, name):
    
    #stem = '.'.join(name.split('.')[:-1])

    inputs  = target_dir + '/output/' + title + '/' + name + '/contigs.fa' 
    outputs = target_dir + '/output/' + title + '/' + name + '/prokka/' + name + '.gff'
    options = {'nodes': 1, 'cores': 8, 'memory': '4g', 'walltime': '04:00:00', 'account': 'clinicalmicrobio'} # initially 2 hours
    spec = f"""


cd {target_dir}/output/{title}/{name}


prokka --cpu 8 --outdir prokka --prefix {name} contigs.fa 2> stderr.txt

cp prokka/*.gff annotation.gff

"""
    
    return inputs, outputs , options, spec





def roary(target_dir, title, gffs):
    # target_dir:   Der hvor den skal gemme outputtet.
    # gffs:         En liste med fulde stier til de .gff-filer som skal analyseres.

    hours = -(-len(gffs)//100)*10 # 10 timer for hver 100 filer
    ram = -(-len(gffs)//100)*8 # 8 for hver 100 filer

    inputs = [target_dir + '/output/' + title + '/' + i for i in gffs]
    outputs = [target_dir + '/output/' + title + '/roary/core_gene_alignment.aln'] # Denne fil skal bruges til at lave træet, så det er den vigtigste. Og så også en liste over alle .gff-filer som er brugt.

    newline_for_f_string_workaround = '\n'
    options = {'nodes': 1, 'cores': 16, 'memory': f'{ram}g', 'walltime': f'{hours}:00:00', 'queue': 'normal', 'account': 'ClinicalMicrobio'}
    spec = f'''


cd {target_dir}/output/{title}

roary -f roary -e -v -r -p 16 {' '.join(gffs)} 2> roary_stderr.txt

#mail -s 'roary done {title}' kobel@pm.me <<< 'Sent from the compare pipeline'"""

echo JOBID $SLURM_JOBID
jobinfo $SLURM_JOBID
'''
    return (inputs, outputs, options, spec)


def fasttree(target_dir, title):
    inputs = target_dir + '/output/' + title + '/roary/core_gene_alignment.aln'
    outputs = target_dir + '/output/' + title + '/fasttree/tree.newick'
    options = {'nodes': 1, 'cores': 8, 'memory': f'8g', 'walltime': f'2:00:00', 'queue': 'normal', 'account': 'ClinicalMicrobio'}
    spec = f"""
cd {target_dir}/output/{title}
mkdir fasttree
cd fasttree



FastTree -nt -gtr ../roary/core_gene_alignment.aln > tree.newick 2> stderr.txt


"""
    return inputs, outputs, options, spec


def abricate():
    # Depends on core gene alignment.
    pass

def quicktree(target_dir, title, names):
    pass


def roary_plots():
    pass

def mailzip():
    spec = """mail -s 'clinmicpipe done {group_name}' -a pangenome_matrix.png kobel@pm.me <<< 'Sent from workflow_templates.py'"""



