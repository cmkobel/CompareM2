
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
 

def prokka(target_dir, title, name):
    
    #stem = '.'.join(name.split('.')[:-1])

    inputs  = target_dir + '/output/' + title + '/' + name + '/contigs.fa' 
    outputs = target_dir + '/output/' + title + '/' + name + '/prokka'
    options = {'nodes': 1, 'cores': 8, 'memory': '4g', 'walltime': '04:00:00', 'account': 'clinicalmicrobio'} # initially 2 hours
    spec = f"""


cd {target_dir}/output/{title}/{name}
touch canttouchthis

prokka --cpu 8 --outdir prokka contigs.fa > stdout.txt 2> stderr.txt

"""
    
    return inputs, outputs , options, spec



def roary(dir, gffs):
    # dir:      Der hvor den skal gemme outputtet.
    # gffs:     En liste med fulde stier til de .gff-filer som skal analyseres.

    hours = -(-len(gffs)//100)*10 # 10 timer for hver 100 filer
    ram = -(-len(gffs)//100)*8 # 8 for hver 100 filer

    inputs = gffs
    outputs = [dir + 'core_gene_alignment.aln', dir + 'gffs.txt'] # Denne fil skal bruges til at lave træet, så det er den vigtigste. Og så også en liste over alle .gff-filer som er brugt.

    newline_for_f_string_workaround = '\n'
    options = {'nodes': 1, 'cores': 8, 'memory': f'{ram}g', 'walltime': f'{hours}:00:00', 'queue': 'normal', 'account': 'ClinicalMicrobio'}
    spec = f'''
roary -f {dir} -e -v -p 8 {' '.join(gffs)} # Easier to spot number of gffs with a newline.

echo "{newline_for_f_string_workaround.join(gffs)}" > {dir}gffs.txt
echo JOBID $SLURM_JOBID
jobinfo $SLURM_JOBID
'''
    #print('\n\n\n\nthisisthespec\n\n', spec) # for debug only.
    return (inputs, outputs, options, spec)