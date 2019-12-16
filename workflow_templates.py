
from gwf import *

DEBUG_STATUS = False

def debug(title = ''):
	if DEBUG_STATUS:
		return '2> ' + str(title).strip() + '_stderr.txt'
	else:
		return ""




def stem(input):
    """ Routine that finds the stem of a file name (removes extension) """
    stem = '.'.join(str(input).split('.')[:-1])
    return stem



def initialize(title, source_dir, target_dir):
    """ Creates the output/{title} directory"""
    inputs = ""#source_dir
    outputs = target_dir + "/output/" + title #+ "/initialized.txt"
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '0:02:00',  'account': 'clinicalmicrobio'}
    spec = f"""

mkdir -p {target_dir}/output/{title}
cd {target_dir}/output/{title}

#echo "started3" >> initialized.txt {debug('init')}


"""
    return inputs, outputs, options, spec


def copy(source, target_dir, title, name):
    """  Copies the contigs to folders in the output directory and converts everything to fasta"""
    inputs = source
    outputs = target_dir + '/output/' + title + '/' + name + '/' + 'contigs.fa'
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '0:05:00',  'account': 'clinicalmicrobio'}
    spec = f"""

mkdir -p {target_dir + '/output/' + title + '/' + name}
any2fasta "{source}" > {target_dir + '/output/' + title + '/' + name}/'contigs.fa'

"""
    return inputs, outputs, options, spec
 

def kraken2(target_dir, title, name):
    inputs = target_dir + '/output/' + title + '/' + name + '/contigs.fa'
    outputs = target_dir + '/output/' + title + '/kraken2/' + name + '_report.txt'
    options = {'nodes': 1, 'cores': 1, 'memory': '8g', 'walltime': '1:00:00', 'account': 'clinicalmicrobio'}
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


prokka --cpu 8 --force --outdir prokka --prefix {name} contigs.fa {debug('prokka')}

cp prokka/*.gff annotation.gff

"""
    
    return inputs, outputs , options, spec


# MLST: Multi Locus Sequence Typing
def mlst(target_dir, title, contigs):
    inputs = [target_dir + '/output/' + title + '/' + i for i in contigs]
    outputs = target_dir + '/output/' + title + '/mlst.tsv' # Denne fil skal bruges til at lave træet, så det er den vigtigste. Og så også en liste over alle .gff-filer som er brugt.
    
    options = {'nodes': 1, 'cores': 2, 'memory': '4g', 'walltime': '01:00:00',  'account': 'clinicalmicrobio'}
    spec = f'''


cd {target_dir}/output/{title}
#mkdir mlst


mlst {' '.join(contigs)} > mlst.tsv {debug('mlst')}

'''
    return (inputs, outputs, options, spec)	



# Roary: The pan genome pipeline
def roary(target_dir, title, gffs):
    # target_dir:   Der hvor den skal gemme outputtet.
    # gffs:         En liste med fulde stier til de .gff-filer som skal analyseres.

    hours = -(-len(gffs)//100)*10 # 10 timer for hver 100 filer
    ram = -(-len(gffs)//100)*8 # 8 for hver 100 filer

    inputs = [target_dir + '/output/' + title + '/' + i for i in gffs]
    outputs = [target_dir + '/output/' + title + '/roary/core_gene_alignment.aln', # Denne fil skal bruges til at lave træet, så det er den vigtigste. Og så også en liste over alle .gff-filer som er brugt.
               target_dir + '/output/' + title + '/roary/gene_presence_absence.csv']
    newline_for_f_string_workaround = '\n'
    options = {'nodes': 1, 'cores': 16, 'memory': f'{ram}g', 'walltime': f'{hours}:00:00', 'account': 'ClinicalMicrobio'}
    spec = f'''


cd {target_dir}/output/{title}

roary -f roary -e -v -r -p 16 {' '.join(gffs)} 2> roary_stderr.txt



echo JOBID $SLURM_JOBID
jobinfo $SLURM_JOBID
'''
    return (inputs, outputs, options, spec)


def fasttree(target_dir, title):
    # todo: time and mem should depend on number of isolates
    inputs = target_dir + '/output/' + title + '/roary/core_gene_alignment.aln'
    outputs = [target_dir + '/output/' + title + '/fasttree/tree.newick',
               target_dir + '/output/' + title + '/fasttree/tree.pdf']
    options = {'nodes': 1, 'cores': 8, 'memory': '8g', 'walltime': '6:00:00', 'account': 'clinicalmicrobio'}
    spec = f"""
cd {target_dir}/output/{title}
mkdir -p fasttree
cd fasttree




FastTree -nt -gtr ../roary/core_gene_alignment.aln > tree.newick {debug('ft')}



Rscript /project/ClinicalMicrobio/faststorage/compare/scripts/R/ape_newick2pdf.r tree.newick "{title} core genome"  {debug('R')}


"""
    return inputs, outputs, options, spec


def abricate():
    # Depends on core gene alignment.
    pass

def quicktree(target_dir, title, names):
    pass



def roary_plots(target_dir, title):
    script_file = '/faststorage/project/ClinicalMicrobio/compare/scripts/py/roary_plots.py'
    #tree_file = 

    inputs = [target_dir + '/output/' + title + '/fasttree/tree.newick',
              target_dir + '/output/' + title + '/roary/gene_presence_absence.csv']
    outputs = [target_dir + '/output/' + title + '/roary_plots/pangenome_matrix.png',
               target_dir + '/output/' + title + '/roary_plots/pangenome_matrix_alternative.svg',
               target_dir + '/output/' + title + '/roary_plots/pangenome_matrix_alternative.pdf']
    #
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '00:10:00', 'account': 'clinicalmicrobio'}
    spec = f'''
cd {target_dir}/output/{title}
mkdir -p roary_plots
cd roary_plots

python {script_file} ../fasttree/tree.newick ../roary/gene_presence_absence.csv 2> stderr.out

perl {target_dir}/scripts/perl/roary2svg.pl ../roary/gene_presence_absence.csv > pangenome_matrix_alternative.svg 2> 2_stderr.out
cairosvg pangenome_matrix_alternative.svg -o pangenome_matrix_alternative.pdf


#mail -s "compare done {title}" -a pangenome_matrix.png $COMPARE_DEFAULT_EMAIL_ADDRESS <<< "Sent from the compare pipeline" &
# find out how to attach multiple files
'''
    return inputs, outputs, options, spec



def panito(target_dir, title):
    inputs = [target_dir + '/output/' + title + '/roary/core_gene_alignment.aln',
              target_dir + '/output/' + title + '/fasttree/tree.newick']
    outputs = target_dir + '/output/' + title + '/panito/ani.pdf'
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '00:10:00', 'account': 'clinicalmicrobio'}
    spec = f'''

cd {target_dir}/output/{title}
mkdir -p panito
cd panito

panito {inputs[0]} > ani.tsv {debug('panito')}
Rscript /project/ClinicalMicrobio/faststorage/compare/scripts/R/aniplot.r ani.tsv {inputs[1]} {debug('r_panito')}


'''
    return inputs, outputs, options, spec




def send_mail(target_dir, title, names):
    inputs = [target_dir + '/output/' + title + '/roary/summary_statistics.txt',
              target_dir + '/output/' + title + '/panito/ani.pdf',
              target_dir + '/output/' + title + '/roary_plots/pangenome_matrix.png',
              target_dir + '/output/' + title + '/roary_plots/pangenome_matrix_alternative.pdf',
              target_dir + '/output/' + title + '/fasttree/tree.newick',
              target_dir + '/output/' + title + '/fasttree/tree.pdf',
              target_dir + '/output/' + title + '/mlst.tsv'
              
                     
        ]

    newline = '\n'
    outputs = target_dir + '/output/' + title + '/mailsente' # If it doesn't have an arbitrary output, the first job (init) will be run
    
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '00:10:00', 'account': 'clinicalmicrobio'}
    spec = f"""


cd {target_dir}/output/{title}
# touch mailsent
# collect mail content

echo -e "Assembly Comparator results for {title}\n" > mail.txt

echo -e "list of assemblies:" >> mail.txt
echo -e "{str(newline).join(names)}" >> mail.txt

echo -e "\n" >> mail.txt

echo "Roary summary statistics:" >> mail.txt
cat roary/summary_statistics.txt >> mail.txt

echo -e "\n" >> mail.txt

echo -e "MLST results:" >> mail.txt
sed 's/\<contigs.fa\>//g' mlst.tsv >> mail.txt

echo -e "\n" >> mail.txt

echo -e "A few small output files from the pipeline has been attached in the zip-file" >> mail.txt

echo -e "To access the full analysis, please visit /project/ClinicalMicrobio/faststorage/compare/{title} on GenomeDK." >> mail.txt



zip -j {title}.zip {' '.join(inputs)}

mail -s "compare done: {title}" -a {title}.zip -q mail.txt kobel@pm.me <<< "" 
#mail -s "compare done: {title}" {' '.join(['-a ' + i for i in inputs])} -q mail.txt kobel@pm.me <<< "" 

rm {title}.zip



#mail -s "evolve compare" {' '.join(['-a ' + i for i in inputs])} kobel@pm.me <<< "please work, from the pipeline" {debug('mail_1')}

#mail -s "comparator done: {title}" {' '.join(['-a ' + i for i in inputs])} $COMPARE_DEFAULT_EMAIL_ADDRESS <<< "To access the full analysis, please visit /project/ClinicalMicrobio/faststorage/compare/{title} on GenomeDK." {debug('mail_2')}
#touch mail_sent_{title}


    """
    return inputs, outputs, options, spec




