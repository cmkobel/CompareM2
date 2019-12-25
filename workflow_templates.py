
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


def kraken2_table(target_dir, title, names):
    inputs = [target_dir + '/output/' + title + '/kraken2/' + name + '_report.txt' for name in names]
    outputs = target_dir + '/output/' + title + '/kraken2-table.txt'
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '00:10:00', 'account': 'clinicalmicrobio'}

    command = '''for f in *_report.txt; do echo ${f::-11} >> ../kraken2-table.txt; cat $f | awk '$4 ~ "^S$" {printf("%6.2f%% %s %s %s %s\\n", $1, $6, $7, $8, $9)}' | head -n 3 >> ../kraken2-table.txt; echo >> ../kraken2-table.txt; done'''
    
    spec = f"""
cd {target_dir}/output/{title}/kraken2

{command}

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
def roary(target_dir, title, gffs, blastp_identity = 95, allow_paralogs = False):
    # target_dir:   Der hvor den skal gemme outputtet.
    # gffs:         En liste med fulde stier til de .gff-filer som skal analyseres.

    hours = -(-len(gffs)//100)*10 # 10 timer for hver 100 filer
    ram = -(-len(gffs)//100)*8 # 8 for hver 100 filer

    if allow_paralogs == True:
        ap_string = "-ap"
    else:
        ap_string = ""

    inputs = [target_dir + '/output/' + title + '/' + i for i in gffs]
    outputs = [target_dir + '/output/' + title + '/roary/core_gene_alignment.aln', # Denne fil skal bruges til at lave træet, så det er den vigtigste. Og så også en liste over alle .gff-filer som er brugt.
               target_dir + '/output/' + title + '/roary/gene_presence_absence.csv',
               target_dir + '/output/' + title + '/roary/Rplots.pdf']
    newline_for_f_string_workaround = '\n'
    options = {'nodes': 1, 'cores': 16, 'memory': f'{ram}g', 'walltime': f'{hours}:00:00', 'account': 'ClinicalMicrobio'}
    spec = f'''


cd {target_dir}/output/{title}

echo "blastp_identity (-i) = {str(blastp_identity)}" > roary_thresholds.txt
echo "allow paralogs (-ap): {allow_paralogs}" >> roary_thresholds.txt

roary -f roary -e -v -r -p 16 -i {int(blastp_identity)} {ap_string} {' '.join(gffs)} {debug('roary')}



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

touch in_roary_plots

python {script_file} ../fasttree/tree.newick ../roary/gene_presence_absence.csv 2> stderr.out

perl {target_dir}/scripts/perl/roary2svg.pl ../roary/gene_presence_absence.csv > pangenome_matrix_alternative.svg 2> 2_stderr.out
cairosvg pangenome_matrix_alternative.svg -o pangenome_matrix_alternative.pdf 2> 3_stderr.out

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
              target_dir + '/output/' + title + '/mlst.tsv',
              target_dir + '/output/' + title + '/roary/Rplots.pdf',
              target_dir + '/output/' + title + '/kraken2-table.txt']

    newline = '\n'
    outputs = target_dir + '/output/' + title + '/mailsente' # If it doesn't have an arbitrary output, the first job (init) will be run
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '00:10:00', 'account': 'clinicalmicrobio'}
    awk_command = """awk '{printf("%s %s (%s)\\n", $1, $2, $3)}'"""
    spec = f"""


cd {target_dir}/output/{title}

# collect mail content


echo -e "Assembly Comparator results for {title}" > mail.txt

echo -e "\n" >> mail.txt

echo -e "MLST results:" >> mail.txt
cat mlst.tsv | sed 's/\/contigs.fa//g' | {awk_command} | column -t >> mail.txt

echo -e "\n" >> mail.txt

echo "Roary summary statistics:" >> mail.txt
cat roary/summary_statistics.txt | column -ts $'\t' >> mail.txt

echo -e "\n" >> mail.txt

echo -e "Top 3 kraken results:" >> mail.txt

cat kraken2-table.txt >> mail.txt


echo -e "\n" >> mail.txt


echo -e "A few small output files from the pipeline has been attached in the zip-file" >> mail.txt
echo -e "To access the full analysis, please visit /project/ClinicalMicrobio/faststorage/compare/output/{title} on GenomeDK." >> mail.txt




zip -j {title}.zip {' '.join(inputs)}

mailx -s "comparator done: {title}" kobel@pm.me <<< ""
mailx -s "[comparator] done: {title}" -a {title}.zip -q mail.txt $COMPARATOR_DEFAULT_EMAIL_ADDRESS <<< "" 

rm {title}.zip

echo $(date) >> mailsent



    """
    return inputs, outputs, options, spec




