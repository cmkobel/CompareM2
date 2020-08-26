
from gwf import *
import os


BLASTP = int(os.environ["BLASTP"])

DEBUG_STATUS = False

def debug(title = ''):
	if DEBUG_STATUS:
		return '2> ' + str(title).strip() + '_stderr.txt'
	else:
		return ""


environment = """\
source /faststorage/project/ClinicalMicrobio/compare/environment/etc/profile.d/conda.sh
conda activate comparator"""




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

        {environment}

        mkdir -p {target_dir}/output/{title}
        cd {target_dir}/output/{title}

        #echo "started3" >> initialized.txt {debug('init')}

    """
    
    return inputs, outputs, options, spec


def copy(source, target_dir, title, name):
    """  Copies the contigs to folders in the output directory and converts everything to fasta"""
    inputs = source
    outputs = target_dir + '/output/' + title + '/' + name + '/' + 'contigs.fa'
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '1:00:00',  'account': 'clinicalmicrobio'}
    spec = f"""

        {environment}

        mkdir -p {target_dir + '/output/' + title + '/' + name}
        any2fasta "{source}" | /project/ClinicalMicrobio/faststorage/compare/scripts/py/fasta_shorten_headers.py > {target_dir + '/output/' + title + '/' + name}/contigs.fa


    """

    return inputs, outputs, options, spec
 

def kraken2(target_dir, title, name):
    inputs = target_dir + '/output/' + title + '/' + name + '/contigs.fa'
    outputs = target_dir + '/output/' + title + '/kraken2/' + name + '_report.txt'
    options = {'nodes': 1, 'cores': 1, 'memory': '8g', 'walltime': '1:00:00', 'account': 'clinicalmicrobio'}
    spec = f"""

        {environment}

        cd {target_dir}/output/{title}

        mkdir -p kraken2
        cd kraken2

        cp ../{name}/contigs.fa {name}.fa

        kraken2 --db /project/ClinicalMicrobio/faststorage/database/minikraken_8GB_20200312 --report {name}_report.txt {name}.fa > /dev/null 2> /dev/null

        rm {name}.fa


    """
    return inputs, outputs, options, spec





def summary_kraken(target_dir, title, names):
    #inputs = [target_dir + '/output/' + title + '/kraken2/' + name + '_report.txt' for name in names]
    inputs = []
    for name in names:
        inputs.append(target_dir + '/output/' + title + '/kraken2/' + name + '_report.txt')
        


    outputs = [target_dir + '/output/' + title + '/kraken2-table.txt']
                
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '00:10:00', 'account': 'clinicalmicrobio'}

    command = '''for f in *_report.txt; do echo ${f::-11} >> ../kraken2-table.txt; cat $f | awk -F '\\t' '$4 ~ "(^S$)|(U)" {gsub(/^[ \\t]+/, "", $6); printf("%6.2f%% %s\\n", $1, $6)}' | sort -gr  | head -n 3 >> ../kraken2-table.txt; echo >> ../kraken2-table.txt; done'''
    #kraken_reads_top_command = """awk '$4 ~ "^S$" {printf("%05.2f\\t%s %s %s %s\\n", $1, $6, $7, $8, $9)}'"""
    
    spec = f"""

        {environment}

        cd {target_dir}/output/{title}/kraken2

        # delete possibly old file.
        touch ../kraken2-table.txt
        rm ../kraken2-table.txt

        {command}

    """

    return inputs, outputs, options, spec

def summary_abricate(target_dir, title, names):
    #inputs = [target_dir + '/output/' + title + '/kraken2/' + name + '_report.txt' for name in names]
    inputs = []
    for name in names:
        inputs.append(target_dir + '/output/' + title + '/abricate/isolates/ncbi_' + name + '.tab')


    outputs = [target_dir + '/output/' + title + '/abricate/abricate_resfinder_summary.tab',
               target_dir + '/output/' + title + '/abricate/abricate_vfdb_summary.tab']
                
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '01:00:00', 'account': 'clinicalmicrobio'}

    
    spec = f"""

        {environment}

        cd {target_dir}/output/{title}/abricate

        #abricate --nopath *.tab --summary > ../amr_virulence_summary.tab

        abricate --nopath isolates/plasmidfinder_*.tab --summary > abricate_plasmidfinder_summary.tab
        abricate --nopath isolates/ecoli_vf_*.tab --summary > abricate_ecoli_vf_summary.tab
        abricate --nopath isolates/ncbi_*.tab --summary > abricate_ncbi_summary.tab
        abricate --nopath isolates/resfinder_*.tab --summary > abricate_resfinder_summary.tab
        abricate --nopath isolates/argannot_*.tab --summary > abricate_argannot_summary.tab
        abricate --nopath isolates/vfdb_*.tab --summary > abricate_vfdb_summary.tab
        abricate --nopath isolates/megares_*.tab --summary > abricate_megares_summary.tab
        abricate --nopath isolates/ecoh_*.tab --summary > abricate_ecoh_summary.tab
        abricate --nopath isolates/card_*.tab --summary > abricate_card_summary.tab

        #cp abricate_ncbi_summary.tab ../amr_virulence_summary.tab



        # #collect all abricate results
        # List=$(ls *.tab)
        # arr=($List)
        # first=${{arr[1]}}
        # 
        # cat $first | grep -E "^#" > abricate_all.tsv
        # cat *.tab | grep -vE "^#" >> abricate_all.tsv


    """

    return inputs, outputs, options, spec
    


def abricate(target_dir, title, name):
    inputs = target_dir + '/output/' + title + '/' + name + '/contigs.fa'
    outputs = [target_dir + '/output/' + title + '/abricate/isolates/ncbi_' + name + '.tab']
               #target_dir + '/output/' + title + '/abricate/' + name + '.tab']
    options = {'nodes': 1, 'cores': 1, 'memory': '8g', 'walltime': '1:00:00', 'account': 'clinicalmicrobio'}
    spec = f"""

        {environment}


        cd {target_dir}/output/{title}



        mkdir -p abricate/isolates
        cd abricate/isolates
        cp ../../{name}/contigs.fa {name}.fa




        echo "starting plasmidfinder abrication"
        abricate --db plasmidfinder {name}.fa > plasmidfinder_{name}.tab

        echo "starting ecoli_vf abrication"
        abricate --db ecoli_vf {name}.fa > ecoli_vf_{name}.tab

        echo "starting ncbi abrication"
        abricate --db ncbi {name}.fa > ncbi_{name}.tab

        echo "starting resfinder abrication"
        abricate --db resfinder {name}.fa > resfinder_{name}.tab

        echo "starting argannot abrication"
        abricate --db argannot {name}.fa > argannot_{name}.tab

        echo "starting vfdb abrication"
        abricate --db vfdb {name}.fa > vfdb_{name}.tab

        echo "starting megares abrication"
        abricate --db megares {name}.fa > megares_{name}.tab

        echo "starting ecoh abrication"
        abricate --db ecoh {name}.fa > ecoh_{name}.tab

        echo "starting card abrication"
        abricate --db card {name}.fa > card_{name}.tab



        rm {name}.fa

    """

    return inputs, outputs, options, spec



def prokka(target_dir, title, name):
    
    #stem = '.'.join(name.split('.')[:-1])

    inputs  = target_dir + '/output/' + title + '/' + name + '/contigs.fa' 
    outputs = target_dir + '/output/' + title + '/' + name + '/' + name + '.gff'
    options = {'nodes': 1, 'cores': 8, 'memory': '4g', 'walltime': '1-12:00:00', 'account': 'clinicalmicrobio'} # initially 2 hours
    spec = f"""

        {environment}


        cd {target_dir}/output/{title}/{name}


        # hash tables

        # Generate hash key from assembly
        hash=$(cat contigs.fa | sha256sum | awk '{{print $1}}')

        # Set up directories
        hash_base="/faststorage/project/ClinicalMicrobio/database/prokka_hash/keys"
        hash_key_dir="${{hash_base}}/${{hash}}"

        # Check if the key exists
        if [[ -d "${{hash_key_dir}}" ]]; then
            echo "Hash key exists"
            echo $hash
            mkdir -p prokka
            
            # log usage
            echo -e "copying from ${{hash}}" > prokka/prokka_hash.txt
            echo -e "copy\t$(pwd)/prokka\t${{hash}}\t{name}\t$(date +%F_%H-%M-%S)" | tee "${{hash_key_dir}}"/usage_log.tab
            
            cp "${{hash_key_dir}}/prokka."* prokka

            # Touch it all to update the modified dates
            touch prokka/*
            
            echo -e "${{hash}}" > prokka/hash.txt

        else
            
            prokka --cpu 8 --force --outdir prokka --prefix prokka contigs.fa 

            if [[ -d "${{hash_base}}" ]]; then
                mkdir -p "${{hash_key_dir}}"

                cp prokka/prokka.* "${{hash_key_dir}}"

                echo -e "{name}\t${{hash}}\tsha256sum\tasscom\t$(date)" >> ${{hash_key_dir}}/index.tab
                touch ${{hash_key_dir}}/{name}.sample_name

                cat ${{hash_key_dir}}/index.tab >> ${{hash_base}}/index.tab
            fi

        fi

        cp prokka/prokka.gff {name}.gff

    """
            
    return inputs, outputs, options, spec


# MLST: Multi Locus Sequence Typing
def mlst(target_dir, title, contigs):
    inputs = [target_dir + '/output/' + title + '/' + i for i in contigs]
    outputs = target_dir + '/output/' + title + '/mlst.tsv' # Denne fil skal bruges til at lave træet, så det er den vigtigste. Og så også en liste over alle .gff-filer som er brugt.
    
    options = {'nodes': 1, 'cores': 2, 'memory': '4g', 'walltime': '01:00:00',  'account': 'clinicalmicrobio'}
    spec = f'''

        {environment}

        cd {target_dir}/output/{title}
        #mkdir mlst


        mlst {' '.join(contigs)} > mlst.tsv {debug('mlst')}

    '''

    return (inputs, outputs, options, spec)	



# Roary: The pan genome pipeline
def roary(target_dir, title, gffs, blastp_identity = 95, allow_paralogs = False):
    # target_dir:   Der hvor den skal gemme outputtet.
    # gffs:         En liste med fulde stier til de .gff-filer som skal analyseres.

    hours = -(-len(gffs)//25)*12 # 12 timer for hver 100 filer
    ram = -(-len(gffs)//25)*8 # 8 for hver 100 filer

    if allow_paralogs == True:
        ap_string = "-ap"
    else:
        ap_string = ""

    inputs = [target_dir + '/output/' + title + '/' + i for i in gffs]
    outputs = [target_dir + '/output/' + title + '/roary/core_gene_alignment.aln', # Denne fil skal bruges til at lave træet, så det er den vigtigste. Og så også en liste over alle .gff-filer som er brugt.
               target_dir + '/output/' + title + '/core_gene_alignment.fasta',
               target_dir + '/output/' + title + '/cg_snp_dists.tab',
               target_dir + '/output/' + title + '/roary/gene_presence_absence.csv',
               target_dir + '/output/' + title + '/roary/Rplots.pdf',
               target_dir + '/output/' + title + '/roary_thresholds.txt',
               target_dir + '/output/' + title + '/roary/blastp_'  + str(BLASTP) + '.setting']
    newline_for_f_string_workaround = '\n'
    options = {'nodes': 1, 'cores': 16, 'memory': f'{ram}g', 'walltime': f'{hours}:00:00', 'account': 'ClinicalMicrobio'}
    spec = f'''

        {environment}


        cd {target_dir}/output/{title}

        echo "blastp_identity (--blastp) = {str(blastp_identity)}" > roary_thresholds.txt

        roary -f roary -e -v -r -p 16 -g 100000 -i {int(blastp_identity)} {ap_string} {' '.join(gffs)} {debug('roary')}

        snp-dists roary/core_gene_alignment.aln > cg_snp_dists.tab

        cp roary/core_gene_alignment.aln core_gene_alignment.fasta


        touch roary/blastp_{str(BLASTP)}.setting


        echo JOBID $SLURM_JOBID
        jobinfo $SLURM_JOBID

    '''

    return (inputs, outputs, options, spec)


def fasttree(target_dir, title, n):
    # todo: time and mem should depend on number of isolates
    inputs = target_dir + '/output/' + title + '/roary/core_gene_alignment.aln'
    outputs = [target_dir + '/output/' + title + '/fasttree/tree.newick',
               target_dir + '/output/' + title + '/fasttree/tree.pdf']
    options = {'nodes': 1, 'cores': 8, 'memory': f'{round(n*2, 0)}g', 'walltime': f'{n*60}:00', 'account': 'clinicalmicrobio'}
    spec = f"""

        {environment}


        cd {target_dir}/output/{title}
        mkdir -p fasttree
        cd fasttree


        FastTree -nt -gtr ../roary/core_gene_alignment.aln > tree.newick 2> ft.log #{debug('ft')}


        Rscript /project/ClinicalMicrobio/faststorage/compare/scripts/R/ape_newick2pdf.r tree.newick "{title} core genome"  {debug('R')} || touch tree.pdf

    """

    return inputs, outputs, options, spec




def iqtree(target_dir, title, n):
    # todo: time and mem should depend on number of isolates
    inputs = target_dir + '/output/' + title + '/roary/core_gene_alignment.aln'
    outputs = [target_dir + '/output/' + title + '/_iqtree/run/core_gene_alignment.fasta.log',
               target_dir + '/output/' + title + '/_iqtree/run/core_gene_alignment.fasta.treefile',


               target_dir + '/output/' + title + '/_iqtree/' + title + '.newick',
               target_dir + '/output/' + title + '/_iqtree/iqtree_ml_bootstrap.pdf']
    options = {'nodes': 1, 'cores': 8, 'memory': f'{round(n*2, 0)}g', 'walltime': f'{n*1}:00:00', 'account': 'clinicalmicrobio'}
    spec = f"""

        {environment}

        cd {target_dir}/output/{title}
        mkdir -p _iqtree/run
        cd _iqtree/run

        cp ../../roary/core_gene_alignment.aln core_gene_alignment.fasta

        iqtree -s core_gene_alignment.fasta -bb 1000 -nt 8 -redo

        cp core_gene_alignment.fasta.treefile ../{title}.newick
        cd ..
        Rscript ../../../scripts/R/gg_newick2graphics.r {title}.newick run/core_gene_alignment.fasta.mldist {title}

    """

    return inputs, outputs, options, spec








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

        {environment}

        cd {target_dir}/output/{title}
        mkdir -p roary_plots
        cd roary_plots



        python {script_file} ../fasttree/tree.newick ../roary/gene_presence_absence.csv 

        perl {target_dir}/scripts/perl/roary2svg.pl ../roary/gene_presence_absence.csv > pangenome_matrix_alternative.svg
        cairosvg pangenome_matrix_alternative.svg -o pangenome_matrix_alternative.pdf

    '''

    return inputs, outputs, options, spec



def panito(target_dir, title):
    inputs = [target_dir + '/output/' + title + '/roary/core_gene_alignment.aln',
              target_dir + '/output/' + title + '/fasttree/tree.newick']
    outputs = target_dir + '/output/' + title + '/panito/ani.pdf'
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '00:10:00', 'account': 'clinicalmicrobio'}
    spec = f'''

        {environment}


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
              target_dir + '/output/' + title + '/kraken2-table.txt',
              target_dir + '/output/' + title + '/abricate/abricate_resfinder_summary.tab',
              target_dir + '/output/' + title + '/abricate/abricate_vfdb_summary.tab',
              target_dir + '/output/' + title + '/roary_thresholds.txt',
              target_dir + '/output/' + title + '/cg_snp_dists.tab',
              target_dir + '/output/' + title + '/core_gene_alignment.fasta',
              target_dir + '/output/' + title + '/_iqtree/iqtree_ml_bootstrap.pdf']
    #for name in names:
    #    inputs.append(target_dir + '/output/' + title + '/abricate/' + name)

    newline = '\n'
    outputs = target_dir + '/output/' + title + '/mailsent_FORCEAGAIN' # If it doesn't have an arbitrary output, the first job (init) will be run
    options = {'nodes': 1, 'cores': 1, 'memory': '1g', 'walltime': '00:10:00', 'account': 'clinicalmicrobio'}
    awk_command = """awk '{printf("%s %s (%s)\\n", $1, $2, $3)}'"""
    spec = f"""

        {environment}




        cd {target_dir}/output/{title}



        echo "this is the email address collected in assemblycomparator through finger"
        echo $COMPARATOR_DEFAULT_EMAIL_ADDRESS

        # collect mail content


        echo -e "Assembly Comparator results for {title}" > mail.txt

        echo -e "\n" >> mail.txt

        echo -e "MLST results:" >> mail.txt
        cat mlst.tsv | sed 's/\/contigs.fa//g' | {awk_command} | column -t >> mail.txt

        echo -e "\n" >> mail.txt


        echo "Roary summary statistics:" >> mail.txt
        cat roary/summary_statistics.txt | column -ts $'\t' >> mail.txt
        cat roary_thresholds.txt | grep "blastp" >> mail.txt

        echo -e "\n" >> mail.txt

        echo -e "Top 3 kraken results:" >> mail.txt

        cat kraken2-table.txt >> mail.txt


        echo -e "\n" >> mail.txt


        echo -e "A few small output files from the pipeline have been attached in the zip-file" >> mail.txt
        echo -e "To access the full analysis, please visit /project/ClinicalMicrobio/faststorage/compare/output/{title} on GenomeDK." >> mail.txt




        zip -j {title}.zip {' '.join(inputs)}


        mail -v -s "[comparator] done: {title}" -a {title}.zip -q mail.txt $COMPARATOR_DEFAULT_EMAIL_ADDRESS <<< "" 



        rm {title}.zip

        echo $(date) $COMPARATOR_DEFAULT_EMAIL_ADDRESS >> mailsent



    """

    #inputs.append(target_dir + '/output/' + title + '/completed_abricate')
    return inputs, outputs, options, spec




