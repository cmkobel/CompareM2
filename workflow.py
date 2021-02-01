#!/usr/bin/env python

import os
from os.path import isfile, isdir, join, dirname, basename
from gwf import *
from workflow_templates import *


BLASTP = int(os.environ["BLASTP"])
print('BLASTP:', BLASTP)



gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
})


def sanify(input, allow_hyphens = False):
    """ Makes sure that the name of the gwf target is not illegal. """
    output = []
    
    additional_characters = [95]
    if allow_hyphens:
        additional_characters.append(45)

    for i in str(input):
        
        ascii = ord(i)
        if (ascii >= 48 and ascii <= 57) or (ascii >= 65 and ascii <= 90) or (ascii >= 97 and ascii <= 122) or ascii in additional_characters:
            output.append(i)
        else:
            output.append('_')

    return ''.join(output)


#Todo replace - with _ in things in the names of jobs

source_dir = os.getcwd()
target_dir = '/project/ClinicalMicrobio/faststorage/compare'
unsane_title = basename(os.getcwd()) # Thought as whatever folder the user is calling kmacompare from
title = sanify(unsane_title)
error_file = target_dir + '/' + title + '/stdout.txt'
#print(' workflow title:', title)
#print(' workflow source_dir:', source_dir)
#print(' workflow target_dir:', target_dir)
#print()

fasta_files = []


fasta_files = [f for f in os.listdir(source_dir) if isfile(join(f))]
n_samples = len(fasta_files)

# This is taken care of from the bash script
#print('These are the files considered:')
#for i in fasta_files:
#    print('', i)
#print()

# todo: sanify title, at least, test spaces in input: see 'test 2'

for file in os.listdir(source_dir):
    if isfile(file):
        with open(file, 'r') as opened_file:
            for line in opened_file:
                
                if line[0] != '>':
                    print(f'Warning: the file: {file} doesn\'t look like a fasta file. Consider its inclusion.')
                
                break


# Initialize
#print(dir(gwf))
gwf.target_from_template(sanify('cmp_init_' + title), initialize(title, source_dir, target_dir))
    
names = []
for raw_name in fasta_files:
    name = sanify(stem(raw_name), allow_hyphens = True) # the name should accept hyphens.
    names.append(name)
    
    #todo: slå copy og prokka sammen?

    # submit copy job
    gwf.target_from_template(sanify('cmp_copy_' + title + '_' + name), copy(source = source_dir + '/' + raw_name,
                                                                    target_dir = target_dir,
                                                                    title = title,
                                                                    name = name))
        

    # submit kraken2
    gwf.target_from_template(sanify('cmp_kraken2_' + name), kraken2(target_dir, title, name))

    # submit abricate
    gwf.target_from_template(sanify('cmp_abricate_' + name), abricate(target_dir, title, name))


    # submit prokka
    gwf.target_from_template(sanify('cmp_prokka_' + title + '_' + name), prokka(target_dir, title, name))



# run as group (list of inputs) from now on (names)
#gwf.target_from_template(sanify('cmp_kraken2_abricateall_' + title), summary_tables(target_dir, title, names))


gwf.target_from_template(sanify('cmp_kraken_sum_' + title), summary_kraken(target_dir, title, names))
gwf.target_from_template(sanify('cmp_abricate_sum_' + title), summary_abricate(target_dir, title, names))



contigs = [name + '/contigs.fa' for name in names]
gwf.target_from_template(sanify('cmp_mlst_' + title), mlst(target_dir, title, contigs))

named_contigs = [name + '/' + name + '.fa' for name in names]
gwf.target_from_template(sanify('cmp_mashtree_' + title), mashtree(target_dir, title, named_contigs))




# submit roary
#annotations = [name + '/prokka/' + name + '.gff' for name in names]
annotations = [name + '/' + name + '.gff' for name in names]

gwf.target_from_template(sanify('cmp_roary_' + str(BLASTP) + '_' + title), roary(target_dir, title, annotations, blastp_identity = BLASTP, allow_paralogs = False))


# submit fasttree
gwf.target_from_template(sanify('cmp_fasttree_' + title), fasttree(target_dir, title, len(names)))

# submit iqtree
gwf.target_from_template(sanify('cmp_iqtree_' + title), iqtree(target_dir, title, len(names)))




# submit roary plots
gwf.target_from_template(sanify('cmp_plot_roary_' + title), roary_plots(target_dir, title))

# submit panito
gwf.target_from_template(sanify('cmp_panito_' + title), panito(target_dir, title))


# send a mail
gwf.target_from_template(sanify('cmp_mail_' + title), send_mail(target_dir, title, names))



