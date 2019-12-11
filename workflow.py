#!/usr/bin/env python

import os
from os.path import isfile, isdir, join, dirname, basename
from gwf import *
from workflow_templates import *


gwf = Workflow()


#Todo replace - with _ in things in the names of jobs

source_dir = os.getcwd()
target_dir = '/project/ClinicalMicrobio/faststorage/compare'
title = basename(os.getcwd()) # Thought as whatever folder the user is calling kmacompare from
error_file = target_dir + '/' + title + '/stdout.txt'
print(' workflow title:', title)
print(' workflow source_dir:', source_dir)
print(' workflow target_dir:', target_dir)
print()

fasta_files = []


fasta_files = [f for f in os.listdir(source_dir) if isfile(join(f))]

# This is taken care of from the bash script
#print('These are the files considered:')
#for i in fasta_files:
#    print('', i)
#print()

for file in os.listdir(source_dir):
    if isfile(file):
        with open(file, 'r') as opened_file:
            for line in opened_file:
                
                if line[0] != '>':
                    print(f'Warning: the file: {file} doesn\'t look like a fasta file. Consider removing it.')
                
                break


# Initialize
#print(dir(gwf))
gwf.target_from_template('cmp_init_' + title, initialize(title, source_dir, target_dir))
        
names = []
for raw_name in fasta_files:
    
    name = stem(raw_name.replace(' ', '_'))
    names.append(name)
    
    # submit copy job
    gwf.target_from_template('cmp_copy_' + title + '_' + name, copy(source = source_dir + '/' + raw_name,
                                                                    target_dir = target_dir + '/output/' + title + '/' + name,
                                                                    target_file = 'contigs.fa'))
        

    # submit kraken2
    gwf.target_from_template('cmp_kraken2_' + name, kraken2(target_dir, title, name))


    # submit prokka
    gwf.target_from_template('cmp_prokka_' + title + '_' + name, prokka(target_dir, title, name))



# run as group (list of inputs) from now on (names)


contigs = [name + '/contigs.fa' for name in names]
gwf.target_from_template('cmp_mlst_' + title, mlst(target_dir, title, contigs))


# submit roary
annotations = [name + '/prokka/' + name + '.gff' for name in names]
gwf.target_from_template('cmp_roary_' + title, roary(target_dir, title, annotations))


# submit fasttree
gwf.target_from_template('cmp_fasttree_' + title, fasttree(target_dir, title))


# submit roary plots
gwf.target_from_template('cmp_roary_plots_' + title, roary_plots(target_dir, title))

# submit panito
gwf.target_from_template('cmp_panito_' + title, panito(target_dir, title))



