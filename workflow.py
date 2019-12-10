#!/usr/bin/env python

import os
from os.path import isfile, isdir, join, dirname, basename
from gwf import *
from workflow_templates import *


gwf = Workflow()


source_dir = os.getcwd()
target_dir = '/project/ClinicalMicrobio/faststorage/compare'
title = basename(os.getcwd()) # Thought as whatever folder the user is calling kmacompare from
error_file = target_dir + '/' + title + '/stdout.txt'
print('title:', title)
print('source_dir:', source_dir)
print('target_dir:', target_dir)
print()

fasta_files = []


fasta_files = [f for f in os.listdir(source_dir) if isfile(join(f))]

print('These are the files considered:')
for i in fasta_files:
    print('', i)
print()

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
    
    gwf.target_from_template('cmp_copy_' + title + '_' + name, copy(source = source_dir + '/' + raw_name,
                                                                    target_dir = target_dir + '/output/' + title + '/' + name,
                                                                    target_file = 'contigs.fa'))
        


    # submit prokka
    gwf.target_from_template('cmp_prokka_' + title + '_' + name, prokka(target_dir, title, name, error_file))


