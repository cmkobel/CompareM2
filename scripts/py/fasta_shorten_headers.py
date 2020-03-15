#!/usr/bin/env python3

# This program shortens the headers of a fasta file
# Pipe fasta file with STDIN, outputs with STDOUT

import sys



LENGTH = 30


for line in sys.stdin:
    if line[0] == '>':
        print(line[:LENGTH], end = '\n') 
    else:
        print(line, end = '') # Aleady contains a newline in the end.
