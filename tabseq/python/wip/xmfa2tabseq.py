#!/usr/bin/env python3

import sys
import re

# Original from biopypthon
XMFA_HEADER_REGEX = re.compile(r'> (?P<id>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[\+-]) (?P<name>.*)')


# Used in gBGC pipeline earlier when binning was a thing.
XMFA_HEADER_REGEX = re.compile(r"> *(?P<gene>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) isolate (?P<seqid>\d+)-(?P<strainid>[0-9a-zA-Z]+) bin (?P<bin>\d+)") 

# should fit Maria's raw data
XMFA_HEADER_REGEX = re.compile(r"> *(?P<sample>[0-9a-zA-Z-]+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) (?P<gene>[0-9a-zA-Z]+)") 


def eprint(string, *args, **kwargs):
    print(string, *args, **kwargs, file = sys.stderr)

eprint('This is the xmfa2tabseq.py conversion script.')
eprint('Please make sure that the headers in your .xmfa file fits the regex pattern given in XMFA_HEADER_REGEX ')
eprint(' Current XMFA_HEADER_REGEX:\n  ', XMFA_HEADER_REGEX)



eprint('Parsing xmfa from stdin...')

print('#species', 'sample', 'part', 'comment', 'sequence', sep = '\t')
# Parse header

current_sequence = ''
num = 0
first_line = True

for line in sys.stdin:
    elif line[0] == '>': # New header commences

        # Check to see if something needs to be written:
        if current_sequence != '':
            print('NA', current_sample, current_gene, f"strand:{current_header['strand']};start:{current_header['start']};end:{current_header['end']}", \
                #current_sequence[:10] + '...' + current_sequence[-10:], \
                current_sequence, \
                sep = '\t')

            current_sequence = ''
            num += 1




        m = re.match(XMFA_HEADER_REGEX, line)
        

        # These keys are strictly necessary
        current_gene = m.group('gene')
        current_sample = m.group('sample')
        current_line = line.strip()
        
        # These keys are not strictly necessary
        current_header = {} 
        keys = ('gene', 'sample', 'start', 'end', 'strand')
        for key in keys:
            try: 
                value = str(m.group(key))
                if key == 'start': 
                    value = int(value) 
                    # Convert to zero based counting 
                    if value > 0: 
                        value -= 1 

                if key == "end": 
                    value = int(value) 
                current_header[key] = value 
            except IndexError: 
                current_header[key] = 'NA'
                # This will occur if we're asking for a group that 
                # doesn't exist. It's fine as long as it is not the gene, which we already checked.
                # pass


        # This is solely for debugging:
        if first_line == True:
            eprint(f" The first header \"{line.strip()}\" is parsed as follows:")
            string = ''
            for key in keys:
                string += f"  {key}: \"{current_header[key]}\"\n"

            eprint(string, end = '')
            eprint('Writing to stdout...')
            first_line = False







    elif line[0] == '#' or line[0] == '=': # comment or alignment end line, skip
        continue

    else: # DNA
        current_sequence += line.strip()


eprint('Wrote', num, 'sequences.\n')