#!/usr/bin/env python3

import sys
import argparse

__author__ = "Carl Mathias Kobel"
__version__ = "0.1.1" # Not tested enough


def check_non_negative_int(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue



# argparser stuff
parser = argparse.ArgumentParser(prog = "fasta2tabseq.py", description='Convert .fasta to .tabseq, from stdin to stdout. The fasta record is written to the part column.')
parser.add_argument('--version', action='version', version=f"tabseq2fasta.py v{__version__}")

# Options to modify columns
parser.add_argument('--fill_sample', nargs='?', help='overwrite the sample data with this value', default = "", metavar = "string", type = str) 
#parser.add_argument('--fill_part', nargs='?', help='overwrite the part data with this value', default = False) # disabled, the record is saved in the part field.
parser.add_argument('--fill_comment', nargs='?', help='overwrite the comment data with this value', default = "", metavar = "string", type = str)

# Options to remove data
# group = parser.add_argument_group('optionally, omit data from being written to the fasta records')
# group.add_argument("--clear_sample", help = "clear sample data", action = "store_true")
# group.add_argument("--clear_part", help = "clear part data", action = "store_true")
# group.add_argument("--clear_comment", help = "clear comment data", action = "store_true")
# group.add_argument("--clear_all", help = "shortcut for the above clear-arguments", action = "store_true")
# 
# # Option to enumerate the records, and even specify the number of padded zeros
# parser.add_argument("--enumerate", help = "prefix the sequences with numbers, optionally specify the number of padded zeros, defaults to 3 padded zeros", nargs = "?", type = check_non_negative_int, metavar = "number of zeros", default = False, const = 3)
# 
# # 
# # 
# parser.add_argument("--separator", help = "choose a symbol to separate the fasta records, defaults to \"|\"", default = "|")



args = parser.parse_args()




def eprint(*args, **kwargs):
    print(*args, **kwargs, file = sys.stderr)
 
eprint(f"""
    compatible with multi line fasta files.
    Parsing fasta from stdin...
""")

# TODO: Implement --fill_sample



# Print the header
print('#sample', 'part', 'comment', 'sequence', sep = '\t')


def write(fasta_header, dna):
    print(f"{args.fill_sample}\t{fasta_header}\t{args.fill_comment}\t{dna}")


fasta_header = ""
dna = ""
init = True
for line in sys.stdin: #[i for i in sys.stdin] + [">"]:


    if line[0] == ">":
        
        # Output contents
        if not init:
            write(fasta_header, dna)
        
        init = False

        dna = "" # Reset buffer

        # Get ready for new DNA
        fasta_header = line[1:].strip() #.replace("\t", " ").replace(" ", "_")
        
    elif line[0] == "\n" or line[0] == "#":
        continue
    else:
        
        dna += line.strip()

write(fasta_header, dna)


