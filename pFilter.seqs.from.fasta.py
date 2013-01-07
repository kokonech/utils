descriptionText = '''
Program to obtain the desired sequences from a multifasta file.
I create it to obtain a cleaner genome from hg19 (without random chromosomes and staff like that)
'''

import sys
import argparse
import os.path

import validateArguments

from Bio import SeqIO


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("fileIn", type=validateArguments.existsFile, help="Input file")
    parser.add_argument("-o", "--fileOut", dest="o", type=argparse.FileType("w"), default="-",help="Output file. If not used output will be stdout. Default \"-\" (stdout)", metavar="fileOut")
    parser.add_argument("chrs", nargs="+", help="Sequence names")

    parser.add_argument("-q", "--quiet",action="store_false", dest="verbose", default=True,help="don't print status messages to stderr")


    args = parser.parse_args()

    seqiter = SeqIO.parse(open(args.fileIn), "fasta")
    SeqIO.write((seq for seq in seqiter if seq.id in args.chrs), args.o, "fasta")


