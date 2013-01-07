import re
import sys
import argparse
import os.path
import time
import tempfile
import subprocess
import os

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

def analyzeSeqs(seqFile, outBlast, outNoHits, db_name, numThreads):
    

    cline = NcbiblastnCommandline(query=seqFile, db=db_name, task="megablast", evalue=0.001, outfmt=6,out=outBlast,num_threads=numThreads)
    
    print time.ctime()
    print cline
    cline()
    print time.ctime()
    
    if outNoHits:
        (fd, tmpIdfileName) = tempfile.mkstemp()
        print "Name of the tmp file for IDs:", tmpIdfileName
        getIDs = "awk 'BEGIN{OFS=\"\\t\"}NR % 2 == 1{id=split($1,a,\">\");print a[2]}' " + seqFile + " | sed 's/ /_/g' | sort -k1,1 > " + tmpIdfileName
        print getIDs
        os.system(getIDs)
        getNoHits = "cut -f 1 %s | sort | uniq | comm -31 - %s > %s"%(outBlast,tmpIdfileName,outNoHits)
        print getNoHits
        os.system(getNoHits)
        os.remove(tmpIdfileName)


            
if __name__ == "__main__":
    
    descriptionText = "Program to execute BLAST against a set of sequences provided in a FASTA file"
    
    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("seqFile", help="Input FASTA file")
    parser.add_argument("dbName", help="BLAST database name")
    parser.add_argument("-o", help="BLAST Output")
    parser.add_argument("-n","--noHits", dest="n", default=None, help="File for the seqs with no hit in the database. If not defined, they will not be output")
    parser.add_argument("-t","--numberThreads", dest="t", default=1, help="Number of threads for blasting")



    args = parser.parse_args()


    res = analyzeSeqs(args.seqFile, args.o,args.n,args.dbName,args.t)






