import re
import sys
import argparse
import os.path
import time
import tempfile
import os

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

SEQ_BLASTED = 0
SEQ_NOT_FOUND = -1

MAX_BUF_SIZE = 100

def parseBlastResult(fileName):
    
    handle = open(fileName)
    blast_records = NCBIXML.parse(handle)
    
    results = []

    for record in blast_records:
        rec_id = str(record.query)
    
        if len(record.alignments) == 0:
            results.append( (rec_id, "-", 0, "-") )
            continue

        for algn in record.alignments:

            evalue = algn.hsps[0].expect
        
            score = 0
            ids = []
            
            for hsp in algn.hsps:
                score += hsp.bits
                ids.append(hsp.identities / float(hsp.align_length))
            
            max_identity = int(max(ids)*100)
            seq_id = algn.hit_id

            results.append( (rec_id, seq_id, max_identity, algn.hit_def ) )
            
    return results
  
def analyzeSeq(seqBuf, outputPath, db_name ):
    
    (fd, qname) = tempfile.mkstemp()
    tmp = os.fdopen(fd, "w")

    for seq in seqBuf:
        seq_id = seq[0]
        seq_data = seq[1]
        #print "Saving ", seq_id, seq_data
        tmp.write(">%s\n%s\n" % (seq_id,seq_data) )
    
    tmp.close()

    #print qname

    cline = NcbiblastnCommandline(query=qname, db=db_name, task="megablast", evalue=0.001, outfmt=5, out=outputPath)
    print cline
    try: 
        cline()
    except:
        os.remove(outputPath)
           
            
if __name__ == "__main__":
    
    descriptionText = "The script reads each sequence from an input .csv file and blasts it for specified organism. The accession id is then extracted from blast result."
    
    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("seqFile", help="Input file with list of sequences")
    parser.add_argument("dbName", help="BLAST database name")

    args = parser.parse_args()

    seqFileName = args.seqFile
    dbName = args.dbName
    print "Input sequences file is ", seqFileName
      
    seqFile = open(seqFileName)  
    sequences = {}
    missing = []

    for line in seqFile:
        items = line.strip().split(",")
        seq_id = items[0]
        try:
            seq = items[1].split("\"")[1]
            sequences[seq_id] = seq    
        except:
            print "No sequence available, skipping ", seq_id
        #print seq_id, seq
    
    
    
    workDirName = "local_analysis_%s_%d" % (seqFileName, MAX_BUF_SIZE)
    
    if not os.path.exists(workDirName):
        os.makedirs(workDirName)
    print "Workdir name is ", workDirName
    
    outputFileName = workDirName + "/" + "resume.txt"
    output = open(outputFileName, "w")
    output.write("#Seq id\tGenbank id\tUrl\tMax identity(%)\tDescription\n")
    
    total = len(sequences)
    bunch_count = 1
    i = 0

    buf = []

    for seq_id in sorted(sequences):
        outputPath = workDirName + "/" + "bunch_" + str(bunch_count) +".xml";
        #print i, seq_id
        if not os.path.exists(outputPath):
            i += 1
            buf.append ( (seq_id, sequences[seq_id]) )
            if len(buf) ==  MAX_BUF_SIZE or i == total:
                print time.ctime()
                res = analyzeSeq(buf, outputPath, dbName)
                print time.ctime()
                buf = []
            else:
                continue
        else:
            i += MAX_BUF_SIZE
        
        bunch_count += 1
        res = parseBlastResult(outputPath)

        for r in res:
            #print r
            output.write( "%s\t%s\t%d\t%s\n" % ( r[0], r[1], r[2], r[3]  )  ) 
        print "Analyzed %d out of %d sequences..." % ( i , total)
        print "\n\n"

    print "Missing sequences:"
    print missing




