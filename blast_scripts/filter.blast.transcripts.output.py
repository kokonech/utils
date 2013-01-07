descriptionText = '''
Program to filter the output the table produced by blasting microarray probes against transcripts (done via blastFromFasta.py. For example see "/data/medips/microarrays/probes.blast.64/transcripts.blast")

The idea is to get a singel gene for each probe, discarding those hits that do not suffice the Quality that the users want
'''
#! /usr/bin/python 

import glob
import os
import argparse
import sys
import commands
import validateArguments


class RowBlast:

    def __init__(self,l,posProbeID=1,posIDorig=2,posIDblasted=1,infoOrig=[1],infoBlasted=[1,2,3]):
        """
        posProbeID: Position of the Probe ID (found in the 1st field) (1 is the first position)
        posIDorig: Position of the original ID which we want to aggregate (found in the first field) (1 is the first position)
        posIDblasted: Position of the blasted ID which we want to aggregate (found in the second field) (1 is the first position)
        infoOrig: List with the positions of the additional information from the probe annotation one wants to output (found in the first field) (1 is the first position)
        infoBlasted: List with the positions of the additional information from the BLAST annotation one wants to output (found in the second field) (1 is the first position)
        Example:
            A_24_P66027|APOBEC3B    APOBEC3B|ENSG00000179750|protein_coding|APOBEC3B-001|ENST00000333467|protein_coding|
        we would select posProbeID=1, posIDorig=2 and posIDblasted=1
        
        infoOrig = [1] and infoBlasted=[1,2,3] would print:
            A_24_P66027 APOBEC3B    ENSG00000179750 protein_coding
        while infoOrig = [1,2] and infoBlasted=[2,1,3] would print:
            A_24_P66027 APOBEC3B    ENSG00000179750 APOBEC3B    protein_coding
        """

        parts = [p for p in l.strip().split()]
        #print parts
        #print parts[0].split("|")
        #print parts[1].split("|")
        #print parts[1].split("|")[4]
        self.probeID = parts[0].split("|")[posProbeID-1]
        self.origID = parts[1].split("|")[posIDorig-1]
        self.blastedID = parts[1].split("|")[posIDblasted-1]
        self.percOverlap = float(parts[2])
        self.bpOverlap = int(parts[3])
        #print "\t".join([parts[0].split("|")[i-1] for i in infoOrig])
        #print infoBlasted
        #print "\t".join([parts[1].split("|")[i-1] for i in infoBlasted])
        self.info = "%s\t%s"%("\t".join([parts[0].split("|")[i-1] for i in infoOrig]), "\t".join([parts[1].split("|")[i-1] for i in infoBlasted]))


    def __str__(self):
        return "%s"%(self.info)

    def passQC(self,percOverlap=95.0,bpOverlap=55):
        return self.percOverlap >= percOverlap and self.bpOverlap >= bpOverlap

    def isControversial(self,r):
        return self.blastedID != r.blastedID

    def toFile(self,f):
        f.write("%s\n"%(str(self)))

def main():
    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fileIn", type=argparse.FileType("r"), help="Blast input file. See \"/data/medips/microarrays/probes.blast.64/transcripts.blast\" as a example")
    parser.add_argument("-o", "--fileOut", dest="o", type=argparse.FileType("w"), default="-",help="Output file. If not used output will be stdout. Default \"-\" (stdout)", metavar="FILE")
    parser.add_argument("-d", "--fileDiscarded", dest="d", type=argparse.FileType("w"), default=None,help="File to save the probes with no hit that fulfils the QC requirements. If not used they are ignored", metavar="FILE")
    parser.add_argument("-c", "--fileControversial", dest="c", type=argparse.FileType("w"), default=None,help="File to save the probes that have two or more good hits for different genes. If not used they are ignored", metavar="FILE")
    parser.add_argument("-p", "--percOverlap", dest="percOverlap", type=float, default=95.0,help="Minimum percentage of overlap between the probe and the transcript (default=95.0)", metavar="PERCENTAGE")
    parser.add_argument("--posProbeID", dest="posProbeID", type=int, default=1,help="Position where the unique ID for the probes in the array can be found. For example in \"ITGB3BP|ENST00000283568|A_23_P23765|3966238\" 3 or 4 could be used (default 1)", metavar="INT")
    parser.add_argument("-b", "--bpOverlap", dest="bpOverlap", type=float, default=55,help="Minimum base pair overlap between the probe and the transcript (default=55)", metavar="PERCENTAGE")
    parser.add_argument('--infoOrig', metavar='INT', type=int, nargs='+', default=[1],help="List with the positions of the additional information from the probe annotation one wants to output (found in the first field) (1 is the first position)")
    parser.add_argument('--infoBlasted', metavar='INT', type=int, nargs='+', default=[1,2,3], help="List with the positions of the additional information from the BLAST annotation one wants to output (found in the second field) (1 is the first position)")

    parser.add_argument("-q", "--quiet",action="store_false", dest="verbose", default=True,help="don't print status messages to stderr")
    args = parser.parse_args()

    numSingleGeneHits = 0
    numControversialProbes = 0
    numDiscardedProbes = 0
    curProbeID = None
    curBestHit = None
    l = args.fileIn.readline()
    while l:
        r = RowBlast(l,posProbeID=args.posProbeID,infoOrig=args.infoOrig,infoBlasted=args.infoBlasted)
        #print l.strip()
        #print r
        if not curProbeID: # if it is the FIRST probe
            curProbeID = r.probeID
            curBestHit = None # To store the best hit and compare the other hits with them
            curControversial = [] # To store de hits with good hits to different genes
        elif curProbeID != r.probeID: # if we find hits from a NEW probe we will store the previous one accordingly
            
            if not curBestHit: # There was no hit that passed the QC
                numDiscardedProbes += 1
                if args.d: # If we want to store them
                    args.d.write("%s\n"%(curProbeID))

            elif len(curControversial) == 0: # The probe has a single gene as a hit
                numSingleGeneHits +=1
                curBestHit.toFile(args.o)

            elif args.c: # There are controversial hits for the same probe and we want to store them
                numControversialProbes += 1
                for c in curControversial:
                    c.toFile(args.c)

            curProbeID = r.probeID
            curBestHit = None # To store the best hit and compare the other hits with them
            curControversial = [] # To store de hits with good hits to different genes
        #print curProbeID, curBestHit, "|".join([str(c) for c in curControversial])
        if r.passQC(percOverlap=args.percOverlap,bpOverlap=args.bpOverlap): # Skip reads that do not fulfil the QC requirements
            #print "Yes QC"
            if curProbeID == r.probeID: # if we find hits from the SAME probe
                #print "hits from the SAME probe"
                if not curBestHit:
                    #print "Fist hit for the probe -> new best hit"
                    curBestHit = r
                elif r.isControversial(curBestHit):
                    if len(curControversial) == 0:
                        curControversial.append(curBestHit)
                        curControversial.append(r)
                    else:
                        newGene = False
                        for contro in curControversial:
                            if not r.isControversial(contro): 
                                newGene = True
                                break
                        if newGene:
                            curControversial.append(r)
            else:
                #print "hits from different probe. Change"
                pass
        #print "curBestHit:", curBestHit
        #raw_input()


        l = args.fileIn.readline()

    # For the last probe
    if not curBestHit: # There was no hit that passed the QC
        numDiscardedProbes += 1
        if args.d: # If we want to store them
            args.d.write("%s\n"%(curProbeID))

    elif len(curControversial) == 0: # The probe has a single gene as a hit
        numSingleGeneHits +=1
        curBestHit.toFile(args.o)

    elif args.c: # There are controversial hits for the same probe and we want to store them
        numControversialProbes += 1
        for c in curControversial:
            c.toFile(args.c)


    if args.verbose:
        print "Discarded probes: %d"%(numDiscardedProbes)
        print "Controversial probes: %d"%(numControversialProbes)
        print "Probes with hits to one single gene: %d"%(numSingleGeneHits)

if __name__ == "__main__":
    main()
