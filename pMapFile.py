descriptionText = '''
Program to substitute a column in a tab file by its mapping given in another
tab file.

The user can select the columns to be mapped.

Special situations are taken like this:

* If there is no mapping available, the user can choose:
    - Save the non-mapped keys in a text file (default. If not file supplied
      then stderr will be used).
    - Printing the line keeping the original key or ignoring the line.    
    - Skip those lines

* If there is more than one mapping possibility the user can decide to:
    - Select the first mapping available
    - Select a random mapping
    - Replicate the line as many times as possible mapping are (default)

'''
#! /usr/bin/python 

import glob
import os
import argparse
import sys
import commands
import validateArguments

def dicFromFile(fileName,cKey,cMapped,sep,ifDup="replicate"):
    """First column in fileName has number 1"""
    
    if ifDup not in ["replicate","keepFirst","random"]:
        print "ERROR! dicFromFile -> ifDup shold be in replicate, keepFirst, random"
        sys.exit(-1)
    
    with open(fileName) as f:
        d = {}
        l = f.readline()
        while l:
            parts = [p.replace("\n","") for p in l.split(sep)]
            if not d.has_key(parts[cKey-1]): d[parts[cKey-1]] = [parts[cMapped-1]]
            else: d[parts[cKey-1]].append(parts[cMapped-1])
            l = f.readline()
    
    for k in d:
        if ifDup == "replicate": 
            pass # Keep all the mapping possibilities
        elif ifDup == "keepFirst": 
            d[k] = [d[k][0]] # Keep the first one
        elif ifDup == "random": 
            import random
            #if len(d[k]) > 3:
                #print random.choice(range(len(d[k])))
                #raw_input()
            d[k] = [d[k][random.choice(range(len(d[k])))]] # Select a random one

    
    return d

def mapLine(l,d,cChanged,sep,ifNotMapping="keepOriginal",fileNotMapped=None):
    """First column in fileName has number 1"""
    ret = []
    if ifNotMapping not in ["keepOriginal","skip","sendToFile"]:
        print "ERROR! dicFromFile -> ifDup shold be in keepOriginal, skip, sendToFile"
        sys.exit(-1)
    if ifNotMapping == "sendToFile" and not fileNotMapped:
        print "ERROR! dicFromFile -> selected sendToFile but not file provided!"
        sys.exit(-1)


    parts = [p.replace("\n","") for p in l.split(sep)]
    #print l
    #print parts[cChanged-1]
    #raw_input()
    #print d[parts[cChanged-1]]
    #print all([i.isspace() for i in d[parts[cChanged-1]]])
    #print all([len(i) == 0 for i in d[parts[cChanged-1]]])
    #raw_input()
    if d.has_key(parts[cChanged-1]) and not all([i.isspace() for i in d[parts[cChanged-1]]]) and not all([len(i) == 0 for i in d[parts[cChanged-1]]]):
        #print parts[cChanged-1], d[parts[cChanged-1]]
        #raw_input()
        for mapping in d[parts[cChanged-1]]:
            
            if cChanged==1:
                curLine="%s"%(sep.join([mapping,sep.join(parts[cChanged:])]))
                ret.append(curLine)
            elif cChanged==len(parts):
                #print sep.join([sep.join(parts[0:cChanged-1])])
                #print d[parts[cChanged-1]]
                #raw_input()
                curLine="%s"%(sep.join([sep.join(parts[0:cChanged-1]), mapping]))
                ret.append(curLine)
            else: # If it is the first column we skeep the first  
                curLine="%s"%(sep.join([sep.join(parts[0:cChanged-1]),mapping,sep.join(parts[cChanged:])]))
                ret.append(curLine)
    else:
        if ifNotMapping == "skip":
            pass
        elif ifNotMapping == "keepOriginal":
            curLine="%s"%(sep.join(parts))
            ret.append(curLine)
        elif ifNotMapping == "sendToFile":
            curLine="%s"%(sep.join(parts))
            #print curLine
            #raw_input()
            fileNotMapped.write("%s\n"%curLine)

    return ret


def main():
    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fileIn", type=validateArguments.existsFile, help="Input file")
    parser.add_argument("fileMap", type=validateArguments.existsFile, help="Mapping file")
    parser.add_argument("-o", "--fileOut", dest="o", type=argparse.FileType("w"), default="-",help="Output file. If not used output will be stdouti. Default \"-\" (stdout)", metavar="fileOut")
    parser.add_argument("-s", "--sep", dest="sep", default="\t", type=validateArguments.singleCharacter, help="Separator character", metavar="chr")
    parser.add_argument("--columnChanged", dest="cChanged", type=int, default=1, help="Number of the column in fileIn to be changed (default=1)",metavar="INT")
    parser.add_argument("--columnKey", dest="cKey", type=int, default=1, help="Number of the column in fileMap used as key (default=1)",metavar="INT")
    parser.add_argument("--columnMapped", dest="cMapped", type=int, default=2, help="Number of the column in fileMap to be used as a mapping (default=2)",metavar="INT")
    parser.add_argument("--ifDuplicates", dest="ifDup", type=str, choices=["replicate","keepFirst","random"], default="replicate", help="What to do in case that there is more than 1 mapping possibility. Values = replicate, keepFirst, random. See description for more info (default = replicate)",metavar="solution")
    parser.add_argument("--ifNotMapping", dest="ifNotMap", type=str, choices=["keepOriginal","skip", "sendToFile"], default="keepOriginal", help="What to do in case that there is NO mapping. Values = keepOriginal, skip, sendToFile. See description for more info (default = sendToFile)",metavar="solution")
    parser.add_argument("--fileNotMapped", dest="fileNotMapped", type=argparse.FileType("w"), default=sys.stderr, help="Output file for not mapped rows. By default it will be stderr", metavar="fileOut")
    parser.add_argument("-q", "--quiet",action="store_false", dest="verbose", default=True,help="don't print status messages to stderr")
    args = parser.parse_args()
    
    #if not args.o: fileOut = sys.stdout
    #else: fileOut = open(args.o,"w")
    fileOut = args.o
    d = dicFromFile(args.fileMap,cKey=args.cKey,cMapped=args.cMapped,sep=args.sep,ifDup=args.ifDup)
    
    with open(args.fileIn) as fIn:
        #with open(args.o,'w') as o
        notMapped = 0 
        mapped = 0
        totalMappings = 0
        l = fIn.readline()
        while l:
            newLines = mapLine(l,d,cChanged=args.cChanged,sep=args.sep,ifNotMapping=args.ifNotMap,fileNotMapped=args.fileNotMapped)
            #print newLines
            #raw_input()
            if len(newLines) == 0 or all([i.isspace() for i in newLines]) or all([len(i) == 0 for i in newLines]):
                notMapped += 1
            else:
                mapped += 1
                for newLine in newLines:
                    totalMappings += 1
                    fileOut.write("%s\n"%newLine)

            l = fIn.readline()
    if args.verbose: sys.stderr.write("Number of mapped elements: %d\n"%(mapped))
    if args.verbose: sys.stderr.write("Number of non mapped elements: %d\n"%(notMapped))
    if args.verbose: sys.stderr.write("Number of total mappings: %d\n"%(totalMappings))

if __name__ == "__main__":
    main()

