descriptionText = '''
Program to aggregate the values of a tab file.

One can select:
    - which column will be used as a key for the aggregation
    - which column will be used for obtaining the values
    - the function for the aggregation

Example:
    Input:
        A   1
        A   3
        C   3
    Would produce:
        A   2
        C   3
'''
#! /usr/bin/python 

import glob
import os
import argparse
import sys
import commands
import numpy
import validateArguments

def dicFromFile(fileName,cKey,cMapped,sep="\t"):
    """First column in fileName has number 1"""
    
    
    with open(fileName) as f:
        d = {}
        l = f.readline()
        while l:
            parts = [p.replace("\n","") for p in l.split(sep)]
            #print parts
            if not d.has_key(parts[cKey-1]): d[parts[cKey-1]] = [parts[cMapped-1]]
            else: d[parts[cKey-1]].append(parts[cMapped-1])
            l = f.readline()
   

    
    return d

def absMax(l):
    """Return the maximum in absolute value of a list"""
    return l[numpy.argmax(numpy.abs(l))]



def aggregate(d,aggMethod,discardDiscordant=False):
    acceptedMethods = ["mean","median","max","min","absMax"]
    if not aggMethod in acceptedMethods:
        raise AttributeError("aggMethod not in %s"%",".join(acceptedMethods))
    #method = decideAggMethod(aggMethod,discarDiscordant)

    if aggMethod == "mean": method = numpy.mean
    if aggMethod == "median": method = numpy.median
    if aggMethod == "max": method = numpy.max
    if aggMethod == "min": method = numpy.min
    if aggMethod == "absMax": method = absMax

    #r = [[k,method([float(i) for i in d[k]])] for k in d if not discardDiscordant or all([float(j) >= 0 for j in d[[k]]] or all([float(j) <= 0 for j in d[[k]]])]
    r = [[k,method([float(i) for i in d[k]])] for k in d if not discardDiscordant or all([float(j) >= 0 for j in d[k]]) or all([float(j) <= 0 for j in d[k]])]

    return r
        


def main():
    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fileIn", type=validateArguments.existsFile, help="Input File")
    parser.add_argument("-o", "--fileOut", dest="o", type=argparse.FileType("w"), default="-",help="Output file. If not used output will be stdout. Default \"-\" (stdout)", metavar="fileOut")
    parser.add_argument("-s", "--sep", dest="sep", default="\t", type=validateArguments.singleCharacter, help="Separator character", metavar="chr")
    parser.add_argument("--columnKey", dest="cKey", type=int, default=1, help="Number of the column in fileIn used as key (default=1)",metavar="INT")
    parser.add_argument("--columnAggregated", dest="cAggregated", type=int, default=2, help="Number of the column in fileIn to be aggregated (default=2)",metavar="INT")
    parser.add_argument("--aggMethod", dest="aggMethod", type=str, choices=["mean","median","max","min","absMax"], default="mean", help="Method to be used when aggregating rows",metavar="Aggregation Method")
    parser.add_argument("--discardDiscordant", action="store_true", default=False,help="Ignore IDs which have positive AND negative values. Useful for filtering inconsistent probes in a microarray for example.")
    parser.add_argument("-q", "--quiet",action="store_false", dest="verbose", default=True,help="don't print status messages to stderr")
    args = parser.parse_args()
    
    d = dicFromFile(args.fileIn,args.cKey,args.cAggregated,args.sep)

    res = aggregate(d, args.aggMethod, discardDiscordant=args.discardDiscordant)
    for i in res:
        args.o.write("%s\t%s\n"%(i[0],i[1])) 

if __name__ == "__main__":
    main()


