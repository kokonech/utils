descriptionText = '''
Program that receives two or more files with one column with the names of the features sorted by a score (the score should not be in the files) and returns thei similari
The program applies the Canberra similarity to the lists. ("Algebraic Comparison of Partial Lists in Bioinformatics" - http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0036540)
The implementation of the canberra similarity was found in the function "canberraq" of the mlpy package (v2.1 - https://launchpad.net/ubuntu/lucid/+source/mlpy/2.1.0~dfsg1-2/+files/mlpy_2.1.0~dfsg1.orig.tar.gz, in the current one, v3.5, there was no possibility of having list with differen number of element)
'''

import sys
import argparse
import os.path

import validateArguments

import numpy as np
import mlpy_v2

def readFile(fName):
    f = open(fName)
    r = [l.strip() for l in f.readlines()]
    return r


def createAlphabet(lists):
    """Used when no alphabet is provided.
    It constructs the union of all the elements on the lists in "lists" and assign them an index for later reference"""
    s = set()

    alphabet = {}
    for l in lists:
        s= s.union(set(l))
    for i, feature in enumerate(s):
        alphabet[feature] = i
    return alphabet

def composeAlphabetPositions(fName):
    """It creates the dicc with de index of the features for later reference"""
    alphabet = {}
    l = readFile(fName)
    for i, feature in enumerate(l):
        alphabet[feature] = i
    return alphabet

def formatLists(lists,alphabet):
    formattedLists = []
    for l in lists:
        curArray = np.zeros(len(alphabet),dtype=int) -1 # Fill all the elements with -1 (they are called as not present unless they are found afterwards)
        for rank, feature in enumerate(l):
            curArray[alphabet[feature]] = rank
        formattedLists.append(curArray)

    return np.array(formattedLists)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("fileIn1", type=validateArguments.existsFile, help="First input file")
    parser.add_argument("fileInList", type=validateArguments.existsFile, nargs="+", help="One or more input files")
    parser.add_argument("-a","--alphabet", type=validateArguments.existsFile, default=None,help="Alphabet File")
    parser.add_argument("-o", "--fileOut", dest="o", type=argparse.FileType("w"), default="-",help="Output file. If not used output will be stdout. Default \"-\" (stdout)", metavar="fileOut")
    parser.add_argument("-c","--core", dest="core", action='store_true', help="Appl")

    parser.add_argument("-q", "--quiet",action="store_false", dest="verbose", default=True,help="don't print status messages to stderr")


    args = parser.parse_args()

    files = [args.fileIn1]
    files.extend(args.fileInList)
    #args.fileInList.extend([args.fileIn1]) # Putting together the input files
    #features = [readFile(f) for f in args.fileInList]
    features = [readFile(f) for f in files]

#    for i in range(len(files)):
#        print "Len %d:\t%d\t%d"%(i+1, len(features[i]),len(set(features[i])))

    if args.alphabet:
        alphabet = composeAlphabetPositions(args.alphabet) # Created from input file. Dicc with feature -> index. Used to have the same order in the different lists
    else:
        alphabet = createAlphabet(features) # Created with the union of the elements in all the lists.

    formattedLists = formatLists(features,alphabet)

#    for i in range(len(formattedLists)):
#        print "List", i+1
#        for j in range(len(formattedLists[i])):
#            print formattedLists[i][j]


    r = mlpy_v2.canberraq(formattedLists,dist=True,complete=(not args.core))
    print r