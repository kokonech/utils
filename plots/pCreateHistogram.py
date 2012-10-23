descriptionText = '''
Program to create a histogram from a text file.
The user can select the column from which the histogram will be computed
'''

import sys
import argparse
import os.path

import validateArguments

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator

def getValues(fileIn, column=2, sep="\t",verbose=False):
    with open(fileIn) as f:
        values = np.array([float(l.split(sep)[column-1]) for l in f.readlines()])
    
    if verbose:
        print "Number of read elements: %d"%len(values)

    #for i in values:
        #print i
        #raw_input()
    return values




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("fileIn", type=validateArguments.existsFile, help="Input file")
    parser.add_argument("-o", "--fileOut", help="Output file for the generated image. If not passed the image will be displayed")
    parser.add_argument("-c", "--column", dest="c", type=int, default="2", help="Column where the values for the histogram are found (default=2)",metavar="INT")
    parser.add_argument("-n", "--numBins", dest="n", type=int, nargs="*", default="15", help="Number of bins for the histogram (default=15)",metavar="INT")
    parser.add_argument("-s", dest="s", default="\t", type=str, help="Field Separator (default=\\t)",metavar="CHAR")
    parser.add_argument("--title", dest="title", type=str, help="Title of the figure",metavar="STR")
    parser.add_argument("--xLabel", dest="xLabel", type=str, help="Label for the X axis",metavar="STR")
    parser.add_argument("--yLabel", dest="yLabel", type=str, help="Label for the Y axis",metavar="STR")
    parser.add_argument("--numXticks", dest="numXticks", type=int, help="Number of ticks in the X axis",metavar="INT")

    parser.add_argument("-q", "--quiet",action="store_false", dest="verbose", default=True,help="don't print status messages to stderr")


    args = parser.parse_args()

    if len(args.n) == 1: args.n = args.n[0]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    #mu, sigma = 100, 15
    #x = mu + sigma*np.random.randn(10000)

    x = getValues(args.fileIn, column=args.c, sep=args.s, verbose=args.verbose)
    
    # the histogram of the data
    n, bins, patches = ax.hist(x, bins=args.n, normed=False, facecolor='green', alpha=0.75,cumulative=False,weights=(np.zeros_like(x) + 1. / x.size))

    if args.title:
        ax.set_title(args.title)
    #plt.axis([0, 300, 0, 0.004])
    #ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: ('%.2f')%(y*1e-3)))
    if args.numXticks:
        ax.xaxis.set_major_locator(MaxNLocator(args.numXticks))

    if args.xLabel:
        ax.set_xlabel(args.xLabel)
    if args.yLabel:
        ax.set_ylabel(args.yLabel)

    ax.grid(True)
    fig.canvas.draw()
    
    if args.fileOut:
        plt.savefig(args.fileOut,dpi=300)
    else:
        plt.show()

