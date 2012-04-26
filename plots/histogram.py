"""
argv[1] fileIn
argv[2] fileOut (png/pdf)
argv[3] title
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys



x = []

file = open(sys.argv[1])

for line in file:
    #print line
    val = int(line)
    #if (val > 3000): 
        #val = 3000
    x.append ( val )

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title(sys.argv[3])

# the histogram of the data
n, bins, patches = ax.hist(x, 50, normed=0, facecolor='green', alpha=0.75)

# hist uses np.histogram under the hood to create 'n' and 'bins'.
# np.histogram returns the bin edges, so there will be 50 probability
# density values in n, 51 bin edges in bins and 50 patches.  To get
# everything lined up, we'll compute the bin centers
#bincenters = 0.5*(bins[1:]+bins[:-1])
# add a 'best fit' line for the normal PDF
#y = mlab.normpdf( bincenters, mu, sigma)
#l = ax.plot(bincenters, 'r--', linewidth=1)

#print bins,n,patches


#ax.set_xlim(0, 3000)
ax.set_xlabel("Length of DMRs")
ax.set_ylabel("Number of DMRs")
ax.grid(True)

plt.savefig(sys.argv[2])

