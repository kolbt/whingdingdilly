'''
#                           This is an 80 character line                       #
This file reads in processed data from text files and outputs whether or not the
given system undergoes phase separation.  This is currently computed by the size
and lifetime of the largest cluster, but, could also feasibly arise from the
percentage of clustered particles in the system (not necessarily one cluster).

This file:
    1.  Reads in post-processed text data
    2.  Examines the largest cluster size
    3.  Writes to a phase separation file the values: peA, peB, xA, ps
'''

import sys
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

phaseTxt = sys.argv[1]

# Import data into arrays
peA, \
peB, \
xA, \
ps = np.loadtxt(phaseTxt, skiprows=1, unpack=True)

xA /= 100.0

for i in range(0, len(peA)):
    if ps[i] == 1:                                      # is phase sep.
        plt.scatter(peA[i], xA[i], facecolor='k', edgecolor='k')
    else:                                               # is not phase sep.
        plt.scatter(peA[i], xA[i], facecolor='none', edgecolor='k')

plt.title(r'Active/Passive Phase Separation $\epsilon=1$')
plt.xlabel(r'$Pe_{Fast}$')
plt.ylabel(r'$x_{Fast}$')
plt.savefig('phase-seperation-plot.png', dpi=1000)
