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
from textwrap import wrap

phaseTxt = sys.argv[1]
catesTxt = '/Users/kolbt/Desktop/compiled/whingdingdilly/phase_diagrammer/cates_phase_data.txt'

# Import data into arrays
peA, \
peB, \
xA, \
ps = np.loadtxt(phaseTxt, skiprows=1, unpack=True)

CpeA, \
CpeB, \
CxA, \
Cps = np.loadtxt(catesTxt, skiprows=1, unpack=True)

xA /= 100.0
CxA /= 100.0

# It'd be best to sort entries and use same index

# for i in range(0, len(peA)):
#     # Find matching Cates value
#     for j in range(0, len(Cps)):
#         if CpeA[j] == peA[i]:
#             for k in range(0, len(Cps)):
#                 if CxA[k] == xA[i]:
#                     print(Cps[k])
#                     if Cps[k] == ps[i] == 1:
#                         plt.scatter(peA[i], xA[i], facecolor='k', edgecolor='k')
#                     elif Cps[k] == ps[i] == 0:
#                         plt.scatter(peA[i], xA[i], facecolor='none', edgecolor='k')
#                     elif Cps[k] == 0 and ps[i] == 1:
#                         plt.scatter(peA[i], xA[i], facecolor='r', edgecolor='k')
#                     else:
#                         plt.scatter(peA[i], xA[i], facecolor='b', edgecolor='k')
#                     break

# Plot Cates data
for i in range(0, len(CpeA)):
    if Cps[i] == 1:
        black = plt.scatter(CpeA[i], CxA[i], facecolor='k', edgecolor='k')
    else:
        empty = plt.scatter(CpeA[i], CxA[i], facecolor='none', edgecolor='k')

# Plot my data on top of it (hacky but it works)
for i in range(0, len(peA)):
    if ps[i] == 1:
        blue = plt.scatter(peA[i], xA[i], facecolor='b', edgecolor='k')

legend = plt.legend(bbox_to_anchor=(0.0, 1.0),
                    loc='center left',
                    ncol=3,
                    handles=[empty, black, blue],
                    labels=['Gas', 'HS: Phase Separated', r'$\epsilon=1:$ Phase Separated'],
                    framealpha=1.0)
plt.title(r'Active/Passive Phase Separation', y=1.03)
plt.xlabel(r'$Pe_{Fast}$')
plt.ylabel(r'$x_{Fast}$')
plt.tight_layout()
plt.savefig('phase-seperation-plot.png', dpi=1000)
