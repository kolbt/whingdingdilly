'''
#                           This is an 80 character line                       #
INTENT: take text file data and plot it.  Saving both the text files and the
plots allows for simple replotting without recalculating.
Note: I'll need to generate another file which compares between simulations
Data:   PeA PeB PeR xA Ep sigT sigAA sigAB sigBB fT fAA fAB fBB phiEff
( 1.) Read data into arrays
( 2.) Put it ALL in one plot
( 3.) Do a little dance
'''

# Import the goods
import sys
import matplotlib.pyplot as plt
import numpy as np

# File to pull data from
inTxt = "effective_particle_radius.txt"

# Import data into arrays
peA,\
peB,\
peR,\
xA,\
eps,\
modeALL,\
modeAA,\
modeAB,\
modeBB,\
fALL,\
fAA,\
fAB,\
fBB,\
phiEff = np.loadtxt(inTxt, skiprows=1, unpack=True)

# Everything is going on one plot
''' Plot data
x-axis: activity Ratio
y-axis: effective area fraction
color: activity of fast species
shape: particle fraction
'''
minALL = min(modeALL)
minAA = min(modeAA)
minAB = min(modeAB)
minBB = min(modeBB)
maxALL = max(modeALL)
maxAA = max(modeAA)
maxAB = max(modeAB)
maxBB = max(modeBB)

minVal = min(minALL, minAA, minAB, minBB)
maxVal = max(maxALL, maxAA, maxAB, maxBB)

# Let's sort the array according to PeA
for i in range(0, len(peR)):
    for j in range(0, len(peR)):
        # If peA value is smaller AND index is larger
        if peR[i] < peR[j] and i > j:
            # Store in temporary value
            tmp1 = peR[j]
            # Swap them
            peR[j] = peR[i]
            peR[i] = tmp1
            # Swap corresponding phi indices
            tmp2 = modeALL[j]
            modeALL[j] = modeALL[i]
            modeALL[i] = tmp2
            tmp2 = modeAA[j]
            modeAA[j] = modeAA[i]
            modeAA[i] = tmp2
            tmp2 = modeAB[j]
            modeAB[j] = modeAB[i]
            modeAB[i] = tmp2
            tmp2 = modeBB[j]
            modeBB[j] = modeBB[i]
            modeBB[i] = tmp2

plt.plot(peR, modeALL, c='r', zorder=1, label='All')
plt.plot(peR, modeAA, c='g', zorder=1, label='AA')
plt.plot(peR, modeAB, c='b', zorder=1, label='AB')
plt.plot(peR, modeBB, c='k', zorder=1, label='BB')
plt.scatter(peR, modeALL, c='r', zorder=2)
plt.scatter(peR, modeAA, c='g', zorder=2)
plt.scatter(peR, modeAB, c='b', zorder=2)
plt.scatter(peR, modeBB, c='k', zorder=2)
plt.legend()
plt.xlim(min(peR), max(peR))
plt.ylim(minVal, maxVal)
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.625), title=r'$(x_{slow}, Pe_{fast})$')
plt.xlabel(r'Activity Ratio $(Pe_{R})$')
plt.ylabel(r'$\sigma_{Effective}$')
#plt.tight_layout()
plt.savefig('diameter-trends.png', bbox_inches='tight', dpi=1000)
