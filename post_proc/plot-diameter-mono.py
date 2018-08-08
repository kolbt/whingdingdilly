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

# Let's sort the array according to PeA
for i in range(0, len(peA)):
    for j in range(0, len(peA)):
        # If peA value is smaller AND index is larger
        if peA[i] < peA[j] and i > j:
            # Store in temporary value
            tmp1 = peA[j]
            # Swap them
            peA[j] = peA[i]
            peA[i] = tmp1
            # Swap corresponding phi indices
            tmp2 = phiEff[j]
            phiEff[j] = phiEff[i]
            phiEff[i] = tmp2
            # Swap corresponding sigma indices
            tmp3 = modeALL[j]
            modeALL[j] = modeALL[i]
            modeALL[i] = tmp3

plt.plot(peA, phiEff, c='r', zorder=1)
plt.scatter(peA, phiEff, c='k', zorder=2)
plt.xlim(min(peA), max(peA))
plt.ylim(min(phiEff), max(phiEff))
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.625), title=r'$(x_{slow}, Pe_{fast})$')
plt.xlabel(r'Activity $(Pe)$')
plt.ylabel(r'$\phi_{Effective}$')
#plt.tight_layout()
plt.savefig('phi-trends.png', bbox_inches='tight', dpi=1000)
plt.close()

plt.plot(peA, modeALL, c='r', zorder=1)
plt.scatter(peA, modeALL, c='k', zorder=2)
plt.xlim(min(peA), max(peA))
plt.ylim(min(modeALL), max(modeALL))
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.625), title=r'$(x_{slow}, Pe_{fast})$')
plt.xlabel(r'Activity $(Pe)$')
plt.ylabel(r'$\sigma_{Effective}$')
#plt.tight_layout()
plt.savefig('diameter-trends.png', bbox_inches='tight', dpi=1000)