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
''' Plot stats
x-axis: activity Ratio
y-axis: effective area fraction
color: activity of fast species
shape: particle fraction
'''

plt.scatter(peR, phiEff)
plt.ylim(0.45, 0.60)
plt.xlabel(r'$Pe_{Ratio}$')
plt.ylabel(r'$\phi_{Effective}$')
plt.show()






