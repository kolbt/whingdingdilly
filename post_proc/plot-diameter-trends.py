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

# Get distinct PeR values
peRs = []
for i in range(0, len(peR)):
    if not (peR[i] in peRs):
        peRs.append(peR[i])

# Get distinct PeFast values
peFast = []
for i in range(0, len(peB)):
    if not (peB[i] in peFast):
        peFast.append(peB[i])
# colors = np.arange(0, len(peFast))

# Get distinct xSlow values
xSlow = []
for i in range(0, len(xA)):
    if not (xA[i] in xSlow):
        xSlow.append(xA[i])
# shapes = np.arange(0, len(xSlow))

# # Loop through PeR
# for h in range(0, len(peR)):
#     # Get index for PeFast/ set color
#     color = 0
#     for i in range (0, len(peFast)):
#         if peB[h] == peFast[i]:
#             color = i
#     # Get index for xSlow/ set shape
#     shape = 0
#     for i in range(0, len(xSlow)):
#         if xA[h] == xSlow[i]:
#             shape = i
#     plt.scatter(peR[h], phiEff[h], c=color, s=shape)

# Need to connect lines that have same PeFast and xSlow
lines = np.zeros( (len(peFast) * len(xSlow), len(peRs) ), dtype=np.float64)

# Loop through all data points
for i in range(0, phiEff):
    for i in range(0, )



# #for loop for color
# for i in range(0, len(peFast)):
#     #for loop for shape
#     for j in range(0, len(xSlow)):
#         #plot the point
#         plt.scatter(peR[g], phiEff[h], c=colors[i], s=shapes[j])

plt.ylim(0.45, 0.60)
plt.xlabel(r'$Pe_{Ratio}$')
plt.ylabel(r'$\phi_{Effective}$')
plt.show()






