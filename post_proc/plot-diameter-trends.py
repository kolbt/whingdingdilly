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
import colorsys

# Color generator
def getColorScheme(numColors):
    '''Generate a colorscheme'''
    colors = []
    r = 0
    g = 0
    b = 0
    step = float(1.0 / numColors)
    for i in range(0, numColors):
        colors.append((r, g, b))
        r += step
        g += step
        b += step
    return colors

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
peRs = sorted(peRs)

# Get distinct PeFast values
peFast = []
for i in range(0, len(peB)):
    if not (peB[i] in peFast):
        peFast.append(peB[i])
peFast = sorted(peFast)

# Get distinct xSlow values
xSlow = []
for i in range(0, len(xA)):
    if not (xA[i] in xSlow):
        xSlow.append(xA[i])
xSlow = sorted(xSlow)

numLines = (len(xSlow) * len(peFast))
linepeR = [ [] for x in xrange(0, numLines) ]
linephiEff = [ [] for x in xrange(0, numLines) ]
key = np.zeros(numLines, dtype=np.ndarray)

dataNum = len(phiEff)
count = 0
for i in range(0, len(xSlow)):
    constxSlow = xSlow[i]
    for j in range(0, len(peFast)):
        constpeFast = peFast[j]
        key[count] = (constxSlow, constpeFast)
        for k in range(0, len(peRs)):
            tmpPeR = peRs[k]
            for l in range(0, dataNum):
                if xA[l] == constxSlow and peB[l] == constpeFast and peR[l] == tmpPeR:
                    linepeR[count].append(peR[l])
                    linephiEff[count].append(phiEff[l])
        count += 1

# Get colorscheme of appropriate length
colorNum = len(peFast)
color = getColorScheme(colorNum)
# Initialize symbol list
symbols =['o', 'v', 's', 'D', '*', 'h', 'p', '8', '+']
# Get shift list
shift = np.zeros(len(xSlow), dtype=np.float64)
step = 0.015
shift -= step * 0.5 * len(xSlow)
for i in range(0, len(shift)):
    shift[i] += step * i

for i in range(0, numLines):
    myLabel = "({}, {})".format(key[i][0]/100.0, key[i][1])
    # if for color (pefast)
    for j in range(0, len(peFast)):
        if key[i][1] == peFast[j]:
            cInd = j
    # if for symbol (xSlow)
    for j in range(0, len(xSlow)):
        if key[i][0] == xSlow[j]:
            sInd = j
    # Plot it
    plt.plot(linepeR[i] + shift[sInd], linephiEff[i], c=color[cInd], marker=symbols[sInd], label=myLabel)

#plt.ylim(0.45, 0.60)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.625), title=r'$(x_{slow}, Pe_{fast})$')
plt.xlabel(r'$Pe_{Ratio}$')
plt.ylabel(r'$\phi_{Effective}$')
#plt.tight_layout()
plt.savefig('diameter-trends.png', bbox_inches='tight', dpi=1000)






