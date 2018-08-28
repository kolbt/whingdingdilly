'''
#                           This is an 80 character line                       #
This file reads in processed data from text files and outputs whether or not the
given system undergoes phase separation.  This is currently computed by the size
and lifetime of the largest cluster, but, could also feasibly arise from the
percentage of clustered particles in the system (not necessarily one cluster).

This file:
    1.  Read in SA txt (l_clust and phi)
    2.  Compute effective activity
    3.  Use Cates formula to check agreement
    4.  Plot Pe_intended vs x_A
    5.  Use filled unfilled to indicate agreement
'''

# Imports and loading the .gsd file
import sys
import numpy as np
import math
import matplotlib.pyplot as plt

# Function that'll grab my parameters from the filenames
def getFromTxt(fname, first, last):
    """Takes a string, text before and after desired text, outs text between"""
    start = fname.index( first ) + len( first )
    end = fname.index( last, start )
    myTxt = fname[start:end]
    return float(myTxt)

sigma = 1.0
# The old model
epsilon = 1.0

# Grab the command line arguments
txtFiles = sys.argv[2:]     # start at 2 to avoid early arguments

# Parse out the activities and fractions from the filenames
peA = np.zeros_like(txtFiles, dtype=np.int)
peB = np.zeros_like(txtFiles, dtype=np.int)
xA = np.zeros_like(txtFiles, dtype=np.float64)
ep = np.zeros_like(txtFiles, dtype=np.int)

try:
    for i in range(0, len(txtFiles)):
        peA[i] = getFromTxt(txtFiles[i], "pa", "_pb")
        peB[i] = getFromTxt(txtFiles[i], "pb", "_xa")
        xA[i] = getFromTxt(txtFiles[i], "xa", "_ep")
        ep[i] = getFromTxt(txtFiles[i], "ep", ".txt")
except:
    for i in range(0, len(txtFiles)):
        peA[i] = getFromTxt(txtFiles[i], "pa", "_pb")
        peB[i] = getFromTxt(txtFiles[i], "pb", "_xa")
        xA[i] = getFromTxt(txtFiles[i], "xa", ".txt")
        ep[i] = 1

try:
    peR = peA.astype(float) / peB.astype(float)         # Compute activity ratio
except:
    peR = np.zeros(len(txtFiles))

# see file format in write-phase-txt.py
phase_file = "phase_sep_w_diameters.txt"
f = open(phase_file, 'w')
f.write(('Act_A').center(10) + ' ' + \
        ('Act_B').center(10) + ' ' + \
        ('Frac_A').center(10) + ' ' + \
        ('Sig_All').center(10) + ' ' + \
        ('Sig_AA').center(10) + ' ' + \
        ('Sig_AB').center(10) + ' ' + \
        ('Sig_BB').center(10) + ' ' + \
        ('Phi_Eff').center(10) + ' ' + \
        ('Phase_Sep').center(10) + \
        '\n')
f.close()

# Loop through each data series
for i in range(0, len(txtFiles)):

    # Import data into arrays
    tst, \
    gasA, \
    gasB, \
    gasTot, \
    denA, \
    denB, \
    denTot, \
    lgClust, \
    MCS, \
    sigALL, \
    sigAA, \
    sigAB, \
    sigBB, \
    phiEff, \
    lC_Area, \
    totC_Area, \
    lC_density, \
    denDen, \
    gasDen = np.loadtxt(txtFiles[i], skiprows=1, unpack=True)

    # Requirement to be consider phase separated
    partNum = gasTot[0]     # everything starts in a gas
    frames = len(tst)
    sizeMin = partNum * 0.25  # 40% of particles in single cluster
    timeMin = frames * 0.50  # cluster present for half of all frames

    # See if data satisfies phase separated requirements
    count = 0
    phaseSep = 0
    ALL = 0
    AA = 0
    AB = 0
    BB = 0
    phiAvg = 0
    # Get last 10% of simulation
    numAvg = (0.10 * len(lgClust))
    avgTime = len(lgClust) - numAvg

    for j in range(0, len(lgClust)):
        # Average over last
        if j >= avgTime:
            ALL += sigALL[j]
            AA += sigAA[j]
            AB += sigAB[j]
            BB += sigBB[j]
            phiAvg += phiEff[j]
        if lgClust[j] >= sizeMin:
            count += 1

    # Average diameter values
    ALL /= numAvg
    AA /= numAvg
    AB /= numAvg
    BB /= numAvg
    phiAvg /= numAvg

    if count >= timeMin:
        phaseSep = 1
        f = open(phase_file, 'a')
        f.write(str(peA[i]).center(10) + ' ' + \
                str(peB[i]).center(10) + ' ' + \
                str(xA[i]).center(10) + ' ' + \
                '{0:.4f}'.format(ALL).center(10) + ' ' + \
                '{0:.4f}'.format(AA).center(10) + ' ' + \
                '{0:.4f}'.format(AB).center(10) + ' ' + \
                '{0:.4f}'.format(BB).center(10) + ' ' + \
                '{0:.3f}'.format(phiAvg).center(10) + ' ' + \
                str(phaseSep).center(10) + \
                '\n')
    else:
        f = open(phase_file, 'a')
        f.write(str(peA[i]).center(10) + ' ' + \
                str(peB[i]).center(10) + ' ' + \
                str(xA[i]).center(10) + ' ' + \
                '{0:.4f}'.format(ALL).center(10) + ' ' + \
                '{0:.4f}'.format(AA).center(10) + ' ' + \
                '{0:.4f}'.format(AB).center(10) + ' ' + \
                '{0:.4f}'.format(BB).center(10) + ' ' + \
                '{0:.3f}'.format(phiAvg).center(10) + ' ' + \
                str(phaseSep).center(10) + \
                '\n')
    f.close()

# At this point we have a textfile with params, diameters, and ps

# Now let's take that data and compare to Cates theory
kappa = 4.05

def effectiveActivity(peIntended, sigEff):
    '''Compute the effective activity'''
    peEffective = peIntended * (sigEff**4)
    return peEffective

def catesTheory(phiEff, peEff):
    '''This uses the Cates theory to provide the minimum active
        fraction of particles required for phase separation'''
    minFrac = (3 * (np.pi ** 2) * kappa) / (4 * phiEff * peEff)
    return minFrac

# Title plot by effective Pe being used (plot using intended Pe)
peA, \
peB, \
xA, \
ALL, \
AA, \
AB, \
BB, \
phiEff, \
phaseSep = np.loadtxt(phase_file, skiprows=1, unpack=True)

# Extract monodisperse values for diameter and phi
distPeAs = []
resultSig = []
resultPhi = []
for i in range(0, len(peA)):
    if xA[i] == 100.0:
        distPeAs.append(peA[i])
        resultSig.append(ALL[i])
        resultPhi.append(phiEff[i])

print("Activities: {}").format(distPeAs)
print("Corresponding Sigmas: {}").format(resultSig)
print("Corresponding Phis: {}").format(resultPhi)

xA /= 100.0

for i in range(0, len(peA)):

    # Grab effective phi and sigma, compute Cates theory
    for j in range(0, len(distPeAs)):
        if peA[i] == distPeAs[j]:
            myPeEff = effectiveActivity(peA[i], resultSig[j])
            minFracEff = catesTheory(resultPhi[j], myPeEff)

    # Intended theory comparison
    minFracInt = catesTheory(0.60, peA[i])

    # Plot the intended and effective theory
    plt.scatter(peA[i], minFracInt, facecolor='g')
    plt.scatter(peA[i], minFracEff, facecolor='c')

    if xA[i] >= minFracEff:     # is phase separated (theory)
        if phaseSep[i] == 1:    # agreement!
            psAgree = plt.scatter(peA[i], xA[i], facecolor='k', edgecolor='k')
        else:                   # theory: yes, simulation: no
            psTheory = plt.scatter(peA[i], xA[i], facecolor='r', edgecolor='k')

    else:                       # is gas (theory)
        if phaseSep[i] == 0:    # agreement!
            gasAgree = plt.scatter(peA[i], xA[i], facecolor='w', edgecolor='k')
        else:                   # theory: no, simulation: yes
            gasTheory = plt.scatter(peA[i], xA[i], facecolor='b', edgecolor='k')

# Ensure each point type exists
psAgree = plt.scatter(-10, -10, facecolor='k', edgecolor='k')
psTheory = plt.scatter(-10, -10, facecolor='r', edgecolor='k')
gasAgree = plt.scatter(-10, -10, facecolor='w', edgecolor='k')
gasTheory = plt.scatter(-10, -10, facecolor='b', edgecolor='k')

plt.xlim(min(peA), max(peA))
plt.ylim(min(xA), max(xA))
legend = plt.legend(bbox_to_anchor=(0.0, 1.0),
                    loc='center left',
                    ncol=4,
                    handles=[psAgree, psTheory, gasAgree, gasTheory],
                    labels=['PS Both', 'PS Theory Only', 'Gas Both', 'Gas Theory Only'],
                    framealpha=1.0)
plt.title(r'Active/Passive Phase Separation', y=1.03)
plt.xlabel(r'$Pe_{Fast}$')
plt.ylabel(r'$x_{Fast}$')
plt.tight_layout()
plt.savefig('renormalize_params.png', dpi=1000)
plt.close()
