import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#################################################################################
### Our theory function... ######################################################
# increasing the gf cutoff gives fewer phase separated pixels ###################
def evalTheory( peA, peB, xA, gf ):
    if ((xA * peA) - (xA * peB) + peB) > ((3 * ((np.pi)**2) * kappa)/(4 * phi)):
        if (3 * ((np.pi)**2) * kappa)/(4 * phi * (xA * (peA - peB) + peB)) < gf:
            return 1
        else:
            return 0    # if there is greater than gf% gas it is a gas
    else:
        return 0
#################################################################################

# read in your arguments first
file = str(sys.argv[1])
pb = float(sys.argv[2])
gf = float(sys.argv[3])
gf /= 100
kappa = float(sys.argv[4])
kappa /= 1000

# theory stuff is gonna go up top
phi = 0.6
#kappa = 2.225

# array to hold theory values
theory_diag = np.zeros((11,16), dtype=np.int8)
fvar = 1.0                      # variable for xa
pvar = 0                        # variable for peA
for i in range(0, 16):
    for j in range(0, 11):
        theory_diag[j][i] = evalTheory(pvar, pb, fvar, gf)
        fvar -= 0.1
    fvar = 1.0
    pvar += 10

# start with 1D array
data = np.loadtxt(file, dtype=np.int8)

# change this array so that it is 11(rows) x 16(columns)
phase_diag = np.zeros((11,16), dtype=np.int8)

# for loop the fuck outta this
r=10
c=0
for i in range(0,len(data)):
    phase_diag[r][c] = data[i]
    r -= 1
    if r == -1:
        r = 10
        c += 1

# compare your arrays to see how well we did!
count = 0
for i in range(0, 16):
    for j in range(0, 11):
        if theory_diag[j][i] != phase_diag[j][i]:
            count += 1
# this gives the number of mismatches
#print(count)
sys.exit(count)

# WORD. Now plot that shit
#plt.imshow(theory_diag, cmap='Blues')
#plt.savefig('py_theory_test.png', dpi=200)
#plt.close()
#plt.imshow(phase_diag, cmap='Blues')
#plt.savefig('myphase_'+str(pb)+'.png', dpi=200)
#plt.close()
