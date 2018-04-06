import numpy as np
import matplotlib.pyplot as plt

threeEtaPiSigma = 1.0
sigma = 1.0
epsilon = 1.0

def computeLJ(r):
    potentialLJ = 4*epsilon*(((sigma/r)**12)-((sigma/r)**6))
    return potentialLJ

def computeNegDeriv(x1, x2, y1, y2):
    negDeriv = -(y2 - y1)/(x2 - x1)
    return negDeriv


def computeLJForce(r):
    forceLJ = 48*epsilon*(((sigma**12)/(r**13))-((sigma**6)/(r**7)))
    return forceLJ

rangeR = np.arange(0.00000, 5.00000, 0.00001)
potLJ = np.zeros_like(rangeR)
forLJ = np.zeros_like(rangeR)
nDeriv = np.zeros((len(rangeR)), dtype=np.float32)

for i in range(0, len(rangeR)):
    potLJ[i] = computeLJ(rangeR[i])
    forLJ[i] = computeLJForce(rangeR[i])

myEvals = np.zeros((4,2), dtype=np.float32)
for i in range(1, len(rangeR)):
    nDeriv[i-1] = computeNegDeriv(rangeR[i-1], rangeR[i], potLJ[i-1], potLJ[i])
    if (999 < nDeriv[i-1] < 1001):
        myEvals[0][0] = rangeR[i]
        myEvals[0][1] = nDeriv[i-1]
    if (499 < nDeriv[i-1] < 501):
        myEvals[1][0] = rangeR[i]
        myEvals[1][1] = nDeriv[i-1]
    if (rangeR[i] == 1.00000):
        myEvals[2][0] = rangeR[i]
        myEvals[2][1] = nDeriv[i-1]
    if (nDeriv[i-1] < 0) and (nDeriv[i-2] > 0):
        myEvals[3][0] = rangeR[i]
        myEvals[3][1] = nDeriv[i-1]


print(myEvals)

plt.plot(rangeR, potLJ)
plt.xlim(0,5)
plt.ylim(-1.1,2)
plt.xlabel(r'$\frac{r}{\sigma}$')
plt.ylabel(r'$U_{Lennard-Jones}$')
plt.show()

#plt.plot(rangeR, forLJ)
plt.plot(rangeR, nDeriv)

plt.scatter(myEvals[0][0], myEvals[0][1])
plt.text(myEvals[0][0] + 0.01, myEvals[0][1] - 25,
         ("{0:.2f}".format(myEvals[0][0]), "{0:.0f}".format(myEvals[0][1])))

plt.scatter(myEvals[1][0], myEvals[1][1])
plt.text(myEvals[1][0] + 0.01, myEvals[1][1],
         ("{0:.2f}".format(myEvals[1][0]), "{0:.0f}".format(myEvals[1][1])))

plt.scatter(myEvals[2][0], myEvals[2][1])
plt.text(myEvals[2][0] - 0.03, myEvals[2][1] + 25,
         ("{0:.2f}".format(myEvals[2][0]), "{0:.0f}".format(myEvals[2][1])))

plt.scatter(myEvals[3][0], myEvals[3][1])
plt.text(myEvals[3][0] - 0.045, myEvals[3][1] + 25,
         ("{0:.2f}".format(myEvals[3][0]), "{0:.0f}".format(myEvals[3][1])))

plt.xlim(0.77,1.14)
plt.ylim(-1,1010)
plt.xlabel(r'$\frac{r}{\sigma}$')
plt.ylabel(r'$F_{Lennard-Jones}$')
plt.show()

#def effectiveSigma(Fp):
#    "Read in the force, compute effective diameter (r)"
#    
#    return effSigma



# Use this to compute the effective particle diameter
#def effectiveSigma(Fp, epsilon):
#    effSigma = epsilon / Fp
#    return effSigma
#
#Pe = np.arange(1.0, 500.0, 1.0)
#effSig = np.zeros_like(Pe)
#for i in range(0, len(Pe)):
#    effSig[i] = effectiveSigma(Pe[i], 1.0)
#
#plt.plot(Pe, effSig)
#plt.show()
