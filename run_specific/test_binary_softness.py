import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

fActive = 24.0
threeEtaPiSigma = 1.0
sigma = 1.0
epsilon = 2.0*fActive / 24.0

def computeLJPotential(r):
    "Given a distance, computes the potential"
    potentialLJ = 4*epsilon*(((sigma/r)**12)-((sigma/r)**6))
    return potentialLJ

def computeLJForce(r):
    "Given a distance, computes the force"
    forceLJ = 4*epsilon*((12*(sigma**12)*(r**-13))-(6*(sigma**12)*(r**-7)))
    return forceLJ

def functionToSolve(r, f):
    return (2*(sigma**6)*(r**-13))-(r**-7)-(f/(24*epsilon*(sigma**6)))

def computeEffectiveSigma(f):
    "Given a force, give the center-to-center distance"
    r_guess = 0.8
    effectiveSigma = fsolve(functionToSolve, r_guess, args=f)
    return effectiveSigma


print(computeEffectiveSigma(0))
print(computeEffectiveSigma(fActive))
print(computeEffectiveSigma(2*fActive))

distRange = np.arange(0.00000, 5.00000, 0.00001)
forcRange = np.arange(0.00, 1000.00, 0.01)
potLJ = np.zeros_like(distRange)
forLJ = np.zeros_like(distRange)
efSig = np.zeros_like(forcRange)

for i in range(0, len(distRange)):
    potLJ[i] = computeLJPotential(distRange[i])
    forLJ[i] = computeLJForce(distRange[i])
    
for i in range(0, len(forcRange)):
    efSig[i] = computeEffectiveSigma(forcRange[i])

plt.plot(distRange, forLJ)
plt.xlim(0.77,1.14)
plt.ylim(-1,1010)
plt.tight_layout()
plt.show()

plt.plot(forcRange, efSig)
plt.plot(forLJ, distRange)
plt.xlim(0, 1000)
plt.ylim(0.7, 1.0)
plt.tight_layout()
plt.show()
