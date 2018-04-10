import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

fActive = 500.0
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

#print(computeEffectiveSigma(0))
#print(computeEffectiveSigma(fActive))
#print(computeEffectiveSigma(2*fActive))

mylab = r'$\epsilon=\frac{2F_{A}}{24.0}$'
for j in range(0,3):
    if j==1:
        epsilon = fActive/24.0
        mylab = r'$\epsilon=\frac{F_{A}}{24.0}$'
    if j==2:
        epsilon = 1.0
        mylab = r'$\epsilon=1$'
    
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

    #plt.plot(distRange, forLJ)
    #plt.xlim(0.77,1.14)
    #plt.ylim(-1,1010)
    #plt.tight_layout()
    #plt.show()
    #
    #plt.plot(forcRange, efSig)
    fA1 = fActive
    eff1 = computeEffectiveSigma(fA1)[0]
    fA2 = 2*fActive
    eff2 = computeEffectiveSigma(fA2)[0]

    plt.plot(forLJ, distRange, label=mylab)
    plt.scatter(fA1, eff1, c='k')
    plt.text(fA1 - 90, eff1 - 0.021, ("{0:.0f}".format(fA1), "{0:.3f}".format(eff1)))
    plt.scatter(fA2, eff2, c='g')
    plt.text(fA2 - 215, eff2 - 0.015, ("{0:.0f}".format(fA2), "{0:.3f}".format(eff2)))

    plt.xlim(0, 1000)
    plt.ylim(0.76, 1.122)
    plt.xlabel('Force', fontsize=20)
    plt.ylabel('Effective Diameter', fontsize=20)
    plt.tight_layout()

plt.legend(loc=3, fontsize=20)
plt.savefig('various_epsilons.png', dpi=1000)
