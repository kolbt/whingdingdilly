import numpy as np
import matplotlib.pyplot as plt

threeEtaPiSigma = 1.0
epsilon = 1.0

# Use this to compute the effective particle diameter
def effectiveSigma(Fp, epsilon):
    effSigma = epsilon / Fp
    return effSigma

Pe = np.arange(1.0, 500.0, 1.0)
effSig = np.zeros_like(Pe)
for i in range(0, len(Pe)):
    effSig[i] = effectiveSigma(Pe[i], 1.0)

plt.plot(Pe, effSig)
plt.show()
