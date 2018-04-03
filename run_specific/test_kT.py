import numpy as np
import matplotlib.pyplot as plt

threeEtaPiSigma = 1.0
epsilon = 1.0
sigma = 1.0
tau_lj = ((sigma**2)*threeEtaPiSigma) / (epsilon)
tau_lj = 1
F_a = 24 * epsilon / sigma
velocity = F_a * threeEtaPiSigma

def computekT(activity, velocity, sigma, threeEtaPiSigma):
    "This takes the activity I want and outputs the corresponding temperature"
    kT = float((velocity * sigma * threeEtaPiSigma) / activity)
    return kT

# This gives a temperature of 298, almost 0 activity
print(computekT(0.0805, velocity, sigma, threeEtaPiSigma))
print((threeEtaPiSigma*(sigma**2)))

evals = np.arange(10,500,0.01)
outs = np.zeros_like(evals)
for i in range(0, len(evals)):
    outs[i] = computekT(evals[i], velocity, sigma, threeEtaPiSigma)

plt.plot(evals, outs)
plt.show()
