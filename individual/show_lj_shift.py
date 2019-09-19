import numpy as np
import matplotlib.pyplot as plt

def ljPotential(r, sigma, epsilon):
    u = 4. * epsilon * ( ((sigma/r)**12) - ((sigma/r)**6) )
    return u

def ljForce(r, sigma, epsilon):
    f = 24. * epsilon * ( (2*((sigma**12)/(r**13))) - ((sigma**6)/(r**7)) )
    return f

def compEps(pe):
    epsB = 10.
    eps = ((4. * pe) / 24.) + epsB
    return eps

def compSig(pe):
    m = -0.095771
    b = 0.307
    sig = (pe**m) * np.exp(b)
    return sig

# Get input values for my functions
rs = np.arange(0., 2., 0.0001)
pes = np.arange(50., 550., 50.)
rSoft = 2**(1./6.)

fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Plot the soft data
softPots = ljPotential(rs, 1.0, 1.0)
softForce = ljForce(rs, 1.0, 1.0)
ax[0].plot(rs, softPots, c='k')
ax[1].plot(rs, softForce, c='k')

for i in pes:
    # Grab epsilon and sigma
    inSig = compSig(i)
    inEps = compEps(i)
    # Compute the potentials and forces for each system
    hardPots = ljPotential(rs, inSig, inEps)
    hardForce = ljForce(rs, inSig, inEps)
    ax[0].plot(rs, hardPots, c=plt.cm.jet(i / 500.))
    ax[1].plot(rs, hardForce, c=plt.cm.jet(i / 500.), label=str(i))

potMax = ljPotential(0.7, compSig(500.), compEps(500.))
forMax = ljForce(0.7, compSig(500.), compEps(500.))

ax[0].set_xlim(0.7, 2.0)
ax[0].set_ylim(-inEps * 1.1, potMax)
ax[1].set_xlim(0.7, rSoft)
ax[1].set_ylim(0., forMax)

# Is there any significance to the crossing point?
#ax[0].axvline(x=(compSig(100) * (2.**(1./6.)) ))
#ax[1].axvline(x=(compSig(100) * (2.**(1./6.)) ))
#ax[0].axvline(x=(compSig(100) ))
#ax[1].axvline(x=(compSig(100) ))

ax[0].set_ylabel(r'$U_{LJ}$')
ax[1].set_ylabel(r'$F_{LJ}$')
ax[0].set_xlabel(r'$r(i,j)$')
ax[1].set_xlabel(r'$r(i,j)$')
cbar = ax[1].legend()
cbar.set_title('Activity')
plt.tight_layout()
plt.savefig('smallHS.png', dpi=1000)
