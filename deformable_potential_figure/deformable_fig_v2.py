'''
A schematic figure to illustrate how Pe (and F_act) sets deformability:
    -LJ potential
    -Overlay 2 Forces (Strong under week)
    -Corresponds to collision angle
'''

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from PIL import Image
import matplotlib.patches as patches
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
fsize = 9
plt.rcParams.update({'font.size': fsize})
params = {'legend.fontsize': fsize,
          'axes.labelsize': fsize,
          'axes.titlesize': fsize,
          'xtick.labelsize': fsize,
          'ytick.labelsize': fsize}
plt.rcParams.update(params)

eps = 0.1
sigma = 1.
def ljPotential(r, eps=0.1, sigma=1.):
    div = (sigma/r)
    U = ( 4. * eps * ((div)**12 - (div)**6) ) + eps
    return U

def ljForce(r, eps=0.1, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
    
def convergeConstPeEps(pe, eps):
    r = 1.112
    while ljForce(r, eps) < pe:
        r -= 0.0001
    return r

# Compute the weak and strong collision force
peWeak = 10.
peMid = 50.
peStrong = 150.
rWeak = convergeConstPeEps(peWeak, eps)
rMid = convergeConstPeEps(peMid, eps)
rStrong = convergeConstPeEps(peStrong, eps)

# Plot the figure
width = 3. + (3./8.)        # single-column figure width
fig = plt.figure(figsize=(width, width))
ax = []
#ax.append(fig.add_subplot(131))                     # left column
#ax.append(fig.add_subplot(332, projection='3d'))    # right top
#ax.append(fig.add_subplot(324, projection='3d'))    # right mid
#ax.append(fig.add_subplot(326, projection='3d'))    # right bottom

ax.append(fig.add_subplot(111))                                 # potential figure
sWidth = 0.3
sHeight = sWidth * (1./2.)
ax.append(plt.axes([0.55, 0.26, sWidth, sHeight], projection='3d'))   # weak pe spheres
ax.append(plt.axes([0.325, 0.425, sWidth, sHeight], projection='3d'))  # mid pe spheres
ax.append(plt.axes([0.225, 0.675, sWidth, sHeight], projection='3d'))   # strong pe spheres
ax[1].patch.set_alpha(0.0)
ax[2].patch.set_alpha(0.0)
ax[3].patch.set_alpha(0.0)

# Plot LJ potential
dist = np.arange(0.0001, ((2.**(1./6.))*sigma)*2., 0.001)
ax[0].plot(dist, ljPotential(dist, eps=eps), c='k', lw=1.5, label='LJ-Potential', zorder=0)

# Overlay the weak activity range
shift = 1.
weakRange = np.arange(rWeak, ((2.**(1./6.))*sigma)*2., 0.001)
ax[0].plot(weakRange, ljPotential(weakRange, eps=eps) + (1.5 * shift), c='r', lw=1.25, label='Weak', zorder=0)
ax[0].scatter(rWeak, ljPotential(rWeak, eps=eps) + (1.5 * shift), c='r', zorder=1)

# Overlay the middle activity range
midRange = np.arange(rMid, ((2.**(1./6.))*sigma)*2., 0.001)
ax[0].plot(midRange, ljPotential(midRange, eps=eps) + (1. * shift), c='b', lw=1.25, label='Mid', zorder=0)
ax[0].scatter(rMid, ljPotential(rMid, eps=eps) + (1. * shift), c='b', zorder=1)

# Overlay the strong activity range
strongRange = np.arange(rStrong, ((2.**(1./6.))*sigma)*2., 0.001)
ax[0].plot(strongRange, ljPotential(strongRange, eps=eps) + (0.5 * shift), c='g', lw=1.25, label='Strong', zorder=0)
ax[0].scatter(rStrong, ljPotential(rStrong, eps=eps) + (0.5 * shift), c='g', zorder=1)

# Limits
ax[0].set_xlim(rStrong - 0.05, (2.**(1./6.))*sigma)
ax[0].set_ylim(0., 10.)
ax[0].set_xlabel(r'Interparticle distance $(\delta_{i,j})$', fontsize=fsize)
ax[0].set_ylabel(r'Lennard-Jones potential $(U_{LJ})$', fontsize=fsize)
#ax[0].legend()

# Plot the overlap of spheres
# For wire mesh
backu, backv = np.mgrid[1*np.pi:2*np.pi:10j, 0:np.pi:10j]
backx = np.cos(backu)*np.sin(backv)
backy = np.sin(backu)*np.sin(backv)
backz = np.cos(backv)
frontu, frontv = np.mgrid[0*np.pi:1*np.pi:10j, 0:np.pi:10j]
frontx = np.cos(frontu)*np.sin(frontv)
fronty = np.sin(frontu)*np.sin(frontv)
frontz = np.cos(frontv)
# For solid sphere
uS, vS = np.mgrid[0:2*np.pi:1000j, 0:np.pi:500j]
xS = np.cos(uS)*np.sin(vS)
yS = np.sin(uS)*np.sin(vS)
zS = np.cos(vS)
backAlph = 0.3
frontAlph = 0.5
ax[1].plot_wireframe(backx - rWeak, backy, backz, color="#808080", alpha=backAlph)
ax[1].plot_wireframe(backx + rWeak, backy, backz, color="#808080", alpha=backAlph)
ax[1].plot_surface((xS*rWeak) - rWeak, yS*rWeak, zS*rWeak, color="r")
ax[1].plot_surface((xS*rWeak) + rWeak, yS*rWeak, zS*rWeak, color="r")
ax[1].plot_wireframe(frontx - rWeak, fronty, frontz, color="#808080", alpha=frontAlph)
ax[1].plot_wireframe(frontx + rWeak, fronty, frontz, color="#808080", alpha=frontAlph)
ax[1].set_axis_off()
ax[1].view_init(0, 90)
ax[1].set_xlim(-2., 2.)
ax[1].set_ylim(-1., 1.)
ax[1].set_zlim(-1., 1.)
ax[1].dist = 6.

ax[2].plot_wireframe(backx - rMid, backy, backz, color="#808080", alpha=backAlph)
ax[2].plot_wireframe(backx + rMid, backy, backz, color="#808080", alpha=backAlph)
ax[2].plot_surface((xS*rMid) - rMid, yS*rMid, zS*rMid, color="b")
ax[2].plot_surface((xS*rMid) + rMid, yS*rMid, zS*rMid, color="b")
ax[2].plot_wireframe(frontx - rMid, fronty, frontz, color="#808080", alpha=frontAlph)
ax[2].plot_wireframe(frontx + rMid, fronty, frontz, color="#808080", alpha=frontAlph)
ax[2].set_axis_off()
ax[2].view_init(0, 90)
ax[2].set_xlim(-2., 2.)
ax[2].set_ylim(-1., 1.)
ax[2].set_zlim(-1., 1.)
ax[2].dist = 6.

ax[3].plot_wireframe(backx - rStrong, backy, backz, color="#808080", alpha=backAlph)
ax[3].plot_wireframe(backx + rStrong, backy, backz, color="#808080", alpha=backAlph)
ax[3].plot_surface((xS*rStrong) - rStrong, yS*rStrong, zS*rStrong, color="g")
ax[3].plot_surface((xS*rStrong) + rStrong, yS*rStrong, zS*rStrong, color="g")
ax[3].plot_wireframe(frontx - rStrong, fronty, frontz, color="#808080", alpha=frontAlph)
ax[3].plot_wireframe(frontx + rStrong, fronty, frontz, color="#808080", alpha=frontAlph)
ax[3].set_axis_off()
ax[3].view_init(0, 90)
ax[3].set_xlim(-2., 2.)
#ax[3].set_ylim(-1.5, 1.5)
#ax[3].set_zlim(-1.5, 1.5)
ax[3].set_ylim(-1., 1.)
ax[3].set_zlim(-1., 1.)
ax[3].dist = 6.

ax[0].text(0.75, 0.75, r'$Pe=$'+"{0:g}".format(peWeak), color='r', transform=ax[0].transAxes, fontsize=fsize)
ax[0].text(0.75, 0.825, r'$Pe=$'+"{0:g}".format(peMid), color='b', transform=ax[0].transAxes, fontsize=fsize)
ax[0].text(0.75, 0.9, r'$Pe=$'+"{0:g}".format(peStrong), color='g', transform=ax[0].transAxes, fontsize=fsize)

# Set tick parameters
ax[0].tick_params(axis='both', direction='in', labelsize=fsize)

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

ax[0].xaxis.set_minor_locator(MultipleLocator(0.05))
ax[0].yaxis.set_minor_locator(MultipleLocator(1))
ax[0].tick_params(axis='both', which='minor', length=2, direction='in')

labels = ax[0].get_xticks().tolist()
print(labels)
for i in range(0, len(labels)):
    if labels[i] == 1.0:
        labels[i] = r'$\sigma$'
print(labels)
ax[0].set_xticklabels(labels)

# Let's add spatial heatmaps to this figure
ld_img = [] # list to hold images
# Path to images
imPath = '/Users/kolbt/Desktop/soft_figures/method_schematic'
imPath = '/Users/kolbt/Desktop/compiled/whingdingdilly/ipython/clusters_soft'
# Image file names
imgs = ['spatial_delta_pa10.0_pb0_xa100.0_frame0600.png',
        'spatial_delta_pa50.0_pb0_xa100.0_frame0600.png',
        'spatial_delta_pa150.0_pb0_xa100.0_frame0600.png']
# The height/width
dim = 6600
for i in imgs:
    im = Image.open(imPath + '/' + i)
    # (left, upper, right, lower)
    # image is 9600 x 7200
    left = 700
    upper = 250
    im1 = im.crop((left, upper, left+dim, upper+dim))
    ld_img.append(im1)
    
# Add an axes to the right of the plot for the heatmaps
imdim = 0.25
base = 0.1125
buff = 0.008
left = 0.92
bottom, width, height = base, imdim, imdim
hm1 = fig.add_axes([left, bottom, width, height])
hm1.imshow(ld_img[0])
hm1.set_axis_off()
hm1.set_aspect('equal')

bottom, width, height = base + imdim + buff, imdim, imdim
hm2 = fig.add_axes([left, bottom, width, height])
hm2.imshow(ld_img[1])
hm2.set_axis_off()
hm2.set_aspect('equal')

bottom, width, height = base + (2.*imdim) + (2.*buff), imdim, imdim
hm3 = fig.add_axes([left, bottom, width, height])
hm3.imshow(ld_img[2])
hm3.set_axis_off()
hm3.set_aspect('equal')

hm1.add_patch(patches.Rectangle((0, 0), 1, 1, linewidth=2.5, edgecolor='r', facecolor='none', transform=hm1.transAxes))
hm2.add_patch(patches.Rectangle((0, 0), 1, 1, linewidth=2.5, edgecolor='b', facecolor='none', transform=hm2.transAxes))
hm3.add_patch(patches.Rectangle((0, 0), 1, 1, linewidth=2.5, edgecolor='g', facecolor='none', transform=hm3.transAxes))

cbax = fig.add_axes([0.55, base, 0.765, 0.765])
divider = make_axes_locatable(cbax)
cax = divider.append_axes("left", size="1%", pad=0.0)
cmap = mpl.cm.jet_r
norm = mpl.colors.Normalize(vmin=0.6, vmax=1.0)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cb1 = fig.colorbar(sm, ax=cax, orientation='vertical')
cb1.set_label(r'$\delta$', fontsize=fsize)
cbax.axis('off')
cax.axis('off')
    
#plt.savefig("particle_deformation_eps" + str(eps) + ".png", dpi=2000, bbox_inches='tight', pad_inches=0)
plt.savefig("particle_deformation_eps" + str(eps) + ".png", dpi=2000, bbox_inches="tight")
plt.close()
