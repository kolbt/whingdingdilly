'''
#                           This is an 80 character line                       #
Read in .gsd file, output a 3 x 2 png at each timestep:
    -original
    -orientation (colorbar)
    -magnitude of direcional active force (heatmap)
    -tracking individual particles
    -particle radii (heatmap)
    -nearest neighbor (heatmap)
'''

import sys

# Run locally
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')
sys.path.append('/Users/kolbt/Desktop/compiled/gsd/build')
# Run on the cpu
sys.path.append('/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build')
# Run on the gpu
sys.path.append('/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build')
sys.path.append('/nas/longleaf/home/kolbt/programs/gsd/build')

import hoomd
from hoomd import md
from hoomd import deprecated

import gsd
from gsd import hoomd
from gsd import pygsd

import freud
from freud import parallel
from freud import box
from freud import density
from freud import cluster

import numpy as np
import math
import random
from scipy import stats

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
plt.rcParams.update({'font.size': 4})
plt.rcParams["xtick.major.size"] = 1
plt.rcParams["xtick.major.width"] = 0.5
plt.rcParams["xtick.minor.size"] = 1
plt.rcParams["xtick.minor.width"] = 0.5
plt.rcParams["ytick.major.size"] = 1
plt.rcParams["ytick.major.width"] = 0.5
plt.rcParams["ytick.minor.size"] = 1
plt.rcParams["ytick.minor.width"] = 0.5

def getDistance(point1, point2x, point2y):
    """Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance
    
def computeTauPerTstep(epsilon, mindt=0.000001):
    '''Read in epsilon, output tauBrownian per timestep'''
    kBT = 1.0
    tstepPerTau = float(epsilon / (kBT * mindt))
    return 1. / tstepPerTau

def roundUp(n, decimals=0):
    '''Round up size of bins to account for floating point inaccuracy'''
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier
    
def getNBins(length, minSz=(2**(1./6.))):
    "Given box size, return number of bins"
    initGuess = int(length) + 1
    nBins = initGuess
    # This loop only exits on function return
    while True:
        if length / nBins > minSz:
            return nBins
        else:
            nBins -= 1

# Functions vector quantities and angles
def quatToVector(quat, type, peZero, peOne):
    "Takes quaternion, returns orientation vector"
    if type == 0:
        mag = peZero
    else:
        mag = peOne
    x = quat[1] * mag
    y = quat[2] * mag
    act_vec = (x, y)
    return act_vec

def getMagnitude(vecF):
    "Take force vector, output magnitude"
    x = vecF[0]
    y = vecF[1]
    magF = np.sqrt((x**2)+(y**2))
    return magF

def quatToAngle(quat):
    "Take vector, output angle between [-pi, pi]"
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y, x)
    return rad

# Get infile and open
inFile = str(sys.argv[1])
f = hoomd.open(name=inFile, mode='rb')

# I NEED TO GRAB THE ACTIVITES AND POSITION FROM THE FILENAME
peIn = 150.
peOut = 500.

start = 0                   # first frame to process
dumps = int(f.__len__())    # get number of timesteps dumped
end = dumps                 # final frame to process
start = dumps - 2

box_data = np.zeros((1), dtype=np.ndarray)  # box dimension holder
r_cut = 2**(1./6.)                          # potential cutoff
tauPerDT = computeTauPerTstep(epsilon=1.)   # brownian time per timestep

# Set the colormap
myCols = plt.cm.viridis
fast = '#d8b365'
slow = '#5ab4ac'
colorsList = [fast, slow]
my_cmap = colors.ListedColormap(colorsList)
trackList = ['#D3D3D3']

myShrink = 0.65                             # shrink the colorbars
padCbar = 5

# Access and read .gsd data
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    first_tstep = snap.configuration.step
    box_data = snap.configuration.box
    # Box dimensions
    x_box = box_data[0]
    y_box = box_data[1]
    l_box = x_box
    hx_box = x_box / 2.0
    hy_box = y_box / 2.0
    h_box = hx_box
    a_box = x_box * y_box
    # Get the number particles
    pos = snap.particles.position
    partNum = len(pos)
    # Get 5 random particles to track
    track = []
    for j in xrange(5):
        track.append(random.randint(0, partNum - 1))
        col = plt.cm.jet(track[-1]/float(partNum))
        rgb = col[:3]
        trackList.append(colors.rgb2hex(rgb))
    track_cmap = colors.ListedColormap(trackList)

    # Compute mesh
    nBins = (getNBins(x_box, r_cut))
    sizeBin = roundUp((x_box / nBins), 6)
    # Check to see the size of the bin in use
    print("Length of box: {}").format(x_box)
    print("{} bins of size {}").format(nBins, sizeBin)

    for j in range(start, end):
        # Set the system snapshot
        snap = t[j]
        # Instantiate empty system mesh
        binParts = [[[] for b in range(nBins)] for a in range(nBins)]
        binForce = np.zeros((nBins, nBins, 2), dtype=np.float64)
        binFMag = np.zeros((nBins, nBins), dtype=np.float64)

        # Easier accessors
        pos = snap.particles.position               # position
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid                 # type
        ori = snap.particles.orientation            # orientation
        ang = np.array(list(map(quatToAngle, ori))) # convert to [-pi, pi]
        tst = snap.configuration.step               # timestep
        tst -= first_tstep                          # normalize by first timestep
        tst *= tauPerDT                             # convert to Brownian time

        # We want these to have default values for the gas
        effSigma = [1.] * partNum
        nearNeigh = [0] * partNum

        # Put particles in their respective bins
        for k in range(0, partNum):
            # Get mesh indices
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Append particle id to appropriate bin
            binParts[x_ind][y_ind].append(k)
            # Compute binned active force (INSIDE IS ALWAYS TYPE 0)
            vec = quatToVector(ori[k], typ[k], peIn, peOut)
            binForce[x_ind][y_ind][0] += vec[0]
            binForce[x_ind][y_ind][1] += vec[1]
        
        # Get the magnitude of the directional binned active force
        for k in range(0, nBins):
            for l in range(0, nBins):
                binFMag[k][l] = getMagnitude(binForce[k][l])
            

        # Compute distance, each pair will be counted twice
        for k in range(0, partNum):
            # Get mesh indices
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Get index of surrounding bins
            l_bin = x_ind - 1  # index of left bins
            r_bin = x_ind + 1  # index of right bins
            b_bin = y_ind - 1  # index of bottom bins
            t_bin = y_ind + 1  # index of top bins
            if r_bin == nBins:
                r_bin -= nBins  # adjust if wrapped
            if t_bin == nBins:
                t_bin -= nBins  # adjust if wrapped
            h_list = [l_bin, x_ind, r_bin]  # list of horizontal bin indices
            v_list = [b_bin, y_ind, t_bin]  # list of vertical bin indices

            # Loop through all bins
            for h in range(0, len(h_list)):
                for v in range(0, len(v_list)):
                    # Take care of periodic wrapping for position
                    wrapX = 0.0
                    wrapY = 0.0
                    if h == 0 and h_list[h] == -1:
                        wrapX -= l_box
                    if h == 2 and h_list[h] == 0:
                        wrapX += l_box
                    if v == 0 and v_list[v] == -1:
                        wrapY -= l_box
                    if v == 2 and v_list[v] == 0:
                        wrapY += l_box
                    # Compute distance between particles: effective radius
                    for b in range(0, len(binParts[h_list[h]][v_list[v]])):
                        ref = binParts[h_list[h]][v_list[v]][b]
                        r = getDistance(pos[k],
                                        pos[ref][0] + wrapX,
                                        pos[ref][1] + wrapY)
                        r = round(r, 4)  # round value to 4 decimal places
                        # Store the effective diameter for area calculation
                        if 0.1 < r < effSigma[k]:
                            effSigma[k] = r
                    # Compute distance between particles: nearest neighbors
                    for b in range(0, len(binParts[h_list[h]][v_list[v]])):
                        ref = binParts[h_list[h]][v_list[v]][b]
                        r = getDistance(pos[k],
                                        pos[ref][0] + wrapX,
                                        pos[ref][1] + wrapY)
                        r = round(r, 4) # round value to 4 decimal places
                        # If particles are near enough, increment neighbor count
                        if 0.1 < r < (r_cut * effSigma[k]):
                            nearNeigh[k] += 1


#        fig, ax = plt.subplots(2, 3)
#        # Plot the original (ax[0, 0])
#        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
#                                                        np.zeros_like(effSigma),
#                                                        offsets=xy, units='xy',
#                                                        cmap=my_cmap,
#                                                        transOffset=ax[0, 0].transData)
#        coll.set_array(np.ravel(typ))
#        ax[0, 0].add_collection(coll)
#        ax[0, 0].set_xlim(-h_box, h_box)
#        ax[0, 0].set_ylim(-h_box, h_box)
#        ax[0, 0].set_aspect('equal')
#        ax[0, 0].tick_params(axis='both', which='both',
#                             bottom=False, top=False, left=False, right=False,
#                             labelbottom=False, labeltop=False, labelleft=False,
#                             labelright=False)
#
#        # Plot the orientation (ax[0, 1])
#        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
#                                                        np.zeros_like(effSigma),
#                                                        offsets=xy, units='xy',
#                                                        cmap=plt.cm.hsv,
#                                                        transOffset=ax[0, 1].transData)
#        coll.set_array(np.ravel(ang))
#        cMin = min(ang)
#        cMax = max(ang)
#        coll.set_clim([cMin, cMax])
#        ax[0, 1].add_collection(coll)
#
##        divider = make_axes_locatable(ax[0, 1])
##        cax = divider.append_axes("right", "5%", pad="3%")
##        plt.colorbar(coll, cax=cax)
#
#        cbar = fig.colorbar(coll, ax=[ax[0, 1]], shrink=myShrink)
#        cbar.set_label(r'Orientation', labelpad=padCbar, rotation=270)
#        ax[0, 1].set_xlim(-h_box, h_box)
#        ax[0, 1].set_ylim(-h_box, h_box)
#        ax[0, 1].set_aspect('equal')
#        ax[0, 1].tick_params(axis='both', which='both',
#                             bottom=False, top=False, left=False, right=False,
#                             labelbottom=False, labeltop=False, labelleft=False,
#                             labelright=False)
#
#        # Plot the imshow magnitude of the binned active force (ax[0, 2])
#        coll = ax[0, 2].imshow(binFMag.T,
#                               extent=(0,nBins,0,nBins),
#                               origin='lower')
##        coll.set_clim([minCol, 6])
#        cbar = fig.colorbar(coll, ax=[ax[0, 2]], shrink=myShrink)
#        cbar.set_label(r'Binned Active Force', labelpad=padCbar, rotation=270)
#        ax[0, 2].set_xticks(())
#        ax[0, 2].set_yticks(())
#        ax[0, 2].set_aspect('equal')
#
#        # Plot the individual particle data (ax[1, 0])
#        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
#                                                        np.zeros_like(effSigma),
#                                                        offsets=xy, units='xy',
#                                                        cmap=track_cmap,
#                                                        transOffset=ax[1, 0].transData)
#        count = 1
#        trackCols = []
#        for k in xrange(partNum):
#            if k in track:
#                trackCols.append(count)
#                count += 1
#            else:
#                trackCols.append(0)
#        coll.set_array(np.ravel(trackCols))
#        ax[1, 0].add_collection(coll)
#        ax[1, 0].set_xlim(-h_box, h_box)
#        ax[1, 0].set_ylim(-h_box, h_box)
#        ax[1, 0].set_aspect('equal')
#        ax[1, 0].tick_params(axis='both', which='both',
#                             bottom=False, top=False, left=False, right=False,
#                             labelbottom=False, labeltop=False, labelleft=False,
#                             labelright=False)
#        count = 1
#        for k in track:
#            ax[1, 0].scatter(pos[k][0], pos[k][1], c=trackList[count])
#            count += 1
#
#        # Plot the radii heatmap (ax[1, 1])
#        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
#                                                        np.zeros_like(effSigma),
#                                                        offsets=xy, units='xy',
#                                                        cmap=myCols,
#                                                        transOffset=ax[1, 1].transData)
#        coll.set_array(np.ravel(effSigma))
#        minCol = min(effSigma)
#        coll.set_clim([minCol, 1.])
#        ax[1, 1].add_collection(coll)
#        cbar = fig.colorbar(coll, ax=[ax[1, 1]], shrink=myShrink, pad=1.)
#        cbar.set_label(r'Effective Diameter', labelpad=padCbar, rotation=270)
#        ax[1, 1].set_xlim(-h_box, h_box)
#        ax[1, 1].set_ylim(-h_box, h_box)
#        ax[1, 1].set_aspect('equal')
#        ax[1, 1].tick_params(axis='both', which='both',
#                             bottom=False, top=False, left=False, right=False,
#                             labelbottom=False, labeltop=False, labelleft=False,
#                             labelright=False)
#
#        # Plot the nearest neighbor heatmap (ax[1, 2])
#        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
#                                                        np.zeros_like(effSigma),
#                                                        offsets=xy, units='xy',
#                                                        cmap=myCols,
#                                                        transOffset=ax[1, 2].transData)
#        coll.set_array(np.ravel(nearNeigh))
#        minCol = min(nearNeigh)
#        coll.set_clim([minCol, 6])
#        ax[1, 2].add_collection(coll)
#        cbar = fig.colorbar(coll, ax=[ax[1, 2]], shrink=myShrink, pad=1.)
#        cbar.set_label(r'Nearest neighbors', labelpad=padCbar, rotation=270)
#        ax[1, 2].set_xlim(-h_box, h_box)
#        ax[1, 2].set_ylim(-h_box, h_box)
#        ax[1, 2].set_aspect('equal')
#        ax[1, 2].tick_params(axis='both', which='both',
#                             bottom=False, top=False, left=False, right=False,
#                             labelbottom=False, labeltop=False, labelleft=False,
#                             labelright=False)
#
#        # Save the file
#        pad = str(j).zfill(4)
#        fig.tight_layout()
##        plt.subplots_adjust(wspace=0.1, hspace=0.1)
#        plt.savefig('my_test_frame' + str(pad) + '.jpg', dpi=1000)



        # Try to plot with gridspec
        fig = plt.figure()
        gs = gridspec.Gridspec(2, 5)
        

        # Plot the original (ax[0, 0])
        ax0 = fig.add_subplot(gs[0, 0])
        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
                                                        np.zeros_like(effSigma),
                                                        offsets=xy, units='xy',
                                                        cmap=my_cmap,
                                                        transOffset=ax0.transData)
        coll.set_array(np.ravel(typ))
        ax0.add_collection(coll)
        ax0.set_xlim(-h_box, h_box)
        ax0.set_ylim(-h_box, h_box)
        ax0.set_aspect('equal')
        ax0.tick_params(axis='both', which='both',
                             bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labeltop=False, labelleft=False,
                             labelright=False)
        
        # Plot the orientation (ax[0, 1])
        ax1 = fig.add_subplot(gs[0, 1])
        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
                                                        np.zeros_like(effSigma),
                                                        offsets=xy, units='xy',
                                                        cmap=plt.cm.hsv,
                                                        transOffset=ax1.transData)
        coll.set_array(np.ravel(ang))
        cMin = min(ang)
        cMax = max(ang)
        coll.set_clim([cMin, cMax])
        ax1.add_collection(coll)
        # Set up the colorbar axis for orientation
        c1 = fig.add_subplot(gs[0, 2])
        cbar = plt.colorbar(coll, ax=c1, shrink=myShrink)
        cbar.set_label(r'Orientation', labelpad=padCbar, rotation=270)
        ax1.set_xlim(-h_box, h_box)
        ax1.set_ylim(-h_box, h_box)
        ax1.set_aspect('equal')
        ax1.tick_params(axis='both', which='both',
                             bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labeltop=False, labelleft=False,
                             labelright=False)
        
        # Plot the imshow magnitude of the binned active force (ax[0, 2])
        coll = ax[0, 2].imshow(binFMag.T,
                               extent=(0,nBins,0,nBins),
                               origin='lower')
#        coll.set_clim([minCol, 6])
        cbar = fig.colorbar(coll, ax=[ax[0, 2]], shrink=myShrink)
        cbar.set_label(r'Binned Active Force', labelpad=padCbar, rotation=270)
        ax[0, 2].set_xticks(())
        ax[0, 2].set_yticks(())
        ax[0, 2].set_aspect('equal')

        # Plot the individual particle data (ax[1, 0])
        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
                                                        np.zeros_like(effSigma),
                                                        offsets=xy, units='xy',
                                                        cmap=track_cmap,
                                                        transOffset=ax[1, 0].transData)
        count = 1
        trackCols = []
        for k in xrange(partNum):
            if k in track:
                trackCols.append(count)
                count += 1
            else:
                trackCols.append(0)
        coll.set_array(np.ravel(trackCols))
        ax[1, 0].add_collection(coll)
        ax[1, 0].set_xlim(-h_box, h_box)
        ax[1, 0].set_ylim(-h_box, h_box)
        ax[1, 0].set_aspect('equal')
        ax[1, 0].tick_params(axis='both', which='both',
                             bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labeltop=False, labelleft=False,
                             labelright=False)
        count = 1
        for k in track:
            ax[1, 0].scatter(pos[k][0], pos[k][1], c=trackList[count])
            count += 1
        
        # Plot the radii heatmap (ax[1, 1])
        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
                                                        np.zeros_like(effSigma),
                                                        offsets=xy, units='xy',
                                                        cmap=myCols,
                                                        transOffset=ax[1, 1].transData)
        coll.set_array(np.ravel(effSigma))
        minCol = min(effSigma)
        coll.set_clim([minCol, 1.])
        ax[1, 1].add_collection(coll)
        cbar = fig.colorbar(coll, ax=[ax[1, 1]], shrink=myShrink, pad=1.)
        cbar.set_label(r'Effective Diameter', labelpad=padCbar, rotation=270)
        ax[1, 1].set_xlim(-h_box, h_box)
        ax[1, 1].set_ylim(-h_box, h_box)
        ax[1, 1].set_aspect('equal')
        ax[1, 1].tick_params(axis='both', which='both',
                             bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labeltop=False, labelleft=False,
                             labelright=False)

        # Plot the nearest neighbor heatmap (ax[1, 2])
        coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
                                                        np.zeros_like(effSigma),
                                                        offsets=xy, units='xy',
                                                        cmap=myCols,
                                                        transOffset=ax[1, 2].transData)
        coll.set_array(np.ravel(nearNeigh))
        minCol = min(nearNeigh)
        coll.set_clim([minCol, 6])
        ax[1, 2].add_collection(coll)
        cbar = fig.colorbar(coll, ax=[ax[1, 2]], shrink=myShrink, pad=1.)
        cbar.set_label(r'Nearest neighbors', labelpad=padCbar, rotation=270)
        ax[1, 2].set_xlim(-h_box, h_box)
        ax[1, 2].set_ylim(-h_box, h_box)
        ax[1, 2].set_aspect('equal')
        ax[1, 2].tick_params(axis='both', which='both',
                             bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labeltop=False, labelleft=False,
                             labelright=False)
        
        # Save the file
        pad = str(j).zfill(4)
        fig.tight_layout()
#        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        plt.savefig('my_test_frame' + str(pad) + '.jpg', dpi=1000)
