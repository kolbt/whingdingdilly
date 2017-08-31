'''
This file will compute the contact force given:
    1. Reference particle
        -position (x and y)
    2. Contact particles (variable length arg)
ALTERNATIVELY, TAKE ONE PARTICLE AT A TIME SIMPLY
ADD TO THE CURRENT INDEX IN OUTSIDE LOOP
        -position (x and y)
        -net force vector (x and y)
It will return one vector value:
    1.) The vector contact force from surrounding
    particles (x and y)
'''

import numpy as np

def mag(x, y):
    magnitude = np.sqrt((x**2)+(y**2))
    return magnitude

def mechPressure(refpos_x, refpos_y, x_j, y_j, fx_j, fy_j):

    # initialize fc
    fc = np.zeros((3), dtype=np.float32)

    # compute Rij
    R_ijx = refpos_x - x_j
    R_ijy = refpos_y - x_y

    # compute angle between Fj and Rij
    theta = np.arccos( (( fx_j*R_ijx ) + ( fy_j*R_ijy ))\
                      / (mag(fx_j,fy_j) *  mag(R_ijx,R_ijy)) )

    # if gate theta
    if np.absolute(theta) >= (np.pi / 2):
        return fc
    
    # compute angle between x and R
    phi = arctan(R_ijy/R_ijx)

    # get the magnitude of the force vector
    Fj = mag(fx_j, fy_j)
    # in direction of ref particle
    Fc_mag = Fj * np.cos(theta)

    # get x and y comps of fc
    fc[0] = Fc_mag * np.cos(phi)
    fc[1] = Fc_mag * np.sin(phi)
    return fc
