''' 
Things this file should do
   1. reads in current:
        q_clust array
        particle_positions
        particle_types
   2. compute this tsteps:
        dx
        dy
        dz
   3. check to see if we've crossed a bound
   4. gives 
        msd_total tot/a/b
        msd_gas tot/a/b
        msd_liq tot/a/b
        "return msd_tot, msd_gas, msd_liq ..."
'''
import numpy as np

def msd( l_box, pos, prev_pos, disp_x=0, disp_y=0, disp_z=0):
    
    part_num = len(pos)
    msd_val = np.zeros((part_num), dtype=np.float64)
    for b in range(0, part_num):
        
        # Compute displacement
        dx = pos[b][0] - prev_pos[b][0]
        dy = pos[b][1] - prev_pos[b][1]
        dz = pos[b][2] - prev_pos[b][2]
    
        # Take care of periodicity
        if dx < -50:
            dx += l_box
        if dx > 50:
            dx -= l_box
        disp_x[b] += dx
        
        if dy < -50:
            dy += l_box
        if dy > 50:
            dy -= l_box
        disp_y[b] += dy
        
        if dz < -50:
            dz += l_box
        if dz > 50:
            dz -= l_box
        disp_z[b] += dz

        # Compute mean square displacement
        msd_val[b] = np.sqrt(((disp_x[b])**2) + ((disp_y[b])**2) + ((disp_z[b])**2))
            
    return msd_val
