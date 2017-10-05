import sys
import os

runfor = 100
dump_freq = 10
# tau = sigma^2 / diffusion coefficient
tau = 1
# dt = 2E-5 * tau, or, x * sigma^2 / Diffusion coeff
my_dt = 0.000001 * tau
# run for 100 tau, 100 * sigma^2 / Diffusion coeff
sim_length = runfor * tau
# compute number of tsteps to achieve this
tsteps = sim_length / my_dt
# calculate number of tsteps which are dumped
dumps = tsteps/dump_freq

import numpy as np

pow = np.log10(1/my_dt)
one_length = int(18*(pow-2)+29)             # gives length of array w/ values below 1 tau
num = tau / my_dt                           # this is 1 tau in terms of tsteps
spacer = tau / (5*my_dt)                    # 1/5th of tau, the spacer
gr_one_len = (tsteps - num)/spacer          # gives length of remaining array (tau > 1)
ar_tot_len = int(gr_one_len + one_length)

#get tsteps for msd calculations, needs to be in tau
msd_dumps = np.zeros((ar_tot_len), dtype=np.float64)
jumper = 5
value_to_dump = 15
count = 10

for iii in range(0,len(msd_dumps)):
    if iii <= 10:
        msd_dumps[iii] = iii
    elif value_to_dump * my_dt >= 1:
        msd_dumps[iii] = num
        num += spacer
    elif count == 95:
        msd_dumps[iii] = value_to_dump
        jumper *= 10
        value_to_dump += jumper
        count = 10
    else:
        msd_dumps[iii] = value_to_dump
        value_to_dump += jumper
        count += 5

#print(msd_dumps*my_dt)
ten_size = 0
for jjj in range(0,len(msd_dumps)):
    if msd_dumps[jjj]*my_dt <= 10:
        ten_size += 1

msd_ten = np.zeros((ten_size), dtype=np.float64)
for hhh in range(0,len(msd_ten)):
    msd_ten[hhh] = msd_dumps[hhh]

#print(msd_ten*my_dt)
msd_dumps += 110000
last_ten = 10 * tau / my_dt
msd_ten += tsteps - last_ten + 110000

#print(msd_ten)
#print(tsteps)
