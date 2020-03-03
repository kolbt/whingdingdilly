'''
#                           This is an 80 character line                       #

This file will run in analyze.sh and will get the relevant parameters from the
source .gsd file names.

'''

import os
import sys

# Functions to grab parameters from files
def checkFile(fname, string):
    for i in range(0, len(fname)):
        if fname[i] == string[0]:
#             print"{} matches {}".format(fname[i], string[0])
            for j in range(1, len(string)):
                if fname[i + j] == string[j]:
#                     print"{} matches {}".format(fname[i+j], string[j])
                    if j == (len(string) - 1):
#                         print"Final match!"
                        return True
                else:
                    break
    return False
    
def txtValue(fname, string):
    out = ""
    index = 0
    for i in range(0, len(fname)):
        if fname[i] == string[0]:
            for j in range(1, len(string)):
                if (i + j) > (len(fname) - 1):
                    break
                elif fname[i + j] == string[j]:
                    if j == (len(string) - 1):
                        # Last index of search string
                        index = i + j
                else:
                    break
                        
    # First index of value
    index += 1
    mybool = True
    while mybool:
        if fname[index].isdigit():
            out = out + fname[index]
            index += 1
        elif fname[index] == ".":
            if fname[index+1].isdigit():
                out = out + fname[index]
                index += 1
            else:
                mybool = False
        else:
            mybool = False
    return float(out)
    
def getdtau(fname, string="dtau"):
    out = ""
    index = 0
    for i in range(0, len(fname)):
        if fname[i] == string[0]:
            for j in range(1, len(string)):
                if (i + j) > (len(fname) - 1):
                    break
                elif fname[i + j] == string[j]:
                    if j == (len(string) - 1):
                        # Last index of search string
                        index = i + j
                else:
                    break
                        
    # First index of value
    index += 1
    mybool = True
    while mybool:
        # We are at . before the file extension
        if fname[index + 1] == 'g':
            return float(out)
        out = out + fname[index]
        index += 1

# Read in the file name
file = str(sys.argv[1])

if checkFile(file, "pe"):
    pe = txtValue(file, "pe")
else:
    pe = 0
print(pe)

if checkFile(file, "pa"):
    paList = txtValue(file, "pa")
else:
    pa = 0
print(pa)

if checkFile(file, "pb"):
    pbList = txtValue(file, "pb")
else:
    pb = 0
print(pb)
    
if checkFile(file, "xa"):
    xa = txtValue(file, "xa")
else:
    xa = 100.
print(xa)

if checkFile(file, "ep"):
    if checkFile(file, "eps"):
        ep = txtValue(file, "eps")
    else:
        ep = txtValue(file, "ep")
else:
    ep = 1.
print(ep)
    
if checkFile(file, "phi"):
    phi = txtValue(file, "phi")
else:
    phi = 60.
print(phi)

if checkFile(file, "cluster"):
    cluster = 1
else:
    cluster = 0
print(cluster)

if checkFile(file, "dtau"):
    dtau = getdtau(file)
else:
    dtau = 0.000001
print(dtau)
