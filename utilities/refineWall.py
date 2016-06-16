#!/usr/bin/python

import os
import sys


number = int(sys.argv[1])
gamma    =  float(sys.argv[2])

patches = sys.argv[3:]

#################




gamma_sum = 0
for k in range(1,number):
    gamma_sum += gamma**k

# scaled to one
dX = [1/gamma_sum]
for k in range(1,number):
    dX.append(dX[0]*gamma**k)
dX = dX[::-1]
edgeV_lst = []
for k in range(number-1):
    eV = 1 - dX[k] / sum(dX[k:])
    edgeV_lst.append(eV)
    
    
for patch in patches:
    for edgeV in edgeV_lst: 
        os.system("refineWallLayer -overwrite {patch} {edgeV}"\
                    .format(patch=patch, edgeV=edgeV))
                
                
                
                
                
                
