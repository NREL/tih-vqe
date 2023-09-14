#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,sys
import numpy as np

import helpers


def gencommute(m1,m2):
    #generate commute
    return np.all(np.matmul(m1,m2) == np.matmul(m2,m1))

def qcommute(p1,p2):
    #qubit wise commute
    p1 = list(p1)
    p2 = list(p2)
    combo = list(zip(p1,p2))
    return np.all([gencommute(basisgates[x[0]],basisgates[x[1]]) for x in combo])


'''

check if any bases were missed

'''

#formal check of commutivity
I = np.array([[1,0],[0,1]])
X = np.array([[0,1],[1,0]])
Y = np.array([[0,np.complex(0,-1)],[np.complex(0,1),0]])
Z = np.array([[1,0],[0,-1]])
basisgates = {'I':I,'X':X,'Y':Y,'Z':Z}

oplist = [2]
idists = [0]  #need idist for each op in oplist
for i in range(len(oplist)):
    op = oplist[i]
    idist = idists[i]

    allbases = helpers.readfile('basesop'+str(op)+'-'+str(idist)+'_original.dat')
    minbases = helpers.readfile('basesop'+str(op)+'-'+str(idist)+'.dat')
    allbases = [x.split()[0] for x in allbases if x != '\n']
    minbases = [x.split()[0] for x in minbases if x != '\n']

    allcovered = True
    for b in allbases: #loop over all original bases
        flag = False
        if list(set(b)) == ['I']:
            continue
        for ab in minbases: #loop over greedy reduced bases
            if qcommute(b,ab): #check if any cover b
                flag = True #if not, then flag stays false
                break
        if not flag:
            allcovered = False
            print('{:} not covered!'.format(b))

    if allcovered:
        print('Looks good!')
