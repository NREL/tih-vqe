#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os,sys
import itertools
import numpy as np
from shutil import copy2
import argparse

import helpers


def get_dir_files(a_dir = './'):
    #return list of all files within a_dir
    a_dir = a_dir.strip('/') + '/'    #ensure ends in slash
    filenames = [f for f in os.listdir(a_dir) if os.path.isfile(os.path.join(a_dir, f))]
    return filenames

def set_cover(universe, subsets):
    """Find a family of subsets that covers the universal set"""
    elements = set(e for s in subsets for e in s[0])
    # Check the subsets cover the universe
    if elements != universe:
        return None
    covered = set()
    cover = []
    # Greedily add the subsets with the most uncovered points
    while covered != elements:
        subset = max(subsets, key=lambda s: len(s[0] - covered))
        cover.append(subset)
        covered |= subset[0]
        print(len(covered),len(elements))
    return cover

def covered(b1,b2,Nb):
    #returns True if b1 "covered" by b2
    if all([True if b1[i] == 'I' or b1[i] == b2[i] else False for i in range(Nb)]):
        return True
    else:
        return False

def parse_args():
    #get arguments from command line
    parser = argparse.ArgumentParser(description="set cover script")
    parser.add_argument("-f", type=str, default=None, help='file to read')
    args = parser.parse_args()
    return args

args = parse_args()
filename = args.f

'''

reduces original set of pauli strings down using commutation rules using greedy algorithm

'''

oplist = [12]
onlyupdatejson = False
for op in oplist:
    mydir = 'basis-sets-parity-bonds/'
    basisfile = '../'+mydir + 'allbasissets.json'
    if os.path.isfile('allbasissets-greedy.json'):
        allbasissetsgreedy = helpers.readjson('allbasissets-greedy.json')
    else:
        allbasissetsgreedy = {}
    if os.path.isfile(basisfile):
        allbasissets = helpers.readjson(basisfile)
    else:
        print('need {:}'.format(basisfile))
        sys.exit()
    if not str(op) in allbasissetsgreedy:
        allbasissetsgreedy[str(op)] = {}
        #allbasissetsgreedy[str(op)] = allbasissets[str(op)]

    basesfiles = get_dir_files('../'+mydir)
    basesfiles = [x for x in basesfiles if 'basesop'+str(op)+'-' in x and '_original.dat' in x]
    print('files to process:')
    print(basesfiles)
    for bf in basesfiles:
        if filename is not None:
            if bf != filename:
                continue
        copy2('../'+mydir+bf,bf)
        newname = bf.replace('_original','')
        for i in allbasissets[str(op)]:
            if allbasissets[str(op)][i] == bf:
                allbasissetsgreedy[str(op)][i] = newname
        if onlyupdatejson:  #if run this in multiple places and need to combine later
            helpers.writejson('allbasissets-greedy.json',allbasissetsgreedy)
            continue

        bases = helpers.readfile(bf)
        bases = [x.split()[0] for x in bases if x != '\n']

        testbases = [x for x in bases]
        print(len(testbases))
        testbases2 = [list(x) for x in testbases]
        #replace pauli with prime number so simply checking if 2 indices commute
        #since we know I commutes with anything but others dont qwc commute with each other without extra work
        testbases2 = np.array([[0 if y == 'I' else 2 if y == 'X' else 3 if y == 'Y' else 5 if y == 'Z' else 9 for y in x] for x in testbases2])
        Nb = len(testbases[0])
        #all possible combos of pauli strings without an I of same string length as set to cover
        filledforI2 = np.array(list(itertools.product([2,3,5],repeat=Nb)))

        universe = set(range(len(testbases)))
        subsets = []

        total = len(filledforI2)  #number of all strings
        print('total permutations of length {:}: {:}'.format(Nb,total))
        Nblock = 5000  #chunk to process
        if len(filledforI2) < Nblock:  #reduce chunk if last one
            Nblock = len(filledforI2)
        currspot = 0
        #aggregates how many pauli strings each possible pauli string (without any I) can cover
        while currspot < total:
            repeatbases = np.repeat(testbases2[:,:,np.newaxis], Nblock, axis=2)
            cancover = np.sum(repeatbases % filledforI2[currspot:currspot+Nblock,:].T, axis = 1)
            subsets.extend([np.where(cancover[:,i] == 0)[0].tolist() for i in range(Nblock)])
            currspot += Nblock
            if total - currspot < Nblock:
                Nblock = total - currspot
            print('currspot',currspot)

        # same result as while loop but can require a ton of memory
        #    repeatbases = np.repeat(testbases2[:,:,np.newaxis], len(filledforI2), axis=2)
        #    cancover = np.sum(repeatbases % filledforI2.T, axis = 1)
        #    subsets = [np.where(cancover[:,i] == 0)[0].tolist() for i in range(len(filledforI2))]

        subsets = [(set(x),xx) for xx,x in enumerate(subsets)]
        print('cover')
        cover = set_cover(universe, subsets)
        bases = [filledforI2[x[1]] for x in cover]
        bases = [''.join(['X' if x == 2 else 'Y' if x == 3 else 'Z' for x in y]) for y in bases]
        print()
        print(bases)
        print(len(cover))
        print()
        bases = set(bases)
        print(bases)
        baseswrite = [x+'\n' for x in bases]
        helpers.writefile(newname,baseswrite)
        helpers.writejson('allbasissets-greedy.json',allbasissetsgreedy)

