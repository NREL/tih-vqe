#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,sys
import numpy as np

import helpers
import vqeops


def pauli_operator_to_dict(pauli_operator):
    d = pauli_operator.to_dict()
    paulis = d['paulis']
    paulis_dict = {}

    for x in paulis:
        label = x['label']
        coeff = x['coeff']['real']
        paulis_dict[label] = coeff
    return paulis_dict

def get_dir_files(a_dir = './'):
    #return list of all files within a_dir
    a_dir = a_dir.strip('/') + '/'    #ensure ends in slash
    filenames = [f for f in os.listdir(a_dir) if os.path.isfile(os.path.join(a_dir, f))]
    return filenames

'''

generates pauli string list for each bond length of op in vqeops

'''

oplist = [2]
useonedist = False
idist = 13

dists = np.arange(0.5,3.6,0.1)
dists = [round(x,2) for x in dists]

mydir = 'basis-sets-parity-bonds/'
basisfile = mydir + 'allbasissets.json'
if os.path.isfile(basisfile):
    allbasissets = helpers.readjson(basisfile)
else:
    allbasissets = {}

#makes list of pauli strings needed for each op and bond length
#stores it if it's unique
#allbasissets.json stores name of file containing this data for each op and bond length
for op in oplist:
    if not str(op) in allbasissets:
        allbasissets[str(op)] = {}
        existingbases = []
        existingbasesname = []
        allfiles = get_dir_files('basis-sets-parity-bonds/')
        filenames = [x for x in allfiles if 'basesop'+str(op)+'-' in x and not '_original' in x]

        #update allbasissets if data files already present from previous run
        if len(filenames) > 0:
            presentindices = list(sorted([int(x.split('_')[0].split('-')[1]) for x in filenames]))
            ispot = 0
            for i in range(len(dists)):
                curri = presentindices[ispot]
                currfile = 'basesop'+str(op)+'-'+str(curri)+'_original.dat'
                testfile = 'basesop'+str(op)+'-'+str(i)+'_original.dat'
                if testfile in filenames:
                    allbasissets[str(op)][str(i)] = testfile
                    if i > 0:  #so lags behind 1 index
                        ispot += 1
                else:
                    allbasissets[str(op)][str(i)] = currfile
            helpers.writejson(basisfile,allbasissets)

            for i in allbasissets[str(op)]:
                if not allbasissets[str(op)][i] in existingbasesname:
                    temp = helpers.readfile(mydir + allbasissets[str(op)][i])
                    temp = [x.split()[0] for x in temp if x != '\n']
                    existingbases.append(temp)
                    existingbasesname.append('basesop'+str(op)+'-'+str(i)+'_original.dat')
    else:
        #to cover case where extra bond lengths are added later
        allfiles = get_dir_files('basis-sets-parity-bonds/')
        filenames = [x for x in allfiles if 'basesop'+str(op)+'-' in x and not '_original' in x]
        if len(filenames) > 0:
            presentindices = list(sorted([int(x.split('_')[0].split('-')[1]) for x in filenames]))
            ispot = 0
            for i in range(len(dists)):
                curri = presentindices[ispot]
                currfile = 'basesop'+str(op)+'-'+str(curri)+'_original.dat'
                testfile = 'basesop'+str(op)+'-'+str(i)+'_original.dat'
                if testfile in filenames:
                    allbasissets[str(op)][str(i)] = testfile
                    ispot += 1
                else:
                    allbasissets[str(op)][str(i)] = currfile
            helpers.writejson(basisfile,allbasissets)

        existingbases = []
        existingbasesname = []
        for i in allbasissets[str(op)]:
            if not allbasissets[str(op)][i] in existingbasesname:
                temp = helpers.readfile(mydir + allbasissets[str(op)][i])
                temp = [x.split()[0] for x in temp if x != '\n']
                existingbases.append(temp)
                existingbasesname.append('basesop'+str(op)+'-'+str(i)+'_original.dat')

    for i in range(len(dists)):
        if useonedist:
            if i != idist:
                continue
        dist = dists[i]
        print(dist)
        if str(i) in allbasissets[str(op)]:
            print('already done')
            continue

        print('generating H...')
        H, num_particles, num_spin_orbitals, shift = vqeops.get_qubit_op(dist,op)
        print('done!')

        pauli_dict = pauli_operator_to_dict(H)

        bases = pauli_dict.keys()
        print('length of all bases: {:}'.format(len(bases)))
        print('Nspinorbs',num_spin_orbitals)
        baseswrite = [x+'\n' for x in bases]

        covered = False
        for j in range(len(existingbases)):
            if bases == existingbases[j]:
                allbasissets[str(op)][str(i)] = existingbasesname[j]
                covered = True
                break
        if not covered:
            allbasissets[str(op)][str(i)] = 'basesop'+str(op)+'-'+str(i)+'_original.dat'
            helpers.writefile(mydir + '/basesop' + str(op) + '-'+str(i)+'_original.dat',baseswrite)
            existingbases.append(bases)
            existingbasesname.append('basesop' + str(op) + '-'+str(i)+'_original.dat')
        helpers.writejson(basisfile,allbasissets)
