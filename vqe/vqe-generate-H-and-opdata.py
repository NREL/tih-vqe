#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os,sys
import numpy as np
import pickle
from qiskit.aqua.algorithms import NumPyEigensolver

import helpers
import vqeops

'''

generates H and exact energy data and stores in pkl for parallel use

'''

oplist = [2]
useonedist = False
store = True

dists = np.arange(0.5,3.6,0.1)
dists = [round(x,2) for x in dists]

dirname = 'opdata/'
filename = 'opdata-bonds.json'
flag = True

if os.path.isfile(dirname+filename) and flag:
    data = helpers.readjson(dirname+filename)
else:
    data = {str(x):{str(y):{'bond_length':0,'num_particles':0,'num_spin_orbitals':0,'shift':0,'complete':False,'reference_energy':0} for y in range(len(dists))} for x in oplist}

for op in oplist:
    if not str(op) in data:
        data[str(op)] = {str(y):{'bond_length':0,'num_particles':0,'num_spin_orbitals':0,'shift':0,'complete':False,'reference_energy':0} for y in range(len(dists))}
    for d in range(len(dists)):
        if not str(d) in data[str(op)]:
            data[str(op)][str(d)] = {'bond_length':0,'num_particles':0,'num_spin_orbitals':0,'shift':0,'complete':False,'reference_energy':0}

for op in oplist:
    print('op',op)
    for dd,dist in enumerate(dists):
        if useonedist:
            if dist != 1.8:
                continue
        print(dd)
        H, num_particles, num_spin_orbitals, shift = vqeops.get_qubit_op(dist, op)

        if store:
            with open('allH/op' + str(op) + 'H-d'+str(dd)+'.pkl','wb') as output:
                pickle.dump(H,output,pickle.HIGHEST_PROTOCOL)
                print('done pickle-ing')

        data[str(op)][str(dd)]['bond_length'] = dist
        data[str(op)][str(dd)]['num_particles'] = num_particles
        data[str(op)][str(dd)]['num_spin_orbitals'] = num_spin_orbitals
        data[str(op)][str(dd)]['shift'] = shift
        data[str(op)][str(dd)]['complete'] = True

        if op <= 14:
            exact_result = NumPyEigensolver(H).run()
            reference_energy = min(np.real(exact_result.eigenvalues))+shift
        else:
            reference_energy = 0
        data[str(op)][str(dd)]['reference_energy'] = reference_energy

        helpers.writejson(dirname+filename,data)


