#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from qiskit.aqua.algorithms import NumPyEigensolver

import vqeops
import helpers


# script to generate exact diag results for different ops and plot results
# show effect of removing different orbitals on exact results

#see vqeops.py for oplist references
#LiH low/high spin removing -3,-2: no kink
oplist = [2,22]
labels = ['S = 1','S = 3']
#LiH low/high spin removing -2,-1: get kink
oplist = [23,24]
labels = ['S = 1','S = 3']

#TiH low/med/high spin removing -4,-3,-2,-1: get kink
oplist = [19,20,21]
labels = ['S = 2','S = 4','S = 6']
#TiH low/med/high spin removing -6,-5,-4,-3: get kink, but smaller
oplist = [9,10,11]
labels = ['S = 2','S = 4','S = 6']


#store results in json
regen = False
distances = [0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
             2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]

filename = 'exact.json'
if os.path.isfile(filename):
    data = helpers.readjson(filename)
else:
    data = {}

#use qiskit and numpy to get exact result
exact_energies = np.zeros([len(distances),len(oplist)])
for idist,dist in enumerate(distances):
    for iop,op in enumerate(oplist):
        print(op,idist)
        if regen:
            qubitOp, num_particles, num_spin_orbitals, shift = vqeops.get_qubit_op(dist,op)
            result = NumPyEigensolver(qubitOp).run()
            energy = np.real(result.eigenvalues)[0] + shift
            exact_energies[idist,iop] = energy
if not regen:
    exact_energies = np.array([data[str(op)] for op in oplist]).T

if regen:
    for iop,op in enumerate(oplist):
        data[str(op)] = exact_energies[:,iop].tolist()
        helpers.writejson(filename,data)

#min across spins
minvals = np.min(exact_energies,axis=1)

#plot
colors = ['dimgray','dodgerblue','orange']
fig = plt.figure(figsize=(5,4),dpi=200)
ax = fig.add_subplot(111)
lines = []
for iop in range(len(oplist)):
    line1, = ax.plot(distances, exact_energies[:,iop]-np.min(exact_energies[-1,:]),color=colors[iop],linestyle='-',marker = 'o',zorder=3)
    lines.append(line1)
line1, = ax.plot(distances, minvals-np.min(exact_energies[-1,:]),color='k',linestyle='-',linewidth=1,marker = None,zorder=3)
ax.legend(tuple(lines),tuple(labels))
ax.axhline(y=0,color='k',linestyle='-',zorder=2)
ax.set_xlabel('Bond length ($\mathrm{\AA}$)',fontsize=14)
ax.set_ylabel('E - E$_\mathrm{dissc}$ (Ha)',fontsize=14)
[ax.spines[x].set_linewidth(2) for x in ax.spines]
ax.tick_params(direction='out',width=2,length=4)
ax.grid(color='gainsboro',linestyle='-',zorder=1,axis='both')
ax.set_ylim([-0.2,0.4])
plt.tight_layout()
plt.show()
