#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os,sys
import numpy as np

from qiskit import assemble
from qiskit.quantum_info import Statevector


def vqe(initial_points,allcircs,parameters,var_form,bases,pauli_dict,pauli_order,shift,backend,eigdict,bitstringorder,mapping):
    #do vqe measurement
    allcounts = {''.join(x):{} for x in bases}
    for bb,b in enumerate(bases):
        if bb == 0 and list(set(b)) == ['I']:  #only true if b is in all 'I' basis
            continue

        circparams = allcircs[bb]._parameter_table.get_keys()
        subparameters = [(x,xx) for xx,x in enumerate(parameters) if x in circparams]
        Np = len(subparameters)
        assembledjob = assemble(allcircs[bb].bind_parameters({subparameters[i][0]:initial_points[subparameters[i][1]] for i in range(Np)}))
        job = backend.run(assembledjob)
        result = job.result()
        statevector = result.get_statevector(allcircs[bb])
        newsv = Statevector(statevector)
        counts = newsv.probabilities_dict()

        allcounts[''.join(b)] = counts

    #aggregate energies
    paulistringenergies = {}
    for p in pauli_dict:
        if list(set(p)) == ['I']:
            continue
        if mapping[p] in allcounts:
            temp = [allcounts[mapping[p]][x] if x in allcounts[mapping[p]] else 0 for x in bitstringorder]
        else:
            temp = [0]*len(bitstringorder)
        paulistringenergies[p] = np.dot(eigdict[p],temp)*pauli_dict[p]

    all_H_energies = np.array([paulistringenergies[x] for x in pauli_order if list(set(x)) != ['I']])
    return all_H_energies

