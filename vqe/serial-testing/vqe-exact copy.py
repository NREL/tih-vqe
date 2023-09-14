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
oplist = [11]

distances = [2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]

data = {}

#use qiskit and numpy to get exact result
exact_energies = np.zeros([len(distances),len(oplist)])
for idist,dist in enumerate(distances):
    for iop,op in enumerate(oplist):
        print(op,idist)

        qubitOp, num_particles, num_spin_orbitals, shift = vqeops.get_qubit_op(dist,op)
        result = NumPyEigensolver(qubitOp).run()
        energy = np.real(result.eigenvalues)[0] + shift
        exact_energies[idist,iop] = energy
        print(energy)

