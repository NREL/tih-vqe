#!/usr/bin/env python3 -W ignore::DeprecationWarning
# -*- coding: utf-8 -*-


import os,sys
import numpy as np
from random import random
from scipy.optimize import minimize
import itertools

from qiskit import ClassicalRegister, transpile, assemble
from qiskit.quantum_info import Statevector
from qiskit.aqua.algorithms import NumPyEigensolver
from qiskit import BasicAer
from qiskit.chemistry.components.variational_forms import UCCSD
from qiskit.chemistry.components.initial_states import HartreeFock
from qiskit.circuit import Parameter as qp

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

def addmeasure_to_vqe_circuit(circ, bases, creg, simulator):
    totaladded = 0
    for bb,b in enumerate(bases):
        if b == 'Z':
            totaladded += 1
        elif b == 'X':
            circ.h(bb)
            totaladded += 2
        elif b == 'Y':
            circ.sdg(bb)
            circ.h(bb)
            totaladded += 3
        if simulator == 'statevector' and b != 'I':
            totaladded -= 1
    return circ, totaladded

def removemeasure_to_vqe_circuit(circ, totaladded):
    if totaladded > 0:
        del circ.data[-1*totaladded:]
    return circ

def getexact(op,runexact=True):
    reference_energy = 0
    if runexact:
        exact_result = NumPyEigensolver(H).run()
        reference_energy = min(np.real(exact_result.eigenvalues))+shift
    return reference_energy

def setupansatz(num_spin_orbitals,num_particles,ucctype):
    initial_state = HartreeFock(
        num_spin_orbitals,
        num_particles,
        qubit_mapping='parity')

    var_form = UCCSD(excitation_type=ucctype,
        num_orbitals=num_spin_orbitals,
        num_particles=num_particles,
        initial_state=initial_state,
        qubit_mapping='parity')
    return var_form

def geteigarrays(pauli_dict):
    eigI = np.array([1,1])
    eigother = np.array([1,-1])
    eigs = {'I':eigI,'X':eigother,'Y':eigother,'Z':eigother}
    def calceig(pauli):
        pauli = [x for x in pauli]
        if len(pauli) <= 1:
            print('Must be at least 2 long')
            sys.exit()
        else:
            E = np.kron(eigs[pauli[0]],eigs[pauli[1]])
            for i in range(2,len(pauli)):
                E = np.kron(E,eigs[pauli[i]])
            return E
    eigdict = {x:calceig(x[::-1]) for x in pauli_dict}   #flip order
    bitstringorder = [str(x)[::-1] for x in ["".join(seq) for seq in itertools.product("01", repeat=var_form.num_qubits)]]
    return eigdict,bitstringorder

def readbases(op,ijob):
    allbasissets = helpers.readjson('../basis-sets-parity-bonds/allbasissets.json')
    basesfile = allbasissets[str(op)][str(ijob)]
    #tempfile = helpers.readfile('../basis-sets-parity-greedy-bonds/' + basesfile)
    #bases = [x.split() for x in tempfile]
    #bases = [x[0] for x in bases if x != []]
    #minbases = [list(x) for x in bases]
    minbases = []  #do not need this, just getting 1 circuit without basis gates

    basesfile = basesfile.replace('.dat','')
    if not '_original' in basesfile:
        basesfile = basesfile + '_original.dat'
    else:
        basesfile += '.dat'
    tempfile = helpers.readfile('../basis-sets-parity-greedy-bonds/' + basesfile)
    bases = [x.split() for x in tempfile]
    bases = [x[0] for x in bases if x != []]
    allbases = [list(x) for x in bases]

    Nbs = len(allbases[0])   #get length of all bases
    return minbases,allbases,Nbs

def domapping(op,minbases,allbases):
    def gencommute(m1,m2):
        return np.all(np.matmul(m1,m2) == np.matmul(m2,m1))

    def qcommute(p1,p2):
        p1 = list(p1)
        p2 = list(p2)
        combo = list(zip(p1,p2))
        return np.all([gencommute(basisgates[x[0]],basisgates[x[1]]) for x in combo])

    I = np.array([[1,0],[0,1]])
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,np.complex(0,-1)],[np.complex(0,1),0]])
    Z = np.array([[1,0],[0,-1]])
    basisgates = {'I':I,'X':X,'Y':Y,'Z':Z}

    allbases = [''.join(x) for x in allbases]
    minbases = [''.join(x) for x in minbases]
    # if rank == 0:
    #     print('all',allbases)
    #     print('min',minbases)

    mapping = {}
    for b in allbases:
        if list(set(b)) == ['I']:
            mapping[b] = b
            continue
        foundcover = False
        for ab in minbases:
            if qcommute(b,ab):
                mapping[b] = ab
                foundcover = True
                break
        if not foundcover:
            print('{:} was not covered! Check pauli commutation!'.format(b))
            sys.exit()
    return mapping

def vqe(initial_points,allcircs,parameters,var_form,bases,pauli_dict,pauli_order,shift,backend,eigdict,bitstringorder,mapping):

    allcounts = {''.join(x):{} for x in bases}
    for bb,b in enumerate(bases):
        if list(set(b)) != ['I']:  #only true if b is in all 'I' basis
            continue

        circparams = allcircs[bb]._parameter_table.get_keys()
        subparameters = [(x,xx) for xx,x in enumerate(parameters) if x in circparams]
        Np = len(subparameters)
        assembledjob = assemble(allcircs[bb].bind_parameters({subparameters[i][0]:initial_points[subparameters[i][1]] for i in range(Np)}))

        Nsqg = 0
        Ndqg = 0
        for g in assembledjob.to_dict()['experiments'][0]['instructions']:
            qg = g['qubits']
            if len(qg) == 1:
                Nsqg += 1
            elif len(qg) == 2:
                Ndqg += 1
            else:
                print('Gate does not have length 1 or 2!')
                sys.exit()
        print('Nsqg',Nsqg)
        print('Ndqg',Ndqg)

    return


oplist = [2]
ucctypes = ['s','sd']
ijobs = [13]

dists = np.arange(0.5,3.6,0.1)
dists = [round(x,2) for x in dists]

for op in oplist:
    for iucc,ucctype in enumerate(ucctypes):
        print(op,ucctype)
        for ijob in ijobs:

            dist = dists[ijob]

            runexact = False
            H, num_particles, num_spin_orbitals, shift = vqeops.get_qubit_op(dist,op)

            pauli_dict = pauli_operator_to_dict(H)  #this step could be combined with the previous line easily here
            pauli_order = list(pauli_dict.keys())
            reference_energy = getexact(op,runexact)

            var_form = setupansatz(num_spin_orbitals,num_particles,ucctype)
            if iucc == 0:
                print('number qubits',var_form.num_qubits)
                print('number spin orbitals',num_spin_orbitals)

            eigdict,bitstringorder = geteigarrays(pauli_dict)

            minbases,allbases,Nbs = readbases(op,ijob)
            #skip since just want original circuit
            #mapping = domapping(op,minbases,allbases)
            mapping = {}
            bases = [x for x in allbases]   #change to all bases so have identity circuit

            simulator = 'statevector'
            backend = BasicAer.get_backend('statevector_simulator')

            initial_points = [0 for x in range(var_form.num_parameters)]
            pnames = list(range(var_form.num_parameters))
            pnames = [str(x) for x in pnames]
            parameters = [qp(x) for x in pnames]

            #print('making circuits...')
            circ = var_form.construct_circuit(parameters)
            c = ClassicalRegister(var_form._num_qubits)
            allcircs = []
            for bb,b in enumerate(bases):
                if list(set(b)) != ['I']:
                    continue
                circ, totaladded = addmeasure_to_vqe_circuit(circ, b[::-1], c, simulator)
                transpiled_circ = transpile(circ,backend)
                circ = removemeasure_to_vqe_circuit(circ, totaladded)
                allcircs.append(transpiled_circ)

            vqe(initial_points,allcircs,parameters,var_form,bases,pauli_dict,pauli_order,shift,backend,eigdict,bitstringorder,mapping)
            print()



