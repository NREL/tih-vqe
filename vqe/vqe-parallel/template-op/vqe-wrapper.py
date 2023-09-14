#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os,sys
import numpy as np
from random import random
from scipy.optimize import minimize
import itertools
from mpi4py import MPI
import argparse
import pickle
from datetime import datetime

from qiskit import ClassicalRegister, transpile
from qiskit.chemistry.components.variational_forms import UCCSD
from qiskit.chemistry.components.initial_states import HartreeFock
from qiskit import BasicAer
from qiskit.circuit import Parameter as qp

import vqe_inner_funcs
import helpers
import vqeops


def logprint(*text):
    #print to log file
    printtext = ''
    for t in text:
        printtext += str(t) + ' '
    print(printtext)
    sys.stdout.flush()

def pauli_operator_to_dict(pauli_operator):
    #make dictionary of coefficients
    d = pauli_operator.to_dict()
    paulis = d['paulis']
    paulis_dict = {}

    for x in paulis:
        label = x['label']
        coeff = x['coeff']['real']
        paulis_dict[label] = coeff
    return paulis_dict

def addmeasure_to_vqe_circuit(circ, bases, creg, simulator):
    #add measurement bases gates
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
    #remove measurement bases gates
    if totaladded > 0:
        del circ.data[-1*totaladded:]
    return circ

def setupansatz(ucctype,num_spin_orbitals,num_particles):
    #make UCCS or UCCSD ansatz
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
    #determines bitstrings to use when calculating energy
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

def readbases(op,ijob,rank):
    #read in set of measurement bases
    if rank == 0:
        logprint('Trying to read bases set')
    allbasissets = helpers.readjson('basis-sets-parity-greedy-bonds/allbasissets-greedy.json')
    basesfile = allbasissets[str(op)][str(ijob)]
    tempfile = helpers.readfile('basis-sets-parity-greedy-bonds/' + basesfile)
    bases = [x.split() for x in tempfile]
    bases = [x[0] for x in bases if x != []]
    minbases = [list(x) for x in bases]

    basesfile = basesfile.replace('.dat','')
    tempfile = helpers.readfile('basis-sets-parity-greedy-bonds/' + basesfile + '_original.dat')
    bases = [x.split() for x in tempfile]
    bases = [x[0] for x in bases if x != []]
    allbases = [list(x) for x in bases]

    if rank == 0:
        logprint('Successful reading of greedy bases!')
    Nbs = len(minbases[0])   #get length of all bases
    return minbases,allbases,Nbs

def domapping(op,minbases,allbases):
    #determine which bases can be used to cover other bases in energy summation
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

def lookforparams(initial_points,ucctype,ijob,rank):
    #try to find existing params to use as starting guess
    paramsfile = 'params-'+str(ucctype)+'-'+str(ijob)
    if os.path.isfile(paramsfile):
        params = helpers.readfile(paramsfile)
        params = [x.split()[0] for x in params if x != '\n']
        params = [float(x) for x in params]
        if len(params) == len(initial_points):
            if rank == 0:
                logprint('Found params from old job, using them')
                Nr = 8
                printparams(params,Nr)
            return params
        else:
            print('Could not read in params file, length did not match needed num_params')
            return initial_points
    else:
        print('Could not find {:} file, using defaults'.format(paramsfile))
        return initial_points

def printparams(initial_points, Nr = 8):
    #print formatted list of params to log file
    for i in range(0,len(initial_points),Nr):
        tempdata = initial_points[i:i+Nr]    #works even if last row has less than Nr elements
        datastring = ' '
        for j in tempdata:
            datastring += '{:9.6f} '.format(j)
        logprint(datastring)

def parallel_vqe(initial_points,allcircs,parameters,var_form,bases,pauli_dict,
                 pauli_order,shift,backend,eigdict,bitstringorder,mapping,rank,stopflagfn):
    #run vqe function in parallel then aggregate results
    if rank == 0:
        srt = datetime.now()
    stopflagfn = comm.bcast(stopflagfn,root=0)

    finalE = 0
    all_vqe_results = np.zeros(len(pauli_dict)-1)
    if rank == 0:
        pauli_order = list(pauli_dict.keys())
    else:
        pauli_order = ['' for x in range(len(pauli_dict))]
    if stopflagfn == 0:
        combo = comm.bcast([initial_points,pauli_order],root=0)  #for some reason multiple bcasts overwrite the previous ones, so combine calls
        initial_points = combo[0]
        pauli_order = combo[1]

        if bases != []:
            #this line is where multi node jobs hang before qiskit tweaked
            local_vqe_result = vqe_inner_funcs.vqe(initial_points, allcircs,parameters,var_form,bases,pauli_dict,pauli_order,shift,backend,eigdict,bitstringorder,mapping)
        else:
            local_vqe_result = np.zeros(len(pauli_dict)-1)   #all energies are zero, do not include identity yet
        comm.Reduce(local_vqe_result,all_vqe_results,op=MPI.SUM,root=0)
        if rank == 0:
            logprint('compiling...')
            all_vqe_results_val = sum(all_vqe_results)
            Istring = 'I'*var_form.num_qubits
            logprint('terms',all_vqe_results_val,pauli_dict[Istring],shift)
            finalE = all_vqe_results_val+pauli_dict[Istring]+shift
            logprint('params')
            Nr = 8   #number of parameters to print per column
            printparams(initial_points,Nr)
            logprint('finalE',finalE)

    elif stopflagfn == 1:
        return stopflagfn

    if rank == 0:
        end = datetime.now()
        diff = end-srt
        logprint('steptime: {:.4f} seconds'.format(diff.total_seconds()))
        logprint('====================================================')

    return finalE

def parse_args():
    #get arguments from command line
    parser = argparse.ArgumentParser(description="script to run custom vqe")
    parser.add_argument("op", type=int, help="which operator function to use")
    parser.add_argument("ijob",type=int,help='job index (bond length index), used to generate parameters')
    parser.add_argument("ucctype",type=str, choices=['s','sd'],help='specify UCCS or UCCSD')
    args = parser.parse_args()
    return args


#parallelization info
comm = MPI.COMM_WORLD
size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

args = parse_args()
op = args.op
ijob = args.ijob
ucctype = args.ucctype

if rank == 0:
    jobstart = datetime.now()
    logprint('Op',op)
    logprint('Starting vqe-wrapper')

dist = 1.5 #default
genflag = False

#read in all pre-prepared info and print to log file
if os.path.isfile('opdata/opdata-bonds.json'):
    opdata = helpers.readjson('opdata/opdata-bonds.json')
    if rank == 0:
        logprint('Read opdata')
else:
    logprint('Could not find opdata-bonds.json!')
try:
    num_particles = opdata[str(op)][str(ijob)]['num_particles']
    num_spin_orbitals = opdata[str(op)][str(ijob)]['num_spin_orbitals']
    shift = opdata[str(op)][str(ijob)]['shift']
    dist = opdata[str(op)][str(ijob)]['bond_length']
    reference_energy = opdata[str(op)][str(ijob)]['reference_energy']
    if rank == 0:
        logprint('Bond index',ijob)
        logprint('Bond length',dist)
        logprint('num_particles',num_particles)
        logprint('num_spin_orbitals',num_spin_orbitals)
        logprint('shift',shift)
        logprint('reference_energy',reference_energy)
except:
    logprint('problem with setting opdata!')
    genflag = True
dirpath = 'allH/'
try:
    with open(dirpath + 'op' + str(op) + 'H-d'+str(ijob)+'-' + str(rank) + '.pkl','rb') as input:
        H = pickle.load(input)
    if rank == 0:
        logprint('Read H.pkl')
except:
    logprint('problem with ' + str(rank) + ' H reading. File: ' + dirpath + 'op' + str(op) + 'H-d'+str(ijob)+'-' + str(rank) + '.pkl')
    genflag = True
if genflag:
    H, num_particles, num_spin_orbitals, shift = vqeops.get_qubit_op(dist,op)

#do calculation setup
pauli_dict = pauli_operator_to_dict(H)
pauli_order = list(pauli_dict.keys())

var_form = setupansatz(ucctype,num_spin_orbitals,num_particles)

eigdict,bitstringorder = geteigarrays(pauli_dict)

minbases,allbases,Nbs = readbases(op,ijob,rank)
#next 3 lines split up bases by processor
basesindex = list(range(len(minbases)))
basesindex = np.array(basesindex[::size])+rank
bases = [minbases[x] for x in basesindex if x < len(minbases)]
###
mapping = domapping(op,minbases,allbases)

simulator = 'statevector'
backend = BasicAer.get_backend('statevector_simulator')

pnames = list(range(var_form.num_parameters))
pnames = [str(x) for x in pnames]
parameters = [qp(x) for x in pnames]

#make circuits
if rank == 0:
    logprint('making circuits...')
circ = var_form.construct_circuit(parameters)
c = ClassicalRegister(var_form._num_qubits)
allcircs = []
for bb,b in enumerate(bases):
    circ, totaladded = addmeasure_to_vqe_circuit(circ, b[::-1], c, simulator)
    transpiled_circ = transpile(circ,backend)
    circ = removemeasure_to_vqe_circuit(circ, totaladded)
    allcircs.append(transpiled_circ)

#start parallel_vqe part
if rank == 0:
    logprint('starting vqe...')
if rank == 0:
    stopflag = 0
    initial_points = [0 for x in range(var_form.num_parameters)]
    initial_points = lookforparams(initial_points,ucctype,ijob,rank)

    eps = 1e-5
    logprint('eps: {:}\n'.format(eps))
    vqe_result = minimize(parallel_vqe,initial_points,args=(allcircs,parameters,var_form,bases,pauli_dict,pauli_order,shift,backend,eigdict,bitstringorder,mapping,rank,stopflag),
                      method="SLSQP",tol=1e-6,bounds=var_form._bounds,options={'disp':True,'maxiter':1000,'eps':eps,'iprint':2})

    logprint('The estimated ground state energy from VQE algorithm is: {}'.format(vqe_result.fun))
    logprint('finalparams')
    for x in vqe_result.x:
        logprint(x)
    stopflag = 1
    #final iter
    tmp = parallel_vqe(vqe_result.x,allcircs,parameters,var_form,bases,pauli_dict,pauli_order,shift,backend,eigdict,bitstringorder,mapping,rank,stopflag)   #I believe this allows everything to terminate
else:
    #stop condition for other ranks
    stopflag = 0
    initial_points = [0.0 for x in range(var_form.num_parameters)]
    while stopflag == 0:
        stopflag = parallel_vqe(initial_points,allcircs,parameters,var_form,bases,pauli_dict,pauli_order,shift,backend,eigdict,bitstringorder,mapping,rank,stopflag)

ha2ryd = 2.0
ryd2ev = 13.605698066
if rank == 0:
    logprint()
    logprint('The exact ground state energy is: {}'.format(reference_energy))
    logprint()
    jobend = datetime.now()
    jobtime = jobend-jobstart
    logprint('Job runtime (s)',jobtime.total_seconds())
    logprint()

