#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os,sys
import numpy as np

from qiskit.aqua.operators import Z2Symmetries
from qiskit.chemistry.drivers import PySCFDriver, UnitsType
from qiskit.chemistry import FermionicOperator


def get_qubit_op(dist, op):
    if op == 1:
        catbasis = 'sto3g'
        driver = PySCFDriver(atom="H 0 0 0; H 0 0 "+str(dist), unit=UnitsType.ANGSTROM,
                         charge=0, spin=0, basis={'H':catbasis})
        molecule = driver.run()
        repulsion_energy = molecule.nuclear_repulsion_energy
        num_particles = molecule.num_alpha+molecule.num_beta
        num_spin_orbitals = molecule.num_orbitals*2
        ferOp = FermionicOperator(h1=molecule.one_body_integrals, h2=molecule.two_body_integrals)
        qubitOp = ferOp.mapping(map_type='parity', threshold=0.00000001)
        qubitOp = Z2Symmetries.two_qubit_reduction(qubitOp, num_particles)
        shift = repulsion_energy
    else:
        #freeze core orbitals (like 1s for Li and 1s2s2p for Na)
        if op == 2:
            cation, spin, catbasis, freeze_list, remove_list = ['Li',0,'sto3g',[0],[-3,-2]] #remove orbitals that don't contribute to bonding
        elif op == 3:
            cation, spin, catbasis, freeze_list, remove_list = ['Li',0,'sto3g',[0],[]] #keep all unoccupied orbitals
        elif op == 4:
            cation, spin, catbasis, freeze_list, remove_list = ['Na',0,'sto3g',[0,1,2,3,4],[]]
        elif op == 5:
            cation, spin, catbasis, freeze_list, remove_list = ['K',0,'sto3g',[0,1,2,3,4,5,6,7,8],[]]

        #remove a lot of TiH least important orbitals
        elif op == 6:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',1,'sto3g',[0,1,2,3,4,5,6,7,8],[-8,-7,-6,-5,-4,-3]]
        elif op == 7:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',3,'sto3g',[0,1,2,3,4,5,6,7,8],[-8,-7,-6,-5,-4,-3]]
        elif op == 8:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',5,'sto3g',[0,1,2,3,4,5,6,7,8],[-8,-7,-6,-5,-4,-3]]

        #remove smaller set of TiH least important orbitals
        elif op == 9:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',1,'sto3g',[0,1,2,3,4,5,6,7,8],[-6,-5,-4,-3]]
        elif op == 10:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',3,'sto3g',[0,1,2,3,4,5,6,7,8],[-6,-5,-4,-3]]
        elif op == 11:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',5,'sto3g',[0,1,2,3,4,5,6,7,8],[-6,-5,-4,-3]]

        #removed smallest set of orbitals that do not significantly contribute to bonding, using MO coeffs
        elif op == 12:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',1,'sto3g',[0,1,2,3,4,5,6,7,8],[-4,-3]]
        elif op == 13:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',3,'sto3g',[0,1,2,3,4,5,6,7,8],[-4,-3]]
        elif op == 14:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',5,'sto3g',[0,1,2,3,4,5,6,7,8],[-4,-3]]

        #larger basis set
        elif op == 15:
            cation, spin, catbasis, freeze_list, remove_list = ['Li',0,'6-31g',[0],[-3,-2]]
        elif op == 16:
            cation, spin, catbasis, freeze_list, remove_list = ['Li',0,'6-31g',[0],[-7,-6,-3,-2]]
          #just used for cost estimate
        elif op == 17:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',3,'6-31g',[0,1,2,3,4,5,6,7,8],[]]
          #just used for cost estimate
        elif op == 18:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',3,'6-31g',[0,1,2,3,4,5,6,7,8],[-13,-12,-11,-8,-7,-5,-4,-3,-2]]


        #next 6 added to test removing different unoccupied orbitals
        #TiH remove arbitrary 4 orbitals that do contribute to bonding
        elif op == 19:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',1,'sto3g',[0,1,2,3,4,5,6,7,8],[-4,-3,-2,-1]]
        elif op == 20:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',3,'sto3g',[0,1,2,3,4,5,6,7,8],[-4,-3,-2,-1]]
        elif op == 21:
            cation, spin, catbasis, freeze_list, remove_list = ['Ti',5,'sto3g',[0,1,2,3,4,5,6,7,8],[-4,-3,-2,-1]]

        #LiH equivalents to TiH bonding orbital ops above
        elif op == 22:
            cation, spin, catbasis, freeze_list, remove_list = ['Li',2,'sto3g',[0],[-3,-2]] #high spin case of op2
          #remove 2 orbitals that do contribute to bonding (well, -1 does)
        elif op == 23:
            cation, spin, catbasis, freeze_list, remove_list = ['Li',0,'sto3g',[0],[-2,-1]]
        elif op == 24:
            cation, spin, catbasis, freeze_list, remove_list = ['Li',2,'sto3g',[0],[-2,-1]]

        else:
            raise ValueError('Could not find op {:}, add to vqeops.py!!!'.format(op))

        qubitOp, num_particles, num_spin_orbitals, shift = qubit_op_frzrm(cation, catbasis, dist, spin, freeze_list, remove_list)
    return qubitOp, num_particles, num_spin_orbitals, shift

#base functions used by each material
def pyscfmolc(cation, catbasis, dist, spin):
    driver = PySCFDriver(atom=str(cation)+" 0 0 0; H 0 0 " + str(dist), unit=UnitsType.ANGSTROM,
                         charge=0, spin=spin, basis={cation:catbasis,'H':'sto3g'})
    molecule = driver.run()
    return molecule

def qubit_op_frzrm(cation, catbasis, dist, spin, freeze_list, remove_list):
    #qiskit op with frozen/removed orbitals
    molecule = pyscfmolc(cation, catbasis, dist, spin)

    repulsion_energy = molecule.nuclear_repulsion_energy
    num_particles = molecule.num_alpha + molecule.num_beta
    num_spin_orbitals = molecule.num_orbitals*2
    remove_list = [x % molecule.num_orbitals for x in remove_list]
    freeze_list = [x % molecule.num_orbitals for x in freeze_list]
    remove_list = [x - len(freeze_list) for x in remove_list]
    remove_list += [x + molecule.num_orbitals - len(freeze_list)  for x in remove_list]
    freeze_list += [x + molecule.num_orbitals for x in freeze_list]
    ferOp = FermionicOperator(h1=molecule.one_body_integrals, h2=molecule.two_body_integrals)
    ferOp, energy_shift = ferOp.fermion_mode_freezing(freeze_list)
    num_spin_orbitals -= len(freeze_list)
    num_particles -= len(freeze_list)
    ferOp = ferOp.fermion_mode_elimination(remove_list)
    num_spin_orbitals -= len(remove_list)
    qubitOp = ferOp.mapping(map_type='parity', threshold=0.00000001)
    qubitOp = Z2Symmetries.two_qubit_reduction(qubitOp, num_particles)
    shift = energy_shift+repulsion_energy
    return qubitOp, num_particles, num_spin_orbitals, shift

def get_qubit_op1(dist, catbasis = 'sto3g', spin = 0):
    #for H2
    driver = PySCFDriver(atom="H 0 0 0; H 0 0 "+str(dist), unit=UnitsType.ANGSTROM,
                         charge=0, spin=spin, basis={'H':catbasis})
    molecule = driver.run()
    repulsion_energy = molecule.nuclear_repulsion_energy
    num_particles = molecule.num_alpha+molecule.num_beta
    num_spin_orbitals = molecule.num_orbitals*2
    ferOp = FermionicOperator(h1=molecule.one_body_integrals, h2=molecule.two_body_integrals)
    qubitOp = ferOp.mapping(map_type='parity', threshold=0.00000001)
    qubitOp = Z2Symmetries.two_qubit_reduction(qubitOp, num_particles)
    shift = repulsion_energy
    return qubitOp, num_particles, num_spin_orbitals, shift


