#!/usr/bin/env python

'''
Creates initial configuration for a system containing branched linear chains.

'''

import sys
import math
import numpy as np
from configuration import Configuration
from config_io import *

#-------------------------------------------------------------------------------

config = Configuration()

#Add simulation box and boundary condition. The box size can be updated later if
#necessary.
config.add_simbox(20.0, 20.0, 20.0, 0)

#Molecule details
num_mols = 1
config.add_molecule_type('BTLBRS', num_mols)

#Number of backbone atoms
na_bbone = int(sys.argv[1])
#Number of atoms on each side chain
na_sc = int(sys.argv[2])
#Number of atoms between consecutive branch points
na_sp = 0
#Number of side chains growing from a branch point
f = 1
#Number of branch points
n_bp = 1 + (na_bbone-1)//(na_sp+1)
#Number of side chains
n_sc = n_bp*f
#Total number of side chain atoms
nat_sc = n_sc*na_sc
#Total number of atoms in this molecule
natm = na_bbone + nat_sc
#Total number of bonds in this molecule
nbnd = (na_bbone - 1) + nat_sc

#Atom types
#All atoms are point particles with unit mass
atm_t_bb = config.add_atom_type('BB', na_bbone, 0, 1.0)
atm_t_sc = config.add_atom_type('SC', nat_sc, 0, 1.0)

#Vdw interaction
eps = 1.0; sigma = 2.0; rcut = 2.0**(1.0/6)*sigma
config.add_ia_vdw('BB', 'BB', 'lj', np.array([eps, sigma, rcut]))
config.add_ia_vdw('BB', 'SC', 'lj', np.array([eps, sigma, rcut]))
config.add_ia_vdw('SC', 'SC', 'lj', np.array([eps, sigma, rcut]))

#config.add_ia_vdw('BB', 'SC', 'tab', np.array([50.0, 0.0, 1, 0]))
#config.add_ia_vdw('SC', 'SC', 'tab', np.array([0.0, 0.0, 0]))

#Bond types
bnd_t_bb = config.add_bond_type('kg', np.array([7.5, 3.0, eps, sigma]))
bnd_t_bs = config.add_bond_type('kg', np.array([7.5, 3.0, eps, sigma]))
bnd_t_ss = config.add_bond_type('kg', np.array([7.5, 3.0, eps, sigma]))

#Angle types
#ang_t_bb = config.add_angle_type('cos', np.array([2.0]))

len_bond = 1.94 #0.97*sigma #Equilibrium bond length
theta_bb = (45.0/180.0)*math.pi
sep = 1.9 # Must be <= len_bond

#Add atoms, bonds, etc
#Loop over all molecules (here only a single molecule type)
iatm = 0; ibnd = 0; iang = 0 #Index of atom, bond, etc.
for imol in range(1, num_mols+1):
    #Backbone
    #Put first atom at the origin
    iatm = config.add_atom(atm_t_bb, 0.0, np.zeros((3,)))
    iatm_beg = iatm
    #Put second atom along x-axis
    iatm = config.append_atom_bonded(atm_t_bb, 0.0, len_bond, 'alignx', iatm)
    #Put the rest of the backbone atoms using efrc links
    for i in range(3,na_bbone+1):
        iatm = config.append_atom_bonded(atm_t_bb, 0.0, len_bond, 'alignx',
                iatm, im2=iatm-1, theta=theta_bb, sep=sep)

    #Backbone bonds (=na_bbone-1)
    for i in range(1, na_bbone):
        ibnd = config.add_bond(bnd_t_bb, iatm_beg+i-1, iatm_beg+i)
        if i == 1:
            ibnd_beg = ibnd

    #Set backbone branch
    config.set_branch(-1, iatm_beg, iatm_beg+na_bbone-1)

    #Side chains (arranged as efjc)
    for i in range(1, na_bbone+1, na_sp+1):
        ibp = iatm_beg + i - 1 #Branch point
        for i_sc in range(1, f+1):
            for ia_sc in range(1, na_sc+1):
                if ia_sc == 1:
                    iatm = config.append_atom_bonded(atm_t_sc, 0.0, len_bond, 'efjc',
                            ibp, sep=sep)
                    ibnd = config.add_bond(bnd_t_bs, ibp, iatm)
                    ia_br_beg = iatm
                else:
                    iatm = config.append_atom_bonded(atm_t_sc, 0.0, len_bond, 'efjc',
                            iatm, sep=sep)
                    ibnd = config.add_bond(bnd_t_ss, iatm-1, iatm)
            #Set branch
            config.set_branch(ibp, ia_br_beg, na_sc)

    #Add molecule
    config.add_molecule(1, natm, iatm_beg, nbnd, ibnd_beg, 0, 0, 0, 0)
#config.update_simbox()

#Read in coordinates from external file
#natm = len(config.atoms)
#fn_coords = sys.argv[3]
#with open(fn_coords, 'r') as fh:
#    lines = fh.readlines()
#for i in range(natm):
#    words = lines[i+2].rstrip('\n').split()
#    config.atoms[i+1]['coords'][0] = float(words[0])
#    config.atoms[i+1]['coords'][1] = float(words[1])
#    config.atoms[i+1]['coords'][2] = float(words[2])
#End of reading coordinates

#write_ldf(config, 'tfr.txt', title='test_frame')
write_cfg(config, 'bb-%d-%d.cfg'%(na_bbone, na_sc), title='bb_cfg')
