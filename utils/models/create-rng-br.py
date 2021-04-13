#!/usr/bin/env python
'''
Creates initial configuration for a system containing branched rings.

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
config.add_simbox(200.0, 200.0, 200.0, 0)

#Molecule type
num_mols = 1
config.add_molecule_type('RNGBRSH', num_mols)

#Number of backbone atoms
na_bbone = 8 #int(sys.argv[1])
#Number of atoms on each side chain
na_sc = 2 #int(sys.argv[2])
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
nbnd = na_bbone + nat_sc

#Atom types
#All atoms are point particles with unit mass
atm_t_bb = config.add_atom_type('BB', na_bbone, 0, 1.0)
atm_t_sc = config.add_atom_type('SC', nat_sc, 0, 1.0)

#Vdw interaction
eps = 1.0; sigma = 2.0; rcut = 2.0**(1.0/6)*sigma
config.add_ia_vdw('BB', 'BB', 'lj', np.array([eps, sigma, rcut]))
config.add_ia_vdw('BB', 'SC', 'lj', np.array([eps, sigma, rcut]))
config.add_ia_vdw('SC', 'SC', 'lj', np.array([eps, sigma, rcut]))

#Bond types
bnd_t_bb = config.add_bond_type('kg', np.array([7.5, 3.0, eps, sigma]))
bnd_t_bs = config.add_bond_type('kg', np.array([7.5, 3.0, eps, sigma]))
bnd_t_ss = config.add_bond_type('kg', np.array([7.5, 3.0, eps, sigma]))

#Angle types (No angles)
#ang_t_bb = config.add_angle_type('cos', np.array([2.0]))

len_bond = 1.94 #0.97*sigma #Equilibrium bond length
theta_bb = (45.0/180.0)*math.pi
sep = 1.9 # Must be <= len_bond

#Loop over all molecules (here only a single molecule type)
iatm = 0; ibnd = 0; iang = 0 #Index of atom, bond, etc.
#Build the ring out of a polygon in the x-y plane.
theta = 2*math.pi/na_bbone
half_theta = theta/2
r = 0.5*len_bond/math.sin(half_theta)

for imol in range(1, num_mols+1):
    #Backbone
    for i in range(na_bbone):
        cos_theta = math.cos(i*theta+half_theta)
        sin_theta = math.sin(i*theta+half_theta)
        ri = np.array([r*cos_theta, r*sin_theta, 0.0])
        iatm = config.add_atom(atm_t_bb, 0, ri)
        if i == 0:
            iatm_beg = iatm

    #Backbone bonds 
    for i in range(1, na_bbone+1):
        atm_i = i; atm_j = 1 + i%na_bbone
        ibnd = config.add_bond(bnd_t_bb, iatm_beg+atm_i-1, iatm_beg+atm_j-1)
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

config.update_simbox()

write_ldf(config, 'tfr.txt', title='test_frame')
write_cfg(config, 'tcfg.cfg', title='test_cfg')
