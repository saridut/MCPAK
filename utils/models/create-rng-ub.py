#!/usr/bin/env python

'''
Creates initial configuration for a system containing unbranched rings.

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
config.add_simbox(1.0, 1.0, 1.0, 0)

#Molecule details
num_mols = 1
config.add_molecule_type('RNG-UB', num_mols)

#Total number of atoms
natm = 8 #int(sys.argv[1])

#Total number of bonds
nbnd = natm

#Total number of angles. If no angles set num_angles to zero.
nang = natm

#Atom types
# natm atoms as point particles with mass 1.0
atm_t_bb = config.add_atom_type('M', natm, 0, 1.0)

#Vdw interaction
eps = 1.0; sigma = 2.0; rcut = 2.0**(1.0/6)*sigma
config.add_ia_vdw('M', 'M', 'lj', np.array([eps, sigma, rcut]))

#Bond types
bnd_t_bb = config.add_bond_type('kg', np.array([7.5, 3.0, eps, sigma]))

#Angle types
ang_t_bb = config.add_angle_type('cos', np.array([2.0]))

#Add atoms, bonds, and angles
len_bond = 2.0 #0.97*sigma #Equilibrium bond length

#Loop over all molecules (here only a single molecule type)
iatm = 0; ibnd = 0; iang = 0 #Index of atom, bond, etc.
#Build the ring out of a polygon in the x-y plane.
theta = 2*math.pi/natm
half_theta = theta/2
r = 0.5*len_bond/math.sin(half_theta)

for imol in range(1, num_mols+1):
    for i in range(natm):
        cos_theta = math.cos(i*theta+half_theta)
        sin_theta = math.sin(i*theta+half_theta)
        ri = np.array([r*cos_theta, r*sin_theta, 0.0])
        iatm = config.add_atom(atm_t_bb, 0, ri)
        if i == 0:
            iatm_beg = iatm

    for i in range(1,nbnd+1):
        atm_i = i; atm_j = 1 + i%natm
        ibnd = config.add_bond(bnd_t_bb, iatm_beg+atm_i-1, iatm_beg+atm_j-1)
        if i == 1:
            ibnd_beg = ibnd

    for i in range(1,nang+1):
        atm_i = i; atm_j = 1 + i%natm; atm_k = 1 + (i+1)%natm
        iang = config.add_angle(ang_t_bb, iatm_beg+atm_i-1, iatm_beg+atm_j-1, 
                iatm_beg+atm_k-1)
        if i == 1:
            iang_beg = iang

    #Add molecule
    config.add_molecule(1, natm, iatm_beg, nbnd, ibnd_beg, nang, iang_beg, 0, 0)

config.update_simbox()

write_ldf(config, 'tfr.txt', title='test_frame')
write_cfg(config, 'tcfg.cfg', title='test_cfg')
