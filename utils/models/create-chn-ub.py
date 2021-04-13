#!/usr/bin/env python

'''
Creates initial configuration for a system containing unbranched linear chains.

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
config.add_molecule_type('CHN-UB', num_mols)

#Number of atoms
natm = 10 #int(sys.argv[1])

#Number of bonds
nbnd = natm - 1

#Number of angles. If no angles set num_angles to zero.
nang = natm - 2 if natm > 2 else 0

#Atom types
# natm atoms as point particles with mass 1.0
atm_t_bb = config.add_atom_type('M', natm, 0, 1.0)

#Vdw interaction
#LJ
eps = 1.0; sigma = 2.0; rcut = 2.5*sigma
config.add_ia_vdw('M', 'M', 'lj', np.array([eps, sigma, rcut]))

#Bond types
bnd_t_bb = config.add_bond_type('kg', np.array([7.5, 3.0, eps, sigma]))

#Angle types
ang_t_bb = config.add_angle_type('cos', np.array([2.0]))

#Add atoms, bonds, and angles
len_bond = 0.97*sigma #Equilibrium bond length
theta = (45.0/180.0)*math.pi
sep = len_bond# Must be <= len_bond

#Loop over all molecules (here only a single molecule type)
iatm = 0; ibnd = 0; iang = 0 #Index of atom, bond, etc.
for imol in range(1, num_mols+1):
    #First atom of molecule imol
    iatm = config.add_atom( atm_t_bb, 0.0, np.array([len_bond, 0.0, 0.0]) )
    iatm_beg = iatm
    #Second atom of molecule imol
    iatm = config.append_atom_bonded(1, 0.0, len_bond, 'alignx', iatm)

    for i in range(3,natm+1):
        iatm = config.append_atom_bonded(atm_t_bb, 0.0, len_bond, 'efjc',
                iatm, im2=iatm-1, theta=theta, sep=sep)
    
    for i in range(1, nbnd+1):
        ibnd = config.add_bond(bnd_t_bb, iatm_beg+i-1, iatm_beg+i)
        if i == 1:
            ibnd_beg = ibnd
    
    for i in range(1, nang+1):
        iang = config.add_angle(ang_t_bb, iatm_beg+i-1, iatm_beg+i, iatm_beg+i+1)
        if i == 1:
            iang_beg = iang
    
    config.add_molecule(1, natm, iatm_beg, nbnd, ibnd_beg, nang, iang_beg, 0, 0)

#Tether atom 1 with a harmonic bond to the origin
config.add_tether('harm', np.array([50.0, len_bond]), 1, np.zeros((3,)) )

#External field: stretching along x-axis with force fex
fex = 1.0
config.add_external('xfrc', np.array([1, natm, fex]) )

config.update_simbox()

#write_ldf(config, 'tfr.txt', title='test_frame')
#write_cfg(config, 'chn-%d.cfg'%natm, title='CHN-LJ-%d'%natm)
write_cfg(config, 'tcfg.cfg', title='test_cfg')
