#!/usr/bin/env python

import numpy as np

#-------------------------------------------------------------------------------

def write_xyz(config, fn, title=''):

    na = len(config.atoms)

    with open(fn,'w') as fh:
        fh.write(str(na) + '\n')
        fh.write(title + '\n')

        for i in range(1, na+1):
            coords = config.atoms[i]['coords']
            fh.write('  '.join(['% 15.6f'%x for x in coords]) + '\n')

#-------------------------------------------------------------------------------

def write_ldf(config, fn, title=''):

    with open(fn,'w') as fh:
        fh.write('#' + title + '\n')

        fh.write('%d atoms\n'%(len(config.atoms)))
        fh.write('%d atom types\n'%(len(config.atom_types)))

        if len(config.bonds) > 0:
            fh.write('%d bonds\n'%(len(config.bonds)))
            fh.write('%d bond types\n'%(len(config.bond_types)))

        if len(config.angles) > 0:
            fh.write('%d angles\n'%(len(config.angles)))
            fh.write('%d angle types\n'%(len(config.angle_types)))

        if len(config.dihedrals) > 0:
            fh.write('%d dihedrals\n'%(len(config.dihedrals)))
            fh.write('%d dihedral types\n'%(len(config.dihedral_types)))

        fh.write('%f  %f  xlo xhi\n'%(config.simbox[0,0], config.simbox[1,0]))
        fh.write('%f  %f  ylo yhi\n'%(config.simbox[0,1], config.simbox[1,1]))
        fh.write('%f  %f  zlo zhi\n'%(config.simbox[0,2], config.simbox[1,2]))
        fh.write('0  0  0  xy xz yz\n') #Tilt factors

        fh.write('\n')
        fh.write('Atoms # molecular\n')
        fh.write('\n')

        for imol in range(1, len(config.molecules)+1):
            natm = config.molecules[imol]['natm']
            iatm_beg = config.molecules[imol]['iatm_beg']
            for i in range(1, natm+1):
                at = config.atoms[iatm_beg+i-1]['type']
                coords = config.atoms[iatm_beg+i-1]['coords']
                buf = '%d  %d  %d  '%(i, imol, at) \
                    + '  '.join( ['% .8f'%x for x in coords] )\
                    + '\n'
                fh.write(buf)

        if len(config.bonds) > 0:
            fh.write('\n')
            fh.write('Bonds\n')
            fh.write('\n')
            for i in range(1, len(config.bonds)+1):
                bt = config.bonds[i]['type']
                atm_i = config.bonds[i]['atm_i']
                atm_j = config.bonds[i]['atm_j']
                buf = '%d  %d  %d  %d\n'%(i, bt, atm_i, atm_j)
                fh.write(buf)

        if len(config.angles) > 0:
            fh.write('\n')
            fh.write('Angles\n')
            fh.write('\n')
            for i in range(1, len(config.angles)+1):
                ant   = config.angles[i]['type']
                atm_i = config.angles[i]['atm_i']
                atm_j = config.angles[i]['atm_j']
                atm_k = config.angles[i]['atm_k']
                buf = '%d  %d  %d  %d  %d\n'%(i, ant, atm_i, atm_j, atm_k)
                fh.write(buf)

        if len(config.dihedrals) > 0:
            fh.write('\n')
            fh.write('Dihedrals\n')
            fh.write('\n')
            for i in range(1, len(config.dihedrals)+1):
                dt    = config.dihedrals[i]['type']
                atm_i = config.dihedrals[i]['atm_i']
                atm_j = config.dihedrals[i]['atm_j']
                atm_k = config.dihedrals[i]['atm_k']
                atm_l = config.dihedrals[i]['atm_l']
                buf = '%d  %d  %d  %d  %d  %d\n'%(i, dt, atm_i, atm_j, atm_k, atm_l)
                fh.write(buf)

#-------------------------------------------------------------------------------

def write_cfg(config, fn, title='Title'):
    '''
    Writes configuration to a config(.cfg) file.
    '''

    with open(fn,'w') as fh:
        fh.write('#' + title + '\n')
        fh.write('\nversion 2.0\n')

        fh.write('\nSIMBOX\n')
        #Box cell vectors: a, b, & c
        lx = config.simbox[1,0] - config.simbox[0,0]
        ly = config.simbox[1,1] - config.simbox[0,2]
        lz = config.simbox[1,2] - config.simbox[0,2]

        a = np.array([lx, 0, 0])
        b = np.array([0, ly, 0])
        c = np.array([0, 0, lz])

        fh.write('% .8f  % .8f  % .8f\n'%(a[0], a[1], a[2]))
        fh.write('% .8f  % .8f  % .8f\n'%(b[0], b[1], b[2]))
        fh.write('% .8f  % .8f  % .8f\n'%(c[0], c[1], c[2]))

        fh.write('IMCON  %d\n'%config.imcon)

        if len(config.atoms) > 0:
            fh.write('\nATOM_TYPES  %d\n'%len(config.atom_types))
            for i in range(1, len(config.atom_types)+1):
                fh.write('%s  %f  %d  %d\n'%(config.atom_types[i]['name'],
                    config.atom_types[i]['mass'], config.atom_types[i]['style'],
                    config.atom_types[i]['natm']))

        if len(config.atoms) > 0:
            fh.write('\nATOMS  %d\n'%len(config.atoms))
            for i in range(1, len(config.atoms)+1):
                buf = '%d  '%config.atoms[i]['type']
                buf += '%f  '%config.atoms[i]['charge']
                buf += '  '.join(['% .14E'%x for x in config.atoms[i]['coords']])
                if 'orientation' in config.atoms[i]:
                    buf += '  '.join(['%f'%x for x in config.atoms[i]['orientation']])
                fh.write(buf + '\n')

        if len(config.bonds) > 0:
            fh.write('\nBOND_TYPES  %d\n'%len(config.bond_types))
            for i in range(1, len(config.bond_types)+1):
                buf = '%s  '%config.bond_types[i]['style']
                buf += '%d  '%config.bond_types[i]['params'].size
                buf += '  '.join(['%f'%x for x in config.bond_types[i]['params']])
                fh.write(buf + '\n')

        if len(config.bonds) > 0:
            fh.write('\nBONDS  %d\n'%len(config.bonds))
            for i in range(1, len(config.bonds)+1):
                fh.write('%d  %d  %d\n'%(config.bonds[i]['type'],
                    config.bonds[i]['atm_i'], config.bonds[i]['atm_j']))

        if len(config.angles) > 0:
            fh.write('\nANGLE_TYPES  %d\n'%len(config.angle_types))
            for i in range(1, len(config.angle_types)+1):
                buf = '%s  '%config.angle_types[i]['style']
                buf += '%d  '%config.angle_types[i]['params'].size
                buf += '  '.join(['%f'%x for x in config.angle_types[i]['params']])
                fh.write(buf + '\n')

        if len(config.angles) > 0:
            fh.write('\nANGLES  %d\n'%len(config.angles))
            for i in range(1, len(config.angles)+1):
                fh.write('%d   %d   %d   %d\n'%(config.angles[i]['type'],
                    config.angles[i]['atm_i'], config.angles[i]['atm_j'],
                    config.angles[i]['atm_k']))

        if len(config.dihedrals) > 0:
            fh.write('\nDIHEDRAL_TYPES  %d\n'%len(config.dihedral_types))
            for i in range(1, len(config.dihedral_types)+1):
                buf = '%s  '%config.dihedral_types[i]['style']
                buf += '%d  '%config.dihedral_types[i]['params'].size
                buf += '  '.join(['%f'%x for x in config.dihedral_types[i]['params']])
                fh.write(buf + '\n')

        if len(config.dihedrals) > 0:
            fh.write('\nDIHEDRALS  %d\n'%len(config.dihedrals))
            for i in range(1, len(config.dihedrals)+1):
                fh.write('%d  %d  %d  %d  %d\n'%(config.dihedrals[i]['type'],
                    config.dihedrals[i]['atm_i'], config.dihedrals[i]['atm_j'],
                    config.dihedrals[i]['atm_k'], config.dihedrals[i]['atm_l']))

        if len(config.branches) > 0:
            fh.write('\nBRANCHES  %d\n'%len(config.branches))
            for i in range(1, len(config.branches)+1):
                fh.write('%d  %d  %d\n'%(config.branches[i]['iatm_teth'],
                    config.branches[i]['na_br'], config.branches[i]['iatm_beg']))

        if len(config.molecules) > 0:
            fh.write('\nMOLECULE_TYPES  %d\n'%len(config.molecule_types))
            for i in range(1, len(config.molecule_types)+1):
                fh.write('%s  %d\n'%(config.molecule_types[i]['name'],
                    config.molecule_types[i]['nmol']))

        if len(config.molecules) > 0:
            fh.write('\nMOLECULES  %d\n'%len(config.molecules))
            for i in range(1, len(config.molecules)+1):
                buf =  '%d  '%config.molecules[i]['type']
                buf += '%d  '%config.molecules[i]['natm']
                buf += '%d  '%config.molecules[i]['iatm_beg']
                buf += '%d  '%config.molecules[i]['nbnd']
                buf += '%d  '%config.molecules[i]['ibnd_beg']
                buf += '%d  '%config.molecules[i]['nang']
                buf += '%d  '%config.molecules[i]['iang_beg']
                buf += '%d  '%config.molecules[i]['ndhd']
                buf += '%d'%config.molecules[i]['idhd_beg']
                fh.write(buf+'\n')

        if len(config.tethers) > 0:
            fh.write('\nTETHERS  %d\n'%len(config.tethers))
            for i in range(1, len(config.tethers)+1):
                buf = '%s  '%config.tethers[i]['style']
                buf += '%d  '%config.tethers[i]['params'].size
                buf += '  '.join(['%f'%x for x in config.tethers[i]['params']])
                buf += '  %d  '%config.tethers[i]['iatm']
                buf += '  '.join(['%f'%x for x in config.tethers[i]['tp']])
                fh.write(buf + '\n')

        if len(config.vdw) > 0:
            fh.write('\nVDW  %d\n'%len(config.vdw))
            for i in range(1, len(config.vdw)+1):
                buf = '%d  '%config.vdw[i]['type_i']
                buf += '%d  '%config.vdw[i]['type_j']
                buf += '%s  '%config.vdw[i]['style']
                if config.vdw[i]['style'] != 'none':
                    buf += '%d  '%config.vdw[i]['params'].size
                    buf += '  '.join(['%f'%x for x in config.vdw[i]['params']])
                fh.write(buf + '\n')

        if len(config.external) > 0:
            fh.write('\nEXTERNAL  %d\n'%len(config.external))
            for i in range(1, len(config.external)+1):
                buf = '%s  '%config.external[i]['style']
                buf += '%d  '%config.external[i]['params'].size
                buf += '  '.join(['%f'%x for x in config.external[i]['params']])
                fh.write(buf + '\n')

#-------------------------------------------------------------------------------
