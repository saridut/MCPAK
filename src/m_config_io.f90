module m_config_io

use m_precision
use m_strings
use m_globals

implicit none

contains

!******************************************************************************

subroutine read_dump(fn)
    !! Reads from DUMP file

    character(len=*), intent(in) :: fn
    integer :: i, it, npar
    integer :: fu, ierr

    open(newunit=fu, file=fn, access='stream', form='unformatted', &
        action='read', status='old')

    read(fu) leql, nblks, nts

    read(fu) smbx_a, smbx_b, smbx_c, imcon

    read(fu) num_atom_types
    allocate(atom_specs (num_atom_types))
    allocate(atom_pop   (num_atom_types))
    do it = 1, num_atom_types
        read(fu) atom_specs(it)%name, atom_specs(it)%mass, &
            atom_specs(it)%style, atom_pop(it)
    end do

    read(fu) num_atoms
    allocate(atoms       (num_atoms))
    allocate(charge      (num_atoms))
    allocate(coordinates(3,num_atoms))
    read(fu) atoms, charge, coordinates

    read(fu) num_bond_types
    if (num_bond_types > 0) then
        allocate(bond_specs(num_bond_types))
        do it = 1, num_bond_types
            read(fu) bond_specs(it)%style, npar
            allocate( bond_specs(it)%params(npar) )
            read(fu) bond_specs(it)%params(1:npar)
        end do
    end if

    read(fu) num_bonds
    if (num_bonds > 0) then
        allocate(bonds(3,num_bonds))
        read(fu) bonds
    end if

    read(fu) num_angle_types
    if (num_angle_types > 0) then
        allocate(angle_specs(num_angle_types))
        do it = 1, num_angle_types
            read(fu) angle_specs(it)%style, npar
            allocate( angle_specs(it)%params(npar) )
            read(fu) angle_specs(it)%params(1:npar)
        end do
    end if

    read(fu) num_angles
    if (num_angles > 0) then
        allocate(angles(4,num_angles))
        read(fu) angles
    end if

    read(fu) num_dihedral_types
    if (num_dihedral_types > 0) then
        allocate(dihedral_specs(num_dihedral_types))
        do it = 1, num_dihedral_types
            read(fu) dihedral_specs(it)%style, npar
            allocate( dihedral_specs(it)%params(npar) )
            read(fu) dihedral_specs(it)%params(1:npar)
        end do
    end if

    read(fu) num_dihedrals
    if (num_dihedrals > 0) then
        allocate(dihedrals(5,num_dihedrals))
        read(fu) dihedrals
    end if

    read(fu) num_branches
    if (num_branches > 0) then
        allocate(branches(3,num_branches))
        read(fu) branches
    end if

    read(fu) num_molecule_types
    allocate(molecule_names(num_molecule_types))
    allocate(molecule_pop  (num_molecule_types))
    read(fu) molecule_names, molecule_pop
    read(fu) num_molecules
    allocate(molecules(9,num_molecules))
    read(fu) molecules

    read(fu) num_tethers
    if (num_tethers > 0) then
        allocate(tethers(num_tethers))
        do i = 1, num_tethers
            read(fu) tethers(i)%style, npar
            allocate( tethers(i)%params(npar) )
            read(fu) tethers(i)%params(1:npar), tethers(i)%atm, tethers(i)%point
        end do
    end if

    read(fu) num_vdw_types
    if (num_vdw_types > 0) then
        allocate(vdw_specs(num_vdw_types))
        allocate(vdw_pairs(2,num_vdw_types))
        do it = 1, num_vdw_types
            read(fu) vdw_pairs(:,it), vdw_specs(it)%style, npar
            allocate( vdw_specs(it)%params(npar) )
            read(fu) vdw_specs(it)%params(1:npar)
        end do
    end if

    read(fu) num_externals
    if (num_externals > 0) then
        allocate( extrn_fields(num_externals) )
        do it = 1, num_externals
            read(fu) extrn_fields(it)%style, npar
            allocate( extrn_fields(it)%params(npar) )
            read(fu) extrn_fields(it)%params(1:npar)
        end do
    end if

    close(fu)

    end subroutine

!******************************************************************************

subroutine write_dump(fn)
    !! Writes to DUMP file

    character(len=*), intent(in) :: fn
    integer :: i, it, npar
    integer :: fu

    open(newunit=fu, file=fn, access='stream', form='unformatted', &
        action='write', status='replace')

    write(fu) leql, nblks, nts

    write(fu) smbx_a, smbx_b, smbx_c, imcon

    write(fu) num_atom_types
    do it = 1, num_atom_types
        write(fu) atom_specs(it)%name, atom_specs(it)%mass, &
            atom_specs(it)%style, atom_pop(it)
    end do
    write(fu) num_atoms, atoms, charge, coordinates

    write(fu) num_bond_types
    if (num_bond_types > 0) then
        do it = 1, num_bond_types
            npar = size(bond_specs(it)%params)
            write(fu) bond_specs(it)%style, npar, bond_specs(it)%params(1:npar)
        end do
    end if
    write(fu) num_bonds
    if (num_bonds > 0) write(fu) bonds

    write(fu) num_angle_types
    if (num_angle_types > 0) then
        do it = 1, num_angle_types
            npar = size(angle_specs(it)%params)
            write(fu) angle_specs(it)%style, npar, angle_specs(it)%params(1:npar)
        end do
    end if
    write(fu) num_angles
    if (num_angles > 0) write(fu) angles

    write(fu) num_dihedral_types
    if (num_dihedral_types > 0) then
        do it = 1, num_dihedral_types
            npar = size(dihedral_specs(it)%params)
            write(fu) dihedral_specs(it)%style, npar, &
                dihedral_specs(it)%params(1:npar)
        end do
    end if
    write(fu) num_dihedrals
    if (num_dihedrals > 0) write(fu) dihedrals

    write(fu) num_branches
    if (num_branches > 0) write(fu) branches

    write(fu) num_molecule_types, molecule_names, molecule_pop
    write(fu) num_molecules, molecules

    write(fu) num_tethers
    if (num_tethers > 0) then
        do i = 1, num_tethers
            npar = size(tethers(i)%params)
            write(fu) tethers(i)%style, npar, tethers(i)%params(1:npar), &
                tethers(i)%atm, tethers(i)%point
        end do
    end if

    write(fu) num_vdw_types
    if (num_vdw_types > 0) then
        do it = 1, num_vdw_types
            npar = size(vdw_specs(it)%params)
            write(fu) vdw_pairs(:,it), vdw_specs(it)%style, npar, &
                vdw_specs(it)%params(1:npar)
        end do
    end if

    write(fu) num_externals
    if (num_externals > 0) then
        do it = 1, num_externals
            npar = size(extrn_fields(it)%params)
            write(fu) extrn_fields(it)%style, npar, &
                extrn_fields(it)%params(1:npar)
        end do
    end if

    close(fu)

    end subroutine

!******************************************************************************

subroutine read_config(fn)
    !! Read from CONFIG file

    character(len=*), intent(in) :: fn
    character(len=mxrdln) :: line
    character(len=:), allocatable :: word
    integer :: i, it, npar, ibr, jt, at_i, at_j
    integer :: fu, ios, ierr

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '#', ios)
        if (ios /= 0) exit

        line = adjustl(line)

        if (str_startswith(line, 'SIMBOX')) then
            !Read simulation box size into lattice vectors defining the box
            read(fu,*) smbx_a(1:3)
            read(fu,*) smbx_b(1:3)
            read(fu,*) smbx_c(1:3)
        end if

        if (str_startswith(line, 'IMCON')) then
            !Read box boundary condition
            call str_split(line, ' ', word)
            imcon = str_to_i(line)
        end if

        if (str_startswith(line, 'ATOM_TYPES')) then
            call str_split(line, ' ', word)
            num_atom_types = str_to_i(line)

            allocate(atom_specs (num_atom_types))
            allocate(atom_pop   (num_atom_types))

            do it = 1, num_atom_types
                read(fu,*) atom_specs(it)%name, atom_specs(it)%mass, &
                    atom_specs(it)%style, atom_pop(it)
            end do
        end if

        if (str_startswith(line, 'ATOMS')) then
            call str_split(line, ' ', word)
            num_atoms = str_to_i(line)

            allocate(atoms (num_atoms))
            allocate(charge (num_atoms))
            allocate(coordinates(3,num_atoms))

            do i = 1, num_atoms
                read(fu,*) atoms(i), charge(i), coordinates(:,i)
            end do
        end if

        if (str_startswith(line, 'BOND_TYPES')) then
            call str_split(line, ' ', word)
            num_bond_types = str_to_i(line)

            allocate(bond_specs(num_bond_types))

            do it = 1, num_bond_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                bond_specs(it)%style = word
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( bond_specs(it)%params(npar) )
                if (npar > 0) then
                    read(line, *) bond_specs(it)%params
                end if
            end do
        end if

        if (str_startswith(line, 'BONDS')) then
            call str_split(line, ' ', word)
            num_bonds = str_to_i(line)
            allocate(bonds(3,num_bonds))
            do i = 1, num_bonds
                read(fu,*) bonds(:,i)
            end do
        end if

        if (str_startswith(line, 'ANGLE_TYPES')) then
            call str_split(line, ' ', word)
            num_angle_types = str_to_i(line)

            allocate(angle_specs(num_angle_types))

            do it = 1, num_angle_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                angle_specs(it)%style = word
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( angle_specs(it)%params(npar) )
                if (npar > 0) then
                    read(line, *) angle_specs(it)%params
                end if
            end do
        end if

        if (str_startswith(line, 'ANGLES')) then
            call str_split(line, ' ', word)
            num_angles = str_to_i(line)
            allocate(angles(4,num_angles))
            do i = 1, num_angles
                read(fu,*) angles(:,i)
            end do
        end if

        if (str_startswith(line, 'DIHEDRAL_TYPES')) then
            call str_split(line, ' ', word)
            num_dihedral_types = str_to_i(line)

            allocate(dihedral_specs(num_dihedral_types))

            do it = 1, num_dihedral_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                dihedral_specs(it)%style = word
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( dihedral_specs(it)%params(npar) )
                if (npar > 0) then
                    read(line, *) dihedral_specs(it)%params
                end if
            end do
        end if

        if (str_startswith(line, 'DIHEDRALS')) then
            call str_split(line, ' ', word)
            num_dihedrals = str_to_i(line)
            allocate(dihedrals(5,num_dihedrals))
            do i = 1, num_dihedrals
                read(fu,*) dihedrals(:,i)
            end do
        end if

        if (str_startswith(line, 'BRANCHES')) then
            call str_split(line, ' ', word)
            num_branches = str_to_i(line)
            allocate(branches(3,num_branches))
            do ibr = 1, num_branches
                read(fu,*) branches(:,ibr)
            end do
        end if

        if (str_startswith(line, 'MOLECULE_TYPES')) then
            call str_split(line, ' ', word)
            num_molecule_types = str_to_i(line)

            allocate(molecule_names(num_molecule_types))
            allocate(molecule_pop(num_molecule_types))

            do it = 1, num_molecule_types
                read(fu, *) molecule_names(it), molecule_pop(it)
            end do
        end if

        if (str_startswith(line, 'MOLECULES')) then
            call str_split(line, ' ', word)
            num_molecules = str_to_i(line)
            allocate(molecules(9,num_molecules))
            do i = 1, num_molecules
                read(fu,*) molecules(:,i)
            end do
        end if

        if (str_startswith(line, 'TETHERS')) then
            call str_split(line, ' ', word)
            num_tethers = str_to_i(line)

            allocate(tethers(num_tethers))

            do i = 1, num_tethers
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                tethers(i)%style = word
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( tethers(i)%params(npar) )
                if (npar > 0) then
                    read(line, *) tethers(i)%params, tethers(i)%atm, &
                        tethers(i)%point
                else
                    read(line, *) tethers(i)%atm, tethers(i)%point
                end if
            end do
        end if

        if (str_startswith(line, 'VDW')) then
            call str_split(line, ' ', word)
            num_vdw_types = str_to_i(line)

            allocate(vdw_specs(num_vdw_types))
            allocate(vdw_pairs(2,num_vdw_types))

            do it = 1, num_vdw_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                at_i = str_to_i(word)
                call str_split(line, ' ', word)
                at_j = str_to_i(word)
                if (at_i < at_j) then
                    jt = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
                    vdw_pairs(1,jt) = at_j
                    vdw_pairs(2,jt) = at_i
                else
                    jt = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
                    vdw_pairs(1,jt) = at_i
                    vdw_pairs(2,jt) = at_j
                end if
                call str_split(line, ' ', word)
                vdw_specs(jt)%style = word
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( vdw_specs(jt)%params(npar) )
                if (npar > 0) then
                    read(line, *) vdw_specs(jt)%params
                end if
            end do
        end if

        if (str_startswith(line, 'EXTERNAL')) then
            call str_split(line, ' ', word)
            num_externals = str_to_i(line)

            allocate(extrn_fields(num_externals))

            do it = 1, num_externals
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                extrn_fields(it)%style = word
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( extrn_fields(it)%params(npar) )
                if (npar > 0) then
                    read(line, *) extrn_fields(it)%params
                end if
            end do
        end if

    end do

    close(fu)

    end subroutine

!******************************************************************************

subroutine read_config_1(fn)
    !! Read config file corresponding to version mcpak-0.1.

    character(len=*), intent(in) :: fn
    character(len=mxrdln) :: line
    character(len=:), allocatable :: word
    integer :: i, it, npar, ibr, jt, at_i, at_j, styl
    integer :: fu, ios, ierr
    character(len=8), dimension(0:3) :: bond_styles
    character(len=8), dimension(0:2) :: angle_styles
    character(len=8), dimension(0:2) :: vdw_styles

    bond_styles = ['none', 'harm', 'fene', 'kg  ']
    angle_styles = ['none', 'cos ', 'harm']
    vdw_styles = ['none', 'lj  ', 'hs  ']

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '#', ios)
        if (ios /= 0) exit

        line = adjustl(line)

        if (str_startswith(line, 'SIMBOX')) then
            !Read simulation box size into lattice vectors defining the box
            read(fu,*) smbx_a(1:3)
            read(fu,*) smbx_b(1:3)
            read(fu,*) smbx_c(1:3)
        end if

        if (str_startswith(line, 'IMCON')) then
            !Read box boundary condition
            call str_split(line, ' ', word)
            imcon = str_to_i(line)
        end if

        if (str_startswith(line, 'ATOM_TYPES')) then
            call str_split(line, ' ', word)
            num_atom_types = str_to_i(line)

            allocate(atom_specs (num_atom_types))
            allocate(atom_pop   (num_atom_types))

            do it = 1, num_atom_types
                read(fu,*) atom_specs(it)%name, atom_pop(it), &
                    atom_specs(it)%style, atom_specs(it)%mass
                     
            end do
        end if

        if (str_startswith(line, 'ATOMS')) then
            call str_split(line, ' ', word)
            num_atoms = str_to_i(line)

            allocate(atoms (num_atoms))
            allocate(charge (num_atoms))
            allocate(coordinates(3,num_atoms))

            do i = 1, num_atoms
                read(fu,*) atoms(i), charge(i), coordinates(:,i)
            end do
        end if

        if (str_startswith(line, 'BOND_TYPES')) then
            call str_split(line, ' ', word)
            num_bond_types = str_to_i(line)

            allocate(bond_specs(num_bond_types))

            do it = 1, num_bond_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                styl = str_to_i(word)
                bond_specs(it)%style = bond_styles(styl)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( bond_specs(it)%params(npar) )
                if (npar > 0) then
                    read(line, *) bond_specs(it)%params
                end if
            end do
        end if

        if (str_startswith(line, 'BONDS')) then
            call str_split(line, ' ', word)
            num_bonds = str_to_i(line)
            allocate(bonds(3,num_bonds))
            do i = 1, num_bonds
                read(fu,*) bonds(:,i)
            end do
        end if

        if (str_startswith(line, 'ANGLE_TYPES')) then
            call str_split(line, ' ', word)
            num_angle_types = str_to_i(line)

            allocate(angle_specs(num_angle_types))

            do it = 1, num_angle_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                styl = str_to_i(word)
                angle_specs(it)%style = angle_styles(styl)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( angle_specs(it)%params(npar) )
                if (npar > 0) then
                    read(line, *) angle_specs(it)%params
                end if
            end do
        end if

        if (str_startswith(line, 'ANGLES')) then
            call str_split(line, ' ', word)
            num_angles = str_to_i(line)
            allocate(angles(4,num_angles))
            do i = 1, num_angles
                read(fu,*) angles(:,i)
            end do
        end if

        if (str_startswith(line, 'BRANCHES')) then
            call str_split(line, ' ', word)
            num_branches = str_to_i(line)
            allocate(branches(3,num_branches))
            do ibr = 1, num_branches
                read(fu,*) branches(:,ibr)
            end do
        end if

        if (str_startswith(line, 'MOLECULE_TYPES')) then
            call str_split(line, ' ', word)
            num_molecule_types = str_to_i(line)

            allocate(molecule_names(num_molecule_types))
            allocate(molecule_pop(num_molecule_types))

            do it = 1, num_molecule_types
                read(fu, *) molecule_names(it), molecule_pop(it)
            end do
        end if

        if (str_startswith(line, 'MOLECULES')) then
            call str_split(line, ' ', word)
            num_molecules = str_to_i(line)
            allocate(molecules(9,num_molecules))
            do i = 1, num_molecules
                read(fu,*) molecules(:,i)
            end do
        end if

        if (str_startswith(line, 'VDW')) then
            call str_split(line, ' ', word)
            num_vdw_types = str_to_i(line)

            allocate(vdw_specs(num_vdw_types))
            allocate(vdw_pairs(2,num_vdw_types))

            do it = 1, num_vdw_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                at_i = str_to_i(word)
                call str_split(line, ' ', word)
                at_j = str_to_i(word)
                if (at_i < at_j) then
                    jt = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
                    vdw_pairs(1,jt) = at_j
                    vdw_pairs(2,jt) = at_i
                else
                    jt = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
                    vdw_pairs(1,jt) = at_i
                    vdw_pairs(2,jt) = at_j
                end if
                call str_split(line, ' ', word)
                styl = str_to_i(word)
                vdw_specs(it)%style = vdw_styles(styl)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( vdw_specs(jt)%params(npar) )
                if (npar > 0) then
                    read(line, *) vdw_specs(jt)%params
                end if
            end do
        end if

        if (str_startswith(line, 'EXTERNAL')) then
            call str_split(line, ' ', word)
            num_externals = str_to_i(line)

            allocate(extrn_fields(num_externals))

            do it = 1, num_externals
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                extrn_fields(it)%style = word
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                allocate( extrn_fields(it)%params(npar) )
                if (npar > 0) then
                    read(line, *) extrn_fields(it)%params
                end if
            end do
        end if

    end do

    close(fu)

    end subroutine

!*******************************************************************************

subroutine write_config(fn, title)
    !! Write to cfg file

    character(len=*), intent(in) :: fn
    character(len=*), intent(in) :: title
    integer :: fu
    integer :: i, it, npar

    open(newunit=fu, file=fn, action='write', status='replace')

    write(fu, '(a)') '#'//trim(title)
    write(fu, *)
    write(fu, '(a)') 'version 2.0'

    write(fu, *)
    write(fu, '(a)') 'SIMBOX'
    write(fu, '(3(g0.6,2x))') smbx_a
    write(fu, '(3(g0.6,2x))') smbx_b
    write(fu, '(3(g0.6,2x))') smbx_c
    write(fu, '(a,2x,i0)') 'IMCON', imcon

    write(fu, *)
    write(fu, '(a,2x,i0)') 'ATOM_TYPES', num_atom_types
    do it = 1, num_atom_types
        write(fu, '(a,2x,g0.6,2x,i0,2x,i0)') trim(atom_specs(it)%name), &
            atom_specs(it)%mass, atom_specs(it)%style, atom_pop(it)
    end do

    write(fu, *)
    write(fu, '(a,2x,i0)') 'ATOMS', num_atoms
    do i = 1, num_atoms
        write(fu, '(i0,2x,g0.6,2x,*(es22.15,2x))') atoms(i), charge(i), coordinates(:,i)
    end do

    if (num_bonds > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'BOND_TYPES', num_bond_types
        do it = 1, num_bond_types
            npar = size(bond_specs(it)%params)
            write(fu, '(a,2x,i0,2x,*(g0.6,2x))') trim(bond_specs(it)%style), &
                npar, bond_specs(it)%params(1:npar)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'BONDS', num_bonds
        do i = 1, num_bonds
            write(fu, '(*(i0,2x))') bonds(:,i)
        end do
    end if

    if (num_angles > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'ANGLE_TYPES', num_angle_types
        do it = 1, num_angle_types
            npar = size(angle_specs(it)%params)
            write(fu, '(a,2x,i0,2x,*(g0.6,2x))') trim(angle_specs(it)%style), &
                npar, angle_specs(it)%params(1:npar)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'ANGLES', num_angles
        do i = 1, num_angles
            write(fu, '(*(i0,2x))') angles(:,i)
        end do
    end if

    if (num_dihedrals > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'DIHEDRAL_TYPES', num_dihedral_types
        do it = 1, num_dihedral_types
            npar = size(dihedral_specs(it)%params)
            write(fu, '(a,2x,i0,2x,*(g0.6,2x))') trim(dihedral_specs(it)%style), &
                npar, dihedral_specs(it)%params(1:npar)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'DIHEDRALS', num_dihedrals
        do i = 1, num_dihedrals
            write(fu, '(*(i0,2x))') dihedrals(:,i)
        end do
    end if

    if (num_branches > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'BRANCHES', num_branches
        do i = 1, num_branches
            write(fu, '(*(i0,2x))') branches(:,i)
        end do
    end if

    write(fu, *)
    write(fu, '(a,2x,i0)') 'MOLECULE_TYPES', num_molecule_types
    do it = 1, num_molecule_types
        write(fu, '(a,2x,i0)') trim(molecule_names(it)), molecule_pop(it)
    end do

    write(fu, *)
    write(fu, '(a,2x,i0)') 'MOLECULES', num_molecules
    do i = 1, num_molecules
        write(fu, '(*(i0,2x))') molecules(:,i)
    end do

    if (num_tethers > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'TETHERS', num_tethers
        do i = 1, num_tethers
            npar = size(tethers(i)%params)
            write(fu, '(a,2x,i0,2x,*(g0.6,2x))', advance='no') &
                trim(tethers(i)%style), npar, tethers(i)%params(1:npar)
            write(fu, '(i0,2x,3(g0.6,2x))') tethers(i)%atm, tethers(i)%point
        end do
    end if

    if (num_vdw_types > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'VDW', num_vdw_types
        do it = 1, num_vdw_types
            npar = size(vdw_specs(it)%params)
            write(fu,'(2(i0,2x),a,2x,i0,2x,*(g0.6,2x))') vdw_pairs(:,it), &
                trim(vdw_specs(it)%style), npar, vdw_specs(it)%params(1:npar)
        end do
    end if

    if (num_externals > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'EXTERNAL', num_externals
        do it = 1, num_externals
            npar = size( extrn_fields(it)%params )
            write(fu,'(a,2x,i0,2x,*(g0.6,2x))') trim(extrn_fields(it)%style), &
                npar, extrn_fields(it)%params(1:npar)
        end do
    end if

    close(fu)

    end subroutine

!******************************************************************************

subroutine write_ldf(fn_ld, title)
    !! Write to a LAMMPS data file. 

    character(len=*),       intent(in) :: fn_ld
    character(len=*),       intent(in) :: title
    real(rp), dimension(3) :: tilt_factors
    integer :: fu_ld
    integer :: cntr_atm, iatm_beg, natm
    integer :: i, iatm, imol

    open(newunit=fu_ld, file=fn_ld, action='write')

    !Header
    write(fu_ld,'(a)') '#'//trim(adjustl(title))

    write(fu_ld,'(i0,2x,a)') num_atoms, 'atoms'
    write(fu_ld,'(i0,2x,a)') num_atom_types, 'atom types'

    if (num_bonds > 0) then
        write(fu_ld,'(i0,2x,a)') num_bonds, 'bonds'
        write(fu_ld,'(i0,2x,a)') num_bond_types, 'bond types'
    end if

    if (num_angles > 0) then
        write(fu_ld,'(i0,2x,a)') num_angles, 'angles'
        write(fu_ld,'(i0,2x,a)') num_angle_types, 'angle types'
    end if

    if (num_dihedrals > 0) then
        write(fu_ld,'(i0,2x,a)') num_dihedrals, 'dihedrals'
        write(fu_ld,'(i0,2x,a)') num_dihedral_types, 'dihedral types'
    end if

    !Simulation box
    write(fu_ld,'(a,2x,g0.6,2x,a)') '0.0', smbx_a(1), 'xlo xhi'
    write(fu_ld,'(a,2x,g0.6,2x,a)') '0.0', smbx_b(2), 'ylo yhi'
    write(fu_ld,'(a,2x,g0.6,2x,a)') '0.0', smbx_c(3), 'zlo zhi'

    !Simulation box tilt factors
    tilt_factors = 0.0_rp
    write(fu_ld,'(3(g0.6,2x),a)') tilt_factors, 'xy xz yz'

    !Body: Atoms
    write(fu_ld,*)
    write(fu_ld,'(a)') 'Atoms # angle'
    write(fu_ld,*)

    cntr_atm = 1
    do imol = 1, num_molecules
        natm = molecules(2,imol)
        iatm_beg = molecules(3,imol)
        do i = 1, natm
            iatm = iatm_beg + i -1
            write(fu_ld,'(i0,2x,i0,2x,i0,2x,3(g0.8,2x))') cntr_atm, imol, &
                atoms(iatm), coordinates(:,iatm)
            cntr_atm = cntr_atm + 1
        end do
    end do

    !Body: Bonds
    if (num_bonds > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Bonds'
        write(fu_ld,*)
        do i = 1, num_bonds
            write(fu_ld,'(i0,2x,3(i0,2x))') i, bonds(:,i)
        end do
    end if

    !Body: Angles
    if (num_angles > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Angles'
        write(fu_ld,*)
        do i = 1, num_angles
            write(fu_ld,'(i0,2x,4(i0,2x))') i, angles(:,i)
        end do
    end if

    !Body: Dihedrals
    if (num_dihedrals > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Dihedrals'
        write(fu_ld,*)
        do i = 1, num_dihedrals
            write(fu_ld,'(i0,2x,5(i0,2x))') i, dihedrals(:,i)
        end do
    end if

    close(fu_ld)

    end subroutine

!*******************************************************************************

subroutine write_xyz(fn_xyz, title)
    !! Write to an XYZ file.

    character(len=*),       intent(in) :: fn_xyz
    character(len=*),       intent(in) :: title
    integer :: i, imol, iatm_beg, iatm
    integer :: natm
    integer :: fu_xyz

    open(newunit=fu_xyz, file=fn_xyz, action='write')

    write(fu_xyz, '(i0)') num_atoms
    write(fu_xyz, '(a)') title

    do imol = 1, num_molecules
        natm = molecules(2,imol)
        iatm_beg = molecules(3,imol)
        do i = 1, natm
            iatm = iatm_beg + i -1
            write(fu_xyz,'(3(es14.7,2x))') coordinates(:,iatm)
        end do
    end do

    close(fu_xyz)

    end subroutine

!*******************************************************************************

subroutine read_tabulated(ierr)

    integer, intent (out) :: ierr
    character(len=12) :: fn_vdw = 'tab_vdw.txt'
    character(len=12) :: fn_bond = 'tab_bnd.txt'
    character(len=12) :: fn_angle = 'tab_ang.txt'
    character(len=12) :: fn_dihedral = 'tab_dhd.txt'
    character(len=8)  :: an_i, an_j, an_if, an_jf
    character(len=:), allocatable :: fn
    character(len=mxrdln) :: line
    character(len=:), allocatable :: word
    integer :: fu, num_tab_pot, at_i, at_j, at_if, at_jf, typ
    integer :: i, it, ios, n
    logical :: lfnd

    ierr = 0

    !Read tabulated vdw potentials
    do it = 1, num_vdw_types
        if (trim(vdw_specs(it)%style) /= 'tab') cycle
        at_i = vdw_pairs(1,it); at_j = vdw_pairs(1,it)
        an_i = atom_specs(at_i)%name; an_j = atom_specs(at_j)%name
        lfnd = .false.

        !Open fn_vdw for reading
        open(newunit=fu, file=fn_vdw, action='read', status='old')

        call readline(fu, line, '#', ios)
        line = adjustl(line)
        read(line,*) num_tab_pot  !Number of tabulated potentials
        do i = 1, num_tab_pot
            call readline(fu, line, '#', ios)
            call str_split(line, ' ', word)
            an_if = word
            call str_split(line, ' ', word)
            an_jf = word
            if ( ((trim(an_i)==trim(an_if)).and.(trim(an_j)==trim(an_jf))) .or. &
                 ((trim(an_i)==trim(an_jf)).and.(trim(an_j)==trim(an_if))) ) then
                fn = trim(adjustl(line)); lfnd = .true.
                exit
            end if
        end do
        close(fu)

        if (.not. lfnd) then
            write(*,*) 'Tabulated vdw potential not found'
            write(*,'(a,1x,a,2x,a)') 'Atom names:', an_i, an_j
            ierr = 1; return
        end if

        !Read tabulated vdw potential
        open(newunit=fu, file=fn, action='read', status='old')
        call readline(fu, line, '#', ios)
        read(line,*) n
        vdw_specs(it)%tab_size = n
        allocate(vdw_specs(it)%tab_t(n))
        allocate(vdw_specs(it)%tab_v(n))
        do i = 1, n
            read(fu, *) vdw_specs(it)%tab_t(i), vdw_specs(it)%tab_v(i)
        end do
        close(fu)
    end do

    !Read tabulated bond potentials
    do it = 1, num_bond_types
        if (trim(bond_specs(it)%style) /= 'tab') cycle
        lfnd = .false.

        !Open fn_bond for reading
        open(newunit=fu, file=fn_bond, action='read', status='old')

        call readline(fu, line, '#', ios)
        line = adjustl(line)
        read(line,*) num_tab_pot  !Number of tabulated potentials
        do i = 1, num_tab_pot
            call readline(fu, line, '#', ios)
            call str_split(line, ' ', word)
            typ = str_to_i(word)
            if ( typ == it ) then
                fn = trim(adjustl(line)); lfnd = .true.
                exit
            end if
        end do
        close(fu)

        if (.not. lfnd) then
            write(*,*) 'Tabulated bond potential not found'
            write(*,'(a,1x,i0)') 'Bond type:', it
            ierr = 1; return
        end if

        !Read tabulated bond potential
        open(newunit=fu, file=fn, action='read', status='old')
        call readline(fu, line, '#', ios)
        read(line,*) n
        bond_specs(it)%tab_size = n
        allocate(bond_specs(it)%tab_t(n))
        allocate(bond_specs(it)%tab_v(n))
        do i = 1, n
            read(fu, *) bond_specs(it)%tab_t(i), bond_specs(it)%tab_v(i)
        end do
        close(fu)
    end do

    !Read tabulated angle potentials
    do it = 1, num_angle_types
        if (trim(angle_specs(it)%style) /= 'tab') cycle
        lfnd = .false.

        !Open fn_angle for reading
        open(newunit=fu, file=fn_angle, action='read', status='old')

        call readline(fu, line, '#', ios)
        line = adjustl(line)
        read(line,*) num_tab_pot  !Number of tabulated potentials
        do i = 1, num_tab_pot
            call readline(fu, line, '#', ios)
            call str_split(line, ' ', word)
            typ = str_to_i(word)
            if ( typ == it ) then
                fn = trim(adjustl(line)); lfnd = .true.
                exit
            end if
        end do
        close(fu)

        if (.not. lfnd) then
            write(*,*) 'Tabulated angle potential not found'
            write(*,'(a,1x,i0)') 'Angle type:', it
            ierr = 1; return
        end if

        !Read tabulated angle potential
        open(newunit=fu, file=fn, action='read', status='old')
        call readline(fu, line, '#', ios)
        read(line,*) n
        angle_specs(it)%tab_size = n
        allocate(angle_specs(it)%tab_t(n))
        allocate(angle_specs(it)%tab_v(n))
        do i = 1, n
            read(fu, *) angle_specs(it)%tab_t(i), angle_specs(it)%tab_v(i)
        end do
        close(fu)
    end do


    !Read tabulated dihedral potentials
    if (num_dihedrals > 0) then
        write(*,*) 'Reading tabulated dihedral potential not implemented yet.'
    end if

    end subroutine

!*******************************************************************************

end module m_config_io
