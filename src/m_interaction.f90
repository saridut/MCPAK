module m_interaction
!! This module provides the potentials listed below. Parameters for each
!! potential is given by the vector `params` of maximum length `mxparam`.
!!

use m_precision
use m_constants_math
use m_qsort, only: dqsort
use m_globals
use m_nbr_lists
use m_ia_bond
use m_ia_angle
use m_ia_dihedral
use m_ia_vdw
use m_ia_tether
use m_ia_external

implicit none

real(rp), dimension(:), allocatable :: xbuf
integer,  dimension(:), allocatable :: indx_xsorted

contains

!******************************************************************************

subroutine ia_setup()
    !! Sets up parameters for potentials 

    integer :: it, i

    do it = 1, num_bond_types
        call ia_bond_setup(bond_specs(it))
    end do

    if (num_angles > 0) then
        do it = 1, num_angle_types
            call ia_angle_setup(angle_specs(it))
        end do
    end if

    if (num_dihedrals > 0) then
        do it = 1, num_dihedral_types
            call ia_dihedral_setup(dihedral_specs(it))
        end do
    end if

    if (lvdw) then
        do it = 1, num_vdw_types
            call ia_vdw_setup(vdw_specs(it))
        end do

        allocate(xbuf(num_atoms))
        xbuf = 0.0_rp
        allocate(indx_xsorted(num_atoms))
        indx_xsorted = 0
    end if

    if (num_tethers > 0) then
        do i = 1, num_tethers
            call ia_tether_setup(tethers(i))
        end do
    end if

    if (num_externals > 0) then
        do i = 1, num_externals
            call ia_external_setup(extrn_fields(i))
        end do
    end if

    end subroutine

!******************************************************************************

subroutine ia_calc_energy(ierr)
    !! Calculates total energy

    integer, intent(out) :: ierr

    ierr = 0

    !Calculation for pairwise interactions
    if (lvdw) then
        if (use_verlet_tab) then
            call ia_calc_vdw_energy_vl(ierr)
        else
            call ia_calc_vdw_energy(ierr)
        end if
        if (ierr /= 0) return
    end if

    !Calculation of bonded interactions
    call ia_calc_bond_energy(ierr)
    if (ierr /= 0) return

    !Calculation of angular interactions
    if (num_angles > 0) call ia_calc_angle_energy()

    !Calculation of dihedral interactions
    if (num_dihedrals > 0) call ia_calc_dihedral_energy()

    !Calculation of tether interactions
    if (num_tethers > 0) then
        call ia_calc_tether_energy(ierr)
        if (ierr /= 0) return
    end if

    !Calculation of external interactions
    if (num_externals > 0) call ia_calc_external_energy(ierr)

    end subroutine

!********************************************************************************

subroutine ia_calc_vdw_energy(ierr)
    !! Calculates total energy due to short-ranged non-bonded pairwise interactions.
    !! Uses direct N^2 calculation. Will overwrite `energy_vdw` in module `m_globals`.

    integer, intent(out) :: ierr
    real(rp), dimension(3) :: ri, rj
    real(rp) :: enrg
    integer :: i, j, at_i, at_j, is, js
    integer :: typ

    ierr = 0
    energy_vdw = 0.0_rp

    !Sorting
    xbuf = coordinates(1,:)
    call dqsort(xbuf, indx_xsorted)

    do is = 1, (num_atoms-1)
        i = indx_xsorted(is)
        ri = coordinates(:,i)
        at_i = atoms(i)
        do js = is+1, num_atoms
            j = indx_xsorted(js)
            !Check if atom j is an excluded atom for atom i
            if ( exat_tab%is_in(i,j) ) cycle
            rj = coordinates(:,j)
            at_j = atoms(j)

            if (at_i < at_j) then
                typ = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
            else
                typ = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
            end if

            call ia_get_vdw_energy(ri, rj, vdw_specs(typ), enrg, ierr)

            if (ierr /= 0) return
            energy_vdw = energy_vdw + enrg
        end do
    end do

    end subroutine

!********************************************************************************

subroutine ia_calc_vdw_energy_vl(ierr)
    !! Calculates total energy due to short-ranged pairwise interactions.
    !! Uses Verlet list. The Verlet list takes care of excluded atoms.

    integer, intent(out) :: ierr
    integer, dimension(:), pointer :: nbrs => null()
    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: rj
    real(rp) :: enrg
    integer :: i, j, k
    integer :: at_i, at_j, typ

    ierr = 0; energy_vdw = 0.0_rp

    call verlet_build()

    do i = 1, num_atoms
        ri = coordinates(:,i)
        at_i = atoms(i)

        !Getting list of neighbors of particle i using pointer nbrs
        call verlet_tab%get_row(i, nbrs)

        do k = 1, size(nbrs)
            j = nbrs(k)
            rj = coordinates(:,j)
            at_j = atoms(j)

            if (at_i < at_j) then
                typ = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
            else
                typ = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
            end if

            call ia_get_vdw_energy(ri, rj, vdw_specs(typ), enrg, ierr)

            if (ierr /= 0) return
            energy_vdw = energy_vdw + enrg
        end do
    end do

    energy_vdw = 0.5*energy_vdw !MC verlet table

    end subroutine
    
!********************************************************************************

subroutine ia_calc_atm_vdw_energy(iatm, enrg_atm, ierr)
    !! Calculates the pair energy of a single atom. No Verlet list is used.

    integer, intent(in) :: iatm
    real(rp), intent(out) :: enrg_atm
    integer, intent(out) :: ierr
    real(rp), dimension(3) :: ri, rj
    real(rp) :: enrg
    integer :: j, typ, at_i, at_j

    ierr = 0; enrg_atm = 0.0_rp

    ri = coordinates(:,iatm)
    at_i = atoms(iatm)

    !Pair energy
    do j = 1, num_atoms
        !Ignore self and atoms bonded to iatm
        if ( exat_tab%is_in(iatm,j) ) cycle
        rj = coordinates(:,j)
        at_j = atoms(j)

        if (at_i < at_j) then
            typ = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
        else
            typ = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
        end if

        call ia_get_vdw_energy(ri, rj, vdw_specs(typ), enrg, ierr)

        if (ierr /= 0) return
        enrg_atm = enrg_atm + enrg
    end do

    end subroutine

!********************************************************************************

subroutine ia_calc_atm_vdw_energy_vl(iatm, enrg_atm, ierr)
    !! Calculates the pair energy of a single atom using Verlet list. The Verlet
    !! list takes care of excluded atoms.

    integer, intent(in) :: iatm
    real(rp), intent(out) :: enrg_atm
    integer, intent(out) :: ierr
    integer, dimension(:), pointer :: nbrs => null()
    real(rp), dimension(3) :: ri, rj
    real(rp) :: enrg
    integer :: j, k, typ, at_i, at_j

    ierr = 0; enrg_atm = 0.0_rp

    ri = coordinates(:,iatm)
    at_i = atoms(iatm)

    call verlet_build()

    !Getting list of neighbors of particle i using pointer nbrs
    call verlet_tab%get_row(iatm, nbrs)

    do k = 1, size(nbrs)
        j = nbrs(k)
        rj = coordinates(:,j)
        at_j = atoms(j)

        if (at_i < at_j) then
            typ = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
        else
            typ = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
        end if

        call ia_get_vdw_energy(ri, rj, vdw_specs(typ), enrg, ierr)

        if (ierr /= 0) return
        enrg_atm = enrg_atm + enrg
    end do

    end subroutine

!*******************************************************************************

subroutine ia_calc_bond_energy(ierr)
    !! Calculates total bond energy. Will overwrite `energy_bond` in module
    !! `m_globals`.

    integer, intent(out) :: ierr
    real(rp), dimension(3) :: ri, rj
    real(rp) :: enrg
    integer :: ibond
    integer :: typ, i, j

    ierr = 0
    energy_bond = 0.0_rp

    do ibond = 1, num_bonds
        typ = bonds(1,ibond)
        i = bonds(2,ibond)
        j = bonds(3,ibond)
        ri = coordinates(:,i)
        rj = coordinates(:,j)
        call ia_get_bond_energy(ri, rj, bond_specs(typ), enrg, ierr)
        if (ierr /= 0) return
        energy_bond = energy_bond + enrg
    end do

    end subroutine
    
!*******************************************************************************

subroutine ia_calc_atm_bond_energy(iatm, enrg_atm, ierr)
    !! Calculates the bond energy of a single atom.

    integer, intent(in) :: iatm
    real(rp), intent(out) :: enrg_atm
    integer, intent(out) :: ierr
    integer, dimension(:), pointer :: bnds_iatm => null()
    real(rp), dimension(3) :: ri, rj
    real(rp) :: enrg
    integer :: ibond, typ, i, j, k

    ierr = 0
    enrg_atm = 0.0_rp

    !Bond energy
    !bnds_iatm: Pointer to ids of bonds incident to atom iatm
    call atbo_tab%get_row(iatm, bnds_iatm)
    !Calculate energies due to bonds listed in bnds_iatm
    do k = 1, size(bnds_iatm)
        ibond = bnds_iatm(k)
        typ = bonds(1,ibond)
        i = bonds(2,ibond); j = bonds(3,ibond)
        ri = coordinates(:,i); rj = coordinates(:,j)
        call ia_get_bond_energy(ri, rj, bond_specs(typ), enrg, ierr)
        if (ierr /= 0) return
        enrg_atm = enrg_atm + enrg
    end do

    end subroutine

!********************************************************************************

subroutine ia_calc_angle_energy()
    !! Calculates total angular energy. Will overwrite `energy_angle` in module
    !! `m_globals`.

    real(rp), dimension(3) :: rim1
    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: rip1
    real(rp) :: enrg
    integer :: iangle
    integer :: typ
    integer :: i, im1, ip1

    energy_angle = 0.0_rp

    do iangle = 1, num_angles
        typ = angles(1,iangle)
        im1 = angles(2, iangle)
        i = angles(3, iangle)
        ip1 = angles(4, iangle)
        rim1 = coordinates(:,im1)
        ri = coordinates(:,i)
        rip1 = coordinates(:,ip1)
        call ia_get_angle_energy(rim1, ri, rip1, angle_specs(typ), enrg)
        energy_angle = energy_angle + enrg
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_calc_atm_angle_energy(iatm, enrg_atm)
    !! Calculates the angle energy due to a single atom.

    integer, intent(in) :: iatm
    real(rp), intent(out) :: enrg_atm
    integer, dimension(:), pointer :: angs_iatm => null()
    real(rp), dimension(3) :: ri, rip1, rim1
    real(rp) :: enrg
    integer :: iangle
    integer :: typ
    integer :: i, k, ip1, im1

    enrg_atm = 0.0_rp

    !Angle energy
    !angs_iatm: Pointer to ids of angles incident to atom iatm
    call atan_tab%get_row(iatm, angs_iatm)
    !Calculate energies due to angles listed in angs_iatm
    do k = 1, size(angs_iatm)
        iangle = angs_iatm(k)
        typ = angles(1,iangle)
        im1 = angles(2,iangle); i = angles(3,iangle); ip1 = angles(4,iangle)
        rim1 = coordinates(:,im1)
        ri = coordinates(:,i)
        rip1 = coordinates(:,ip1)
        call ia_get_angle_energy(rim1, ri, rip1, angle_specs(typ), enrg)
        enrg_atm = enrg_atm + enrg
    end do

    end subroutine

!********************************************************************************

subroutine ia_calc_dihedral_energy()
    !! Calculates total dihedral energy. Will overwrite `energy_dihedral` in module
    !! `m_globals`.

    real(rp), dimension(3) :: ri, rj, rk, rl
    real(rp) :: enrg
    integer :: idhd
    integer :: typ
    integer :: i, j, k, l

    energy_dihedral = 0.0_rp

    do idhd = 1, num_dihedrals
        typ = dihedrals(1,idhd)
        i = dihedrals(2, idhd)
        j = dihedrals(3, idhd)
        k = dihedrals(4, idhd)
        l = dihedrals(5, idhd)

        ri = coordinates(:,i)
        rj = coordinates(:,j)
        rk = coordinates(:,k)
        rl = coordinates(:,l)

        call ia_get_dihedral_energy(ri, rj, rk, rl, dihedral_specs(typ), enrg)
        energy_dihedral = energy_dihedral + enrg
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_calc_atm_dihedral_energy(iatm, enrg_atm)
    !! Calculates the dihedral energy due to a single atom.

    integer, intent(in) :: iatm
    real(rp), intent(out) :: enrg_atm
    integer, dimension(:), pointer :: dhds_iatm => null()
    real(rp), dimension(3) :: ri, rj, rk, rl
    real(rp) :: enrg
    integer :: idhd
    integer :: typ
    integer :: i, j, k, l, m

    enrg_atm = 0.0_rp

    !Dihedral energy
    !dhds_iatm: Pointer to ids of dihedrals incident to atom iatm
    call atdh_tab%get_row(iatm, dhds_iatm)

    !Calculate energies due to dihedrals listed in dhds_iatm
    do m = 1, size(dhds_iatm)
        idhd = dhds_iatm(k)
        typ = dihedrals(1,idhd)

        i = dihedrals(2, idhd)
        j = dihedrals(3, idhd)
        k = dihedrals(4, idhd)
        l = dihedrals(5, idhd)

        ri = coordinates(:,i)
        rj = coordinates(:,j)
        rk = coordinates(:,k)
        rl = coordinates(:,l)

        call ia_get_dihedral_energy(ri, rj, rk, rl, dihedral_specs(typ), enrg)
        enrg_atm = enrg_atm + enrg
    end do

    end subroutine

!********************************************************************************

subroutine ia_calc_tether_energy(ierr)
    !! Calculates total tether energy. Will overwrite `energy_tether` in module
    !! `m_globals`.

    integer, intent(out) :: ierr
    real(rp), dimension(3) :: r
    real(rp) :: enrg
    integer :: iteth, teth_iatm

    ierr = 0
    energy_tether = 0.0_rp

    do iteth = 1, num_tethers
        teth_iatm = tethers(iteth)%atm !Index of the tethered atom
        r = coordinates(:,teth_iatm)
        call ia_get_tether_energy(r, tethers(iteth), enrg, ierr)
        if (ierr /= 0) return
        energy_tether = energy_tether + enrg
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_calc_atm_tether_energy(iatm, enrg_atm, ierr)
    !! Calculates the tether energy due to a single atom.

    integer, intent(in) :: iatm
    real(rp), intent(out) :: enrg_atm
    integer, intent(out) :: ierr
    real(rp), dimension(3) :: r
    real(rp) :: enrg
    integer :: iteth, teth_iatm

    ierr = 0
    enrg_atm = 0.0_rp

    do iteth = 1, num_tethers
        teth_iatm = tethers(iteth)%atm !Index of the tethered atom
        if (iatm /= teth_iatm) cycle
        r = coordinates(:,iatm)
        call ia_get_tether_energy(r, tethers(iteth), enrg, ierr)
        if (ierr /= 0) return
        enrg_atm = enrg_atm + enrg
    end do

    end subroutine

!********************************************************************************

subroutine ia_calc_external_energy(ierr)
    !! Calculates total energy due to external fields. Will overwrite 
    !!`energy_external` in module `m_globals`.

    integer, intent(out) :: ierr
    real(rp), dimension(3) :: r
    real(rp) :: enrg
    integer :: iext

    ierr = 0
    energy_external = 0.0_rp

    do iext = 1, num_externals
        call ia_get_external_energy(extrn_fields(iext), coordinates, enrg, ierr)
        if (ierr /= 0) return
        energy_external = energy_external + enrg
    end do

    end subroutine
    
!******************************************************************************

end module m_interaction
