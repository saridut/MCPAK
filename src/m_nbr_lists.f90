module m_nbr_lists

use m_precision
use m_vector
use m_table
use m_globals

implicit none

private

public :: verlet_init, verlet_build, verlet_delete, &
            atbo_build, atan_build, atdh_build, exat_build
public :: verlet_tab, atbo_tab, atan_tab, atdh_tab, exat_tab

type(itable_t) :: verlet_tab
    !! Verlet table
type(itable_t) :: atbo_tab
    !! Atoms -> bonds table
type(itable_t) :: atan_tab
    !! Atoms -> angles table
type(itable_t) :: atdh_tab
    !! Atoms -> dihedrals table
type(itable_t) :: atat_tab
    !! Atoms -> bonded atoms table (1-ring)
type(itable_t) :: exat_tab
    !! Atoms -> excluded atoms (from vdw calculation) table
real(rp), dimension(:,:), allocatable :: coordinates_save
    !!  (3, *num_atoms*) array
real(rp), dimension(:,:), allocatable :: coordinates_dr
    !!  (3, *num_atoms*) array

real(rp) :: rskin_sq = 0.0_rp
real(rp) :: tskin_sq = 0.0_rp

contains

!******************************************************************************

subroutine verlet_init(rskin, tskin)

    real(rp), intent(in) :: rskin
    real(rp), intent(in) :: tskin

    rskin_sq = rskin**2
    tskin_sq = tskin**2

    call itbl_init(verlet_tab, num_atoms)

    allocate( coordinates_save(3,num_atoms) )
    allocate( coordinates_dr(3,num_atoms)   )

    !Initializing to zero
    coordinates_save = 0.0_rp
    coordinates_dr = 0.0_rp

    end subroutine

!******************************************************************************

subroutine verlet_delete()

    if (allocated(coordinates_save)) deallocate(coordinates_save)
    if (allocated(coordinates_dr)  ) deallocate(coordinates_dr)
    call verlet_tab%delete()

    end subroutine

!******************************************************************************

subroutine verlet_build()

    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: rj
    real(rp), dimension(3) :: rij
    real(rp) :: dr_sq_max
    real(rp) :: rij_sq
    logical, save :: first_call = .true.
    integer :: i, j

    !On first call no check for rebuilding
    if (.not. first_call) then
        !Check whether rebuilding the list is necessary
        coordinates_dr = coordinates - coordinates_save
        dr_sq_max = maxval(sum(coordinates_dr**2, dim=1))
        if ( 4*dr_sq_max < tskin_sq ) then
            return
        !else
        !    print*, nblks, nts, 4*dr_sq_max
        end if
    end if

    first_call = .false.

    !Clear table
    call verlet_tab%clear()

    !Loop over all pairs to build list (different from MD type table)
    do i = 1, num_atoms
        ri = coordinates(:,i)
        do j = 1, num_atoms
            if ( exat_tab%is_in(i,j) ) cycle
            rj = coordinates(:,j)
            rij = rj - ri
            rij_sq = sum(rij**2)
            if (rij_sq < rskin_sq) call verlet_tab%append(i,j)
        end do
    end do

    !Release additional memory
    call verlet_tab%shrink_to_fit()

    !Back up positions
    coordinates_save = coordinates

    end subroutine

!********************************************************************************

subroutine atbo_build()

    type(ivector_t), dimension(:), allocatable :: buf_map
    integer :: iatm, jatm, ibnd
    integer :: i, j

    !Build table as list of lists
    allocate(buf_map(num_atoms))
    do iatm = 1, num_atoms
        call ivector_init(buf_map(iatm))
    end do

    do ibnd = 1, num_bonds
        iatm = bonds(2,ibnd)
        jatm = bonds(3,ibnd)
        call buf_map(iatm)%append(ibnd)
        call buf_map(jatm)%append(ibnd)
    end do

    !Sort in ascending order
    do iatm = 1, num_atoms
        if (buf_map(iatm)%len > 1) call buf_map(iatm)%sort()
    end do

    !Copy over to table
    call itbl_init(atbo_tab, num_atoms)
    do iatm = 1, num_atoms
        do j = 1, buf_map(iatm)%len
            jatm = buf_map(iatm)%get_val(j)
            call atbo_tab%append(iatm, jatm)
        end do
    end do

    !Delete list of lists table
    do iatm = 1, num_atoms
        call buf_map(iatm)%delete()
    end do
    deallocate(buf_map)

    !Release additional memory
    call atbo_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

subroutine atan_build()

    type(ivector_t), dimension(:), allocatable :: buf_map
    integer :: iatm, jatm, katm, iang
    integer :: i, j

    !Build table as list of lists
    allocate(buf_map(num_atoms))
    do iatm = 1, num_atoms
        call ivector_init(buf_map(iatm))
    end do

    do iang = 1, num_angles
        iatm = angles(2,iang)
        jatm = angles(3,iang)
        katm = angles(4,iang)
        call buf_map(iatm)%append(iang)
        call buf_map(jatm)%append(iang)
        call buf_map(katm)%append(iang)
    end do

    !Sort in ascending order
    do iatm = 1, num_atoms
        if (buf_map(iatm)%len > 1) call buf_map(iatm)%sort()
    end do

    !Copy over to table
    call itbl_init(atan_tab, num_atoms)
    do iatm = 1, num_atoms
        do j = 1, buf_map(iatm)%len
            jatm = buf_map(iatm)%get_val(j)
            call atan_tab%append(iatm, jatm)
        end do
    end do

    !Delete list of lists table
    do iatm = 1, num_atoms
        call buf_map(iatm)%delete()
    end do
    deallocate(buf_map)

    !Release additional memory
    call atan_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

subroutine atdh_build()

    type(ivector_t), dimension(:), allocatable :: buf_map
    integer :: iatm, jatm, katm, latm, idhd
    integer :: i, j

    !Build table as list of lists
    allocate(buf_map(num_atoms))
    do iatm = 1, num_atoms
        call ivector_init(buf_map(iatm))
    end do

    do idhd = 1, num_dihedrals
        iatm = dihedrals(2,idhd)
        jatm = dihedrals(3,idhd)
        katm = dihedrals(4,idhd)
        latm = dihedrals(5,idhd)
        call buf_map(iatm)%append(idhd)
        call buf_map(jatm)%append(idhd)
        call buf_map(katm)%append(idhd)
        call buf_map(latm)%append(idhd)
    end do

    !Sort in ascending order
    do iatm = 1, num_atoms
        if (buf_map(iatm)%len > 1) call buf_map(iatm)%sort()
    end do

    !Copy over to table
    call itbl_init(atdh_tab, num_atoms)
    do iatm = 1, num_atoms
        do j = 1, buf_map(iatm)%len
            jatm = buf_map(iatm)%get_val(j)
            call atdh_tab%append(iatm, jatm)
        end do
    end do

    !Delete list of lists table
    do iatm = 1, num_atoms
        call buf_map(iatm)%delete()
    end do
    deallocate(buf_map)

    !Release additional memory
    call atdh_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

subroutine atat_build()

    integer, dimension(:), pointer :: inc_bonds => null()
    type(ivector_t) :: excl_atms
    integer :: iatm, jatm, jbnd, jbnd_atm1, jbnd_atm2
    integer :: i, j

    !Initialize table
    call itbl_init(atat_tab, num_atoms)
    !List of excluded atoms
    call ivector_init(excl_atms)

    do iatm = 1, num_atoms
        call excl_atms%clear()
        !Get incident bonds
        call atbo_tab%get_row(iatm, inc_bonds)
        !Add atoms from each of the incident bonds, excluding iatom
        do j = 1, size(inc_bonds)
            jbnd = inc_bonds(j)
            jbnd_atm1 = bonds(2,jbnd)
            jbnd_atm2 = bonds(3,jbnd)
            if (jbnd_atm1 /= iatm) call excl_atms%append(jbnd_atm1)
            if (jbnd_atm2 /= iatm) call excl_atms%append(jbnd_atm2)
        end do
        !Sort and remove duplicates. These are the atoms in the 1-ring.
        call excl_atms%unique()

        !Add to atat_tab. Note that iatm does not appear in the atat table
        !for iatm.
        do j = 1, excl_atms%len
            jatm = excl_atms%get_val(j)
            call atat_tab%append(iatm, jatm)
        end do
    end do

    !Release additional memory
    call atat_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

subroutine exat_build()

    integer, dimension(:), pointer :: nbr_atms => null()
    integer, dimension(:), pointer :: nbr2_atms => null()
    integer, dimension(:), pointer :: nbr3_atms => null()
    type(ivector_t) :: excl_atms
    integer :: iatm, jatm, j2atm, j3atm
    integer :: i, j, j2, j3

    !Build atat_tab. This is the 1-ring neighborhood of each atom.
    call atat_build()
    call atat_tab%shrink_to_fit()

    !Debug statements
    !print*, 'ATAT_TAB'
    !call atat_tab%print()

    !Initialize table
    call itbl_init(exat_tab, num_atoms)
    !List of excluded atoms
    call ivector_init(excl_atms)

    do iatm = 1, num_atoms
        call excl_atms%clear()
        if (excluded_atoms > 0) then
            !Get atoms in 1-ring neighborhood
            call atat_tab%get_row(iatm, nbr_atms)
            !Add atoms to excluded atoms list. First adding iatm itself.
            call excl_atms%append(iatm)
            do j = 1, size(nbr_atms)
                jatm = nbr_atms(j)
                call excl_atms%append(jatm)
                !Atoms for second ring neighbors
                if (excluded_atoms > 1) then
                    call atat_tab%get_row(jatm, nbr2_atms)
                    do j2 = 1, size(nbr2_atms)
                        j2atm = nbr2_atms(j2)
                        call excl_atms%append(j2atm)
                        !Atoms for third ring neighbors
                        if (excluded_atoms > 2) then
                            call atat_tab%get_row(j2atm, nbr3_atms)
                            do j3 = 1, size(nbr3_atms)
                                j3atm = nbr3_atms(j3)
                                call excl_atms%append(j3atm)
                            end do
                            !Finished adding all third ring neighbors
                        end if
                    end do
                    !Finished adding all second ring neighbors
                end if
                !Finished adding all first ring neighbors
            end do
        end if

        !Sort and remove duplicates
        call excl_atms%unique()

        !Add to exat_tab. Note that iatom appears in the excluded atoms list for
        !iatom.
        do j = 1, excl_atms%len
            jatm = excl_atms%get_val(j)
            call exat_tab%append(iatm, jatm)
        end do
    end do

    !Release additional memory
    call exat_tab%shrink_to_fit()

    !atat_tab no longer required
    call atat_tab%delete()

    end subroutine

!********************************************************************************

end module m_nbr_lists
