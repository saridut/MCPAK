module m_mc_moves

use m_precision
use m_constants_math
use m_ran_num
use m_globals
use m_nbr_lists
use m_interaction

implicit none

private
public :: mcm_dspl, mcm_crnk, mcm_sc_crnk, mcm_pivot, mcm_dbl_pivot, &
        mcm_sc_pivot, mcm_res_pivot, coordinates_pvt, to_com

!Memory buffer for pivot moves (allocated & deallocated in module m_mc_solver)
real(rp), dimension(:,:), allocatable :: coordinates_pvt

contains
 
!******************************************************************************

subroutine mcm_dspl()
    !! Performs displacement move.

    real(rp), dimension(3) :: coords
    real(rp), dimension(3) :: rdspl
    real(rp) :: enrg_vdw, enrg_bnd, enrg_ang, enrg_dhd, enrg_teth
    real(rp) :: enrg_vdw_, enrg_bnd_, enrg_ang_, enrg_dhd_, enrg_teth_
    real(rp) :: energy_external_
    real(rp) :: enrg_dif
    integer  :: iatm
    integer  :: natmpts, naccept
    integer  :: iatmpt, iaccept
    integer  :: ierr

    !Number of attempts
    natmpts = num_mc_moves(5)

    !Initialize number of accepts to zero
    naccept = 0

    do iatmpt = 1, natmpts
        !Pick an atom
        iatm = get_iuniform(1, num_atoms+1)
        !Choose a random unit vector
        call ransphere(rdspl)
        !Scale the random unit vector to generate a displacement vector
        rdspl = rdspl*get_uniform(0.0_rp, mcm_mxdspl)

        !Get energy based on old configuration. For a good configuration, ierr
        !will return 0.
        enrg_vdw_ = 0.0_rp; enrg_bnd_ = 0.0_rp
        enrg_ang_ = 0.0_rp; enrg_dhd_ = 0.0_rp; enrg_teth_ = 0.0_rp
        energy_external_ = energy_external

        if (lvdw) then
            if (use_verlet_tab) then
                call ia_calc_atm_vdw_energy_vl(iatm, enrg_vdw_, ierr)
            else
                call ia_calc_atm_vdw_energy(iatm, enrg_vdw_, ierr)
            end if
        end if
        if (num_bonds > 0) call ia_calc_atm_bond_energy(iatm, enrg_bnd_, ierr)
        if (num_angles > 0) call ia_calc_atm_angle_energy(iatm, enrg_ang_)
        if (num_dihedrals > 0) call ia_calc_atm_dihedral_energy(iatm, enrg_dhd_)
        if (num_tethers > 0) call ia_calc_atm_tether_energy(iatm, enrg_teth_, ierr)
        if (num_externals > 0) call ia_calc_external_energy(ierr)

        !Save atom position
        coords = coordinates(:,iatm)

        !Displace atom
        coordinates(:,iatm) = coordinates(:,iatm) + rdspl

        !Get energy based on new configuration
        enrg_vdw = 0.0_rp; enrg_bnd = 0.0_rp
        enrg_ang = 0.0_rp; enrg_dhd = 0.0_rp; enrg_teth = 0.0_rp

        if (lvdw) then
            if (use_verlet_tab) then
                call ia_calc_atm_vdw_energy_vl(iatm, enrg_vdw, ierr)
            else
                call ia_calc_atm_vdw_energy(iatm, enrg_vdw, ierr)
            end if
        end if

        if ((ierr==0) .and. (num_bonds > 0)) then
            call ia_calc_atm_bond_energy(iatm, enrg_bnd, ierr)
        end if
        if ((ierr==0) .and. (num_angles > 0)) then
            call ia_calc_atm_angle_energy(iatm, enrg_ang)
        end if
        if ((ierr==0) .and. (num_dihedrals > 0)) then
            call ia_calc_atm_dihedral_energy(iatm, enrg_dhd)
        end if
        if ((ierr==0) .and. (num_tethers > 0)) then
            call ia_calc_atm_tether_energy(iatm, enrg_teth, ierr)
        end if
        if ((ierr==0) .and. (num_externals > 0)) then
            call ia_calc_external_energy(ierr)
        end if

        if (ierr /= 0) then
            enrg_dif = huge(0.0_rp)
        else
            enrg_dif = (enrg_vdw - enrg_vdw_) + (enrg_ang - enrg_ang_) &
                + (enrg_ang - enrg_ang_) + (enrg_dhd - enrg_dhd_)      &
                + (enrg_teth - enrg_teth_)                             &
                + (energy_external - energy_external_)
        end if

        iaccept = metro_crit(enrg_dif)
        naccept = naccept + iaccept
        if (iaccept == 0) then
            !If move is not accepted, revert position
            coordinates(:,iatm) =  coords
            energy_external = energy_external_
        else
            !If move is accepted, update energy
            energy_bond = energy_bond + enrg_bnd - enrg_bnd_
            energy_angle = energy_angle + enrg_ang - enrg_ang_
            energy_dihedral = energy_dihedral + enrg_dhd - enrg_dhd_
            energy_vdw = energy_vdw + enrg_vdw - enrg_vdw_
            energy_tether = energy_tether + enrg_teth - enrg_teth_
            if (num_tethers == 0) call to_com()
        end if
    end do

    mc_rec(1,5) = mc_rec(1,5) + natmpts
    mc_rec(2,5) = mc_rec(2,5) + naccept

    end subroutine

!******************************************************************************

subroutine mcm_crnk(is_ring)
    !! Performs crankshaft move (including end-bond rotation) on an unbranched
    !! chain or ring. End-bond rotation will not be performed for rings.

    logical, intent(in) :: is_ring
    real(rp), dimension(3,3) :: rotmat
    real(rp), dimension(3) :: coords
    real(rp), dimension(3) :: axis
    real(rp), dimension(3) :: r_pvt
    real(rp), dimension(3) :: ristar
    real(rp) :: axis_len
    real(rp) :: angle
    real(rp) :: enrg_vdw, enrg_ang, enrg_dhd, enrg_teth
    real(rp) :: enrg_vdw_, enrg_ang_, enrg_dhd_, enrg_teth_
    real(rp) :: energy_external_
    real(rp) :: enrg_dif
    integer  :: iatm
    integer  :: natmpts, naccept
    integer  :: iatmpt, iaccept
    integer  :: crnk_beg, crnk_end
    integer  :: ierr
    integer  :: na_bbone

    !No branches: Total number of atoms same as the number of backbone atoms.
    na_bbone = num_atoms

    !There must be two or more atoms
    if (na_bbone < 2) return

    !Number of attempts
    natmpts = num_mc_moves(2)

    !Initialize number of accepts to zero
    naccept = 0

    do iatmpt = 1, natmpts
        !Pick an atom
        if (is_ring) then
            !For rings: No ends
            iatm = get_iuniform(1, na_bbone+1)
            if (iatm == 1) then
                crnk_beg = na_bbone; crnk_end = 2
            else if (iatm == na_bbone) then
                crnk_beg = na_bbone - 1; crnk_end = 1
            else
                crnk_beg = iatm - 1; crnk_end = iatm + 1
            end if
            axis = coordinates(:,crnk_end) - coordinates(:,crnk_beg)
        else
            !For chains
            if (num_tethers == 0) then
                iatm = get_iuniform(1, na_bbone+1)
            else
                iatm = get_iuniform(2, na_bbone+1) !<-- Tether constraints here
            end if
            if (iatm == 1) then
                !Ends are rotated about a random axis
                call ransphere(axis)
                crnk_beg = 2
            else if (iatm == na_bbone) then
                !Ends are rotated about a random axis
                call ransphere(axis)
                crnk_beg = na_bbone - 1
            else
                !Non-ends are rotated about an axis determined by adjacent atom
                !positions
                crnk_beg = iatm - 1; crnk_end = iatm + 1
                axis = coordinates(:,crnk_end) - coordinates(:,crnk_beg)
            end if
        end if

        axis_len = norm2(axis)
        if (axis_len < 1.0E-8_rp) then
            !If axis_len is too small reject the move
            cycle
        else
            !Get unit vector along the axis
            axis = axis/axis_len
        end if

        angle = get_uniform(-math_pi, math_pi)
        call get_rotmat(axis, angle, rotmat)

        r_pvt = coordinates(:,crnk_beg) !Pivot point

        !Get energy based on old configuration. For a good configuration, ierr
        !will return 0.
        enrg_vdw_ = 0.0_rp; enrg_ang_ = 0.0_rp
        enrg_dhd_ = 0.0_rp; enrg_teth_ = 0.0_rp
        energy_external_ = energy_external

        if (lvdw) then
            if (use_verlet_tab) then
                call ia_calc_atm_vdw_energy_vl(iatm, enrg_vdw_, ierr)
            else
                call ia_calc_atm_vdw_energy(iatm, enrg_vdw_, ierr)
            end if
        end if
        if (num_angles > 0) call ia_calc_atm_angle_energy(iatm, enrg_ang_)
        if (num_dihedrals > 0) call ia_calc_atm_dihedral_energy(iatm, enrg_dhd_)
        if (num_tethers > 0) call ia_calc_atm_tether_energy(iatm, enrg_teth_, ierr)
        if (num_externals > 0) call ia_calc_external_energy(ierr)

        !Save atom position
        coords = coordinates(:,iatm)

        !Rotate atom
        ristar = coordinates(:,iatm) - r_pvt
        coordinates(:,iatm) = r_pvt + rotmat(:,1)*ristar(1) &
                + rotmat(:,2)*ristar(2) + rotmat(:,3)*ristar(3)

        !Get energy based on new configuration
        enrg_vdw = 0.0_rp; enrg_ang = 0.0_rp
        enrg_dhd = 0.0_rp; enrg_teth = 0.0_rp

        if (lvdw) then
            if (use_verlet_tab) then
                call ia_calc_atm_vdw_energy_vl(iatm, enrg_vdw, ierr)
            else
                call ia_calc_atm_vdw_energy(iatm, enrg_vdw, ierr)
            end if
        end if

        if ((ierr==0) .and. (num_angles > 0)) then
            call ia_calc_atm_angle_energy(iatm, enrg_ang)
        end if
        if ((ierr==0) .and. (num_dihedrals > 0)) then
            call ia_calc_atm_dihedral_energy(iatm, enrg_dhd)
        end if
        if ((ierr==0) .and. (num_tethers > 0)) then
            call ia_calc_atm_tether_energy(iatm, enrg_teth, ierr)
        end if
        if ((ierr==0) .and. (num_externals > 0)) then
            call ia_calc_external_energy(ierr)
        end if

        if (ierr /= 0) then
            enrg_dif = huge(0.0_rp)
        else
            enrg_dif = (enrg_vdw - enrg_vdw_) + (enrg_ang - enrg_ang_) &
                + (enrg_dhd - enrg_dhd_) + (enrg_teth - enrg_teth_)    &
                + (energy_external - energy_external_)
        end if

        iaccept = metro_crit(enrg_dif)
        naccept = naccept + iaccept
        if (iaccept == 0) then
            !If move is not accepted, revert position
            coordinates(:,iatm) =  coords
            energy_external = energy_external_
        else
            !If move is accepted, update energy
            energy_angle = energy_angle + enrg_ang - enrg_ang_
            energy_dihedral = energy_dihedral + enrg_dhd - enrg_dhd_
            energy_vdw = energy_vdw + enrg_vdw - enrg_vdw_
            energy_tether = energy_tether + enrg_teth - enrg_teth_
            if (num_tethers == 0) call to_com()
        end if
    end do

    mc_rec(1,2) = mc_rec(1,2) + natmpts
    mc_rec(2,2) = mc_rec(2,2) + naccept

    end subroutine

!******************************************************************************

subroutine mcm_sc_crnk()
    !! Performs side chain crankshaft move. 
     
    real(rp), dimension(3,3) :: rotmat
    real(rp), dimension(3) :: coords
    real(rp), dimension(3) :: axis
    real(rp), dimension(3) :: r_pvt
    real(rp), dimension(3) :: ristar
    real(rp) :: axis_len
    real(rp) :: angle
    real(rp) :: enrg_vdw, enrg_ang, enrg_dhd, enrg_teth
    real(rp) :: enrg_vdw_, enrg_ang_, enrg_dhd_, enrg_teth_
    real(rp) :: energy_external_
    real(rp) :: enrg_dif
    integer  :: iatm
    integer  :: natmpts, naccept
    integer  :: iatmpt, iaccept
    integer  :: crnk_beg, crnk_end
    integer  :: ierr
    integer  :: na_br, ia_br_beg, ia_br, ibr

    !Number of attempts
    natmpts = num_mc_moves(2)
    !Initialize number of accepts to zero
    naccept = 0

    do iatmpt = 1, natmpts
        !Pick a side chain
        ibr = get_iuniform(2, num_branches+1) !Ignore the backbone
        na_br = branches(2,ibr)
        ia_br_beg = branches(3,ibr)

        !Pick an atom on this side chain and crankshaft axis
        if (na_br == 1) then
            !There is only one atom on this side chain
            iatm = ia_br_beg
            call ransphere(axis)
            crnk_beg = branches(1,ibr)
        else
            !There are multiple atoms on this side chain
            iatm = get_iuniform(ia_br_beg, ia_br_beg+na_br)
            if (iatm == ia_br_beg) then
                !First side chain atom
                crnk_beg = branches(1,ibr)
                crnk_end = iatm + 1
                axis = coordinates(:,crnk_end) - coordinates(:,crnk_beg)
            else if (iatm == (ia_br_beg+na_br-1)) then
                !Last side chain atom
                call ransphere(axis)
                crnk_beg = iatm - 1
            else
                !Intermediate side chain atom
                crnk_beg = iatm - 1; crnk_end = iatm + 1
                axis = coordinates(:,crnk_end) - coordinates(:,crnk_beg)
            end if
        end if

        axis_len = norm2(axis)
        if (axis_len < 1.0E-8_rp) then
            !If axis_len is too small reject the move
            cycle
        else
            !Get unit vector along the axis
            axis = axis/axis_len
        end if

        angle = get_uniform(-math_pi, math_pi)
        call get_rotmat(axis, angle, rotmat)

        r_pvt = coordinates(:,crnk_beg) !Pivot point

        !Get energy based on old configuration. For a good configuration, ierr
        !will return 0.
        enrg_vdw_ = 0.0_rp; enrg_ang_ = 0.0_rp
        enrg_dhd_ = 0.0_rp; enrg_teth_ = 0.0_rp
        energy_external_ = energy_external

        if (lvdw) then
            if (use_verlet_tab) then
                call ia_calc_atm_vdw_energy_vl(iatm, enrg_vdw_, ierr)
            else
                call ia_calc_atm_vdw_energy(iatm, enrg_vdw_, ierr)
            end if
        end if
        if (num_angles > 0) call ia_calc_atm_angle_energy(iatm, enrg_ang_)
        if (num_dihedrals > 0) call ia_calc_atm_dihedral_energy(iatm, enrg_dhd_)
        if (num_tethers > 0) call ia_calc_atm_tether_energy(iatm, enrg_teth_, ierr)
        if (num_externals > 0) call ia_calc_external_energy(ierr)

        !Save atom position
        coords = coordinates(:,iatm)

        !Rotate atom
        ristar = coordinates(:,iatm) - r_pvt
        coordinates(:,iatm) = r_pvt + rotmat(:,1)*ristar(1) &
                + rotmat(:,2)*ristar(2) + rotmat(:,3)*ristar(3)
    
        !Get energy based on new configuration
        enrg_vdw = 0.0_rp; enrg_ang = 0.0_rp
        enrg_dhd = 0.0_rp; enrg_teth = 0.0_rp

        if (lvdw) then
            if (use_verlet_tab) then
                call ia_calc_atm_vdw_energy_vl(iatm, enrg_vdw, ierr)
            else
                call ia_calc_atm_vdw_energy(iatm, enrg_vdw, ierr)
            end if
        end if

        if ((ierr==0) .and. (num_angles > 0)) then
            call ia_calc_atm_angle_energy(iatm, enrg_ang)
        end if
        if ((ierr==0) .and. (num_dihedrals > 0)) then
            call ia_calc_atm_dihedral_energy(iatm, enrg_dhd)
        end if
        if ((ierr==0) .and. (num_tethers > 0)) then
            call ia_calc_atm_tether_energy(iatm, enrg_teth, ierr)
        end if
        if ((ierr==0) .and. (num_externals > 0)) then
            call ia_calc_external_energy(ierr)
        end if

        if (ierr /= 0) then
            enrg_dif = huge(0.0_rp)
        else
            enrg_dif = (enrg_vdw - enrg_vdw_) + (enrg_ang - enrg_ang_) &
                + (enrg_dhd - enrg_dhd_) + (enrg_teth - enrg_teth_)    &
                + (energy_external - energy_external_)
        end if

        iaccept = metro_crit(enrg_dif)
        naccept = naccept + iaccept
        if (iaccept == 0) then
            !If move is not accepted, revert position
            coordinates(:,iatm) =  coords
            energy_external = energy_external_
        else
            !If move is accepted, update energy
            energy_angle = energy_angle + enrg_ang - enrg_ang_
            energy_dihedral = energy_dihedral + enrg_dhd - enrg_dhd_
            energy_vdw = energy_vdw + enrg_vdw - enrg_vdw_
            energy_tether = energy_tether + enrg_teth - enrg_teth_
            if (num_tethers == 0) call to_com()
        end if
    end do

    mc_rec(1,2) = mc_rec(1,2) + natmpts
    mc_rec(2,2) = mc_rec(2,2) + naccept

    end subroutine

!******************************************************************************

subroutine mcm_pivot()
    !! Performs a pivot move. Does not use a Verlet list.
    !I have left some print statements for debugging purposes.

    integer :: iaccept
    real(rp), dimension(3,3) :: rotmat
    real(rp), dimension(3) :: axis
    real(rp), dimension(3) :: r_pvt
    real(rp), dimension(3) :: ristar
    real(rp) :: angle
    real(rp) :: enrg_ang_, enrg_ang, enrg_dhd_, enrg_dhd
    real(rp) :: energy_vdw_, energy_external_
    real(rp) :: enrg_dif
    integer  :: pvt
    integer  :: ierr
    integer  :: na_bbone, na_br, tp
    integer  :: i, jbeg, jend, ibr, ia_br, ia_br_beg 

    iaccept = 0; ierr = 0
    coordinates_pvt = coordinates 

    if (num_branches == 0) then
        na_bbone = num_atoms
    else
        na_bbone = branches(2,1)
    end if

    call ransphere(axis)
    angle = get_uniform(-math_pi, math_pi)
    call get_rotmat(axis, angle, rotmat)

    !Pick a pivot point on the backbone. The section to be rotated cannot
    !have any tethering.
    if (num_tethers == 0) then
        !For a free molecule, rotate the shorter part
        pvt = get_iuniform(2, na_bbone)
        if (pvt < na_bbone/2 ) then
            jbeg = 1; jend = pvt - 1     !Will be rotated
        else
            jbeg = pvt + 1; jend = na_bbone !Will be rotated
        end if
    else
        !Change this part to suit the tether location. Here assuming that atom 1
        !is tethered to a wall, rotate the section away from the wall.
        pvt = get_iuniform(1, na_bbone)
        jbeg = pvt + 1; jend = na_bbone !Will be rotated
    end if
    r_pvt = coordinates(:,pvt)
    !print*, 'pvt ', pvt
    !print*, 'jbeg:jend ', jbeg, jend

    !Energy before rotation. This includes only the contribution from the pivot
    !point atom. The total vdw and total external energies are already known
    !from the last iteration.
    enrg_ang_ = 0.0_rp; enrg_dhd_ = 0.0_rp
    energy_vdw_ = energy_vdw           !Caching
    energy_external_ = energy_external !Caching
    if (num_angles > 0) call ia_calc_atm_angle_energy(pvt, enrg_ang_)
    if (num_dihedrals > 0) call ia_calc_atm_dihedral_energy(pvt, enrg_dhd_)

    !Rotate backbone
    do i = jbeg, jend
       ristar = coordinates(:,i) - r_pvt
       coordinates(:,i) = r_pvt + rotmat(:,1)*ristar(1) &
                       + rotmat(:,2)*ristar(2)          &
                       + rotmat(:,3)*ristar(3)
    end do 
    !Rotate side chains
    if (num_branches > 0) then
        do ibr = 2, num_branches !Backbone is the first branch
            tp = branches(1,ibr) !Index of tether point
            if ( (tp < jbeg) .or. (tp > jend) ) cycle
            na_br = branches(2,ibr) !Number of atoms in branch
            ia_br_beg = branches(3,ibr) !Index of the beginning atom in branch
            do ia_br = 1, na_br
                i = ia_br_beg + ia_br - 1
                ristar = coordinates(:,i) - r_pvt
                coordinates(:,i) = r_pvt + rotmat(:,1)*ristar(1) &
                        + rotmat(:,2)*ristar(2) + rotmat(:,3)*ristar(3)
            end do
        end do 
    end if

    !Total energy after rotation
    enrg_ang = 0.0_rp; enrg_dhd = 0.0_rp

    call ia_calc_external_energy(ierr)
    if ((ierr == 0) .and. (lvdw)) then
        call ia_calc_vdw_energy(ierr)
    end if
    if ((ierr == 0) .and. (num_angles > 0)) then
        call ia_calc_atm_angle_energy(pvt, enrg_ang)
    end if
    if ((ierr == 0) .and. (num_dihedrals > 0)) then
        call ia_calc_atm_dihedral_energy(pvt, enrg_dhd)
    end if

    !Energy difference
    if (ierr /= 0) then
        enrg_dif = huge(0.0_rp)
    else
        enrg_dif = (energy_vdw - energy_vdw_) + (enrg_ang - enrg_ang_) &
            + (enrg_dhd - enrg_dhd_) + (energy_external - energy_external_)
    end if

    iaccept = metro_crit(enrg_dif)
    if (iaccept == 0) then
        !If the move is rejected
        coordinates = coordinates_pvt
        energy_vdw = energy_vdw_
        energy_external = energy_external_
    else
        energy_angle = energy_angle + enrg_ang - enrg_ang_
        energy_dihedral = energy_dihedral + enrg_dhd - enrg_dhd_
        if (num_tethers == 0) call to_com()
    end if

    mc_rec(1,1) = mc_rec(1,1) + 1
    mc_rec(2,1) = mc_rec(2,1) + iaccept

    end subroutine

!******************************************************************************

subroutine mcm_res_pivot()
    !!Performs a restricted pivot move on the backbone. Does not use a Verlet list.
    !!Unlike usual pivot, the rotation axis for is chosen along the
    !!bond incident to the pivot point. This move may be useful for molecules that
    !!are restricted to be in a stretched conformation or if angular constraints
    !!need to be respected.
    !I have left some print statements for debugging purposes.

    integer :: iaccept
    real(rp), dimension(3,3) :: rotmat
    real(rp), dimension(3) :: axis
    real(rp), dimension(3) :: r_pvt
    real(rp), dimension(3) :: ristar
    real(rp) :: angle
    real(rp) :: enrg_ang_, enrg_ang, enrg_dhd_, enrg_dhd
    real(rp) :: energy_vdw_, energy_external_
    real(rp) :: enrg_dif
    integer  :: pvt
    integer  :: ierr
    integer  :: na_bbone, na_br, tp
    integer  :: i, jbeg, jend, ibr, ia_br, ia_br_beg 

    iaccept = 0; ierr = 0
    coordinates_pvt = coordinates 

    if (num_branches == 0) then
        na_bbone = num_atoms
    else
        na_bbone = branches(2,1)
    end if

    !Pick a pivot point on the backbone. The section to be rotated cannot
    !have any tethering.
    if (num_tethers == 0) then
        !For a free molecule, rotate the shorter part
        pvt = get_iuniform(2, na_bbone)
        if (pvt < na_bbone/2 ) then
            jbeg = 1; jend = pvt - 1     !Will be rotated
            axis = coordinates(:,pvt) - coordinates(:,pvt+1)
            axis = axis/norm2(axis)
        else
            jbeg = pvt + 1; jend = na_bbone !Will be rotated
            axis = coordinates(:,pvt) - coordinates(:,pvt-1)
            axis = axis/norm2(axis)
        end if
    else
        !Change this part to suit the tether location. Here, assuming that atom 1
        !is tethered to a wall, rotate the section away from the wall.
        pvt = get_iuniform(1, na_bbone)
        jbeg = pvt + 1; jend = na_bbone !Will be rotated
        axis = coordinates(:,pvt) - coordinates(:,pvt-1)
        axis = axis/norm2(axis)
    end if
    r_pvt = coordinates(:,pvt)
    !print*, 'pvt ', pvt
    !print*, 'jbeg:jend ', jbeg, jend

    !Choose the angle of rotation
    angle = get_uniform(-math_pi, math_pi)
    call get_rotmat(axis, angle, rotmat)

    !Energy before rotation. This includes only the contribution from the pivot
    !point atom. The total vdw and total external energies are already known
    !from the last iteration.
    enrg_ang_ = 0.0_rp; enrg_dhd_ = 0.0_rp
    energy_vdw_ = energy_vdw           !Caching
    energy_external_ = energy_external !Caching
    if (num_angles > 0) call ia_calc_atm_angle_energy(pvt, enrg_ang_)
    if (num_dihedrals > 0) call ia_calc_atm_dihedral_energy(pvt, enrg_dhd_)

    !Rotate backbone
    do i = jbeg, jend
       ristar = coordinates(:,i) - r_pvt
       coordinates(:,i) = r_pvt + rotmat(:,1)*ristar(1) &
                       + rotmat(:,2)*ristar(2)          &
                       + rotmat(:,3)*ristar(3)
    end do 
    !Rotate side chains
    if (num_branches > 0) then
        do ibr = 2, num_branches !Backbone is the first branch
            tp = branches(1,ibr) !Index of tether point
            if ( (tp < jbeg) .or. (tp > jend) ) cycle
            na_br = branches(2,ibr) !Number of atoms in branch
            ia_br_beg = branches(3,ibr) !Index of the beginning atom in branch
            do ia_br = 1, na_br
                i = ia_br_beg + ia_br - 1
                ristar = coordinates(:,i) - r_pvt
                coordinates(:,i) = r_pvt + rotmat(:,1)*ristar(1) &
                        + rotmat(:,2)*ristar(2) + rotmat(:,3)*ristar(3)
            end do
        end do 
    end if

    !Total energy after rotation
    enrg_ang = 0.0_rp; enrg_dhd = 0.0_rp

    call ia_calc_external_energy(ierr)
    if ((ierr == 0) .and. (lvdw)) then
        call ia_calc_vdw_energy(ierr)
    end if
    if ((ierr == 0) .and. (num_angles > 0)) then
        call ia_calc_atm_angle_energy(pvt, enrg_ang)
    end if
    if ((ierr == 0) .and. (num_dihedrals > 0)) then
        call ia_calc_atm_dihedral_energy(pvt, enrg_dhd)
    end if

    !Energy difference
    if (ierr /= 0) then
        enrg_dif = huge(0.0_rp)
    else
        enrg_dif = (energy_vdw - energy_vdw_) + (enrg_ang - enrg_ang_) &
            + (enrg_dhd - enrg_dhd_) + (energy_external - energy_external_)
    end if

    iaccept = metro_crit(enrg_dif)
    if (iaccept == 0) then
        !If the move is rejected
        coordinates = coordinates_pvt
        energy_vdw = energy_vdw_
        energy_external = energy_external_
    else
        energy_angle = energy_angle + enrg_ang - enrg_ang_
        energy_dihedral = energy_dihedral + enrg_dhd - enrg_dhd_
        if (num_tethers == 0) call to_com()
    end if

    mc_rec(1,4) = mc_rec(1,4) + 1
    mc_rec(2,4) = mc_rec(2,4) + iaccept

    end subroutine

!******************************************************************************

subroutine mcm_dbl_pivot()
    !!Performs double pivot move (for rings). Does not use Verlet list.
    !It is assumed that there are no tethers. If there are tethers, it is better
    !to incorporate them via energetic penalty.
    !I have left some print statements for debugging purposes.

    integer :: iaccept
    real(rp), dimension(3,3) :: rotmat
    real(rp), dimension(3) :: axis
    real(rp), dimension(3) :: r_pvt
    real(rp), dimension(3) :: ristar
    real(rp) :: angle
    real(rp) :: enrg_ang_, enrg_ang, enrg_dhd_, enrg_dhd
    real(rp) :: energy_vdw_, energy_external_
    real(rp) :: enrg, enrg_dif
    real(rp) :: axis_len
    integer  :: ierr
    integer  :: na_bbone, na_br, tp
    integer  :: ax_beg, ax_end, tmp, mstps, jstps
    integer  :: i, j, jj, jbeg, ibr, ia_br, ia_br_beg

    iaccept = 0; ierr = 0
    coordinates_pvt = coordinates 

    if (num_branches == 0) then
        na_bbone = num_atoms
    else
        na_bbone = branches(2,1)
    end if

    ax_beg = get_iuniform(1, na_bbone+1)
!   print*, 'ax_beg ', ax_beg
    do
        mstps = get_iuniform(2, na_bbone-1)
        ax_end = mod(ax_beg+mstps-1, na_bbone) + 1
        axis = coordinates(:,ax_end) - coordinates(:,ax_beg)
        axis_len = norm2(axis)
        !Check if axis_len is too small
        if (axis_len > 1.0E-8_rp) exit
    end do
!   print*, 'mstps: ax_end ', mstps, ax_end
    !Swap ax_beg and ax_end such that ax_end > ax_beg
    if (ax_end < ax_beg) then
        tmp = ax_end; ax_end = ax_beg; ax_beg = tmp
    end if
!   print*, 'ax_beg:ax_end ', ax_beg, ax_end

    !Get unit vector along the axis
    axis = axis/axis_len

    angle = get_uniform(-math_pi, math_pi)
    call get_rotmat(axis, angle, rotmat)
    r_pvt = coordinates(:,ax_beg) !First pivot point

    !Find the shorter part of the molecule to rotate. Note that the atoms ax_end
    !and ax_beg (i.e. atoms on the rotation axis) are not to be rotated.
    if ((ax_end-ax_beg) < na_bbone/2 ) then
        !Will be rotated
        jbeg = ax_beg + 1; jstps = (ax_end-1) - (ax_beg+1)     
    else
        !Will be rotated
        jbeg = ax_end + 1; jstps = na_bbone - (ax_end+1) + (ax_beg-1) 
    end if
!   print*, 'jbeg:jstps ', jbeg, jstps

    !Energy before rotation. This includes only the contribution from the pivot
    !point atom. The total vdw and total external energies are already known
    !from the last iteration.
    enrg_ang_ = 0.0_rp; enrg_dhd_ = 0.0_rp
    energy_vdw_ = energy_vdw           !Caching
    energy_external_ = energy_external !Caching
    if (num_angles > 0) then
        call ia_calc_atm_angle_energy(ax_beg, enrg)
        enrg_ang_ = enrg
        call ia_calc_atm_angle_energy(ax_end, enrg)
        enrg_ang_ = enrg_ang_ + enrg
    end if
    if (num_dihedrals > 0) then
        call ia_calc_atm_dihedral_energy(ax_beg, enrg)
        enrg_dhd_ = enrg
        call ia_calc_atm_dihedral_energy(ax_end, enrg)
        enrg_dhd_ = enrg_dhd_ + enrg
    end if

    !Rotate backbone
    do jj = 0, jstps
        j = mod(jbeg+jj-1, na_bbone) + 1
        ristar = coordinates(:,j) - r_pvt
        coordinates(:,j) = r_pvt + rotmat(:,1)*ristar(1) &
                        + rotmat(:,2)*ristar(2)          &
                        + rotmat(:,3)*ristar(3)
        !Rotate side chains
        if (num_branches > 0) then
            do ibr = 2, num_branches !Backbone is the first branch
                tp = branches(1,ibr) !Index of tether point
                if (tp /= j) cycle
                na_br = branches(2,ibr) !Number of atoms in branch
                ia_br_beg = branches(3,ibr) !Index of the beginning atom in branch
                do ia_br = 1, na_br
                    i = ia_br_beg + ia_br - 1
                    ristar = coordinates(:,i) - r_pvt
                    coordinates(:,i) = r_pvt + rotmat(:,1)*ristar(1) &
                            + rotmat(:,2)*ristar(2) + rotmat(:,3)*ristar(3)
                end do
            end do 
        end if
    end do 

    !Total energy after rotation
    enrg_ang = 0.0_rp; enrg_dhd = 0.0_rp

    call ia_calc_external_energy(ierr)

    if ((ierr == 0) .and. (lvdw)) then
        call ia_calc_vdw_energy(ierr)
    end if

    if ((ierr == 0) .and. (num_angles > 0)) then
        call ia_calc_atm_angle_energy(ax_beg, enrg)
        enrg_ang = enrg
        call ia_calc_atm_angle_energy(ax_end, enrg)
        enrg_ang = enrg_ang + enrg
    end if

    if ((ierr == 0) .and. (num_dihedrals > 0)) then
        call ia_calc_atm_dihedral_energy(ax_beg, enrg)
        enrg_dhd = enrg
        call ia_calc_atm_dihedral_energy(ax_end, enrg)
        enrg_dhd = enrg_dhd + enrg
    end if

    !Energy difference
    if (ierr /= 0) then
        enrg_dif = huge(0.0_rp)
    else
        enrg_dif = (energy_vdw - energy_vdw_) + (enrg_ang - enrg_ang_) &
            + (enrg_dhd - enrg_dhd_) + (energy_external - energy_external_)
    end if

    iaccept = metro_crit(enrg_dif)
    if (iaccept == 0) then
        !If the move is rejected
        coordinates = coordinates_pvt
        energy_vdw = energy_vdw_
        energy_external = energy_external_
    else
        energy_angle = energy_angle + enrg_ang - enrg_ang_
        energy_dihedral = energy_dihedral + enrg_dhd - enrg_dhd_
        call to_com()
    end if

    mc_rec(1,1) = mc_rec(1,1) + 1
    mc_rec(2,1) = mc_rec(2,1) + iaccept

    end subroutine

!******************************************************************************

subroutine mcm_sc_pivot()
    !! Performs pivot move on side chains. Does not use Verlet list.

    real(rp), dimension(3,3) :: rotmat
    real(rp), dimension(3) :: axis
    real(rp), dimension(3) :: r_pvt
    real(rp), dimension(3) :: ristar
    real(rp) :: angle
    real(rp) :: enrg_ang_, enrg_ang, enrg_dhd_, enrg_dhd
    real(rp) :: energy_vdw_, energy_external_
    real(rp) :: enrg_dif
    integer  :: pvt
    integer  :: ierr
    integer  :: natmpts, naccept
    integer  :: iatmpt, iaccept
    integer  :: i, jbeg, jend, ibr, ia_br, ia_br_beg, na_br 

    !Number of attempts
    natmpts = num_mc_moves(3)
    !Initialize number of accepts to zero
    naccept = 0

    do iatmpt = 1, natmpts
        coordinates_pvt = coordinates
        !Pick a side chain
        ibr = get_iuniform(2, num_branches) !Ignore the backbone
        na_br = branches(2,ibr)
        ia_br_beg = branches(3,ibr)

        call ransphere(axis)
        angle = get_uniform(-math_pi, math_pi)
        call get_rotmat(axis, angle, rotmat)

        !Pick a pivot point on this side chain
        if (na_br == 1) then
            !There is only one atom on this side chain. Pivot point
            !is the backbone atom to which this side chain is tethered to.
            pvt = branches(1,ibr)
            jbeg = ia_br_beg; jend = ia_br_beg + na_br - 1
        else
            !There are multiple atoms on this side chain
            pvt = get_iuniform(0, na_br)
            if (pvt == 0) then
                !Pivot point is the backbone atom to which this side chain
                !is tethered to.
                pvt = branches(1,ibr)
                jbeg = ia_br_beg; jend = ia_br_beg + na_br - 1
            else
                pvt = pvt + ia_br_beg - 1
                jbeg = pvt + 1; jend = ia_br_beg + na_br - 1
            end if
        end if

        r_pvt = coordinates(:,pvt) !Pivot point

        !Energy before rotation. This includes only the contribution from the pivot
        !point atom. The total vdw and total external energies are already known
        !from the last iteration.
        enrg_ang_ = 0.0_rp; enrg_dhd_ = 0.0_rp
        energy_vdw_ = energy_vdw           !Caching
        energy_external_ = energy_external !Caching
        if (num_angles > 0) call ia_calc_atm_angle_energy(pvt, enrg_ang_)
        if (num_dihedrals > 0) call ia_calc_atm_dihedral_energy(pvt, enrg_dhd_)

        !Rotate 
        do i = jbeg, jend
           ristar = coordinates(:,i) - r_pvt
           coordinates(:,i) = r_pvt + rotmat(:,1)*ristar(1) &
                           + rotmat(:,2)*ristar(2)          &
                           + rotmat(:,3)*ristar(3)
        end do

        !Total energy after rotation
        enrg_ang = 0.0_rp; enrg_dhd = 0.0_rp

        call ia_calc_external_energy(ierr)
        if ((ierr == 0) .and. (lvdw)) then
            call ia_calc_vdw_energy(ierr)
        end if
        if ((ierr == 0) .and. (num_angles > 0)) then
            call ia_calc_atm_angle_energy(pvt, enrg_ang)
        end if
        if ((ierr == 0) .and. (num_dihedrals > 0)) then
            call ia_calc_atm_dihedral_energy(pvt, enrg_dhd)
        end if

        !Energy difference
        if (ierr /= 0) then
            enrg_dif = huge(0.0_rp)
        else
            enrg_dif = (energy_vdw - energy_vdw_) + (enrg_ang - enrg_ang_) &
                + (enrg_dhd - enrg_dhd_) + (energy_external - energy_external_)
        end if

        iaccept = metro_crit(enrg_dif)
        naccept = naccept + iaccept
        if (iaccept == 0) then
            !If move is not accepted, revert position
            coordinates = coordinates_pvt
            energy_vdw = energy_vdw_
            energy_external = energy_external_
        else
            !If move is accepted, update energy
            energy_angle = energy_angle + enrg_ang - enrg_ang_
            energy_dihedral = energy_dihedral + enrg_dhd - enrg_dhd_
            if (num_tethers == 0) call to_com()
        end if
    end do

    mc_rec(1,3) = mc_rec(1,3) + natmpts
    mc_rec(2,3) = mc_rec(2,3) + naccept

    end subroutine

!******************************************************************************

function metro_crit(enrg_dif) result (res)
    !! Applies Metropolis criterion to a given energy difference. Returns zero
    !! for rejection and one for acceptance.

    real(rp), intent(in) :: enrg_dif
    integer  :: res
    real(rp) :: ran
    real(rp) :: prob

    res = 0

    if (enrg_dif <= 0.0_rp) then
        res = 1
    else if ((enrg_dif >= 0.0_rp) .and. (enrg_dif <= 40.0_rp)) then
        !Avoid considering very high energy differences. Here, any move with
        !a energy difference of more than 40 is rejected.
        ran = get_uniform(0.0_rp, 1.0_rp)
        prob = exp(-enrg_dif)
        if (ran <= prob) then
            res = 1
        end if
    end if

    end function

!******************************************************************************

subroutine get_rotmat(axis, alpha, rotmat)
    !! Returns the rotation matrix for axis-angle

    real(rp), dimension(3,3), intent(out):: rotmat
    real(rp), dimension(3), intent(in) :: axis
    real(rp), intent(in) :: alpha
    real(rp)  :: sinalpha
    real(rp)  :: cosalpha
    real(rp)  :: icosalpha
    real(rp)  :: icosxy
    real(rp)  :: icosyz
    real(rp)  :: icoszx
    real(rp)  :: sinx
    real(rp)  :: siny
    real(rp)  :: sinz
    
    sinalpha  = sin(alpha)
    cosalpha  = cos(alpha)
    icosalpha = 1.0_rp - cosalpha
    
    icosxy = axis(1)*axis(2)*icosalpha 
    icosyz = axis(2)*axis(3)*icosalpha 
    icoszx = axis(3)*axis(1)*icosalpha 
    sinx   = axis(1)*sinalpha
    siny   = axis(2)*sinalpha
    sinz   = axis(3)*sinalpha
    
    rotmat(1,1) =  cosalpha + icosalpha*axis(1)*axis(1)
    rotmat(2,1) =  sinz + icosxy
    rotmat(3,1) =  -siny + icoszx
    rotmat(1,2) =  icosxy - sinz
    rotmat(2,2) =  cosalpha + icosalpha*axis(2)*axis(2)
    rotmat(3,2) =  sinx + icosyz
    rotmat(1,3) =  siny + icoszx
    rotmat(2,3) =  -sinx + icosyz
    rotmat(3,3) =  cosalpha + icosalpha*axis(3)*axis(3)
    
    end subroutine

!******************************************************************************

subroutine to_com()
    !!Brings the center-of-mass of the molecule to the origin.

    real(rp), dimension(3) :: cm
    integer :: i

    cm = sum(coordinates,2)/num_atoms

    do i = 1, num_atoms
        coordinates(:,i) = coordinates(:,i) - cm
    end do

    end subroutine

!******************************************************************************

end module m_mc_moves
