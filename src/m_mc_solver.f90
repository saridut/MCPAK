module m_mc_solver

use m_precision
use m_constants_math
use m_strings
use m_globals
use m_nbr_lists
use m_mc_moves
use m_interaction, only: ia_calc_energy
use m_stats_io
use m_config_io, only: write_dump
use m_logger

implicit none

private

public :: mcs_init, mcs_run, mcs_finish

contains
 
!******************************************************************************

subroutine mcs_init()

    !Number of moves per MC cycle & necessary memory buffer allocation

    if (use_mc_dspl) num_mc_moves(5) = num_atoms

    if (use_mc_pvt) then
        num_mc_moves(1) = 1
        allocate(coordinates_pvt(3, num_atoms))
    end if

    if (use_mc_crnk) then
        if (num_branches == 0) then
            num_mc_moves(2) = num_atoms
        else
            num_mc_moves(2) = num_atoms - branches(2,1)
        end if
    end if

    if (use_mc_sc_pvt) then
        if (num_branches == 0) then
            num_mc_moves(3) = 0
        else
            num_mc_moves(3) = num_branches
        end if
    end if

    if (use_mc_res_pvt) then
        num_mc_moves(4) = 1
        if (.not. (allocated(coordinates_pvt)) ) then
            allocate(coordinates_pvt(3, num_atoms))
        end if
    end if

    !Initialize and build atom -> bond, atom -> atom and next nearest bonded
    !neighbor tables
    call atbo_build()
    !Debug statements
    !write(*,*) 'ATBO_TAB'
    !call atbo_tab%print()

    if (num_angles > 0) call atan_build()
    if (num_dihedrals > 0) call atdh_build()
    call exat_build()
    !Debug statements
    !write(*,*) 'EXAT_TAB'
    !call exat_tab%print()

    !Initialize verlet list
    if (use_verlet_tab) call verlet_init(rcutoff+tskin, tskin)

    end subroutine

!******************************************************************************

subroutine mcs_finish()

    if (allocated(coordinates_pvt)) deallocate(coordinates_pvt)
    call atbo_tab%delete()
    if (num_angles > 0) call atan_tab%delete()
    if (num_dihedrals > 0) call atdh_tab%delete()
    call exat_tab%delete()
    if (use_verlet_tab) call verlet_delete()

    end subroutine

!******************************************************************************

subroutine mcs_run(ierr)

    integer, intent(out) :: ierr
    logical :: is_ring
    integer :: na_bbone

    ierr = 0

    !Bring the center-of-mass of the molecule to the origin
    call to_com()

    !Is this a ring molecule?
    if (num_branches == 0) then
        !Unbranched molecule
        if (num_atoms == num_bonds) then
            !Unbranched ring
            is_ring = .true.
        else
            !Unbranched chain
            is_ring = .false.
        end if
    else
        !Branched molecule
        na_bbone = branches(2,1)
        if ( (bonds(2,na_bbone)==na_bbone) .and. (bonds(2,na_bbone)==1) ) then
            !Branched ring
            is_ring = .true.
        else
            !Branched chain
            is_ring = .false.
        end if
    end if

    !Calculate total energy
    call ia_calc_energy(ierr)
    if (ierr /= 0) then
        call logger%log_msg('<mcs_run>: Bad initial configuration')
        return
    end if

    !Equilibration run
    if (leql) then
        do while (nts < nts_eql)
            call mcs_mc_cycle(is_ring)
            nts = nts + 1

            !Logging
            if (mod(nts,nts_log) == 0) then
                call logger%log_msg('nts: '//str_from_num(nts))
            end if

            !Dump revive file
            if (mod(nts,nts_dump)==0) then
                call write_dump(fn_revive//trim(adjustl(job_tag)))
            end if

            !Equilibration stats
            if (write_eql_stats) then
                if (mod(nts,nts_eql_samp)==0) call stats_write_eql(.false.)
            end if
        end do
        call logger%log_msg('Equilibration completed')
        call stats_finish()
        call write_dump(fn_revive//trim(adjustl(job_tag)))
        leql = .false.
        nts = 0
        call stats_init()
    end if

    !Production run
    if (.not. leql) then
        do while (nblks < nblks_sim)
            do while (nts < block_size)
                call mcs_mc_cycle(is_ring)
                nts = nts + 1
                call stats_collect()

                !Logging
                if (mod(nts,nts_log) == 0) then
                    call logger%log_msg('nblks: '//str_from_num(nblks)&
                        &//' nts: '//str_from_num(nts))
                end if

                !Dump revive file
                if (mod(nts,nts_dump)==0) then
                    call write_dump(fn_revive//trim(adjustl(job_tag)))
                end if
            end do
            nblks = nblks + 1
            !Write stats
            call stats_write(.false.)
            !Write traj
            if (write_traj) call traj%append_frame(nblks, coordinates)
            nts = 0
            !Dump revive file after a block as well
            if (mod(nts,nts_dump)==0) then
                call write_dump(fn_revive//trim(adjustl(job_tag)))
            end if
        end do
    end if

    end subroutine

!******************************************************************************

subroutine mcs_mc_cycle(is_ring)

    logical, intent(in) :: is_ring

    !Displacement
    if (use_mc_dspl) call mcm_dspl()

    !Crankshaft
    if (use_mc_crnk) then
        if (num_branches == 0) then
            !Unbranched molecule: crankshaft on backbone
            call mcm_crnk(is_ring)
        else
            !Branched molecule: Crankshaft on side chains only
            call mcm_sc_crnk()
        end if
    end if

    !Side chain pivot
    if ( (use_mc_sc_pvt) .and. (num_branches > 0) ) call mcm_sc_pivot()

    !Backbone pivot
    if (use_mc_pvt) then
        if (is_ring) then
            call mcm_dbl_pivot()
        else
            call mcm_pivot()
        end if
    end if

    !Backbone restricted pivot (not for rings)
    if ( (use_mc_res_pvt) .and. (.not. is_ring) ) call mcm_res_pivot()

    end subroutine

!******************************************************************************

end module m_mc_solver
