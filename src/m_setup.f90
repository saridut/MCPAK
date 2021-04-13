module m_setup

use m_precision
use m_ran_num
use m_globals
use m_config_io
use m_stats_io, only: stats_init, stats_finish
use m_mc_solver, only: mcs_init, mcs_finish
use m_interaction, only: ia_setup
use m_logger

implicit none

contains

!*******************************************************************************

subroutine setup()

    integer :: ierr
    logical :: lexists

    if ( lrevive ) then
        !Restarting simulation
        !Read revive file
        call read_dump(fn_revive//trim(adjustl(job_tag)))
        !Read tabulated potentials, if any
        call read_tabulated(ierr)


        !Append (write) to existing (new) trajectory file
        if (write_traj) then
            !Check if trajectory file exists
            inquire(file=fn_traj//trim(adjustl(job_tag)), exist=lexists)
            if (lexists) then
                !Open existing file for appending
                call traj%open(fn_traj//trim(adjustl(job_tag)), 'rw')
            else
                !Create new file
                call traj%create(fn_traj//trim(adjustl(job_tag)), num_atoms)
            end if
        end if
    else 
        !New simulation
        !Read initial configuration file
        call read_config(fn_cfg//trim(adjustl(job_tag)))
        !Read tabulated potentials, if any
        call read_tabulated(ierr)
        leql = .true.
        nblks = 0; nts = 0

        !Create new trajectory file
        if (write_traj) then
            call traj%create(fn_traj//trim(adjustl(job_tag)), num_atoms)
        end if
    end if

    !Initialize random number generator
    if (read_seed) then
        call init_stream('random_seed.txt'//trim(adjustl(job_tag)))
    else 
        call init_stream('')
    end if
    if (write_seed) call save_seed('random_seed.txt'//trim(adjustl(job_tag)))

    !Initialize stats collection
    call stats_init()

    !Set up MC solver
    call mcs_init()

    !Set up interactions
    if (num_vdw_types == 0) lvdw = .false.
    call ia_setup()

    end subroutine

!*******************************************************************************

subroutine config_clear()
    !! Clears out all configuration related variables in module `m_globals`.

    smbx_a = 0.0_rp
    smbx_b = 0.0_rp
    smbx_c = 0.0_rp
    imcon = 0

    num_atom_types = 0; num_atoms = 0
    if (allocated(atom_specs))  deallocate(atom_specs)
    if (allocated(atom_pop  ))  deallocate(atom_pop)
    if (allocated(atoms      )) deallocate(atoms)
    if (allocated(charge     )) deallocate(charge)
    if (allocated(coordinates)) deallocate(coordinates)

    num_bond_types = 0; num_bonds = 0
    if (allocated(bond_specs)) deallocate(bond_specs)
    if (allocated(bonds      )) deallocate(bonds)

    num_angle_types = 0; num_angles = 0
    if (allocated(angle_specs)) deallocate(angle_specs)
    if (allocated(angles      )) deallocate(angles)

    num_dihedral_types = 0; num_dihedrals = 0
    if (allocated(dihedral_specs)) deallocate(dihedral_specs)
    if (allocated(dihedrals      )) deallocate(dihedrals)

    num_branches = 0
    if (allocated(branches)) deallocate(branches)

    num_molecule_types = 0; num_molecules = 0
    if (allocated(molecule_names)) deallocate(molecule_names)
    if (allocated(molecule_pop  )) deallocate(molecule_pop)
    if (allocated(molecules     )) deallocate(molecules)

    num_tethers = 0
    if (allocated(tethers)) deallocate(tethers)

    num_vdw_types = 0
    if (allocated(vdw_specs)) deallocate(vdw_specs)
    if (allocated(vdw_pairs )) deallocate(vdw_pairs)

    num_externals = 0
    if (allocated(extrn_fields)) deallocate(extrn_fields)

    end subroutine

!*******************************************************************************

subroutine finish()

    logical :: lopen

    call traj%close()
    call mcs_finish()
    call stats_finish()
    call config_clear()
    call logger%finish()

    end subroutine

!*******************************************************************************

end module m_setup
