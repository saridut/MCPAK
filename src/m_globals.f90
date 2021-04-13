module m_globals

use m_precision
use m_ia_types
use m_trajectory

implicit none

!General
integer, parameter :: mxrdln = 1024
    !! Maximum length of character string for input line buffer.
integer, parameter :: mx_mcmtyp = 5
    !! Maximum number of MC move types.
    !! Type 1: backbone pivot, type 2: side chain crankshaft, type 3: side chain
    !! pivot, type 4: backbone restricted pivot, type 5: displacement

!Simulation box
real(rp), dimension(3) :: smbx_a
    !! Lattice vector of simulation box along x-axis.
real(rp), dimension(3) :: smbx_b
    !! Lattice vector of simulation box along y-axis.
real(rp), dimension(3) :: smbx_c
    !! Lattice vector of simulation box along z-axis.
integer :: imcon = 0
    !! Flag specifying boundary condition on the simulation box.

!Particle configuration: Atoms
integer :: num_atom_types = 0
    !! Number of *atom_type*s
type(atm_specs_t), dimension(:), allocatable :: atom_specs
    !! (*num_atom_types*,) array.
integer, dimension(:), allocatable :: atom_pop
    !! (*num_atom_types*,) array. Population of atoms of each type.
integer :: num_atoms = 0
    !!  Number of atoms
integer , dimension(:), allocatable:: atoms
    !! (*num_atoms*,) array. For atom *i*, its type *at = atoms(i)*,
    !! charge *charge(i)*, position *coordinates(:,i)*,
    !! orientation (if the style requires) *orientations(:,i)*. 
real(rp), dimension(:), allocatable :: charge
    !! (*num_atoms*,) array.
real(rp), dimension(:,:), allocatable :: coordinates
    !! (3, *num_atoms*) array
real(rp), dimension(:,:), allocatable :: orientation
    !! (4, *num_atoms*) array

!Particle configuration: Bonds
integer :: num_bond_types = 0
    !!  Number of *bond_type*s
type(ia_specs_t), dimension(:), allocatable :: bond_specs
     !!  (*num_bond_types*,) array
integer :: num_bonds = 0
    !!  Total number of bonds.
integer, dimension(:,:), allocatable :: bonds
    !! (3, *num_bonds*) array. Bond *i* is of type *bt = bonds(1,i)*,  directed from
    !! atom *bonds(2,i)* to *bonds(3,i)*. 

!Particle configuration: Angles
integer :: num_angle_types = 0
    !!  Number of *angle_type*s
type(ia_specs_t), dimension(:), allocatable :: angle_specs
     !!  (*num_angle_types*,) array
integer :: num_angles = 0
    !!  Number of angles
integer, dimension(:,:), allocatable :: angles
    !! (4, *num_angles*) array. Angle *i* is of type *ant = angles(1,i)*, incident
    !! to atoms *angles(2,i)*, *angles(3,i)*, and *angles(4,i)*.

!Particle configuration: Dihedrals
integer :: num_dihedral_types = 0
    !!  Number of *dihedral_type*s
type(ia_specs_t), dimension(:), allocatable :: dihedral_specs
     !!  (*num_dihedral_types*,) array
integer :: num_dihedrals = 0
    !!  Number of dihedrals
integer, dimension(:,:), allocatable :: dihedrals
    !! (5, *num_dihedrals*) array. Dihedral *i* is of type *dt = dihedrals(1,i)*, incident
    !! to atoms *dihedrals(2,i)*, *dihedrals(3,i)*, *dihedrals(4,i)*, and *dihedrals(5,i)*.

!Particle configuration: Branches
integer :: num_branches = 0
    !! Total number of branches (including the backbone)
integer , dimension(:,:), allocatable:: branches
    !! (3,*num_branches*) array. Branch *i* is tethered to atom *branches(1,i)*,
    !! contains *branches(2,i)* atoms, with the beginning atom index *branches(3,i)*.

!Particle configuration: Molecules
integer :: num_molecule_types = 0
    !!  Number of *molecule_type*s
character(len=8), dimension(:), allocatable :: molecule_names
    !! (*num_molecule_types*,) array
integer, dimension(:), allocatable :: molecule_pop
    !! (*num_molecule_types*,) array
integer :: num_molecules = 0
    !!  Number of molecules
integer , dimension(:,:), allocatable:: molecules
    !! (9,*num_molecules*) array. For molecule *i*, its type *mt = molecules(1,i)*, 
    !! containing *molecules(2,i)* atoms with beginning index *molecules(3,i)*,
    !! *molecules(4,i)* bonds with beginning index *molecules(5,i)*,
    !! *molecules(6,i)* angles with beginning index *molecules(7,i)*, and
    !! *molecules(8,i)* dihedrals with beginning index *molecules(9,i)*.

!Particle configuration: Tethers
integer :: num_tethers = 0
    !!  Number of tethers
type(tether_t), dimension(:), allocatable :: tethers
     !!  (*num_tethers*,) array

!Particle configuration: VDW (pair) interactions
integer :: num_vdw_types = 0
    !!  Number of *vdw_type*s
integer, dimension(:,:), allocatable :: vdw_pairs
    !!  (2, *num_vdw_types*) array. Stores atom type of interacting pairs, such
    !! that at_i >= at_j.
type(ia_specs_t), dimension(:), allocatable :: vdw_specs
     !!  (*num_vdw_types*,) array.

!Particle configuration: External field
integer :: num_externals = 0
    !!  Number of external fields
type(extrn_field_t), dimension(:), allocatable :: extrn_fields
    !!  (*num_externals*,) array

!End of configuration related globals

!Variables controlling runtime behavior
logical :: leql = .true.
    !! Is the system equilibrating? {T, F}
logical :: lrevive = .false.
    !! Is this a restart run? {T, F}.
integer :: nts
    !! Counter for MC cycles
integer :: nblks
    !! Counter for MC blocks
integer :: nts_log = 1
    !! Interval for logging (in MC cycles)
integer :: nts_dump = 1
    !! Interval for dumping to revive file (in MC cycles)
integer :: nts_eql = 0
    !! Number of MC cycles for equilibration
integer :: nts_eql_samp = 1
    !! Sampling interval during equilibration (in MC cycles)
integer :: block_size = 100
    !! Number of MC cycles in each block
integer :: nblks_sim = 0
!! Total number of blocks in production run
logical :: use_mc_pvt = .false.
    !! Use MC pivot moves? {T, F}
logical :: use_mc_res_pvt = .false.
    !! Use MC restricted pivot moves? {T, F}
logical :: use_mc_dspl = .false.
    !! Use MC displacement moves? {T, F}.
logical :: use_mc_crnk = .false.
    !! Use MC crankshaft moves? {T, F}. If side chains are present, crankshaft
    !! moves will operate only on side chain atoms. Otherwise, they will be
    !! performed on the backbone atoms.
logical :: use_mc_sc_pvt = .false.
    !! Use MC pivot (only side chains) moves? {T, F}
logical :: use_verlet_tab = .false.
    !! Use Verlet neighbor table? {T, F}
real(rp) :: mcm_mxdspl = 0.0_rp
    !!  Maximum displacement for MC displacement move.
real(rp) :: rcutoff = 0.0_rp
    !!  Maximum cutoff for short-ranged forces
real(rp) :: tskin = 0.0_rp
    !!  Thickness of the skin sphere (same for all)

!End of Variables controlling runtime behavior

!Variables for I/O
character(len=:), allocatable :: fn_cfg
    !! Name of the file containing the initial configuration
character(len=:), allocatable :: fn_revive
    !! Name of the revive file
character(len=:), allocatable :: fn_traj
    !! Name of the trajectory file
character(len=:), allocatable :: fn_stats
    !! Name of the statistics file
type(trajectory_t) :: traj
    !! Trajectory object
character(len=8) :: job_tag = ''
    !! A tag useful for array jobs, available only as a command line argument
logical :: read_seed = .false.
    !! {T, F}
    !!
    !! Whether to initialize the random number generator by reading a seed from
    !! a file. If `read_seed` == T, the seed will be read from a file
    !! 'random_seed.txt'
logical :: write_seed = .false.
    !! {T, F}
    !!
    !!  Whether to write the random number generator seed. If
    !!  `write_seed` == T the seed will be written to a file named
    !!  'random_seed.txt'
logical :: write_eql_stats = .false.
    !! During equilibration, should the statistics file be written? {T, F}
logical :: write_traj = .false.
    !! Should the trajectory be written to file? {T, F}

!End of variables for I/O

!Miscellaneous variables
real(rp) :: energy_bond = 0.0_rp
    !! Bond energy
real(rp) :: energy_angle = 0.0_rp
    !! Angle energy
real(rp) :: energy_dihedral = 0.0_rp
    !! Dihedral energy
real(rp) :: energy_vdw = 0.0_rp
    !! vdw interaction energy
real(rp) :: energy_tether = 0.0_rp
    !! Energetic contribution from tethers
real(rp) :: energy_external = 0.0_rp
    !! Energetic contribution from external fields
integer, dimension(mx_mcmtyp) :: num_mc_moves = 0
    !! Number of MC moves of each type in an MC cycle. 
integer(ip_long), dimension(2,mx_mcmtyp) :: mc_rec = 0
    !! Number of attempts (first element) & accepts (second element) of MC
    !! moves. See above for type details.
integer :: excluded_atoms = 0
    !! Control for excluded atoms in vdw calculation.
    !! No exclusion: 0, exclude 1-ring bonded neighbors: 1,
    !! exclude 2-ring bonded neighbors: 2, exclude 3-ring bonded neighbors: 3.
logical :: lvdw = .false.
    !! Whether to calculate VDW interactions
logical :: lelectrostatics = .false.
    !! Whether to calculate electrostatic interactions

!End of miscellaneous variables

end module m_globals
