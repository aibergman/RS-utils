module system_data
  implicit none
  public

  ! Define a real kind with high precision
  integer, parameter :: rp = selected_real_kind(15, 307)

  !----------------------------------------------------------------
  ! Derived type for the lattice file
  !----------------------------------------------------------------
  type :: lattice_t
    ! Declare integer variables
    integer :: nbulk_bulk
    integer :: nbulk_org
    integer :: ntot
    integer :: ntype
    integer :: nbas
    integer :: nrec

    ! Declare allocatable integer arrays
    integer, allocatable :: izp(:)
    integer, allocatable :: no(:)
    integer, allocatable :: iu(:)
    integer, allocatable :: ib(:)
    integer, allocatable :: irec(:)

    ! Declare allocatable real arrays and variables
    real(rp), allocatable :: a(:, :)   ! Lattice vectors
    real(rp), allocatable :: crd(:, :) ! Basis vectors
    real(rp), allocatable :: ct(:)
    real(rp), allocatable :: ws(:)
    real(rp) :: r2

  end type lattice_t
  !----------------------------------------------------------------
  ! Derived type for the &element block
  !----------------------------------------------------------------
  type :: element_t
    character(len=8) :: symbol       ! Element symbol
    integer :: atomic_number         ! Atomic number
    integer :: core                  ! Core electrons
    integer :: valence               ! Valence electrons
    integer :: f_core                ! Core electrons in f-orbitals
    integer :: num_quant_s           ! Number of s-orbitals
    integer :: num_quant_p           ! Number of p-orbitals
    integer :: num_quant_d           ! Number of d-orbitals
  end type element_t

  !----------------------------------------------------------------
  ! Derived type for the potential block
  !----------------------------------------------------------------
  type :: potential_t
    integer :: lmax                  ! Maximum angular momentum
    real(rp) :: pl(3, 2)             ! Potential parameters
    real(rp) :: ql(3, 3, 2)          ! Quadrupole moments
    real(rp) :: center_band(3, 2)    ! Center of the band
    real(rp) :: shifted_band(3, 2)   ! Shifted band parameters
    real(rp) :: width_band(3, 2)     ! Width of the band
    real(rp) :: ws_r                 ! Wigner-Seitz radius
    real(rp) :: sumec                ! Sum of core energies
    real(rp) :: sumev                ! Sum of valence energies
    real(rp) :: etot                 ! Total energy
    real(rp) :: utot                 ! Total potential energy
    real(rp) :: ekin                 ! Kinetic energy
    real(rp) :: rhoeps               ! Electron density energy
    real(rp) :: vmad                 ! Madelung potential
    real(rp) :: enu(3, 2)            ! Orthogonal potential parameters
    real(rp) :: c(3, 2)              ! Coefficients for potential
    real(rp) :: srdel(3, 2)          ! Screening parameters
    real(rp) :: qpar(3, 2)           ! Quadrupole parameters
    real(rp) :: ppar(3, 2)           ! Potential parameters
    real(rp) :: obar(3, 2)           ! Overlap parameters
    real(rp) :: vl(3, 2)             ! Local potential parameters
  end type potential_t

  !-------------------------------------------------------------------
  ! Derived type for the &lattice namelist in input.nml
  !-------------------------------------------------------------------
  type :: lattice_inp
    integer :: ntype                 ! Number of types
    integer :: rc                    ! Radius cut
    real(kind=8) :: alat             ! Lattice parameter
    character(len=8) :: crystal_sym  ! Crystal symmetry ('bcc', 'fcc', 'file', etc.)
    real(kind=8) :: wav              ! Average WS radius
    logical*1 :: pbc                 ! Periodic BC
    integer :: ndim                  ! Rough atom estimate
    real(rp), allocatable :: ct(:)   ! Cutoff radius per site
    real(rp) :: r2                   ! Max cutoff radius (squared)
  end type lattice_inp

  !-------------------------------------------------------------------
  ! Derived type for the &atoms namelist in input.nml
  !-------------------------------------------------------------------
  type :: atoms_inp
    character(len=256) :: database   ! Path to the database
    character(len=8), dimension(:), allocatable :: label     ! Array of atom labels (here fixed to one element)
    integer, dimension(:), allocatable :: atomic_number     ! Array of atom labels (here fixed to one element)
  end type atoms_inp

  type :: lmto_t
    integer :: konfig(0:3, 1)
    integer :: nclass, lmx(1), nl, nsp
    real(kind=8) :: pnu(0:3-1, 2, 1)
    real(kind=8) :: qnu(0:2, 0:3-1, 2, 1)
    real(kind=8) :: z(1)
    real(kind=8) :: avw
    real(kind=8) :: wsr(1)
    real(kind=8) :: pp(6, 0:3-1, 2, 1)
    real(kind=8) :: pdum(0:3-1, 2, 1)
    logical :: lafm, lfree
  end type lmto_t
end module system_data
