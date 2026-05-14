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
    integer :: num_quant_f           ! Number of f-orbitals
  end type element_t

  !----------------------------------------------------------------
  ! Derived type for the potential block
  !----------------------------------------------------------------
  type :: potential_t
    integer :: lmax                  ! Maximum angular momentum
    real(rp), allocatable :: pl(:, :)             ! (l, spin) Potential parameters
    real(rp), allocatable :: ql(:, :, :)          ! (moment, l, spin) Moments
    real(rp), allocatable :: center_band(:, :)    ! (l, spin) Center of the band
    real(rp), allocatable :: shifted_band(:, :)   ! (l, spin) Shifted band parameters
    real(rp), allocatable :: width_band(:, :)     ! (l, spin) Width of the band
    real(rp) :: ws_r                 ! Wigner-Seitz radius
    real(rp) :: sumec                ! Sum of core energies
    real(rp) :: sumev                ! Sum of valence energies
    real(rp) :: etot                 ! Total energy
    real(rp) :: utot                 ! Total potential energy
    real(rp) :: ekin                 ! Kinetic energy
    real(rp) :: rhoeps               ! Electron density energy
    real(rp) :: vmad                 ! Madelung potential
    real(rp), allocatable :: enu(:, :)            ! (l, spin) Orthogonal potential parameters
    real(rp), allocatable :: c(:, :)              ! (l, spin) Coefficients for potential
    real(rp), allocatable :: srdel(:, :)          ! (l, spin) Screening parameters
    real(rp), allocatable :: qpar(:, :)           ! (l, spin) Quadrupole parameters
    real(rp), allocatable :: ppar(:, :)           ! (l, spin) Potential parameters
    real(rp), allocatable :: obar(:, :)           ! (l, spin) Overlap parameters
    real(rp), allocatable :: vl(:, :)             ! (l, spin) Local potential parameters
    real(rp), allocatable :: screening_alpha(:)   ! Screening alpha per l
    real(rp), allocatable :: screening_alpha_dot(:) ! Screening alpha dot per l
    real(rp), dimension(2) :: vmtz = 0.0_rp
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
    logical*1 :: b1                  ! Boundary condition 1
    logical*1 :: b2                  ! Boundary condition 2
    logical*1 :: b3                  ! Boundary condition 3
    integer :: n1                    ! Grid points 1
    integer :: n2                    ! Grid points 2
    integer :: n3                    ! Grid points 3
    integer :: ndim                  ! Rough atom estimate
    real(rp), allocatable :: ct(:)   ! Cutoff radius per site
    real(rp) :: r2                   ! Max cutoff radius (squared)
    integer :: lmax = 2              ! Target maximum angular momentum for RS output
    character(len=16) :: strux_backend = ''
    character(len=16) :: screening = ''
    logical :: strux_want_sdot = .false.
    real(rp) :: strux_solve_scale = 1.0_rp
    real(rp), allocatable :: screening_alpha(:)
    real(rp), allocatable :: screening_sigma(:)
    logical :: has_strux_backend = .false.
    logical :: has_screening = .false.
    logical :: has_strux_want_sdot = .false.
    logical :: has_strux_solve_scale = .false.
    logical :: has_screening_alpha = .false.
    logical :: has_screening_sigma = .false.
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
    real(kind=8), allocatable :: pnu(:, :, :)
    real(kind=8), allocatable :: qnu(:, :, :, :)
    real(kind=8) :: z(1)
    real(kind=8) :: avw
    real(kind=8) :: wsr(1)
    real(kind=8), allocatable :: pp(:, :, :, :)
    real(kind=8), allocatable :: pdum(:, :, :)
    logical :: lafm, lfree
  end type lmto_t
end module system_data
