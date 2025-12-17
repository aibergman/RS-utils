  module rsutil_lib
    
  use system_data
  use rsutil_config
  use input_io
  use lattice_io
  use spheres
  use potpar_wrapper

  implicit none

  ! Derived type variables for the input.nml data.
  type(lattice_inp) :: lat_input
  type(atoms_inp)   :: atm_input
  type(lattice_t)   :: lat_data
  type(potential_t) :: pot_data

  real(rp), dimension(:,:), allocatable :: new_basis
  character(len=8), dimension(:), allocatable :: new_label
  integer, dimension(:), allocatable :: new_atomic_number

  integer :: i

  contains 
subroutine process_lattice_impl()
   implicit none
  !-----------------------------------------------------------------
  ! Read the lattice namelist from lattice.nml using lattice_io.
  !-----------------------------------------------------------------
  call read_lattice_namelist("lattice.nml", lat_data)

  !-----------------------------------------------------------------
  ! Read the input namelists from input.nml using input_io.
  ! This subroutine reads the &lattice and &atoms groups and sets
  ! the ntype variable.
  !-----------------------------------------------------------------
  call read_input_namelists("input.nml", lat_input, atm_input, lat_data)

  !-----------------------------------------------------------------
  ! Set PBC based on cluster type from config
  !-----------------------------------------------------------------
  if (trim(cluster_type) == "free") then
     lat_input%pbc = .false.
  else
     lat_input%pbc = .true.
  end if

  !-----------------------------------------------------------------
  ! Initialize the spheres module.
  ! This subroutine sets up necessary data structures and parameters
  ! for the spheres module.
  !-----------------------------------------------------------------
  call init_spheres()

  !-----------------------------------------------------------------
  ! Update the structure based on the input data.
  ! This subroutine modifies the lattice and atomic data according
  ! to the input parameters.
  !-----------------------------------------------------------------
  call update_structure(lat_input, atm_input, lat_data)

  !-----------------------------------------------------------------
  ! Expand and fill arrays with new lattice data.
  ! This subroutine processes the lattice data to generate new basis,
  ! labels, and atomic numbers.
  !-----------------------------------------------------------------
  print *,'Lattice parameter:' , lat_input%alat
  call expand_and_fill_arrays(lat_data, lat_input%alat)

  !-----------------------------------------------------------------
  ! Copy the updated data to the lattice_t type.
  ! This subroutine transfers the processed data into the lat_data
  ! variable of type lattice_t.
  !-----------------------------------------------------------------
  call copy_to_lattice_t(lat_data)

  !-----------------------------------------------------------------
  ! Write the updated lattice namelist to lattice_new.nml.
  ! This subroutine outputs the modified lattice data to a new
  ! namelist file.
  !-----------------------------------------------------------------
  call write_lattice_namelist("lattice_new.nml", lat_data)

  !-----------------------------------------------------------------
  ! Update the input namelists with new data.
  ! This subroutine writes the updated &lattice and &atoms groups
  ! to a new input namelist file.
  !-----------------------------------------------------------------
  call update_namelists("input.nml", "input_new.nml", lat_input, atm_input)

  call update_atoms(pot_data, atm_input, lat_data, lat_input)

end subroutine process_lattice_impl

end module rsutil_lib
