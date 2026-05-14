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

  call apply_setup_geometry_transforms(lat_input, lat_data)
  if (use_2d_mode) then
     lat_input%b3 = .false.
     lat_input%n3 = 1
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

  if (lmax_override >= 0) then
    call update_atoms(pot_data, atm_input, lat_data, lat_input, lmax_override)
  else
    call update_atoms(pot_data, atm_input, lat_data, lat_input, 2)
  end if

end subroutine process_lattice_impl

subroutine apply_setup_geometry_transforms(lat_inp, lat_dat)
  implicit none
  type(lattice_inp), intent(inout) :: lat_inp
  type(lattice_t), intent(inout) :: lat_dat
  real(rp) :: sigma_l_norm, sigma_b_norm, det_latt

  if (alat_scale_factor <= 0.0_rp) error stop "alat-scale must be > 0"
  if (lattice_dist_sigma < 0.0_rp) error stop "lattice-dist must be >= 0"
  if (basis_dist_sigma < 0.0_rp) error stop "basis-dist must be >= 0"
  if (z_iso_factor == 0.0_rp) error stop "z-iso-factor must be > 0 when provided"

  if (use_fixed_seed) call set_seed_from_integer(random_seed_value)

  lat_inp%alat = lat_inp%alat * alat_scale_factor
  if (lat_inp%alat <= 0.0_rp) error stop "Scaled alat became non-positive"

  if (use_2d_mode) call apply_2d_cell_mode(lat_inp, lat_dat)

  sigma_l_norm = lattice_dist_sigma / lat_inp%alat
  sigma_b_norm = basis_dist_sigma / lat_inp%alat

  if (sigma_l_norm > 0.0_rp) then
    call add_gaussian_noise_matrix(lat_dat%a, sigma_l_norm)
    det_latt = lat_dat%a(1,1)*(lat_dat%a(2,2)*lat_dat%a(3,3)-lat_dat%a(2,3)*lat_dat%a(3,2)) - &
               lat_dat%a(1,2)*(lat_dat%a(2,1)*lat_dat%a(3,3)-lat_dat%a(2,3)*lat_dat%a(3,1)) + &
               lat_dat%a(1,3)*(lat_dat%a(2,1)*lat_dat%a(3,2)-lat_dat%a(2,2)*lat_dat%a(3,1))
    if (abs(det_latt) < 1.0e-8_rp) error stop "Distorted lattice is near-singular"
  end if

  if (sigma_b_norm > 0.0_rp) then
    do i = 1, size(lat_dat%crd, 2)
      call add_gaussian_noise_vector(lat_dat%crd(:, i), sigma_b_norm)
    end do
  end if
end subroutine apply_setup_geometry_transforms

subroutine apply_2d_cell_mode(lat_inp, lat_dat)
  implicit none
  type(lattice_inp), intent(in) :: lat_inp
  type(lattice_t), intent(inout) :: lat_dat
  integer :: nsites, i, n_bottom, n_top
  real(rp) :: zmin, zmax, zspan, lz_ang, lz_norm, area_xy, nrm
  real(rp) :: normal(3)
  real(rp), parameter :: ztol = 1.0e-6_rp

  nsites = size(lat_dat%crd, 2)
  if (nsites <= 0) error stop "2d mode requires at least one site"

  zmin = minval(lat_dat%crd(3, :))
  zmax = maxval(lat_dat%crd(3, :))
  zspan = zmax - zmin

  if (z_iso_factor > 0.0_rp) then
    n_bottom = 0
    n_top = 0
    do i = 1, nsites
      if (abs(lat_dat%crd(3, i) - zmin) <= ztol) n_bottom = n_bottom + 1
      if (abs(lat_dat%crd(3, i) - zmax) <= ztol) n_top = n_top + 1
    end do
    if (n_bottom <= 0 .or. n_top <= 0) error stop "2d layer detection failed"
    area_xy = sqrt(sum(cross3(lat_dat%a(:, 1), lat_dat%a(:, 2))**2))
    if (area_xy <= 0.0_rp) error stop "Invalid in-plane area in 2d mode"
    lz_norm = zspan + z_iso_factor*sqrt(area_xy/real(n_bottom, rp)) + &
              z_iso_factor*sqrt(area_xy/real(n_top, rp))
  else
    lz_ang = zspan*lat_inp%alat + lat_inp%alat
    lz_norm = lz_ang/lat_inp%alat
  end if

  if (lz_norm <= 0.0_rp) error stop "Computed 2d z-length is non-positive"

  normal = cross3(lat_dat%a(:, 1), lat_dat%a(:, 2))
  nrm = sqrt(sum(normal**2))
  if (nrm <= 1.0e-12_rp) error stop "Cannot build 2d normal: a1 and a2 are collinear"
  normal = normal / nrm
  if (normal(3) < 0.0_rp) normal = -normal

  lat_dat%a(:, 3) = normal * lz_norm
end subroutine apply_2d_cell_mode

subroutine set_seed_from_integer(seed_value)
  implicit none
  integer, intent(in) :: seed_value
  integer :: n, i
  integer, allocatable :: seed(:)
  call random_seed(size=n)
  allocate(seed(n))
  do i = 1, n
    seed(i) = mod(abs(seed_value) + 104729*i, huge(1) - 1)
    if (seed(i) == 0) seed(i) = i
  end do
  call random_seed(put=seed)
  deallocate(seed)
end subroutine set_seed_from_integer

subroutine add_gaussian_noise_matrix(mat, sigma)
  implicit none
  real(rp), intent(inout) :: mat(:, :)
  real(rp), intent(in) :: sigma
  integer :: i, j
  do j = 1, size(mat, 2)
    do i = 1, size(mat, 1)
      mat(i, j) = mat(i, j) + sigma*gaussian_rand()
    end do
  end do
end subroutine add_gaussian_noise_matrix

subroutine add_gaussian_noise_vector(vec, sigma)
  implicit none
  real(rp), intent(inout) :: vec(:)
  real(rp), intent(in) :: sigma
  integer :: i
  do i = 1, size(vec)
    vec(i) = vec(i) + sigma*gaussian_rand()
  end do
end subroutine add_gaussian_noise_vector

real(rp) function gaussian_rand()
  implicit none
  real(rp) :: u1, u2
  call random_number(u1)
  call random_number(u2)
  u1 = max(u1, 1.0e-12_rp)
  gaussian_rand = sqrt(-2.0_rp*log(u1))*cos(2.0_rp*acos(-1.0_rp)*u2)
end function gaussian_rand

pure function cross3(v1, v2) result(v3)
  implicit none
  real(rp), intent(in) :: v1(3), v2(3)
  real(rp) :: v3(3)
  v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
  v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
  v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
end function cross3

real(rp) function det3(amat)
  implicit none
  real(rp), intent(in) :: amat(3, 3)
  det3 = amat(1,1)*(amat(2,2)*amat(3,3)-amat(2,3)*amat(3,2)) - &
         amat(1,2)*(amat(2,1)*amat(3,3)-amat(2,3)*amat(3,1)) + &
         amat(1,3)*(amat(2,1)*amat(3,2)-amat(2,2)*amat(3,1))
end function det3


end module rsutil_lib
