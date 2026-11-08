module rsutil_config
  ! This module holds global configuration flags for rsutil
  ! It's small and has no dependencies, avoiding circular imports
  use system_data
  implicit none

  character(len=20) :: cluster_type = "periodic"
  logical :: enable_empty_spheres = .false.
  character(len=16) :: empty_sphere_backend = "xband"
  real(rp) :: min_empty_sphere_radius = 1.2_rp
  real(rp) :: max_empty_sphere_radius = 2.5_rp
  integer :: lmax_override = -1
  real(rp) :: lattice_dist_sigma = 0.0_rp
  real(rp) :: basis_dist_sigma = 0.0_rp
  real(rp) :: alat_scale_factor = 1.0_rp
  integer :: random_seed_value = 0
  logical :: use_fixed_seed = .false.
  logical :: use_2d_mode = .false.
  real(rp) :: z_iso_factor = -1.0_rp
  integer :: rec_steps = 20
  character(len=16) :: rec_method = "block"

end module rsutil_config
