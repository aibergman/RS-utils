module rsutil_config
  ! This module holds global configuration flags for rsutil
  ! It's small and has no dependencies, avoiding circular imports
  use system_data
  implicit none

  character(len=20) :: cluster_type = "periodic"
  logical :: enable_empty_spheres = .false.
  real(rp) :: min_empty_sphere_radius = 1.2_rp
  real(rp) :: max_empty_sphere_radius = 2.5_rp
  integer :: lmax_override = -1

end module rsutil_config
