module rsutil_c_interface
  ! This module provides the C interface for process_lattice.
  ! It is NOT wrapped by f90wrap, only called directly from ctypes.
  use iso_c_binding
  use rsutil_config
  use rsutil_lib, only: process_lattice_impl
  implicit none

contains

subroutine process_lattice(cluster_type_str, empty_spheres_flag, min_radius, max_radius, lmax_flag) bind(c,name='process_lattice')
   implicit none
   character(kind=c_char), intent(in) :: cluster_type_str(*)
   integer(c_int), intent(in), value :: empty_spheres_flag
   real(c_double), intent(in), value :: min_radius
   real(c_double), intent(in), value :: max_radius
   integer(c_int), intent(in), value :: lmax_flag
   character(len=20) :: cluster_type
   integer :: i, slen

   ! Convert C string to Fortran string
   i = 1
   do while (cluster_type_str(i) /= c_null_char .and. i < 20)
      cluster_type(i:i) = cluster_type_str(i)
      i = i + 1
   end do
   slen = i - 1
   cluster_type = cluster_type(1:slen)

   ! Store config variables
   enable_empty_spheres = (empty_spheres_flag == 1)
   min_empty_sphere_radius = real(min_radius)
   max_empty_sphere_radius = real(max_radius)
   lmax_override = int(lmax_flag)

   ! Print settings
   if (trim(cluster_type) == "free") then
      print *, "Setting cluster type to FREE (pbc = .false.)"
   else
      print *, "Setting cluster type to PERIODIC (pbc = .true.)"
   end if
   
   if (enable_empty_spheres) then
      print *, "Empty spheres ENABLED"
      print *, "  Min radius = ", min_empty_sphere_radius
      print *, "  Max radius = ", max_empty_sphere_radius
   else
      print *, "Empty spheres DISABLED"
   end if

   ! Call the main processing routine
   call process_lattice_impl()

end subroutine process_lattice

end module rsutil_c_interface
