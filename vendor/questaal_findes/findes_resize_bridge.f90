! Compatibility entry points required by the revisioned FINDES source.
subroutine sclwsr(opt, nbas, nbasp, nspec, alat, plat, basp, ips, modep, &
                  slabl, z, lock, volfac, wsmax, omax1, omax2, wsr)
  use iso_fortran_env, only: real64
  use questaal_sphere_resize, only: resize_questaal_spheres
  implicit none
  integer, intent(in) :: opt, nbas, nbasp, nspec, modep(3), ips(nbas)
  integer, intent(inout) :: lock(nspec)
  real(real64), intent(in) :: alat, plat(3,3), basp(3,nbasp), z(nspec)
  real(real64), intent(in) :: wsmax, omax1(3), omax2(3)
  real(real64), intent(inout) :: volfac, wsr(nspec)
  character(8), intent(in) :: slabl(*)
  real(real64) :: actual, target

  if (opt /= 30) error stop 'standalone FINDES bridge supports opt=30 only'
  if (nbasp < nbas) error stop 'invalid padded basis size in FINDES bridge'
  target = volfac
  call resize_questaal_spheres(plat, alat, basp(:,1:nbas), ips, z, wsr, lock, &
                               target, wsmax, omax1, omax2, actual)
  ! FINDES compares the packing fraction returned by VOLSPh to its requested
  ! target after this call. Keep its target intact; the original application
  ! sclwsr also returns an achieved fraction, but that would make this caller
  ! treat every constrained resize as if it had met the requested target.
end subroutine sclwsr

double precision function volsph(nspec, nrspec, wsr)
  use iso_fortran_env, only: real64
  implicit none
  integer, intent(in) :: nspec, nrspec(nspec)
  real(real64), intent(in) :: wsr(nspec)
  volsph = 4.1887902047863909846_real64*sum(real(nrspec,real64)*wsr**3)
end function volsph
