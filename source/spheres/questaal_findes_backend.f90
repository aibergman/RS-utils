! RS-utils adapter for the pinned, standalone Questaal FINDES source.
module questaal_findes_backend
  use iso_fortran_env, only: real64
  use structures, only: str_ctrl, str_lat, str_site
  implicit none
  private
  public :: add_questaal_empty_spheres

  interface
    subroutine findes(mode,sincl,nescut,mxesspec,s_ctrl,s_lat,s_site,alat,nbas,nclass,nrclas,mxclas, &
                      nsymop,symopm,symopv,ishorec,ishores,lock,wsr,z)
      use structures, only: str_ctrl, str_lat, str_site
      implicit none
      integer :: mode,nescut,mxesspec,nbas,nclass,nrclas(*),mxclas,nsymop,ishorec(3),ishores(3),lock(*)
      real(8) :: alat,symopm(9,*),symopv(3,*),wsr(*),z(*)
      character(*) :: sincl(2)
      type(str_ctrl) :: s_ctrl
      type(str_lat) :: s_lat
      type(str_site) :: s_site(*)
    end subroutine findes
  end interface

contains

  subroutine add_questaal_empty_spheres(lattice, alat, basis, atomic_number, radii, min_radius, max_radius, &
                                         basis_out, species_out, radii_out, status)
    real(real64), intent(in) :: lattice(3,3), alat, basis(:,:), atomic_number(:), radii(:)
    real(real64), intent(in) :: min_radius, max_radius
    real(real64), allocatable, intent(out) :: basis_out(:,:), radii_out(:)
    integer, allocatable, intent(out) :: species_out(:)
    integer, intent(out) :: status
    type(str_ctrl) :: ctrl
    type(str_lat) :: lat
    type(str_site) :: site(1)
    integer, allocatable :: nrclas(:), lock(:)
    real(real64), allocatable :: wsr(:), z(:)
    real(real64) :: symopm(9,1), symopv(3,1)
    character(8), allocatable :: labels(:)
    character(8) :: sincl(2)
    integer :: nbas, nclass, mxclas, i

    status = 1
    nbas = size(basis,2)
    if (size(basis,1) /= 3 .or. size(atomic_number) /= nbas .or. size(radii) /= nbas .or. &
        nbas < 1 .or. nbas > 200 .or. alat <= 0.0_real64) return
    ! RS-utils currently presents each site as its own sphere class. Keep that
    ! mapping; FINDES appends new classes and positions to these arrays.
    nclass = nbas
    mxclas = min(200, max(nbas + 64, 2*nbas))
    allocate(ctrl%ipc(nbas), ctrl%ips(nbas), ctrl%clabl(mxclas), ctrl%spid(nbas))
    allocate(lat%pos(3,nbas), nrclas(mxclas), lock(mxclas), wsr(mxclas), z(mxclas), labels(mxclas))
    ctrl%nesabc = [24,24,24]
    ctrl%nspec = nbas
    ctrl%rmines = min_radius
    ctrl%rmaxes = max_radius
    ctrl%sclwsr = 1.0_real64
    ! Match the overlapping-sphere packing used by the XBand finder. With
    ! both limits at zero, FINDES shrinks the atomic spheres to tangency and
    ! can leave a large part of the cell unfilled even after adding ES.
    ctrl%omax1 = 0.14_real64
    ctrl%omax2 = 0.26_real64
    ctrl%wsrmax = 0.0_real64
    lat%plat = lattice
    lat%plat0 = lattice
    lat%pos = basis
    ctrl%ipc = [(i,i=1,nbas)]
    ctrl%ips = ctrl%ipc
    nrclas = 0; nrclas(1:nbas) = 1
    lock = 0
    wsr = 0.0_real64; wsr(1:nbas) = radii
    z = 0.0_real64; z(1:nbas) = atomic_number
    do i=1,mxclas
      write(labels(i),'("S",i3.3)') i
    end do
    ctrl%clabl = labels
    ctrl%spid(1:nbas) = labels(1:nbas)
    symopm = 0.0_real64
    symopm(1,1)=1.0_real64; symopm(5,1)=1.0_real64; symopm(9,1)=1.0_real64
    symopv = 0.0_real64
    sincl = ' '
    call findes(1,sincl,mxclas-nbas,mxclas,ctrl,lat,site,alat,nbas,nclass,nrclas,mxclas,1, &
                symopm,symopv,[0,0,0],[0,0,0],lock,wsr,z)
    if (nbas > size(lat%pos,2) .or. nbas > size(ctrl%ipc)) return
    allocate(basis_out(3,nbas), species_out(nbas), radii_out(nbas))
    basis_out = lat%pos(:,1:nbas)
    species_out = ctrl%ipc(1:nbas)
    radii_out = wsr(species_out)
    status = 0
  end subroutine add_questaal_empty_spheres

end module questaal_findes_backend
