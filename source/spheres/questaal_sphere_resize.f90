! Standalone adaptation of Questaal's sclwsr/sclws2 radius policy.
! The geometric neighbor list is built locally instead of calling pairc;
! the iterative locking and constrained-volume scaling follows sclws2.
module questaal_sphere_resize
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: resize_questaal_spheres

  real(real64), parameter :: pi43 = 4.1887902047863909846_real64
  real(real64), parameter :: tiny = 1.0e-5_real64

contains

  subroutine resize_questaal_spheres(lattice, alat, basis, species, nuclear_charge, &
                                      radii, lock, target_fraction, max_radius, &
                                      max_overlap_distance, max_overlap_radius, actual_fraction)
    real(real64), intent(in) :: lattice(3,3), alat, basis(:,:)
    integer, intent(in) :: species(:)
    real(real64), intent(in) :: nuclear_charge(:), max_radius
    real(real64), intent(in) :: max_overlap_distance(3), max_overlap_radius(3)
    real(real64), intent(inout) :: radii(:)
    integer, intent(inout) :: lock(:)
    real(real64), intent(in) :: target_fraction
    real(real64), intent(out) :: actual_fraction

    integer, allocatable :: n_per_species(:), lock_local(:), lock_saved(:)
    real(real64), allocatable :: radii_initial(:), neighbor_distance(:)
    integer, allocatable :: pair_i(:), pair_j(:), pair_shift(:,:)
    integer :: nbas, nspec, i, is, n_pairs
    real(real64) :: cell_volume, empty_target, ratio, neighbor_cutoff
    logical :: has_empty

    nbas = size(basis,2)
    nspec = size(radii)
    if (size(species) /= nbas .or. size(lock) /= nspec .or. &
        size(nuclear_charge) /= nspec .or. size(basis,1) /= 3 .or. &
        size(lattice,1) /= 3 .or. size(lattice,2) /= 3) error stop 'invalid sphere-resize dimensions'
    if (any(species < 1) .or. any(species > nspec)) error stop 'invalid sphere species index'
    if (alat <= 0.0_real64 .or. target_fraction <= 0.0_real64) error stop 'invalid sphere-resize scale'

    allocate(n_per_species(nspec), radii_initial(nspec), lock_local(nspec), lock_saved(nspec))
    n_per_species = 0
    do i = 1, nbas
      n_per_species(species(i)) = n_per_species(species(i)) + 1
    end do
    radii_initial = radii
    lock_local = lock
    cell_volume = abs(det3(lattice))*alat**3
    if (cell_volume <= tiny) error stop 'singular cell in sphere resizing'

    ! Include at least the nearest periodic images even when initial Wigner-
    ! Seitz radii are much smaller than half the lattice spacing. Those image
    ! pairs cap radius growth in sclws2; omitting them lets a sphere expand
    ! through its own periodic copies before FINDES tries adding interstitials.
    neighbor_cutoff = max(2.0_real64*maxval(radii), &
                          alat*sum(sqrt(sum(lattice*lattice,dim=1))))
    call build_neighbors(lattice, alat, basis, neighbor_cutoff, &
                         pair_i, pair_j, pair_shift, neighbor_distance)
    n_pairs = size(pair_i)
    has_empty = any(nuclear_charge == 0.0_real64)

    ! Questaal calls sclwsr with opt=30: first resize occupied spheres while
    ! empty-sphere radii are zero, then resize empty spheres with atom radii
    ! held fixed.
    if (has_empty) then
      lock_saved = lock_local
      do is = 1, nspec
        if (nuclear_charge(is) == 0.0_real64) then
          lock_local(is) = 2
          radii(is) = 0.0_real64
        end if
      end do
      empty_target = target_fraction
      call sclws2_local(basis, species, nuclear_charge, n_per_species, &
                        pair_i, pair_j, pair_shift, neighbor_distance, max_radius, &
                        max_overlap_distance, max_overlap_radius, lock_local, &
                        cell_volume, empty_target, radii)
      lock_local = lock_saved
      do is = 1, nspec
        if (nuclear_charge(is) == 0.0_real64) radii(is) = radii_initial(is)
      end do
      do is = 1, nspec
        if (nuclear_charge(is) /= 0.0_real64) lock_local(is) = 2
      end do
      ratio = sphere_volume(radii, n_per_species)/cell_volume
      if (ratio < target_fraction) then
        empty_target = target_fraction
        call sclws2_local(basis, species, nuclear_charge, n_per_species, &
                          pair_i, pair_j, pair_shift, neighbor_distance, max_radius, &
                          max_overlap_distance, max_overlap_radius, lock_local, &
                          cell_volume, empty_target, radii)
      end if
    else
      empty_target = target_fraction
      call sclws2_local(basis, species, nuclear_charge, n_per_species, &
                        pair_i, pair_j, pair_shift, neighbor_distance, max_radius, &
                        max_overlap_distance, max_overlap_radius, lock_local, &
                        cell_volume, empty_target, radii)
    end if

    actual_fraction = sphere_volume(radii, n_per_species)/cell_volume
    lock = lock_local
  end subroutine resize_questaal_spheres

  subroutine sclws2_local(basis, species, nuclear_charge, n_per_species, &
                           pair_i, pair_j, pair_shift, distance, max_radius, omax1, omax2, &
                           lock, cell_volume, target_fraction, radii)
    real(real64), intent(in) :: basis(:,:), nuclear_charge(:)
    integer, intent(in) :: species(:), n_per_species(:), pair_i(:), pair_j(:), pair_shift(:,:)
    real(real64), intent(in) :: distance(:), max_radius, omax1(3), omax2(3), cell_volume
    integer, intent(inout) :: lock(:)
    real(real64), intent(inout) :: target_fraction, radii(:)

    real(real64) :: dovl1(3), dovl2(3), dscl(size(radii)), wsrs(size(radii))
    real(real64) :: avw, gamma, gw, amax, amax1, amax2, amax3, amax4
    real(real64) :: wsi, wsj, r, opo1, omo2, a, b, bmax, vola, volb
    real(real64) :: dm(0:3), vconst, rr, ss, tt, pp, qq, dd, uu, vv, x
    real(real64) :: ratio
    integer :: nspec, nloop, i, j, is, js, ip, q
    logical :: finished

    nspec = size(radii)
    gamma = 0.0_real64
    avw = (cell_volume/real(max(1,size(basis,2)),real64)*3.0_real64/(4.0_real64*pi43))**(1.0_real64/3.0_real64)
    gw = avw*gamma
    wsrs = radii

    do nloop = 1, nspec+1
      amax = 9.0e9_real64
      do is = 1, nspec
        if (lock(is) /= 2) lock(is) = 0
      end do

      ! Lock species which already meet an overlap or global-radius limit.
      do q = 1, size(pair_i)
        i = pair_i(q); j = pair_j(q)
        is = species(i); js = species(j)
        if (is == js .and. i == j .and. all(pair_shift(:,q) == 0)) cycle
        ip = overlap_class(nuclear_charge(is), nuclear_charge(js))
        wsi = radii(is); wsj = radii(js); r = distance(q)
        opo1 = 1.0_real64 + omax1(ip)
        omo2 = 1.0_real64 - omax2(ip)
        if (abs(opo1*r-wsi-wsj) < tiny .or. abs(r-omo2*wsi-wsj) < tiny .or. &
            abs(r-wsi-omo2*wsj) < tiny) lock(is) = max(lock(is),1)
        if (max_radius > 0.0_real64 .and. abs(wsi-max_radius) < tiny) lock(is)=max(lock(is),1)
      end do

      ! Release a tentative lock if any constraint is exceeded.
      do q = 1, size(pair_i)
        i = pair_i(q); j = pair_j(q)
        is = species(i); js = species(j)
        if (lock(is) == 2) cycle
        if (is == js .and. i == j .and. all(pair_shift(:,q) == 0)) cycle
        ip = overlap_class(nuclear_charge(is), nuclear_charge(js))
        wsi = radii(is); wsj = radii(js); r = distance(q)
        opo1 = 1.0_real64 + omax1(ip)
        omo2 = 1.0_real64 - omax2(ip)
        if (opo1*r-wsi-wsj < -tiny .or. r-omo2*wsi-wsj < -tiny .or. &
            r-wsi-omo2*wsj < -tiny) lock(is)=0
        if (max_radius > 0.0_real64 .and. wsi-max_radius > tiny) lock(is)=0
      end do

      ! Largest allowed scaling factor for all currently floating species.
      do is = 1, nspec
        if (lock(is) /= 0) cycle
        do q = 1, size(pair_i)
          i = pair_i(q); j = pair_j(q)
          if (species(i) /= is) cycle
          js = species(j)
          if (is == js .and. i == j .and. all(pair_shift(:,q) == 0)) cycle
          ip = overlap_class(nuclear_charge(is), nuclear_charge(js))
          r = distance(q); wsi = radii(is); wsj = radii(js)
          opo1 = 1.0_real64 + omax1(ip)
          omo2 = 1.0_real64 - omax2(ip)
          amax1 = huge(1.0_real64); amax2 = amax1; amax3 = amax1; amax4 = amax1
          if (lock(js) /= 0) then
            amax1 = (opo1*r-wsj+gw)/(wsi+gw)
            if (omo2 > 0.0_real64) amax2 = (r-wsj+omo2*gw)/(omo2*wsi+omo2*gw)
            amax3 = (r-omo2*wsj+gw)/(wsi+gw)
          else
            amax1 = (opo1*r+2.0_real64*gw)/(wsi+wsj+2.0_real64*gw)
            if (wsj+omo2*wsi+(1.0_real64+omo2)*gw > 0.0_real64) &
              amax2 = (r+(1.0_real64+omo2)*gw)/(wsj+omo2*wsi+(1.0_real64+omo2)*gw)
            if (wsi+omo2*wsj+(1.0_real64+omo2)*gw > 0.0_real64) &
              amax3 = (r+(1.0_real64+omo2)*gw)/(wsi+omo2*wsj+(1.0_real64+omo2)*gw)
          end if
          if (max_radius > 0.0_real64 .and. wsi > 0.0_real64) amax4=max_radius/wsi
          amax = min(amax,amax1,amax2,amax3,amax4)
        end do
      end do

      if (amax == 9.0e9_real64 .or. amax == huge(1.0_real64)) exit
      bmax = (1.0_real64-1.0_real64/amax)*gw
      vola = 0.0_real64; volb = 0.0_real64
      do i = 1, size(species)
        is = species(i)
        if (lock(is) == 0) then
          volb = volb + (amax*radii(is)+bmax)**3
        else
          vola = vola + radii(is)**3
        end if
      end do
      vola = vola*pi43; volb = volb*pi43

      if (cell_volume*target_fraction < vola+volb) then
        dm = 0.0_real64
        do i = 1, size(species)
          is = species(i)
          if (lock(is) /= 0) cycle
          a = radii(is)+gw
          if (abs(gamma) > 1.0_real64) then
            b = radii(is)
          else
            b = -gw
          end if
          dm(0)=dm(0)+b**3
          dm(1)=dm(1)+3.0_real64*a*b*b
          dm(2)=dm(2)+3.0_real64*a*a*b
          dm(3)=dm(3)+a**3
        end do
        vconst = vola-pi43*dm(3)
        if (vconst*(1.0_real64-tiny) > cell_volume*target_fraction) exit
        if (abs(dm(3)) > tiny) then
          rr=dm(2)/dm(3); ss=dm(1)/dm(3)
          tt=(dm(0)-(cell_volume*target_fraction-vola)/pi43)/dm(3)
          pp=ss-rr*rr/3.0_real64
          qq=2.0_real64*rr**3/27.0_real64-rr*ss/3.0_real64+tt
          dd=pp**3/27.0_real64+qq*qq/4.0_real64
          if (dd >= 0.0_real64) then
            uu=sign(abs(sqrt(dd)-qq/2.0_real64)**(1.0_real64/3.0_real64),sqrt(dd)-qq/2.0_real64)
            if (abs(uu) > tiny) then
              vv=-pp/uu/3.0_real64
              x=uu+vv-rr/3.0_real64
              if (abs(gamma) > 1.0_real64) then
                amax=x+1.0_real64; bmax=x*gw/amax
              else
                amax=x; bmax=(1.0_real64-1.0_real64/amax)*gw
              end if
            end if
          end if
        end if
      end if

      finished = .true.
      do is = 1, nspec
        dscl(is)=1.0_real64
        if (lock(is) == 0) then
          if (radii(is) > tiny) then
            dscl(is)=amax+bmax/radii(is)
            radii(is)=dscl(is)*radii(is)
          else
            radii(is)=max(0.0_real64,amax*radii(is)+bmax)
          end if
          finished=.false.
        end if
      end do
      ratio=sphere_volume(radii,n_per_species)/cell_volume
      finished=finished .or. abs(ratio-target_fraction)<tiny .or. nloop==nspec+1
      if (finished) then
        dovl1=0.0_real64; dovl2=0.0_real64
        do q = 1, size(pair_i)
          i=pair_i(q); j=pair_j(q); is=species(i); js=species(j)
          if (is == js .and. i == j .and. all(pair_shift(:,q)==0)) cycle
          ip=overlap_class(nuclear_charge(is),nuclear_charge(js))
          wsi=radii(is); wsj=radii(js); r=distance(q)
          if (r /= 0.0_real64) dovl1(ip)=max(dovl1(ip),(wsi+wsj-r)/r)
          if (wsi /= 0.0_real64) dovl2(ip)=max(dovl2(ip),(wsi+wsj-r)/wsi)
          if (wsj /= 0.0_real64) dovl2(ip)=max(dovl2(ip),(wsi+wsj-r)/wsj)
        end do
        target_fraction=ratio
        return
      end if
    end do
  end subroutine sclws2_local

  subroutine build_neighbors(lattice, alat, basis, cutoff, pair_i, pair_j, pair_shift, distance)
    real(real64), intent(in) :: lattice(3,3), alat, basis(:,:), cutoff
    integer, allocatable, intent(out) :: pair_i(:), pair_j(:), pair_shift(:,:)
    real(real64), allocatable, intent(out) :: distance(:)
    integer :: nmax(3), i, j, a, b, c, n, cap, k
    integer, allocatable :: pi(:), pj(:), ps(:,:)
    real(real64), allocatable :: pd(:)
    real(real64) :: qlat(3,3), det, delta(3), r

    call inverse3(lattice, qlat, det)
    if (abs(det) < tiny) error stop 'singular lattice in neighbor generation'
    do k=1,3
      nmax(k)=ceiling(max(0.0_real64,cutoff/alat*sqrt(sum(qlat(:,k)**2))))+1
    end do
    cap=max(128,size(basis,2)**2*product(2*nmax+1))
    allocate(pi(cap),pj(cap),ps(3,cap),pd(cap))
    n=0
    do i=1,size(basis,2)
      do j=1,size(basis,2)
        do a=-nmax(1),nmax(1)
          do b=-nmax(2),nmax(2)
            do c=-nmax(3),nmax(3)
              if (i==j .and. a==0 .and. b==0 .and. c==0) cycle
              delta=basis(:,i)-basis(:,j)-matmul(lattice,[real(a,real64),real(b,real64),real(c,real64)])
              r=sqrt(sum(delta*delta))*alat
              if (r>cutoff+tiny) cycle
              n=n+1
              if (n>cap) error stop 'neighbor list capacity estimate failed'
              pi(n)=i; pj(n)=j; ps(:,n)=[a,b,c]; pd(n)=r
            end do
          end do
        end do
      end do
    end do
    allocate(pair_i(n),pair_j(n),pair_shift(3,n),distance(n))
    pair_i=pi(:n); pair_j=pj(:n); pair_shift=ps(:,:n); distance=pd(:n)
  end subroutine build_neighbors

  integer function overlap_class(zi,zj) result(k)
    real(real64), intent(in) :: zi,zj
    k=2
    if (nint(zi)/=0 .and. nint(zj)/=0) k=1
    if (nint(zi)==0 .and. nint(zj)==0) k=3
  end function overlap_class

  real(real64) function sphere_volume(radii, counts) result(volume)
    real(real64), intent(in) :: radii(:)
    integer, intent(in) :: counts(:)
    volume=pi43*sum(real(counts,real64)*radii**3)
  end function sphere_volume

  real(real64) function det3(a) result(det)
    real(real64), intent(in) :: a(3,3)
    det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
       -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
       +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  end function det3

  subroutine inverse3(a, inverse, determinant)
    real(real64), intent(in) :: a(3,3)
    real(real64), intent(out) :: inverse(3,3), determinant
    determinant=det3(a)
    if (abs(determinant)<tiny) then
      inverse=0.0_real64
      return
    end if
    inverse(1,1)= (a(2,2)*a(3,3)-a(2,3)*a(3,2))/determinant
    inverse(1,2)= (a(1,3)*a(3,2)-a(1,2)*a(3,3))/determinant
    inverse(1,3)= (a(1,2)*a(2,3)-a(1,3)*a(2,2))/determinant
    inverse(2,1)= (a(2,3)*a(3,1)-a(2,1)*a(3,3))/determinant
    inverse(2,2)= (a(1,1)*a(3,3)-a(1,3)*a(3,1))/determinant
    inverse(2,3)= (a(1,3)*a(2,1)-a(1,1)*a(2,3))/determinant
    inverse(3,1)= (a(2,1)*a(3,2)-a(2,2)*a(3,1))/determinant
    inverse(3,2)= (a(1,2)*a(3,1)-a(1,1)*a(3,2))/determinant
    inverse(3,3)= (a(1,1)*a(2,2)-a(1,2)*a(2,1))/determinant
  end subroutine inverse3
end module questaal_sphere_resize
