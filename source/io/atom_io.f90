!-------------------------------------------------------------------
! Module defining the derived types for the namelist blocks,
! the namelist declarations, and a routine to read the data.
!-------------------------------------------------------------------
module atom_io
   use system_data
    implicit none

    ! Screening parameters
    !> Original screening (From Jepsen)
    real(rp), dimension(4), parameter :: qm_canonical = [.348485d0, .053030d0, .010714d0, .002000d0]
    !> Screening from LMTO47 (RWS = average WS for all)
    real(rp), dimension(4), parameter :: qm_lmto47 = [0.3500000d0, 0.0571667d0, 0.0168070d0, 0.002000d0]
    ! !----------------------------------------------------------------
    ! ! Derived type for the &element block
    ! !----------------------------------------------------------------
    ! type :: element_t
    !    character(len=8) :: symbol
    !    integer :: atomic_number
    !    integer :: core
    !    integer :: valence
    !    integer :: f_core
    !    integer :: num_quant_s
    !    integer :: num_quant_p
    !    integer :: num_quant_d
    ! end type element_t
  
    ! type :: potential_t
    !     integer :: lmax
    !     real(rp) :: pl(3,2)
    !     ! Array dimensions as in your example: for ql, the dimensions are (3,3,2)
    !     real(rp) :: ql(3,3,2)
    !     real(rp) :: center_band(3,2)
    !     real(rp) :: shifted_band(3,2)
    !     real(rp) :: width_band(3,2)
    !     real(rp) :: ws_r
    !     real(rp) :: sumec
    !     real(rp) :: sumev
    !     real(rp) :: etot
    !     real(rp) :: utot
    !     real(rp) :: ekin
    !     real(rp) :: rhoeps
    !     real(rp) :: vmad
    !     ! Orthogonal potential parameters
    !     real(rp) :: enu(3,2)
    !     real(rp) :: c(3,2)
    !     real(rp) :: srdel(3,2)
    !     real(rp) :: qpar(3,2)
    !     real(rp) :: ppar(3,2)
    !     real(rp) :: obar(3,2)
    !     real(rp) :: vl(3,2)
    ! end type potential_t

contains

    !----------------------------------------------------------------
    ! Subroutine to read the namelist file and assign its contents
    ! to the derived type variables.
    !----------------------------------------------------------------
    subroutine read_atom(filename, element_data, par_data)
      character(len=*), intent(in) :: filename
      type(element_t), intent(out) :: element_data
      type(potential_t), intent(out) :: par_data
      integer :: unit, ioerr

      ! Namelist declarations for the &element and &par groups.
      ! integer, parameter :: rp = selected_real_kind(15, 307)

      ! Including namelists from RS-LMTO
      include "../include/potential.f90"
      include "../include/element.f90"

      ! Open the file for reading.
      open(newunit=unit, file=filename, status='old', action='read', iostat=ioerr)
      if (ioerr /= 0) then
         print *, "Error opening file: ", trim(filename)
         stop
      end if

      ! Read the namelist groups. The names in the file must match the
      ! names declared in the namelist statements above.
      read(unit, nml=element)
      read(unit, nml=par)
      close(unit)

      !----------------------------------------------------------------
      ! Assign the read values to the derived type components.
      !----------------------------------------------------------------
      element_data%symbol       = symbol
      element_data%atomic_number  = atomic_number
      element_data%core         = core
      element_data%valence      = valence
      element_data%f_core       = f_core
      element_data%num_quant_s  = num_quant_s
      element_data%num_quant_p  = num_quant_p
      element_data%num_quant_d  = num_quant_d

      par_data%lmax = lmax
      if (par_data%lmax < 0 .or. par_data%lmax > 3) then
         print *, "Error: read_atom found unsupported lmax=", par_data%lmax, " (supported 0..3)"
         stop 1
      end if
      call allocate_potential_arrays(par_data, par_data%lmax)
      par_data%pl           = pl(1:par_data%lmax+1, :)
      par_data%ql           = ql(:, 1:par_data%lmax+1, :)
      par_data%center_band  = center_band(1:par_data%lmax+1, :)
      par_data%shifted_band = shifted_band(1:par_data%lmax+1, :)
      par_data%width_band   = width_band(1:par_data%lmax+1, :)
      par_data%ws_r         = ws_r
      par_data%sumec        = sumec
      par_data%sumev        = sumev
      par_data%etot         = etot
      par_data%utot         = utot
      par_data%ekin         = ekin
      par_data%rhoeps       = rhoeps
      par_data%vmad         = vmad
      par_data%enu          = enu(1:par_data%lmax+1, :)
      par_data%c            = c(1:par_data%lmax+1, :)
      par_data%obar         = obar(1:par_data%lmax+1, :)
      par_data%srdel        = srdel(1:par_data%lmax+1, :)
      par_data%qpar         = qpar(1:par_data%lmax+1, :)
      par_data%ppar         = ppar(1:par_data%lmax+1, :)
      par_data%vl           = vl(1:par_data%lmax+1, :)


    end subroutine read_atom

    subroutine write_atom(filename, element, par)
      character(len=*), intent(in) :: filename
      type(element_t), intent(in) :: element
      type(potential_t), intent(in) :: par
      integer :: unit, i, j, k, ioerr, nl
  
      ! Open the output file with a new unit.
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ioerr)
      if (ioerr /= 0) then
         print *, "Error opening file: ", filename
         stop
      end if
      nl = par%lmax + 1
  
      !---------------------------
      ! Write the &element namelist
      !---------------------------
      write(unit,*) "&element"
      write(unit, '(4X,A,A,A)') "symbol='", trim(element%symbol) // "',"
      write(unit, '(4X,A,I0,A)') "atomic_number=", element%atomic_number, ","
      write(unit, '(4X,A,I0,A)') "core=", element%core, ","
      write(unit, '(4X,A,I0,A)') "valence=", element%valence, ","
      write(unit, '(4X,A,I0,A)') "f_core=", element%f_core, ","
      write(unit, '(4X,A,I0,A)') "num_quant_s=", element%num_quant_s, ","
      write(unit, '(4X,A,I0,A)') "num_quant_p=", element%num_quant_p, ","
      write(unit, '(4X,A,I0)') "num_quant_d=", element%num_quant_d
      write(unit,*) "/"
      write(unit,*)
  
      !---------------------------
      ! Write the &par namelist
      !---------------------------
      write(unit,*) "&par"
      ! Write ql arrays:
      write(unit, '(4X,A)') "! Charge (Q0)"
      do k = 1, 2
         write(unit, '(4X,A,I0,A)', advance='no') "ql(1, :, ", k, ") = "
         do i = 1, nl
            if (i < nl) then
               write(unit, '(F12.10,A)', advance='no') par%ql(1,i,k), ", "
            else
               write(unit, '(F12.10)') par%ql(1,i,k)
            end if
         end do
      end do
  
      write(unit, '(4X,A)') "! Second-order moments (Q2)"
      do k = 1, 2
         write(unit, '(4X,A,I0,A)', advance='no') "ql(3, :, ", k, ") = "
         do i = 1, nl
            if (i < nl) then
               write(unit, '(F12.10,A)', advance='no') par%ql(3,i,k), ", "
            else
               write(unit, '(F12.10)') par%ql(3,i,k)
            end if
         end do
      end do
  
      ! Write the pl array (dimensions 3x2)
      write(unit, '(4X,A)') "! Potential functions (PL)"
      call write_2d_array(unit, "pl", par%pl)
  
      ! Write center_band and width_band arrays.
      write(unit, '(4X,A)') "! TB potential parameters"
      call write_2d_array(unit, "center_band", par%center_band)
      call write_2d_array(unit, "width_band", par%width_band)
  
      ! Write scalar parameters.
      write(unit, '(4X,A)') "! Wigner-Seitz radius"
      write(unit, '(4X,A,F12.6)') "ws_r = ", par%ws_r

      write(unit, '(4X,A,I0)') "lmax = ", par%lmax
  
      ! write(unit, '(4X,A,F20.12)') "sumec= ", par%sumec
      ! write(unit, '(4X,A,F20.12)') "sumev= ", par%sumev
      ! write(unit, '(4X,A,F20.12)') "etot= ", par%etot
      ! write(unit, '(4X,A,F20.12)') "utot= ", par%utot
      ! write(unit, '(4X,A,F20.12)') "ekin= ", par%ekin
      ! write(unit, '(4X,A,F20.12)') "rhoeps= ", par%rhoeps
      ! write(unit, '(4X,A,F20.12)') "vmad= ", par%vmad
  
      ! Write Orthogonal potential parameters via helper subroutine.
      write(unit, '(4X,A)') "! ORT potential parameters"
      call write_2d_array(unit, "enu", par%enu)
      call write_2d_array(unit, "c",   par%c)
      call write_2d_array(unit, "srdel", par%srdel)
      call write_2d_array(unit, "qpar", par%qpar)
      call write_2d_array(unit, "ppar", par%ppar)
      call write_2d_array(unit, "vl", par%vl)
  
      write(unit,*) "/"
  
      close(unit)
    end subroutine write_atom
  
    !-------------------------------------------------------------------
    ! Helper subroutine to write a (nl)x2 real array.
    !-------------------------------------------------------------------
    subroutine write_2d_array(unit, name, arr)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real(kind=8), intent(in) :: arr(:, :)
      integer :: i, j
  
      do j = 1, size(arr, 2)
        write(unit, '(4X,A,A,I0,A)', advance='no') trim(name), "(:, ", j, ") = "
        do i = 1, size(arr, 1)
           if (i < size(arr, 1)) then
              write(unit, '(G20.14,A)', advance='no') arr(i,j), ", "
           else
              write(unit, '(G20.14)') arr(i,j)
           end if
        end do
     end do
    end subroutine write_2d_array



subroutine screen_ort(potential, wsm)
    ! Port of predls()
    implicit none

    type(potential_t), intent(inout) :: potential   !< Potential parameters and other info
    real(rp) :: wsm !< Average WS radius
    real(rp) :: x, y, wow
    real(rp), dimension(:, :), allocatable :: dele, qi, c, enu, del, q
    real(rp), dimension(4) :: qm
    integer :: i, j, lmax

    wow = wsm/potential%ws_r
    lmax = potential%lmax
    allocate(dele(potential%lmax + 1, 2), qi(potential%lmax + 1, 2))
    allocate(c(potential%lmax + 1, 2), enu(potential%lmax + 1, 2))
    allocate(del(potential%lmax + 1, 2), q(potential%lmax + 1, 2))
    enu(:, :) = potential%enu(:, :)
    c(:, :) = potential%c(:, :)
    qm(:) = qm_canonical(:)
    q(:, :) = potential%qpar(:, :)
    del(:, :) = potential%srdel(:, :)

    do J = 1, 2
       do I = 1, lmax + 1
        DELE(I, J) = DEL(I, J)*WOW**(0.5d0 - I)
        QI(I, J) = Q(I, J)*WOW**(1 - 2*I)
          X = 1 - ((QI(I, J) - QM(I))*(C(I, J) - ENU(I, J))/DELE(I, J)/DELE(I, J))
          Y = (QI(I, J) - QM(I))/ &
              (((C(I, J) - ENU(I, J))*(QI(I, J) - QM(I))) - (DELE(I, J)*DELE(I, J)))
          potential%center_band(i, j) = (C(I, J) - ENU(I, J))*X + ENU(I, J)
          potential%shifted_band(i, j) = (C(I, J) - ENU(I, J))*X
          potential%width_band(i, j) = DELE(I, J)*X
          potential%obar(i, j) = Y
       end do
    end do
 end subroutine screen_ort

 subroutine allocate_potential_arrays(par_data, lmax)
   type(potential_t), intent(inout) :: par_data
   integer, intent(in) :: lmax
   integer :: nl

   nl = lmax + 1
   if (allocated(par_data%pl)) deallocate(par_data%pl)
   if (allocated(par_data%ql)) deallocate(par_data%ql)
   if (allocated(par_data%center_band)) deallocate(par_data%center_band)
   if (allocated(par_data%shifted_band)) deallocate(par_data%shifted_band)
   if (allocated(par_data%width_band)) deallocate(par_data%width_band)
   if (allocated(par_data%enu)) deallocate(par_data%enu)
   if (allocated(par_data%c)) deallocate(par_data%c)
   if (allocated(par_data%srdel)) deallocate(par_data%srdel)
   if (allocated(par_data%qpar)) deallocate(par_data%qpar)
   if (allocated(par_data%ppar)) deallocate(par_data%ppar)
   if (allocated(par_data%obar)) deallocate(par_data%obar)
   if (allocated(par_data%vl)) deallocate(par_data%vl)

   allocate(par_data%pl(nl, 2))
   allocate(par_data%ql(3, nl, 2))
   allocate(par_data%center_band(nl, 2))
   allocate(par_data%shifted_band(nl, 2))
   allocate(par_data%width_band(nl, 2))
   allocate(par_data%enu(nl, 2))
   allocate(par_data%c(nl, 2))
   allocate(par_data%srdel(nl, 2))
   allocate(par_data%qpar(nl, 2))
   allocate(par_data%ppar(nl, 2))
   allocate(par_data%obar(nl, 2))
   allocate(par_data%vl(nl, 2))
 end subroutine allocate_potential_arrays

end module atom_io
