module lattice_io
  use system_data
  implicit none
  public


    ! Declare integer variables
    integer :: nbulk_bulk
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

    real(rp), allocatable :: a(:, :)   ! Lattice vectors
    real(rp), allocatable :: crd(:, :) ! Basis vectors
    real(rp), allocatable :: ct(:)
    real(rp), allocatable :: ws(:)
    real(rp) :: r2

    ! Namelist declaration
    namelist /lattice/ nbulk_bulk, ntot, ntype, nbas, nrec, a, crd, &
            izp, no, iu, ib, irec, ct, r2

  interface write_array
    module procedure write_array_real
    module procedure write_array_int
  end interface write_array

    interface expand_array
      module procedure expand_array_int
      module procedure expand_array_real
    end interface expand_array

contains

  ! Routine to read the lattice namelist
  subroutine read_lattice_namelist(filename, lattice_struct)
    implicit none
    character(len=*), intent(in) :: filename
    type(lattice_t), intent(out) :: lattice_struct

    ! Local variables
    integer :: read_status


    ! Attempt to read namelist without allocating arrays
    open(unit=10, file=filename, status='old', action='read')
    read(10, nml=lattice, iostat=read_status)

    ! Check the read status
    if (read_status /= 0) then
       ! Partial read: Use available variables to determine dimensions
       ! print *, "Determining array sizes from available data..."
       if (nrec > 0) then
          allocate(a(3, 3), crd(3, nbulk_bulk))
          allocate(izp(nrec), no(nrec), irec(nrec))
          allocate (ib(ntot), iu(ntot), ct(ntype))
          allocate (ws(ntot))
          ! Sett WS to two since it is not read
          ws = 2.0d0
       else
          print *, "Error: 'nrec' not defined in the input file."
          stop
       end if
    end if

    ! Re-read the namelist after allocating arrays
    rewind(10)  ! Rewind the file for a fresh read
    read(10, nml=lattice, iostat=read_status)

    ! Check for errors during the second read
    if (read_status /= 0) then
       print *, "Error reading the namelist even after allocation. Status:", read_status
       stop
    end if

    ! Close the file
    close(10)

    print *, "Namelist read successfully."

    call copy_to_lattice_t(lattice_struct)

  end subroutine read_lattice_namelist

  function itoa(i)
    implicit none
    integer, intent(in) :: i
    character(len=10) :: itoa
    write(itoa, '(I0)') i
  end function itoa


  subroutine write_lattice_namelist(filename, lattice_struct)
    implicit none
    character(len=*), intent(in) :: filename
    type(lattice_t), intent(in) :: lattice_struct

    ! Local variables
    integer :: i, j
  
    ! Copy the lattice data to the module variables
    call copy_from_lattice_t(lattice_struct)

    ! Open the file for writing
    open(unit=20, file=filename, status='replace', action='write')
  
    ! Begin writing the namelist
    write(20, '(A)') "&lattice"
    write(20, '(A, I4)') "nbulk_bulk = ", nbulk_bulk
    write(20, '(A, I4)') "ntot = ", ntot
    write(20, '(A, I4)') "ntype = ", ntype
    write(20, '(A, I4)') "nbas = ", nbas
    write(20, '(A, I4)') "nrec = ", nrec
  
    ! Write the matrix `a`
    call write_array_2d(20, a, 3, 3, 'a')
  
    ! Write the matrix `crd`
    call write_array_2d(20, crd, 3, nrec, 'crd')
  
    ! Write the array `izp`
    call write_array(20, izp, nrec, "izp")
  
    ! Write the array `no`
    call write_array(20, no, nrec, "no")
  
    ! Write the array `iu`
    call write_array(20, iu, ntot, "iu")
  
    ! Write the array `ib`
    call write_array(20, ib, ntot, "ib")
  
    ! Write the array `irec`
    call write_array(20, irec, nrec, "irec")
  
    ! Write the array `ct`
    call write_array(20, ct, ntype, "ct")
  
    ! Write the scalar `r2`
    write(20, '(A, F10.4)') "r2 =", r2
  
    ! End the namelist
    write(20, '(A)') "/"
    ! End the namelist
    write(20, '(A)') "/"
  
    ! Close the file
    close(20)
  
    write(*,*) "Lattice namelist written to file:", filename
  
  end subroutine write_lattice_namelist



  subroutine write_array_real(unit, array, length, label)
    implicit none
    integer, intent(in) :: unit, length
    real*8, intent(in) :: array(:)
    character(len=*), intent(in) :: label
    integer :: i

    do i = 1, length
      write(unit, '(A, F8.2)') trim(label) // "(" // trim(adjustl(itoa(i))) // ") =", array(i)
    end do
    write(unit, '(A)') "!"
  end subroutine write_array_real

  subroutine write_array_int(unit, array, length, label)
    implicit none
    integer, intent(in) :: unit, length
    integer, intent(in) :: array(:)
    character(len=*), intent(in) :: label
    integer :: i

    do i = 1, length
      write(unit, '(A, I8)') trim(label) // "(" // trim(adjustl(itoa(i))) // ") =", array(i)
    end do
    write(unit, '(A)') "!"
  end subroutine write_array_int

  subroutine write_array_2d(unit, array, dim1, dim2, label)
    implicit none
    integer, intent(in) :: unit, dim1, dim2
    real*8, intent(in) :: array(:, :)
    character(len=*), intent(in) :: label
    integer :: i, j

    do j = 1, dim2
      write(unit, '(A, 2(F12.8,","), F12.8)') trim(label) // "(:, " // trim(adjustl(itoa(j))) // ") = ", array(:, j)
    end do
    write(unit, '(A)') "!"
  end subroutine write_array_2d

  subroutine copy_to_lattice_t(lattice_struct)
    implicit none
    type(lattice_t), intent(out) :: lattice_struct

    lattice_struct%nbulk_bulk = nbulk_bulk
    lattice_struct%ntot = ntot
    lattice_struct%ntype = ntype
    lattice_struct%nbas = nbas
    lattice_struct%nrec = nrec

    allocate(lattice_struct%a(size(a, 1), size(a, 2)))
    lattice_struct%a = a

    allocate(lattice_struct%crd(size(crd, 1), size(crd, 2)))
    lattice_struct%crd = crd

    allocate(lattice_struct%izp(size(izp)))
    lattice_struct%izp = izp

    allocate(lattice_struct%ws(size(ws)))
    lattice_struct%ws = ws

    allocate(lattice_struct%no(size(no)))
    lattice_struct%no = no

    allocate(lattice_struct%iu(size(iu)))
    lattice_struct%iu = iu

    allocate(lattice_struct%ib(size(ib)))
    lattice_struct%ib = ib

    allocate(lattice_struct%irec(size(irec)))
    lattice_struct%irec = irec

    allocate(lattice_struct%ct(size(ct)))
    lattice_struct%ct = ct

    lattice_struct%r2 = r2
  end subroutine copy_to_lattice_t

  subroutine copy_from_lattice_t(lattice_struct)
    implicit none
    type(lattice_t), intent(in) :: lattice_struct

    nbulk_bulk = lattice_struct%nbulk_bulk
    ntot = lattice_struct%ntot
    ntype = lattice_struct%ntype
    nbas = lattice_struct%nbas
    nrec = lattice_struct%nrec

    if (allocated(a)) deallocate(a)
    allocate(a(size(lattice_struct%a, 1), size(lattice_struct%a, 2)))
    a = lattice_struct%a

    if (allocated(crd)) deallocate(crd)
    allocate(crd(size(lattice_struct%crd, 1), size(lattice_struct%crd, 2)))
    crd = lattice_struct%crd

    if (allocated(izp)) deallocate(izp)
    allocate(izp(size(lattice_struct%izp)))
    izp = lattice_struct%izp

    if (allocated(ws)) deallocate(ws)
    allocate(ws(size(lattice_struct%ws)))
    ws = lattice_struct%ws

    if (allocated(no)) deallocate(no)
    allocate(no(size(lattice_struct%no)))
    no = lattice_struct%no

    if (allocated(iu)) deallocate(iu)
    allocate(iu(size(lattice_struct%iu)))
    iu = lattice_struct%iu

    if (allocated(ib)) deallocate(ib)
    allocate(ib(size(lattice_struct%ib)))
    ib = lattice_struct%ib

    if (allocated(irec)) deallocate(irec)
    allocate(irec(size(lattice_struct%irec)))
    irec = lattice_struct%irec

    if (allocated(ct)) deallocate(ct)
    allocate(ct(size(lattice_struct%ct)))
    ct = lattice_struct%ct

    r2 = lattice_struct%r2
  end subroutine copy_from_lattice_t

    subroutine initialize_arrays_to_nbulk()
      implicit none
      integer :: i

      ! Set ntot, ntype, nbas, and nrec to nbulk_bulk
      ntot = nbulk_bulk
      ntype = nbulk_bulk
      nbas = nbulk_bulk
      nrec = nbulk_bulk

      ! Allocate arrays
      if (allocated(izp)) deallocate(izp)
      if (allocated(no)) deallocate(no)
      if (allocated(iu)) deallocate(iu)
      if (allocated(ib)) deallocate(ib)
      if (allocated(irec)) deallocate(irec)
      if (allocated(irec)) deallocate(ct)
      allocate(izp(nbulk_bulk), no(nbulk_bulk), iu(nbulk_bulk), ib(nbulk_bulk), irec(nbulk_bulk), ct(nbulk_bulk))

      ! Fill arrays with running numbers
      do i = 1, nbulk_bulk
        izp(i) = i
        no(i) = i
        iu(i) = i
        ib(i) = i
        irec(i) = i
        ct(i) = 4.0d0
      end do
    end subroutine initialize_arrays_to_nbulk

    ! Declare allocatable real arrays and variables
    subroutine expand_and_fill_arrays(lat_data)
      implicit none

      type(lattice_t), intent(in) :: lat_data

      integer :: i, old_size

      ! Store the old size
      old_size = size(izp)

      ! Set ntot, ntype, nbas, and nrec to nbulk_bulk
      nbulk_bulk =  lat_data%nbulk_bulk
      ntot =  lat_data%nbulk_bulk
      ntype = lat_data%nbulk_bulk
      nbas =  lat_data%nbulk_bulk
      nrec =  lat_data%nbulk_bulk

      ! Update lattice (redundant)
      if (allocated(a)) deallocate(a)
      allocate(a(size(lat_data%a, 1), size(lat_data%a, 2)))
      a = lat_data%a

      ! Update basis 
      if (allocated(crd)) deallocate(crd)
      allocate(crd(size(lat_data%crd, 1), size(lat_data%crd, 2)))
      crd = lat_data%crd

      ! Expand arrays to new dimensions
      call expand_array(izp, lat_data%nbulk_bulk)
      call expand_array(no, lat_data%nbulk_bulk)
      call expand_array(iu, lat_data%nbulk_bulk)
      call expand_array(ib, lat_data%nbulk_bulk)
      call expand_array(irec, lat_data%nbulk_bulk)
      call expand_array(ct, lat_data%nbulk_bulk)
      call expand_array(ws, lat_data%nbulk_bulk)

      ! Fill the new entries with running numbers
      do i = 1, lat_data%nbulk_bulk
        izp(i) = i
        no(i) = i
        iu(i) = i
        ib(i) = i
        irec(i) = i
        ! Copy ws even if not included in namelist
        ws(i) = lat_data%ws(i)
        ! Set a default ct scaling here. Actual scaling = 1.5*1.88 due to bohr vs Å
        ct(i) = sqrt(2.01d0) * lat_data%ws(i)
      end do

      ! Calculate r2 as max(ct)^2
      r2 = (1.01d0*maxval(ct))**2
    end subroutine expand_and_fill_arrays

    subroutine expand_array_int(array, new_size)
      implicit none
      integer, intent(inout), allocatable :: array(:)
      integer, intent(in) :: new_size
      integer, allocatable :: temp(:)

      ! Allocate temporary array and copy old values
      allocate(temp(size(array)))
      temp = array

      ! Deallocate old array and allocate new array with new size
      deallocate(array)
      allocate(array(new_size))

      ! Copy back old values to the new array
      array(1:size(temp)) = temp

      ! Deallocate temporary array
      deallocate(temp)
    end subroutine expand_array_int

    subroutine expand_array_real(array, new_size)
      implicit none
      real(rp), intent(inout), allocatable :: array(:)
      integer, intent(in) :: new_size
      real(rp), allocatable :: temp(:)

      ! Allocate temporary array and copy old values
      allocate(temp(size(array)))
      temp = array

      ! Deallocate old array and allocate new array with new size
      deallocate(array)
      allocate(array(new_size))

      ! Copy back old values to the new array
      array(1:size(temp)) = temp

      ! Deallocate temporary array
      deallocate(temp)
    end subroutine expand_array_real

end module lattice_io