module input_io
    use system_data
    implicit none

contains

    !-------------------------------------------------------------------
    ! Subroutine: read_namelists
    !
    ! Opens the file, reads the namelist groups &lattice and &atoms,
    ! and then assigns the values to the derived type variables.
    !-------------------------------------------------------------------
    subroutine read_input_namelists(filename, lat, atm, lat_lat)
        character(len=*), intent(in)  :: filename
        type(lattice_inp), intent(out) :: lat
        type(atoms_inp), intent(out) :: atm
        type(lattice_t), intent(in)   :: lat_lat

        ! Local variables
        integer :: unit, ioerr

        !-------------------------------------------------------------------
        ! Temporary variables for reading the lattice namelist.
        ! The names here must match the variable names in the input file.
        !-------------------------------------------------------------------
        integer :: rc, ndim=10000
        real(kind=8) :: alat
        character(len=8) :: crystal_sym
        real(kind=8) :: wav
        logical :: pbc=.true., b1=.true., b2=.true., b3=.true.
        integer :: n1=10, n2=10, n3=10, i, pos, j
        ! Additions from previous lattice.nml entries
        integer :: ntype = 1
        real(rp), allocatable :: ct(:)
        real(rp) :: r2
        character(len=16) :: strux_backend = ''
        logical :: strux_want_sdot = .false.
        real(rp) :: strux_solve_scale = 1.0_rp
        logical :: has_strux_backend, has_strux_want_sdot, has_strux_solve_scale
        character(len=1024) :: scan_line
        integer :: scan_ios

        namelist /lattice/ rc, alat, crystal_sym, wav, pbc, b1, b2, b3, n1, n2, n3, ndim, &
            ntype, ct, r2, strux_backend, strux_want_sdot, strux_solve_scale

        !-------------------------------------------------------------------
        ! Temporary variables for reading the atoms namelist.
        !-------------------------------------------------------------------
        character(len=256) :: database   ! directory of the database
        character(len=8), allocatable :: label(:)
        namelist /atoms/ database, label

        ! Local variables
        integer :: read_status

        ! Open the file for reading.
        open(newunit=unit, file=filename, status='old', action='read', iostat=ioerr)
        if (ioerr /= 0) then
            print *, "Error opening file: ", trim(filename)
            stop 1
        end if

        ! Set ntype to zero to track if it was read from the file.
        ntype = lat_lat%ntype
        r2 = lat_lat%r2
        has_strux_backend = .false.
        has_strux_want_sdot = .false.
        has_strux_solve_scale = .false.
        do
          read(unit, '(A)', iostat=scan_ios) scan_line
          if (scan_ios /= 0) exit
          scan_line = lower_string(scan_line)
          if (index(scan_line, "strux_backend") > 0) has_strux_backend = .true.
          if (index(scan_line, "strux_want_sdot") > 0) has_strux_want_sdot = .true.
          if (index(scan_line, "strux_solve_scale") > 0) has_strux_solve_scale = .true.
        end do
        rewind(unit)

        ! Read the &lattice namelist.
        allocate(ct(0))
        read(unit, nml=lattice,  iostat=read_status)
        ! If reading files, allocate ct array with ntype elements.
        if (read_status /= 0) then

           ! Allocate ct if needed
           deallocate(ct)
           allocate(ct(ntype))

           ! Re-read the namelist after allocating arrays
           rewind(unit)  ! Rewind the file for a fresh read
           read(unit, nml=lattice)
        end if

        ! Catch if ct was not defined in input.nml
            if (size(ct) /= ntype) then
                  deallocate(ct)
            allocate(ct(ntype))
            if(allocated(lat_lat%ct)) then
               ct = lat_lat%ct
            else
               ct = alat * 0.8_rp
            end if
        end if

        ! Check if r2 makes sense
         if (r2 == 0.0_rp) then
            r2 = (1.005d0*maxval(ct))**2
         end if

        ! Re-allocate atom label array
        allocate(atm%label(ntype))
        allocate(atm%atomic_number(ntype))
        allocate(label(ntype))

        ! Read the &atoms namelist.
        read(unit, nml=atoms)
        close(unit)

        ! Transfer the read values into the derived type instances.
        lat%rc = rc
        lat%alat = alat
        lat%crystal_sym = crystal_sym
        lat%wav = wav
        lat%pbc = pbc
        lat%b1 = b1
        lat%b2 = b2
        lat%b3 = b3
        lat%n1 = n1
        lat%n2 = n2
        lat%n3 = n3
        lat%ntype = ntype
        lat%ndim = ndim
        lat%r2 = r2
        lat%strux_backend = strux_backend
        lat%has_strux_backend = has_strux_backend
        lat%strux_want_sdot = strux_want_sdot
        lat%has_strux_want_sdot = has_strux_want_sdot
        lat%strux_solve_scale = strux_solve_scale
        lat%has_strux_solve_scale = has_strux_solve_scale
        allocate(lat%ct(lat%ntype))
        lat%ct = ct

        atm%database = database
         ! Strip any digits from the label
         do i = 1, ntype
             label(i) = adjustl(label(i))
             label(i) = trim(adjustl(label(i)))
             !      ! Loop to find the position where the number starts
             !  do j = len_trim(label(i)), 1, -1
             !     if (ichar(label(i)(j:j)) < ichar('0') .or. ichar(label(i)(j:j)) > ichar('9')) then
             !        pos = j
             !        exit
             !     end if
             !  end do
             ! if (pos > 0) then
             !    label(i) = label(i)(1:pos)  ! Take the substring before the digit
             ! end if
             print *,'Label:', label(i)
         end do
        atm%label = label
        atm%atomic_number = get_atomic_numbers(label, ntype)

    end subroutine read_input_namelists

      
        !-------------------------------------------------------------------
        ! Function: int_to_string
        ! Converts an integer to a character string.
        !-------------------------------------------------------------------
        function int_to_string(i) result(str)
          implicit none
          integer, intent(in) :: i
          character(len=16) :: str
          write(str, '(I0)') i
        end function int_to_string
      
        !-------------------------------------------------------------------
        ! Subroutine: write_lattice_block
        !
        ! Writes the updated &lattice namelist block to the output unit.
        !-------------------------------------------------------------------
        subroutine write_lattice_block(unit, lat)
          implicit none
          integer, intent(in) :: unit
          type(lattice_inp), intent(in) :: lat
          integer :: itype, ndim

          write(unit, '(A)') "&lattice"
          write(unit, '(A, I0)') "rc = ", lat%rc
          write(unit, '(A, F12.6)') "alat = ", lat%alat
          write(unit, '(A, A)') "crystal_sym = ", "'" // trim(lat%crystal_sym) // "'"
          write(unit, '(A, F12.10)') "wav = ", lat%wav
          write(unit, '(A, L)') "pbc = ", lat%pbc
          if (lat%pbc) then
            write(unit, '(A, L)') "b1 = ", lat%b1
            write(unit, '(A, L)') "b2 = ", lat%b2
            write(unit, '(A, L)') "b3 = ", lat%b3
            write(unit, '(A, I0)') "n1 = ", lat%n1
            write(unit, '(A, I0)') "n2 = ", lat%n2
            write(unit, '(A, I0)') "n3 = ", lat%n3
            ndim = lat%n1 * lat%n2 * lat%n3 * lat%ntype
            write(unit, '(A, I0)') "ndim = ", ndim
          else
            write(unit, '(A, I0)') "ndim = ", lat%ndim
          end if
          write(unit, '(A, I0)') "ntype = ", lat%ntype
          write(unit, '(A, F10.4)') "ct(:) = ", lat%ct(1)
          write(unit, '(A, F12.6)') "r2 = ", lat%r2
          if (lat%has_strux_backend) write(unit, '(A, A, A)') "strux_backend = '", trim(lat%strux_backend), "'"
          if (lat%has_strux_want_sdot) write(unit, '(A, L1)') "strux_want_sdot = ", lat%strux_want_sdot
          if (lat%has_strux_solve_scale) write(unit, '(A, ES16.8)') "strux_solve_scale = ", lat%strux_solve_scale
          write(unit, '(A)') "/"
        end subroutine write_lattice_block

        function lower_string(s) result(out)
          character(len=*), intent(in) :: s
          character(len=len(s)) :: out
          integer :: ii, ia
          out = s
          do ii = 1, len(s)
            ia = iachar(s(ii:ii))
            if (ia >= iachar('A') .and. ia <= iachar('Z')) out(ii:ii) = achar(ia + 32)
          end do
        end function lower_string
      
        !-------------------------------------------------------------------
        ! Subroutine: write_atoms_block
        !
        ! Writes the updated &atoms namelist block to the output unit.
        !-------------------------------------------------------------------
        subroutine write_atoms_block(unit, atm)
          implicit none
          integer, intent(in) :: unit
          type(atoms_inp), intent(in) :: atm
          integer :: i, n_labels
          character(len=16) :: idx_str
      
          n_labels = size(atm%label)
          write(unit, '(A)') "&atoms"
          write(unit, '(A, A)') "database = ", "'" // trim(atm%database) // "'"
          do i = 1, n_labels
             idx_str = int_to_string(i)
             write(unit, '(A)') "label(" // trim(idx_str) // ") = '" // trim(atm%label(i)) // "'"
          end do
          write(unit, '(A)') "/"
        end subroutine write_atoms_block
      
        !-------------------------------------------------------------------
        ! Subroutine: update_namelists
        !
        ! Reads the input file line by line and writes it to the output file,
        ! replacing the &lattice and &atoms blocks with updated versions.
        !-------------------------------------------------------------------
        subroutine update_namelists(in_file, out_file, new_lat, new_atm, rec_steps, rec_method)
          implicit none
          character(len=*), intent(in) :: in_file
          character(len=*), intent(in) :: out_file
          type(lattice_inp), intent(in) :: new_lat
          type(atoms_inp), intent(in)   :: new_atm
          integer, intent(in) :: rec_steps
          character(len=*), intent(in) :: rec_method
      
          integer :: in_unit, out_unit, ios
          character(len=1024) :: line, lower_line
          logical :: control_found, wrote_lld, wrote_recur
      
          ! Open the input and output files.
          open(newunit=in_unit, file=in_file, status='old', action='read', iostat=ios)
          if (ios /= 0) then
            print *, "Error opening input file: ", trim(in_file)
            call write_lattice_atoms_file(out_file, new_lat, new_atm)
          end if
          open(newunit=out_unit, file=out_file, status='replace', action='write', iostat=ios)
          if (ios /= 0) then
            print *, "Error opening output file: ", trim(out_file)
            stop
          end if
      
          control_found = .false.

          do
            read(in_unit, '(A)', iostat=ios) line
            if (ios < 0) exit  ! End-of-file reached.
            if (ios > 0) then
               print *, "Error reading input file."
               stop
            end if
      
            ! Check if this line begins the &lattice block.
            if (index(trim(adjustl(line)), "&lattice") == 1) then
               ! Write the updated lattice block.
               call write_lattice_block(out_unit, new_lat)
               ! Skip the old lattice block until the "/" marker is found.
               do
                  read(in_unit, '(A)', iostat=ios) line
                  if (ios < 0) exit
                  if (trim(line) == "/") exit
               end do
               cycle  ! Continue with next line after the block.
            end if
      
            ! Check if this line begins the &atoms block.
            if (index(trim(adjustl(line)), "&atoms") == 1) then
               call write_atoms_block(out_unit, new_atm)
               do
                  read(in_unit, '(A)', iostat=ios) line
                  if (ios < 0) exit
                  if (trim(line) == "/") exit
               end do
               cycle
            end if

            ! Check if this line begins the &control block.
            lower_line = lower_string(trim(adjustl(line)))
            if (index(lower_line, "&control") == 1) then
               control_found = .true.
               wrote_lld = .false.
               wrote_recur = .false.
               write(out_unit, '(A)') "&control"
               do
                  read(in_unit, '(A)', iostat=ios) line
                  if (ios < 0) exit
                  lower_line = lower_string(trim(adjustl(line)))
                  if (trim(adjustl(line)) == "/") then
                     if (.not. wrote_lld) write(out_unit, '(A, I0)') "lld = ", rec_steps
                     if (.not. wrote_recur) write(out_unit, '(A, A, A)') "recur = '", trim(rec_method), "'"
                     write(out_unit, '(A)') "/"
                     exit
                  end if
                  if (index(lower_line, "lld") == 1) then
                     write(out_unit, '(A, I0)') "lld = ", rec_steps
                     wrote_lld = .true.
                     cycle
                  end if
                  if (index(lower_line, "recur") == 1) then
                     write(out_unit, '(A, A, A)') "recur = '", trim(rec_method), "'"
                     wrote_recur = .true.
                     cycle
                  end if
                  write(out_unit, '(A)') trim(line)
               end do
               cycle
            end if
      
            ! Otherwise, copy the line unchanged.
            write(out_unit, '(A)') trim(line)
          end do

          if (.not. control_found) then
             write(out_unit, '(A)') ""
             write(out_unit, '(A)') "&control"
             write(out_unit, '(A, A, A)') "recur = '", trim(rec_method), "'"
             write(out_unit, '(A, I0)') "lld = ", rec_steps
             write(out_unit, '(A)') "/"
          end if
      
          close(in_unit)
          close(out_unit)


        write(*,*) "Input namelist written to file:", out_file
        end subroutine update_namelists
      

  !-------------------------------------------------------------------
  ! Subroutine: write_lattice_atoms_file
  !
  ! Writes only the &lattice and &atoms namelist blocks to the given file.
  ! This routine can be used to create a new input file if none exists.
  !-------------------------------------------------------------------
        subroutine write_lattice_atoms_file(out_file, new_lat, new_atm)
            implicit none
            character(len=*), intent(in) :: out_file
            type(lattice_inp), intent(in) :: new_lat
            type(atoms_inp), intent(in) :: new_atm
            integer :: out_unit, ios
        
            ! Open the output file (replace if it already exists).
            open(newunit=out_unit, file=out_file, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
              print *, "Error opening output file: ", trim(out_file)
              stop 1
            end if
        
            ! Write only the &lattice and &atoms blocks.
            call write_lattice_block(out_unit, new_lat)
            call write_atoms_block(out_unit, new_atm)
        
            close(out_unit)
          end subroutine write_lattice_atoms_file

    !-------------------------------------------------------------------
    ! Function: get_atomic_numbers
    !
    ! Takes an array of atom labels and returns an array of atomic numbers.
    !-------------------------------------------------------------------
    function get_atomic_numbers(labels, ntype) result(atomic_numbers)
      implicit none
      character(len=*), intent(in) :: labels(:)
      integer, intent(in) :: ntype
      integer :: atomic_numbers(ntype)
      integer :: i, j, pos
      character(len=2) :: sym

      ! Define a lookup table of element symbols (all two characters).
      character(len=2), parameter :: element_symbols(119) = (/ 'E ',    &
           'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
           'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
           'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
           'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
           'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
           'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
           'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
           'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
           'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', &
           'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og' /)

      ! Define a corresponding lookup table of atomic numbers.
      integer, parameter :: atomic_numbers_table(119) = (/ &
           0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
           11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
           21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
           31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
           41, 42, 43, 44, 45, 46, 47, 48, 49, 50, &
           51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
           61, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
           71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
           81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
           91, 92, 93, 94, 95, 96, 97, 98, 99, 100, &
           101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
           111, 112, 113, 114, 115, 116, 117, 118 /)

      do i = 1, ntype
         ! Extract the first two characters from the label.
         sym = adjustl(labels(i)(1:2))
         ! Remove digits from label
         pos = scan(sym, '0123456789') 
         if ( pos > 1) then
            sym = sym(1:pos-1)
         endif
         ! Trim and pad the symbol to length 2.
         sym = trim(sym)
         if (len(sym) < 2) then
            sym = sym // ' '
         endif

         atomic_numbers(i) = -1
         do j = 1, 118
            if (trim(element_symbols(j)) == sym) then
               atomic_numbers(i) = atomic_numbers_table(j)
               exit
            end if
         end do
         if (atomic_numbers(i) < 0) then
            print *, "Unknown element: ", sym
            stop 1
         end if
      end do

    end function get_atomic_numbers

    pure integer function expected_highest_f_electrons(Z) result(nf)
  implicit none
  integer, intent(in) :: Z

  if (Z < 58) then
     ! For elements before the lanthanide series, no f electrons are occupied.
     nf = 0
  else if (Z >= 58 .and. Z <= 71) then
     ! Lanthanides: highest occupied shell is 4f; occupancy increases linearly.
     nf = Z - 57
  else if (Z >= 72 .and. Z <= 90) then
     ! 5d transition metals: the 4f shell is fully occupied.
     nf = 14
  else if (Z >= 91 .and. Z <= 103) then
     ! Actinides: the highest occupied shell is 5f (the 4f remains fully occupied as core).
     nf = Z - 89
  else
     ! For elements outside these ranges, return 0.
     nf = 0
  end if

end function expected_highest_f_electrons

end module input_io
