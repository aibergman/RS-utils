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
    subroutine read_input_namelists(filename, lat, atm, ntype)
        character(len=*), intent(in)  :: filename
        type(lattice_inp), intent(out) :: lat
        type(atoms_inp), intent(out) :: atm
        integer, intent(in) :: ntype

        ! Local variables
        integer :: unit, ioerr

        !-------------------------------------------------------------------
        ! Temporary variables for reading the lattice namelist.
        ! The names here must match the variable names in the input file.
        !-------------------------------------------------------------------
        integer :: rc
        real(kind=8) :: alat
        character(len=8) :: crystal_sym
        real(kind=8) :: wav
        logical :: pbc, b1, b2, b3
        integer :: n1, n2, n3

        namelist /lattice/ rc, alat, crystal_sym, wav, pbc, b1, b2, b3, n1, n2, n3

        !-------------------------------------------------------------------
        ! Temporary variables for reading the atoms namelist.
        !-------------------------------------------------------------------
        character(len=256) :: database   ! directory of the database
        character(len=8), allocatable :: label(:)
        namelist /atoms/ database, label

        ! Open the file for reading.
        open(newunit=unit, file=filename, status='old', action='read', iostat=ioerr)
        if (ioerr /= 0) then
            print *, "Error opening file: ", trim(filename)
            stop 1
        end if

        ! Re-allocate atom label array
        allocate(atm%label(ntype))
        allocate(atm%atomic_number(ntype))
        allocate(label(ntype))

        ! Read the &lattice namelist.
        read(unit, nml=lattice)
        ! Read the &atoms namelist.
        read(unit, nml=atoms)
        close(unit)

        ! Transfer the read values into the derived type instances.
        lat%rc = rc
        lat%alat = alat
        lat%crystal_sym = crystal_sym
        lat%wav = wav

        atm%database = database
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
          write(unit, '(A)') "&lattice"
          write(unit, '(A, I0)') "rc = ", lat%rc
          write(unit, '(A, F12.6)') "alat = ", lat%alat
          write(unit, '(A, A)') "crystal_sym = ", "'" // trim(lat%crystal_sym) // "'"
          write(unit, '(A, F12.10)') "wav = ", lat%wav
          write(unit, '(A)') "/"
        end subroutine write_lattice_block
      
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
        subroutine update_namelists(in_file, out_file, new_lat, new_atm)
          implicit none
          character(len=*), intent(in) :: in_file
          character(len=*), intent(in) :: out_file
          type(lattice_inp), intent(in) :: new_lat
          type(atoms_inp), intent(in)   :: new_atm
      
          integer :: in_unit, out_unit, ios
          character(len=1024) :: line
      
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
      
            ! Otherwise, copy the line unchanged.
            write(out_unit, '(A)') trim(line)
          end do
      
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
      integer :: i, j
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
         ! If the first character is a digit, take only the first character.
         if (index('0123456789', sym(1:1)) > 0) then
            sym = adjustl(labels(i)(1:1))
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
