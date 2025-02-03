module element_info_mod
  implicit none
  private
  public :: get_element_info

  ! Number of elements in our lookup (119 entries).
  integer, parameter :: nElements = 105

  ! The lookup table of element symbols.
  ! Each element symbol is stored as a two-character string.
  character(len=2), parameter :: element_symbols(nElements) = (/ &
       'E ', 'H ', 'He','Li','Be','B ', 'C ', 'N ', 'O ', 'F ', &
       'Ne','Na','Mg','Al','Si','P ', 'S ', 'Cl','Ar','K ', &
       'Ca','Sc','Ti','V ', 'Cr','Mn','Fe','Co','Ni','Cu', &
       'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ', &
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In', &
       'Sn','Sb','Te','I ', 'Xe','Cs','Ba','La','Ce','Pr', &
       'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm', &
       'Yb','Lu','Hf','Ta','W ', 'Re','Os','Ir','Pt','Au', &
       'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac', &
       'Th','Pa','U ', 'Np','Pu','Am','Cm','Bk','Cf','Es', &
       'Fm','Md','No','Lr','Rf'/)
       !'Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt', &
       !'Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' /)

  ! The corresponding lookup table for the IZCOR values.
  ! This table is taken from the original data statement:
  !   data IZCOR/0, 0, 0, 6*2, 4, 4, 6*10, 12, 12, 12*18, 4*28, 30, 30, &
  !             12*36, 4*46, 48, 48, 16*54, 10*68, 4*78, 80, 80, 18*86/
  integer, parameter :: izcor(nElements) = (/ &
       0, 0, 0, &           ! 1-3: 0,0,0
       2,2,2,2,2,2, &        ! 4-9: six copies of 2
       4,4, &                ! 10-11: two copies of 4
       10,10,10,10,10,10, &   ! 12-17: six copies of 10
       12,12, &              ! 18-19: two copies of 12
       18,18,18,18,18,18,18,18,18,18,18,18, & ! 20-31: twelve copies of 18
       28,28,28,28, &        ! 32-35: four copies of 28
       30,30, &              ! 36-37: two copies of 30
       36,36,36,36,36,36,36,36,36,36,36,36, & ! 38-49: twelve copies of 36
       46,46,46,46, &        ! 50-53: four copies of 46
       48,48, &              ! 54-55: two copies of 48
       54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54, & ! 56-71: sixteen copies of 54
       68,68,68,68,68,68,68,68,68,68, & ! 72-81: ten copies of 68
       78,78,78,78, &        ! 82-85: four copies of 78
       80,80, &              ! 86-87: two copies of 80
       86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86 /) ! 88-105: eighteen copies of 86

contains

  !-------------------------------------------------------------------
  ! Subroutine: get_element_info
  !
  ! Given an element symbol TXT (2 characters) and an atomic index IZ,
  ! this routine returns the (zero-based) atomic number in IZ and the
  ! corresponding lookup value (IZC) from our izcor array.
  !
  ! If IZ is provided as a nonnegative number, the routine simply returns:
  !    IZC = izcor(IZ+1)
  !
  ! If IZ is negative, the routine searches for TXT in our lookup table.
  !-------------------------------------------------------------------
  subroutine get_element_info(TXT, IZ, IZC)
    implicit none
    character(len=*), intent(in) :: TXT  ! Input element symbol (at least 2 characters)
    integer, intent(inout) :: IZ         ! Atomic index (if < 0, then determined from TXT)
    integer, intent(out) :: IZC          ! Corresponding izcor value
    integer :: i
    character(len=2) :: symbol

    if (IZ >= 0) then
      if (IZ < nElements) then
         IZC = izcor(IZ + 1)
         return
      else
         print *, "Error: atomic index out of range in get_element_info."
         stop 1
      end if
    end if

    ! Use the first two characters of TXT.
    symbol = adjustl(TXT(1:2))
    ! If the second character is an underscore, replace it with a space.
    if (symbol(2:2) == '_') symbol(2:2) = ' '
    do i = 1, nElements
       if (trim(element_symbols(i)) == trim(symbol)) then
          IZ = i - 1
          IZC = izcor(i)
          return
       end if
    end do

    print *, "Error: Element not found: ", TXT
    stop 1
  end subroutine get_element_info

end module element_info_mod
