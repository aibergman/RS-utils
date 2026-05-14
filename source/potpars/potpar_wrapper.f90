module potpar_wrapper
  use atom_io
  use lmto_lib
  use system_data
  use input_io, only : expected_highest_f_electrons
  implicit none
contains

  subroutine update_atoms(pot_data, atm_input, lat_data, lat_input)
    implicit none

    type(potential_t), intent(out) :: pot_data
    type(atoms_inp), intent(in) :: atm_input
    type(lattice_t), intent(in) :: lat_data
    type(lattice_inp), intent(in) :: lat_input

    ! Local variables
    integer :: I
    type(element_t) :: element
    character(len=20) :: filename
    
    do I = 1, lat_data%nbulk_bulk
      call process_one_atom(pot_data, atm_input, lat_data, lat_input, element, I)
      !  write(*,*) atm_input%label(I)
      !  write(*,*) element%symbol
      filename = trim(adjustl(atm_input%label(I))) // ".nml"
      call write_atom(filename, element, pot_data)
    end do

  end subroutine update_atoms

  subroutine process_one_atom(pot_data, atm_input, lat_data, lat_input, element, iatom)
  !subroutine process_one_atom(pot_data, z, wsr, wavg)
    implicit none
    type(potential_t), intent(out) :: pot_data
    type(atoms_inp), intent(in) :: atm_input
    type(lattice_t), intent(in) :: lat_data
    type(lattice_inp), intent(in) :: lat_input
    type(element_t), intent(out) :: element
    integer, intent(in) :: iatom

    ! Local variables
    type(lmto_t) :: lm

    print *,'-------------------------------------------------'
    print *,'Processing atom', iatom
    print *,'Atomic number:', atm_input%atomic_number(iatom)
    print *,'Wigner-Seitz radius:', lat_data%ws(iatom)
    print *,'Label:', atm_input%label(iatom)
    print *,'-------------------------------------------------'

    lm%z(1) = 1.0d0*atm_input%atomic_number(iatom)
    lm%wsr(1) = lat_data%ws(iatom)
    lm%avw = lat_input%wav

    ! Initialization (one atom at a time, configurable up to Lmax=3)
    lm%nclass = 1
    lm%nl = lat_input%lmax + 1
    if (lat_input%lmax < 0 .or. lat_input%lmax > 3) then
      print *, "Error: requested lmax=", lat_input%lmax, " but this phase supports 0..3."
      stop 1
    end if
    lm%nsp = 2
    lm%lafm = .false.
    lm%lfree = .false.
    allocate(lm%pnu(0:lm%nl-1, lm%nsp, lm%nclass))
    allocate(lm%qnu(0:2, 0:lm%nl-1, lm%nsp, lm%nclass))
    allocate(lm%pp(6, 0:lm%nl-1, lm%nsp, lm%nclass))
    allocate(lm%pdum(0:lm%nl-1, lm%nsp, lm%nclass))

    ! Transform and compute moments (defcnf & defmom)
    call perform_transformations(lm%konfig, lm%lmx, lm%nclass, lm%nl, lm%nsp, lm%pnu, lm%qnu, lm%z, lm%lafm, lm%lfree)

    ! Compute potential parameters (defpp and related processing)
    call compute_potential(lm%avw, lm%nclass, lm%nl, lm%nsp, lm%pp, lm%wsr, lm%z, lm%pdum)

    ! Store the potential parameters in a potential_t instance
    call fill_potential_t(pot_data, element, atm_input, lm, iatom)

    ! Calculate TB parameters from ORT
    call screen_ort(pot_data, lm%avw)

    ! Print the results
    ! call print_results(lm%pnu, lm%qnu, lm%pp)

    ! Write namlist atom file
    ! call write_atom('atom_out.nml', element, pot_data)
  end subroutine process_one_atom
    

  subroutine get_user_input(z, wsr)
    implicit none
    ! Reads user input for the atom number (z) and the Wigner-Seitz radius (wsr),
    ! then sets avw equal to wsr.
    real(kind=8), intent(out) :: z(1)
    real(kind=8), intent(out) :: wsr(1)

    print *, 'Enter the atom number (z):'
    read(*,*) z(1)
    print *, 'Enter the Wigner-Seitz radius (wsr):'
    read(*,*) wsr(1)
  end subroutine get_user_input

  subroutine perform_transformations(konfig, lmx, nclass, nl, nsp, pnu, qnu, z, lafm, lfree)
    implicit none
    ! Calls defcnf and defmom from lm47_utils and adjusts the quantum numbers.
    integer, intent(inout) :: konfig(0:3, 1)
    integer, intent(inout) :: lmx(1)
    integer, intent(in) :: nclass, nl, nsp
    real(kind=8), intent(inout) :: pnu(0:, :, :)
    real(kind=8), intent(inout) :: qnu(0:, 0:, :, :)
    real(kind=8), intent(in) :: z(1)
    logical, intent(inout) :: lafm, lfree

    call defcnf(konfig, lmx(1), z(1))
    lmx(1) = min(nl - 1, 3)
    ! print *, 'Calling defmom', lfree
    ! lfree = .true.
    call defmom(konfig, lafm, lfree, lmx, nclass, nl, nsp, pnu, qnu, z)
    ! Adding finite values to second order moments for stability:
    qnu(2, :, :, :) = qnu(0, :, :, :) * 0.005d0
  end subroutine perform_transformations

  subroutine compute_potential(avw, nclass, nl, nsp, pp, wsr, z, pdum)
    implicit none
    ! Calls defpp from lm47_utils and then reorders and transforms the potential parameters.
    real(kind=8), intent(in) :: avw
    integer, intent(in) :: nclass, nl, nsp
    real(kind=8), intent(inout) :: pp(:, 0:, :, :)
    real(kind=8), intent(in) :: wsr(1), z(1)
    real(kind=8), intent(inout) :: pdum(0:, :, :)
    ! print *, 'Calling defpp'
    call defpp(avw, nclass, nl, nsp, pp, wsr, z)
    ! Transform from LMTO47 to RS-LMTO convention:
    pp(4, :, :, :) = 1.0d0 / sqrt(pp(4, :, :, :))
    ! Swap the 4th and 5th components:
    pdum = pp(4, :, :, :)
    pp(4, :, :, :) = pp(5, :, :, :)
    pp(5, :, :, :) = pdum
  end subroutine compute_potential

  subroutine print_results(pnu, qnu, pp)
    ! Prints the computed moments and potential parameters.
    real(kind=8), intent(in) :: pnu(0:, :, :)
    real(kind=8), intent(in) :: qnu(0:, 0:, :, :)
    real(kind=8), intent(in) :: pp(:, 0:, :, :)

    print *, 'The result of defmom is:'
    print '(a, 6f10.6)', 'PNU up:', pnu(:, 1, 1)
    print '(a, 6f10.6)', 'PNU dw:', pnu(:, 2, 1)
    print '(a, 6f10.6)', 'QNU s up:', qnu(:, 0, 1, 1)
    print '(a, 6f10.6)', 'QNU p up:', qnu(:, 1, 1, 1)
    print '(a, 6f10.6)', 'QNU d up:', qnu(:, 2, 1, 1)
    print '(a, 6f10.6)', 'QNU s dw:', qnu(:, 0, 2, 1)
    print '(a, 6f10.6)', 'QNU p dw:', qnu(:, 1, 2, 1)
    print '(a, 6f10.6)', 'QNU d dw:', qnu(:, 2, 2, 1)
    print *, ''
    print *, 'The result of defpp is:'
    print *, '            E_nu        C      SqDel     Qpar       Ppar      d_nu'
    print '(a, 6f10.6)', 'PP s up:', pp(:, 0, 1, 1)
    print '(a, 6f10.6)', 'PP p up:', pp(:, 1, 1, 1)
    print '(a, 6f10.6)', 'PP d up:', pp(:, 2, 1, 1)
    print '(a, 6f10.6)', 'PP s dw:', pp(:, 0, 2, 1)
    print '(a, 6f10.6)', 'PP p dw:', pp(:, 1, 2, 1)
    print '(a, 6f10.6)', 'PP d dw:', pp(:, 2, 2, 1)
  end subroutine print_results

  subroutine fill_potential_t(potential, element, atom, lm, iatom)
    ! Fills a potential_t instance with data from potpar_wrapper
    type(potential_t), intent(out) :: potential
    type(element_t), intent(out) :: element
    type(atoms_inp), intent(in) :: atom
    type(lmto_t), intent(in) :: lm
    integer, intent(in) :: iatom

    integer :: i, j, k, nf

    potential%lmax = lm%nl - 1
    call allocate_potential_arrays(potential, potential%lmax)

    ! Potential function
    potential%pl = lm%pnu(:,:,1)

    ! Charge and band moments
    potential%ql = lm%qnu(:,:,:,1)

    ! Wigner-Seitz radius
    potential%ws_r = lm%wsr(1)

    ! ORT potential parameters
    potential%enu = lm%pp(1, :, :, 1)
    potential%c = lm%pp(2, :, :, 1)
    potential%srdel = lm%pp(3, :, :, 1)
    potential%qpar = lm%pp(4, :, :, 1)
    potential%ppar = lm%pp(5, :, :, 1)
    ! Actually vl is not pp(6) but it is not an input anyway
    potential%vl =  lm%pp(6, :, :, 1)
    !potential%dnu = lm%pp(6, :, :, 1)

    element%symbol = atom%label(iatom)
    element%atomic_number = atom%atomic_number(iatom)
    element%valence = nint(1000.0d0*sum(potential%ql(1,:,:)))/1000.0d0
    ! print *,'AB valence:', potential%ql(1,:,:)
    ! print *,'AB sum valence:', sum(potential%ql(1,:,:))
    element%core = element%atomic_number - element%valence
    !print *,'AB core:', element%atomic_number, element%valence

    element%num_quant_s = int(potential%pl(1,1))
    element%num_quant_p = 0
    element%num_quant_d = 0
    if (size(potential%pl, 1) >= 2) element%num_quant_p = int(potential%pl(2,1))
    if (size(potential%pl, 1) >= 3) element%num_quant_d = int(potential%pl(3,1))
    ! Add small values to moments for empty spheres
    if(element%atomic_number==0)  then
      call random_number(potential%ql)
      potential%ql = potential%ql * 1.0d-6
    end if

    element%f_core = expected_highest_f_electrons(element%atomic_number)

    !!! ! Pretty print the filled values
    !!! print *, 'Potential parameters:'
    !!! print *, 'Lmax:', potential%lmax
    !!! print *, 'PNU:'
    !!! do i = 1, size(potential%pl, 1)
    !!!   print '(a, 3f10.6)', '  ', potential%pl(i, :)
    !!! end do
    !!! print *, 'QNU:'
    !!! do i = 1, size(potential%ql, 1)
    !!!   do j = 1, size(potential%ql, 2)
    !!!   print '(a, 3f10.6)', '  ', potential%ql(i, j, :)
    !!!   end do
    !!! end do
    !!! print *, 'ENU:'
    !!! do i = 1, size(potential%enu, 1)
    !!!   print '(a, 3f10.6)', '  ', potential%enu(i, :)
    !!! end do
    !!! print *, 'C:'
    !!! do i = 1, size(potential%c, 1)
    !!!   print '(a, 3f10.6)', '  ', potential%c(i, :)
    !!! end do
    !!! print *, 'SRDEL:'
    !!! do i = 1, size(potential%srdel, 1)
    !!!   print '(a, 3f10.6)', '  ', potential%srdel(i, :)
    !!! end do
    !!! print *, 'QPAR:'
    !!! do i = 1, size(potential%qpar, 1)
    !!!   print '(a, 3f10.6)', '  ', potential%qpar(i, :)
    !!! end do
    !!! print *, 'PPAR:'
    !!! do i = 1, size(potential%ppar, 1)
    !!!   print '(a, 3f10.6)', '  ', potential%ppar(i, :)
    !!! end do

    !!! ! Pretty print the filled values for element
    !!! print *, 'Element properties:'
    !!! print *, 'Symbol:', element%symbol
    !!! print *, 'Atomic number:', element%atomic_number
    !!! print *, 'Valence:', element%valence
    !!! print *, 'Core:', element%core
    !!! print *, 'Number of s quantum states:', element%num_quant_s
    !!! print *, 'Number of p quantum states:', element%num_quant_p
    !!! print *, 'Number of d quantum states:', element%num_quant_d
    !!! print *, 'f core:', element%f_core
  end subroutine fill_potential_t

end module potpar_wrapper
