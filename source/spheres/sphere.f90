module spheres
    use system_data
    implicit none

    ! Parameters
    integer, parameter :: NTMAX = 200
    integer, parameter :: NRMAX = 300
    integer, parameter :: NRAD1 = 251
    integer, parameter :: MAXDIM = 1000000
    integer, parameter :: MAX2SORT = NTMAX
    real*8, parameter :: ANGSTROM_TO_BOHR = 1.8897259886

    ! Common variables
    real*8 :: AV(3), BV(3), CV(3), QLAT(3, 3), RMAXES, RMINES
    integer :: N1ES, N2ES, N3ES, NATES
    common / EMPTYS / NATES, N1ES, N2ES, N3ES, RMINES, RMAXES
    common / L2LAT / AV, BV, CV, QLAT

    ! Shared variables
    real*8 :: ALAT, DPAS, R0ACT
    real*8 :: BAS(3, NTMAX)
    real*8 :: BASNEW(3, MAX2SORT)
    real*8 :: CONC(NTMAX)
    real*8 :: RWS(NTMAX)
    real*8 :: SES(MAX2SORT)
    real*8 :: Z(NTMAX)

    logical :: DB
    ! real*8, external :: DNEVMOD
    character(100) :: FFF
    integer :: IMQ(NTMAX), IMT(NTMAX)
    integer :: IQA(48), ITOP(48), NG
    integer :: ISNEW(MAX2SORT)
    integer :: NATOM, NATOMNEW, NDNEV
    integer :: NEED_ES 
    integer :: NM, NQ, NRADMAX, NSORT, NT
    integer :: W(MAXDIM)
    character(4) :: TXTT(NTMAX)

contains

    subroutine init_spheres()
        implicit none

        ! Initialization of parameters
        R0ACT = 1.D-2
        NRADMAX = NRMAX - 1
        DPAS = 0.05D0
        NDNEV = 14
        ! DB = .true.
        DB = .false.
        N1ES = 24
        N2ES = 24
        N3ES = 24

    end subroutine init_spheres


    subroutine parse_structure(lat_inp, at_inp, lat_data)
        use rsutil_config, only: enable_empty_spheres, min_empty_sphere_radius, max_empty_sphere_radius
        implicit none

        ! Arguments
        type(lattice_inp), intent(in) :: lat_inp
        type(atoms_inp), intent(in) :: at_inp
        type(lattice_t), intent(in) :: lat_data
        ! real*8, intent(in) :: a(3, 3)
        ! real*8, intent(in) :: crd(3, nbulk_bulk)
        ! integer, intent(in) :: nbulk_bulk


        ! Local variables
        integer :: ISP1, ISP2, I, J, IDUM, INI

        ! Set NEED_ES based on enable_empty_spheres flag
        ! NEED_ES = 0 means calculate empty spheres, NEED_ES = 1 means skip them
        if (enable_empty_spheres) then
            NEED_ES = 0
        else
            NEED_ES = 1
        end if
        
        ! Old defaults: RMINES=1.2D0, RMAXES=2.5D0
        RMINES = min_empty_sphere_radius
        RMAXES = max_empty_sphere_radius

        ALAT = lat_inp%alat*ANGSTROM_TO_BOHR
        if (DB) print *, ' ALAT:', ALAT

        ! As usual, consider transpose..
        AV = lat_data%a(:, 1)
        BV = lat_data%a(:, 2)
        CV = lat_data%a(:, 3)

        if (DB) print *, 'AV:', (AV(I), I=1, 3)
        if (DB) print *, 'BV:', (BV(I), I=1, 3)
        if (DB) print *, 'CV:', (CV(I), I=1, 3)

        NQ = lat_data%nbulk_bulk
        if (DB) print *, 'N*:', NQ

        do I = 1, NQ
            IMQ(I) = I
            BAS(:, I) = lat_data%crd(:, I)
            if (DB) print *, I, "BAS", BAS(1, I), BAS(2, I), BAS(3, I), "IMQ", IMQ(I)
        end do

        NM = NQ
        if (DB) print *, 'NM:', NM

        NSORT = NQ
        if (DB) print *, 'NSORT:', NSORT

        NT = NSORT
        do I = 1, NSORT
            Z(I) = at_inp%atomic_number(I)
            TXTT(I) = at_inp%label(I)
            CONC(I) = 1.0D0
            IMT(I) = I
            if (DB) print *, I, "Z,TXTT,CONC,IMT:", Z(I), TXTT(I), CONC(I), IMT(I)
        end do

        ! Hard-wired number of symmetry operations to one
        NG = 1
        ITOP(NG) = 1
        IQA(NG) = 1

        if (DB) then
            print *, '===================================================='
            print*, 'PARAMETERS:', NT, NQ, NM, 'NT,NQ,NM,-- ALAT:', ALAT
            print*, 'Z(i),TXTT(i),CONC(i),IMT(i)'
            do I = 1, NT
                print*, Z(I), TXTT(I), CONC(I), IMT(I)
            end do
            print*, 'qx(i),qy(i),qz(i),IMQ(i)'
            do I = 1, NQ
                print*, (BAS(J, I), J=1, 3), IMQ(I)
            end do
            print*, 'Lattice vectors:'
            print*, AV
            print*, BV
            print*, CV
            print *, '===================================================='
        end if

        NG = 1

        ! Very important row below.
        if (1 == 2) write (*, *) IDUM
    end subroutine parse_structure


    subroutine perform_calculations()
        use rsutil_config, only: enable_empty_spheres, empty_sphere_backend, min_empty_sphere_radius, max_empty_sphere_radius
        use questaal_findes_backend, only: add_questaal_empty_spheres
        implicit none

        ! Local variables
        integer :: I, J, IPNT, IT, ITYPE, ISORT, IATOM, IST
        integer :: NSORTES

        real*8 :: AMTC(NRMAX, NTMAX) 
        real*8 :: R0(NTMAX)
        real*8 :: R0SITE(MAX2SORT)
        real*8 :: WS
        real*8 :: WSREST(MAX2SORT)
        real*8 :: RHOSITE(NRAD1, MAX2SORT)
        real*8 :: RAD, RAT(NRAD1, MAX2SORT)
        real*8 :: RINT, TAUES(3, MAX2SORT)
        real*8 :: ZN, ZZ(MAX2SORT)

        real*8 :: DNEVMOD
        real*8, allocatable :: quest_basis(:,:), quest_radii(:)
        integer, allocatable :: quest_species(:)
        integer :: quest_status
        open(unit=16, file='spheres.log', status='unknown', action='write')
        call WKINIT(MAXDIM)
        call GETGBASIS()

        R0 = 9999D0
        AMTC = 9999D0
        RHOSITE(:, :) = 9999D0

        write(16, *) "TXTT:", TXTT(1:NT)
        call READRO(NT, AMTC, Z, R0, TXTT, RHOSITE, NRADMAX, NRAD1)

        do IT = 1, NT
            do I = 1, NRAD1
                RAT(I, IT) = R0(IT) * exp(dble(I - 1) * DPAS)
            end do
        end do

        do I = 1, NM
            WSREST(I) = 0
            ZZ(I) = 0
            R0SITE(I) = R0ACT
        end do

        RHOSITE(:, :) = 0.0D0

        do ITYPE = 1, NT
            ISORT = IMQ(IMT(ITYPE))
            ZN = Z(ITYPE)
            call DEFWSR(WS, ZN)
            if (DB) write(16,*) 'type:', ITYPE, ' Z=', ZN, ' Rad=', WS, ' conc=', CONC(ITYPE)
            WSREST(ISORT) = WSREST(ISORT) + WS * CONC(ITYPE)
            ZZ(ISORT) = ZZ(ISORT) + ZN * CONC(ITYPE)
            do IPNT = 1, NRAD1
                RAD = R0ACT * exp(dble(IPNT - 1) * DPAS)
                if (RAD < RAT(NRAD1 - 2, ITYPE)) then
                    RINT = max(0.D0, DNEVMOD(RAT(1, ITYPE), AMTC(1, ITYPE), RAD, NRAD1, NDNEV))
                else
                    RINT = 0
                end if
                RHOSITE(IPNT, ISORT) = RHOSITE(IPNT, ISORT) + CONC(ITYPE) * RINT
            end do
        end do

        do I = 1, NM
            if (DB) print*, ' Site:', I, '  RaD:', WSREST(I), ' Z:', ZZ(I)
        end do

        NSORT = NM
        NATOM = NQ

        if (DB) then
            print*, 'av=', AV
            print*, 'bv=', BV
            print*, 'cv=', CV
            print*, 'imq', (IMQ(I), I=1, NATOM)
            print*, 'BAS'
            do I = 1, NATOM
                print*, (BAS(J, I), J=1, 3)
            end do
            do IST = 1, NSORT
                print*, 'R0S', R0SITE(IST), ' ZZ', ZZ(IST), WSREST(IST)
            end do
            print*, 'dpas,alat,natom,nsort:', DPAS, ALAT, NATOM, NSORT
        end if

        call DEFMTR(W, NRAD1, AV, BV, CV, IMQ, BAS, RHOSITE, R0SITE, DPAS, ALAT, NATOM, NSORT, ZZ, RWS, NEED_ES, NRAD1, WSREST)

        NATOMNEW = NATOM
        
        if (NEED_ES == 0 .and. enable_empty_spheres) then
            if (trim(empty_sphere_backend) == "questaal") then
                call add_questaal_empty_spheres(reshape([AV,BV,CV],[3,3]), ALAT, BAS(:,1:NATOM), Z(1:NATOM), &
                    RWS(1:NATOM), RMINES, RMAXES, quest_basis, quest_species, quest_radii, quest_status)
                if (quest_status /= 0) error stop 'Questaal FINDES backend rejected this lattice'
                NATOMNEW = size(quest_radii)
                BASNEW(:,1:NATOMNEW) = quest_basis
                ISNEW(1:NATOMNEW) = quest_species
                do I = 1, NATOMNEW
                    IMQ(I) = quest_species(I)
                    RWS(I) = quest_radii(I)
                end do
                NSORTES = NATOMNEW-NATOM
                if (NSORTES > 0) SES(1:NSORTES) = quest_radii(NATOM+1:NATOMNEW)
                deallocate(quest_basis, quest_species, quest_radii)
            else
                call EMPTYSPHERES(AV, BV, CV, BAS, IMQ, ALAT, NATOM, NSORT, RWS, NSORTES, TAUES, SES, MAX2SORT, ITOP, IQA, NG, &
                   ISNEW, BASNEW, NATOMNEW)
                do I = NATOM+1, NATOMNEW
                    RWS(I) = SES(ISNEW(I)-NM)
                end do
            end if
        else
            do IATOM = 1, NATOM
                ISNEW(IATOM) = IMQ(IATOM)
                do I = 1, 3
                    BASNEW(I, IATOM) = BAS(I, IATOM)
                end do
            end do
        end if
        close(16)
    end subroutine perform_calculations

    subroutine export_structure(basis, radii)
        use rsutil_config, only: empty_sphere_backend
        implicit none

        ! Arguments
        real*8, intent(out), allocatable, dimension(:,:) :: basis
        real*8, intent(out), allocatable, dimension(:) :: radii


        ! Local variables
        integer :: I, J, ISR

        allocate(basis(3, NATOMNEW))
        allocate(radii(NATOMNEW))

        do I = 1, NATOMNEW
            basis(:, I) = BASNEW(:, I)
            if (trim(empty_sphere_backend) == "questaal") then
                radii(I) = RWS(I)
            else
                radii(I) = RWS(IMQ(I))
            end if
        end do

    end subroutine export_structure

    subroutine write_output()
        use rsutil_config, only: empty_sphere_backend
        implicit none

        ! Local variables
        integer :: I, J, ISR

        open (22, FILE='radii.out')
        write (22, *) NATOMNEW

        do I = 1, NATOMNEW
            ISR = ISNEW(I)
            if (I <= NATOM) then
                if (trim(empty_sphere_backend) == "questaal") then
                    write (22, 99001) ISNEW(I), (BASNEW(J, I), J=1, 3), RWS(I), ' / A/'
                else
                    write (22, 99001) ISNEW(I), (BASNEW(J, I), J=1, 3), RWS(IMQ(I)), ' / A/'
                end if
            else if (trim(empty_sphere_backend) == "questaal") then
                write (22, 99001) ISR, (BASNEW(J, I), J=1, 3), RWS(I), ' / E/'
            else
                write (22, 99001) ISR, (BASNEW(J, I), J=1, 3), SES(ISR - NM), ' / E/'
            end if
        end do
        close (22)

99001 format(i4, 2X, 3F21.15, 3X, f21.15, 3X, a)
    end subroutine write_output

    subroutine update_structure(lat_input, atm_input, lat_data)
        implicit none

        ! Arguments
        type(lattice_inp), intent(inout) :: lat_input
        type(atoms_inp), intent(inout) :: atm_input
        type(lattice_t), intent(inout) :: lat_data
        !real*8, intent(in) :: crd(3, nbulk_bulk)
        !integer, intent(in) :: nbulk_bulk

        ! Local variables
        integer :: I, J
        ! real*8, allocatable :: temp_basis(:,:)
        ! real*8, allocatable :: temp_ws(:)
        ! integer :: temp_nbulk_bulk
        ! character(len=8), allocatable :: temp_label(:)
        ! integer, allocatable :: temp_atomic_number(:)


        ! Update structure
        call parse_structure(lat_input, atm_input, lat_data)

        call perform_calculations()

        call transfer_to_struct(lat_input, atm_input, lat_data)

        call write_output()

    end subroutine update_structure

    subroutine transfer_to_struct(lat_input, atm_input, lat_data)
        implicit none

        ! Arguments
        type(lattice_inp), intent(inout) :: lat_input
        type(atoms_inp), intent(inout) :: atm_input
        type(lattice_t), intent(inout) :: lat_data

        ! Local variables
        integer :: I
        real*8, allocatable :: temp_basis(:,:)
        real*8, allocatable :: temp_ws(:)
        integer :: temp_nbulk_bulk
        character(len=8), allocatable :: temp_label(:)
        integer, allocatable :: temp_atomic_number(:)

        ! Check if structure has grown, then update data
        !if (NATOMNEW > lat_data%nbulk_bulk) then
            ! Allocate temporary arrays
            allocate(temp_basis(3, lat_data%nbulk_bulk))
            allocate(temp_label(lat_data%nbulk_bulk))
            allocate(temp_atomic_number(lat_data%nbulk_bulk))
            allocate(temp_ws(lat_data%nbulk_bulk))

            ! Copy data from lat_data to temporary arrays
            temp_nbulk_bulk = lat_data%nbulk_bulk
            do I = 1, temp_nbulk_bulk
                temp_basis(:, I) = lat_data%crd(:, I)
                temp_label(I) = atm_input%label(I)
                temp_atomic_number(I) = atm_input%atomic_number(I)
                temp_ws(I) = lat_data%ws(I)
            end do

            ! Deallocate existing arrays
            deallocate(lat_data%crd, atm_input%label, atm_input%atomic_number, lat_data%ws)

            ! Expand dimensions of arrays
            allocate(lat_data%crd(3, NATOMNEW))
            allocate(atm_input%label(NATOMNEW))
            allocate(atm_input%atomic_number(NATOMNEW))
            allocate(lat_data%ws(NATOMNEW))

            ! Fill arrays with temporary data
            do I = 1, temp_nbulk_bulk
                lat_data%crd(:, I) = temp_basis(:, I)
                atm_input%label(I) = temp_label(I)
                atm_input%atomic_number(I) = temp_atomic_number(I)
                lat_data%ws(I) = RWS(I)
            end do

            ! print *,' Test label :',atm_input%label(temp_nbulk_bulk)
            ! Fill arrays with new data
            do I = temp_nbulk_bulk + 1, NATOMNEW
                lat_data%crd(:, I) = BASNEW(:, I)
                write(atm_input%label(I), '(A,I0)') 'E', I-temp_nbulk_bulk
                atm_input%atomic_number(I) = 0
                lat_data%ws(I) = RWS(I)
            end do
            lat_data%nbulk_bulk = NATOMNEW
         
        !else
        !    ! No new atoms: Update RWS only
        !    do I = 1, lat_data%nbulk_bulk
        !        lat_data%ws(I) = RWS(I)
        !    end do
        !end if

        lat_input%ntype = NATOMNEW
        lat_input%wav = sum(lat_data%ws(1:NATOMNEW)) / NATOMNEW / ANGSTROM_TO_BOHR

    end subroutine transfer_to_struct

end module spheres
