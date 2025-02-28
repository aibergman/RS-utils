subroutine radii
    !*==aa0001.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    integer NTMAX, NRMAX, NRAD1, MAXDIM, MAX2SORT
    parameter(NTMAX=200, NRMAX=300, NRAD1=251, MAXDIM=1000000, MAX2SORT=NTMAX)
    !
    ! COMMON variables
    !
    real*8 AV(3), BV(3), CV(3), QLAT(3, 3), RMAXES, &
        RMINES
    integer N1ES, N2ES, N3ES, NATES
    common / EMPTYS / NATES, N1ES, N2ES, N3ES, RMINES, RMAXES
    common / L2LAT / AV, BV, CV, QLAT
    !
    ! Local variables
    !
    real*8 ALAT, AMTC(NRMAX, NTMAX), BAS(3, NTMAX), &
        BASNEW(3, MAX2SORT), CONC(NTMAX), DPAS, R0(NTMAX), R0ACT, &
        R0SITE(MAX2SORT), RAD, RAT(NRAD1, MAX2SORT), RHOSITE(NRAD1, &
        MAX2SORT), RINT, RWS(NTMAX), SES(MAX2SORT), TAUES(3, MAX2SORT), &
        WS, WSREST(MAX2SORT), Z(NTMAX), ZN, ZZ(MAX2SORT)
    logical DB
    real*8 DNEVMOD
    character(100) FFF
    integer I, IATOM, IDUM, IMQ(NTMAX), IMT(NTMAX), INI, IPNT, IQA( &
        48), ISNEW(MAX2SORT), ISORT, ISP1, ISP2, ISR, IST, IT, ITOP(48), &
        ITYPE, J, NATOM, NATOMNEW, NDNEV, NEED_ES, NG, NM, NQ, NRADMAX, &
        NSORT, NSORTES, NT, W(MAXDIM)
    character(4) TXTT(NTMAX)
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    ! ==================================================
    !
    !  parameters which can be read from input
    !
    ! - basic translation vectors
    R0ACT = 1.D-2
    NRADMAX = NRMAX - 1
    !exp. pass
    DPAS = 0.05D0
    !interpolation
    NDNEV = 14
    !
    !debug
    DB = .true.
    !cell division
    N1ES = 24
    N2ES = 24
    N3ES = 24
    !
    !      rmines=1.2D0              !minimal ES
    !      rmaxes=2.5D0              !maximal ES
    !      need_es=0                 !0 - activate search of es
    !
    ! ==================================================
    !
    !       READING INPUT DATA:
    !
    INI = 23
    call GETARG(1, FFF)
    if (FFF /= ' ') then
        open (INI, FILE=FFF, STATUS='old')
    else
        INI = 5
        print*,'READING FROM THE STANDART INPUT:'
    end if
    !----------------------------------------------------
    read (INI, *) NEED_ES
    read (INI, *) RMINES
    read (INI, *) RMAXES
    !----------------------------------------------------
    do I = 1, 9
        read (INI, *)
    end do
    read (INI, *) ISP1, ISP2
    if (DB) write(16, *) 'ISP1, ISP2: ', ISP1, ISP2
    read (INI, *)
    read (INI, *)
    read (INI, *)
    read (INI, *) ALAT
    if (DB) write(16, *) 'ALAT:', ALAT
    read (INI, *)
    read (INI, *)
    read (INI, *)
    read (INI, *)
    read (INI, *)
    read (INI, *)
    read (INI, *)
    read (INI, *) (AV(I), I=1, 3)
    if (DB) write(16, *) 'AV:',(AV(I), I=1, 3)
    read (INI, *) (BV(I), I=1, 3)
    read (INI, *) (CV(I), I=1, 3)
    if (DB) print *, 'BV:',(BV(I), I=1, 3)
    if (DB) print *, 'CV:',(CV(I), I=1, 3)
    read (INI, *)
    read (INI, *) NQ
    if (DB) print *,'N*:', NQ
    read (INI, *)
    if (DB) print *,'NATOM:', NQ
    do I = 1, NQ
        read (INI, *) IDUM, IMQ(I), (BAS(J, I), J=1, 3)
        if (DB) print *, I, "BAS", BAS(1, I), BAS(2, I), BAS(3, I), "IMQ", IMQ(I)
        !,RWS(IMQ(i)),nl,
    end do
    read (INI, *)
    read (INI, *) NM
    if (DB) print *, 'NM:', NM
    read (INI, *)
    do I = 1, NM
        read (INI, *)
    end do
    read (INI, *)
    read (INI, *) NSORT
    if (DB) print *,'NSORT:', NSORT
    read (INI, *)
    NT = NSORT
    do I = 1, NSORT
        read (INI, *) IDUM, Z(I), TXTT(I), IDUM, CONC(I), IMT(I)
        if (DB) print *, I, "Z,TXTT,CONC,IMT:",  Z(I), TXTT(I), CONC(I), IMT(I)
    end do
    ! Remove symmetry analysis for now
    NG = 1
    ITOP(NG) = 1
    IQA(NG) = 1
    ! read (INI, *) NG
    ! read (INI, *) (ITOP(I), I=1, NG)
    ! if (DB) print *,'ITOP:', (ITOP(I), I=1, NG)
    ! read (INI, *) (IQA(I), I=1, NG)
    ! if (DB) print *,'IQA:', (IQA(I), I=1, NG)
    if (DB) then
        print*,'PARAMETERS:', NT, NQ, NM, 'NT,NQ,NM,-- ALAT:', &
            ALAT
        print*,'Z(i),TXTT(i),CONC(i),IMT(i)'
        do I = 1, NT
            print*,Z(I), TXTT(I), CONC(I), IMT(I)
        end do
        print*,'qx(i),qy(i),qz(i),IMQ(i)'
        do I = 1, NQ
            print*,(BAS(J, I), J=1, 3), IMQ(I)
        end do
        !
        print*,'Lattice vectors:'
        print*,AV
        print*,BV
        print*,CV
    end if
    !
    NG = 1
    !
    !     ==================================================
    !
    !
    call WKINIT(MAXDIM)
    !
    call GETGBASIS
    !
    !
    ! amtc,dpas,r0,rom,rop,r0,rop,rom,      dummy initialisation
    !
    R0 = 9999D0
    AMTC = 9999D0
    RHOSITE(:, :) = 9999D0
    !
    call READRO(NT, AMTC, Z, R0, TXTT, RHOSITE, NRADMAX, NRAD1)
    !
    ! Producing new radial mesh:
    !
    do IT = 1, NT
        do I = 1, NRAD1
            RAT(I, IT) = R0(IT) * exp(dble(I - 1) * DPAS)
        end do
    end do
    !
    do I = 1, NM
        WSREST(I) = 0
        ZZ(I) = 0
        R0SITE(I) = R0ACT
    end do
    !
    RHOSITE(:, :) = 0.0D0
    !
    do ITYPE = 1, NT
        !
        ! this is cycle over atomic types
        ! IMT(itype) - number of inequivalent lattice position, which is
        ! occupied by this atom
        !
        ISORT = IMQ(IMT(ITYPE))
        ZN = Z(ITYPE)
        call DEFWSR(WS, ZN)
        if (DB) print*,'type:', ITYPE, ' Z=', ZN, ' Rad=', WS, ' conc=&
           &', CONC(ITYPE)
        WSREST(ISORT) = WSREST(ISORT) + WS * CONC(ITYPE)
        ZZ(ISORT) = ZZ(ISORT) + ZN * CONC(ITYPE)
        do IPNT = 1, NRAD1
            RAD = R0ACT * exp(dble(IPNT - 1) * DPAS)
            if (RAD < RAT(NRAD1 - 2, ITYPE)) then
                RINT = max(0.D0, DNEVMOD(RAT(1, ITYPE), AMTC(1, ITYPE), &
                                         RAD, NRAD1, NDNEV))
            else
                RINT = 0
            end if
            RHOSITE(IPNT, ISORT) = RHOSITE(IPNT, ISORT) + CONC(ITYPE) &
                                   * RINT
        end do
    end do
    !
    do I = 1, NM
        if (DB) print*,' Site:', I, '  RaD:', WSREST(I), ' Z:', ZZ(&
           &I)
    end do
    !
    !      stop
    !
    NSORT = NM
    NATOM = NQ
    !
    !
    if (DB) then
        print*,'av=', AV
        print*,'bv=', BV
        print*,'cv=', CV
        print*,'imq', (IMQ(I), I=1, NATOM)
        print*,'BAS'
        do I = 1, NATOM
            print*,(BAS(J, I), J=1, 3)
        end do
        do IST = 1, NSORT
            print*,'R0S', R0SITE(IST), ' ZZ', ZZ(IST), WSREST(IST)
        end do
        print*,'dpas,alat,natom,nsort:', DPAS, ALAT, NATOM, NSORT
    end if
    !
    call DEFMTR(W, NRAD1, AV, BV, CV, IMQ, BAS, RHOSITE, R0SITE, &
                DPAS, ALAT, NATOM, NSORT, ZZ, RWS, NEED_ES, NRAD1, WSREST)
    !
    NATOMNEW = NATOM
    if (NEED_ES == 0) then
        call EMPTYSPHERES(AV, BV, CV, BAS, IMQ, ALAT, NATOM, NSORT, &
                          RWS, NSORTES, TAUES, SES, MAX2SORT, ITOP, IQA, NG, ISNEW, &
                          BASNEW, NATOMNEW)
        !ccccccccccccccccccccccccccccccccccccccccccccccccccc
    else
        do IATOM = 1, NATOM
            ISNEW(IATOM) = IMQ(IATOM)
            do I = 1, 3
                BASNEW(I, IATOM) = BAS(I, IATOM)
            end do
        end do
        !cccccccccccccccccccccccccccccccccccccccccccccccccc
    end if
    !
    !
    ! ==================================================
    !
    !       WRITING OUTPUT DATA:
    !
    !
    open (22, FILE='radii.out')
    write (22, *) NATOMNEW
    !
    do I = 1, NATOMNEW
        ISR = ISNEW(I)
        if (I <= NATOM) then
            write (22, 99001) ISNEW(I), (BASNEW(J, I), J=1, 3), &
                RWS(IMQ(I)), ' / A/'
        else
            write (22, 99001) ISR, (BASNEW(J, I), J=1, 3), SES(ISR - &
                                                               NM), ' / E/'
        end if
    end do
    close (22)
99001 format(i4, 2X, 3F21.15, 3X, f21.15, 3X, a)
end subroutine radii

!*==readro.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
subroutine READRO(NSORT, RO, Z, R0, TXT, ROM, NRADMAX, NRAD1)
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    integer NRAD1, NRADMAX, NSORT
    real*8 R0(*), RO(NRADMAX + 1, *), ROM(NRAD1, &
                                                   *), Z(*)
    character(4) TXT(*)
    !
    ! Local variables
    !
    integer I, ISORT, ISR, IZC, IZN, IZNUC
    !
    !*** End of declarations rewritten by SPAG
    !
    open (8, STATUS='scratch')
    do ISORT = 1, NSORT
        if (Z(ISORT) < 0) then
            IZNUC = -10
            write(16, *) 'CALING QZZZZZZC'
            call QZZC(TXT(ISORT), IZNUC, IZC)
            Z(ISORT) = IZNUC
        end if
        if (Z(ISORT) < 0.3D0) then
            ! emty sphere
            do I = 1, NRAD1
                RO(I, ISORT) = 0.D0
                ROM(I, ISORT) = 0.D0
                R0(ISORT) = 1D-5
            end do
        else
            do ISR = 1, ISORT - 1
                if (nint(Z(ISR)) == nint(Z(ISORT))) then
                    do I = 1, NRAD1
                        RO(I, ISORT) = RO(I, ISR)
                        ROM(I, ISORT) = ROM(I, ISR)
                    end do
                    R0(ISORT) = R0(ISR)
                    goto 100
                end if
            end do
            !          write(buf,2)nint(z(isort)),nint(zc)
            write (16, '(/'' For atom '',a4,''    Z='',f4.0/)') TXT(ISORT), Z&
               & (ISORT)
            IZN = nint(Z(ISORT))
            IZC = 0
            call RHFDS(IZN, IZC, RO(1, ISORT), ROM(1, ISORT), R0( &
                       ISORT))
            !
        end if
100 end do
    close (8)
end subroutine READRO
!*==qzzc.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
subroutine QZZC(TXT, IZ, IZC)
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    integer IZ, IZC
    character(2) TXT
    !
    ! Local variables
    !
    character(210) ELEMENTS
    integer I, IZCOR(105)
    character(2) TX
    !
    !*** End of declarations rewritten by SPAG
    !
    data IZCOR/0, 0, 0, 6*2, 4, 4, 6*10, 12, 12, 12*18, 4* &
        28, 30, 30, 12*36, 4*46, 48, 48, 16*54, 10*68, 4*78, 80, &
        80, 18*86/
    ELEMENTS = 'E H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFe&
       &CoNi'//'CuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLa&
       &'//'CePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn'//'&
       &FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLwKu'
    if (IZ >= 0) then
        IZC = IZCOR(IZ + 1)
        return
    end if
    TX = TXT
    if (TX(2:2) == '_') TX(2:2) = ' '
    do I = 1, 105
        !        print 1,elements(i*2-1:i*2),i-1,izcor(i)
        if (ELEMENTS(I * 2 - 1:I * 2) == TXT(1:2)) then
            IZ = I - 1
            IZC = IZCOR(I)
            return
        end if
    end do
    write (16, *) TXT
    stop 'Error : Such Chemical element is ABSENT'
end subroutine QZZC
!*==shortn.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
subroutine SHORTN(P, P1)
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    real*8 EPS
    parameter(EPS=1.D-7)
    !
    ! COMMON variables
    !
    real*8 PLAT(3, 3), QLAT(3, 3)
    common / L2LAT / PLAT, QLAT
    !
    ! Dummy arguments
    !
    real*8 P(3), P1(3)
    !
    ! Local variables
    !
    integer I
    real*8 X(3)
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    do I = 1, 3
        X(I) = P(1) * QLAT(1, I) + P(2) * QLAT(2, I) + P(3) * QLAT( &
               3, I)
        X(I) = X(I) - nint(X(I) - EPS)
        !$$$        if(x(i).lt.-0.5d0+eps)x(i)=x(i)+1.d0
    end do
    do I = 1, 3
        P1(I) = X(1) * PLAT(I, 1) + X(2) * PLAT(I, 2) + X(3) &
                * PLAT(I, 3)
    end do
end subroutine SHORTN
!*==vec.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
subroutine VEC(V0, V1, V2, G1)
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    real*8 G1(3, 3), V0(3), V1(3), V2(3)
    !
    ! Local variables
    !
    integer I, J
    !
    !*** End of declarations rewritten by SPAG
    !
    do I = 1, 3
        V0(I) = V1(I)
        do J = 1, 3
            V0(I) = V0(I) + V2(J) * G1(J, I)
        end do
    end do
    call SHORTN(V0, V0)
end subroutine VEC
!*==getgbasis.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
subroutine GETGBASIS
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    real*8 AVC(3), BVC(3), CVC(3), FVC(3), HVC(3), &
        PVC(3)
    common / L2LAT / AVC, BVC, CVC, FVC, PVC, HVC
    !
    ! Local variables
    !
    integer I
    real*8 VOLCELL, W1
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    call CROSS(FVC, BVC, CVC)
    call CROSS(PVC, CVC, AVC)
    call CROSS(HVC, AVC, BVC)
    W1 = AVC(1) * FVC(1) + AVC(2) * FVC(2) + AVC(3) * FVC(3)
    VOLCELL = abs(W1)
    !
    if (VOLCELL < 1.D-5) then
        print*,AVC
        print*,BVC
        print*,CVC
        call ENDJOB(10, 'GETGBASIS: Lattice vectors are complanar!')
    end if
    do I = 1, 3
        FVC(I) = FVC(I) / W1
        PVC(I) = PVC(I) / W1
        HVC(I) = HVC(I) / W1
    end do
end subroutine GETGBASIS
