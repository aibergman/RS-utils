!*==emptyspheres.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
subroutine EMPTYSPHERES(AVC, BVC, CVC, BAS, IS, ALAT, NATOM, &
                        NSORT, S, NSORTES, TAUES, SES, MAX2SORT, ITOP, IQA, NG, ISNEW, &
                        BASNEW, NATOMNEW)
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    integer NPATB1
    parameter(NPATB1=6000)
    !
    ! Dummy arguments
    !
    real(8) ALAT
    integer MAX2SORT, NATOM, NATOMNEW, NG, NSORT, NSORTES
    real(8) AVC(3), BAS(3, *), BASNEW(*), BVC(3), CVC(3), &
        S(*), SES(*), TAUES(3, *)
    integer IQA(*), IS(*), ISNEW(*), ITOP(*)
    !
    ! Local variables
    !
    real(8) ANG(3), GEN(3, 3, 48), S2NEW(:), SNEW(:), TAU(:, :), &
        VC(3, 48), VV(3)
    integer I, IER, IG, IGA, IGC, ISB(:), ISKIP(:), I_ISB, I_ISKIP, &
        I_NHSORTES, I_S2NEW, I_SNEW, I_TAU, J, MAXNPAT, MAXNPATG, &
        NHSORTES(:)
    !
    !*** End of declarations rewritten by SPAG
    !
    ALLOCATABLESNEW, S2NEW, TAU, NHSORTES, ISKIP, ISB
    !
    MAXNPAT = MAX2SORT
    MAXNPATG = NPATB1
    !
    do I = 1, 3
        VV(I) = -BAS(I, 1)
    end do
    !
    do IGC = 1, NG
        IG = ITOP(IGC)
        IGA = IG
        if (IGA > 32) IGA = IGA - 32
        call EILANG(ANG, IGA)
        call TURNM(ANG, GEN(1, 1, IGC))
        !                              !DEF TURN MATRIX
        if (IG /= IGA) then
            !                              !ADD INVERSION IF NEED
            do I = 1, 3
                do J = 1, 3
                    GEN(I, J, IGC) = -GEN(I, J, IGC)
                end do
            end do
        end if
        call VEC(VC(1, IGC), BAS(1, IQA(IGC)), VV, GEN(1, 1, IGC))
                print*,Igc,ig,(vc(i,igc),i=1,3)
    end do
    !
    do I = 1, NG
        IG = ITOP(I)
        !
    end do
    !
    call DEFRR(I_SNEW, MAXNPAT)
    call DEFRR(I_S2NEW, MAXNPAT)
    call DEFI(I_ISB, MAXNPATG)
    call DEFRR(I_TAU, MAXNPATG * 3)
    call DEFI(I_NHSORTES, MAXNPAT)
    call DEFI(I_ISKIP, MAXNPAT)
    !
    allocate (SNEW(MAXNPAT), S2NEW(MAXNPAT), TAU(3, MAXNPATG))
    allocate (NHSORTES(MAXNPAT), ISKIP(MAXNPAT), ISB(MAXNPATG))
    !
    call EMPTY(IER, AVC, BVC, CVC, BAS, IS, ALAT, NATOM, NSORT, S, &
               BASNEW, ISNEW, SNEW, S2NEW, ISB, TAU, MAXNPAT, MAXNPATG, NSORTES, &
               TAUES, SES, NHSORTES, ISKIP, GEN, VC, NG, NATOMNEW)
    !
    !     nsortes=iter-1
    write (16, *) ' IRR. POSITIONS:'
    do I = 1, NSORTES
        write (16, '(3f21.15)') (TAUES(J, I), J=1, 3)
    end do
    write (16, *) 'Spheres done!'
    !      stop
end subroutine EMPTYSPHERES
!*==empty.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
subroutine EMPTY(IER, AVC, BVC, CVC, BAS, IS, ALAT, NATOM, NSORT, &
                 S, BASNEW, ISNEW, SNEW, S2NEW, ISB, TAU, MAXNPAT, MAXNPATG, &
                 NSORTES, TAUES, SES, NHSORTES, ISKIP, GM, VC, NG, NATOMNEW)
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    real(8) BOA, COA, RMAXES, RMINES, U(3, 3), UM1(3, 3)
    integer IDATASTR, N1ES, N2ES, N3ES, NATES
    common / EMPTYS / NATES, N1ES, N2ES, N3ES, RMINES, RMAXES
    common / TAU2 / BOA, COA
    common / U / U, UM1, IDATASTR
    !
    ! Dummy arguments
    !
    real(8) ALAT
    integer IER, MAXNPAT, MAXNPATG, NATOM, NATOMNEW, NG, NSORT, &
        NSORTES
    real(8) AVC(3), BAS(3, *), BASNEW(3, *), BVC(3), CVC(3), &
        GM(3, 3, 48), S(*), S2NEW(*), SES(*), SNEW(*), &
        TAU(3, *), TAUES(3, *), VC(3, 48)
    integer IS(*), ISB(*), ISKIP(*), ISNEW(*), NHSORTES( &
        *)
    !
    ! Local variables
    !
    real(8) ANG(3), CS, CVOL, DEMIN, EPS_ANG, EPS_ANGMAX, PA, PB, PC, &
        PI43, RAD, RADNEW, SM, SMN, SSM, SVOLA, SVOLE, VE(3, 48), &
        VSHIFT(3), VV(3)
    real(8) DPI, TRNT
    integer I, IATOM, IG, IIN, IN, INVMAX, IODN, ISORT, ITAV, ITER, &
        IUN, J, K, NATB, NATS, NSORTNEW
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
    do I = 1, 3
        do J = 1, 3
            U(I, J) = 0
            UM1(I, J) = 0
            GM(I, J, 1) = 0
        end do
        U(I, I) = 1
        UM1(I, I) = 1
        GM(I, I, 1) = 1
        VSHIFT(I) = 0
        VC(I, 1) = 0
    end do
    !      it(1)=1
    !     ng=1
    BOA = 1
    COA = 1
    !---------------------------------------
    IUN = 16
    NATES = 0
    SSM = 0
    do I = 1, NSORT
        SNEW(I) = S(I) / ALAT
        SSM = max(SSM, SNEW(I))
        S2NEW(I) = SNEW(I)**2
    end do
    INVMAX = 0
    if (RMAXES < 0) INVMAX = 1
    if (RMINES < 0) INVMAX = -INVMAX
    RMINES = abs(RMINES)
    RMAXES = abs(RMAXES)
    SM = RMAXES / ALAT + SSM
    SMN = RMINES / ALAT
    do IATOM = 1, NATOM
        ISNEW(IATOM) = IS(IATOM)
        do I = 1, 3
            BASNEW(I, IATOM) = BAS(I, IATOM)
        end do
    end do
    NATOMNEW = NATOM
    NSORTNEW = NSORT
    ITER = 0
    write(IUN, *) 'BOA, COA:', BOA, COA, ALAT
    write (IUN, *) '  Searching of Empty spheres  '
    write (*, *) '  Searching of Empty spheres  '
    write (IUN, '(a,3i3)') ' Mesh: ', N1ES, N2ES, N3ES
    write (IUN, 99001) 'Smin, Smax (a.u.): ', RMINES, RMAXES
100 continue
    ITER = ITER + 1
    write (IUN, *) '  Iteration ', ITER
    !
    call GETEPOS(MAXNPATG, AVC, BVC, CVC, BASNEW, NATOMNEW, ISNEW, &
                 TAU, NATB, ISB, SNEW, S2NEW, N1ES, N2ES, N3ES, PA, PB, PC, RAD, &
                 IER, SM)
    !
    write (IUN, 99001) ' Maximum possible sphere is ', RAD * ALAT
    write (IUN, 99001) ' Position in a,b,c:', PA, PB, PC
    if (RAD < SMN) then
        write (IUN, *) 'skipped'
        write(IUN, *) 'RAD=', RAD * ALAT, ' RMINES=', RMINES
        goto 300
    end if
    !
    if (IER == 0) then
        ISKIP(ITER) = 0
        if (RAD > RMAXES / ALAT) then
            if (INVMAX == 0) then
                write (IUN, 99001) 'Contracted to be ', RMAXES
                ! rad is used later
                RAD = RMAXES / ALAT
            else
                ISKIP(ITER) = INVMAX
                write (IUN, 99001) 'Sphere is too large', RAD * ALAT
                write (IUN, 99001) 'Small auxiliary sphere is inserted instead'
            end if
        end if
        !
        ! multiplication :
        VE(1, 1) = PA * AVC(1) + PB * BVC(1) + PC * CVC(1) &
                   - VSHIFT(1)
        VE(2, 1) = PA * AVC(2) + PB * BVC(2) + PC * CVC(2) &
                   - VSHIFT(2)
        VE(3, 1) = PA * AVC(3) + PB * BVC(3) + PC * CVC(3) &
                   - VSHIFT(3)
        !
        ITAV = 0
150     continue
        ITAV = ITAV + 1
        IN = 1
        EPS_ANGMAX = 0.D0
        call SHORTN(VE, VE)
        !       print*,(ve(i,1),i=1,3)
        do IG = 1, NG
            call VEC(VV, VC(1, IG), VE, GM(1, 1, IG))
            !         print*,vv
            do IIN = 1, IN
                ANG(1) = VE(1, IIN) - VV(1)
                ANG(2) = VE(2, IIN) - VV(2)
                ANG(3) = VE(3, IIN) - VV(3)
                call SHORTN(ANG, ANG)
                EPS_ANG = ANG(1)**2 + ANG(2)**2 + ANG(3)**2
                if (EPS_ANG < 1.D-5) then
                    EPS_ANGMAX = max(EPS_ANG, EPS_ANGMAX)
                    !              print*,iin,eps_ang
                    goto 200
                end if
            end do
            IN = IN + 1
            VE(1, IN) = VV(1)
            VE(2, IN) = VV(2)
            VE(3, IN) = VV(3)
            !         print*,'in=',in,vv
200     end do
        VV(1) = 0
        VV(2) = 0
        VV(3) = 0
        IODN = 0
        DEMIN = 1.D10
        do IIN = 2, IN
            ANG(1) = VE(1, IIN) - VE(1, 1)
            ANG(2) = VE(2, IIN) - VE(2, 1)
            ANG(3) = VE(3, IIN) - VE(3, 1)
            call SHORTN(ANG, ANG)
            EPS_ANG = ANG(1)**2 + ANG(2)**2 + ANG(3)**2
            EPS_ANG = sqrt(EPS_ANG)
            if (DEMIN > EPS_ANG) DEMIN = EPS_ANG
            if (EPS_ANG < RAD * 1.5D0) then
                VV(1) = VV(1) + ANG(1)
                VV(2) = VV(2) + ANG(2)
                VV(3) = VV(3) + ANG(3)
                IODN = IODN + 1
            end if
            ! iin
        end do
        if (IODN /= 0) then
            write (IUN, '(a,i3,2f21.15)') 'a bit shifted position', &
                IODN, real(DEMIN), real(RAD)
            write (IUN, 99001) 'OLD:', VE(1, 1), VE(2, 1) / BOA, &
                VE(3, 1) / COA
            VE(1, 1) = VE(1, 1) + VV(1) / (IODN + 1)
            VE(2, 1) = VE(2, 1) + VV(2) / (IODN + 1)
            VE(3, 1) = VE(3, 1) + VV(3) / (IODN + 1)
            write (IUN, 99001) 'new:', VE(1, 1), VE(2, 1) / BOA, &
                VE(3, 1) / COA
            call REGET_RE(TAU, ISB, SNEW, NATB, VSHIFT, VE(1, 1), &
                          RADNEW)
            write (IUN, 99001) 'rad,radtst', RAD * ALAT, RADNEW * &
                ALAT
            RAD = RADNEW
            ! to avoide infinite cycle
            if (ITAV <= 5) goto 150
            write (IUN, *) 'can not find averaged position'
            stop
        elseif (DEMIN < 2.D0 * RAD) then
            ! overlapping spheres
            write (IUN, *) 'radius is decreased from', RAD, ' to', &
                DEMIN / 2.D0
            ! decrease radius to touching
            RAD = DEMIN / 2.D0
            if (RAD < SMN) then
                write (IUN, *) 'sphere is less then rmines'
                write (IUN, *) 'rad=', RAD * ALAT, ' rmines=', &
                    RMINES
                ! will be skipped later
                ISKIP(ITER) = -1
            end if
        end if
        write (IUN, '(a,i4,a)') 'it will be ', IN, ' ES per unit cell:'
        if (NATOMNEW + IN > MAXNPAT) then
            print*,'EMPTY: too many atoms', NATOMNEW + IN, ' maxnpat=', MAX&
               &NPAT
            print*,'try to increase array sizes in empty.f'
            stop
        end if
        do IATOM = NATOMNEW + 1, NATOMNEW + IN
            ISNEW(IATOM) = NSORTNEW + 1
            do I = 1, 3
                BASNEW(I, IATOM) = VE(I, IATOM - NATOMNEW) + VSHIFT(I)
            end do
            write (IUN, '(2x,3f15.8)') (BASNEW(I, IATOM), I=1, 3)
        end do
        !
        do K = 1, 3
            VV(K) = 0
            do J = 1, 3
                VV(K) = VV(K) + VE(J, 1) * UM1(J, K)
            end do
            TAUES(K, ITER) = VV(K)
        end do
        !      write(iun,'(i3,3f15.8)')iter,vv
        NSORTNEW = NSORTNEW + 1
        NATOMNEW = NATOMNEW + IN
        NHSORTES(ITER) = IN
        if (ISKIP(ITER) == 0) then
            SNEW(NSORTNEW) = min(RAD, RMAXES / ALAT)
        else
            SNEW(NSORTNEW) = RMINES / ALAT
        end if
        S2NEW(NSORTNEW) = SNEW(NSORTNEW)**2
        SES(ITER) = SNEW(NSORTNEW) * ALAT
        write(IUN, *) IN, ' Empty spheres added, S=', SNEW(NSORTNEW) * ALAT
        NATES = NATES + IN
        goto 100
    end if
    !
300 continue
    NSORTES = ITER - 1
    if (INVMAX >= 0) then
        ! check if some spheres should be skipped due to small radius
        do ISORT = 1, NSORTES
            if (ISKIP(ISORT) < 0) INVMAX = ISKIP(ISORT)
        end do
    end if
    if (INVMAX < 0) then
        ! skip auxiliary spheres
        !$$$        write(iun,*)'nsortes,natomnew',nsortes,natomnew
        !$$$        do isort=1,nsortes
        !$$$          write(iun,'(i3,3f21.15,2i4,f21.15)')isort,(taues(k,isort),
        !$$$     $         ,nhsortes(isort),iskip(isort),ses(isort)
        !$$$        enddo
        !$$$        do iatom=1,natomnew
        !$$$          write(iun,'(i3,3f21.15)')iatom,(basnew(k,iatom),k=1,3)
        !$$$        enddo
        !$$$        write(iun,*)'rmines',rmines
        IATOM = 1
        ISORT = 1
        do while (ISORT <= NSORTES)
            write (IUN, *) 'isort,nsortes,s', ISORT, NSORTES, SES(ISORT)
            NATS = NHSORTES(ISORT)
            if (ISKIP(ISORT) < 0) then
                write (IUN, 99001) 'sphere is skipped', (TAUES(K, ISORT) &
                                                         , K=1, 3)
                do I = ISORT + 1, NSORTES
                    NHSORTES(I - 1) = NHSORTES(I)
                    SES(I - 1) = SES(I)
                    ISKIP(I - 1) = ISKIP(I)
                    do K = 1, 3
                        TAUES(K, I - 1) = TAUES(K, I)
                    end do
                end do
                do I = IATOM + NATS, NATOMNEW
                    do K = 1, 3
                        BASNEW(K, I - NATS) = BASNEW(K, I)
                    end do
                end do
                NSORTES = NSORTES - 1
                NATOMNEW = NATOMNEW - NATS
            else
                IATOM = IATOM + NATS
                ISORT = ISORT + 1
            end if
        end do
        !$$$        write(iun,*)'nsortes,natomnew',nsortes,natomnew
        !$$$        do isort=1,nsortes
        !$$$          write(iun,'(i3,3f21.15,i4,f21.15)')isort,(taues(k,isort),k
        !$$$     $         ,nhsortes(isort),ses(isort)
        !$$$        enddo
        !$$$        do iatom=1,natomnew
        !$$$          write(iun,'(i3,3f21.15)')iatom,(basnew(k,iatom),k=1,3)
        !$$$        enddo
    end if
    !$$$      if(nsortes.gt.0)then
    !$$$        open(24,status='scratch',form='unformatted')
    !$$$        do iatom=natom+1,natomnew
    !$$$          write(24)(basnew(i,iatom),i=1,3)
    !$$$c$$$        write(iun,'(2x,3f15.8,i4)')(basnew(i,iatom),i=1,3),isnew
    !$$$        enddo
    !$$$      endif
    ! MT -> ASA spheres
    PI43 = 4.D0 * DPI() / 3.D0
    SVOLA = 0.D0
    do IATOM = 1, NATOM
        SVOLA = SVOLA + PI43 * S(IS(IATOM))**3
    end do
    SVOLE = 0.D0
    do ISORT = 1, NSORTES
        SVOLE = SVOLE + NHSORTES(ISORT) * PI43 * SES(ISORT)**3
    end do
    CVOL = TRNT(AVC, BVC, CVC) * ALAT**3
    CS = CVOL / (SVOLA + SVOLE)
    SVOLA = SVOLA * CS
    SVOLE = SVOLE * CS
    CS = abs(CS)**(1.D0 / 3.D0)
    do ISORT = 1, NSORT
        S(ISORT) = CS * S(ISORT)
    end do
    do ISORT = 1, NSORTES
        SES(ISORT) = CS * SES(ISORT)
    end do
    write (IUN, 99001) 'Vat, Ves, V', SVOLA, SVOLE, SVOLA + SVOLE
    write (IUN, 99001) 'Vat/V, Ves/V', SVOLA / CVOL, SVOLE / CVOL
99001 format(1X, a, 3F21.15)
end subroutine EMPTY
!*==getepos.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
subroutine GETEPOS(NPAT, A, B, C, BAS, NATOM, IS, TAU, NATB, ISB, &
                   S, S2, N1, N2, N3, PAE, PBE, PCE, RAD, IER, SM)
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    integer IER, N1, N2, N3, NATB, NATOM, NPAT
    real(8) PAE, PBE, PCE, RAD, SM
    real(8) A(3), B(3), BAS(3, *), C(3), S(*), S2(*), &
        TAU(3, NPAT)
    integer IS(NPAT), ISB(NPAT)
    !
    ! Local variables
    !
    real(8) ANA, ANB, ANC, DD, DDMAX, DIAG, DIST2, DMAX, DMIN, DX, DY, &
        DZ, PA, PB, PC, X, Y, Z
    integer I, IA, IATOM, IB, IC, J, K
    !
    integer, parameter :: IUN = 16
    !*** End of declarations rewritten by SPAG
    !
    !
    DIAG = max((A(1) + B(1) + C(1))**2 + (A(2) + B(2) + C(2))**2 &
               + (A(3) + B(3) + C(3))**2, (A(1) - B(1) + C(1))**2 &
               + (A(2) - B(2) + C(2))**2 + (A(3) - B(3) + C(3))**2, &
               (A(1) + B(1) - C(1))**2 + (A(2) + B(2) - C(2))**2 &
               + (A(3) + B(3) - C(3))**2, (A(1) - B(1) - C(1))**2 &
               + (A(2) - B(2) - C(2))**2 + (A(3) - B(3) - C(3))**2)
    !
    IER = 0
    DIAG = sqrt(DIAG) / 2
    DDMAX = (SM + DIAG)**2
    !
    NATB = 0
    do I = -3, 3
        do J = -3, 3
            do K = -3, 3
                DX = I * A(1) + J * B(1) + K * C(1)
                DY = I * A(2) + J * B(2) + K * C(2)
                DZ = I * A(3) + J * B(3) + K * C(3)
                do IATOM = 1, NATOM
                    X = DX + BAS(1, IATOM)
                    Y = DY + BAS(2, IATOM)
                    Z = DZ + BAS(3, IATOM)
                    !              print*,X,Y,Z
                    DD = X**2 + Y**2 + Z**2
                    if (DD < DDMAX) then
                        NATB = NATB + 1
                        TAU(1, NATB) = X
                        TAU(2, NATB) = Y
                        TAU(3, NATB) = Z
                        if (NATB > NPAT) then
                            IER = 1
                            write(IUN, *) 'IER=', IER
                            write(IUN, *) 'GETEPOS: too many atoms are generated ', NATB, ' npat=', NPAT
                            write(IUN, *) 'try to increase npatb1 in empty.f'
                            close(IUN)
                            return
                        end if
                        ISB(NATB) = IS(IATOM)
                    end if
                    ! iatom
                end do
                ! k
            end do
            ! j
        end do
        ! i
    end do
    !
    PAE = 0
    PBE = 0
    PCE = 0
    !
    ANA = 0.5D0 / N1
    ANB = 0.5D0 / N2
    ANC = 0.5D0 / N3
    DMAX = 0
    !      npnt=0
    do IA = -N1, N1 - 1, 2
        PA = ANA * IA
        do IB = -N2, N2 - 1, 2
            PB = ANB * IB
            do IC = -N3, N3 - 1, 2
                !            npnt=npnt+1
                PC = ANC * IC
                !$$$            lg=ia.eq.-24.and.ib.eq.0.and.ic.eq.-72
                X = PA * A(1) + PB * B(1) + PC * C(1)
                Y = PA * A(2) + PB * B(2) + PC * C(2)
                Z = PA * A(3) + PB * B(3) + PC * C(3)
                !           if(lg)print*,'XYZ:',x,y,z
                DMIN = 1.D10
                do I = 1, NATB
                    DIST2 = (X - TAU(1, I))**2 + (Y - TAU(2, I))**2 + (Z - TAU( &
                                                                       3, I))**2
                    !             if(lg.and.isb(i).eq.8)print*,'TTT:',tau(1,i),tau(2,i)
                    !    $             ,tau(3,i)
                    !$$$              print *,'ia,ib,ic,i',ia,ib,ic,i,natb
                    !$$$              print *,'dist2,isb(i),sq',dist2,isb(i),s2(isb(i))
                    !$$$                PRINT*,I,DIST2,ISB(I),S2(ISB(I))
                    !$$$                print*, x,y,z,ia,ib,ic
                    !$$$                print*,tau(1,i),tau(2,i),tau(3,i)
                    if (DIST2 < S2(ISB(I))) goto 20
                    DD = sqrt(DIST2) - S(ISB(I))
                    DMIN = min(DD, DMIN)
                    ! i
                end do
                !            print*,DMIN
                if (DMIN > DMAX) then
                    DMAX = DMIN
                    PAE = PA
                    PBE = PB
                    PCE = PC
                    !$$$              iae=ia
                    !$$$              ibe=ib
                    !$$$              ice=ic
                end if
                ! ic
20          end do
            ! ib
        end do
        ! ia
    end do
    !      print*,iter,'MAX is:',DMAX,iae,ibe,ice,pae,pbe,pce
    !      print*,iter,'npnt=',npnt
    RAD = DMAX
end subroutine GETEPOS
!*==reget_re.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
subroutine REGET_RE(TAU, ISB, S, NATB, VSHIFT, VS, RAD)
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    integer NATB
    real(8) RAD
    integer ISB(*)
    real(8) S(*), TAU(3, *), VS(*), VSHIFT(*)
    !
    ! Local variables
    !
    real(8) D, DMIN, V(3)
    integer I, K
    !
    !*** End of declarations rewritten by SPAG
    !
    do K = 1, 3
        V(K) = VS(K) + VSHIFT(K)
    end do
    DMIN = 1.D10
    do I = 1, NATB
        D = (V(1) - TAU(1, I))**2 + (V(2) - TAU(2, I))**2 + &
            (V(3) - TAU(3, I))**2
        D = sqrt(D) - S(ISB(I))
        if (D < 0.D0) then
            print*,'REGET_RE: error'
            print*,I, D + S(ISB(I)), ISB(I), S(ISB(I))
            print*,(V(K), K=1, 3)
            print*,TAU(1, I), TAU(2, I), TAU(3, I)
            return
        end if
        DMIN = min(D, DMIN)
        ! natb
    end do
    RAD = DMIN
end subroutine REGET_RE
