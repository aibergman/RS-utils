!*==defmtr.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine DEFMTR(W, NRAD1, AV, BV, CV, IS, BAS, RO, R0, DPAS, &
                        ALAT, NATOM, NSORT, Z, SMT, IMT, NRADMAX, WSREST)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 ALAT, DPAS
          integer IMT, NATOM, NRAD1, NRADMAX, NSORT
          real*8 AV(3), BAS(3, *), BV(3), CV(3), R0(*), &
              RO(NRADMAX, *), SMT(*), WSREST(*), Z(*)
          integer IS(*), W(*)
!
! Local variables
!
          real*8 A(:), AVW, BOUND(3, 3), DV(:), ORIGIN(3), &
              PI, PLAT(3, 3), POT(:, :), R(:), VOLCEL
          real*8 DPI
          integer I, J, LTMAX(3), NDEL(3), NR(:), NRMAX, OA, OBAS, ODV, &
              ONR, OPOT, OR
          save A, AVW, BOUND, DV, I, J, LTMAX, NDEL, NRMAX, OA, OBAS, ODV, &
              ONR, OPOT, OR, ORIGIN, PI, PLAT, POT, R, VOLCEL
!
!*** End of declarations rewritten by SPAG
!
          ALLOCATABLEPOT, DV, R, A, NR
!
          PI = DPI()
          NRMAX = NRAD1 + 1
          call DEFDR(OA, NSORT)
          call DEFDR(OPOT, NSORT * NRMAX)
          call DEFI(ONR, NSORT)
          call DEFDR(OR, NRMAX)
          call DEFDR(ODV, NSORT * NRMAX)
!      call hpot(W(oa),r0,nsort,W(onr),nrad1,W(opot),Z,RO,dpas,
!     .     W(or),W(odv),nradmax)
!
          allocate (POT(NRMAX, NSORT))
          allocate (DV(NRAD1), R(NRAD1), A(NSORT), NR(NSORT))
!
          call HPOT(A, R0, NSORT, NR, NRAD1, POT, Z, RO, DPAS, R, DV, &
                    NRADMAX)
!
          call RLSE(OR)
!
          VOLCEL = AV(1) * (BV(2) * CV(3) - BV(3) * CV(2)) + AV(2) &
                   * (BV(3) * CV(1) - BV(1) * CV(3)) + AV(3) * (BV(1) * CV(2) &
                                                                - BV(2) * CV(1))
          VOLCEL = abs(VOLCEL)
!    average radii:
          AVW = (VOLCEL / NATOM / (4.D0 * PI) * 3.D0)**(1.D0 / 3.D0) &
                * ALAT
          PLAT(1, 1) = AV(1)
          PLAT(1, 2) = BV(1)
          PLAT(1, 3) = CV(1)
          PLAT(2, 1) = AV(2)
          PLAT(2, 2) = BV(2)
          PLAT(2, 3) = CV(2)
          PLAT(3, 1) = AV(3)
          PLAT(3, 2) = BV(3)
          PLAT(3, 3) = CV(3)
!
          do I = 1, 3
              do J = 1, 3
                  BOUND(I, J) = PLAT(I, J)
              end do
              LTMAX(I) = 2
              NDEL(I) = 0
              ORIGIN(I) = 0
          end do
          call DEFDR(OBAS, NATOM * 3)
!
          call HRTREE(W, ALAT, AVW, BAS, BOUND, IS, LTMAX, NATOM, NSORT, &
                      NDEL, ORIGIN, PLAT, SMT, Z, WSREST, A, R0, NR, POT, NRMAX, IMT)
!
      end subroutine DEFMTR
!*==hpot.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      subroutine HPOT(A, B, NCLASS, NR, NRAD1, POT, Z, RO, DPAS, R, DV, &
                      NRADMAX)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          real*8 HUGE
          parameter(HUGE=1.D6)
!
! Dummy arguments
!
          real*8 DPAS
          integer NCLASS, NRAD1, NRADMAX
          real*8 A(*), B(*), DV(NRAD1), POT(NRAD1 + 1, &
                                                     *), R(NRAD1), RO(NRADMAX, *), Z(*)
          integer NR(*)
!
! Local variables
!
          integer I, ICL, NP, NRMAX
          real*8 ZZ
          save I, ICL, NP, NRMAX, ZZ
!
!*** End of declarations rewritten by SPAG
!
!- Calculate H. potential
! ----------------------------------------------------------------------
!i Inputs:
!i   clabl :name of the different inequivalent atom
!i   nclass:number of classes, atoms in same class are symmetry-related
!i   nrmax :maximum number of mesh points
!i   nsp   :=1 spin degenerate, =2 non-degenerate
!i   z     :nuclear charge
!o Outputs:
!o   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!o   b     :                 -//-
!o   ierr  :0 if no error occurs, ic if read error for class ic
!o   nr    :number of mesh points
!o   pot   :spherical Hartree potential
!w   v     :spherical potential (electronic contribution)
!r Remarks:
! ----------------------------------------------------------------------
!
          ! - radial mesh
          NRMAX = NRAD1 + 1
          do ICL = 1, NCLASS
              do I = 1, NRAD1
                  R(I) = B(ICL) * exp(dble(I - 1) * DPAS)
              end do
              NP = NRAD1
              call POTS(POT(2, ICL), RO(1, ICL), DV, R, DPAS, Z(ICL), &
                        NP)
              ZZ = 2 * Z(ICL)
              A(ICL) = DPAS
              NR(ICL) = NRMAX
              do I = 2, NRMAX
                  POT(I, ICL) = POT(I, ICL) - ZZ / R(I - 1)
              end do
              POT(1, ICL) = -HUGE
          end do
      end subroutine HPOT
!*==pots.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine POTS(DV, D, DP, DR, DPAS, Z, NP)
!
! iHTEgPiPOBAHiE pOTEHciAlA pO 4 TO~KAM
! DV - POTENTIAL   D CH.DENCITY  DP WORK ARRAY  DR RADIAL MESH
! DPAS EXP. PASS
! Z ATOM NUMBER     NP - NUMBER OF POINTS
! **********************************************************************
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 DPAS, Z
          integer NP
          real*8 D(*), DP(*), DR(*), DV(*)
!
! Local variables
!
          real*8 DAS, DLO, DLO2
          integer I, J, K
          save DAS, DLO, DLO2, I, J, K
!
!*** End of declarations rewritten by SPAG
!
!
          DAS = DPAS / 24.D0
          do I = 1, NP
              DV(I) = D(I) * DR(I)
          end do
          DLO = exp(DPAS)
          DLO2 = DLO * DLO
          DP(2) = DR(1) * (D(2) - D(1) * DLO2) / (12.D0 * (DLO - 1.D0))
          DP(1) = DV(1) / 3.D0 - DP(2) / DLO2
          DP(2) = DV(2) / 3.D0 - DP(2) * DLO2
          J = NP - 1
          do I = 3, J
              DP(I) = DP(I - 1) + DAS * (13.D0 * (DV(I) + DV(I - 1)) &
                                         - (DV(I - 2) + DV(I + 1)))
          end do
          DP(NP) = DP(J)
          DV(J) = DP(J)
          DV(NP) = DP(J)
          do I = 3, J
              K = NP + 1 - I
              DV(K) = DV(K + 1) / DLO + DAS * (13.D0 * (DP(K + 1) / DLO + DP &
                                                        (K)) - (DP(K + 2) / DLO2 + DP(K - 1) * DLO))
          end do
          DV(1) = DV(3) / DLO2 + DPAS * (DP(1) + 4.D0 * DP(2) / DLO + &
                                         DP(3) / DLO2) / 3.D0
          do I = 1, NP
              DV(I) = DV(I) / DR(I) * 2
          end do
      end subroutine POTS
!*==fmesh.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      real*8 function FMESH(A, B, F, NR, R)
!- Computes the value of f at r for given values on a mesh
! ----------------------------------------------------------------------
!i Inputs:
!i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!i   b     :                 -//-
!i   f     :function defined on a mesh
!i   nr    :number of mesh points
!i   r     :radial distance
!o Outputs:
!o   fmesh :interpolated value of f at r
!r Remarks:
!r   if r is insite the mesh a quadratic fit is used
!r   if r is outside the mesh an exponential fit is used
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          integer NP
          parameter(NP=4)
!
! Dummy arguments
!
          real*8 A, B, R
          integer NR
          real*8 F(*)
!
! Local variables
!
          real*8 DELSQF, DI3INT
          integer IS, NSTART
          real*8 XX
          save IS, NSTART, XX
          external DELSQF, DI3INT
!
!*** End of declarations rewritten by SPAG
!
          if (R >= B) then
              XX = log(R / B) / A + 2
              IS = IDNINT(XX) - 1
          else
              FMESH = F(2)
              return
          end if
          if (IS <= NR - 2) then
              IS = MAX0(1, IS)
              FMESH = DI3INT(IS, F(IS), XX)
          else
              NSTART = NR - NP + 1
              FMESH = DELSQF(NP, NSTART, F(NSTART), XX)
          end if
      end function FMESH
!*==iprint.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      integer function IPRINT()
!- Gets last integer of print priority stack
! ----------------------------------------------------------------------
!o Outputs:
!o   iprint:defines the verbosity level
!r Remarks:
!r   verbosity:   0  nearly nothing is printed
!r               10  very terse
!r               20  terse
!r               30  normal
!r               40  verbose
!r               50  very verbose
!r               60  highest verbosity
!r              100  low-level debugging
!r              110  intermediate-level debugging
!r              120  high-level debugging
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
          IPRINT = 30
      end function IPRINT
!*==potxn.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      subroutine POTXN(A, ALAT, B, BAS, IAX, ICLASS, NPR, NR, NRMAX, &
                       PLAT, POT, POTL, XN)
!-Calculates Hartree potential at point xn
! ----------------------------------------------------------------------
!i Inputs:
!i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!i   alat  :length scale
!i   b     :                 -//-
!i   bas   :basis vectors (scaled by alat)
!i   iax   :information about positions around a specified pair of atom
!i   iclass:the jth atom belongs to class iclass(j)
!i   npr   :number of neighbors
!i   nr    :number of mesh points
!i   nrmax :maximum number of mesh points
!i   plat  :primitive lattice vectors (scaled by alat)
!i   plat  :primitive translation vectors in real space
!i   pot   :spherical Hartree potential
!i   xn    :cartesian coordinates where potential is calculated
!o Outputs:
!i   potl  :ovelapping Hartree potential
!r Remarks:
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 ALAT, POTL
          integer NPR, NRMAX
          real*8 A(*), B(*), BAS(3, *), PLAT(3, 3), &
              POT(NRMAX, *), XN(3)
          integer IAX(5, *), ICLASS(*), NR(*)
!
! Local variables
!
          real*8 D(3), R
          real*8 DNRM23, FMESH
          integer JBAS, JC, JPR, K
          save D, JBAS, JC, JPR, K, R
          external DNRM23 !, FMESH
!
!*** End of declarations rewritten by SPAG
!
          POTL = 0.D0
          do JPR = 1, NPR
              JBAS = IAX(2, JPR)
              JC = ICLASS(JBAS)
              do K = 1, 3
                  D(K) = BAS(K, JBAS) - XN(K) + PLAT(K, 1) * IAX(3, JPR) &
                         + PLAT(K, 2) * IAX(4, JPR) + PLAT(K, 3) * IAX(5, JPR)
              end do
              R = sqrt(DNRM23(D)) * ALAT
              POTL = POTL + FMESH(A(JC), B(JC), POT(1, JC), NR(JC), &
                                  R)
          end do
!
      end subroutine POTXN
!*==di3int.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      real*8 function DI3INT(IX, YA, X)
!- Interpolates y = f(x) for given xa and ya=f(xa)
! ----------------------------------------------------------------------
!i Inputs:
!i   ix    :xa(1) is integer and equals ix
!i          and xa(1),xa(2) and xa(3) differ exactly by 1.d0
!i   ya    :value of f at xa
!i   x     :x-value at which f is interpolated
!o Outputs:
!o   di3int:interpolated value
!r Remarks:
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer IX
          real*8 X
          real*8 YA(3)
!
! Local variables
!
          real*8 XA(3)
!
!*** End of declarations rewritten by SPAG
!
          XA(1) = dble(IX)
          XA(2) = dble(IX + 1)
          XA(3) = dble(IX + 2)
!
          DI3INT = 0.5D0 * (X - XA(2)) * ((X - XA(3)) * YA(1) + &
                                          (X - XA(1)) * YA(3)) - (X - XA(1)) * (X - XA(3)) * YA(2)
      end function DI3INT
!*==i_shell.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      subroutine I_SHELL(M, N, IARRAY)
!- shell sort of a array of integer vectors
! ----------------------------------------------------------------------
!i Inputs:
!i   m     :number of components in iarray
!i   n     :number of elements in iarray
!i   iarray:array to be sorted
!o Outputs:
!o   iarray:array to be sorted
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer M, N
          integer IARRAY(M, 0:N - 1)
!
! Local variables
!
          integer I, IT, J, K, L, LOGNB2, MM, MMM, N2, NN
!
!*** End of declarations rewritten by SPAG
!
          LOGNB2 = int(log(FLOAT(N + 1)) * 1.4426950)
          N2 = N
          do NN = 1, LOGNB2
              N2 = N2 / 2
              K = N - N2
              do J = 1, K
                  I = J - 1
20                continue
                  L = I + N2
                  do MM = 1, M
                      if (IARRAY(MM, L) < IARRAY(MM, I)) then
                          do MMM = 1, M
                              IT = IARRAY(MMM, I)
                              IARRAY(MMM, I) = IARRAY(MMM, L)
                              IARRAY(MMM, L) = IT
                          end do
                          I = I - N2
                          if (I < 0) exit
                          goto 20
                      elseif (IARRAY(MM, L) /= IARRAY(MM, I)) then
                          exit
                      end if
                  end do
              end do
          end do
!
      end subroutine I_SHELL
!*==hrtree.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      subroutine HRTREE(W, ALAT, AVW, BAS, BOUND, ICLASS, LTMAX, NBAS, &
                        NCLASS, NDEL, ORIGIN, PLAT, WSR, Z, WSREST, A, B, NR, POT, NRMAX, &
                        IMT)
!- Calculate Hartree potential on a mesh and determine muffin-tin radia
! ----------------------------------------------------------------------
!i Inputs:
!i   alat  :length scale
!i   avw   :average Wigner Seitz radius
!i   bas   :basis vectors (scaled by alat)
!i   bound :two vectors spanning the plane (scaled by alat)
!i   clabl :name of the different inequivalent atom
!i   iclass:the jth atom belongs to class iclass(j)
!i   ltmax :ltmax(i)= limit in i-dirction for unit cells considered
!i   nbas  :number of atoms in the basis
!i   nclass:number of classes, atoms in same class are symmetry-related
!i   ndel  :ndel(i)=number of mesh points along the bound(i) vector
!i   nrclas:number of atoms in the i-th class
!i   origin:origin of the plane (scaled by alat)
!i   plat  :primitive lattice vectors (scaled by alat)
!i   z     :nuclear charge
!o Outputs:
!io  wsr   :Wigner-Seitz sphere radius (in atomic units)
!r Remarks:
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          real*8 FACN1
          parameter(FACN1=1.5D0)
!
! Dummy arguments
!
          real*8 ALAT, AVW
          integer IMT, NBAS, NCLASS, NRMAX
          real*8 A(*), B(*), BAS(3, *), BOUND(3, 3), &
              ORIGIN(*), PLAT(3, 3), POT(NRMAX, *), WSR(*), WSREST(* &
                                                                   ), Z(*)
          integer ICLASS(*), LTMAX(3), NDEL(*), NR(*), W(*)
!
! Local variables
!
          real*8 DSCL(:, :), WMAX
          integer I, IC, J, NEIGHM, ODSCL, OIAX1, OLOCK, ONPR1
          integer IDAMAX, IPRINT
          logical lock(:)
          save DSCL, I, IC, J, LOCK, NEIGHM, ODSCL, OIAX1, OLOCK, ONPR1, &
              WMAX
          ! external DEFDR, DEFI, DEFWSR, IPRINT, NGHBR1, POTMAX, POTSUM, &
          external DEFDR, DEFI, DEFWSR, NGHBR1, POTMAX, POTSUM, &
              RLSE
!
!*** End of declarations rewritten by SPAG
!
          ALLOCATABLEDSCL, LOCK
!
          if (IPRINT() > 90) then
              print*,'******************************************'
              print*,'alat', ALAT, '  avw', AVW
              print*,' === bas. v. ==='
              do J = 1, NBAS
                  print 99001, (BAS(I, J), I=1, 3)
              end do
              print*,' === bound ==='
              do J = 1, 3
                  print 99001, (BOUND(I, J), I=1, 3)
              end do
              print*,' === ltmax ==='
              print 99002, LTMAX
              print*,' === ndel ==='
              print 99002, (NDEL(I), I=1, 3)
              print*,' === origin ==='
              print 99001, (ORIGIN(I), I=1, 3)
              print*,' === plat ==='
              do J = 1, 3
                  print 99001, (PLAT(I, J), I=1, 3)
              end do
              print*,'******************************************'
          end if
! --- Calculate Hartree potential on a mesh
          if (NDEL(1) * NDEL(2) /= 0) call POTSUM(A, ALAT, B, BAS, &
                                                  BOUND, ICLASS, LTMAX, NBAS, NDEL, NR, NRMAX, ORIGIN, PLAT, POT)
! --- Determine maximum of Hartree potential
          do IC = 1, NCLASS
              call DEFWSR(WSREST(IC), Z(IC))
              WSR(IC) = WSREST(IC)
          end do
          WMAX = WSREST(IDAMAX(NCLASS, WSREST, 1)) / AVW
          I = 2 * int((2.D0 * FACN1 * WMAX + 1.D0)**3)
          NEIGHM = max(2 * int((2.D0 * FACN1 * WMAX + 1.D0)**3), &
                       50)
          call DEFI(OIAX1, 5 * NEIGHM * NCLASS)
          call DEFI(ONPR1, NCLASS)
          call NGHBR1(W, ALAT, BAS, FACN1, W(OIAX1), ICLASS, NBAS, NCLASS, &
                      NEIGHM, W(ONPR1), PLAT, WSREST)
!
          call POTMAX(W, A, ALAT, B, BAS, W(OIAX1), ICLASS, NBAS, NCLASS, &
                      W(ONPR1), NR, NRMAX, PLAT, POT, WMAX, WSR, WSREST, Z)
          call DEFDR(ODSCL, NCLASS * (NCLASS + 2))
          call DEFDR(OLOCK, NCLASS)
!     call blowup(alat,bas,1.d0,0.d0,w(oiax1),iclass,nbas,
!    $           nclass,w(onpr1),0.d0,0.d0,plat,wsr,W(odscl),W(olock))
!     if(imt.ne.0)call blowup(alat,bas,1.d0,0.d0,w(oiax1),iclass,nbas,
!    $         nclass,w(onpr1),0.4d0,0.8d0,plat,wsr,W(odscl),W(olock))
!
          allocate (DSCL(NCLASS, 0:NCLASS + 1))
          allocate (lock(NCLASS))
!
          call BLOWUP(ALAT, BAS, 1.D0, 0.D0, W(OIAX1), ICLASS, NBAS, &
                      NCLASS, W(ONPR1), 0.D0, 0.D0, PLAT, WSR, DSCL, LOCK)
          if (IMT /= 0) call BLOWUP(ALAT, BAS, 1.D0, 0.D0, W(OIAX1), &
                                    !ICLASS, NBAS, NCLASS, W(ONPR1), 0.1D0, 0.3D0, PLAT, WSR, DSCL, &
                                    ICLASS, NBAS, NCLASS, W(ONPR1), 0.4D0, 0.8D0, PLAT, WSR, DSCL, &
                                    LOCK)
!
!
          call RLSE(OIAX1)
99001     format(20F21.15)
99002     format(20I10)
!
      end subroutine HRTREE
!*==potsum.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      subroutine POTSUM(A, ALAT, B, BAS, BOUND, ICLASS, LTMAX, NBAS, &
                        NDEL, NR, NRMAX, ORIGIN, PLAT, POT)
!- Calculates the overlapping Hartree potential
! ----------------------------------------------------------------------
!i Inputs:
!i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!i   b     :                 -//-
!i   alat  :length scale
!i   bas   :basis vectors (scaled by alat)
!i   bound :two vectors spanning the plane (scaled by alat)
!i   iclass:the jth atom belongs to class iclass(j)
!i   ltmax :ltmax(i)= limit in i-dirction for unit cells considered
!i   nbas  :number of atoms in the basis
!i   ndel  :ndel(i)=number of mesh points along the bound(i) vector
!i   nr    :number of mesh points
!i   nrmax :maximum number of mesh points
!i   origin:origin of the plane (scaled by alat)
!i   plat  :primitive lattice vectors (scaled by alat)
!i   pot   :spherical Hartree potential
!o Output written to file POT
!r Remarks:
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 ALAT
          integer NBAS, NRMAX
          real*8 A(*), B(*), BAS(3, *), BOUND(3, 3), &
              ORIGIN(3), PLAT(3, 3), POT(NRMAX, *)
          integer ICLASS(*), LTMAX(3), NDEL(*), NR(*)
!
! Local variables
!
          real*8 D01, D02, D03, D11, D12, D13, D21, D22, D23, D31, &
              D32, D33, DN1, DN2, FAC1, FAC2, POTL, RAD, XN1, XN2, XN3
          real*8 FMESH
          integer I, I1, I2, IBAS, IC, J, J1, J2, J3
          save D01, D02, D03, D11, D12, D13, D21, D22, D23, D31, D32, D33, &
              DN1, DN2, FAC1, FAC2, I, I1, I2, IBAS, IC, J, J1, J2, J3, POTL, &
              RAD, XN1, XN2, XN3
!
!*** End of declarations rewritten by SPAG
!
          if (NDEL(1) /= 1) FAC1 = 1.D0 / dble(NDEL(1) - 1)
          if (NDEL(2) /= 1) FAC2 = 1.D0 / dble(NDEL(2) - 1)
          write (16, 99001) ORIGIN, ((BOUND(I, J), I=1, 3), J=1, 2)
          write (16, 99002)
!
          do I1 = 0, NDEL(1) - 1
              DN1 = I1 * FAC1
              do I2 = 0, NDEL(2) - 1
                  DN2 = I2 * FAC2
                  XN1 = ORIGIN(1) + BOUND(1, 1) * DN1 + BOUND(1, 2) * DN2
                  XN2 = ORIGIN(2) + BOUND(2, 1) * DN1 + BOUND(2, 2) * DN2
                  XN3 = ORIGIN(3) + BOUND(3, 1) * DN1 + BOUND(3, 2) * DN2
                  POTL = 0.D0
                  do IBAS = 1, NBAS
                      IC = ICLASS(IBAS)
! --------  ltmax(i) lattice translations in i-direction considered
                      D01 = BAS(1, IBAS) - XN1
                      D02 = BAS(2, IBAS) - XN2
                      D03 = BAS(3, IBAS) - XN3
                      do J1 = -LTMAX(1), LTMAX(1)
                          D11 = D01 + PLAT(1, 1) * J1
                          D12 = D02 + PLAT(2, 1) * J1
                          D13 = D03 + PLAT(3, 1) * J1
                          do J2 = -LTMAX(2), LTMAX(2)
                              D21 = D11 + PLAT(1, 2) * J2
                              D22 = D12 + PLAT(2, 2) * J2
                              D23 = D13 + PLAT(3, 2) * J2
                              do J3 = -LTMAX(3), LTMAX(3)
                                  D31 = D21 + PLAT(1, 3) * J3
                                  D32 = D22 + PLAT(2, 3) * J3
                                  D33 = D23 + PLAT(3, 3) * J3
                                  RAD = sqrt(D31 * D31 + D32 * D32 + D33 * D33) * ALAT
                                  POTL = POTL + FMESH(A(IC), B(IC), POT(1, IC), NR(IC), &
                                                      RAD)
                              end do
                          end do
                      end do
                  end do
              end do
          end do
!
99001     format(//25('-')//, 'make plot for plane', //, 'ORIGIN:',           &
                         &        3F21.15/'R1    :', 3F21.15, /'R2    :', 3F21.15//, 25('-'))
99002     format('Begin to make POT ...')
      end subroutine POTSUM
!*==nghbr1.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine NGHBR1(W, ALAT, BAS, FACN1, IAX, ICLASS, NBAS, NCLASS, &
                        NEIGHM, NPR, PLAT, WSR)
!- Create a table of all neighbors within facn1*(wsr(1)+wsr(2))
! ----------------------------------------------------------------------
!i Inputs:
!i   alat  :length scale
!i   bas   :basis vectors (scaled by alat)
!i   clabl :name of the different inequivalent atom
!i   facn1 :see remarks
!i   iclass:the jth atom belongs to class iclass(j)
!i   nbas  :number of atoms in the basis
!i   nclass:number of classes, atoms in same class are symmetry-related
!i   neighm:maximum number of neighbors
!i   plat  :primitive lattice vectors (scaled by alat)
!i   wsr   :Wigner-Seitz sphere radius (in atomic units)
!o Outputs:
!i   iax   :see remarks
!i   npr   :number of neighbors around each atom
!r Remarks:
!r   Creates a neighour list for a specified atom, generating iax
!r   which contains all neigbors which fulfill:
!r
!r        distance(i,j) <= (wsr(i)+wsr(j))*facn1
!r
!r    iax(1): not used
!r    iax(2): ibas = atom in cluster
!r    iax(3): i
!r    iax(4): j
!r    iax(5): k
!r   To be sure that at least one pair is found, fac is increased
!r   by 1.2 until a neighbor is found.
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 ALAT, FACN1
          integer NBAS, NCLASS, NEIGHM
          real*8 BAS(3, *), PLAT(3, 3), WSR(*)
          integer IAX(5, NCLASS, *), ICLASS(*), NPR(*), W(*)
!
! Local variables
!
          real*8 D2, DR(3), FAC, WI, WJ, WJPWI, WJPWI2
          real*8 DRR2
          integer I, I1, I2, I3, IBAS, IPR, IPR5, J, JBAS, JC, K, OIWK
          integer ICLBAS, IPRINT
          character(72) MESSG
          save D2, DR, FAC, I, I1, I2, I3, IBAS, IPR, IPR5, J, JBAS, JC, K, &
              MESSG, OIWK, WI, WJ, WJPWI, WJPWI2
!
!*** End of declarations rewritten by SPAG
!
          call DEFI(OIWK, NEIGHM)
          do JC = 1, NCLASS
              JBAS = ICLBAS(JC, ICLASS, NBAS, 1)
              WJ = WSR(JC)
              IPR = 0
              FAC = FACN1
              do while (IPR < 2)
                  IPR = 0
                  do IBAS = 1, NBAS
                      WI = WSR(ICLASS(IBAS))
                      WJPWI = (WJ + WI) / ALAT * FAC
                      call LATLIM(PLAT, WJPWI, I1, I2, I3)
                      WJPWI2 = WJPWI * WJPWI
! --------- Sweep lattice translations to find all neighbors
! --------- within wjpwi
                      do I = -I1, I1
                          do J = -I2, I2
                              do K = -I3, I3
                                  D2 = DRR2(PLAT, BAS(1, JBAS), BAS(1, IBAS), I, J, K, DR)
                                  if (D2 <= WJPWI2) then
                                      if (IPR >= NEIGHM) then
                                          write (MESSG, 99002) NEIGHM, I1, I2, I3
                                          call ERRMSG(MESSG, 4)
                                      end if
                                      IPR5 = IPR * 5
                                      W(OIWK + IPR5 + 0) = 10000 * D2
                                      W(OIWK + IPR5 + 1) = IBAS
                                      W(OIWK + IPR5 + 2) = I
                                      W(OIWK + IPR5 + 3) = J
                                      W(OIWK + IPR5 + 4) = K
                                      IPR = IPR + 1
                                  end if
                              end do
                          end do
                      end do
                  end do
                  FAC = FAC * 1.2D0
              end do
!
              call I_SHELL(5, IPR, W(OIWK))
              do I = 1, IPR
                  call INCOPY(5, W(OIWK + 5 * (I - 1)), 1, IAX(1, JC, I), &
                              1)
                  IAX(1, JC, I) = JBAS
              end do
              NPR(JC) = IPR
!
!
! ----- Printout
              if (IPRINT() >= 80) then
!
                  do IPR = 1, NPR(JC)
                      IBAS = IAX(2, JC, IPR)
                      I = IAX(3, JC, IPR)
                      J = IAX(4, JC, IPR)
                      K = IAX(5, JC, IPR)
                      D2 = DRR2(PLAT, BAS(1, JBAS), BAS(1, IBAS), I, J, K, DR)
                  end do
                  write (16, *)
              end if
!
          end do
          if (IPRINT() >= 70) write (16, 99001) NEIGHM, (NPR(JC), JC=1, &
                                                        NCLASS)
!
99001     format(/' NGHBR1: neighm=', i4, '> npr=', 50(11I4, /26X))
99002     format(' NGHBR1: too many pairs, neighm=', i3, 6X, 'i1,i2,i3 =', 3I3,&
                         &        '$')
      end subroutine NGHBR1
!*==potmax.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      subroutine POTMAX(W, A, ALAT, B, BAS, IAX1, ICLASS, NBAS, NCLASS, &
                        NPR1, NR, NRMAX, PLAT, POT, WMAX, WSR, WSREST, Z)
!- Calculates the overlapping Hartree potential and finds it's maximum
! ----------------------------------------------------------------------
!i Inputs:
!i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!i   alat  :length scale
!i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!i   bas   :basis vectors (scaled by alat)
!i   clabl :name of the different inequivalent atom
!i   iax1  :information about positions around a specified atom
!i   iclass:the jth atom belongs to class iclass(j)
!i   nbas  :number of atoms in the basis
!i   nclass:number of classes, atoms in same class are symmetry-related
!i   npr1  :number of neighbors
!i   nr    :number of mesh points
!i   nrmax :maximum number of mesh points
!i   plat  :primitive lattice vectors (scaled by alat)
!i   pot   :spherical Hartree potential
!i   wmax  :largest wsr/avw
!i   wsrest:estimation of Wigner-Seitz sphere radius (in atomic units)
!i   z     :nuclear charge
!o Outputs:
!o   wsr   :Wigner-Seitz sphere radius (in atomic units)
!r Remarks:
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          integer NSTEP
          parameter(NSTEP=20)
          real*8 FACN2, TINY, TOL
          parameter(FACN2=3.D0, TINY=1.D-5, TOL=1.D-2)
!
! Dummy arguments
!
          real*8 ALAT, WMAX
          integer NBAS, NCLASS, NRMAX
          real*8 A(*), B(*), BAS(3, *), PLAT(3, 3), &
              POT(NRMAX, *), WSR(*), WSREST(*), Z(*)
          integer IAX1(5, NCLASS, *), ICLASS(*), NPR1(*), NR(*), &
              W(*)
!
! Local variables
!
          real*8 DLAMCH, DRR2
          real*8 DR(3), POTL, RIK, RIKO, RPMX, STEP, X(:), &
              XN(3), Y(:)
          integer I, IBAS, IC, IMIN, ISTEP, K, K1, K2, K3, KBAS, KC, KCO, &
              KPR, NEIGHM, NPR2, OIAX2
          integer ICLBAS, IPRINT
          save DR, I, IBAS, IC, IMIN, ISTEP, K, K1, K2, K3, KBAS, KC, KCO, &
              KPR, NEIGHM, NPR2, OIAX2, POTL, RIK, RIKO, RPMX, STEP, X, XN, Y
          external DEFI, DLAMCH, DSCAL, NGHBR2, POTXN
          !external DEFI, DLAMCH, DRR2, DSCAL, ICLBAS, IPRINT, NGHBR2, POTXN
!
!*** End of declarations rewritten by SPAG
!
          ALLOCATABLEx, y
          allocate (X(0:NSTEP), Y(0:NSTEP))
          NEIGHM = MAX0(int(4.D0 * (FACN2 * WMAX + 1.D0)**3), 50)
          call DEFI(OIAX2, 5 * NEIGHM)
!
          do IC = 1, NCLASS
              if (IDNINT(Z(IC)) /= 0) then
                  if (IPRINT() >= 40) write (16, 99001)
                  IBAS = ICLBAS(IC, ICLASS, NBAS, 1)
!
                  WSR(IC) = DLAMCH('o')
                  RIKO = -1.D0
                  KCO = -1
                  do KPR = 2, NPR1(IC)
                      KBAS = IAX1(2, IC, KPR)
                      K1 = IAX1(3, IC, KPR)
                      K2 = IAX1(4, IC, KPR)
                      K3 = IAX1(5, IC, KPR)
                      KC = ICLASS(KBAS)
                      RIK = sqrt(DRR2(PLAT, BAS(1, IBAS), BAS(1, KBAS), K1, K2, &
                                      K3, DR)) * ALAT
                      if (IDNINT(Z(KC)) /= 0 .and. (abs(RIK - RIKO) &
                                                    > TINY .or. KCO /= KC)) then
!
                          RIKO = RIK
                          KCO = KC
                          STEP = RIK / NSTEP
                          call NGHBR2(ALAT, BAS, FACN2, W(OIAX2), IBAS, ICLASS, &
                                      K1, K2, K3, KBAS, NBAS, NEIGHM, NPR2, PLAT, WSREST)
                          call DSCAL(3, 1.D0 / DFLOAT(NSTEP), DR, 1)
! --------  Find aproximate position of first maximun
                          do ISTEP = 0, NSTEP
                              do K = 1, 3
                                  XN(K) = BAS(K, IBAS) + ISTEP * DR(K)
                              end do
                              call POTXN(A, ALAT, B, BAS, W(OIAX2), ICLASS, NPR2, NR, &
                                         NRMAX, PLAT, POT, POTL, XN)
                              X(ISTEP) = dble(ISTEP) * STEP
                              Y(ISTEP) = max(POTL, -1.D5)
                              if (ISTEP /= 0) then
                                  if (ISTEP >= 2 .and. Y(ISTEP) < Y(ISTEP - 1)) goto &
                                      5
                              end if
                          end do
                          call ERRMSG('bug in potmax.$', 4)
!
! --------- Now a single maximum is supposed between istep-2 and istep
5                         continue
                          do K = 1, 3
                              X(K) = X(ISTEP - 3 + K)
                              Y(K) = Y(ISTEP - 3 + K)
                          end do
                          RPMX = 0.D0
                          IMIN = (ISTEP - 2) * 2
                          do while (X(3) - X(1) > TOL)
                              X(5) = X(3)
                              X(3) = X(2)
                              Y(5) = Y(3)
                              Y(3) = Y(2)
                              STEP = STEP * 0.5D0
                              call DSCAL(3, 0.5D0, DR, 1)
                              do I = 1, 3, 2
                                  do K = 1, 3
                                      XN(K) = BAS(K, IBAS) + (IMIN + I) * DR(K)
                                  end do
                                  call POTXN(A, ALAT, B, BAS, W(OIAX2), ICLASS, NPR2, NR, &
                                             NRMAX, PLAT, POT, Y(I + 1), XN)
                                  X(I + 1) = (IMIN + I) * STEP
                              end do
                              do I = 2, 4
                                  if (Y(I + 1) < Y(I)) then
                                      do K = 1, 3
                                          Y(K) = Y(I - 2 + K)
                                          X(K) = X(I - 2 + K)
                                      end do
                                      exit
                                  end if
                              end do
                              IMIN = 2 * (IMIN + I - 2)
                          end do
                          RPMX = 0.5D0 * (Y(1) * (X(2) + X(3)) - 2.D0 * Y(2) &
                                          * (X(1) + X(3)) + Y(3) * (X(1) + X(2))) / (Y(1) &
                                                                                     - Y(2) - Y(2) + Y(3))
                          WSR(IC) = min(WSR(IC), RPMX)
                      end if
                  end do
                  if (IPRINT() >= 40) write (16, 99002)
!        if (iprint().ge.20) write(6,303)clabl(ic),wsr(ic)
              end if
          end do
!
99001     format(/' POTMAX: ATOM1   ATOM2   DIST      VMAX   R(VMAX)', /, 8X,&
                         &        42('-'))
99002     format(8X, 42('-'))
!
      end subroutine POTMAX
!*==nghbr2.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine NGHBR2(ALAT, BAS, FACN2, IAX, IBAS, ICLASS, K1, K2, &
                        K3, KBAS, NBAS, NEIGHM, NPR, PLAT, WSR)
!- Create a table of all neighbors around two atoms
! ----------------------------------------------------------------------
!i Inputs:
!i   alat  :length scale
!i   bas   :basis vectors (scaled by alat)
!i   clabl :name of the different inequivalent atom
!i   facn2 :see remarks
!i   ibas  :index to first atom
!i   iclass:the jth atom belongs to class iclass(j)
!i   k1-3  :second atom is translated by k1*plat1+k2*plat2+p3*plat3
!i   kbas  :index to second atom
!i   nbas  :number of atoms in the basis
!i   neighm:maximum number of neighbors
!i   plat  :primitive lattice vectors (scaled by alat)
!i   wsr   :Wigner-Seitz sphere radius (in atomic units)
!o Outputs:
!i   iax   :see remarks
!i   npr   :number of neighbors
!r Remarks:
!r   Creates a neighour list for a pair of atoms ibas and kbas,
!r   generating iax which contains all neigbors jbas within
!r   facn2*wsr(jbas) around ibas or kbas , i.e. those which fulfill:
!r
!r        distance(i,j) <= wsr(j)*facn2
!r    or  distance(k,j) <= wsr(j)*facn2
!r
!r    iax(1): not used
!r    iax(2): ibas = atom in cluster
!r    iax(3): i
!r    iax(4): j
!r    iax(5): k
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 ALAT, FACN2
          integer IBAS, K1, K2, K3, KBAS, NBAS, NEIGHM, NPR
          real*8 BAS(3, *), PLAT(3, 3), WSR(*)
          integer IAX(5, *), ICLASS(*)
!
! Local variables
!
          real*8 B(3, 2), D1, D2, DB11, DB12, DB21, DB22, DB31, &
              DB32, DR11, DR12, DR13, DR21, DR22, DR23, DR31, DR32, DR33, DRB11, &
              DRB12, DRB21, DRB22, DRB31, DRB32, WJF, WJF2
          integer I, I1, I2, I3, II, IPR, J, JBAS, JC, K
          integer IPRINT
          character(72) MESSG
          save B, D1, D2, DB11, DB12, DB21, DB22, DB31, DB32, DR11, DR12, &
              DR13, DR21, DR22, DR23, DR31, DR32, DR33, DRB11, DRB12, DRB21, &
              DRB22, DRB31, DRB32, I, I1, I2, I3, II, IPR, J, JBAS, JC, K, &
              MESSG, WJF, WJF2
          external DCOPY, LATLIM
          ! external DCOPY, IPRINT, LATLIM
!
!*** End of declarations rewritten by SPAG
!
          call DCOPY(3, BAS(1, IBAS), 1, B(1, 1), 1)
          B(1, 2) = BAS(1, KBAS) + K1 * PLAT(1, 1) + K2 * PLAT(1, 2) &
                    + K3 * PLAT(1, 3)
          B(2, 2) = BAS(2, KBAS) + K1 * PLAT(2, 1) + K2 * PLAT(2, 2) &
                    + K3 * PLAT(2, 3)
          B(3, 2) = BAS(3, KBAS) + K1 * PLAT(3, 1) + K2 * PLAT(3, 2) &
                    + K3 * PLAT(3, 3)
          IPR = 0
          do JBAS = 1, NBAS
              JC = ICLASS(JBAS)
              WJF = WSR(JC) / ALAT * FACN2
              call LATLIM(PLAT, WJF, I1, I2, I3)
              WJF2 = WJF * WJF
              DB11 = BAS(1, JBAS) - B(1, 1)
              DB21 = BAS(2, JBAS) - B(2, 1)
              DB31 = BAS(3, JBAS) - B(3, 1)
              DB12 = BAS(1, JBAS) - B(1, 2)
              DB22 = BAS(2, JBAS) - B(2, 2)
              DB32 = BAS(3, JBAS) - B(3, 2)
              do I = -I1, I1
                  DR11 = PLAT(1, 1) * I
                  DR21 = PLAT(2, 1) * I
                  DR31 = PLAT(3, 1) * I
                  do J = -I2, I2
                      DR12 = DR11 + PLAT(1, 2) * J
                      DR22 = DR21 + PLAT(2, 2) * J
                      DR32 = DR31 + PLAT(3, 2) * J
                      do K = -I3, I3
                          DR13 = DR12 + PLAT(1, 3) * K
                          DR23 = DR22 + PLAT(2, 3) * K
                          DR33 = DR32 + PLAT(3, 3) * K
                          DRB11 = DR13 + DB11
                          DRB21 = DR23 + DB21
                          DRB31 = DR33 + DB31
                          DRB12 = DR13 + DB12
                          DRB22 = DR23 + DB22
                          DRB32 = DR33 + DB32
                          D1 = DRB11 * DRB11 + DRB21 * DRB21 + DRB31 * DRB31
                          D2 = DRB12 * DRB12 + DRB22 * DRB22 + DRB32 * DRB32
                          if (D1 <= WJF2 .or. D2 <= WJF2) then
                              IPR = IPR + 1
                              if (IPR >= 2 * NEIGHM) then
                                  write (MESSG, 99002) 2 * NEIGHM, I1, I2, I3
                                  call ERRMSG(MESSG, 4)
                              end if
                              IAX(2, IPR) = JBAS
                              IAX(3, IPR) = I
                              IAX(4, IPR) = J
                              IAX(5, IPR) = K
                              II = 0
                              if (D1 <= WJF2) II = II + 1
                              if (D2 <= WJF2) II = II + 2
! ------------- Printout
                              if (IPRINT() >= 80) write (16, 99001) IPR, D1, D2, JBAS, I, &
                                  J, K, II
                          end if
                      end do
                  end do
              end do
          end do
!
          NPR = IPR
!
99001     format(i3, 2F21.15, 4I5, 2X, i3)
99002     format(' NGHBR2: too many pairs, neighm=', i3, 6X, 'i1,i2,i3 =', 3I3,&
                         &        '$')
      end subroutine NGHBR2
!*==blowup.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      subroutine BLOWUP(ALAT, BAS, FACVOL, GAMMA, IAX1, ICLASS, NBAS, &
                        NCLASS, NPR1, OMMAX1, OMMAX2, PLAT, WSR, DSCL, LOCK)
!- Blows up the spheres
! ----------------------------------------------------------------------
!i Inputs:
!i   alat  :length scale
!i   bas   :basis vectors (scaled by alat)
!i   clabl :name of the different inequivalent atom
!i   iax1  :information about positions around a specified atom
!i   iclass:the jth atom belongs to class iclass(j)
!i   nbas  :number of atoms in the basis
!i   nclass:number of classes, atoms in same class are symmetry-related
!i   npr1  :number of neighbors
!i   facvol:sum of sphere volumes = facvol * cell volume
!i   gamma :scaling is r -> a(r+b) with a*b=gamma*(a-1)*avw
!i   ommax1:maximum overlap divided by distance (s1+s2-d)/d  < ommax1
!i   ommax2:maximum overlap divided by  radius  (s1+s2-d)/s1 < ommax2
!i   plat  :primitive lattice vectors (scaled by alat)
!io  wsr   :Wigner-Seitz sphere radius (in atomic units)
!o Outputs:
!io  wsr   :Wigner-Seitz sphere radius (in atomic units)
!r Remarks:
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          real*8 FPI3, TINY
          parameter(FPI3=4.18879020478639053D0, TINY=1.D-4)
!
! Dummy arguments
!
          real*8 ALAT, FACVOL, GAMMA, OMMAX1, OMMAX2
          integer NBAS, NCLASS
          real*8 BAS(3, *), DSCL(NCLASS, 0:NCLASS + 1), PLAT(3, &
                                                                      3), WSR(*)
          integer IAX1(5, NCLASS, *), ICLASS(*), NPR1(*)
          logical lock(*)
!
! Local variables
!
          real*8 A, AMAX, AMAX1, AMAX2, AMAX3, AVW, B, BMAX, D, DM &
              (0:3), DR(3), DSCLMX, GW, OMO2, OPO1, P, Q, R, RIK, RIKO, S, T, &
              TMO2, U, V, VOL, VOLA, VOLB, VOLSPH, WSRI, WSRK, X
          real*8 DDET33, DLAMCH, DRR2
          logical FIN
          integer IBAS, IC, ICMAX, ICMIN, ILOOP, K1, K2, K3, KBAS, KC, KCO, &
              KPR, LINE, NLOOP
          integer ICLBAS, IPRINT
          character(72) MESSG
          save A, AMAX, AMAX1, AMAX2, AMAX3, AVW, B, BMAX, D, DM, DR, &
              DSCLMX, FIN, GW, IBAS, IC, ICMAX, ICMIN, ILOOP, K1, K2, K3, KBAS, &
              KC, KCO, KPR, LINE, MESSG, NLOOP, OMO2, OPO1, P, Q, R, RIK, RIKO, &
              S, T, TMO2, U, V, VOL, VOLA, VOLB, VOLSPH, WSRI, WSRK, X
          external ERRMSG
          ! external DLAMCH, DRR2, ERRMSG, ICLBAS, IPRINT
!
!*** End of declarations rewritten by SPAG
!
          OPO1 = 1.D0 + OMMAX1
          OMO2 = 1.D0 - OMMAX2
          TMO2 = 2.D0 - OMMAX2
          VOL = abs(DDET33(PLAT)) * ALAT**3
          AVW = (VOL / FPI3 / NBAS)**(1.D0 / 3.D0)
          GW = AVW * GAMMA
!
          VOLSPH = 0.D0
          do IBAS = 1, NBAS
              VOLSPH = VOLSPH + FPI3 * WSR(ICLASS(IBAS))**3
          end do
          if (IPRINT() >= 30) write (16, 99009) VOL, VOLSPH
!
!     call dcopy(nclass,1.d0,0,dscl(1,0),1)
          DSCL(1:NCLASS, 0) = 1.D0
!
          do NLOOP = 1, NCLASS + 1
              AMAX = DLAMCH('o')
! ----- lock radii of those spheres whith maximum allowed overlap
              lock(1:NCLASS) = .false.
              do IC = 1, NCLASS
                  IBAS = ICLBAS(IC, ICLASS, NBAS, 1)
                  WSRI = WSR(IC)
                  do KPR = 2, NPR1(IC)
                      KBAS = IAX1(2, IC, KPR)
                      K1 = IAX1(3, IC, KPR)
                      K2 = IAX1(4, IC, KPR)
                      K3 = IAX1(5, IC, KPR)
                      KC = ICLASS(KBAS)
                      WSRK = WSR(KC)
                      RIK = sqrt(DRR2(PLAT, BAS(1, IBAS), BAS(1, KBAS), K1, K2, K3, &
                                      DR))
                      RIK = RIK * ALAT
                      if (abs(OPO1 * RIK - WSRI - WSRK) < TINY) lock(IC) = .true.
                      if (abs(RIK - OMO2 * WSRI - WSRK) < TINY) lock(IC) = .true.
                      if (abs(RIK - WSRI - OMO2 * WSRK) < TINY) lock(IC) = .true.
                  end do
              end do
!
              do IC = 1, NCLASS
                  if (.not. lock(IC)) then
                      IBAS = ICLBAS(IC, ICLASS, NBAS, 1)
                      RIKO = -1.D0
                      KCO = -1
                      WSRI = WSR(IC)
                      do KPR = 2, NPR1(IC)
                          KBAS = IAX1(2, IC, KPR)
                          K1 = IAX1(3, IC, KPR)
                          K2 = IAX1(4, IC, KPR)
                          K3 = IAX1(5, IC, KPR)
                          KC = ICLASS(KBAS)
                          RIK = sqrt(DRR2(PLAT, BAS(1, IBAS), BAS(1, KBAS), K1, K2, &
                                          K3, DR))
                          RIK = RIK * ALAT
                          if (abs(RIK - RIKO) > TINY .or. KCO /= KC) then
                              WSRK = WSR(KC)
                              RIKO = RIK
                              KCO = KC
                              if (lock(KC)) then
                                  AMAX1 = (OPO1 * RIK - WSRK + GW) / (WSRI + GW)
                                  if (OMO2 > 0.D0) AMAX2 = (RIK - WSRK + OMO2 * GW) &
                                                           / (WSRI + GW) / OMO2
                                  AMAX3 = (RIK - OMO2 * WSRK + GW) / (WSRI + GW)
                              else
                                  AMAX1 = (OPO1 * RIK + GW + GW) / (WSRI + WSRK + GW + &
                                                                    GW)
                                  if (WSRK + OMO2 * WSRI + TMO2 * GW > 0.D0) AMAX2 = &
                                      (RIK + TMO2 * GW) / (WSRK + OMO2 * WSRI + TMO2 * GW)
                                  if (WSRI + OMO2 * WSRK + TMO2 * GW > 0.D0) AMAX3 = &
                                      (RIK + TMO2 * GW) / (WSRI + OMO2 * WSRK + TMO2 * GW)
                              end if
                              AMAX = min(AMAX, AMAX1, AMAX2, AMAX3)
                          end if
                      end do
                  end if
              end do
              BMAX = (1.D0 - 1.D0 / AMAX) * GW
!
              VOLA = 0.D0
              VOLB = 0.D0
              do IBAS = 1, NBAS
                  IC = ICLASS(IBAS)
                  if (.not. lock(IC)) then
                      VOLB = VOLB + (AMAX * WSR(IC) + BMAX)**3
                  else
                      VOLA = VOLA + WSR(IC)**3
                  end if
              end do
              VOLA = VOLA * FPI3
              VOLB = VOLB * FPI3
!
              if (VOL * FACVOL < VOLA + VOLB) then
!
                  DM(:) = 0.0D0
                  do IBAS = 1, NBAS
                      IC = ICLASS(IBAS)
                      if (.not. lock(IC)) then
                          A = WSR(IC) + GW
! ----------- for numerical reasons distinguish cases
                          if (abs(GAMMA) > 1.D0) then
                              B = WSR(IC)
                          else
                              B = -GW
                          end if
                          DM(0) = DM(0) + B * B * B
                          DM(1) = DM(1) + 3.D0 * A * B * B
                          DM(2) = DM(2) + 3.D0 * A * A * B
                          DM(3) = DM(3) + A * A * A
                      end if
                  end do
                  if (abs(DM(3)) > TINY) then
                      R = DM(2) / DM(3)
                      S = DM(1) / DM(3)
                      T = (DM(0) - (VOL * FACVOL - VOLA) / FPI3) / DM(3)
                      P = S - R * R / 3.D0
                      Q = 2.D0 * R * R * R / 27.D0 - R * S / 3.D0 + T
                      D = P * P * P / 27.D0 + Q * Q / 4.D0
                      U = (sqrt(D) - Q / 2.D0)**(1.D0 / 3.D0)
                      V = -P / U / 3.D0
                      X = U + V - R / 3.D0
                      if (abs(GAMMA) > 1.D0) then
                          AMAX = X + 1.D0
                          BMAX = X * GW / AMAX
                      else
                          AMAX = X
                          BMAX = (1.D0 - 1.D0 / AMAX) * GW
                      end if
                      if (IPRINT() >= 100) then
                          write (16, 99008) 'R S T', R, S, T
                          write (16, 99008) 'P Q  ', P, Q
                          write (16, 99008) '  D  ', D
                          write (16, 99008) ' U V ', U, V
                          write (16, 99008) ' AMAX', AMAX
                          write (16, 99008) ' BMAX', BMAX
                          write (16, 99008) ' -------------------------'
                      end if
                  end if
              end if
!
              FIN = .true.
              DSCLMX = 1.D0
              do IC = 1, NCLASS
                  DSCL(IC, NLOOP) = 1.D0
                  if (.not. lock(IC)) then
                      DSCLMX = AMAX + BMAX / WSR(IC)
                      WSR(IC) = DSCLMX * WSR(IC)
                      DSCL(IC, 0) = DSCLMX * DSCL(IC, 0)
                      DSCL(IC, NLOOP) = DSCLMX
                      FIN = .false.
                  end if
              end do
              FIN = FIN .or. abs(DSCLMX - 1.D0) < TINY / 256
!
              VOLSPH = 0.D0
              do IBAS = 1, NBAS
                  VOLSPH = VOLSPH + FPI3 * WSR(ICLASS(IBAS))**3
              end do
!
!
              if (FIN .or. NLOOP == NCLASS + 1) then
                  if (IPRINT() >= 30) then
                      write (16, 99001) OMMAX1, OMMAX2
                      do LINE = 0, (NCLASS - 1) / 6
                          ICMIN = 1 + 6 * LINE
                          ICMAX = MIN0(NCLASS, 6 + 6 * LINE)
!              write(6,301)(clabl(ic),ic=icmin,icmax)
                          write (16, 99007) ('============', IC=ICMIN, ICMAX + 1)
                          write (16, 99002) (WSR(IC) / DSCL(IC, 0), IC=ICMIN, &
                                            ICMAX)
                          write (16, 99007) ('------------', IC=ICMIN, ICMAX + 1)
                          do ILOOP = 1, NLOOP
                              write (16, 99003) ILOOP, (DSCL(IC, ILOOP), IC=ICMIN, &
                                                       ICMAX)
                          end do
                          write (16, 99007) ('------------', IC=ICMIN, ICMAX + 1)
                          write (16, 99004) (WSR(IC), IC=ICMIN, ICMAX)
                          write (16, 99005) (DSCL(IC, 0), IC=ICMIN, ICMAX)
                          write (16, 99006) (WSR(IC) - WSR(IC) / DSCL(IC, 0), &
                                            IC=ICMIN, ICMAX)
                          write (16, 99007) ('------------', IC=ICMIN, ICMAX + 1)
                      end do
                      write (16, 99009) VOL, VOLSPH
                  end if
                  if (abs(OMMAX1 * OMMAX2) > TINY .and. abs(VOLSPH - VOL * &
                                                            FACVOL) > TINY) then
                      write (MESSG, 99010)
                      call ERRMSG(MESSG, 0)
                  end if
                  return
              end if
          end do
!
99001     format(/' BLOWUP: scale radii;   OMMAX1=', f21.15, ', OMMAX2=',    &
                         &        f21.15)
99002     format(6X, 'OLD WSR:    ', 6F21.15)
99003     format(6X, 'LOOP:', i3, '  * ', 6F21.15)
99004     format(6X, 'NEW WSR:    ', 6F21.15)
99005     format(6X, 'W_NEW/W_OLD:', 6F21.15)
99006     format(6X, 'W_NEW-W_OLD:', 6F21.15)
99007     format(6X, a12, 6A9)
99008     format(6X, a, 3F21.15)
99009     format(/' CELL VOLUME=', f21.15, '   SUM OF SPHERE VOLUMES=',      &
                         &        f21.15)
99010     format(' BLOWUP: impossible to reach VOL, increase OMMAX.$')
!
      end subroutine BLOWUP
!*==delsqf.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      real*8 function DELSQF(NDATA, IX, YA, X)
!- exponential least square fit (  a * exp(bx) )
! ----------------------------------------------------------------------
!i Inputs:
!i   ndata :number of data points
!i   ix    :xa(1) is integer and equals ix
!i          and xa(1),xa(2),xa(3),... differ exactly by 1.d0
!i   ya    :value of f at xa
!i   x     :x-value at which f is interpolated
!o Outputs:
!o   delsqf:interpolated value
!r Remarks:
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer IX, NDATA
          real*8 X
          real*8 YA(*)
!
! Local variables
!
          real*8 A, B, SS, SX, SXX, SXY, SY, XX, YY
          integer IDATA
          save A, B, IDATA, SS, SX, SXX, SXY, SY, XX, YY
!
!*** End of declarations rewritten by SPAG
!
          SX = 0.D0
          SY = 0.D0
          SXY = 0.D0
          SXX = 0.D0
!
          do IDATA = 1, NDATA
              XX = dble(IX + IDATA - 1)
              if (YA(IDATA) <= 0.D0) then
                  DELSQF = 0.D0
                  return
              end if
              YY = log(YA(IDATA))
              SX = SX + XX
              SY = SY + YY
              SXX = SXX + XX * XX
              SXY = SXY + XX * YY
          end do
          SS = dble(NDATA)
!
          B = (SXY - SX * SY / SS) / (SXX - SX * SX / SS)
          A = exp((SY - B * SX) / SS)
!
          DELSQF = A * exp(B * X)
!
      end function DELSQF
!*==drr2.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      real*8 function DRR2(PLAT, BAS1, BAS2, I, J, K, DR)
!- Calculates the vector connecting two sites in a solid
! ----------------------------------------------------------------
!i Inputs:
!i   plat  :primitive lattice vectors
!i   bas1  :basis vector of first site
!i   bas2  :basis vector of second site
!i   i,j,k :the number of primitive lattice vectors separating sites
!o Outputs:
!o   dr    :connecting vector bas2 - bas1
!o   drr2  :square of the length of this vector
!r Remarks:
!r   Using the TB package and a table of indices iax, the connecting
!r   vector and the square of the distance is obtained by:
!r   rsqr=drr2(plat,bas(1,iax(1)),bas(1,iax(2)),iax(3),iax(4),iax(5),dr)
! ----------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer I, J, K
          real*8 BAS1(3), BAS2(3), DR(3), PLAT(3, 3)
!
! Local variables
!
          integer IX
!
!*** End of declarations rewritten by SPAG
!
          DRR2 = 0.D0
          do IX = 1, 3
              DR(IX) = BAS2(IX) - BAS1(IX) + PLAT(IX, 1) * I + PLAT(IX, 2) &
                       * J + PLAT(IX, 3) * K
              DRR2 = DRR2 + DR(IX) * DR(IX)
          end do
      end function DRR2
!*==iclbas.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      integer function ICLBAS(IC, ICLASS, NBAS, NRBAS)
!- Returns an index to iclbas atom in basis given class
! ----------------------------------------------------------------------
!i Inputs:
!i   ic    :class label
!i   iclass:the jth atom belongs to class iclass(j)
!i   nbas  :number of atoms in the basis
!i   nrbas :the nrbas-th basis atom of class ic is seeked
!o Outputs:
!o   iclbas:the iclbas-th atom belongs to class ic
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer IC, NBAS, NRBAS
          integer ICLASS(*)
!
! Local variables
!
          integer IBAS, N
          character(72) MESSG
          save IBAS, MESSG, N
          external ERRMSG
!
!*** End of declarations rewritten by SPAG
!
          ICLBAS = 1
          N = 0
          do IBAS = 1, NBAS
              if (ICLASS(IBAS) == IC) N = N + 1
              if (N == NRBAS) then
                  ICLBAS = IBAS
                  return
              end if
          end do
          write (MESSG, 99001) IC, NRBAS, N
          call ERRMSG(MESSG, 4)
!
99001     format(' ICLBAS: class', i3, ', nrbas=', i3, ', but only', i3,        &
                         &        ' basis atoms exist.$')
      end function ICLBAS
!*==latlim.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      subroutine LATLIM(PQLAT, RMAXS, I1, I2, I3)
!- Set limits in X Y Z direction
! ----------------------------------------------------------------------
!i Inputs:
!i   pqlat :primitive lattice vectors (real or reciprocal space)
!i   rmaxs :maximum length of connecting vector
!o Outputs:
!o   i1,i2,i3:all connecting vectors lie within these multiples of the
!o            lattice vectors.
!r Remarks:
!r   This routine only returns the integer part
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer I1, I2, I3
          real*8 RMAXS
          real*8 PQLAT(3, 3)
!
! Local variables
!
          real*8 DET, QPLAT(3, 3), X1, X2, X3
          real*8 DLAMCH, DNRM23
          integer IPRINT
          save DET, QPLAT, X1, X2, X3
          external DINV33
          ! external DINV33, DLAMCH, DNRM23, IPRINT
!
!*** End of declarations rewritten by SPAG
!
          call DINV33(PQLAT, 1, QPLAT, DET)
          X1 = RMAXS * sqrt(DNRM23(QPLAT(1, 1)))
          X2 = RMAXS * sqrt(DNRM23(QPLAT(1, 2)))
          X3 = RMAXS * sqrt(DNRM23(QPLAT(1, 3)))
!
          I1 = 1 + int(X1 - DLAMCH('p'))
          I2 = 1 + int(X2 - DLAMCH('p'))
          I3 = 1 + int(X3 - DLAMCH('p'))
          if (IPRINT() >= 80) write (16, 99001) RMAXS, X1, X2, X3, I1, I2, &
              I3
!
99001     format(/' LATLIM: rmaxs=', f21.15, ', x1,x2,x3=', 3F21.15,          &
                         &        ', i1,i2,i3=', 3I2)
      end subroutine LATLIM
!*==defwsr.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine DEFWSR(WSR, Z)
!- Returns default value of radii for given nuclear charge
! ----------------------------------------------------------------------
!i Inputs:
!i   z     :nuclear charge
!o Outputs:
!o   wsr   :Wigner-Seitz sphere radius (in atomic units)
!r Remarks:
!r These are equilibrium average Wigner-Seitz radii for closed-packed
!r structure of pure element.
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 WSR, Z
!
! Local variables
!
          real*8 CPWSR(0:100)
          integer IZ
          character(72) MESSG
          save IZ, MESSG
          external ERRMSG
!
!*** End of declarations rewritten by SPAG
!
          data CPWSR/2.00D0, 1.39D0, 2.55D0, 3.04D0, 2.27D0, 1.96D0, &
              1.66D0, 1.90D0, 1.90D0, 2.17D0, 2.89D0, 3.76D0, 3.25D0, 2.95D0, &
              2.63D0, 2.56D0, 2.70D0, 2.85D0, 3.71D0, 4.66D0, 3.88D0, 3.31D0, &
              2.99D0, 2.76D0, 2.64D0, 2.57D0, 2.52D0, 2.52D0, 2.55D0, 2.62D0, &
              2.78D0, 2.75D0, 2.79D0, 2.83D0, 2.94D0, 3.13D0, 4.32D0, 4.95D0, &
              4.22D0, 3.61D0, 3.28D0, 3.03D0, 2.91D0, 2.82D0, 2.77D0, 2.78D0, &
              2.84D0, 2.95D0, 3.14D0, 3.30D0, 3.45D0, 3.30D0, 3.31D0, 3.50D0, &
              4.31D0, 5.30D0, 4.20D0, 3.91D0, 3.80D0, 3.75D0, 3.70D0, 3.65D0, &
              3.60D0, 3.55D0, 3.52D0, 3.61D0, 3.67D0, 3.70D0, 3.73D0, 3.75D0, &
              3.56D0, 3.44D0, 3.23D0, 3.04D0, 2.93D0, 2.86D0, 2.82D0, 2.83D0, &
              2.88D0, 2.98D0, 3.27D0, 3.57D0, 3.62D0, 3.37D0, 3.46D0, 3.63D0, &
              4.44D0, 5.81D0, 4.30D0, 3.84D0, 3.52D0, 3.32D0, 3.13D0, 3.02D0, &
              2.96D0, 2.93D0, 2.93D0, 2.95D0, 2.99D0, 3.05D0, 3.17D0/
!
          IZ = IDNINT(Z)
          if (IZ < 0 .or. IZ > 100) then
              write (MESSG, 99001) IZ
              call ERRMSG(MESSG, 0)
              WSR = CPWSR(100)
          else
              WSR = CPWSR(IZ)
          end if
99001     format('DEFWSR: bad nuclear charge, Z=', i3)
      end subroutine DEFWSR
