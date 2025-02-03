!*==dneville.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      function DNEVILLE(X, Y, T, X1, N)
!
!     Neville -Algorithmus
!
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer N
          real*8 X1
          real*8 DNEVILLE
          real*8 T(N), X(N), Y(N)
!
! Local variables
!
          integer I, J
!
!*** End of declarations rewritten by SPAG
!
          do I = 1, N
              T(I) = Y(I)
          end do
!
! It is not clear for me, what does it mean t(0)
!
          do I = 1, N - 1
              ! <=== My changes!!!
              do J = N, I + 1, -1
                  T(J) = T(J) + (X1 - X(J)) * (T(J) - T(J - 1)) / (X(J) &
                                                                   - X(J - I))
                  ! j
              end do
              ! i
          end do
          DNEVILLE = T(N)
      end function DNEVILLE
!*==dnevmod.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      function DNEVMOD(X, Y, X1, N, M)
!
!     Neville -Algorithmus (modified)
!
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          integer MMAX
          parameter(MMAX=20)
!
! Dummy arguments
!
          integer M, N
          real*8 X1
          real*8 DNEVMOD
          real*8 X(N), Y(N)
!
! Local variables
!
          real*8 DNEVILLE
          integer J, JJ, MM
          real*8 T(MMAX)
!
!*** End of declarations rewritten by SPAG
!
          if (M > MMAX) stop ' DNEVMOD: m>mmax'
          call BISEC(J, N, X1, X)
          MM = min(M, 2 * J, 2 * (N - J))
          MM = max(MM, 8)
          JJ = min(max(J - MM / 2, 1), N - MM + 1)
!     print*,'JJ=',jj,x(jj),' N=',N
          DNEVMOD = DNEVILLE(X(JJ), Y(JJ), T, X1, MM)
      end function DNEVMOD
!*==bisec.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine BISEC(K, N, R, RMESH)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer K, N
          real*8 R
          real*8 RMESH(*)
!
! Local variables
!
          real*8 HR
          integer IB, IE, IH
!
!*** End of declarations rewritten by SPAG
!
          K = 0
          if (R < RMESH(1) .or. R > RMESH(N)) return
!
          IB = 1
          IE = N
!
100       continue
          IH = IB + (IE - IB) / 2
!
          HR = RMESH(IH)
          if (IE - IB == 1) then
              K = IB
              return
          end if
          if (R > HR) then
              IB = IH
              goto 100
          elseif (R < HR) then
              IE = IH
              goto 100
          end if
          K = IH
      end subroutine BISEC
!*==qd.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine QD(H, Y, Z, NDIM)
!
!       THE SUBROUTINE CALCULATES NDIM INTEGRALS
!       BETWEEN X1 AND XN OF FUNCTION Z ( N.LE.NDIM
!        XN=X1+(N-1)*H ) FOR THE SPECIAL CASE:
!       NDIM=2*ND+1 Y(1)=0 . Z(2*I+1) IS CALCULATED
!       BY THE SIMPSON'S METHOD AND Z(2*I+2)=Z(2*I+1)+
!       DELT(2*I+2) . TO INCREASE A ACCURACY IT IS
!       USED DELT(2*I+2)=DELT4(2*I+2)*(Z(2*I+3)-Z(2*I+1))/
!       (DELT4(2*I+2)+DELT4(2*I+3)) WHERE DELT4(2*I+2)
!       IS CALCULATED BY USING A CUBIC INTERPOLATION
!       BETWEEN Y(2*I),Y(2*I+1),Y(2*I+2),Y(2*I+3).
!
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 H
          integer NDIM
          real*8 Y(*), Z(*)
!
! Local variables
!
          real*8 C24, C3, DELT, DELT3, DELT4, DELT5, SUM2, Y1, Y2, &
              Y3, Y4, Y5
          integer I
!
!*** End of declarations rewritten by SPAG
!
          if ((NDIM / 2 + NDIM / 2) == NDIM) stop 'THE SUBPROGRAM QD WORKS&
         & ONLY IF NDIM=2*ND+1'
!
          C3 = H / 3.D0
          C24 = H / 24.D0
          Z(1) = 0.D0
          if (NDIM /= 1) then
              Y2 = Y(1)
              Y3 = Y(2)
              Y4 = Y(3)
              Z(2) = C24 * (10 * Y2 + 16 * Y3 - 2 * Y4)
              SUM2 = Y2 + 4 * Y3 + Y4
              Z(3) = C3 * SUM2
              if (NDIM /= 3) then
                  Y5 = Y(4)
!
! the main loop of a integration
!
                  do I = 4, NDIM - 1, 2
                      Y1 = Y3
                      Y2 = Y4
                      Y3 = Y5
                      Y4 = Y(I + 1)
                      DELT3 = Y2 + 4 * Y3 + Y4
                      SUM2 = SUM2 + DELT3
                      Z(I + 1) = C3 * SUM2
                      DELT4 = -Y1 + 13 * Y2 + 13 * Y3 - Y4
                      if (I < NDIM - 1) then
                          Y5 = Y(I + 2)
                          DELT5 = -Y1 + 12 * Y2 + 26 * Y3 + 12 * Y4 - Y5
                      else
                          DELT5 = DELT3 * 8
                      end if
                      if (DELT5 > 1.D-20) then
                          DELT = 8 * DELT4 * DELT3 / DELT5
                          Z(I) = Z(I - 1) + C24 * DELT
                      else
                          Z(I) = Z(I - 1) + C24 * DELT4
                      end if
                  end do
              end if
          end if
      end subroutine QD
!*==intlog.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine INTLOG(R1, DX, Y, Y0, S, YINT, RI)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 DX, R1, S, Y0, YINT
          real*8 RI(*), Y(*)
!
! Local variables
!
          real*8 A, ALF, BET, DEL, F0, F1, F2, FM1, FM2, FSR1, GAM, &
              P, PNMA, PS, RALF, SS, SUM, SUM1, SUM2, X0, XNMA, XS
          real*8 FI
          integer I, I0, I1, I2, IM1, IMA, N, NMA
!
!*** End of declarations rewritten by SPAG
!
          FI(P) = F0 * P + .5D0 * ALF * P**2 + BET / 3.D0 * P**3 + GAM / &
                  4.D0 * P**4 + DEL / 5.D0 * P**5
!  PERFORMS THE INTEGRAL OF Y IN THE INTERVAL (0,S)
!  USING SIMPSON'S RULE
! LOGARITHMIC NET RI(I)=R1*EXP((I-1)*DX)
!  IF (Y0.EQ.0.) Y(R)=A*R**ALF IS ASSUMED AND CALCULATED FROM
!  THE FIRST 2 NET POINTS
!     WRITE(6,*) S,R1,(Y(I),I=1,10)
          XS = log(S / R1)
          SS = 1.D0 + XS / DX + .00001D0
          I2 = 1 + SS
          I1 = I2 - 1
          I0 = I1 - 1
          IM1 = I1 - 2
          N = I1 / 4
          NMA = 4 * N + 1
          IMA = NMA - 3
          SUM = (Y(1) * R1 - Y(NMA) * RI(NMA)) * 7.D0
!      print*,'START'
          do I = 2, IMA, 4
              SUM = SUM + 32.D0 * Y(I) * RI(I) + 12.D0 * Y(I + 1) * RI(I + 1) + &
                    32.D0 * Y(I + 2) * RI(I + 2) + 14.D0 * Y(I + 3) * RI(I + 3)
          end do
          SUM = SUM * DX * 2.D0 / 45.D0
!      print*,'CONT'
          if (Y0 == 0.) then
              SUM1 = 0.D0
              if (Y(1) /= 0.) then
                  ALF = log(Y(2) / Y(1)) / DX
                  RALF = R1**ALF
                  A = Y(1) / RALF
                  SUM1 = A / (ALF + 1.D0) * RALF * R1
              end if
          else
              FSR1 = 1.D0 / DX * (-11.D0 * Y(1) + 18.D0 * Y(2) - 9.D0 * &
                                  Y(3) + 2.D0 * Y(4)) / 6.D0 / RI(1)
              SUM1 = RI(1) * (Y0 + Y(1))*.5D0 - RI(1)**2 / 12.D0 * &
                     FSR1
          end if
!  19 FORMAT(' SUB. INTLOG: Y(1),Y(2),A,ALF,SUM1 =',1P5E21.10)
!     PRINT 19,Y(1),Y(2),A,ALF,SUM1
!      print*,'C33',I2,Y(I2),RI(I2)
          F2 = Y(I2) * RI(I2)
!      print*,'C--'
          F1 = Y(I1) * RI(I1)
          F0 = Y(I0) * RI(I0)
!      print*,'C--'
!
          FM1 = Y(IM1) * RI(IM1)
          FM2 = Y(IM1 - 1) * RI(IM1 - 1)
          ALF = (2.D0 * FM2 - 16.D0 * FM1 + 16.D0 * F1 - 2.D0 * F2) / 24.D0
          BET = (-FM2 + 16.D0 * FM1 - 30.D0 * F0 + 16.D0 * F1 - F2) &
                / 24.D0
          GAM = (-2.D0 * FM2 + 4.D0 * FM1 - 4.D0 * F1 + 2.D0 * F2) &
                / 24.D0
          DEL = (FM2 - 4.D0 * FM1 + 6.D0 * F0 - 4.D0 * F1 + F2) / 24.D0
          XNMA = (NMA - 1) * DX
          X0 = (I0 - 1) * DX
          PNMA = (XNMA - X0) / DX
          PS = (XS - X0) / DX
          SUM2 = (FI(PS) - FI(PNMA)) * DX
          YINT = SUM + SUM1 + SUM2
!      print*,'C44'
!     PRINT 101,XS,I2,NMA
!  101 FORMAT(' XS,I2,NMA : ',F21.15,2I5)
!  102 FORMAT(' RADII :',//,50(/5E21.10))
!  103 FORMAT(' INTEGRAND :',//,50(/5E21.10))
!     PRINT 102,RI
!     PRINT 103,Y
!     PRINT 104,F0,ALF,BET,GAM
!  104 FORMAT(' F0,ALF,BET,GAM,DEL : ',5E21.10)
!  105 FORMAT(' PNMA,PS ;',2E21.10)
!     PRINT 105,PNMA,PS
!     WRITE(6,100) SUM,SUM1,SUM2,YINT
!  100 FORMAT('SUM,SUM1,SUM2,YINT:',4D18.10)
      end subroutine INTLOG
!*==au2a.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      real*8 function AU2A()
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          real*8 A0
          parameter(A0=0.529177D0)
!
!*** End of declarations rewritten by SPAG
!
          AU2A = A0
!
      end function AU2A
!*==incopy.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine INCOPY(N, X, INCX, Y, INCY)
!- Copies an integer vector, x, to an integer vector, y.
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :lenght of x
!i   x     :vector to copy
!i   incx  :incrementation for x
!i   incy  :incrementation for y
!o Outputs:
!o   y     :result vector is stored in y
!r Remarks:
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer INCX, INCY, N
          integer X(*), Y(*)
!
! Local variables
!
          integer I, IX, IY
!
!*** End of declarations rewritten by SPAG
!
! Passed parameters:
! Local parameters:
!
          if (N <= 0) return
          if (INCX /= 1 .or. INCY /= 1) then
! ----  code for unequal increments or equal increments not equal to 1
              IX = 1
              IY = 1
              if (INCX < 0) IX = (-N + 1) * INCX + 1
              if (INCY < 0) IY = (-N + 1) * INCY + 1
              do I = 1, N
                  Y(IY) = X(IX)
                  IX = IX + INCX
                  IY = IY + INCY
              end do
          else
! ----  code for both increments equal to 1
              do I = 1, N
                  Y(I) = X(I)
              end do
          end if
      end subroutine INCOPY
!*==errmsg.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
! ------ ERRMSG ------
      subroutine ERRMSG(TEXT, I)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer I
          character(*) TEXT
!
! Local variables
!
          character(12) FATAL, WARN
!
!*** End of declarations rewritten by SPAG
!
! i<0 info
! i=0 warning
! i>0 fatal error
          data FATAL/'Fatal error:'/, WARN/'Warning:'/
!
          if (I <= 0) then
              print 99001, WARN, TEXT
!        if(iun.ne.6)write(iun,1000)warn,text
          else
              print 99001, FATAL, TEXT
!        if(iun.ne.6)write(iun,1000)fatal,text
              call ENDJOB(-11, TEXT)
          end if
99001     format(/2X, a12, 1X, a)
      end subroutine ERRMSG
!*==infomsg.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
! ------ INFOMSG ------
      subroutine INFOMSG(TEXT)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          character(*) TEXT
!
!*** End of declarations rewritten by SPAG
!
!
          print 99001, TEXT
!      if(iun.ne.6)write(iun,1000)text
99001     format(2X, a)
      end subroutine INFOMSG
!*==warning.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
! ------ WARNING ------
      subroutine WARNING(A)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          character(*) A
!
!*** End of declarations rewritten by SPAG
!
          print '(1x,a)', A
      end subroutine WARNING
!*==ddet33.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      real*8 function DDET33(A)
!- Determinant of 3x3 matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   a :input matrix
!o Outputs:
!o   ddet33 :determinant
! ----------------------------------------------------------------------
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 A(*)
!
! Local variables
!
          real*8 V(3)
!
!*** End of declarations rewritten by SPAG
!
! Passed parameters:
! Local parameters:
!
          call CROSS(V, A(4), A(7))
          DDET33 = A(1) * V(1) + A(2) * V(2) + A(3) * V(3)
      end function DDET33
!*==trnt.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      real*8 function TRNT(A, B, C)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 A(3), B(3), C(3)
!
!*** End of declarations rewritten by SPAG
!
          TRNT = A(1) * B(2) * C(3) + A(2) * B(3) * C(1) + A(3) &
                 * B(1) * C(2) - A(3) * B(2) * C(1) - A(2) * B(1) * C(3) &
                 - A(1) * B(3) * C(2)
      end function TRNT
!*==eilang.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine EILANG(ANG, IG)
!
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          real*8 P, Z, H, H3, T, T2, T4, T5
          parameter(P=3.141592653589793D0, Z=0.D0, H=P / 2, H3=3 * &
                    H, T=P / 3, T2=T * 2, T4=T * 4, T5=T * 5)
!
! Dummy arguments
!
          integer IG
          real*8 ANG(3)
!
! Local variables
!
          real*8 C(3, 32)
!
!*** End of declarations rewritten by SPAG
!
          data C/Z, Z, Z, Z, P, Z, Z, P, P, P, Z, Z, P, H, H, Z, H, H3, P, &
              H, H3, Z, H, H, H, H, Z, H, H, P, H3, H, P, H3, H, Z, H3, P, Z, H, &
              Z, Z, H3, Z, Z, H, P, Z, Z, H, P, P, H, Z, Z, H, Z, P, H, P, H3, &
              H, H3, H3, H, H, H, H, H, H, H, H3, T, Z, Z, T2, Z, Z, T4, Z, Z, &
              T5, Z, Z, T4, P, Z, T5, P, Z, T, P, Z, T2, P, Z/
          ANG(1) = C(1, IG)
          ANG(2) = C(2, IG)
          ANG(3) = C(3, IG)
      end subroutine EILANG
!*==turnm.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine TURNM(OM, U)
!      INCLUDE 'PREC.H'
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 OM(3), U(3, 3)
!
! Local variables
!
          real*8 COSA, COSB, COSG, SINA, SINB, SING
!
!*** End of declarations rewritten by SPAG
!
          SINA = sin(OM(1))
          SINB = sin(OM(2))
          SING = sin(OM(3))
          COSA = cos(OM(1))
          COSB = cos(OM(2))
          COSG = cos(OM(3))
          U(1, 1) = COSA * COSB * COSG - SINA * SING
          U(1, 2) = -COSA * COSB * SING - SINA * COSG
          U(1, 3) = COSA * SINB
          U(2, 1) = SINA * COSB * COSG + COSA * SING
          U(2, 2) = -SINA * COSB * SING + COSA * COSG
          U(2, 3) = SINA * SINB
          U(3, 1) = -SINB * COSG
          U(3, 2) = SINB * SING
          U(3, 3) = COSB
      end subroutine TURNM
!*==dnrm23.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      real*8 function DNRM23(R)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 R(3)
!
!*** End of declarations rewritten by SPAG
!
!-  'norm'-square
!
          DNRM23 = R(1) * R(1) + R(2) * R(2) + R(3) * R(3)
!
      end function DNRM23
!*==lengths.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      integer function LENGTHS(A)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          character(*) A
!
! Local variables
!
          integer J, L
!
!*** End of declarations rewritten by SPAG
!
! actual length of the string
          L = len(A)
          J = 0
          do while (J < L .and. A(J + 1:) /= ' ')
              J = J + 1
          end do
          LENGTHS = J
      end function LENGTHS
!*==dpi.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
      real*8 function DPI()
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
          logical FIRST
          real*8 PI
          save PI
!
!*** End of declarations rewritten by SPAG
!
          data FIRST/.true./
!
          if (FIRST) then
              FIRST = .false.
              PI = atan(1.D0) * 4.D0
!$$$        print *,pi
          end if
          DPI = PI
      end function DPI
!*==dinv33.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine DINV33(MATRIX, IOPT, INVRSE, DET)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 DET
          integer IOPT
          real*8 INVRSE(3, 3), MATRIX(3, 3)
!
! Local variables
!
          real*8 DDOT
          integer I, J
          real*8 XX
!
!*** End of declarations rewritten by SPAG
!
!- INVERTS 3X3 MATRIX
! ----------------------------------------------------------------
!I INPUTS
!I   INVERSE: INPUT MATRIX
!I   IOPT:  IF 0, USUAL INVERSE
!I             1, TRANSPOSE OF INVERSE
!I             2, 2*PI*INVERSE
!I             3, 2*PI*TRANSPOSE OF INVERSE
!O OUTPUTS
!O   INVERSE, AS MODIFIED ACCORDING TO IOPT
!O   DET:      DETERMINANT, OR DET/2*PI (SIGN OK ??)
!R REMARKS
!R   TO GENERATE RECIPROCAL LATTICE VECTORS, CALL DINV33(PLAT,3,RLAT)
! ----------------------------------------------------------------
!     IMPLICIT NONE
          call CROSS(INVRSE, MATRIX(1, 2), MATRIX(1, 3))
          call CROSS(INVRSE(1, 2), MATRIX(1, 3), MATRIX)
          call CROSS(INVRSE(1, 3), MATRIX, MATRIX(1, 2))
          DET = DDOT(3, MATRIX, 1, INVRSE, 1)
          if (DET == 0.D0) call ENDJOB(10, 'INV33: VANISHING DETERMINANT')
          if (IOPT >= 2) DET = DET / (8 * DATAN(1.D0))
          if (mod(IOPT, 2) == 0) then
              do I = 1, 3
                  do J = I + 1, 3
                      XX = INVRSE(I, J)
                      INVRSE(I, J) = INVRSE(J, I)
                      INVRSE(J, I) = XX
                  end do
              end do
          end if
          call DSCAL(9, 1.D0 / DET, INVRSE, 1)
      end subroutine DINV33
!*==cross.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine CROSS(A, B, C)
!
!  cross product (ax,ay,az)=(bx,by,bz)*(cx,cy,cz)
!
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 A(3), B(3), C(3)
!
!*** End of declarations rewritten by SPAG
!
!
          A(1) = B(2) * C(3) - B(3) * C(2)
          A(2) = B(3) * C(1) - B(1) * C(3)
          A(3) = B(1) * C(2) - B(2) * C(1)
!
      end subroutine CROSS
!*==dscal.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine DSCAL(N, DA, DX, INCX)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 DA
          integer INCX, N
          real*8 DX(*)
!
! Local variables
!
          integer I, M, MP1, NINCX
!
!*** End of declarations rewritten by SPAG
!
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
          if (N <= 0 .or. INCX <= 0) return
          if (INCX == 1) then
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
              M = mod(N, 5)
              if (M /= 0) then
                  do I = 1, M
                      DX(I) = DA * DX(I)
                  end do
                  if (N < 5) return
              end if
              MP1 = M + 1
              do I = MP1, N, 5
                  DX(I) = DA * DX(I)
                  DX(I + 1) = DA * DX(I + 1)
                  DX(I + 2) = DA * DX(I + 2)
                  DX(I + 3) = DA * DX(I + 3)
                  DX(I + 4) = DA * DX(I + 4)
              end do
          else
!
!        code for increment not equal to 1
!
              NINCX = N * INCX
              do I = 1, NINCX, INCX
                  DX(I) = DA * DX(I)
              end do
              return
          end if
      end subroutine DSCAL
!*==ddot.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      real*8 function DDOT(N, DX, INCX, DY, INCY)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer INCX, INCY, N
          real*8 DX(*), DY(*)
!
! Local variables
!
          real*8 DTEMP
          integer I, IX, IY, M, MP1
!
!*** End of declarations rewritten by SPAG
!
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
          DDOT = 0.0D0
          DTEMP = 0.0D0
          if (N <= 0) return
          if (INCX == 1 .and. INCY == 1) then
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
              M = mod(N, 5)
              if (M /= 0) then
                  do I = 1, M
                      DTEMP = DTEMP + DX(I) * DY(I)
                  end do
                  if (N < 5) then
                      DDOT = DTEMP
                      goto 99999
                  end if
              end if
              MP1 = M + 1
              do I = MP1, N, 5
                  DTEMP = DTEMP + DX(I) * DY(I) + DX(I + 1) * DY(I + 1) &
                          + DX(I + 2) * DY(I + 2) + DX(I + 3) * DY(I + 3) + DX(I + &
                                                                               4) * DY(I + 4)
              end do
              DDOT = DTEMP
          else
!
!        code for unequal increments or equal increments
!          not equal to 1
!
              IX = 1
              IY = 1
              if (INCX < 0) IX = (-N + 1) * INCX + 1
              if (INCY < 0) IY = (-N + 1) * INCY + 1
              do I = 1, N
                  DTEMP = DTEMP + DX(IX) * DY(IY)
                  IX = IX + INCX
                  IY = IY + INCY
              end do
              DDOT = DTEMP
              return
          end if
99999     continue
      end function DDOT
!*==endjob.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine ENDJOB(J, TEXT)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer J
          character(*) TEXT
!
!*** End of declarations rewritten by SPAG
!
!-----------------
          write (*, *) ' === STOP === ', J, ' ', TEXT
          stop
      end subroutine ENDJOB
!*==wkinit.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine WKINIT(NSIZE)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! COMMON variables
!
          integer LIMIT, OFREE
          real*8 QDUMMY8
          common / Q_LMTO_ / QDUMMY8, OFREE, LIMIT
!
! Dummy arguments
!
          integer LENG, NSIZE, ONAME
!
! Local variables
!
          integer I8, LENGTH, L_WORD_C, L_WORD_CC, L_WORD_I, L_WORD_R, &
              L_WORD_RR
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
          data L_WORD_I/1/, L_WORD_R/1/, L_WORD_RR/2/, L_WORD_C/ &
              2/, L_WORD_CC/4/
!
! ----- DEFINE STORAGE SIZE ------
!  START OF FIRST ARRAY AND MAX NUMBER TO BE DEFINED:
          LIMIT = NSIZE
          OFREE = 5
          return
! ------ SUBROUTINES TO DEFINE ARRAYS OF VARIOUS TYPES -----
          entry DEFI(ONAME, LENG)
          LENGTH = LENG * L_WORD_I
          I8 = 0
          goto 100
          entry DEFR(ONAME, LENG)
          LENGTH = LENG * L_WORD_R
          I8 = 0
          goto 100
          entry DEFC(ONAME, LENG)
          LENGTH = LENG * L_WORD_C
          I8 = 0
          goto 100
          entry DEFRR(ONAME, LENG)
          LENGTH = LENG * L_WORD_RR
          I8 = 1
          goto 100
          entry DEFDR(ONAME, LENG)
          LENGTH = LENG * L_WORD_RR
          I8 = 1
          goto 100
          entry DEFCC(ONAME, LENG)
          LENGTH = LENG * L_WORD_CC
          I8 = 1
100       continue
          if (LENGTH < 0) call ENDJOB(10, '**** LENGTH OF ARRAY NEGATIVE')
          if (LENGTH == 0) LENGTH = 1
          if (I8 == 1 .and. mod(OFREE, 2) == 0) OFREE = OFREE + 1
          ONAME = OFREE
          OFREE = OFREE + LENGTH
          if (OFREE >= LIMIT) then
              print 99001, 'def0', OFREE, LIMIT
              call ENDJOB(5, ' NO MEMORY')
          end if
          return
!
          entry RLSE(ONAME)
          if (ONAME > LIMIT) call ENDJOB(10, '**** RESET POINTER GT LIMIT'&
         &)
          if (ONAME < 3) call ENDJOB(10, '**** RESET POINTER LT 3')
          OFREE = ONAME
99001     format(/45('*')/'  In ', a, ' memory pool exhausted', /'  N___=', I7,&
                         &        ' NMAX_= ', i7, /45('*'))
      end subroutine WKINIT
!*==rhfds.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
! ! 03-03-94 - PERLOV : NDP cann'n calculate (0.0)**(1./3.)
! !                         changed by (0.0 + 1.d-42)**(1./3.)
! **********************************************************************
!
! HARTREE FOCK DIRAC SLATER          J P DESCLAUX     CEA PARIS 1969
! MODIFIE JUILLET 1970
!
! **********************************************************************
      subroutine RHFDS(IZ, IZCORE, RO, ROM, DR00)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! COMMON variables
!
          character(4) BAR(10), TITRE(30)
          real*8 D(251), DC(251), DCOP, DEN(30), DEXE, DEXV, &
              DFL(30), DGC(251, 30), DP(251), DPAS, DPC(251, 30), DQ(251), &
              DQ1(30), DR(251), DV(251), DVF(251), DVN(251), Z
          integer ICUT, ION, IPRAT, MAG(30), NCORB(30), NEL(30), NES, &
              NITER, NK(30), NMAX(30), NORB, NP, NQL(30), NQN(30), NSTOP, &
              NUC, NWF(30)
          real(8) TEST, TESTE, TESTV, TESTY, TETS
          common DEN, DQ1, DFL, NQN, NQL, NK, NMAX, NEL, NORB, MAG
          common / DEUX / DVN, DVF, D, DC, DGC, DPC
          common / DIRA / DV, DR, DP, DQ, DPAS, Z, NSTOP, NES, TETS, NP, &
              NUC
          common / PS2 / DEXV, DEXE, DCOP, TITRE, BAR, TEST, TESTE, TESTY, &
              TESTV, NITER, ION, ICUT, IPRAT
          common / UKP / NCORB, NWF
!
! Dummy arguments
!
          real*8 DR00
          integer IZ, IZCORE
          real*8 RO(251), ROM(251)
!
! Local variables
!
          real*8 DE, DENER, DVAL
          real(8) EMAX, VAL, VMAX, Y, YMAX, YN
          integer I, ICEL, IM, IMAX, ITER, J, K, KK, L, N, NPCH
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
          NSTOP = 1
100       continue
          if (IZ >= 0) then
              call INSLD(NPCH, IZ, IZCORE)
          else
              call INSLD(NPCH, IZ, IZCORE)
          end if
          ITER = 1
          do I = 1, NP
              do J = 1, NORB
                  DGC(I, J) = 0.D0
                  DPC(I, J) = 0.D0
              end do
          end do
          write (8, 99001) BAR
          write (8, 99002)
          N = -(ION + 1)
200       continue
          do I = 1, NP
              D(I) = 0.D0
          end do
          TETS = TEST
          YMAX = 0.
          VMAX = 0.
          EMAX = 0.
! PE{EHiE uPABHEHiq diPAKA dlq KAvdOj OPbiTAli
! **********************************************************************
          do J = 1, NORB
              DE = DEN(J)
250           continue
              call RESLD(NQN(J), NQL(J), NK(J), IMAX, DEN(J), DFL(J), &
                         DQ1(J), J)
              if (NSTOP == 0) then
                  VAL = abs((DEN(J) - DE) / DE)
                  if (VAL > EMAX) EMAX = VAL
                  NMAX(J) = IMAX
                  do I = 1, NP
                      VAL = DGC(I, J) - DP(I)
                      if (abs(DP(I)) > 1.D0) VAL = VAL / DP(I)
                      if (abs(VAL) >= abs(YMAX)) then
                          YMAX = VAL
                          Y = DP(I)
                          YN = DGC(I, J)
                      end if
                      VAL = DPC(I, J) - DQ(I)
                      if (abs(DQ(I)) > 1.D0) VAL = VAL / DQ(I)
                      if (abs(VAL) >= abs(YMAX)) then
                          YMAX = VAL
                          Y = DQ(I)
                          YN = DPC(I, J)
                      end if
                      DGC(I, J) = DP(I)
                      DPC(I, J) = DQ(I)
                      D(I) = D(I) + NEL(J) * (DP(I) * DP(I) + DQ(I) * DQ(I))
                  end do
              elseif (NSTOP /= 362 .or. ITER >= 10 .or. TETS > TEST) then
                  write (8, 99008) NSTOP, NQN(J), TITRE(J)
                  goto 100
              else
                  TETS = TESTV
                  goto 250
              end if
          end do
          call POTSL(DC, D, DP, DR, DPAS, DEXV, Z, NP, ION, ICUT)
          if (NUC > 0) then
              do I = 1, NUC
                  DC(I) = DC(I) + Z / DR(I) + Z * ((DR(I) / DR(NUC))**2 &
                                                   - 3.D0) / (DR(NUC) + DR(NUC))
              end do
          end if
          do I = 1, NP
              DVAL = abs(DC(I) - DV(I))
              if ((DR(I) * DC(I)) <= N) DVAL = -DVAL / DC(I)
              if (DVAL > VMAX) then
                  VMAX = DVAL
                  J = I
              end if
          end do
          write (8, 99003) ITER, VMAX, DR(J), DV(J), DC(J), EMAX, YMAX, &
              YN, Y
!$$$      call textout(20,20,0,7,40,buf)
          if ( &
              TETS > TEST .or. EMAX > TESTE .or. VMAX > TESTV .or. YMAX > TESTY) &
              then
              ITER = ITER + 1
              if (ITER <= NITER) then
! pOTEHciAl dlq ClEdu`}Ej iTEPAcii
! **********************************************************************
                  DVAL = 1.D0 - DCOP
                  do I = 1, NP
                      DVN(I) = DV(I)
                      DVF(I) = DC(I)
                      DV(I) = DVAL * DV(I) + DCOP * DC(I)
                  end do
                  goto 200
              else
                  write (8, 99004) NITER
                  NSTOP = 2
              end if
          end if
          if (NPCH /= 0) then
              write (7, 99009) BAR
              write (7, 99010) (DR(KK), D(KK), KK=1, 251)
          end if
          write (16, 99005)
          write (8, 99005)
! CPEdHiE zHA~EHiq R
! **********************************************************************
          do I = 1, NP
              DVF(I) = DC(I)
              DQ(I) = 0.
          end do
          DVAL = 0.
          do I = 1, NORB
              IM = NMAX(I)
              DVAL = DVAL + NEL(I) * DEN(I)
              do J = 1, IM
                  DC(J) = DGC(J, I) * DGC(J, I) + DPC(J, I) * DPC(J, I)
              end do
              L = 5
              if (IABS(NK(I)) == 1) L = L - 1
              do J = 1, L
                  DP(J) = DFL(I) + DFL(I)
                  if (J == 1) then
                      N = 4
                  elseif (J == 2) then
                      N = 2
                  elseif (J == 3) then
                      N = 1
                  elseif (J == 4) then
                      N = -1
                  else
                      N = -3
                  end if
                  call SOMM(DR, DC, DQ, DPAS, DP(J), N, IM)
              end do
              write (8, 99006) NQN(I), TITRE(I), DEN(I) * 2, NEL(I), &
                  NCORB(I), (DP(J), J=2, 3)
              write (16, 99006) NQN(I), TITRE(I), DEN(I) * 2, NEL(I), &
                  NCORB(I), (DP(J), J=2, 3)
          end do
! pOlHAq |HEPgiq uCPEdHEHHAq pO CfEPE
! **********************************************************************
          DC(1) = 1.D0
          do I = 1, NP
              DP(I) = D(I) / DR(I)
          end do
          if (NUC > 0) then
              do I = 1, NUC
                  DP(I) = D(I) * (3.D0 - DR(I) * DR(I) / (DR(NUC) * DR( &
                                                          NUC))) / (DR(NUC) + DR(NUC))
              end do
              DC(1) = 4.D0
          end if
          call SOMM(DR, DP, DQ, DPAS, DC(1), 0, NP)
          do I = 1, NP
              DP(I) = D(I) * DVF(I)
              D(I) = D(I) * ((D(I) * DR(I) + 1.D-42)**(1.D0 / 3.D0))
          end do
          DC(2) = 3.D0
          DC(3) = 1.D0
          if (NUC /= 0) DC(3) = 4.D0
          call SOMM(DR, DP, DQ, DPAS, DC(3), 0, NP)
          call SOMM(DR, D, DQ, DPAS, DC(2), -1, NP)
          DC(2) = -3.D0 * DC(2) / (105.27578D0**(1.D0 / 3.D0))
          DC(1) = -Z * DC(1)
          DC(4) = DVAL - DC(3)
          DVAL = DVAL + (DC(1) - DC(3) + (DEXE - DEXV) * DC(2)) / 2.D0
          DC(3) = (DC(3) - DC(1) - DEXV * DC(2)) / 2.D0
          DC(2) = DC(2) * DEXE / 2.D0
!     WRITE(8,1007) DVAL,DC(4),DC(3),DC(2),DC(1)
          if (NORB /= 1) then
!     WRITE(8,1001) BAR
!     WRITE(8,1008)
! iHTEgPAly pEPEKPyTiq
! **********************************************************************
              do I = 2, NORB
                  K = I - 1
                  do J = 1, K
                      if (NQL(I) == NQL(J) .and. NK(I) == NK(J)) then
                          IM = NMAX(J)
                          if (NMAX(I) < IM) IM = NMAX(I)
                          do L = 1, IM
                              DQ(L) = DPC(L, I) * DPC(L, J)
                              DC(L) = DGC(L, I) * DGC(L, J)
                          end do
                          DVAL = DFL(I) + DFL(J)
                          call SOMM(DR, DC, DQ, DPAS, DVAL, 0, IM)
                          if (DVAL > 1.D-3) write (8, 99007) NQN(I), TITRE(I), &
                              NQN(J), TITRE(J), DVAL
                      end if
                  end do
              end do
          end if
!     WRITE(8,4003)
          do J = 1, NORB
              if (NWF(J) /= 0) then
                  DENER = DEN(J) * 2
                  write (41, 99011) NQN(J), NQL(J), DENER, DR(1)
                  write (41, 99012) (DGC(I, J), I=1, 250)
              end if
          end do
!
          ICEL = 0
          do J = 1, NORB
              if (NCORB(J) /= 0) ICEL = ICEL + NEL(J)
          end do
!
          DR00 = DR(1)
          do I = 1, 251
              RO(I) = 0.D0
              ROM(I) = 0.D0
          end do
          do J = 1, NORB
!
              do I = 1, 251
                  RO(I) = RO(I) + NEL(J) * (DGC(I, J) * DGC(I, J) + DPC(I, J) &
                                            * DPC(I, J))
                  ROM(I) = ROM(I) + MAG(J) * (DGC(I, J) * DGC(I, J) + DPC(I, &
                                                                          J) * DPC(I, J))
              end do
          end do
!
99001     format(5X, 10A4/)
99002     format(' ITER', 4X, 'DVMAX', 10X, 'R', 14X, 'VN-1', 13X, 'VN', 10X,       &
                         &        'DEMAX', 6X, 'DPMAX', 9X, 'PN-1', 13X, 'PN')
99003     format(I5, 1PE21.10, 3(E21.10), 2(1PE21.10), 2(E21.10))
99004     format(' Number of iterations is exceeded', i4)
99005     format('Level  Energy      Occup. Val.      R**2          R')
99006     format(I2, A2, 1X, E21.10, 2X, I2, 3X, I2, 2X, E21.10, 2X, E21.10)
99007     format(34X, I1, A2, I3, A2, F21.15)
99008     format('  NSTOP=', I4, '  dlq OPbiTAli   ', I3, A2)
99009     format(10A4)
99010     format(4D15.8)
99011     format(I3, 1X, I3, 1X, G14.7, 1X, G14.7)
99012     format(5(E21.10))
      end subroutine RHFDS
!*==insld.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine INSLD(NPCH, IZ, IZCORE)
!
! ~TEHiE CTAPTOBOgO pOTEHciAlA
! pOClEdHqq KAPTA dOlvHA COdEPvATx 4 zBEzdO~Ki
! **********************************************************************
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! COMMON variables
!
          character(4) BAR(10), TITRE(30)
          real*8 DCOP, DEN(30), DEXE, DEXV, DFL(30), DP(251), &
              DPAS, DQ(251), DQ1(30), DR(251), DV(251), Z
          integer ICUT, IEX, ION, IPRAT, IWAT, MAG(30), NCORB(30), &
              NEL(30), NES, NITER, NK(30), NMAX(30), NORB, NP, NQL(30), &
              NQN(30), NSTOP, NUC, NWF(30)
          real(8) RWAT, TEST, TESTE, TESTV, TESTY, TETS
          common DEN, DQ1, DFL, NQN, NQL, NK, NMAX, NEL, NORB, MAG
          common / DIRA / DV, DR, DP, DQ, DPAS, Z, NSTOP, NES, TETS, NP, &
              NUC
          common / HLS / RWAT, IWAT, IEX
          common / PS2 / DEXV, DEXE, DCOP, TITRE, BAR, TEST, TESTE, TESTY, &
              TESTV, NITER, ION, ICUT, IPRAT
          common / UKP / NCORB, NWF
!
! Dummy arguments
!
          integer IZ, IZCORE, NPCH
!
! Local variables
!
          real*8 D1, DR1, DVAL, DVC
          character(4) ENDA, ITXCH(3, 2), TTIRE(9)
          real(8) FPOT
          integer I, ICH, IDEP, IMCH, ION1, IZ1, J, K, L, NKAI, NUC1
          real(8) R, VAL
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
          data TTIRE/'S   ', 'P*  ', 'P   ', 'D*  ', 'D   ', 'F*  ', 'F   &
         &', 'G*  ', 'G   '/, ENDA/'****'/
          data ITXCH/'BART', 'H-HE', 'DIN ', 'X-AL', 'PHA ', '    '/
!      IF (NSTOP.EQ.0) GO TO 2
          if (NSTOP == 1) then
              DPAS = 0.05D0
              DR1 = 1D-2
              !0.01D0
              NES = 15
              NITER = 50
              TESTE = 5.D-06
              TESTY = 1.D-05
              TESTV = 1.D-05
              TEST = 1.D-07
              NP = 251
              NSTOP = 30
              DEXV = 1.D0
              DEXE = 1.5D0
              DCOP = 0.3D0
              !137.0373D0
              DVC = 137.0359895D0
              IZ1 = 0
              ION1 = 0
              NUC1 = -1
              ION = 0
              IDEP = 0
              ICUT = 0
              IPRAT = 1
              IEX = 0
              I = 0
              J = 0
              L = 0
              K = 0
              NUC = 0
              NPCH = 0
              IWAT = 0
              DVAL = 0.D0
              do I = 1, 7
                  BAR(I) = '----'
              end do
!      READ (5,1000) (BAR(I),I=1,7)
!
! BAR TITRE COdEPvAT        24 CiMBOlA
! **********************************************************************
              if (BAR(1) == ENDA) call EXIT_A
              I = 0
              call DEFORB(NQN, NK, NEL, NCORB, NORB, IZ, IZCORE, MAG)
!
! IZ ATOMHyj HOMEP      ION=IZ-~iClO |lEKTPOHOB
! NORB ~iClO OPbiTAlEj    IDEP dOlvH byTx PABEH 0 ili 1
! IDEP=0 CTAPTOByj pOTHciAl =pOTEHciAl TOMACA-fEPMi
! IDEP=1 CTAPTOByj pOTEHciAl ~iTAETCq C KAPT
! ECli ICUT HOlx UL ON CORRIGE LE POTENTIEL EN -(ION+1)/R
! ECli IPRAT HOlx iCpOlxzyETCq pPOcEduPA pPATTA
! ECli IEX HOlx - ObMEH  BARTH-HEDIN'A
! I ~iClO TO~EK dlq iHTEgPiPOBAHiq (251 ECli I=0
! J ~iClO pOpyTOK pOdOgHATx |HEPgi`(15 ECli J=0
! K ~iClO iTEPAcij (50 ECli K=0
! L=0 CTAHdAPTHO BybPAHHAq TO~HOCTx
! KOHE~HyE PAzMEPy qdPA ECli NUC pOlOviTElEH
! **********************************************************************
              do ICH = 8, 10
                  IMCH = ICH - 7
                  BAR(ICH) = ITXCH(IMCH, IEX + 1)
              end do
              if (NORB <= NSTOP) then
                  if (I > 0) then
                      I = 2 * (I / 2) + 1
                      if (I <= NP) then
                          NP = I
                      else
                          write (8, 99011) I
                          goto 100
                      end if
                  end if
                  if (J > 0) NES = J
                  if (K > 0) NITER = K
                  if (IEX == 0) then
                  end if
!      READ (5,1005) DEXV,DEXE
!
! DEXV KO|fficiEHT pPi ObMEHHOM pOTEHciAlE      DEXV=1. dlq  SLATER
! DEXE KO|fficiEHT pPi ObMEHHOj |HEPgii
! DEXV dOlvEH byTx PABEH 2.*DEXE/3. dlq TEOPEMy BiPiAlA
! **********************************************************************
                  if (L == 0) then
                  end if
!      READ (5,1002) DPAS,DR1,TEST,TESTE,TESTY,TESTV
!
! DPAS |KCp.{Ag    DR1 OpPEdElqET pEPBu` TO~Ku  =DR1/IZ
! TEST TO~HOCTx pO |HEPgii dlq RESLD
! TESTE KPiTEPij CAMOCOglACOBAHiq OdHO|lEKTPOHHyX |HEPgij
! TESTE KPiTEPij CAMOCOglACOBAHiq BOlHOByX fuHKcij
! TESTE KPiTEPij CAMOCOglACOBAHiq pOTEHciAlA
! **********************************************************************
                  if (IPRAT /= 0) DCOP = 0.3D0
!
! VI(N+1)=(1.-DCOP)*VI(N)+DCOP*VF(N)
! **********************************************************************
                  Z = DFLOAT(IZ)
                  if (NUC > 0) then
!
! DVAL ATOMHAq MACCA ECli NUC pOlOviTElEH
! **********************************************************************
                      DVAL = Z * (DVAL**(1.D0 / 3.D0)) * 2.2677D-05 / exp &
                             (4.D0 * DPAS)
                      if (DVAL <= DR1) then
                          DR1 = DVAL
                          NUC = 5
                          goto 20
                      else
                          DVAL = DVAL * exp(4.D0 * DPAS)
                          do I = 6, NP
                              D1 = DR1 * exp(DFLOAT(I - 1) * DPAS)
                              if (D1 >= DVAL) goto 10
                          end do
                          write (8, 99018)
                          goto 100
                      end if
10                    continue
                      NUC = I
                      DR1 = DR1 * DVAL / D1
                  end if
20                continue
                  write (8, 99002) IZ, ION, NITER, TESTE, TESTY, TESTV
                  write (8, 99009) NP, DR1, IZ, DPAS
                  write (8, 99010) TEST, NES
                  if (IEX == 0) write (8, 99001)
                  if (IEX /= 0) write (8, 99013) DEXV, DEXE
                  K = 0
                  DVAL = Z * Z / (DVC * DVC)
                  if (NUC > 0) write (8, 99017)
!
! dAHHyE pO OPbiTAlqM
! **********************************************************************
                  do I = 1, NORB
!
! DEN |HEPgiq OPbiTAli B Ed.lAHdAu (<0)
! NQN glABHOE KBAHTOBOE ~iClO      NK KBAHTOBOE ~iClO KAppA
! NEL zAHqTOCTx OPbiTAli NCORB pPizHAK BAlEHTHOCTi(=1 dlq BAl.
! **********************************************************************
                      K = K + NEL(I)
                      DEN(I) = -Z * Z / (4.D0 * NQN(I) * NQN(I))
                      NQL(I) = IABS(NK(I))
                      if (NK(I) < 0) NQL(I) = NQL(I) - 1
                      if (NUC > 0) then
                          NKAI = IABS(NK(I))
                          DFL(I) = DFLOAT(NKAI)
                      else
                          DFL(I) = NK(I) * NK(I)
                          DFL(I) = sqrt(DFL(I) - DVAL)
                      end if
                      L = 2 * IABS(NK(I))
                      if (NQL(I) < NQN(I) .and. NEL(I) <= L .and. NQN(I) &
                          > 0 .and. NQL(I) <= 4) then
                          J = NQL(I) + IABS(NK(I))
                          TITRE(I) = TTIRE(J)
                          write (8, 99006) NQN(I), TITRE(I), NCORB(I), &
                              NEL(I), DEN(I)
                      else
                          write (8, 99004) DEN(I), NQN(I), NQL(I), J, NEL(I)
                          goto 100
                      end if
                  end do
                  if (K == (IZ - ION)) then
                      if (IPRAT /= 0) then
                          write (8, 99016) DCOP
                      else
                          write (8, 99015)
                      end if
                      if (NUC == NUC1) then
                          if (IZ == IZ1 .and. ION == ION1) goto 50
                          if (IZ == IZ1) goto 40
                      end if
                      DR(1) = DR1 / Z
                      do I = 2, NP
                          DR(I) = DR(1) * exp(DFLOAT(I - 1) * DPAS)
                      end do
                  else
                      write (8, 99005)
                      goto 100
                  end if
!
! CTAPTOByj pOTEHciAl
! **********************************************************************
40                continue
                  VAL = -ION - 1
                  if (IDEP == 1) then
!      READ(5,1004) (DV(I),I=1,NP)
!
! ~TEHiE CTAPTOBOgO pOTEHciAlA   (EN U.A. ET NEGATIF) ECli IDEP=1
! **********************************************************************
                      write (8, 99012) BAR, (DV(I), I=1, NP)
                      DVAL = -Z / DV(1)
                      if (NUC > 0) DVAL = 1.D0
                      do I = 1, NP
                          DV(I) = DV(I) * DVAL / DR(I)
                      end do
                  elseif (IDEP == 0) then
                      if (IZ /= IZ1 .or. ION <= ION1 .or. NUC /= NUC1) then
                          do I = 1, NP
                              R = DR(I)
                              DV(I) = FPOT(R, Z, VAL)
                          end do
                          if (IWAT /= 0) then
                              do I = 1, NP
                                  if (DR(I) > RWAT) then
                                      DV(I) = DV(I) + ION / DR(I)
                                  else
                                      DV(I) = DV(I) + ION / RWAT
                                  end if
                              end do
                          end if
                          if (NUC > 0) then
                              do I = 1, NUC
                                  DV(I) = DV(I) + Z / DR(I) + Z * ((DR(I) &
                                                                    / DR(NUC))**2 - 3.D0) / (DR(NUC) + DR(NUC))
                              end do
                          end if
                          goto 50
                      end if
                  else
                      write (8, 99014)
                      goto 100
                  end if
                  if (ICUT == 0) then
                      do I = 1, NP
                          if ((DR(I) * DV(I)) > VAL) DV(I) = VAL / DR(I)
                      end do
                  end if
                  VAL = Z + DV(1) * DR(1)
                  if (NUC > 0) VAL = Z + DV(NUC) * DR(NUC)
                  if (abs(VAL) >= 0.1D0) then
                      write (8, 99007)
                      goto 100
                  end if
              else
                  write (8, 99003) NORB
                  goto 100
              end if
50            continue
              if (NORB /= 1) then
                  do I = 2, NORB
                      K = I - 1
                      do J = 1, K
                          if (NQN(I) == NQN(J) .and. NK(I) == NK(J)) then
                              write (8, 99008)
                              goto 100
                          end if
                      end do
                  end do
              end if
              IZ1 = IZ
              ION1 = ION
              NUC1 = NUC
              do I = 1, NORB
                  NMAX(I) = NP
                  L = 1
                  J = NQN(I) - NQL(I)
                  if ((J - 2 * (J / 2)) == 0) L = -L
                  NKAI = IABS(NK(I))
                  DQ1(I) = DFLOAT(L * NK(I) / NKAI)
                  if (NUC /= 0) then
                      if (NK(I) < 0) DQ1(I) = DQ1(I) * (NK(I) - DFL(I)) &
                                              * DVC / Z
                  end if
              end do
              return
          end if
100       continue
          write (8, 99019)
!    READ (5,1000) BAR(1)
          if (BAR(1) == ENDA) then
          end if
          NSTOP = 1
99001     format(' Exchange: Barth-Hedin   ')
99002     format(' Atomic number  ', I3, '   Ionicity  ', I2 / 1X,              &
                         &        'Maximal number of iterations ',                          &
                         &        I4/' Precision in energy ', E21.10 / 14X, 'wave function ',   &
                         &        E21.10 / 14X, 'potential', E21.10/)
99003     format(' NORB=', I3, 'too large ************** ')
99004     format(' Input error       ', E21.10, I1, 2I2)
99005     format(' Erroneous number of electrons  **************')
99006     format(7X, I1, A2, 2I8, 1PE21.10)
99007     format(' Error in potential   ')
99008     format(' Bad configuration ')
99009     format(' Integration will be in ', i3,                            &
                         &        ' points'/' the first point is ', f7.4, '/', i2, '  step is ',&
                         &        f7.4,/)
!$$$ 2008 FORMAT ('  iHTEgPiPOBAHiE budET pPOBOdiTxCq B  ',I3,
!$$$     = ' TO4KAX ',/,'    pEPByj 4lEH PABEH    ',F7.4,'/'
!$$$     =,I2,' {Ag PABEH   ',F7.4,/)
99010     format(' In the program RESLD the relative precision in energy', &
                         &        ' is ', 1PE9.2, /' the number of attempts ', i3,/)
          !/,'
!    2IDEP=',I3,2X,'ICUT=',I3,2X,'IPRAT=',I3,2X,'ION=',I3,2X,'IWAT=',
!    1I3,2X,'RWAT=',E10.3,/)
!$$$ 2009 FORMAT ('       pPi  PE{EHii uPABHEHiq diPAKA ',
!$$$     */,'    OTHOCiTElxHAq TO4HOCTx pO |HEPgii ',E21.10,
!$$$     */,'                4iClO pOpyTOK         ',I3,/)
99011     format('  NP=', I3, ' is too large *******')
99012     format(20X, 10A4, //, 5X, 'Starting potential * R'/, 5(2X, F21.15))
99013     format('  Exchange: Slater X-alpha, DEXV=', f21.15, '     DEXE=',  &
                         &        f21.15,/)
99014     format(' Error in  IDEP *********************')
99015     format(' Pratt''s mixing is used  '/)
99016     format(' Potential is mixed with  ADMIX=', E21.10,/)
99017     format(10X, ' Finite nucleus'/)
99018     format(' Error in atomic weight  **************** ')
99019     format(' The next case ')
!      GO TO 1
      end subroutine INSLD
!*==somm.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine SOMM(DR, DP, DQ, DPAS, DA, M, NP)
!
! iHTEgPiPOBAHiE METOdOM CiMpCOHA  (DP+DQ)*DR**M OT 0 dO R=DR(NP)
! DPAS |KCpOHEHciAlxHyj {Ag dlq R->0   (DP+DQ)=CTE*R**DA
! **********************************************************************
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 DA, DPAS
          integer M, NP
          real*8 DP(*), DQ(*), DR(*)
!
! Local variables
!
          real*8 D1, DB, DC, DL
          integer I, MM
!
!*** End of declarations rewritten by SPAG
!
          MM = M + 1
          D1 = DA + DFLOAT(MM)
          DA = 0.D0
          DB = 0.D0
          do I = 1, NP
              DL = DR(I)**MM
              if (I /= 1 .and. I /= NP) then
                  DL = DL + DL
                  if ((I - 2 * (I / 2)) == 0) DL = DL + DL
              end if
              DC = DP(I) * DL
              if (DC < 0) then
                  DB = DB + DC
              elseif (DC /= 0) then
                  DA = DA + DC
              end if
              DC = DQ(I) * DL
              if (DC < 0) then
                  DB = DB + DC
              elseif (DC /= 0) then
                  DA = DA + DC
              end if
          end do
          DA = DPAS * (DA + DB) / 3.D0
          DC = exp(DPAS) - 1.D0
          DB = D1 * (D1 + 1.D0) * DC * exp((D1 - 1.D0) * DPAS)
          DB = DR(1) * (DR(2)**M) / DB
          DC = (DR(1)**MM) * (1.D0 + 1.D0 / (DC * (D1 + 1.D0))) / D1
          DA = DA + DC * (DP(1) + DQ(1)) - DB * (DP(2) + DQ(2))
      end subroutine SOMM
!*==fpot.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      function FPOT(R, Z, WA)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real(8) R, WA
          real*8 Z
          real(8) FPOT
!
! Local variables
!
          real(8) S, WC, WD, WE
!
!*** End of declarations rewritten by SPAG
!
!
! pOTEHciAl TOMACA-fEPMi B TO~KE  R  Z ATOMHyj HOMEP
! WA ~iClO |lEKTPOHOB-Z-1
! **********************************************************************
          S = Z
          WC = sqrt((R * (S + WA)**(1.D0 / 3.D0)) / 0.8853D0)
          WD = WC * (0.60112 * WC + 1.81061D0) + 1.D0
          WE = WC * (WC * (WC * (WC * (0.04793D0 * WC + 0.21465D0) + &
                                 0.77112D0) + 1.39515D0) + 1.81061D0) + 1.D0
          WC = (Z + WA) * (WD / WE)**2 - WA
          FPOT = -WC / R
      end function FPOT
!*==potsl.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine POTSL(DV, D, DP, DR, DPAS, DEXV, Z, NP, ION, ICUT)
!
! iHTEgPiPOBAHiE pOTEHciAlA pO 4 TO~KAM
! DV pOTEHciAl   D plOTHOCTx DP BLOC DE TRAVAIL  DR PAdiAlxHAq {KAlA
! DPAS |KCp.{Ag          DEXV MHOviTElx dlq ObMEHA
! Z ATOMHyj HOMEP     NP ~iClO TO~EK        ION=Z-~iClO |lEKTPOHOB
! SI ICUT EST NUL ON CORRIGE EVENTUELLEMENT LE POTENTIEL EN -(ION+1)/R
! **********************************************************************
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! COMMON variables
!
          integer IEX, IWAT
          real(8) RWAT
          common / HLS / RWAT, IWAT, IEX
!
! Dummy arguments
!
          real*8 DEXV, DPAS, Z
          integer ICUT, ION, NP
          real*8 D(*), DP(*), DR(*), DV(*)
!
! Local variables
!
          real*8 DARS, DAS, DBRS, DCF, DCP, DLO, DLO2, DNY, DR2, &
              DS, DS1, DSF, DSF2, DSF3, DSP, DSP2, DSP3, DVXC
          integer I, J, K
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
          DLO = -DFLOAT(ION + 1)
          if (IEX == 0) then
              do I = 1, NP
!
!     pOTEHciAl pO BARTH-HEDIN
!  *********************************************************************
                  DR2 = DR(I)**2
                  DS1 = (D(I) / (3.D0 * DR2) + 1.D-42)**(1.D0 / 3.D0)
!
                  if (abs(DS1) < 1.D-10) then
                      DVXC = 0.D0
                  else
                      DS = 1.D0 / DS1
                      DSF = DS / 75.D0
                      DSF2 = DSF * DSF
                      DSF3 = DSF2 * DSF
                      DSP = DS / 30.D0
                      DSP2 = DSP * DSP
                      DSP3 = DSP2 * DSP
                      DCF = (1.D0 + DSF3) * log(1.D0 + 1.D0 / DSF) + 0.5D0 * &
                            DSF - DSF2 - 0.3333333333D0
                      DCP = (1.D0 + DSP3) * log(1.D0 + 1.D0 / DSP) + 0.5D0 * &
                            DSP - DSP2 - 0.3333333333D0
                      DNY = 5.1297628D0 * (0.0504D0 * DCP - 0.0254D0 * DCF)
                      DARS = -1.22177412D0 / DS + DNY
                      DBRS = -0.0504D0 * log(1.D0 + 30.D0 / DS) - DNY
                      DVXC = DARS + DBRS
                  end if
                  DV(I) = DV(I) - (Z - 0.5D0 * DR(I) * DVXC)
!  *********************************************************************
                  if (ICUT == 0) then
                      if (DV(I) > DLO) DV(I) = DLO
                  end if
                  DV(I) = DV(I) / DR(I)
              end do
              if (IWAT == 0) return
              do I = 1, NP
                  if (DR(I) > RWAT) then
                      DV(I) = DV(I) + ION / DR(I)
                  else
                      DV(I) = DV(I) + ION / RWAT
                  end if
              end do
              goto 99999
          end if
!
!     pOTEHciAl pO Cl|TTEPu
!
          do I = 1, NP
              DV(I) = DV(I) - (Z + 3.D0 * DEXV * ((DR(I) * D(I) / &
                                                   105.27578D0 + 1.D-42)**(1.D0 / 3.D0)))
              if (ICUT == 0) then
                  if (DV(I) > DLO) DV(I) = DLO
              end if
              DV(I) = DV(I) / DR(I)
          end do
          if (IWAT == 0) return
          do I = 1, NP
              if (DR(I) > RWAT) then
                  DV(I) = DV(I) + ION / DR(I)
              else
                  DV(I) = DV(I) + ION / RWAT
              end if
          end do
          return
99999     continue
      end subroutine POTSL
!*==deforb.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine DEFORB(XN, NK, XZ, IVAL, NORB, IZN, IZC, MAG)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer IZC, IZN, NORB
          integer IVAL(30), MAG(30), NK(30), XN(30), XZ(30)
!
! Local variables
!
          integer I, IXT, IXV, IXZ1(105), IXZ2(105), IXZ3(105), IXZ4( &
              105), IXZ5(105), IZ, IZV, J, JF1, JFUL, JFULR(105), JREL(105), &
              NMAG
          real(8) XJ(35), XL(35), XMAG, XML, XXJ(31), XXL(31), XXN(31), &
              XXZ(31), ZN
!
!*** End of declarations rewritten by SPAG
!
!================================================================
!
!     STORED DATA FOR FULL ORBITALS (AS THEY ARE NUMBERED IN 'ATOMGL'):
!
!     STORED DATA FOR RELATIVISTIC ORBITALS:
!
          data XXN/1.0D0, 2.0D0, 2.0D0, 2.0D0, 3.0D0, 3.0D0, 3.0D0, 4.0D0, &
              3.0D0, 3.0D0, 4.0D0, 4.0D0, 5.0D0, 4.0D0, 4.0D0, 5.0D0, 5.0D0, &
              6.0D0, 4.0D0, 4.0D0, 5.0D0, 5.0D0, 6.0D0, 6.0D0, 7.0D0, 5.0D0, &
              5.0D0, 6.0D0, 6.0D0, 7.0D0, 7.0D0/
          data XXL/0.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, &
              2.0D0, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 2.0D0, 2.0D0, 1.0D0, 1.0D0, &
              0.0D0, 3.0D0, 3.0D0, 2.0D0, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 3.0D0, &
              3.0D0, 2.0D0, 2.0D0, 1.0D0, 1.0D0/
          data XXJ/0.5D0, 0.5D0, 0.5D0, 1.5D0, 0.5D0, 0.5D0, 1.5D0, 0.5D0, &
              1.5D0, 2.5D0, 0.5D0, 1.5D0, 0.5D0, 1.5D0, 2.5D0, 0.5D0, 1.5D0, &
              0.5D0, 2.5D0, 3.5D0, 1.5D0, 2.5D0, 0.5D0, 1.5D0, 0.5D0, 2.5D0, &
              3.5D0, 1.5D0, 2.5D0, 0.5D0, 1.5D0/
          data XXZ/2.0D0, 2.0D0, 2.0D0, 4.0D0, 2.0D0, 2.0D0, 4.0D0, 2.0D0, &
              4.0D0, 6.0D0, 2.0D0, 4.0D0, 2.0D0, 4.0D0, 6.0D0, 2.0D0, 4.0D0, &
              2.0D0, 6.0D0, 8.0D0, 4.0D0, 6.0D0, 2.0D0, 4.0D0, 2.0D0, 6.0D0, &
              8.0D0, 4.0D0, 6.0D0, 2.0D0, 4.0D0/
          data JFULR/0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 6, 7, &
              7, 8, 8, 8, 8, 7, 9, 9, 9, 9, 7, 10, 10, 11, 11, 11, 11, 12, 12, &
              13, 13, 13, 12, 12, 14, 12, 12, 12, 12, 15, 15, 16, 16, 16, 16, &
              17, 17, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, &
              20, 20, 20, 20, 21, 21, 21, 21, 17, 17, 22, 22, 23, 23, 23, 23, &
              24, 24, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, &
              27, 27, 27, 27/
          data JREL/1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7, &
              8, 8, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 12, 12, 12, 12, &
              13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 16, 16, 17, 17, &
              17, 17, 18, 18, 21, 21, 19, 19, 19, 19, 20, 21, 20, 20, 20, 20, &
              20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 23, 23, 24, 24, &
              24, 24, 25, 25, 28, 28, 28, 28, 28, 26, 27, 28, 28, 27, 27, 27, &
              27, 27, 28, 28, 28/
          data IXZ1/1, 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 0, 1, 0, 1, 2, 3, 0, &
              1, 0, 1, 2, 3, 1, 1, 2, 3, 4, 1, 0, 1, 0, 1, 2, 3, 0, 1, 0, 1, 2, &
              1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 0, 0, 1, 3, 4, 5, 0, &
              1, 1, 3, 4, 5, 6, 7, 0, 1, 2, 3, 0, 1, 2, 3, 1, 1, 0, 1, 0, 1, 2, &
              3, 0, 1, 0, 0, 0, 2, 3, 4, 0, 1, 1, 2, 4, 5, 6, 7, 0, 1, 2, 3/
          data IXZ2/23*0, 4, 4*0, 4, 11*0, 4, 4, 0, 4, 4, 4, 4, 16* &
              0, 1, 13*0, 6, 6, 16*0, 1, 1, 8*0/
          data IXZ3/23*0, 1, 4*0, 6, 12*0, 1, 0, 3, 4, 6, 6, 9*0, &
              1, 1, 19*0, 8, 8, 9*0, 1, 2, 1, 1, 1, 12*0/
          data IXZ4/77*0, 4, 4, 26*0/
          data IXZ5/77*0, 5, 6, 26*0/
!================================================================
          ZN = IZN
          IZ = int(ZN + 0.000001D0)
          IZV = IZN - IZC
!------ J-REPRESENTATION ------
          JFUL = JFULR(IZ)
!
          J = JREL(IZ)
          do I = 1, JFUL
              XN(I) = nint(XXN(I))
              XL(I) = XXL(I)
              XJ(I) = XXJ(I)
              XZ(I) = nint(XXZ(I))
          end do
!
          XZ(JFUL + 1) = IXZ1(IZ)
          XZ(JFUL + 2) = IXZ2(IZ)
          XZ(JFUL + 3) = IXZ3(IZ)
          XZ(JFUL + 4) = IXZ4(IZ)
          XZ(JFUL + 5) = IXZ5(IZ)
!
!--------------------------------
          if (JFUL /= J) then
              JF1 = JFUL + 1
              do I = JF1, J
                  XN(I) = nint(XXN(I))
                  XL(I) = XXL(I)
                  XJ(I) = XXJ(I)
!--------------------------------
              end do
          end if
!===============================================================
          IXV = 0
          IXT = 0
          do I = J, 1, -1
              IXT = IXT + XZ(I)
              NK(I) = nint(XL(I))
              MAG(I) = 0
              XMAG = 0
              if (XJ(I) > XL(I)) NK(I) = nint(-XL(I) - 1)
              if (IXV + XZ(I) <= IZV) then
                  if (I >= 31) then
                      write (16, *) '########### POTSL: I = ', I, ' > 30'
                      stop
                  end if
                  IVAL(I) = 1
                  IXV = IXV + XZ(I)
                  if (NK(I) > 1) then
                      XML = XZ(I) + XZ(I + 1)
                      NMAG = NK(I) * 2 + 1
                      XMAG = NMAG - abs(XML - NMAG)
                      MAG(I) = nint(XMAG)
!            XMAG1=(nmag-1)*xmag/(2.*nmag)
!            XMAG2=xmag-xmag1
!            print*,xz(i),xz(i+1)
!            print*, i,nk(i),' l=',nint(xl(i)),' n='
!     $           ,xn(i),' val=',nint(XMAG),xmag1,xmag2
                  end if
              else
                  IVAL(I) = 0
              end if
          end do
          if (IXV /= IZV) stop ' error - I can not make initial data '
          NORB = J
      end subroutine DEFORB
!*==exit_a.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine EXIT_A
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      end subroutine EXIT_A
!*==resld.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine RESLD(NQN, NQL, NK, IMAX, DE, DFL, DQ1, JC1)
!
! PE{EHiE uPABHEHiq diPAKA
! NQN NglABHOE KBAHTOBOE ~iClO   NQL OPbiTAlxHOE KBAHTOBOE ~iClO
! NK KBAHTOBOE ~iClO KAppA    IMAX pOClEdHqq TO~KA B TAbulqcii BOl-
! HOBOj fuHKcii     DE |HEPgiq   DFL pOKAzATElx CTEpEHi B PAzlOvEHii
! BOlHOBOj fuHKcii       DQ1 OTHO{EHiE DP K DQ B HA~AlE KOOPdiHAT
! **********************************************************************
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! COMMON variables
!
          real*8 DB, DEP(5), DEQ(5), DK, DM, DP(251), DPAS, &
              DPNO(4, 30), DQ(251), DQNO(4, 30), DR(251), DSAL, DV(251), &
              DVC, Z
          integer NES, NP, NSTOP, NUC
          real(8) TEST
          common / DIRA / DV, DR, DP, DQ, DPAS, Z, NSTOP, NES, TEST, NP, &
              NUC
          common / PS1 / DEP, DEQ, DB, DVC, DSAL, DK, DM
          common / TROIS / DPNO, DQNO
!
! Dummy arguments
!
          real*8 DE, DFL, DQ1
          integer IMAX, JC1, NK, NQL, NQN
!
! Local variables
!
          real*8 DBE, DD, DKOEF, DPM, DPQ, DQM, DSUM, DVAL
          real(8) ELIM, VAL
          integer I, IES, IMAT, IMM, J, JC, K, LLL, M, ND, NOEUD
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
! DV pOTEHciAl B AT.Ed.(<0)        DR PAdiAlxHAq CETKA
! DP bOlx{Aq KOMpOHEHTA  DQ MAlAq KOMpOHEHTA    DPAS |KCpOHEHc.{Ag
! Z ATOMHyj HOMEP     NSTOP KOHTPOlx ~iClEHHOgO iHTEgPiPOBAHiq
! NES MAKCiMAlxHOE ~iClO iTEPAcij pO OpPEdElEHi` |HEPgii
! TEST TO~HOCTx pOlu~EHiq |HEPgii          NP NMAKCiMAlxHOE ~iClO TO~EK
! KOHE~HyE PAzMEPy qdPA ECli    NUC HE   0
! **********************************************************************
!
! DEP,DEQ pPOizBOdHyE DP i  DQ   DB=ENERGIE/DVC    DVC CKOPOCTx CBETA
!  B AT.Ed.         DSAL=2.*DVC   DK KBAHTOBOE ~iClO KAppA
! DM=|KCpOHEHciAlxHyj {Ag/720., DKOEF=1./720.
! **********************************************************************
          data DKOEF/.1388888888888888D-2/
          JC = 1
          if (JC1 > 0) JC = JC1
          NSTOP = 0
          DVC = 137.0373D0
          DSAL = DVC + DVC
          IMM = 0
          IES = 0
          DK = NK
          LLL = (NQL * (NQL + 1)) / 2
          ND = 0
          NOEUD = NQN - NQL
          if (LLL /= 0) then
              ELIM = DV(1) + LLL / (DR(1) * DR(1))
              do I = 2, NP
                  VAL = DV(I) + LLL / (DR(I) * DR(I))
                  if (VAL <= ELIM) ELIM = VAL
              end do
              if (ELIM < 0) then
                  if (DE <= ELIM) DE = ELIM * 0.5D0
              else
                  NSTOP = 17
! 2*V+L*(L+1)/R**2   BC`du pOlOviTElEH
! **********************************************************************
                  return
              end if
          else
              ELIM = -Z * Z / (1.5 * NQN * NQN)
              if (DE <= ELIM) DE = ELIM * 0.5D0
          end if
100       continue
          if (IMM /= 1) then
              do I = 7, NP, 2
                  IMAT = NP + 1 - I
                  if ((DV(IMAT) + DFLOAT(LLL) / (DR(IMAT) * DR(IMAT)) &
                       - DE) <= 0.D0) exit
              end do
              if (IMAT <= 5) then
                  DE = DE * 0.5D0
                  if (DE < -TEST .and. ND <= NOEUD) goto 100
                  NSTOP = 28
! 2*V+L*(L+1)/R**2-2*E BC`du pOlOviTElEH
! **********************************************************************
                  return
              end if
          end if
! HA~AlxHyE zHA~EHiq dlq iHTEgPiPOBAHiq pO BHuTPEHHEj OblACTi
! **********************************************************************
          DB = DE / DVC
          call INOUH(DP, DQ, DR, DQ1, DFL, DV(1), Z, TEST, NUC, NSTOP, JC)
          if (NSTOP /= 0) goto 99999
!     NSTOP=45
! HET CXOdiMOCTi B HA~AlE KOOPdiHAT
! **********************************************************************
          ND = 1
          do I = 1, 5
              DVAL = DR(I)**DFL
              if (I /= 1) then
                  if (DP(I - 1) /= 0.D0) then
                      if ((DP(I) / DP(I - 1)) <= 0.D0) ND = ND + 1
                  end if
              end if
              DP(I) = DP(I) * DVAL
              DQ(I) = DQ(I) * DVAL
              DEP(I) = DEP(I) * DVAL
              DEQ(I) = DEQ(I) * DVAL
          end do
          K = -1 + 2 * (NOEUD - 2 * (NOEUD / 2))
          if ((DP(1) * DFLOAT(K)) > 0.D0) then
              if ((DFLOAT(K) * DFLOAT(NK) * DQ(1)) >= 0.D0) then
                  DM = DPAS * DKOEF
! IiHTEgPiPOBAHiE pO BHuTPEHHEj OblACTi
! **********************************************************************
                  do I = 6, IMAT
                      DP(I) = DP(I - 1)
                      DQ(I) = DQ(I - 1)
                      call INTH(DP(I), DQ(I), DV(I), DR(I))
                      if (DP(I - 1) /= 0.D0) then
                          if ((DP(I) / DP(I - 1)) <= 0.D0) then
                              ND = ND + 1
                              if (ND > NOEUD) goto 200
                          end if
                      end if
                  end do
                  if (ND == NOEUD) then
! HA~AlHyE zHA~EHiq dlq iHTgPiPOBAHiq pO BHE{HEj OblACTi
! **********************************************************************
                      DQM = DQ(IMAT)
                      DPM = DP(IMAT)
                      if (IMM /= 1) then
                          do I = 1, NP, 2
                              IMAX = NP + 1 - I
                              if (((DV(IMAX) - DE) * DR(IMAX) * DR(IMAX)) &
                                  <= 300.D0) exit
                          end do
                      end if
                      DD = sqrt(-DE * (2.D0 + DB / DVC))
                      DPQ = -DD / (DSAL + DB)
                      DM = -DM
                      do I = 1, 5
                          J = IMAX + 1 - I
                          DP(J) = exp(-DD * DR(J))
                          DEP(I) = -DD * DP(J) * DR(J)
                          DQ(J) = DPQ * DP(J)
                          DEQ(I) = DPQ * DEP(I)
                      end do
                      M = IMAX - 5
! iHTEgPiPOBAHiE pO BHE{HEj OblACTi
!***********************************************************************
                      do I = IMAT, M
                          J = M + IMAT - I
                          DP(J) = DP(J + 1)
                          DQ(J) = DQ(J + 1)
                          call INTH(DP(J), DQ(J), DV(J), DR(J))
                      end do
! C{iBKA bOlx{Oj KOMpOHEHTy
! **********************************************************************
                      DVAL = DPM / DP(IMAT)
                      if (DVAL > 0.D0) then
                          do I = IMAT, IMAX
                              DP(I) = DP(I) * DVAL
                              DQ(I) = DQ(I) * DVAL
                          end do
! By~iClEHiE HOPMy
! **********************************************************************
                          DSUM = 3.D0 * DR(1) * (DP(1)**2 + DQ(1)**2) &
                                 / (DPAS * (DFL + DFL + 1.D0))
                          do I = 3, IMAX, 2
                              DSUM = DSUM + DR(I) * (DP(I)**2 + DQ(I)**2) &
                                     + 4.D0 * DR(I - 1) * (DP(I - 1)**2 + DQ(I - 1) &
                                                           **2) + DR(I - 2) * (DP(I - 2)**2 + DQ(I - 2) &
                                                                               **2)
                          end do
                          DSUM = DPAS * (DSUM + DR(IMAT) * (DQM * DQM - DQ &
                                                            (IMAT) * DQ(IMAT))) * 0.3333333333333333D0
! izMEHEHiE |HEPgii
! **********************************************************************
                          DBE = DP(IMAT) * (DQM - DQ(IMAT)) * DVC / DSUM
                          IMM = 0
                          VAL = abs(DBE / DE)
                          if (VAL <= TEST) then
                              if (JC1 <= 0 .and. IMAX > (-JC1)) then
                                  IMAX = -JC1
                                  DSUM = 3.D0 * DR(1) * (DP(1)**2 + DQ(1)**2) &
                                         / (DPAS * (DFL + DFL + 1.D0))
                                  do I = 3, IMAX, 2
                                      DSUM = DSUM + DR(I) * (DP(I)**2 + DQ(I)**2) &
                                             + 4.D0 * DR(I - 1) * (DP(I - 1)**2 + DQ(I - 1) &
                                                                   **2) + DR(I - 2) * (DP(I - 2)**2 + DQ(I - 2) &
                                                                                       **2)
                                  end do
                                  DSUM = DPAS * (DSUM + DR(IMAT) * (DQM * DQM - DQ &
                                                                    (IMAT) * DQ(IMAT))) * 0.3333333333333333D0
                              end if
                              DSUM = sqrt(DSUM)
                              DQ1 = DQ1 / DSUM
                              do I = 1, IMAX
                                  DP(I) = DP(I) / DSUM
                                  DQ(I) = DQ(I) / DSUM
                              end do
                              do I = 1, 4
                                  DPNO(I, JC) = DPNO(I, JC) / DSUM
                                  DQNO(I, JC) = DQNO(I, JC) / DSUM
                              end do
                              if (IMAX /= NP) then
                                  J = IMAX + 1
                                  do I = J, NP
                                      DP(I) = 0.D0
                                      DQ(I) = 0.D0
                                  end do
                              end if
                              NSTOP = 0
                              goto 99999
                          else
102                           continue
                              DVAL = DE + DBE
                              if (DVAL < 0.D0) then
                                  DE = DVAL
                                  if (VAL <= 0.1D0) IMM = 1
                                  IES = IES + 1
                                  if (IES <= NES) goto 100
                                  NSTOP = 362
! pPEBy{EH ~iClO iTEPAcij
! **********************************************************************
                                  return
                              else
                                  DBE = DBE * 0.5D0
                                  VAL = VAL * 0.5D0
                                  if (VAL > TEST) goto 102
                                  NSTOP = 345
! HulEBAq |HEPgiq
! **********************************************************************
                                  return
                              end if
                          end if
                      else
                          NSTOP = 312
! O{ibKA B zHAKE bOlx{Oj KOMpOHEHTy
! **********************************************************************
                          return
                      end if
                  else
                      DE = 0.8D0 * DE
                      if (DE < -TEST) goto 100
                      NSTOP = 206
! ~iClO uzlOB Cli{KOM MAlO
! **********************************************************************
                      return
                  end if
              end if
          end if
          NSTOP = 53
! O{ibKA B PAzlOvEHii B HEA~AlE KOOPdiHAT
! **********************************************************************
          return
200       continue
          DE = 1.2D0 * DE
          if (DE > ELIM) goto 100
          NSTOP = 210
! ~iClO uzlOB Cli{KOM BEliKO
! **********************************************************************
          return
99999     continue
      end subroutine RESLD
!*==inouh.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine INOUH(DP, DQ, DR, DQ1, DFL, DV, Z, TEST, NUC, NSTOP, &
                       JC)
!
! HA~AlxHyE zHA~EHiq dlq iHTEgPiPOBAHiq pO BHuTPEHHEj OblACTi
! DP bOlx{Aq KOMpOHEHTA   DQ MAlAq KOMpOHEHTA     DR PAdiAlxHAq.CETKA
! DQ1 OTHO{EHiE DP K DQ B HA~AlE KOOPd.DFL pOKAzATElx CTEpEHi glABHOgO
! ~lEHA PAzlOvEHiq B HA~AlE DV pOTEHciAl B pEPBOj TO~KE
! Z ATOMHyj HOMEP      TEST TO~HOCTx
! KOHE~HyE PAzMEPy qdPA ECli    NUC HE 0
! NSTOP KOHTPOlx CXOdiMOCTi pPi PAzlOvEHii B Pqd
! **********************************************************************
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! COMMON variables
!
          real*8 DD, DEP(5), DEQ(5), DK, DM, DPNO(4, 30), &
              DQNO(4, 30), DSAL, DVC
          common / PS1 / DEP, DEQ, DD, DVC, DSAL, DK, DM
          common / TROIS / DPNO, DQNO
!
! Dummy arguments
!
          real*8 DFL, DQ1, DV, Z
          integer JC, NSTOP, NUC
          real(8) TEST
          real*8 DP(*), DQ(*), DR(*)
!
! Local variables
!
          real*8 DBE, DEVA1, DEVA2, DEVA3, DPR, DQR, DSUM, DVAL
          integer I, J, M
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
! DEP,DEQ pPOizBOdHyE DP i  DQ   DD=|HEPgiq/DVC    DVC CKOPOCTx CBETA
!  B AT. Ed.        DSAL=2.*DVC   DK KBAHTOBOE ~iClO KAppA
! DM=|KCp.{Ag/720.
! **********************************************************************
          do I = 1, 10
              DP(I) = 0.D0
              DQ(I) = 0.D0
          end do
          if (NUC <= 0) then
              DVAL = Z / DVC
              DEVA1 = -DVAL
              DEVA2 = DV / DVC + DVAL / DR(1) - DD
              DEVA3 = 0.D0
              if (DK <= 0) then
                  DBE = (DK - DFL) / DVAL
              else
                  DBE = DVAL / (DK + DFL)
              end if
              DQ(10) = DQ1
              DP(10) = DBE * DQ1
          else
              DVAL = DV + Z * (3.D0 - DR(1) * DR(1) / (DR(NUC) * DR(NUC) &
                                                       )) / (DR(NUC) + DR(NUC))
              DEVA1 = 0.D0
              DEVA2 = (DVAL - 3.D0 * Z / (DR(NUC) + DR(NUC))) / DVC - &
                      DD
              DEVA3 = Z / (DR(NUC) * DR(NUC) * DR(NUC) * DSAL)
              if (DK <= 0) then
                  DP(10) = DQ1
              else
                  DQ(10) = DQ1
              end if
          end if
          do I = 1, 5
              DP(I) = DP(10)
              DQ(I) = DQ(10)
              DEP(I) = DP(I) * DFL
              DEQ(I) = DQ(I) * DFL
          end do
          M = 1
100       continue
          DM = M + DFL
          DSUM = DM * DM - DK * DK + DEVA1 * DEVA1
          DQR = (DSAL - DEVA2) * DQ(M + 9) - DEVA3 * DQ(M + 7)
          DPR = DEVA2 * DP(M + 9) + DEVA3 * DP(M + 7)
          DVAL = ((DM - DK) * DQR - DEVA1 * DPR) / DSUM
          DSUM = ((DM + DK) * DPR + DEVA1 * DQR) / DSUM
          J = -1
          do I = 1, 5
              DPR = DR(I)**M
              DQR = DSUM * DPR
              DPR = DVAL * DPR
              if (M /= 1) then
                  if (abs(DPR / DP(I)) <= TEST .and. abs(DQR / DQ(I)) &
                      <= TEST) J = 1
              end if
              DP(I) = DP(I) + DPR
              DQ(I) = DQ(I) + DQR
              DEP(I) = DEP(I) + DPR * DM
              DEQ(I) = DEQ(I) + DQR * DM
          end do
          if (J /= 1) then
              DP(M + 10) = DVAL
              DQ(M + 10) = DSUM
              M = M + 1
              if (M <= 20) goto 100
              NSTOP = 45
          end if
          do I = 1, 4
              DPNO(I, JC) = DP(I + 9)
              DQNO(I, JC) = DQ(I + 9)
          end do
      end subroutine INOUH
!*==inth.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine INTH(DP, DQ, DV, DR)
!
! iHTEgPiPOBAHiE pO METOdu AdAMCA pO 5-TO~KAM bOlx{Oj KOMpOHEHTy  DP i
! MAlOj KOMpOHEHTy DQ   B TO~KE DR  ; dv-pOTEHciAl B |TOj TO~KE
! **********************************************************************
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! COMMON variables
!
          real*8 DB, DEP(5), DEQ(5), DK, DM, DSAL, DVC
          common / PS1 / DEP, DEQ, DB, DVC, DSAL, DK, DM
!
! Dummy arguments
!
          real*8 DP, DQ, DR, DV
!
! Local variables
!
          real*8 DKOEF1, DKOEF2, DPR, DQR, DSUM
          integer I
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
! DEP,DEQ pPOizBOdHyE DP i  DQ   DD=|HEPgiq/DVC    DVC CKOPOCTx CBETA
!  B AT. Ed.        DSAL=2.*DVC   DK KBAHTOBOE ~iClO KAppA
! DM=|KCp.{Ag/720.
! DKOEF1=405./502., DKOEF2=27./502.
! **********************************************************************
          data DKOEF1/.9462151394422310D0/, DKOEF2/ &
              .5378486055776890D-1/
          DPR = DP + DM * ((251.D0 * DEP(1) + 2616.D0 * DEP(3) + 1901.D0 &
                            * DEP(5)) - (1274.D0 * DEP(2) + 2774.D0 * DEP(4)))
          DQR = DQ + DM * ((251.D0 * DEQ(1) + 2616.D0 * DEQ(3) + 1901.D0 &
                            * DEQ(5)) - (1274.D0 * DEQ(2) + 2774.D0 * DEQ(4)))
          do I = 2, 5
              DEP(I - 1) = DEP(I)
              DEQ(I - 1) = DEQ(I)
          end do
          DSUM = (DB - DV / DVC) * DR
          DEP(5) = -DK * DPR + (DSAL * DR + DSUM) * DQR
          DEQ(5) = DK * DQR - DSUM * DPR
          DP = DP + DM * ((106.D0 * DEP(2) + 646.D0 * DEP(4) + 251.D0 * &
                           DEP(5)) - (19.D0 * DEP(1) + 264.D0 * DEP(3)))
          DQ = DQ + DM * ((106.D0 * DEQ(2) + 646.D0 * DEQ(4) + 251.D0 * &
                           DEQ(5)) - (19.D0 * DEQ(1) + 264.D0 * DEQ(3)))
          DP = DKOEF1 * DP + DKOEF2 * DPR
          DQ = DKOEF1 * DQ + DKOEF2 * DQR
          DEP(5) = -DK * DP + (DSAL * DR + DSUM) * DQ
          DEQ(5) = DK * DQ - DSUM * DP
      end subroutine INTH
!*==dlamch.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      real*8 function DLAMCH(CMACH)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          real*8 ONE, ZERO
          parameter(ONE=1.0D+0, ZERO=0.0D+0)
!
! Dummy arguments
!
          character CMACH
!
! Local variables
!
          real*8 BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, &
              RND, SFMIN, SMALL, T
          integer BETA, IMAX, IMIN, IT
          logical FIRST, LRND
          logical LSAME
          save BASE, EMAX, EMIN, EPS, PREC, RMAX, RMIN, RND, SFMIN, T
          external DLAMC2, LSAME
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Save statement ..
!     ..
!     .. Data statements ..
          data FIRST/.true./
!     ..
!     .. Executable Statements ..
!
          if (FIRST) then
              FIRST = .false.
              call DLAMC2(BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX)
              BASE = BETA
              T = IT
              if (LRND) then
                  RND = ONE
                  EPS = (BASE**(1 - IT)) / 2
              else
                  RND = ZERO
                  EPS = BASE**(1 - IT)
              end if
              PREC = EPS * BASE
              EMIN = IMIN
              EMAX = IMAX
              SFMIN = RMIN
              SMALL = ONE / RMAX
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
              if (SMALL >= SFMIN) SFMIN = SMALL * (ONE + EPS)
          end if
!
          if (LSAME(CMACH, 'E')) then
              RMACH = EPS
          elseif (LSAME(CMACH, 'S')) then
              RMACH = SFMIN
          elseif (LSAME(CMACH, 'B')) then
              RMACH = BASE
          elseif (LSAME(CMACH, 'P')) then
              RMACH = PREC
          elseif (LSAME(CMACH, 'N')) then
              RMACH = T
          elseif (LSAME(CMACH, 'R')) then
              RMACH = RND
          elseif (LSAME(CMACH, 'M')) then
              RMACH = EMIN
          elseif (LSAME(CMACH, 'U')) then
              RMACH = RMIN
          elseif (LSAME(CMACH, 'L')) then
              RMACH = EMAX
          elseif (LSAME(CMACH, 'O')) then
              RMACH = RMAX
          end if
!
          DLAMCH = RMACH
!
!     End of DLAMCH
!
      end function DLAMCH
!*==dlamc1.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
      subroutine DLAMC1(BETA, T, RND, IEEE1)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer BETA, T
          logical IEEE1, RND
!
! Local variables
!
          real*8 A, B, C, F, ONE, QTR, SAVEC, T1, T2
          real*8 DLAMC3
          logical FIRST, LIEEE1, LRND
          integer LBETA, LT
          save LBETA, LIEEE1, LRND, LT
          external DLAMC3
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
          data FIRST/.true./
!
          if (FIRST) then
              FIRST = .false.
              ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
              A = 1
              C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
50            continue
              if (C == ONE) then
                  A = 2 * A
                  C = DLAMC3(A, ONE)
                  C = DLAMC3(C, -A)
                  goto 50
              end if
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
              B = 1
              C = DLAMC3(A, B)
!
!+       WHILE( C.EQ.A )LOOP
100           continue
              if (C == A) then
                  B = 2 * B
                  C = DLAMC3(A, B)
                  goto 100
              end if
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
              QTR = ONE / 4
              SAVEC = C
              C = DLAMC3(C, -A)
              LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
              B = LBETA
              F = DLAMC3(B / 2, -B / 100)
              C = DLAMC3(F, A)
              if (C == A) then
                  LRND = .true.
              else
                  LRND = .false.
              end if
              F = DLAMC3(B / 2, B / 100)
              C = DLAMC3(F, A)
              if ((LRND) .and. (C == A)) LRND = .false.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
              T1 = DLAMC3(B / 2, A)
              T2 = DLAMC3(B / 2, SAVEC)
              LIEEE1 = (T1 == A) .and. (T2 > SAVEC) .and. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
              LT = 0
              A = 1
              C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
150           continue
              if (C == ONE) then
                  LT = LT + 1
                  A = A * LBETA
                  C = DLAMC3(A, ONE)
                  C = DLAMC3(C, -A)
                  goto 150
              end if
!+       END WHILE
!
          end if
!
          BETA = LBETA
          T = LT
          RND = LRND
          IEEE1 = LIEEE1
!
!     End of DLAMC1
!
      end subroutine DLAMC1
!*==dlamc2.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
      subroutine DLAMC2(BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer BETA, EMAX, EMIN, T
          real*8 EPS, RMAX, RMIN
          logical RND
!
! Local variables
!
          real*8 A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
              SIXTH, SMALL, THIRD, TWO, ZERO
          real*8 DLAMC3
          logical FIRST, IEEE, IWARN, LIEEE1, LRND
          integer GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, NGNMIN, NGPMIN
          save LBETA, LEMAX, LEMIN, LEPS, LRMAX, LRMIN, LT
          external DLAMC1, DLAMC4, DLAMC5
          !external DLAMC1, DLAMC3, DLAMC4, DLAMC5
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) DOUBLE PRECISION
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) DOUBLE PRECISION
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Save statement ..
!     ..
!     .. Data statements ..
          data FIRST/.true./, IWARN/.false./
!     ..
!     .. Executable Statements ..
!
          if (FIRST) then
              FIRST = .false.
              ZERO = 0
              ONE = 1
              TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
              call DLAMC1(LBETA, LT, LRND, LIEEE1)
!
!        Start to find EPS.
!
              B = LBETA
              A = B**(-LT)
              LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
              B = TWO / 3
              HALF = ONE / 2
              SIXTH = DLAMC3(B, -HALF)
              THIRD = DLAMC3(SIXTH, SIXTH)
              B = DLAMC3(THIRD, -HALF)
              B = DLAMC3(B, SIXTH)
              B = abs(B)
              if (B < LEPS) B = LEPS
!
              LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
50            continue
              if ((LEPS > B) .and. (B > ZERO)) then
                  LEPS = B
                  C = DLAMC3(HALF * LEPS, (TWO**5) * (LEPS**2))
                  C = DLAMC3(HALF, -C)
                  B = DLAMC3(HALF, C)
                  C = DLAMC3(HALF, -B)
                  B = DLAMC3(HALF, C)
                  goto 50
              end if
!+       END WHILE
!
              if (A < LEPS) LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
              RBASE = ONE / LBETA
              SMALL = ONE
              do I = 1, 3
                  SMALL = DLAMC3(SMALL * RBASE, ZERO)
              end do
              A = DLAMC3(ONE, SMALL)
              call DLAMC4(NGPMIN, ONE, LBETA)
              call DLAMC4(NGNMIN, -ONE, LBETA)
              call DLAMC4(GPMIN, A, LBETA)
              call DLAMC4(GNMIN, -A, LBETA)
              IEEE = .false.
!
              if ((NGPMIN == NGNMIN) .and. (GPMIN == GNMIN)) then
                  if (NGPMIN == GPMIN) then
                      LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
                  elseif ((GPMIN - NGPMIN) == 3) then
                      LEMIN = NGPMIN - 1 + LT
                      IEEE = .true.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
                  else
                      LEMIN = min(NGPMIN, GPMIN)
!            ( A guess; no known machine )
                      IWARN = .true.
                  end if
!
              elseif ((NGPMIN == GPMIN) .and. (NGNMIN == GNMIN)) then
                  if (abs(NGPMIN - NGNMIN) == 1) then
                      LEMIN = max(NGPMIN, NGNMIN)
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
                  else
                      LEMIN = min(NGPMIN, NGNMIN)
!            ( A guess; no known machine )
                      IWARN = .true.
                  end if
!
              elseif ((abs(NGPMIN - NGNMIN) == 1) .and. (GPMIN == GNMIN) &
                      ) then
                  if ((GPMIN - min(NGPMIN, NGNMIN)) == 3) then
                      LEMIN = max(NGPMIN, NGNMIN) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
                  else
                      LEMIN = min(NGPMIN, NGNMIN)
!            ( A guess; no known machine )
                      IWARN = .true.
                  end if
!
              else
                  LEMIN = min(NGPMIN, NGNMIN, GPMIN, GNMIN)
!         ( A guess; no known machine )
                  IWARN = .true.
              end if
!**
! Comment out this if block if EMIN is ok
              if (IWARN) then
                  FIRST = .true.
                  write (16, FMT=99001) LEMIN
              end if
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine DLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
              IEEE = IEEE .or. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
              LRMIN = 1
              do I = 1, 1 - LEMIN
                  LRMIN = DLAMC3(LRMIN * RBASE, ZERO)
              end do
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
              call DLAMC5(LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX)
          end if
!
          BETA = LBETA
          T = LT
          RND = LRND
          EPS = LEPS
          EMIN = LEMIN
          RMIN = LRMIN
          EMAX = LEMAX
          RMAX = LRMAX
!
          return
!
99001     format(//' WARNING. The value EMIN may be incorrect:-',          &
                         &        '  EMIN = ', I8,                                           &
                         &        /' If, after inspection, the value EMIN looks',           &
                         &        ' acceptable please comment out ',                        &
                         &        /' the IF block as marked within the code of routine',    &
                         &        ' DLAMC2,', /' otherwise supply EMIN explicitly.',/)
!
!     End of DLAMC2
!
      end subroutine DLAMC2
!*==dlamc3.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
      real*8 function DLAMC3(A, B)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          real*8 A, B
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) DOUBLE PRECISION
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
          DLAMC3 = A + B
!
!
!     End of DLAMC3
!
      end function DLAMC3
!*==dlamc4.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
      subroutine DLAMC4(EMIN, START, BASE)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer BASE, EMIN
          real*8 START
!
! Local variables
!
          real*8 A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
          real*8 DLAMC3
          integer I
          ! external DLAMC3
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAMC4 is a service routine for DLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) DOUBLE PRECISION
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
          A = START
          ONE = 1
          RBASE = ONE / BASE
          ZERO = 0
          EMIN = 1
          B1 = DLAMC3(A * RBASE, ZERO)
          C1 = A
          C2 = A
          D1 = A
          D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
100       continue
          if ((C1 == A) .and. (C2 == A) .and. (D1 == A) .and. (D2 == A)) &
              then
              EMIN = EMIN - 1
              A = B1
              B1 = DLAMC3(A / BASE, ZERO)
              C1 = DLAMC3(B1 * BASE, ZERO)
              D1 = ZERO
              do I = 1, BASE
                  D1 = D1 + B1
              end do
              B2 = DLAMC3(A * RBASE, ZERO)
              C2 = DLAMC3(B2 / RBASE, ZERO)
              D2 = ZERO
              do I = 1, BASE
                  D2 = D2 + B2
              end do
              goto 100
          end if
!+    END WHILE
!
!
!     End of DLAMC4
!
      end subroutine DLAMC4
!*==dlamc5.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
      subroutine DLAMC5(BETA, P, EMIN, IEEE, EMAX, RMAX)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
          real*8 ZERO, ONE
          parameter(ZERO=0.0D0, ONE=1.0D0)
!
! Dummy arguments
!
          integer BETA, EMAX, EMIN, P
          logical IEEE
          real*8 RMAX
!
! Local variables
!
          real*8 DLAMC3
          integer EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
          real*8 OLDY, RECBAS, Y, Z
          ! external DLAMC3
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
          LEXP = 1
          EXBITS = 1
100       continue
          TRY = LEXP * 2
          if (TRY <= (-EMIN)) then
              LEXP = TRY
              EXBITS = EXBITS + 1
              goto 100
          end if
          if (LEXP == -EMIN) then
              UEXP = LEXP
          else
              UEXP = TRY
              EXBITS = EXBITS + 1
          end if
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
          if ((UEXP + EMIN) > (-LEXP - EMIN)) then
              EXPSUM = 2 * LEXP
          else
              EXPSUM = 2 * UEXP
          end if
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
          EMAX = EXPSUM + EMIN - 1
          NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
          if ((mod(NBITS, 2) == 1) .and. (BETA == 2)) EMAX = EMAX - 1
!
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
          if (IEEE) EMAX = EMAX - 1
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
          RECBAS = ONE / BETA
          Z = BETA - ONE
          Y = ZERO
          do I = 1, P
              Z = Z * RECBAS
              if (Y < ONE) OLDY = Y
              Y = DLAMC3(Y, Z)
          end do
          if (Y >= ONE) Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
          do I = 1, EMAX
              Y = DLAMC3(Y * BETA, ZERO)
          end do
!
          RMAX = Y
!
!     End of DLAMC5
!
      end subroutine DLAMC5
!*==lsame.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      logical function LSAME(CA, CB)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          character CA, CB
!
! Local variables
!
          integer INTA, INTB, ZCODE
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     January 31, 1994
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
          LSAME = CA == CB
          if (LSAME) return
!
!     Now test for equivalence if both characters are alphabetic.
!
          ZCODE = ichar('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
          INTA = ichar(CA)
          INTB = ichar(CB)
!
          if (ZCODE == 90 .or. ZCODE == 122) then
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
              if (INTA >= 97 .and. INTA <= 122) INTA = INTA - 32
              if (INTB >= 97 .and. INTB <= 122) INTB = INTB - 32
!
          elseif (ZCODE == 233 .or. ZCODE == 169) then
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
              if (INTA >= 129 .and. INTA <= 137 .or. INTA >= 145 .and. INTA <= 153 .or.&
             &INTA >= 162 .and. INTA <= 169) INTA = INTA + 64
              if (INTB >= 129 .and. INTB <= 137 .or. INTB >= 145 .and. INTB <= 153 .or.&
             &INTB >= 162 .and. INTB <= 169) INTB = INTB + 64
!
          elseif (ZCODE == 218 .or. ZCODE == 250) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
              if (INTA >= 225 .and. INTA <= 250) INTA = INTA - 32
              if (INTB >= 225 .and. INTB <= 250) INTB = INTB - 32
          end if
          LSAME = INTA == INTB
!
!     RETURN
!
!     End of LSAME
!
      end function LSAME
!*==idamax.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      integer function IDAMAX(N, DX, INCX)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer INCX, N
          real*8 DX(*)
!
! Local variables
!
          real*8 DMAX
          integer I, IX
!
!*** End of declarations rewritten by SPAG
!
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
          IDAMAX = 0
          if (N < 1 .or. INCX <= 0) return
          IDAMAX = 1
          if (N == 1) return
          if (INCX == 1) then
!
!        code for increment equal to 1
!
              DMAX = DABS(DX(1))
              do I = 2, N
                  if (DABS(DX(I)) > DMAX) then
                      IDAMAX = I
                      DMAX = DABS(DX(I))
                  end if
              end do
              goto 99999
          end if
!
!        code for increment not equal to 1
!
          IX = 1
          DMAX = DABS(DX(1))
          IX = IX + INCX
          do I = 2, N
              if (DABS(DX(IX)) > DMAX) then
                  IDAMAX = I
                  DMAX = DABS(DX(IX))
              end if
              IX = IX + INCX
          end do
          return
99999     continue
      end function IDAMAX
!*==dcopy.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      subroutine DCOPY(N, DX, INCX, DY, INCY)
          implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
          integer INCX, INCY, N
          real*8 DX(*), DY(*)
!
! Local variables
!
          integer I, IX, IY, M, MP1
!
!*** End of declarations rewritten by SPAG
!
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
          if (N <= 0) return
          if (INCX == 1 .and. INCY == 1) then
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
              M = mod(N, 7)
              if (M /= 0) then
                  do I = 1, M
                      DY(I) = DX(I)
                  end do
                  if (N < 7) return
              end if
              MP1 = M + 1
              do I = MP1, N, 7
                  DY(I) = DX(I)
                  DY(I + 1) = DX(I + 1)
                  DY(I + 2) = DX(I + 2)
                  DY(I + 3) = DX(I + 3)
                  DY(I + 4) = DX(I + 4)
                  DY(I + 5) = DX(I + 5)
                  DY(I + 6) = DX(I + 6)
              end do
          else
!
!        code for unequal increments or equal increments
!          not equal to 1
!
              IX = 1
              IY = 1
              if (INCX < 0) IX = (-N + 1) * INCX + 1
              if (INCY < 0) IY = (-N + 1) * INCY + 1
              do I = 1, N
                  DY(IY) = DX(IX)
                  IX = IX + INCX
                  IY = IY + INCY
              end do
              return
          end if
      end subroutine DCOPY
