
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! #  This is vfzmat_Numerical.f90.  Utility routines.           #
! #  The subroutines in vfzmat_Numerical are listed below with  #
! #     their source of origin (order of appearance).           #
! #                                                             #
! #     vfzmat_Interp                                           #
! #     vfzmat_Interp_2                                         #
! #     GETQUAD2                     M. Christi, 2017           #
! #     BSPLINE, SEVAL, DSEVAL       Public Bspline software    #
! #                                  (Adapted by R. Spurr 2017) #
! #                                                             #
! ###############################################################

Module vfzmat_Numerical_m

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

      USE VLIDORT_PARS_m, only : fpk, zero, one, two, half, quarter, pie

private
!mick mod 9/19/2017 - use new numerical subroutines 
public  vfzmat_Interp, vfzmat_Interp_2, &
        GETQUAD2, BSPLINE, SEVAL, DSEVAL

contains

subroutine vfzmat_Interp &
    ( Max_InAngles, Max_OutAngles, do_Quadrature, dtr, &
      N_OutAngles, N_InAngles, InAngles, InScatMat,    &
      OutAngles, OutCosines, OutWeights, OutScatMat )

!  Stand-alone routine to interpret Scattering Matrix onto Quadrature grid

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensioning

   Integer, intent(in) :: Max_InAngles,  Max_OutAngles

!  Logical for Quadrature

   Logical, intent(in) :: do_Quadrature

!  Angle numbers

   Integer, intent(in) :: N_InAngles,  N_OutAngles

!  Input angles and scattering matrix entries

   Real(dpk), intent(in) :: dtr, InAngles ( Max_InAngles )
   Real(dpk), intent(in) :: InScatMat ( Max_InAngles, 6 )

!  Outputs
!  -------

!   Quadrature

   Real(dpk), intent(inout) :: OutAngles ( Max_OutAngles )
   Real(dpk), intent(inout) :: OutCosines ( Max_OutAngles )
   Real(dpk), intent(inout) :: OutWeights ( Max_OutAngles )

!  Interpolated scattering matrix entries

   Real(dpk), intent(out) :: OutScatMat ( Max_OutAngles, 6 )

!  Local
!  -----

   real(dpk), parameter :: zero = 0.0_dpk, one = 1.0_dpk, mone = -1.0_dpk
   logical    :: loop
   integer    :: n, m, k, k1, k2, kstart
   real(dpk)  :: InCosines(Max_InAngles), x, f1, f2, x1, x2, w1, w2

!  Start of Code
!  -------------

!  initialize

   OutScatMat = zero
   if ( do_Quadrature ) then
      OutAngles  = zero
      OutCosines = zero
      OutWeights = zero
   endif

!  Develop Quadrature

   if ( do_Quadrature ) then
      Call GETQUAD2 ( mone, one, N_OutAngles, OutCosines, OutWeights )
      do k1 = 1, N_OutAngles / 2
         k2 = N_OutAngles + 1 - k1 
         x1 = OutCosines(k1) ; x2 = OutCosines(k2)
         w1 = OutWeights(k1) ; w2 = OutWeights(k2)
         OutCosines(k1) = x2 ; OutCosines(k2) = x1
         OutWeights(k1) = w2 ; OutWeights(k2) = w1
      enddo
      do k = 1, N_OutAngles
         OutAngles(k) = acos ( OutCosines(k) ) / dtr
      enddo
   endif

!  Cosines of input angles

   do k = 1, N_InAngles
      InCosines(k) = cos ( InAngles(k) * dtr )
   enddo

!  interpolate

   kstart = 0
   do n = 1, N_OutAngles
      x = OutCosines(n)
      k = kstart ; loop = .true.
      do while (loop .and.k.lt.N_InAngles)
         k = k + 1 
         loop = (x.lt.InCosines(k))
      enddo
      kstart = k - 1 ; k1 = k - 1 ; k2 = k
      f1 = ( InCosines(k2) - x ) / ( InCosines(k2) - InCosines(k1) ) ; f2 = one - f1
      do m = 1, 6
         OutScatmat(n,m) = f1 * InScatmat(k1,m) + f2 * InScatmat(k2,m)
      enddo
   enddo

!  Done

   return
end subroutine vfzmat_Interp

!

subroutine vfzmat_Interp_2 &
    ( Max_InAngles, Max_OutAngles, do_Quadrature, dtr, &
      N_OutAngles, N_InAngles, InAngles, InScatMat,    &
      OutAngles, OutCosines, OutWeights, OutScatMat )

!  Stand-alone routine to interpret Scattering Matrix onto Quadrature grid
!  This one based on Splines

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensioning

   Integer, intent(in) :: Max_InAngles,  Max_OutAngles

!  Logical for Quadrature

   Logical, intent(in) :: do_Quadrature

!  Angle numbers

   Integer, intent(in) :: N_InAngles,  N_OutAngles

!  Input angles and scattering matrix entries

   Real(dpk), intent(in) :: dtr, InAngles ( Max_InAngles )
   Real(dpk), intent(in) :: InScatMat ( Max_InAngles, 6 )

!  Outputs
!  -------

!   Quadrature

   Real(dpk), intent(inout) :: OutAngles  ( Max_OutAngles )
   Real(dpk), intent(inout) :: OutCosines ( Max_OutAngles )
   Real(dpk), intent(inout) :: OutWeights ( Max_OutAngles )

!  Interpolated scattering matrix entries

   Real(dpk), intent(out) :: OutScatMat ( Max_OutAngles, 6 )

!  Local
!  -----

   real(dpk), parameter :: zero = 0.0_dpk, one = 1.0_dpk, mone = -1.0_dpk
   integer    :: m, k, k1, k2
   real(dpk)  :: InCosines(Max_InAngles), Local_InScatMat(Max_InAngles,6)
   !real(dpk)  :: x1,x2,w1,w2,x,y,y2(Max_InAngles,6),yp1,ypn
   real(dpk)  :: x1,x2,w1,w2,x,y
   real(dpk)  :: bbas( Max_InAngles ),cbas( Max_InAngles ),dbas( Max_InAngles )

!  Start of Code
!  -------------

!  initialize

   OutScatMat = zero
   if ( do_Quadrature ) then
      OutAngles  = zero
      OutCosines = zero
      OutWeights = zero
   endif

!  Develop Quadrature

   if ( do_Quadrature ) then
      Call GETQUAD2 ( mone, one, N_OutAngles, OutCosines, OutWeights )
      do k1 = 1, N_OutAngles / 2
         k2 = N_OutAngles + 1 - k1 
         x1 = OutCosines(k1) ; x2 = OutCosines(k2)
         w1 = OutWeights(k1) ; w2 = OutWeights(k2)
         OutCosines(k1) = x2 ; OutCosines(k2) = x1
         OutWeights(k1) = w2 ; OutWeights(k2) = w1
      enddo
      do k = 1, N_OutAngles
         OutAngles(k) = acos ( OutCosines(k) ) / dtr
      enddo
   endif

!  Cosines of input angles.
!    -- Reverse directions for the Splining (Monotonically increasing)

   do k = 1, N_InAngles
     k1 = N_InAngles + 1 - k 
     InCosines(k1) = cos ( InAngles(k) * dtr )
     Local_Inscatmat(k1,1:6) = Inscatmat(k,1:6)  
   enddo

!  Develop Splines
!    - Set the End-point gradient (YPN) = input gradient at forward-peak
!    - This make a HUGE difference to the accuracy

    do m = 1, 6
       !yp1 = zero
       !ypn = ( Local_InScatMat(N_InAngles,m) -  Local_InScatMat(N_InAngles-1,m) ) / &
       !       ( InCosines(N_InAngles) -  InCosines(N_InAngles-1) )
       !Call vfzmat_SPLINE(Max_InAngles,InCosines,Local_InScatMat(:,m),N_InAngles,yp1,ypn,y2(:,m))
       Call BSpline (Max_InAngles,N_InAngles,InCosines,Local_InScatMat(:,m),bbas,cbas,dbas)
    enddo

!  interpolate

   do k = 1, N_OutAngles
      k1 = N_OutAngles + 1 - k 
      x = OutCosines(k1)
      do m = 1, 6
         !Call vfzmat_SPLINT(Max_InAngles,InCosines,Local_InScatMat(:,m),y2(:,m),N_InAngles,x,y)
         Call Seval (Max_InAngles,N_InAngles,x,InCosines,Local_InScatMat(:,m),bbas,cbas,dbas,y)
         OutScatMat(k1,m) = y
      enddo
   enddo

!  Done

   return
end subroutine vfzmat_Interp_2

!

      SUBROUTINE GETQUAD2(A,B,N,ROOTS,WGTS)

!  Computes N roots and weights for Gauss-Legendre quadrature on the interval (a,b)

      IMPLICIT NONE

!   Precision

      INTEGER, PARAMETER :: DPK = SELECTED_REAL_KIND(15)

!  Limits of interval

      REAL(DPK), INTENT(IN)  :: A, B

!  Dimension

      INTEGER, INTENT(IN) :: N

!  Quadrature roots and weights

      REAL(DPK), INTENT(OUT) :: ROOTS(N), WGTS(N)

!  Local variables

      INTEGER   :: I, M, N2, NM1
      REAL(DPK) :: IR, MR, NR
      REAL(DPK) :: MIDPT, SFAC
      REAL(DPK) :: DLP_DX, LP, LPM1, LPM2, X, XOLD, XX

!  Threshold for Newton's Method

      REAL(DPK), PARAMETER :: QEPS = 1.0D-13

!  Since roots are symmetric about zero on the interval (-1,1), split the interval
!  in half and only work on the lower half of the interval (-1,0).

      N2 = INT((N + 1)/2)
      NR = REAL(N,DPK)

!  Define the shift [midpoint of (a,b)] and scale factor to later move roots from
!  the interval (-1,1) to the interval (a,b)

      MIDPT = HALF*(B + A)
      SFAC  = HALF*(B - A)

      DO M = 1, N2

!  Find current root of the related Nth order Legendre Polynomial on (-1,0) by Newton's
!  Method using two Legendre Polynomial recurrence relations (e.g. see Abramowitz &
!  Stegan (1972))

         !Define starting point [ after Tricomi (1950) ]
         MR = REAL(M,DPK)
         XX = PIE*(MR - QUARTER)/(NR + HALF)
         X  = (ONE - (NR - ONE)/(8.0_DPK*NR**3) &
             - ONE/(384.0_DPK*NR**4)*(39.0_DPK - 28.0_DPK/SIN(XX)**2))*COS(XX)

         !Use Newton's Method
         DO 
            LPM1 = ZERO ; LP = ONE
            DO I = 1, N
               IR = REAL(I,DPK) ; LPM2 = LPM1 ; LPM1 = LP
               LP = ((TWO*IR - ONE)*X*LPM1 - (IR - ONE)*LPM2)/IR
            ENDDO
            DLP_DX = NR*(X*LP - LPM1)/(X**2 - ONE)
            XOLD = X ; X = XOLD - LP/DLP_DX
            IF (ABS(X-XOLD) <= QEPS) EXIT
         ENDDO

!  Shift and scale the current root (and its symmetric counterpart) from the interval (-1,1)
!  to the interval (a,b).  Define their related weights (e.g. see Abramowitz & Stegan (1972)).
!  Note:
!  If (1) N is even or (2) N is odd and M /= N2, then ROOTS(M) and ROOTS(NM1) are unique.
!  If N is odd and M = N2, then M = NM1 and ROOTS(M) = ROOTS(NM1) are one and the same root.

         !On interval lower half: (a,midpt)
         ROOTS(M)   = MIDPT - SFAC*X
         WGTS(M)    = (TWO*SFAC)/((ONE - X**2)*DLP_DX**2)

         !On interval upper half: (midpt,b)
         NM1 = N - M + 1
         ROOTS(NM1) = MIDPT + SFAC*X
         WGTS (NM1) = WGTS(M)

      ENDDO

      END SUBROUTINE GETQUAD2

!

      SUBROUTINE SEVAL (NMAX,N,U,X,Y,B,C,D,VAL)

!------------------------------------------------------------------------
!     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!     BY THE SPLINE SUBROUTINE.
!
!     INPUTS:
!     N       NUMBER OF POINTS OF CURVE Y = F(X)
!     U       ABSCISSA OF POINT TO BE INTERPOLATED
!     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!             OF CURVE F(X)
!     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!             CUBIC SPLINE
!
!     OUTPUTS:
!     VAL   INTERPOLATED VALUE
!             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!------------------------------------------------------------------------

      implicit none
      integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Arguments

      integer  , intent(in) :: NMAX,N
      REAL(dpk), intent(in) :: B(NMAX),C(NMAX),D(NMAX),X(NMAX),Y(NMAX),U
      REAL(dpk), intent(out):: VAL

!  Local
     
      integer   :: I, J, K
      real(dpk) :: DX
      DATA I/1/

!     BINARY SEARCH

      IF (I.GE.N) I = 1
      IF (U.LT.X(I))   GO TO 10
      IF (U.LE.X(I+1)) GO TO 30
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF (U.LT.X(K)) J = K
      IF (U.GE.X(K)) I = K
      IF (J.GT.I+1) GO TO 20

!     SPLINE EVALUATION

   30 DX  = U-X(I)
      VAL = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))

      RETURN
      END SUBROUTINE SEVAL

!

      SUBROUTINE DSEVAL (NMAX,N,U,X,Y,B,C,D,VAL,DVAL)

!------------------------------------------------------------------------
!     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!     BY THE SPLINE SUBROUTINE.
!
!     INPUTS:
!     N       NUMBER OF POINTS OF CURVE Y = F(X)
!     U       ABSCISSA OF POINT TO BE INTERPOLATED
!     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!             OF CURVE F(X)
!     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!             CUBIC SPLINE
!
!     OUTPUTS:
!     VAL   INTERPOLATED VALUE
!             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!------------------------------------------------------------------------

      implicit none
      integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Arguments

      integer  , intent(in) :: NMAX,N
      REAL(dpk), intent(in) :: B(NMAX),C(NMAX),D(NMAX),X(NMAX),Y(NMAX),U
      REAL(dpk), intent(out):: VAL, DVAL

!  Local
     
      integer   :: I, J, K
      real(dpk) :: DX
      DATA I/1/

!     BINARY SEARCH

      IF (I.GE.N) I = 1
      IF (U.LT.X(I))   GO TO 10
      IF (U.LE.X(I+1)) GO TO 30
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF (U.LT.X(K)) J = K
      IF (U.GE.X(K)) I = K
      IF (J.GT.I+1) GO TO 20

!     SPLINE EVALUATION

   30 CONTINUE
      DX  = U-X(I)
      VAL  = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      DVAL = B(I) + DX*(2.0_dpk*C(I)+3.0_dpk*DX*D(I))

      RETURN
      END SUBROUTINE DSEVAL

!

      SUBROUTINE BSPLINE (NMAX,N,X,Y,B,C,D)

!---------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
!     SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
!
!     INPUTS:
!     N       NUMBER OF GIVEN POINTS
!     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!     OUTPUTS:
!     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!
!     REFERENCE:
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!---------------------------------------------------------------------

      implicit none
      integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Arguments

      integer  , intent(in)  :: NMAX,N
      REAL(dpk), intent(in)  :: X(NMAX),Y(NMAX)
      REAL(dpk), intent(out) :: B(NMAX),C(NMAX),D(NMAX)

!  Local
     
      integer   :: NM1, L, I
      real(dpk) :: T

      NM1 = N-1
      IF (N.LT.2) RETURN
      IF (N.LT.3) GO TO 50

!     BUILD THE TRIDIAGONAL SYSTEM
!     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)

      D(1) = X(2)-X(1)
      C(2) = (Y(2)-Y(1))/D(1)
      DO 10 I = 2,NM1
         D(I) = X(I+1)-X(I)
         B(I) = 2.0_dpk*(D(I-1)+D(I))
         C(I+1) = (Y(I+1)-Y(I))/D(I)
         C(I) = C(I+1)-C(I)
   10 CONTINUE

!     CONDITIONS AT LIMITS
!     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES

      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.0_dpk
      C(N) = 0.0_dpk
      IF (N.EQ.3) GO TO 15
      C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))

!     FORWARD ELIMINATION

   15 CONTINUE

      DO 20 I = 2,N
        T = D(I-1)/B(I-1)
        B(I) = B(I)-T*D(I-1)
        C(I) = C(I)-T*C(I-1)
   20 CONTINUE

!     BACK SUBSTITUTION

      C(N) = C(N)/B(N)
      DO 30 L = 1,NM1
        I = N-L
        C(I) = (C(I)-D(I)*C(I+1))/B(I)
   30 CONTINUE

!     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL

      B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.0_dpk*C(N))
      DO 40 I = 1,NM1
        B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.0_dpk*C(I))
        D(I) = (C(I+1)-C(I))/D(I)
        C(I) = 3.0_dpk*C(I)
   40 CONTINUE
      C(N) = 3.0_dpk*C(N)
      D(N) = D(NM1)
      RETURN

!     CAS N = 2

   50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.0_dpk
      D(1) = 0.0_dpk
      B(2) = B(1)
      C(2) = 0.0_dpk
      D(2) = 0.0_dpk

      RETURN
      END SUBROUTINE BSPLINE

!  End module

End Module vfzmat_Numerical_m

