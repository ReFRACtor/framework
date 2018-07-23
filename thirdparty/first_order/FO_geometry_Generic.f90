! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #     FIRST-ORDER SCALAR MODEL (EXACT SINGLE SCATTERING)  #
! #                                                         #
! #  This Version :   1.4 F90                               #
! #  Release Date :   August 2013                           #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3,  29 October  2012, Observation geometry  #
! #   Version 1.4,  31 July     2013, Lattice geometries    #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_geometry_Generic_m

!  Following routines are generic ray-tracing routines

! subroutine FindSun
! subroutine FindSunPaths_D
! subroutine FindSunPaths_T
! subroutine FindSunPath

!  Following is Gaussian-quadrature numerical

!  Subroutine GAULEG_NG

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains


!  General Routines for Sun positioning
!  ++++++++++++++++++++++++++++++++++++

subroutine FindSun(DoNadir,Do_OverheadSun,Radius,SolarDirection,CumAngle,theta_boa_R,theta,stheta,ctheta,DirSun)

!  Find the solar anlge along the LOS path, for given radius and cumulative angle from BOA
!    SolarDirection is defined at BOA, with azimuth relative to the LOS direction.

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs

   logical   , Intent(In)    :: DoNadir,Do_OverheadSun
   real(ffp) , Intent(in)    :: Radius,SolarDirection(3),CumAngle,theta_boa_R

!  Outputs

   real(ffp) , Intent(out)   :: theta,stheta,ctheta
   logical   , Intent(InOut) :: DirSun

!  Local

   real(ffp) :: px(3),b
   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Calculation (Nadir view scenario)

   if ( doNadir ) then
      DirSun = .true.
      theta = theta_boa_R
      ctheta = cos(theta_boa_R)
      stheta = sin(theta_boa_R)
      return
   endif

!  Calculation (overhead sun)

   if ( Do_OverheadSun ) then
      DirSun = .true.
      ctheta = cos(CumAngle)
      stheta = sin(CumAngle)
      theta  = CumAngle
      return
   endif

!  Calculation (General)

   px(1) = - Radius * sin(CumAngle)
   px(2) = zero
   px(3) =   Radius * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta = -b/Radius
   DirSun = ( ctheta.ge.zero )
   stheta = sqrt(one-ctheta*ctheta)
   theta  = acos(ctheta)

!  Done

   return
end subroutine FindSun


subroutine FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                           thstart,sthstart,N,sunpaths)

!  Sunpaths for the Direct-sun illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N
!  Special case = Overhead sun at BOA

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   LOGICAL   , Intent(In)   :: Do_ZeroSunBOA
   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart

!  Output

   real(ffp), Intent(InOut) :: Sunpaths(maxlayers)

!  Local

   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1
   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Layer boundary upper

   N1 = N - 1

!  SBOA condition

   if ( Do_ZeroSunBOA ) then
      sunpaths(n) = radii(n1) - Radstart
      do k = n1, 1, -1
         sunpaths(k) = radii(k-1) - radii(k)
      enddo
      return
   endif

!  First layer

   sth0 = sthstart
   th0  = thstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = asin(sth1)
   ks1  = th0-th1
   sunpaths(n) = sin(ks1)*Radstart/sth1

!  Other layers to TOA

   sth0 = sth1
   th0  = th1
   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = asin(sth1)
      ks1  = th0-th1
      sunpaths(k) = sin(ks1)*radii(k)/sth1
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_D

subroutine FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,thstart,sthstart,N,sunpaths,NT)

!  Sunpaths for the Tangent-height illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart,Pie

!  Output

   INTEGER   , Intent(InOut)   :: NT
   real(ffp), Intent(InOut)    :: Sunpaths(maxlayers)

!  Local

   logical    :: trawl
   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1, tanr

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp
   real(ffp), parameter :: two  = 2.0_ffp

!  Layer boundary upper

   N1 = N - 1

!  tangent height, Find which layer NT

   NT = N
   tanr = sthstart * Radstart
   k = n1 ; trawl = .true.
   do while (k.ge.n1.and.trawl)
      trawl = (radii(k).gt.tanr) ; k = k + 1
   enddo
   nt = k-1 !; write(*,*)n,nt

!  Distances for layers N and below to NT

   if ( nt.gt.n ) then
      th0  = pie - thstart ; sth0 = sthstart
      sth1 = sth0*Radstart/radii(n)
      th1  = asin(sth1) ; ks1  = th0-th1
      sunpaths(n) = two * sin(ks1)*Radstart/sth1
      sth0 = sth1 ; th0 = th1
      do k = n+1,nt-1
        sth1 = sth0*radii(k-1)/radii(k)
        th1  = asin(sth1) ; ks1  = th0-th1
        sunpaths(k) = two * sin(ks1)*radii(k)/sth0
        sth0 = sth1 ; th0 = th1
      enddo
      sth1 = one ; ks1 = 0.5_ffp * pie - th0
      sunpaths(nt) = two * sin(ks1)*radii(nt-1)
   else if ( nt.eq.n ) then
      sunpaths(n) = - two * Radstart * cos(thstart)
   endif

!  Rest of layer n up to the upper boundary

   th0 = pie - thstart ; sth0 = sthstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = asin(sth1) ; ks1  = th0-th1
   sunpaths(n) = sunpaths(n) + sin(ks1)*Radstart/sth1
   sth0 = sth1 ; th0 = th1

!  Trawl up from layers above n, to TOA

   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = asin(sth1)
      ks1  = th0-th1
      sunpaths(k) = sin(ks1)*radii(k)/sth1 
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_T

SUBROUTINE GAULEG_NG(X1,X2,X,W,N,NMAX)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input/Output

      INTEGER  , intent(in)  :: N,NMAX
      REAL(ffp), intent(in)  :: X1, X2
      REAL(ffp), intent(out) :: X(NMAX),W(NMAX)

      INTEGER     :: I, M, J
      REAL(ffp)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(ffp)   :: QUANT, RN, RJ, PIE, ARG

      REAL(ffp), PARAMETER :: EPS = 3.0D-14
      real(ffp), parameter :: zero = 0.0_ffp
      real(ffp), parameter :: one  = 1.0_ffp
      real(ffp), parameter :: two  = 2.0_ffp
      real(ffp), parameter :: half  = 0.5_ffp
      real(ffp), parameter :: qtr   = 0.25_ffp

      M=(N+1)/2
      XM = half * (X2+X1)
      XL = half * (X2-X1)
      RN = real(N,ffp)
      Z1 = zero
      pie = acos(-one)

      DO I=1,M
            arg = ( real(i,ffp) - qtr ) / ( rn + half )
!            Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z  = COS ( pie * arg )
            DO WHILE (ABS(Z-Z1).GT.EPS)
                  P1=one
                  P2=zero
                  DO J=1,N
                        RJ = real(J,ffp)
                        P3=P2
                        P2=P1
                        P1= ( ( two*RJ-one)*Z*P2-(RJ-one)*P3 ) / RJ
                        P1= ( ( two*RJ-one)*Z*P2-(RJ-one)*P3 ) / RJ
                  ENDDO
                  PP=RN*(Z*P1-P2)/(Z*Z-one)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            QUANT = two / ( ( one-Z*Z) * PP *PP )
            X(I)     = XM - XL*Z
            X(N+1-I) = XM + XL*Z
            W(I)     = QUANT * XL
            W(N+1-I) = W(I)
      ENDDO
      RETURN
END SUBROUTINE GAULEG_NG


subroutine FindSunPath  ( x, xboa, rtop, raycon, sundir, sundist, theta0 )

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  I/O

   real(ffp), intent(in)  ::  x, xboa, rtop, raycon, sundir(3)
   real(ffp), intent(out) ::  sundist, theta0

!  Local

   real(ffp) :: xicum, sinx, rad, c0, s0, s1, c1
   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Subroutine for the quick calculation of sunpath from point X on View Path to point at Top of layer

   xicum = xboa - x
   sinx  = sin(x)
   rad   = Raycon / sinx
   c0 = sundir(1) * sin(xicum) - sundir(3) * cos(xicum)
   theta0 = - acos(c0)
   s0 = sqrt(one-c0*c0)
   s1 = s0 * rad / rtop
   c1 = sqrt(one-s1*s1)
   sundist = -rtop * (s1*c0-s0*c1)/s0

!  finish

   return
end subroutine FindSunPath

!  Finish

end module FO_geometry_Generic_m

