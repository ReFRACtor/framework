
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

module vfzmat_Rayleigh_m

   Use vfzmat_Rotation_m
   Use vfzmat_PhasMat_m   , only : vfzmat_PhasMat1

public

contains

subroutine vfzmat_Rayleigh &
   ( max_geoms, max_szas, max_vzas, max_azms,                   & ! input  Dimensions (VLIDORT)
     do_upwelling, do_dnwelling, do_ObsGeoms, Sunlight,         & ! input  Flags
     nstokes, n_geoms, n_szas, n_vzas, n_azms,                  & ! input  Numbers
     offsets, dtr, szas, vzas, azms, obsgeoms, Depol,           & ! Input  Geometries + Depol
     RayFmatrices_up, RayFmatrices_dn, RayZmatrices_up, RayZmatrices_dn, RayCoeffs )

!  Programmed 19 September 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

!  Purpose
!  -------

!  Stand-alone routine to develop F-Matrices, Z-matrices and expansion coefficients for Rayleigh scattering
!    given only the depolarization ratio and input VLIDORT-style geometry

!  Allows for various VLIDORT-based geometry options, upwelling and/or downwelling.

!  Stage 1. (a) Interpolate F-matrices to values implied by geometrical input.
!           (b) Transformation from scattering plane (Fmat) to meridional plane 
!               (Zmat) follows the rotations given in the VLIDORT code.

!  Stage 2. (a) Interpolate F-matrix input to quadrature grid for integration
!           (b) develop Coefficients from integrations using spherical-functions

!  Fmatrix input convention (J is the angle grid on input, N the layer index)
!  ------------------------

!       InFmatrices(N,J,1) = F11(J)
!       InFmatrices(N,J,2) = F22(J)
!       InFmatrices(N,J,3) = F33(J)
!       InFmatrices(N,J,4) = F44(J)
!       InFmatrices(N,J,5) = F12(J)
!       InFmatrices(N,J,6) = F34(J)

!   For Mie scattering, F11 = F22, F33 = F44

!  Convention for using FMatCoeffs output to get VLIDORT Greekmat input
!  --------------------------------------------------------------------

!       GREEKMAT(L,N,1)  --> use + FMatCoeffs(L,N,1)
!       GREEKMAT(L,N,2)  --> use - FMatCoeffs(L,N,5)
!       GREEKMAT(L,N,5)  --> use - FMatCoeffs(L,N,5)
!       GREEKMAT(L,N,6)  --> use + FMatCoeffs(L,N,2)
!       GREEKMAT(L,N,11) --> use + FMatCoeffs(L,N,3)
!       GREEKMAT(L,N,12) --> use - FMatCoeffs(L,N,6)
!       GREEKMAT(L,N,15) --> use + FMatCoeffs(L,N,6)
!       GREEKMAT(L,N,16) --> use + FMatCoeffs(L,N,4)

!    all other entries zero.

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions

   INTEGER  , INTENT(IN) :: max_geoms, max_szas, max_vzas, max_azms

!  Directional flags

   LOGICAL, intent(in)   :: do_upwelling, do_dnwelling

!  Flags for use of observational geometry, Sunlight
 
   LOGICAL, intent(in)   :: do_obsgeoms, Sunlight

!  numbers (general)

   INTEGER, INTENT (IN) :: nstokes

!  Geometry numbers

   INTEGER, INTENT (IN) :: n_geoms, n_szas, n_vzas, n_azms
   INTEGER, INTENT (IN) :: offsets(max_szas,max_vzas)

!  Angles. Convention as for  VLIDORT

   REAL(KIND=dpk), INTENT (IN) :: dtr
   REAL(KIND=dpk), INTENT (IN) :: szas(max_szas)
   REAL(KIND=dpk), INTENT (IN) :: vzas(max_vzas)
   REAL(KIND=dpk), INTENT (IN) :: azms(max_azms)
   REAL(KIND=dpk), INTENT (IN) :: obsgeoms(max_geoms,3)

!  depolarization ratio

   REAL(KIND=dpk), INTENT (IN) :: depol

!  output
!  ------

!  Output Fmatrices (Calculated from Coefficients)

   Real(dpk), INTENT (Out) :: RayFmatrices_up   ( Max_Geoms, 6 )
   Real(dpk), INTENT (Out) :: RayFmatrices_dn   ( Max_Geoms, 6 )

!  Zmatrices

   REAL(KIND=dpk), INTENT (Out) :: RayZmatrices_up(max_geoms,4,4)
   REAL(KIND=dpk), INTENT (Out) :: RayZmatrices_dn(max_geoms,4,4)

!  Fmatrix coefficients

   REAL(KIND=dpk), INTENT (Out) :: RayCoeffs(0:2,6)

!  local variables
!  ---------------

!  Parameters

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk, d_half  = 0.5_dpk, d_two  = 2.0_dpk

!  rotational angles

   REAL(KIND=dpk) :: C1(max_geoms), S1(max_geoms)
   REAL(KIND=dpk) :: C2(max_geoms), S2(max_geoms)

!  Scattering angle cosines

   REAL(KIND=dpk) :: COSSCAT_up(max_geoms)
   REAL(KIND=dpk) :: COSSCAT_dn(max_geoms)

!  Directional sign

   REAL(KIND=dpk) :: vsign

!  Help variables

   INTEGER            :: v, ib, um, ia
   REAL    (KIND=dpk) :: ctheta, stheta, calpha, salpha, cphi
   REAL    (KIND=dpk) :: mu, P2, beta2

!  Debug

   logical, parameter :: check_expansion = .false.

!  Initialize
!  ----------

!  Zero output

   RayFmatrices_up = d_zero
   RayFmatrices_dn = d_zero
   RayZmatrices_up = d_zero
   RayZmatrices_dn = d_zero
   RayCoeffs       = d_zero

!  STAGE 1. Develop Z-matrices
!  ===========================

!  1. Get the scattering angle cosines
!  -----------------------------------

   IF ( .not. Do_Obsgeoms ) THEN
     DO IB = 1, n_szas
       ctheta = cos ( szas(ib) * dtr )
       stheta = sin ( szas(ib) * dtr )
       DO UM = 1, n_vzas
         calpha = cos ( vzas(um) * dtr )
         salpha = sin ( vzas(um) * dtr )
         DO IA = 1, n_azms
           cphi = cos ( azms(ia) * dtr )
           V = OFFSETS(IB,UM) + IA
           COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
           COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
         ENDDO
       ENDDO
     ENDDO
   ELSE
     DO V = 1, n_geoms
       CTHETA = cos(obsgeoms(v,1)*dtr)
       STHETA = sin(obsgeoms(v,1)*dtr)
       CALPHA = cos(obsgeoms(v,2)*dtr)
       SALPHA = sin(obsgeoms(v,2)*dtr)
       CPHI   = cos(obsgeoms(v,3)*dtr)
       COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
       COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
     ENDDO
   ENDIF

!  Set the Rayleigh coefficients

   RayCoeffs      = d_zero
   beta2 = ( d_one - depol ) / ( d_two + depol )
   RayCoeffs(0,1) =  d_one
   RayCoeffs(1,4) =  3.0_dpk * ( d_one - d_two * DEPOL ) / (d_two + DEPOL )
   RayCoeffs(2,1) =  beta2
   RayCoeffs(2,5) =  - SQRT(6.0_dpk) * beta2
   RayCoeffs(2,2) =  6.0_dpk * beta2

!  Construct the F-matrices. Code needs checking.
!   5/25/20. Version 2.8.2 Upgrade. Factor half removed

  if ( do_upwelling ) then
    do V = 1, n_geoms
      mu = COSSCAT_up(v) ; P2 = 1.5_dpk * mu * mu - d_half
!      RayFmatrices_up(v,1) = RayCoeffs(0,1) + d_half * P2 * RayCoeffs(2,1)
      RayFmatrices_up(v,1) = RayCoeffs(0,1) + P2 * RayCoeffs(2,1)
      if ( nstokes .gt. 1 ) then
        RayFmatrices_up(v,2) =   P2 * RayCoeffs(2,2)
        RayFmatrices_up(v,5) = - P2 * RayCoeffs(2,5)
        if ( nstokes.eq.4 ) RayFmatrices_up(v,6) = mu * RayCoeffs(1,4)
      endif
    enddo
  endif

  if ( do_dnwelling ) then
    do V = 1, n_geoms
      mu = COSSCAT_dn(v) ; P2 = 1.5_dpk * mu * mu - d_half
!      RayFmatrices_dn(v,1) = RayCoeffs(0,1) + d_half * P2 * RayCoeffs(2,1)
      RayFmatrices_dn(v,1) = RayCoeffs(0,1) + P2 * RayCoeffs(2,1)
      if ( nstokes .gt. 1 ) then
        RayFmatrices_dn(v,2) = P2 * RayCoeffs(2,2)
        RayFmatrices_dn(v,5) = - P2 * RayCoeffs(2,5)
         if ( nstokes.eq.4 ) RayFmatrices_dn(v,6) = mu * RayCoeffs(1,4)
      endif
    enddo
  endif

!  3. Get the Z-matrices
!  ---------------------

!  upwelling. [ C1/S1/C2/S2 are local, will be overwritten ]

   if ( do_upwelling ) then
      Call vfzmat_Rotation &
       ( max_geoms, max_szas, max_vzas, max_azms, vsign, dtr,   & ! Inputs
         do_ObsGeoms, nstokes, n_geoms, n_szas, n_vzas, n_azms, & ! inputs
         offsets, szas, vzas, azms, obsgeoms,                   & ! Inputs
         C1, S1, C2, S2 )
      Call vfzmat_PhasMat1 &
       ( max_geoms, nstokes, n_geoms, Sunlight, & ! inputs
         C1, S1, C2, S2, RayFmatrices_up,       & ! Inputs
         RayZmatrices_up )
   endif

!  Downwelling

   if ( do_dnwelling ) then
      vsign = + d_one
      Call vfzmat_Rotation &
       ( max_geoms, max_szas, max_vzas, max_azms, vsign, dtr,   & ! Inputs
         do_ObsGeoms, nstokes, n_geoms, n_szas, n_vzas, n_azms, & ! inputs
         offsets, szas, vzas, azms, obsgeoms,                   & ! Inputs
         C1, S1, C2, S2 )
      Call vfzmat_PhasMat1 &
       ( max_geoms, nstokes, n_geoms, Sunlight, & ! inputs
         C1, S1, C2, S2, RayFmatrices_dn,       & ! Inputs
         RayZmatrices_dn )
   endif

!  Done

   return
end subroutine vfzmat_Rayleigh

!  done module

end module vfzmat_Rayleigh_m

