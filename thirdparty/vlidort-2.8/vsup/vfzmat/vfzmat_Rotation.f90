
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

module vfzmat_Rotation_m

public

contains

subroutine vfzmat_Rotation &
   ( max_geoms, max_szas, max_vzas, max_azms, vsign, dtr,   & ! Inputs
     do_ObsGeoms, nstokes, n_geoms, n_szas, n_vzas, n_azms, & ! inputs
     offsets, szas, vzas, azms, obsgeoms,                   & ! Inputs
     C1, S1, C2, S2 )

!  Stand-alone routine to develop Rotation angles for Z-matrices 

!  Programmed 19 September 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

!  Transformation from scattering plan (Fmat) to Meriodional Plane (Zmat)
!   follows the rotations given in the VLIDORT code.

!  Remarks
!  -------

!  1. vsign = +1 for upwelling, -1 for downwelling
!  2. Angles and numbers should be set up beforehand as for VLIDORT.

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input
!  -----

!  dimensions

   INTEGER, INTENT (IN) :: max_geoms, max_szas, max_vzas, max_azms

!  numbers

   INTEGER, INTENT (IN) :: nstokes
   INTEGER, INTENT (IN) :: n_geoms, n_szas, n_vzas, n_azms
   INTEGER, INTENT (IN) :: offsets(max_szas,max_vzas)

!  Flags for use of observational geometry
 
   LOGICAL, intent(in)   :: do_obsgeoms

!  Direction sign

   REAL(KIND=dpk), INTENT (IN) :: vsign, dtr

!  Angles. Convention as for  VLIDORT

   REAL(KIND=dpk), INTENT (IN) :: szas(max_szas)
   REAL(KIND=dpk), INTENT (IN) :: vzas(max_vzas)
   REAL(KIND=dpk), INTENT (IN) :: azms(max_azms)
   REAL(KIND=dpk), INTENT (IN) :: obsgeoms(max_geoms,3)

!  output
!  ------

!  rotational-angle output

   REAL(KIND=dpk), INTENT (Out)   :: C1(max_geoms), S1(max_geoms)
   REAL(KIND=dpk), INTENT (Out)   :: C2(max_geoms), S2(max_geoms)

!  local variables
!  ---------------

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   real(dpk), parameter :: d_half  = 0.5_dpk, d_two  = 2.0_dpk

   INTEGER            :: v, ib, um, ia
   REAL    (KIND=dpk) :: ctheta, stheta, calpha, salpha, cphi
   REAL    (KIND=dpk) :: COSSCAT, SINSCAT, HELP_SINSCAT
   REAL    (KIND=dpk) :: CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

!  Initialization

   C1 = d_zero ; C2 = d_zero
   S1 = d_zero ; S2 = d_zero

!  no Angles if scalar

   if ( nstokes == 1 ) return

!  Lattice geometry
!  ----------------

   IF ( .not. Do_Obsgeoms ) THEN

!  Geometry loops

     DO IB = 1, n_szas
       ctheta = cos ( szas(ib) * dtr )
       stheta = sin ( szas(ib) * dtr )
       DO UM = 1, n_vzas
         calpha = cos ( vzas(um) * dtr )
         salpha = sin ( vzas(um) * dtr )
         DO IA = 1, n_azms
           cphi = cos ( azms(ia) * dtr )
           V = OFFSETS(IB,UM) + IA

!  cosine scatter angle (this is valid only for non-refracting atmosphere
!  VSIGN = -1 for upwelling, +1 for downwelling

           COSSCAT = VSIGN * CTHETA * CALPHA + STHETA * SALPHA * CPHI

!  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

!  a. safety
!    Watch for sin^2(scatter angle) less than zero (machine precision)
!    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

           HELP_SINSCAT = ( D_ONE - COSSCAT * COSSCAT )
           IF ( HELP_SINSCAT.LE.D_ZERO ) THEN
             SINSCAT = 1.0E-12_dpk
           ELSE
             SINSCAT = SQRT ( HELP_SINSCAT )
           ENDIF

!  b. necessary limit analyses - Hovenier limits.
!     R. Spurr and V. Natraj, 17 January 2006

           IF ( ABS(SINSCAT) .LE. 1.0E-12_dpk ) THEN
             CSIG1 = D_ZERO ; CSIG2 = D_ZERO
           ELSE
             IF ( STHETA .EQ. D_ZERO ) THEN
               CSIG1 = -  CPHI
             ELSE
               CSIG1   = ( - VSIGN * CALPHA + CTHETA*COSSCAT ) / SINSCAT / STHETA
             ENDIF
             IF ( SALPHA .EQ. D_ZERO ) THEN
               CSIG2 = -  CPHI
             ELSE
               CSIG2   = ( - CTHETA + VSIGN*CALPHA*COSSCAT ) / SINSCAT / SALPHA
             ENDIF
           ENDIF

!  These lines are necessary to avoid bad values

           IF ( CSIG1 .GT. D_ONE  ) CSIG1 = D_ONE
           IF ( CSIG1 .LT. -D_ONE ) CSIG1 = -D_ONE

           IF ( CSIG2 .GT. D_ONE  ) CSIG2 = D_ONE
           IF ( CSIG2 .LT. -D_ONE ) CSIG2 = -D_ONE

!  output, H/VdM, Eqs. (89)-(94)
!    Rotation sines and cosines. Same for all layers

           CSIG1_2 = D_TWO * CSIG1
           CSIG2_2 = D_TWO * CSIG2
           IF ( ABS(CSIG1-D_ONE).LT.1.0E-12_dpk)THEN
             SSIG1 = D_ZERO
           ELSE
             SSIG1 = SQRT ( D_ONE - CSIG1 * CSIG1 )
           ENDIF
           IF ( ABS(CSIG2-D_ONE).LT.1.0E-12_dpk)THEN
             SSIG2 = D_ZERO
           ELSE
             SSIG2 = SQRT ( D_ONE - CSIG2 * CSIG2 )
           ENDIF

!  For relazm in [180,360), need sign reversal for S1 and S2
!  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

           C1(V) = CSIG1_2 * CSIG1 - D_ONE
           C2(V) = CSIG2_2 * CSIG2 - D_ONE

           IF (AZMS(IA) .LE. 180.0_dpk) THEN
             S1(V) = CSIG1_2 * SSIG1
             S2(V) = CSIG2_2 * SSIG2
           ELSE
             S1(V) = -CSIG1_2 * SSIG1
             S2(V) = -CSIG2_2 * SSIG2
           ENDIF

!  End geometry loops

         ENDDO
       ENDDO
     ENDDO

!  End lattice geometry

   ENDIF

!  Observational geometry
!  ----------------------

   IF ( Do_Obsgeoms ) then

!  Start loop geometry

     DO V = 1, n_geoms

       CTHETA = cos(obsgeoms(v,1)*dtr)
       STHETA = sin(obsgeoms(v,1)*dtr)
       CALPHA = cos(obsgeoms(v,2)*dtr)
       SALPHA = sin(obsgeoms(v,2)*dtr)
       CPHI   = cos(obsgeoms(v,3)*dtr)

!  cosine scatter angle (this is valid only for non-refracting atmosphere
!  VSIGN = -1 for upwelling, +1 for downwelling

       COSSCAT = VSIGN * CTHETA * CALPHA + STHETA * SALPHA * CPHI

!  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

!  a. safety
!    Watch for sin^2(scatter angle) less than zero (machine precision)
!    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

       HELP_SINSCAT = ( D_ONE - COSSCAT * COSSCAT )
       IF ( HELP_SINSCAT.LE.D_ZERO ) THEN
         SINSCAT = 1.0E-12_dpk
       ELSE
         SINSCAT = SQRT ( HELP_SINSCAT )
       ENDIF

!  b. necessary limit analyses - Hovenier limits.
!     R. Spurr and V. Natraj, 17 January 2006

       IF ( ABS(SINSCAT) .LE. 1.0E-12_dpk ) THEN
         CSIG1 = D_ZERO ; CSIG2 = D_ZERO
       ELSE
         IF ( STHETA .EQ. D_ZERO ) THEN
           CSIG1 = -  CPHI
         ELSE
           CSIG1   = ( - VSIGN * CALPHA + CTHETA*COSSCAT ) / SINSCAT / STHETA
         ENDIF
         IF ( SALPHA .EQ. D_ZERO ) THEN
           CSIG2 = -  CPHI
         ELSE
           CSIG2   = ( - CTHETA + VSIGN*CALPHA*COSSCAT ) / SINSCAT / SALPHA
         ENDIF
       ENDIF

!  These lines are necessary to avoid bad values

       IF ( CSIG1 .GT. D_ONE  ) CSIG1 = D_ONE
       IF ( CSIG1 .LT. -D_ONE ) CSIG1 = -D_ONE

       IF ( CSIG2 .GT. D_ONE  ) CSIG2 = D_ONE
       IF ( CSIG2 .LT. -D_ONE ) CSIG2 = -D_ONE

!  output, H/VdM, Eqs. (89)-(94)
!    Rotation sines and cosines. Same for all layers

       CSIG1_2 = D_TWO * CSIG1
       CSIG2_2 = D_TWO * CSIG2
       IF ( ABS(CSIG1-D_ONE).LT.1.0E-12_dpk)THEN
         SSIG1 = D_ZERO
       ELSE
         SSIG1 = SQRT ( D_ONE - CSIG1 * CSIG1 )
       ENDIF
       IF ( ABS(CSIG2-D_ONE).LT.1.0E-12_dpk)THEN
         SSIG2 = D_ZERO
       ELSE
         SSIG2 = SQRT ( D_ONE - CSIG2 * CSIG2 )
       ENDIF

!  For relazm in [180,360), need sign reversal for S1 and S2
!  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

       C1(V) = CSIG1_2 * CSIG1 - D_ONE
       C2(V) = CSIG2_2 * CSIG2 - D_ONE

       IF (Obsgeoms(V,3) .LE. 180.0_dpk) THEN
         S1(V) = CSIG1_2 * SSIG1
         S2(V) = CSIG2_2 * SSIG2
       ELSE
         S1(V) = -CSIG1_2 * SSIG1
         S2(V) = -CSIG2_2 * SSIG2
       ENDIF

!  End geometry loop

     ENDDO

!  End Observational Geometry

   ENDIF

!  Done

   RETURN
END SUBROUTINE vfzmat_Rotation

!  End module

End module vfzmat_Rotation_m

