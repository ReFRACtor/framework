
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

      MODULE VLIDORT_LinOutputs_def_m

!  This module contains the following VLIDORT output structures,
!  with intents :

!          VLIDORT_LinCol    nested in VLIDORT_LinOutputs
!         VLIDORT_LinProf    nested in VLIDORT_LinOutputs
!        VLIDORT_LinAtmos    nested in VLIDORT_LinOutputs
!         VLIDORT_LinSurf    nested in VLIDORT_LinOutputs
!      VLIDORT_LinOutputs    Intent(Out)

!  1/31/21. Version 2.8.3. Add MAXLAYERS to this list (for MSST output)

      USE VLIDORT_PARS_m, Only : fpk, MAXSTOKES,  MAXLAYERS, MAX_SURFACEWFS, MAX_ATMOSWFS, MAX_DIRECTIONS, &
                                 MAXLAYERS, MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_LEVELS, MAX_GEOMETRIES

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinCol

!mick mod 9/19/2017 - integrated diffuse quantities now output fully
!                     separately (renamed variables)

!  Atmospheric column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_COLUMNWF

!  Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_MEANST_DIFFUSE_COLWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_FLUX_DIFFUSE_COLWF

!  Direct Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES ) :: TS_DNMEANST_DIRECT_COLWF

!  Direct Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES ) :: TS_DNFLUX_DIRECT_COLWF

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. Special Fourier-component output for Rayleigh + Planetary-problem situations
!           TOA UPWELLING OUTPUT ONLY.
!      REAL(fpk) :: TS_TOAUP_RAYCOLWF_FOURIER ( MAX_ATMOSWFS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, 0:2 )
! #################### INNOVATIONS 5/5/20 ##############

!  4/26/19 Special Media-property output. -- Introduced by R. Spurr.
!     ** LINEARIZED Output for User-angles streams, also fluxes 

      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_VZANGLES, MAX_ATMOSWFS ) :: TS_ALBMED_USER_COLWF
      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_VZANGLES, MAX_ATMOSWFS ) :: TS_TRNMED_USER_COLWF
      REAL(fpk), dimension ( MAXSTOKES, 2, MAX_ATMOSWFS )                 :: TS_ALBMED_FLUXES_COLWF
      REAL(fpk), dimension ( MAXSTOKES, 2, MAX_ATMOSWFS )                 :: TS_TRNMED_FLUXES_COLWF
      REAL(fpk), dimension ( MAXSTOKES, MAX_SZANGLES, MAX_ATMOSWFS )      :: TS_TRANSBEAM_COLWF

!  4/28/19 Special Planetary Problem LINEARIZED output. -- Introduced by R. Spurr.
      
      REAL(fpk), dimension ( MAXSTOKES, MAX_GEOMETRIES, MAX_ATMOSWFS ) :: TS_PLANETARY_TRANSTERM_COLWF
      REAL(fpk), dimension ( MAX_ATMOSWFS )                            :: TS_PLANETARY_SBTERM_COLWF

!  1/31/21. Version 2.8.3. MSST outputs, Column linearizations.
!    -- Revised conditions to include EITHER UPWELLING OR DOWNWELLING (NOT BOTH)
!    -- For the individual layers, also for the surface terms (upwelling only)
!    -- Additional requirement: LC_LOSTRANS  (Path layer transmittances) from the FOCORR Outgoing results.
!    -- Conditions: Upwelling or Downwelling, Observational geometry, FullRad mode, FOCORR and FOCORR Outgoing

     REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAX_SZANGLES, MAXLAYERS )            :: TS_LC_LOSTRANS
     REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAX_SZANGLES, MAXSTOKES, MAXLAYERS ) :: TS_LC_LAYER_MSSTS
     REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAX_SZANGLES, MAXSTOKES )            :: TS_LC_SURF_MSSTS

      END TYPE VLIDORT_LinCol

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinProf

!mick mod 9/19/2017 - integrated diffuse quantities now output fully
!                     separately (renamed variables)

!  Atmospheric profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_PROFILEWF

!  Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_MEANST_DIFFUSE_PROFWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_FLUX_DIFFUSE_PROFWF

!  Direct Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES ) :: TS_DNMEANST_DIRECT_PROFWF

!  Direct Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_SZANGLES, MAXSTOKES ) :: TS_DNFLUX_DIRECT_PROFWF

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. Special Fourier-component output for Rayleigh + Planetary-problem situations
!           TOA UPWELLING OUTPUT ONLY.
!      REAL(fpk) :: TS_TOAUP_RAYPROFWF_FOURIER ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, 0:2 )
! #################### INNOVATIONS 5/5/20 ##############

!  4/26/19 Special Media-property output. -- Introduced by R. Spurr.
!     ** LINEARIZED Output for User-angles streams, also fluxes 

      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_VZANGLES, MAXLAYERS, MAX_ATMOSWFS ) :: TS_ALBMED_USER_PROFWF
      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_VZANGLES, MAXLAYERS, MAX_ATMOSWFS ) :: TS_TRNMED_USER_PROFWF
      REAL(fpk), dimension ( MAXSTOKES, 2, MAXLAYERS, MAX_ATMOSWFS )                 :: TS_ALBMED_FLUXES_PROFWF
      REAL(fpk), dimension ( MAXSTOKES, 2, MAXLAYERS, MAX_ATMOSWFS )                 :: TS_TRNMED_FLUXES_PROFWF
      REAL(fpk), dimension ( MAXSTOKES, MAX_SZANGLES, MAXLAYERS, MAX_ATMOSWFS )      :: TS_TRANSBEAM_PROFWF

!  4/28/19 Special Planetary Problem LINEARIZED output. -- Introduced by R. Spurr.

      REAL(fpk), dimension ( MAXSTOKES, MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS ) :: TS_PLANETARY_TRANSTERM_PROFWF
      REAL(fpk), dimension ( MAXLAYERS, MAX_ATMOSWFS )                            :: TS_PLANETARY_SBTERM_PROFWF

!  1/31/21. Version 2.8.3. MSST outputs, Profile linearizations.
!    -- Revised conditions to include EITHER UPWELLING OR DOWNWELLING (NOT BOTH)
!    -- For the individual layers, also for the surface terms (upwelling only)
!    -- Additional requirement: LP_LOSTRANS  (Path layer transmittances) from the FOCORR Outgoing results.
!    -- Conditions: Upwelling or Downwelling, Observational geometry, FullRad mode, FOCORR and FOCORR Outgoing

     REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAX_SZANGLES, MAXLAYERS )                       :: TS_LP_LOSTRANS
     REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAXLAYERS, MAX_SZANGLES, MAXSTOKES, MAXLAYERS ) :: TS_LP_LAYER_MSSTS
     REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAXLAYERS, MAX_SZANGLES, MAXSTOKES )            :: TS_LP_SURF_MSSTS

      END TYPE VLIDORT_LinProf

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinAtmos

!  Atmospheric thermal emission: Blackbody weighting functions (BBWFS)
!     26 March 2014, Introduced for Version 2.7

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
        0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)           :: TS_ABBWFS_JACOBIANS
      REAL(fpk), dimension ( MAX_USER_LEVELS, 2, &
        0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)           :: TS_ABBWFS_FLUXES

!  These Jacobians were superseded in Version 2.7 by the above ABBWFS arrays
!      LTE linearization - Temperature weighting functions for BB functions (RT Solutions Use Only)
!      REAL(fpk), dimension ( 0:MAXLAYERS, MAX_USER_LEVELS, &
!        MAX_USER_VZANGLES, MAX_DIRECTIONS ) :: TS_LTE_ATMOSWF

      END TYPE VLIDORT_LinAtmos

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinSurf

!  Surface weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_SURFACEWF

!mick mod 9/19/2017 - integrated diffuse quantities renamed

!  Mean intensity weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_MEANST_DIFFUSE_SURFWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_FLUX_DIFFUSE_SURFWF

!  1/31/21. Version 2.8.3. MSST outputs, Surface linearizations.
!    -- Revised conditions to include EITHER UPWELLING OR DOWNWELLING (NOT BOTH)
!    -- For the individual layers, also for the surface terms (upwelling only)
!    -- Conditions: Upwelling or Downwelling, Observational geometry, FullRad mode, FOCORR and FOCORR Outgoing

     REAL(fpk), DIMENSION ( MAX_SURFACEWFS, MAX_SZANGLES, MAXSTOKES, MAXLAYERS ) :: TS_LS_LAYER_MSSTS
     REAL(fpk), DIMENSION ( MAX_SURFACEWFS, MAX_SZANGLES, MAXSTOKES )            :: TS_LS_SURF_MSSTS

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. Special Fourier-component output for Rayleigh + Planetary-problem situations
!           TOA UPWELLING OUTPUT ONLY. Only 1 Surface WF here --> THE ALBEDO !!!!
!      REAL(fpk) :: TS_TOAUP_RAYSURFWF_FOURIER ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES )
! #################### INNOVATIONS 5/5/20 ##############

!  Surface thermal emission: Blackbody weighting functions (BBWFS)
!     26 March 2014, Introduced for Version 2.7

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
                  MAXSTOKES, MAX_DIRECTIONS )  :: TS_SBBWFS_JACOBIANS
      REAL(fpk), dimension ( MAX_USER_LEVELS, 2, &
                  MAXSTOKES, MAX_DIRECTIONS )  :: TS_SBBWFS_FLUXES

      END TYPE VLIDORT_LinSurf

! #####################################################################
! #####################################################################

      TYPE VLIDORT_LinOutputs

      TYPE(VLIDORT_LinCol)   :: Col
      TYPE(VLIDORT_LinProf)  :: Prof
      TYPE(VLIDORT_LinAtmos) :: Atmos
      TYPE(VLIDORT_LinSurf)  :: Surf

      END TYPE VLIDORT_LinOutputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_LinCol,&
                VLIDORT_LinProf,&
                VLIDORT_LinAtmos,&
                VLIDORT_LinSurf,&
                VLIDORT_LinOutputs

      END MODULE VLIDORT_LinOutputs_def_m
