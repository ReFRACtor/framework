
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

      MODULE VLIDORT_Outputs_def_m

!  This module contains the following structures:

!                 VLIDORT_Main_Outputs  nested in VLIDORT_Outputs
!           VLIDORT_Exception_Handling  nested in VLIDORT_Outputs
!     VLIDORT_Input_Exception_Handling  Intent(Out) from Input settings
!                      VLIDORT_Outputs  Intent(Out)

!  1/31/21. Version 2.8.3. Add MAXLAYERS to this list (for MSST output)

      USE VLIDORT_PARS_m, only : fpk, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, &
                                 MAXLAYERS, MAXBEAMS, MAXSTOKES, MAXSTREAMS, MAXMOMENTS, MAX_DIRECTIONS, MAX_MESSAGES

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Main_Outputs


!  Fourier values

!      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
!        MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS ) :: TS_STOKES_F

!  Fourier-summed values

      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_STOKES

!  Results for mean-value output:

!  Complete Actinic and Regular Fluxes (including Direct terms)
!mick mod 9/19/2017 - integrated diffuse quantities now output fully
!                     separately (renamed variables)

      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_MEANST_DIFFUSE
      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES, MAX_DIRECTIONS ) :: TS_FLUX_DIFFUSE

!  Direct Fluxes only

      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES ) :: TS_DNMEANST_DIRECT
      REAL(fpk), DIMENSION ( MAX_USER_LEVELS, MAX_SZANGLES, &
        MAXSTOKES ) :: TS_DNFLUX_DIRECT

!  4/26/19 Special Media-property output. -- Introduced by R. Spurr.
!     ** Output for User-angles streams, also fluxes 

      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_VZANGLES ) :: TS_ALBMED_USER,   TS_TRNMED_USER 
      REAL(fpk), dimension ( MAXSTOKES, 2 )                 :: TS_ALBMED_FLUXES, TS_TRNMED_FLUXES
      
!  4/28/19 Special Planetary Problem output. -- Introduced by R. Spurr.

      REAL(fpk), dimension ( MAXSTOKES, MAX_GEOMETRIES ) :: TS_PLANETARY_TRANSTERM
      REAL(fpk)                                          :: TS_PLANETARY_SBTERM

!  1/31/21. Version 2.8.3. MSST outputs
!    -- Revised conditions to include EITHER UPWELLING OR DOWNWELLING (NOT BOTH)
!    -- LOSTRANS and PATHGEOMS are necessary to complete the external sphericity correction.

!    For the individual layers, also for the surface terms (upwelling only)
!    Additional requirement: PATHGEOMS (Path SZA/VZA values in Radians ) from the FOCORR Outgoing Geometry.
!    Additional requirement: LOSTRANS  (Path layer transmittances      ) from the FOCORR Outgoing results.
!    Conditions: Upwelling or Downwelling, Observational geometry, FullRad mode, FOCORR and FOCORR Outgoing

     REAL(fpk), DIMENSION ( 2, 0:MAXLAYERS )                     :: TS_PATHGEOMS
     REAL(fpk), DIMENSION ( MAX_SZANGLES, MAXLAYERS )            :: TS_LOSTRANS
     REAL(fpk), DIMENSION ( MAX_SZANGLES, MAXSTOKES, MAXLAYERS ) :: TS_LAYER_MSSTS
     REAL(fpk), DIMENSION ( MAX_SZANGLES, MAXSTOKES )            :: TS_SURF_MSSTS

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. Special Fourier-component output for Rayleigh + Planetary-problem situations
!           TOA UPWELLING OUTPUT ONLY.
!      REAL(fpk), DIMENSION ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, 0:2 ) :: TS_TOAUP_RAYSTOKES_FOURIER
! #################### INNOVATIONS 5/5/20 ##############

!  Contribution functions (TOA Upwelling only)
!    ==> Fourier components of MS component of Diffuse Field (Not included here)
!    ==> Fourier-summed values of Whole Diffuse Field

!      REAL(fpk), DIMENSION ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  ) :: TS_MS_CONTRIBS_F

      REAL(fpk), DIMENSION ( MAX_GEOMETRIES, MAXSTOKES,  MAXLAYERS ) :: TS_CONTRIBS

!  Fourier numbers used (bookkeeping)

      INTEGER, dimension ( MAX_SZANGLES )  :: TS_FOURIER_SAVED

!  Number of geometries (bookkeeping output)

      INTEGER                              :: TS_N_GEOMETRIES

!  Offsets for indexing geometries
!    -- 1/31/21. Version 2.8.3. (2/16/21) Add the Doublet output offset TS_SZD_OFFSETS

      INTEGER :: TS_SZA_OFFSETS ( MAX_SZANGLES )
      INTEGER :: TS_VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )
      INTEGER :: TS_SZD_OFFSETS ( MAX_SZANGLES )

!  Solar Beam Transmittance to BOA
!  rob fix 11/17/2014, for diagnostic use only

      REAL(fpk), dimension ( MAXBEAMS ) :: TS_SOLARBEAM_BOATRANS

!  Contribution function (TOA Upwelling only)

!      REAL(fpk), dimension ( MAX_GEOMETRIES, MAXSTOKES, &
!        MAXLAYERS ) :: TS_SS_CONTRIBS


      END TYPE VLIDORT_Main_Outputs

! #####################################################################
! #####################################################################

type VLIDORT_WLAdjusted_Outputs

!  This is the Water-leaving output after adjustment with the diffuse-flux transmittance
!  ( downwelling at the ocean surface) in order to get proper coupling between atmospheric RT
!  and water-leaving sources. Introduced for Version 2.8.1, 3/18/19 by R. Spurr.
!   -  Output is controlled by the flag DO_WLADJUSTED_OUTPUT
   
!  Isotropic Surface leaving term

      REAL(fpk), dimension ( MAXSTOKES, MAXBEAMS ) :: TS_WLADJUSTED_ISOTROPIC

!  Direct Surface-Leaving term

      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAXBEAMS ) :: TS_WLADJUSTED_DIRECT

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )        :: TS_WLADJUSTED_F_Ords_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_VZANGLES, MAXBEAMS ) :: TS_WLADJUSTED_F_User_0

end type VLIDORT_WLAdjusted_Outputs

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Exception_Handling


!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: TS_STATUS_INPUTCHECK
      INTEGER      :: TS_NCHECKMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: TS_CHECKMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: TS_ACTIONS

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER             :: TS_STATUS_CALCULATION
      CHARACTER (Len=120) :: TS_MESSAGE, TS_TRACE_1, TS_TRACE_2, TS_TRACE_3


      END TYPE VLIDORT_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Input_Exception_Handling


!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: TS_STATUS_INPUTREAD
      INTEGER      :: TS_NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: TS_INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: TS_INPUTACTIONS


      END TYPE VLIDORT_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Outputs

!   -- Add the WLAdjusted_Outputs type structure, Version 2.8.1, 3/18/19

      TYPE(VLIDORT_Main_Outputs)       :: Main
      TYPE(VLIDORT_WLAdjusted_Outputs) :: WLOut
      TYPE(VLIDORT_Exception_Handling) :: Status

      END TYPE VLIDORT_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE
!   -- Add the WLAdjusted_Outputs type structure, Version 2.8.1, 3/18/19

      PRIVATE
      PUBLIC :: VLIDORT_Main_Outputs, &
                VLIDORT_WLAdjusted_Outputs, &
                VLIDORT_Exception_Handling, &
                VLIDORT_Input_Exception_Handling, &
                VLIDORT_Outputs

      END MODULE VLIDORT_Outputs_def_m
