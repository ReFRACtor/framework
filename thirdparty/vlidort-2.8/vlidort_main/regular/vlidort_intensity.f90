
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
! # Subroutines in this Module                                  #
! #                                                             #
! #           VLIDORT_UPUSER_INTENSITY (master)                 #
! #           VLIDORT_DNUSER_INTENSITY (master)                 #
! #                                                             #
! #              GET_TOA_SOURCE                                 #
! #              GET_BOA_SOURCE                                 #
! #                                                             #
! #           VLIDORT_INTEGRATED_OUTPUT (master)                #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. Green's Function changes
!    -- Need additional Taylor routines to be used, for Green's function post-processing
!    -- Whole/part Sourceterm routines completely rewritten       (Green's function included)
!    -- QuadIntens partial-layer subroutines completely rewritten (Green's function included)
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier

!  1/31/21. Version 2.8.3. Other post-processing changes
!    -- Converge routines have been moved to their own module (vlidort_converge.f90)
!    -- MSST option is now included, generates output for sphericity correctio
!    -- DO_DOUBLET_GEOMETRY post-processing operation is in force.
!    -- Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation vs. lattice/doublet geometry choices
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier
!    -- Quadrature and Post-processing sourceterm subroutines moved to their own module.
!    -- (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES) number of STOKES components

      MODULE vlidort_intensity_m

!  1/31/21. Version 2.8.3. 
!    ==> Need Taylor series routines for Green's function post-processing
!    ==> Source term postprocessing routines now have their own module.

      USE vlidort_Taylor_m, only : TAYLOR_SERIES_1, TAYLOR_SERIES_2
      USE vlidort_PostProcessing_m

!  1/31/21. Version 2.8.3.  ==> Only GET_TOA_SOURCE, GET_BOA_SOURCE private. All others now in PostProcessing.

      PRIVATE :: GET_TOA_SOURCE, GET_BOA_SOURCE

!  1/31/21. Version 2.8.3.  ==> Convergence routines no longer in this list

      PUBLIC :: VLIDORT_UPUSER_INTENSITY,  &
                VLIDORT_DNUSER_INTENSITY,  &
                VLIDORT_INTEGRATED_OUTPUT

      CONTAINS

      SUBROUTINE VLIDORT_UPUSER_INTENSITY ( & 
        DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,        & ! Input flags (RT mode)
        DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,       & ! Input flags (sources)
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   DO_INCLUDE_SURFEMISS,       & ! Input flags (Surface)
        DO_DBCORRECTION,     DO_INCLUDE_DIRECTRF,     DO_INCLUDE_DIRECTSL,       & ! Input flags (Surface)
        DO_LAYER_SCATTERING, DO_MSMODE_VLIDORT,       DO_MSMODE_THERMAL,         & ! Input flags (RT mode)
        DO_TOA_CONTRIBS, DO_INCLUDE_BOAFLUX, DO_CLASSICAL_SOLUTION, DO_MSSTS,    & ! Input flags (RT mode)
        FOURIER, IBEAM, ECHT_NSTOKES, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, & ! Input numbers (basic)
        N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, MUELLER_INDEX, TAYLOR_ORDER,   & ! Input bookkeeping + levels
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,            & ! Input partial-layer control
        FLUX_MULTIPLIER, BOAFLUX, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,    & ! Input Flux and quadrature
        T_DELT_DISORDS,  T_DELT_USERM, T_UTUP_USERM, CUMTRANS,                   & ! Input Transmittances
        ALBEDO, BRDF_F, USER_BRDF_F, SURFBB, EMISSIVITY, USER_EMISSIVITY,        & ! Input Surface BRDF/Emiss.
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                   & ! Input Homog. RTE Soln.
        WLOWER, LCON, MCON, T_WLOWER, LAYER_TSUP_UP, LAYER_TSUP_UTUP,            & ! Input RTE PI and thermal
        DELTAU_VERT, PARTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,      & ! Input Beam for Greens 
        T_UTDN_MUBAR, ITRANS_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,       & ! Input Greens function
        USER_DIRECT_BEAM, SL_USERTERM, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2,  & ! Input User solutions
        HMULT_1, HMULT_2, EMULT_UP, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, SIGMA_P, & ! Input multipliers
        PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,                               & ! Output multipliers (Greens)
        STOKES_DOWNSURF, BOA_THTONLY_SOURCE, MS_CONTRIBS_F,                         & ! OUTPUT 1 (Auxiliary)
        STOKES_F, CUMSOURCE_UP, LAYER_MSSTS_F, SURF_MSSTS_F )                         ! OUTPUT 2 (Main)

!  1/31/21. Version 2.8.3. SEVERAL INNOVATIONS and CHANGES
!    -- additional Taylor routines used for Green's function post-processing. Input taylor Order TAYLORM
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F, SURF_MSSTS_F
!    -- DO_CLASSICAL_SOLUTION flag added to control solution method (Greens vs. Classical)
!    -- Whole/part Sourceterm routines completely rewritten (Green's function included)
!    -- Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation vs. lattice/doublet geometry choices
!    -- (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES) number of STOKES components

!  1/31/21. Version 2.8.3. Following additional arrays are Inputs and outputs for the Green's function
!    -- DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, ITRANS_USERM, & ! Input Greens function stuff
!    -- ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, SIGMA_P,                   & ! Input Greens function stuff
!    -- PMULT_UU, PMULT_UD,  UT_PMULT_UU, UT_PMULT_UD,                       & ! Output Greens function multipliers

!   Streamlined for Version 2.8. 7/6/16
!   Version 2.8.1,  3/23/2019. Introduce Control for including BOA illumination

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS,             &
                                 MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_DIRECTIONS, &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, UPIDX

      IMPLICIT NONE

!  INPUTS
!  ======

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_MVOUTPUT

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFEMISS

      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTRF
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTSL

      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           DO_MSMODE_THERMAL

!  1/31/21. Version 2.8.3. Add DO_CLASSICAL_SOLUTION flag

      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_BOAFLUX ! New 3/23/19
      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::           DO_MSSTS

!  numbers
!    -- 1/31/21. Version 2.8.3. (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES)

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           ECHT_NSTOKES, NSTOKES, NLAYERS
      INTEGER, INTENT (IN) ::           NSTREAMS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Level output control + MUELLER_INDEX bookkeeping

      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           LEVELMASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX       ( MAXSTOKES, MAXSTOKES )

!  1/31/21. Version 2.8.3. Input Taylor-series parameter TAYLOR_ORDER

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER        ! New 2.8.3

!  Partial-layer output control

      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  BOA Flux (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) ::  BOAFLUX

!  Flux multipliers and quadrature

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )

!  User-stream and discrete ordinatetransmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM   ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  CUMTRANS       ( MAXLAYERS, MAX_USER_STREAMS )

!  Surface reflectance
!  1/31/21. Version 2.8.3.  ==> BRDF and SLEAVE arrays are defined locally for each Fourier

      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F       ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F  ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT (IN) ::  USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, intent (IN) ::  SL_USERTERM      ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  emissivity

      DOUBLE PRECISION, INTENT (IN) ::  EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SURFBB

!  RTE solutions (Homog.)

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

!  RTE PI solutions

      DOUBLE PRECISION, INTENT (IN) ::  WLOWER  ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON    ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON    ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UP   ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  User stream RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  EMULT_UP    ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Additional Variables for Green's function Post-processing (NEW)
!  ---------------------------------------------------------------------------------

!  Input optical depths

      DOUBLE PRECISION, INTENT (IN)   :: DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN)   :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  User/solar direction sigma factor

      DOUBLE PRECISION, INTENT (IN)   ::  SIGMA_P  ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Initial/Layer transmittance factors for solar beams, and divided by user-cosines

      INTEGER         , INTENT (IN)   :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN)   :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  OUTPUT
!  ======

!  1/31/21. Version 2.8.3. Post-processed Greens function multipliers (whole layer)

      DOUBLE PRECISION, INTENT (INOUT) :: PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (INOUT) :: PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
  
!  1/31/21. Version 2.8.3. Post-processed Greens function multipliers (Partial layer)

      DOUBLE PRECISION, INTENT (INOUT) :: UT_PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (INOUT) :: UT_PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  BOA source material

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_DOWNSURF    ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  Fourier contributions to radiance

      DOUBLE PRECISION, INTENT (INOUT) :: CUMSOURCE_UP ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_F     ( MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, 2 )
      DOUBLE PRECISION, INTENT (INOUT) :: MS_CONTRIBS_F  ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )

!  1/31/21. Version 2.8.3.   ==> Additional layer_mssts and surf_mssts output

      DOUBLE PRECISION, INTENT (INOUT) :: LAYER_MSSTS_F  ( MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION, INTENT (INOUT) :: SURF_MSSTS_F   ( MAX_SZANGLES, MAXSTOKES  )

!  local variables
!  ---------------

!  Help variables

      LOGICAL ::          SFLAG
      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NC, O1
      INTEGER ::          UT, UTA, UM, LUM
      DOUBLE PRECISION :: FINAL_SOURCE

!  Local sources

      DOUBLE PRECISION :: BOA_DIFFUSE_SOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: BOA_DIRECT_SOURCE  ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LAYER_SOURCE       ( MAX_USER_STREAMS, MAXSTOKES )

!  START OF CODE
!  =============

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).
!    -- 1/31/21. Version 2.8.3.  (RTS 2/16/21). Zero with the true NSTOKES value (ECHT_NSTOKES)
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present); rearranged order of loops

      !IF ( DO_USER_STREAMS ) THEN
      !  DO LUM = 1, N_PPSTREAMS
      !    UM = PPSTREAM_MASK(LUM,IBEAM)
      !    DO UTA = 1, N_USER_LEVELS
      !      DO O1 = 1, ECHT_NSTOKES
      !        STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = ZERO
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ENDIF

      IF ( DO_USER_STREAMS ) THEN
        DO O1 = 1, ECHT_NSTOKES
          DO LUM = 1, N_PPSTREAMS
            DO UTA = 1, N_USER_LEVELS
              STOKES_F(UTA,LUM,IBEAM,O1,UPIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. ==> Zero MSST outputs
!rob & mick fix 1/5/2021 - removed IF condition and changed limit from NSTOKES to ECHT_NSTOKES

      !IF ( .not. DO_MSSTS ) then
         DO O1 = 1, ECHT_NSTOKES
            LAYER_MSSTS_F(IBEAM,O1,1:NLAYERS) = ZERO
            SURF_MSSTS_F (IBEAM,O1)           = ZERO
         ENDDO
      !ENDIF

!  Initialize post-processing recursion
!  ====================================

!mick fix 3/31/2015 - moved GET_BOA_SOURCE outside of and before
!  DO_USER_STREAMS if block.  Get the BOA source terms: (1) quad
!  output and (2) user-angle output (diffuse + direct)
      
!   Version 2.8.1,  3/23/2019. Introduce Control for including BOA illumination

!  1/31/21. Version 2.8.3. Use N_PPSTREAMS, PPSTREAM_MASK masking system. Drop DO_OBSERVATION_GEOMETRY

      CALL GET_BOA_SOURCE ( &
        DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_INCLUDE_BOAFLUX,                          & ! Input flags
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS, DO_DBCORRECTION,  & ! Input flags
        DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, & ! Input flags
        FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK,            & ! Input Numbers
        QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX,                                         & ! Input bookkeeping
        BOAFLUX, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! Input surface stuff
        USER_DIRECT_BEAM, SL_USERTERM, SURFBB, EMISSIVITY, USER_EMISSIVITY, & ! Input Direct beam, emiss.
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,              & ! Input Homog solutions
        WLOWER, LCON, MCON, T_DELT_DISORDS, T_WLOWER,                       & ! Input Thermal and PI
        STOKES_DOWNSURF,   BOA_THTONLY_SOURCE,                              & ! BOA outputs
        BOA_DIFFUSE_SOURCE, BOA_DIRECT_SOURCE )                               ! BOA source outputs

!  1/31/21. Version 2.8.3. ==> Set the SURFACE MSST outputs if flagged
!    --  Important Note, Only for observational geometry

      IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY) THEN
         DO O1 = 1, NSTOKES
            SURF_MSSTS_F(IBEAM,O1) = FLUX_MULTIPLIER * BOA_DIFFUSE_SOURCE(IBEAM,O1)
         ENDDO
      ENDIF

!  Set the cumulative source term equal to BOA values
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

      IF ( DO_USER_STREAMS ) THEN
         NC = 0
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            DO O1 = 1, NSTOKES
               CUMSOURCE_UP(UM,O1,NC) = BOA_DIFFUSE_SOURCE(UM,O1) + BOA_DIRECT_SOURCE(UM,O1)
            ENDDO
         ENDDO
      ENDIF

!  Recursion Loop in Source function integration
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART   = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = LEVELMASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            SFLAG = DO_LAYER_SCATTERING(FOURIER,N)
            NC = NLAYERS + 1 - N

!  1/31/21. Version 2.8.3. 
!    -- New routine does classical or Greens, add flag DO_CLASSICAL_SOLUTION, Add Taylor ORder TAYLORM
!    -- New set of arguments needed for Green's function, including SIGMA_P, EMULT_UP
!    -- New outputs PMULT_UU, PMULT_UD from Greens calculation. Argument list reordered.
!    -- Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK). Drop DO_OBSERVATION_GEOMETRY

            CALL WHOLELAYER_STERM_UP ( &
               DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,         & ! Input flags
               DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SFLAG,                       & ! Input flags
               N, FOURIER, IBEAM, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,  & ! Input numbers
               DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,                 & ! Input Beam stuff for Greens
               ITRANS_USERM, T_DELT_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,  & ! Input Green's solution
               K_REAL, K_COMPLEX, LCON, MCON, UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2, & ! Input Homog solutions
               LAYER_TSUP_UP, UPAR_UP_1, UPAR_UP_2, SIGMA_P, EMULT_UP,                & ! Input Partic Integrals
               LAYER_SOURCE, PMULT_UU, PMULT_UD )                                       ! output 

!  1/31/21. Version 2.8.3. ==> For observational geometry, set  the MSST source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY) THEN
              DO O1 = 1, NSTOKES
                LAYER_MSSTS_F(IBEAM,O1,N) = FLUX_MULTIPLIER * LAYER_SOURCE(IBEAM,O1)
              ENDDO
            ENDIF

!  cumulative source
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IBEAM)
               DO O1 = 1, NSTOKES
                  IF ( DO_TOA_CONTRIBS ) THEN
                     MS_CONTRIBS_F(UM,IBEAM,O1,N) = FLUX_MULTIPLIER * CUMTRANS(N,UM) * LAYER_SOURCE(UM,O1)
                  ENDIF
                  !CALL TP7A (N,NC,UM,O1,LAYER_SOURCE,T_DELT_USERM,CUMSOURCE_UP)
                  CUMSOURCE_UP(UM,O1,NC) = LAYER_SOURCE(UM,O1) + T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,O1,NC-1)
               ENDDO
            ENDDO

!  End layer loop, and post processing clause

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

!  1/31/21. Version 2.8.3. 
!    -- New routine does classical or Greens, add flag DO_CLASSICAL_SOLUTION, Add Taylor ORder TAYLORM
!    -- New set of arguments needed for Green's function, including SIGMA_P, EMULT_UP
!    -- New outputs PMULT_UU, PMULT_UD from Greens calculation. Argument list reordered.
!    -- Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK). Drop OBSERVATION_GEOMETRY flag.

            CALL PARTLAYER_STERM_UP ( &
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
              DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SFLAG,                 & ! Input flags
              N, UT, IBEAM, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER, & ! Input numbers
              DELTAU_VERT, PARTAU_VERT, T_UTUP_USERM, ITRANS_USERM,            & ! Input delt
              BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,          & ! Input Beam stuff for Greens
              ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, K_REAL, K_COMPLEX,     & ! input Greens function
              LCON, MCON, UHOM_UPDN, UHOM_UPUP, UT_HMULT_UU, UT_HMULT_UD,      & ! Input Homog solutions
              LAYER_TSUP_UTUP, UPAR_UP_1, UPAR_UP_2, UT_EMULT_UP, SIGMA_P,     & ! Input Partic Integrals
              LAYER_SOURCE, UT_PMULT_UU, UT_PMULT_UD )                           ! output 

!  Final Partial layer source
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IBEAM)
               DO O1 = 1, NSTOKES
                  FINAL_SOURCE = LAYER_SOURCE(UM,O1) + T_UTUP_USERM(UT,UM)*CUMSOURCE_UP(UM,O1,NC)
                  STOKES_F(UTA,LUM,IBEAM,O1,UPIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
                  !CALL TP7B1 (UTA,UM,IBEAM,UT,NC,O1,FLUX_MULTIPLIER,&
                  !            LAYER_SOURCE,T_UTUP_USERM,CUMSOURCE_UP,STOKES_F)
               ENDDO
            ENDDO


          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term
!   Commented out code helps with FD testing
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IBEAM)
                DO O1 = 1, NSTOKES
                  FINAL_SOURCE = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,O1,NC)
                  STOKES_F(UTA,LUM,IBEAM,O1,UPIDX) = FINAL_SOURCE
!                  if (DABS(FINAL_SOURCE).GT.1.0d-10 ) then           ! FD
!                    STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = FINAL_SOURCE   ! FD
!                  endif                                              ! FD
                  !CALL TP7B2 (UTA,UM,IBEAM,NC,O1,FLUX_MULTIPLIER,CUMSOURCE_UP,STOKES_F)
                ENDDO
             ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_UPUSER_INTENSITY

!

      SUBROUTINE VLIDORT_DNUSER_INTENSITY ( &
        DO_USER_STREAMS,     DO_OBSERVATION_GEOMETRY, DO_INCLUDE_TOAFLUX,        & ! Input flags (RT mode)
        DO_SOLAR_SOURCES,    DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,      & ! Input flags (sources)
        DO_LAYER_SCATTERING, DO_MSMODE_VLIDORT, DO_CLASSICAL_SOLUTION, DO_MSSTS, & ! Input flags (RT mode)
        FOURIER, IBEAM, ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS,           & ! Input numbers (basic)
        N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_DN, TAYLOR_ORDER,                  & ! Input bookkeeping + levels
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,            & ! Input partial-layer control
        FLUX_MULTIPLIER, TOAFLUX, T_DELT_USERM, T_UTDN_USERM,                    & ! Input Transmittances, Flux
        K_REAL, K_COMPLEX, LCON, MCON, LAYER_TSUP_DN, LAYER_TSUP_UTDN,           & ! Input RTE Sol + thermal
        DELTAU_VERT, PARTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,      & ! Input Beam for Greens
        T_UTDN_MUBAR, ITRANS_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,    & ! Input Greens function 2   
        UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2, SIGMA_M,                     & ! Input User solutions
        HMULT_1, HMULT_2, EMULT_DN,  UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN,      & ! Input multipliers
        PMULT_DD, PMULT_DU, UT_PMULT_DU, UT_PMULT_DD,                            & ! Output mutlipliers (Greens)
        STOKES_F, CUMSOURCE_DN, LAYER_MSSTS_F  )                                   ! Main output

!  1/31/21. Version 2.8.3. SEVERAL INNOVATIONS and CHANGES
!    -- additional Taylor routines used for Green's function post-processing. Input taylor Order TAYLORM
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- DO_CLASSICAL_SOLUTION flag added to control solution method (Greens vs. Classical)
!    -- Whole/part Sourceterm routines completely rewritten (Green's function included)
!    -- Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation vs. lattice/doublet geometry choices
!    -- (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES) number of STOKES components

!  1/31/21. Version 2.8.3. Following additional arrays are Inputs and outputs for the Green's function
!    -- DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, ITRANS_USERM, & ! Input Greens function stuff
!    -- ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, SIGMA_M,                   & ! Input Greens function stuff
!    -- PMULT_DU, PMULT_DD,  UT_PMULT_DU, UT_PMULT_DD,                       & ! Output Greens function multipliers

!   Streamlined for Version 2.8. 7/6/16
        
!   Version 2.8.1,  3/23/2019. Introduce Control for including TOA illumination

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS,              &
                                 MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_DIRECTIONS, &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, DNIDX

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_TOAFLUX ! New 3/23/19

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT

!  1/31/21. Version 2.8.3. Add DO_CLASSICAL_SOLUTION flag

      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::           DO_MSSTS

!  numbers
!    -- 1/31/21. Version 2.8.3. NLAYERS added as input for MSSTS option
!    -- 1/31/21. Version 2.8.3. (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES)

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           ECHT_NSTOKES, NSTOKES, NLAYERS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Input Taylor-series parameter TAYLOR_ORDER

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER        ! New 2.8.3

!  Level output control

      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           LEVELMASK_DN  ( MAX_USER_LEVELS )

!  Partial-layer output control

      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  TOA Flux (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) ::  TOAFLUX
      
!  Flux multipliers and User-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM   ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  RTE solutions

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON    ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON    ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_DN   ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  User stream RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  EMULT_DN    ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Additional Variables for Green's function Post-processing
!  ---------------------------------------------------------------------------------

!  Input optical depths

      DOUBLE PRECISION, INTENT (IN)   :: DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN)   :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  User/solar direction sigma factor

      DOUBLE PRECISION, INTENT (IN)   :: SIGMA_M  ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Initial/Layer transmittance factors for solar beams, and divided by user-cosines

      INTEGER         , INTENT (IN)   :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN)   :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  OUTPUT
!  ======

!  1/31/21. Version 2.8.3. Post-processed Greens function multipliers (whole layer)

      DOUBLE PRECISION, INTENT (INOUT) :: PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (INOUT) :: PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
  
!  1/31/21. Version 2.8.3. Post-processed Greens function multipliers (Partial layer)

      DOUBLE PRECISION, INTENT (INOUT) :: UT_PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (INOUT) :: UT_PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Fourier contributions to radiance

      DOUBLE PRECISION, INTENT (INOUT) :: CUMSOURCE_DN ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_F     ( MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, 2 )

!  1/31/21. Version 2.8.3.   ==> Additional layer_mssts output

      DOUBLE PRECISION, INTENT (INOUT) :: LAYER_MSSTS_F  ( MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )

!  local variables
!  ---------------

      LOGICAL ::          SFLAG
      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER ::          UT, UTA, UM, NC, O1, LUM
      DOUBLE PRECISION :: TOA_SOURCE   ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: FINAL_SOURCE

!  Zero all Fourier components close to zenith
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).
!    -- 1/31/21. Version 2.8.3.  (RTS 2/16/21). Zero with the true NSTOKES value (ECHT_NSTOKES)
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present); rearranged order of loops

      !IF ( DO_USER_STREAMS ) THEN
      !  DO LUM = 1, N_PPSTREAMS
      !    UM = PPSTREAM_MASK(LUM,IBEAM)
      !    DO UTA = 1, N_USER_LEVELS
      !      DO O1 = 1, ECHT_NSTOKES
      !        STOKES_F(UTA,UM,IBEAM,O1,DNIDX) = ZERO
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ENDIF

      IF ( DO_USER_STREAMS ) THEN
        DO O1 = 1, ECHT_NSTOKES
          DO LUM = 1, N_PPSTREAMS
            DO UTA = 1, N_USER_LEVELS
              STOKES_F(UTA,LUM,IBEAM,O1,DNIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  1/31/21. Version 2.8.3.  ==> Zero MSST outputs
!rob & mick fix 1/5/2021 - removed IF condition and changed limit from NSTOKES to ECHT_NSTOKES

      !IF ( .not. DO_MSSTS ) then
         DO O1 = 1, ECHT_NSTOKES
            LAYER_MSSTS_F(IBEAM,O1,1:NLAYERS) = ZERO
         ENDDO
      !ENDIF

!  Initialize recursion for user-defined stream angles only
!    -- Version 2.8.1,  3/23/2019. Introduce Control for including TOA illumination
!    -- 1/31/21. Version 2.8.3. Use N_PPSTREAMS and PPSTREAM_MASK. Drop OBSERVATION_GEOMETRY flag.

      IF ( DO_USER_STREAMS ) THEN
        CALL GET_TOA_SOURCE (  &
          DO_INCLUDE_TOAFLUX, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, IBEAM, TOAFLUX, TOA_SOURCE )
        NC = 0
        DO LUM = 1, N_PPSTREAMS
           UM = PPSTREAM_MASK(LUM,IBEAM)
           DO O1 = 1, NSTOKES
              CUMSOURCE_DN(UM,O1,NC) = TOA_SOURCE(UM,O1)
           ENDDO
        ENDDO
      ENDIF

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = LEVELMASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            SFLAG = DO_LAYER_SCATTERING(FOURIER,N)
            NC = N

!  1/31/21. Version 2.8.3. 
!    -- New routine does classical or Greens, add flag DO_CLASSICAL_SOLUTION, Add Taylor ORder TAYLORM
!    -- New set of arguments needed for Green's function, including SIGMA_M, EMULT_DN
!    -- New outputs PMULT_DD, PMULT_DU from Greens calculation. Argument list reordered.
!    -- Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK). Drop DO_OBSERVATION_GEOMETRY

            CALL WHOLELAYER_STERM_DN ( &
               DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,         & ! Input flags
               DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SFLAG,                       & ! Input flags
               N, FOURIER, IBEAM, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,  & ! Input numbers
               DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,                 & ! Input Beam stuff for Greens
               ITRANS_USERM, T_DELT_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,  & ! Input Green's solution
               K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, HMULT_1, HMULT_2, & ! Input Homog solutions
               LAYER_TSUP_DN, UPAR_DN_1, UPAR_DN_2, SIGMA_M, EMULT_DN,                & ! Input Partic Integrals
               LAYER_SOURCE, PMULT_DD, PMULT_DU )                                       ! output 

!  Rob debug 3/1/21
!            if ( fourier.eq.2) then
!              DO LUM = 1, N_PPSTREAMS
!                 UM = PPSTREAM_MASK(LUM,IBEAM)
!                 write(400,*)N,NC,UM,LAYER_SOURCE(UM,1),CUMSOURCE_DN(UM,1,NC-1)
!               ENDDO
!            ENDIF

!  1/31/21. Version 2.8.3. ==> For observational geometry, set  the MSST source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY) THEN
              DO O1 = 1, NSTOKES
                LAYER_MSSTS_F(IBEAM,O1,N) = FLUX_MULTIPLIER * LAYER_SOURCE(IBEAM,O1)
              ENDDO
            ENDIF

!  Cumulative sources
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IBEAM)
               DO O1 = 1, NSTOKES
                  CUMSOURCE_DN(UM,O1,NC) = LAYER_SOURCE(UM,O1) + T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,O1,NC-1)
               ENDDO
            ENDDO

!  End layer recursion loop

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

!  1/31/21. Version 2.8.3. 
!    -- New routine does classical or Greens, add flag DO_CLASSICAL_SOLUTION, Add Taylor ORder TAYLORM
!    -- New set of arguments needed for Green's function, including SIGMA_P, EMULT_UP
!    -- New outputs PMULT_UU, PMULT_UD from Greens calculation. Argument list reordered.
!    -- Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK). Drop OBSERVATION_GEOMETRY

            CALL PARTLAYER_STERM_DN ( &
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
              DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SFLAG,                 & ! Input flags
              N, UT, IBEAM, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER, & ! Input numbers
              PARTAU_VERT, T_UTDN_USERM, ITRANS_USERM,                         & ! Input delt/trans
              BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,          & ! Input Beam stuff
              ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, K_REAL, K_COMPLEX,     & ! input Green s function
              LCON, MCON, UHOM_DNDN, UHOM_DNUP, UT_HMULT_DD, UT_HMULT_DU,      & ! Input Homog solutions
              LAYER_TSUP_UTDN, UPAR_DN_1, UPAR_DN_2, UT_EMULT_DN, SIGMA_M,     & ! Input Partic Integrals
              LAYER_SOURCE, UT_PMULT_DU, UT_PMULT_DD )                           ! output 

!  Rob debug 3/1/21. Helped to fix partial-layer Taylor error
!           if ( fourier.lt.6) write(*,*)'gronk partials',Fourier,LAYER_SOURCE(1:4,1)

!  Final Partial layer source
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IBEAM)
               DO O1 = 1, NSTOKES
                  FINAL_SOURCE = LAYER_SOURCE(UM,O1) + T_UTDN_USERM(UT,UM)*CUMSOURCE_DN(UM,O1,NC)
                  STOKES_F(UTA,LUM,IBEAM,O1,DNIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
                  !CALL TP7D1 (UTA,UM,IBEAM,UT,NC,FLUX_MULTIPLIER,&
                  !            LAYER_SOURCE,T_UTDN_USERM,CUMSOURCE_DN,STOKES_F)
               ENDDO
            ENDDO

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IBEAM)
                DO O1 = 1, NSTOKES
                   STOKES_F(UTA,LUM,IBEAM,O1,DNIDX) = FLUX_MULTIPLIER * CUMSOURCE_DN(UM,O1,NC)
                ENDDO
             ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_DNUSER_INTENSITY

!

      SUBROUTINE GET_TOA_SOURCE ( &
          DO_INCLUDE_TOAFLUX, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, IBEAM, TOAFLUX, TOA_SOURCE )

!    Version 2.8.1,  3/23/2019. Introduce Control for including TOA illumination

!    -- 1/31/21. Version 2.8.3. Use N_PPSTREAMS and PPSTREAM_MASK. Drop OBSERVATION_GEOMETRY flag.

      USE VLIDORT_PARS_m, only : MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, ZERO

      IMPLICIT NONE

!  Input

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_TOAFLUX ! New 3/23/19
      INTEGER, INTENT (IN) ::          NSTOKES, IBEAM

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  TOA Flux (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) ::  TOAFLUX
      
!  output

      DOUBLE PRECISION, INTENT (OUT) :: TOA_SOURCE ( MAX_USER_STREAMS, MAXSTOKES )

!  local variables

      INTEGER :: UM, O1, LUM

!  initialise TOA source function
!    - Also include TOA illumination if flagged
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IBEAM)
         DO O1 = 1, NSTOKES
            TOA_SOURCE(UM,O1) = ZERO
         ENDDO
         IF ( DO_INCLUDE_TOAFLUX ) TOA_SOURCE(UM,1) = TOA_SOURCE(UM,1) + TOAFLUX
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE GET_TOA_SOURCE

!

      SUBROUTINE GET_BOA_SOURCE (  &
        DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_INCLUDE_BOAFLUX,                           & ! Input flags
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS, DO_DBCORRECTION,   & ! Input flags
        DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,  & ! Input flags
        FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK,             & ! Input Numbers
        QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX,                                          & ! Input bookkeeping
        BOAFLUX, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! Input surface stuff
        USER_DIRECT_BEAM, SL_USERTERM, SURFBB, EMISSIVITY, USER_EMISSIVITY, & ! Input Direct beam, emiss.
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,              & ! Input Homog solutions
        WLOWER, LCON, MCON, T_DELT_DISORDS, T_WLOWER,                       & ! Input Thermal and PI
        STOKES_DOWNSURF,    BOA_THTONLY_SOURCE,                              & ! BOA outputs
        BOA_DIFFUSE_SOURCE, BOA_DIRECT_SOURCE )                               ! BOA source outputs

!  Bottom of the atmosphere source term, Lambertian or BRDF
!    Version 2.8.1,  3/23/2019. Introduce Control for including BOA illumination
!    Version 2.8.1,  4/29/2019. Separate control for DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL

!mick fix 3/31/2015 - added DO_USER_STREAMS flag

!  1/31/21. Version 2.8.3. BRDF and SLEAVE arrays are defined locally for each Fourier
!    -- Drop MAXMOMENTS from the parameter list.

!  1/31/21. Version 2.8.3. Use N_PPSTREAMS and PPSTREAM_MASK. Drop OBSERVATION_GEOMETRY flag.

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAXBEAMS, &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE

      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags
!    -- 1.31.21. Version 2.8.3. Drop OBSERVATION_GEOMETRY

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_MVOUTPUT

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTRF
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTSL

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_BOAFLUX ! New 3/23/19

      LOGICAL, INTENT (IN) ::           DO_MSMODE_THERMAL
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFEMISS

!  Numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Bookkeeping (quadrature etc)
!    -- 1/31/21. Version 2.8.3. Drop LOCAL_UM_START

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  BOA Flux (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) ::  BOAFLUX
      
!  Surface reflectance
!  1/31/21. Version 2.8.3. BRDF and SLEAVE arrays are defined locally for each Fourier

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F      ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      
      DOUBLE PRECISION, INTENT (IN) ::  USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, intent (IN) ::  SL_USERTERM      ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Emissivity

      DOUBLE PRECISION, INTENT (IN) ::  EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SURFBB

!  Homogeneous solutions

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

!  Particular integrals and thermal

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER   ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Output
!  ======

!  Quad-angle output

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_DOWNSURF    ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  User-angle output

      DOUBLE PRECISION, INTENT (INOUT) :: BOA_DIFFUSE_SOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: BOA_DIRECT_SOURCE  ( MAX_USER_STREAMS, MAXSTOKES )


!  local variables
!  ---------------

      LOGICAL          :: DO_QTHTONLY
      INTEGER ::          N, J, I, IR, IROW, UM, O1, O2, O11, OM, LUM
      INTEGER ::          K, KO1, K0, K1, K2, M, IB
      DOUBLE PRECISION :: REFLEC, S_REFLEC, STOTAL, KMULT
      DOUBLE PRECISION :: SPAR, HOM1, HOM2, SHOM_R, LOCAL_EMISS, FP_SBB
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR
      DOUBLE PRECISION :: THELP ( MAXSTREAMS )

!  Local index

      IB = IBEAM

!  initialise boa source function

!mick fix 3/31/2015 - added DO_USER_STREAMS if condition
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
               BOA_DIFFUSE_SOURCE(UM,O1) = ZERO
               BOA_DIRECT_SOURCE(UM,O1)  = ZERO
            ENDDO
         ENDDO
      ENDIF

!  Version 2.8.1. 3/23/19.  Also include BOA illumination if flagged
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         IF ( DO_INCLUDE_BOAFLUX ) BOA_DIFFUSE_SOURCE(UM,1) = BOA_DIFFUSE_SOURCE(UM,1) + BOAFLUX
      ENDDO

!  Special flag
!   Version 2.8, 7/5/16. Remove DO_QUAD_OUTPUT

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .AND. DO_INCLUDE_MVOUTPUT )
!      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .AND. ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) )

!  Layer index and offset

      N = NLAYERS

!mick thermal fix - added if condition
      !KO1 = K_REAL(N) + 1
      IF ( .NOT. DO_THERMAL_TRANSONLY ) KO1 = K_REAL(N) + 1

!  fourier component

      M = FOURIER

!  1-1 element index

      O11 = 1

!  reflectance from surface
!  ------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN

!mick chg 3/31/2015 - transformed the separate "DO_THERMAL_TRANSONLY" &
!  ".NOT.DO_THERMAL_TRANSONLY" if blocks into a composite IF/ELSEIF block

!  Thermal transmittance solution, build from TOA downwards

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            THELP(I) = ZERO
          ENDDO
          DO K = 1, NLAYERS
            DO I = 1, NSTREAMS
              THELP(I) = THELP(I)*T_DELT_DISORDS(I,K) + T_WLOWER(I,K)
            ENDDO
          ENDDO
          DO I = 1, NSTREAMS
            STOKES_DOWNSURF(I,O11) = QUAD_WEIGHTS(I) * THELP(I)
          ENDDO

!  Full solution: Downward intensity at computational angles (beam/homog
!     --> Develop reflectance integrand  a(j).x(j).I(-j)

        ELSE IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN

!  Downward Stokes vector at surface at computational angles
!    reflectance integrand  a(j).x(j).I(-j)

          DO I = 1, NSTREAMS
            IR = ( I - 1 ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              SPAR = WLOWER(I,O1,N)
              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
                MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
                HOM1 = LXR * T_DELT_EIGEN(K,N)
                HOM2 = MXR
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO
              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K-2
                K1 = KO1 + K0
                K2 = K1  + 1
                LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - &
                          LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + &
                          LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - &
                          MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                HOM1CR = LXR_CR*T_DELT_EIGEN(K1,N) &
                        -LXR_CI*T_DELT_EIGEN(K2,N)
                HOM2CR = MXR_CR
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
              ENDDO
              STOTAL = SPAR + SHOM_R + SHOM_CR
              STOKES_DOWNSURF(I,O1) =  QUAD_STRMWTS(I) * STOTAL
            ENDDO
          ENDDO

        ENDIF

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  ###### Lambertian reflectance (same for all user-streams)
!   @@@ Rob Fix, 2/9/11, Position of DO_QTHTONLY clause is wrong
!  Version 2.8, 7/5/16. Code streamlined

        IF ( DO_LAMBERTIAN_SURFACE ) THEN

          IF ( FOURIER .EQ. 0 ) THEN
            KMULT = SURFACE_FACTOR * ALBEDO
            O1 = 1
            REFLEC = SUM(STOKES_DOWNSURF(1:NSTREAMS,O1))
!            IF ( DO_QTHTONLY ) BOA_THTONLY_SOURCE(1:NSTREAMS,O1) = REFLEC
            REFLEC = KMULT * REFLEC

!mick fix 3/31/2015 - added DO_USER_STREAMS if block
!  1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

            IF ( DO_USER_STREAMS ) THEN            
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                BOA_DIFFUSE_SOURCE(UM,O1) = REFLEC
              ENDDO
            ENDIF
            IF ( DO_QTHTONLY ) BOA_THTONLY_SOURCE(1:NSTREAMS,O1) = REFLEC
          ENDIF

!  ###### bidirectional reflectance
!   @@@ Rob Fix, 2/9/11, DO_QTHTONLY clause is wrong, and wrong position

!  1/31/21. Version 2.8.3. BRDF arrays defined locally, drop M Fourier index
!  1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

        ELSE

          KMULT = SURFACE_FACTOR
          DO LUM = 1, N_PPSTREAMS
             UM = PPSTREAM_MASK(LUM,IB)
             DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                   S_REFLEC = ZERO
                   DO O2 = 1, NSTOKES
                      OM = MUELLER_INDEX(O1,O2)
                      S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * USER_BRDF_F(OM,UM,J)
                   ENDDO
                   REFLEC = REFLEC + S_REFLEC
!                 IF ( DO_QTHTONLY ) BOA_THTONLY_SOURCE(1:NSTREAMS,O1) = REFLEC
                ENDDO
                BOA_DIFFUSE_SOURCE(UM,O1) = KMULT * REFLEC
             ENDDO
          ENDDO

          IF ( DO_QTHTONLY ) THEN
            O11 = 1
            DO I = 1, NSTREAMS
              REFLEC = DOT_PRODUCT(STOKES_DOWNSURF(1:NSTREAMS,O11), BRDF_F(O11,I,1:NSTREAMS) )
              BOA_THTONLY_SOURCE(I,O11) = KMULT * REFLEC
            ENDDO
          ENDIF

        ENDIF

!  Add reflected direct beam if flagged
!   Addition of the DBCORRECTION clause is now required

!mick fix 3/31/2015 - added DO_USER_STREAMS to if condition
!  1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

        IF ( DO_USER_STREAMS .AND. DO_INCLUDE_DIRECTRF .AND. .NOT.DO_DBCORRECTION ) THEN
           DO LUM = 1, N_PPSTREAMS
              UM = PPSTREAM_MASK(LUM,IB)
              DO O1 = 1, NSTOKES
                 BOA_DIRECT_SOURCE(UM,O1) = BOA_DIRECT_SOURCE(UM,O1) + USER_DIRECT_BEAM(UM,IB,O1)
              ENDDO
           ENDDO
        ENDIF

!  Add direct surface-leaving if flagged. New. 2.8.1
!  1/31/21. Version 2.8.3.  Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

        IF ( DO_INCLUDE_DIRECTSL .AND. DO_USER_STREAMS ) THEN
           DO LUM = 1, N_PPSTREAMS
              UM = PPSTREAM_MASK(LUM,IB)
              DO O1 = 1, NSTOKES
                 BOA_DIRECT_SOURCE(UM,O1) = BOA_DIRECT_SOURCE(UM,O1) + SL_USERTERM(LUM,IB,O1)
              ENDDO
           ENDDO
        ENDIF

!  End inclusion of surface terms

      ENDIF

!  Add surface emission term if flagged
!  Version 2.8. 7/5/16. removed "goto 23"

      IF ( DO_INCLUDE_SURFEMISS ) THEN

         FP_SBB = SURFBB
         IF ( .not. DO_MSMODE_THERMAL ) then
!mick fix 3/31/2015 - added DO_USER_STREAMS if block
            IF ( DO_USER_STREAMS ) THEN
               IF ( DO_LAMBERTIAN_SURFACE ) THEN
                  O1 = 1
                  LOCAL_EMISS = FP_SBB * ( ONE - ALBEDO )
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     BOA_DIFFUSE_SOURCE(UM,O1) = BOA_DIFFUSE_SOURCE(UM,O1) + LOCAL_EMISS
                  ENDDO
               ELSE
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     DO O1 = 1, NSTOKES
!  @@@ Rob fix, ordering of USER_EMISSIVITY indices wrong
!                       BOA_DIFFUSE_SOURCE(UM,O1) = BOA_DIFFUSE_SOURCE(UM,O1) + FP_SBB*USER_EMISSIVITY(UM,O1)
                        BOA_DIFFUSE_SOURCE(UM,O1) = BOA_DIFFUSE_SOURCE(UM,O1) + FP_SBB*USER_EMISSIVITY(O1,UM)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDIF   

!  Thermal tranmsittance only

         IF ( DO_QTHTONLY ) THEN
            IF ( DO_LAMBERTIAN_SURFACE ) THEN
               O1 = 1
               LOCAL_EMISS = FP_SBB * ( ONE - ALBEDO )
               DO I = 1, NSTREAMS
                  BOA_THTONLY_SOURCE(I,O1) = BOA_THTONLY_SOURCE(I,O1) + LOCAL_EMISS
               ENDDO
            ELSE
               DO I = 1, NSTREAMS
                  DO O1 = 1, NSTOKES
!  @@@ Rob fix, ordering of EMISSIVITY indices wrong
!                    BOA_THTONLY_SOURCE(I,O1) = BOA_THTONLY_SOURCE(I,O1) + FP_SBB * EMISSIVITY(I,O1)
                     BOA_THTONLY_SOURCE(I,O1) = BOA_THTONLY_SOURCE(I,O1) + FP_SBB * EMISSIVITY(O1,I)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF

!  End surface emission clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE GET_BOA_SOURCE

!

      SUBROUTINE VLIDORT_INTEGRATED_OUTPUT ( &
        DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM,                             & ! Input flags
        DO_CLASSICAL_SOLUTION, DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX,          & ! Input flags
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,          & ! Input flags
        IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, N_DIRECTIONS,         & ! Input numbers
        NSTKS_NSTRMS, LEVELMASK_UP, LEVELMASK_DN, WHICH_DIRECTIONS,             & ! Input Bookkeeping
        PARTLAYERS_LAYERIDX, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,           & ! Input partial layer control
        TAYLOR_ORDER, FLUX_MULTIPLIER, TOAFLUX, BOAFLUX, FLUX_FACTOR,           & ! Flux inputs
        FLUXVEC, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, PARTAU_VERT,         & ! Quadrature inputs
        BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                 & ! Solar beam Param.
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, LOCAL_CSZA,                     & ! Solar beam transmittances
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                         & ! Discrete ordinate transmittances
        ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P,                               & ! Input Greens function
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, LCON, MCON,                    & ! Input RTE
        T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, BVEC, WUPPER, WLOWER,         & ! Input RTE
        T_WLOWER, T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,                    & ! Input Thermal
        MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT,           & ! MAIN output
        UT_GMULT_UP, UT_GMULT_DN )                                                ! Auxiliary Output

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
        
!  Version 2.8. 2017.
!mick mod 9/19/2017 - output flux variables renamed: distinguish between diffuse and direct
!mick fix 9/19/2017 - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to facilitate correction of direct flux
        
!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19

!  1/31/21. Version 2.8.3.
!   -- Add flag DO_CLASSICAL_SOLUTION, controls use of Greens function
!   -- Additional inputs BVEC and INITIAL_TRANS, T_UTDN_MUBAR, TAYLOR ORDER, NSTKS_NSTRMS
!   -- Additional inputs ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P foro the Green's function


      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXBEAMS, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_LEVELS, &
                                 MAX_SZANGLES, MAX_DIRECTIONS, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS,         &
                                 ZERO, HALF, ONE, PI2, PI4, UPIDX, DNIDX

      IMPLICIT NONE

!  input flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_TOAFLUX ! New 3/23/19
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_BOAFLUX ! New 3/23/19
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  3/16/20GF. Greens function, add flag DO_CLASSICAL_SOLUTION

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Input numbers

      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          TAYLOR_ORDER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Optical depths

      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  TOA/BOA Flux (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) :: TOAFLUX
      DOUBLE PRECISION, INTENT (IN) :: BOAFLUX

!  Flux control

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )

!  output levels

      INTEGER, INTENT (IN) ::          LEVELMASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          LEVELMASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS    ( MAX_DIRECTIONS )

!  Partial layer control

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  solar beam parameterization 

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA  ( 0:MAXLAYERS, MAXBEAMS )

!  3/16/20GF. Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  Derived Solar-beam Transmittance at all levels

      DOUBLE PRECISION, INTENT (IN) :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

!  Discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  RTE solutions

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Use BVEC for the classical solution (replaces use of WUPPER)

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER    ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  output
!  ------

!  1/31/21. Version 2.8.3. Local Green functions multipliers for off-grid optical depths

      DOUBLE PRECISION, INTENT (INOUT) :: UT_GMULT_UP ( MAXEVALUES,MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UT_GMULT_DN ( MAXEVALUES,MAX_PARTLAYERS )

!mick mod 9/19/2017 - renamed fluxes for separating diffuse and direct

!  Mean stokes (actinic flux), Regular Flux, Diffuse values

      DOUBLE PRECISION, INTENT (INOUT) :: MEANST_DIFFUSE ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_DIFFUSE   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

!  Direct-beam contributions output separately, 26 May 11

      DOUBLE PRECISION, INTENT (INOUT) :: DNMEANST_DIRECT ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: DNFLUX_DIRECT   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  local variables
!  ---------------

!  Local quadrature output

      DOUBLE PRECISION :: QSTOKES_F ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  1/31/21. Version 2.8.3. Local flag for partial-layer Green's function multiplier creation.

      LOGICAL          :: HAVE_MULT ( MAX_PARTLAYERS )

!  help variables

      INTEGER          :: I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1
      DOUBLE PRECISION :: SUM_MI, SUM_FX
      DOUBLE PRECISION :: DIRECT_TRANS, FTRANS, DNDIRECT_MEANST, DNDIRECT_FLUX
!mick fix 9/19/2017 - added to facilitate correction of direct flux
      DOUBLE PRECISION :: DIRECT_TRANS_SCALED, FTRANS_SCALED, DNDIRECT_MEANST_SCALED, DNDIRECT_FLUX_SCALED

!  Mean intensity and flux
!  -----------------------

!mick fix 1/5/2021 - initialized HAVE_MULT
      HAVE_MULT = .false.

!  Direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

!  Upwelling Stokes output at Quadrature angles

!  1/31/21. Version 2.8.3. QUADINTENS_OFFGRID_DN Call
!    -- Additional I/O arguments associated with the Green's function option
!    -- New Inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER, PARTAU_VERT
!    -- Use Input BVEC instead of WUPPER for the classical solution
!    -- New Inputs BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR (solar beam parameterization)
!    -- New Inputs ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P (Green's function solution)
!    -- New Outputs: Green's multipliers UT_GMULT_UP, UT_GMULT_DN and flag HAVE_MULT

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = LEVELMASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADINTENS_OFFGRID_UP ( &
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,            & ! Input flags
                DO_CLASSICAL_SOLUTION, DO_INCLUDE_BOAFLUX, N, UTA, UT, IBEAM,             & ! Input flags/indices
                NSTOKES, NSTREAMS, NSTKS_NSTRMS, NLAYERS, TAYLOR_ORDER, HAVE_MULT(UT),    & ! Input numbers/Bookkeeping
                FLUX_MULTIPLIER, BOAFLUX, PARTAU_VERT, QUAD_STREAMS,                      & ! Input quad/Fluxes/delt
                BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                   & ! Input Beam for Greens
                T_DELT_DISORDS, T_DISORDS_UTUP, ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, & ! input Disords/Greens
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! Input Homog. solutions
                BVEC, LCON, MCON, T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,              & ! Input thermal and PI solutions
                QSTOKES_F, UT_GMULT_UP, UT_GMULT_DN  )                                      ! OUTPUT
            ELSE
              CALL QUADINTENS_LEVEL_UP ( &
                DO_THERMAL_TRANSONLY, DO_INCLUDE_BOAFLUX,                  & ! Input Flags
                NLEVEL, UTA, NSTOKES, NSTREAMS, NLAYERS,                   & ! Input numbers 
                FLUX_MULTIPLIER, BOAFLUX, QUAD_STREAMS, T_DELT_DISORDS,    & ! Input Flux/Quad/Trans
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,     & ! Input Homog. solutions
                WLOWER, WUPPER, LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,  & ! Input thermal and PI solutions
                QSTOKES_F )                                                  ! OUTPUT
            ENDIF
          ENDDO
        ENDIF

!  Downwelling Stokes output at Quadrature angles
!   * Version 2.8a, Control for TOA isotropic illumination added, 3/23/19

!  1/31/21. Version 2.8.3. QUADINTENS_OFFGRID_DN Call
!    -- Additional I/O arguments associated with the Green's function option
!    -- New Inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER, PARTAU_VERT
!    -- Use Input BVEC instead of WUPPER for the classical solution
!    -- New Inputs BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR (solar beam parameterization)
!    -- New Inputs ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P (Green's function solution)
!    -- New Outputs: Green's multipliers UT_GMULT_UP, UT_GMULT_DN and flag HAVE_MULT

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = LEVELMASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADINTENS_OFFGRID_DN ( &
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,            & ! Input flags
                DO_CLASSICAL_SOLUTION, DO_INCLUDE_TOAFLUX, N, UTA, UT, IBEAM,             & ! Input flags/indices
                NSTOKES, NSTREAMS, NSTKS_NSTRMS, TAYLOR_ORDER, HAVE_MULT(UT),             & ! Input numbers/Bookkeeping
                FLUX_MULTIPLIER, TOAFLUX, PARTAU_VERT, QUAD_STREAMS,                      & ! Input quad/Fluxes/delt
                BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                   & ! Input Beam for Greens
                T_DELT_DISORDS, T_DISORDS_UTDN, ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, & ! Input Disords/Greens
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! Input Homog. solutions
                BVEC, LCON, MCON, T_WLOWER, UT_T_PARTIC,                                  & ! Input thermal and PI solutions
                QSTOKES_F, UT_GMULT_UP, UT_GMULT_DN  )                                      ! OUTPUT
            ELSE
              CALL QUADINTENS_LEVEL_DN ( &
                DO_THERMAL_TRANSONLY, DO_INCLUDE_TOAFLUX,              & ! Input flags
                NLEVEL, UTA, NSTOKES, NSTREAMS, FLUX_MULTIPLIER,       & ! Input numbers and Flux
                TOAFLUX, QUAD_STREAMS, T_DELT_DISORDS,                 & ! Input TOAFlux, stream transmittances
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, & ! Input Homog. solutions
                WLOWER, LCON, MCON, T_WLOWER,                          & ! Input thermal and PI solutions
                QSTOKES_F )                                              ! OUTPUT
            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux output
!  ------------------------------

        IF ( DO_INCLUDE_MVOUTPUT ) THEN

!  Diffuse term output

          DO UTA = 1, N_USER_LEVELS
            DO O1 = 1, NSTOKES
              SUM_MI = ZERO
              SUM_FX = ZERO
              DO I = 1, NSTREAMS
                SUM_MI = SUM_MI + QUAD_WEIGHTS(I) * QSTOKES_F(UTA,I,O1)
                SUM_FX = SUM_FX + QUAD_STRMWTS(I) * QSTOKES_F(UTA,I,O1)
              ENDDO
              MEANST_DIFFUSE(UTA,IBEAM,O1,WDIR) = SUM_MI * HALF
              FLUX_DIFFUSE(UTA,IBEAM,O1,WDIR)   = SUM_FX * PI2
            ENDDO
          ENDDO

!  Nothing to do if no solar sources

!  Version 2.8. remove GOTO 455
!          IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) GO TO 455

!  For the downward direction, add the direct beam contributions
!mick fix 9/19/2017 - use unscaled optical thicknesses for calculation of direct flux
!                     (following LIDORT)

          IF ( DO_INCLUDE_DIRECTBEAM .and. (WDIR .EQ. DNIDX) ) THEN

!  Loop over all the output optical depths

            DO UTA = 1, N_USER_LEVELS

!  For the offgrid values.......
!     .....Only contributions for layers above the PI cutoff

              IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
                UT = PARTLAYERS_OUTINDEX(UTA)
                N  = PARTLAYERS_LAYERIDX(UT)

!                IF ( N .LE. BEAM_CUTOFF(IBEAM) ) THEN
!                  DIRECT_TRANS = INITIAL_TRANS(N,IBEAM) * T_UTDN_MUBAR(UT,IBEAM)
!                  DO O1 = 1, NSTOKES
!
!!  @@@ Rob fix (forgot the Flux factor)
!!                   FTRANS = FLUXVEC(O1) * DIRECT_TRANS
!!                   FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS
!!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!!                   FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * DIRECT_TRANS
!! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
!
!                    FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS
!                    DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
!                    DIRECT_MEANI = FTRANS / PI4
!                    MEAN_DIRECT(UTA,IBEAM,O1) = DIRECT_MEANI
!                    FLUX_DIRECT(UTA,IBEAM,O1) = DIRECT_FLUX
!                    MEAN_STOKES(UTA,IBEAM,O1,WDIR) = MEAN_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_MEANI
!                    FLUX_STOKES(UTA,IBEAM,O1,WDIR) = FLUX_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_FLUX
!                  ENDDO
!                ENDIF

                IF ( N .LE. BEAM_CUTOFF(IBEAM) ) THEN

!  Direct transmittances, scaled and unscaled

                  DIRECT_TRANS        = PARTIALS_SOLARTRANS(UT,IBEAM)
                  DIRECT_TRANS_SCALED = INITIAL_TRANS(N,IBEAM) * T_UTDN_MUBAR(UT,IBEAM)

                  DO O1 = 1, NSTOKES

!  Direct calculation with non-scaled transmittances

                    FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS
                    DNDIRECT_MEANST = FTRANS / PI4
                    DNDIRECT_FLUX   = FTRANS * LOCAL_CSZA(N,IBEAM)
                    DNMEANST_DIRECT(UTA,IBEAM,O1) = DNDIRECT_MEANST
                    DNFLUX_DIRECT(UTA,IBEAM,O1)   = DNDIRECT_FLUX

!  Diffuse calculation

                    FTRANS_SCALED = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS_SCALED
                    DNDIRECT_MEANST_SCALED = FTRANS_SCALED / PI4
                    DNDIRECT_FLUX_SCALED   = FTRANS_SCALED * LOCAL_CSZA(N,IBEAM)
                    MEANST_DIFFUSE(UTA,IBEAM,O1,WDIR) = &
                          MEANST_DIFFUSE(UTA,IBEAM,O1,WDIR) + ( DNDIRECT_MEANST_SCALED - DNDIRECT_MEANST )
                    FLUX_DIFFUSE(UTA,IBEAM,O1,WDIR)   = &
                          FLUX_DIFFUSE(UTA,IBEAM,O1,WDIR)   + ( DNDIRECT_FLUX_SCALED   - DNDIRECT_FLUX )
                  ENDDO
                ENDIF

!  For the on-grid values

              ELSE
                N = LEVELMASK_DN(UTA)

!                IF ( N .LE. BEAM_CUTOFF(IBEAM) ) THEN
!                  IF ( N .EQ. 0 ) THEN
!                    DIRECT_TRANS = ONE
!                  ELSE
!                    DIRECT_TRANS = INITIAL_TRANS(N,IBEAM)*T_DELT_MUBAR(N,IBEAM)
!                  ENDIF
!                  DO O1 = 1, NSTOKES
!
!!  @@@ Rob fix (forgot the Flux factor)
!!                 FTRANS = FLUXVEC(O1) * DIRECT_TRANS
!!                 FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS
!!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!!                  FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * DIRECT_TRANS
!! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
!
!                    FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS
!                    DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
!                    DIRECT_MEANI = FTRANS / PI4
!                    MEAN_DIRECT(UTA,IBEAM,O1) = DIRECT_MEANI
!                    FLUX_DIRECT(UTA,IBEAM,O1) = DIRECT_FLUX
!                    MEAN_STOKES(UTA,IBEAM,O1,WDIR) = MEAN_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_MEANI
!                    FLUX_STOKES(UTA,IBEAM,O1,WDIR) = FLUX_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_FLUX
!                  ENDDO
!                ENDIF

                IF ( N .LE. BEAM_CUTOFF(IBEAM) ) THEN

!  Direct transmittances, scaled and unscaled

                  IF ( N .EQ. 0 ) THEN
                    DIRECT_TRANS = ONE
                    DIRECT_TRANS_SCALED = ONE
                  ELSE
                    DIRECT_TRANS = LEVELS_SOLARTRANS(N,IBEAM)
                    DIRECT_TRANS_SCALED = INITIAL_TRANS(N,IBEAM)*T_DELT_MUBAR(N,IBEAM)
                  ENDIF

                  DO O1 = 1, NSTOKES

!  Direct calculation with non-scaled transmittances

                    FTRANS = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS
                    DNDIRECT_MEANST = FTRANS / PI4
                    DNDIRECT_FLUX   = FTRANS * LOCAL_CSZA(N,IBEAM)
                    DNMEANST_DIRECT(UTA,IBEAM,O1) = DNDIRECT_MEANST
                    DNFLUX_DIRECT(UTA,IBEAM,O1)   = DNDIRECT_FLUX

!  Diffuse calculation

                    FTRANS_SCALED = FLUX_FACTOR * FLUXVEC(O1) * DIRECT_TRANS_SCALED
                    DNDIRECT_MEANST_SCALED = FTRANS_SCALED / PI4
                    DNDIRECT_FLUX_SCALED   = FTRANS_SCALED * LOCAL_CSZA(N,IBEAM)
                    MEANST_DIFFUSE(UTA,IBEAM,O1,WDIR) = &
                          MEANST_DIFFUSE(UTA,IBEAM,O1,WDIR) + ( DNDIRECT_MEANST_SCALED - DNDIRECT_MEANST )
                    FLUX_DIFFUSE(UTA,IBEAM,O1,WDIR)   = &
                          FLUX_DIFFUSE(UTA,IBEAM,O1,WDIR)   + ( DNDIRECT_FLUX_SCALED   - DNDIRECT_FLUX )
                  ENDDO
                ENDIF

              ENDIF
            ENDDO

!  Finish downwelling direct contribution

          ENDIF

!  Finish MV output

        ENDIF

!  Continuation point for avoiding direct beam calculation
! 455    CONTINUE

!  Finish directional loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_INTEGRATED_OUTPUT

!  End module

      END MODULE vlidort_intensity_m
