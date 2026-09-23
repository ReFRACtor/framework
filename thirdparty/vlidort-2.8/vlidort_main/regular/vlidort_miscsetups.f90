
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
! #            VLIDORT_MISCSETUPS (master, calling:)            #
! #              VLIDORT_DELTAMSCALE                            #
! #              VLIDORT_SSALBINIT                              #
! #              VLIDORT_QSPREP                                 #
! #              VLIDORT_PREPTRANS                              #
! #                                                             #
! #            VLIDORT_DIRECTRADIANCE                           #
! #                                                             #
! #            VLIDORT_PIMATRIX_SETUP                           #
! #            VLIDORT_PIMATRIX_SETUP_OMP  (Version 2.7)        #
! #            EMULT_MASTER                                     #
! #            EMULT_MASTER_OBSGEO                              #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. 
!     -- BRDF and SLEAVE arrays are defined locally, each Fourier (DIRECTRADIANCE routine)
!     -- Emult routines (multipliers for the Beam solution) Moved here

      MODULE vlidort_miscsetups_m

      PRIVATE :: VLIDORT_DELTAMSCALE, &
                 VLIDORT_SSALBINIT,   &
                 VLIDORT_QSPREP,      &
                 VLIDORT_PREPTRANS

      PUBLIC :: VLIDORT_MISCSETUPS, &
                VLIDORT_DIRECTRADIANCE, &
                VLIDORT_PIMATRIX_SETUP, &
                VLIDORT_PIMATRIX_SETUP_OMP, &
                EMULT_MASTER, &
                EMULT_MASTER_OBSGEO

      CONTAINS

      SUBROUTINE VLIDORT_MISCSETUPS ( &
        DO_SOLAR_SOURCES, DO_DELTAM_SCALING, DO_PLANE_PARALLEL,        & ! Input flags
        DO_REFRACTIVE_GEOMETRY,                                        & ! Input flags
        DO_PARTLAYERS, DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY,       & ! Input flags
        DO_SOLUTION_SAVING, DO_SPECIALIST_OPTION_3, DO_TOA_CONTRIBS,   & ! Input flags
        NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,     & ! Input Numbers
        NBEAMS, NMOMENTS, N_PARTLAYERS, NLAYERS_CUTOFF,                & ! Input Numbers
        MUELLER_INDEX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Input level/Stokes indices
        PARTLAYERS_OUTFLAG,  PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & ! Input Partlayers
        QUAD_STREAMS, USER_SECANTS, COS_SZANGLES, SUN_SZA_COSINES,     & ! Input streams, SZA cosines
        OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,    & ! Input Optical
        TAUGRID_INPUT, CHAPMAN_FACTORS, PARTIAL_CHAPFACS,              & ! Input Chapman
        DELTAU_VERT, PARTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL,         & ! Output from DELTAMSCALE
        TAUGRID, DELTAU_SLANT, DELTAU_SLANT_UNSCALED,                  & ! Output from DELTAMSCALE
        PARTAU_SLANT_UNSCALED, LEVELS_SOLARTRANS,                      & ! Output from DELTAMSCALE
        PARTIALS_SOLARTRANS, SOLARBEAM_BOATRANS, TRUNC_FACTOR, FAC1,   & ! Output from DELTAMSCALE
        OMEGA_GREEK,                                                   & ! Output from SSALBINIT
        DO_REFLECTED_DIRECTBEAM, BEAM_CUTOFF, TRANS_SOLAR_BEAM,        & ! Output QSPREP
        INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                     & ! output QSPREP
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                & ! Output PREPTRANS (Discrete Ords.)
        T_DELT_MUBAR, T_UTUP_MUBAR, T_UTDN_MUBAR,                      & ! Output PREPTRANS (Solar beams)
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                      & ! Output PREPTRANS (User-vza)
        CUMTRANS, ITRANS_USERM )                                         ! Output PREPTRANS (auxiliary)       

!  Miscellaneous set-ups
!     Delta-M scaling, average secant parameterization, various transmittances

!mick fix 9/19/2017 - the following arguments added to facilitate correction of direct flux:
!                     to input  - DO_SOLAR_SOURCES, PARTIAL_CHAPFACS
!                     to output - DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED, 
!                                 LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS,          &
                                 MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_SZANGLES, &
                                 MAXMOMENTS_INPUT, MAXMOMENTS, MAXSTOKES_SQ


      IMPLICIT NONE

!  Inputs
!  ------

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY

      LOGICAL, INTENT (IN) ::           DO_PARTLAYERS
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY

      LOGICAL, INTENT (IN) ::           DO_SOLUTION_SAVING
      LOGICAL, INTENT (IN) ::           DO_SPECIALIST_OPTION_3
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS

!  Input numbers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      INTEGER, INTENT (IN) ::           NLAYERS_CUTOFF

!  Input numbers and indices

      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

!  Input partial layer control

      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

!  Input streams, SZA cosines

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )

!  input optical

      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Inputs from Chapman factor routine

      DOUBLE PRECISION, INTENT (IN) ::  TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Output from the Delta-M scaling routine (DELTAMSCALE)
!  -----------------------------------------------------

!  Outputs (optical)

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) :: TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_SLANT_UNSCALED ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTAU_SLANT_UNSCALED ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Output truncation factors

      DOUBLE PRECISION, INTENT (OUT) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: FAC1 ( MAXLAYERS )

!  Rob fix 11/17/14. Added Argument for Diagnostic output

      DOUBLE PRECISION, INTENT (OUT) :: LEVELS_SOLARTRANS   ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SOLARBEAM_BOATRANS  ( MAXBEAMS )

!  Output from SSALBINIT.

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Output from the average-secant parameterization routine (QSPREP)
!  ----------------------------------------------------------------

!  output solar beam control

      LOGICAL         , INTENT (OUT) :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      INTEGER         , INTENT (OUT) :: BEAM_CUTOFF ( MAXBEAMS )

!  output average-secant parameterization

      DOUBLE PRECISION, INTENT (OUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )
 
!  output from the Transmittance routine (PREPTRANS)
!  -------------------------------------------------

!  Output transmittances (Discrete Ords.)

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  Output transmittances (Solar beams)

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Output transmittances (User vzangles)

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Auxiliary output

      DOUBLE PRECISION, INTENT (OUT) :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  miscellaneous setup operations, Master routine

!  Performance set-up is in VLIDORT_DERIVE_INPUTS (vlidort_inputs.f)
!      CALL VLIDORT_PERFORMANCE_SETUP

!  Delta-m scaling of input quantities
!    Rob fix 11/17/14 added argument SOLARBEAM_BOATRANS
!mick fix 9/19/2017 - the following arguments added to facilitate correction of direct flux:
!                     (1) DO_SOLAR_SOURCES, N_PARTLAYERS, & PARTIAL_CHAPFACS as inputs
!                     (2) DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED, 
!                         LEVELS_SOLARTRANS, & PARTIALS_SOLARTRANS as outputs

      CALL VLIDORT_DELTAMSCALE ( &
        DO_SOLAR_SOURCES, DO_DELTAM_SCALING, DO_PARTLAYERS,              & ! Input flags
        NSTOKES, NLAYERS, N_PARTLAYERS, NMOMENTS, NBEAMS, N_USER_LEVELS, & ! Input Numbers
        PARTLAYERS_OUTFLAG,  PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,     & ! Input Partlayers
        OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,      & ! Input Optical
        TAUGRID_INPUT, CHAPMAN_FACTORS, PARTIAL_CHAPFACS,                & ! Input Chapman
        DELTAU_VERT, PARTAU_VERT, TAUGRID, OMEGA_TOTAL,                  & ! Output optical
        GREEKMAT_TOTAL, FAC1, TRUNC_FACTOR,                              & ! Output optical
        DELTAU_SLANT, DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED,      & ! Output optical
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, SOLARBEAM_BOATRANS )       ! Output solar trans

!  initialise single scatter albedo terms
!    GREEKMAT x OMEGA

      CALL VLIDORT_SSALBINIT ( &
        NSTOKES, NLAYERS, NMOMENTS, MUELLER_INDEX, & ! Input numbers and indices
        OMEGA_TOTAL, GREEKMAT_TOTAL,               & ! Input optical
        OMEGA_GREEK )                                ! Output

!  Prepare quasi-spherical solar-beam attenuations

      CALL VLIDORT_QSPREP ( &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_SPECIALIST_OPTION_3, & ! Input flags
        NLAYERS, NLAYERS_CUTOFF, NBEAMS, COS_SZANGLES, SUN_SZA_COSINES,    & ! Input numbers, SZA cosines
        DELTAU_VERT, TAUGRID, DELTAU_SLANT,                                & ! Input optical
        DO_REFLECTED_DIRECTBEAM, BEAM_CUTOFF, TRANS_SOLAR_BEAM,            & ! Output solar beam
        INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA )                          ! output solar beam

!  Transmittances and Transmittance factors

      CALL VLIDORT_PREPTRANS ( &
        DO_SOLUTION_SAVING, DO_USER_STREAMS,                         & ! Input flags
        DO_OBSERVATION_GEOMETRY, DO_TOA_CONTRIBS,                    & ! Input flags
        NSTREAMS, NLAYERS, N_PARTLAYERS, NBEAMS, N_USER_STREAMS,     & ! Input numbers 
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, PARTLAYERS_LAYERIDX, & ! Input level control
        QUAD_STREAMS, USER_SECANTS, DELTAU_VERT, PARTAU_VERT,        & ! Input streams and optical
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,                  & ! Input solar beam
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,              & ! Output transmittances (Discrete Ords.)
        T_DELT_MUBAR, T_UTUP_MUBAR, T_UTDN_MUBAR,                    & ! Output transmittances (Solar beams)
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Output transmittances (User-vza)
        CUMTRANS, ITRANS_USERM )                                       ! Output auxiliary       

!  Finish

!      write(*,*)INITIAL_TRANS(1,1), AVERAGE_SECANT(1,1), LOCAL_CSZA(1,1)
!      write(*,*)INITIAL_TRANS(2,1), AVERAGE_SECANT(2,1), LOCAL_CSZA(2,1)
!      write(*,*)INITIAL_TRANS(23,1), AVERAGE_SECANT(23,1), LOCAL_CSZA(23,1)

      RETURN
      END SUBROUTINE VLIDORT_MISCSETUPS

!

      SUBROUTINE VLIDORT_DELTAMSCALE ( &
        DO_SOLAR_SOURCES, DO_DELTAM_SCALING, DO_PARTLAYERS,              & ! Input flags
        NSTOKES, NLAYERS, N_PARTLAYERS, NMOMENTS, NBEAMS, N_USER_LEVELS, & ! Input Numbers
        PARTLAYERS_OUTFLAG,  PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,     & ! Input Partlayers
        OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,      & ! Input Optical
        TAUGRID_INPUT, CHAPMAN_FACTORS, PARTIAL_CHAPFACS,                & ! Input Chapman
        DELTAU_VERT, PARTAU_VERT, TAUGRID, OMEGA_TOTAL,                  & ! Output optical
        GREEKMAT_TOTAL, FAC1, TRUNC_FACTOR,                              & ! Output optical
        DELTAU_SLANT, DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED,      & ! Output optical
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, SOLARBEAM_BOATRANS )       ! Output solar trans

!mick fix 9/19/2017 - the following arguments added to facilitate correction of direct flux:
!                     (1) DO_SOLAR_SOURCES, N_PARTLAYERS, & PARTIAL_CHAPFACS as inputs
!                     (2) DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED, 
!                         LEVELS_SOLARTRANS, & PARTIALS_SOLARTRANS as outputs

!  The deltam-scaling routine

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_USER_LEVELS, MAX_PARTLAYERS, MAXBEAMS, MAXMOMENTS_INPUT, &
                                 MAXMOMENTS, MAXSTOKES_SQ, ZERO, ONE, MAX_TAU_SPATH

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::           DO_PARTLAYERS

!  Input numbers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

!  Input partial layer control

      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

!  input optical

      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Inputs from Chapman factor routine

      DOUBLE PRECISION, INTENT (IN) ::  TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Outputs (optical)

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) :: TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_SLANT_UNSCALED ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTAU_SLANT_UNSCALED ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Output truncation factors

      DOUBLE PRECISION, INTENT (OUT) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: FAC1 ( MAXLAYERS )

!  Rob fix 11/17/14. Added Argument for Diagnostic output

      DOUBLE PRECISION, INTENT (OUT) :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  local variables

      DOUBLE PRECISION :: FDEL, FAC2, DNL1, FDNL1
      DOUBLE PRECISION :: DNM1, DELS, DT, XTD
      INTEGER          :: K, K1, K2, N, N1, L, UT, UTA, NM1, IB

!  Indexing

      INTEGER, DIMENSION(4) :: KTYPE1, KTYPE2

!mick - singularity buster output
      LOGICAL   :: SBUST(6)

!  Setup the two ktype arrays

      KTYPE1 = (/ 1, 6, 11, 16 /)
      KTYPE2 = (/ 2, 5, 12, 15 /)

!mick mod 9/19/2017 - put slant calculations at the bottom of subroutine
!                     (following LIDORT)

!mick fix - initialise
      GREEKMAT_TOTAL = ZERO

!  DELTAM SCALING
!  ==============

      IF ( DO_DELTAM_SCALING ) THEN

!  New section by R. Spurr, RT Solutions Inc.
!   Based in part on the code in VDISORT.

        TAUGRID(0) = ZERO
        NM1  = NMOMENTS+1
        DNM1 = DBLE(2*NM1+1)

!  Scaling for layer input
!  -----------------------

        DO N = 1, NLAYERS

          N1 = N - 1

!  Overall truncation factor

          FDEL = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          FAC2 = ONE - FDEL
          FAC1(N)         = ONE - FDEL * OMEGA_TOTAL_INPUT(N)
          TRUNC_FACTOR(N) = FDEL

!  Scale Greek Matrix entries

          DO L = 0, NMOMENTS
            DNL1  = DBLE(2*L + 1 )
            FDNL1 = FDEL * DNL1
            IF ( NSTOKES .GT. 1 ) THEN
              DO K = 1, 4
                K1 = KTYPE1(K)
                K2 = KTYPE2(K)
                GREEKMAT_TOTAL(L,N,K1) = ( GREEKMAT_TOTAL_INPUT(L,N,K1) - FDNL1 ) / FAC2
                GREEKMAT_TOTAL(L,N,K2) = GREEKMAT_TOTAL_INPUT(L,N,K2) / FAC2
              ENDDO
            ELSE
              GREEKMAT_TOTAL(L,N,1) = ( GREEKMAT_TOTAL_INPUT(L,N,1) - FDNL1 ) / FAC2
            ENDIF
          ENDDO

!  Maintain phase function normalization

          GREEKMAT_TOTAL(0,N,1) = ONE

!  Scale optical depth grid and single scatter albedo

          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N) * FAC1(N)
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N) * FAC2 / FAC1(N)
          TAUGRID(N)     = TAUGRID(N1) + DELTAU_VERT(N)

!  end layer loop

        ENDDO

!  Scaling for user-defined off-grid optical depths
!     (on-grid values have already been scaled)
!mick mod 9/19/2017 - code simplified (following LIDORT)

        IF ( DO_PARTLAYERS ) THEN
          DO UT = 1, N_PARTLAYERS
            N   = PARTLAYERS_LAYERIDX(UT)
            DT  = PARTLAYERS_VALUES(UT)
            XTD = DELTAU_VERT_INPUT(N) * DT
            PARTAU_VERT(UT) = XTD * FAC1(N)
          ENDDO
        ENDIF

!  Old code
!        IF ( DO_PARTLAYERS ) THEN
!          UT = 0
!          DO UTA = 1, N_USER_LEVELS
!            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
!              UT  = UT + 1
!              N   = PARTLAYERS_LAYERIDX(UT)
!              DT  = PARTLAYERS_VALUES(UT)
!              XTD = DELTAU_VERT_INPUT(N) * DT
!              PARTAU_VERT(UT) = XTD * FAC1(N)
!            ENDIF
!          ENDDO
!        ENDIF

!mick mod 9/19/2017 - put scaling of layer path thickness at bottom of subroutine
!                     (following LIDORT)

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

        TAUGRID(0) = ZERO
        DO N = 1, NLAYERS
          FAC1(N)        = ONE
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N)
          TAUGRID(N)     = TAUGRID_INPUT(N)
          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N)
          DO L = 0, NMOMENTS
!mick fix 1/21/2013 -  added if structure and else section
            IF ( NSTOKES .GT. 1 ) THEN
              DO K = 1, 4
                K1 = KTYPE1(K)
                K2 = KTYPE2(K)
                GREEKMAT_TOTAL(L,N,K1) = GREEKMAT_TOTAL_INPUT(L,N,K1)
                GREEKMAT_TOTAL(L,N,K2) = GREEKMAT_TOTAL_INPUT(L,N,K2)
              ENDDO
            ELSE
              GREEKMAT_TOTAL(L,N,1) = GREEKMAT_TOTAL_INPUT(L,N,1)
            ENDIF
          ENDDO
        ENDDO

!  Scaling for user-defined off-grid optical depths
!     (on-grid values have already been scaled)
!mick mod 9/19/2017 - code simplified (following LIDORT)

        IF ( DO_PARTLAYERS ) THEN
          DO UT = 1, N_PARTLAYERS
            N   = PARTLAYERS_LAYERIDX(UT)
            DT  = PARTLAYERS_VALUES(UT)
            XTD = DELTAU_VERT_INPUT(N) * DT
            PARTAU_VERT(UT) = XTD
          ENDDO
        ENDIF

!        IF ( DO_PARTLAYERS ) THEN
!          UT = 0
!          DO UTA = 1, N_USER_LEVELS
!            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
!              UT  = UT + 1
!              N   = PARTLAYERS_LAYERIDX(UT)
!              DT  = PARTLAYERS_VALUES(UT)
!              XTD = DELTAU_VERT_INPUT(N) * DT
!              PARTAU_VERT(UT) = XTD
!            ENDIF
!          ENDDO
!        ENDIF

      ENDIF

!mick fix 2/13/2012 - singularity busters added

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if
!        necessary

      DO N = 1, NLAYERS
        SBUST = .false.

        !Singularity buster for single scatter albedo
        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
          OMEGA_TOTAL(N) = 0.999999999D0
          SBUST(1) = .true.
        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
          OMEGA_TOTAL(N) = 1.0D-9
          SBUST(2) = .true.
        END IF

        !Singularity buster for asymmetry parameter
        !(1) Divide by 2L+1 where L = 1 to get the asym par
        GREEKMAT_TOTAL(1,N,1) = GREEKMAT_TOTAL(1,N,1)/3.0D0
        !(2) Modify the asym par if necessary
        IF (GREEKMAT_TOTAL(1,N,1) > 0.999999999D0) THEN
          GREEKMAT_TOTAL(1,N,1) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (GREEKMAT_TOTAL(1,N,1) < -0.999999999D0) THEN
          GREEKMAT_TOTAL(1,N,1) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((GREEKMAT_TOTAL(1,N,1) >= 0.0D0) .AND. &
                 (GREEKMAT_TOTAL(1,N,1) < 1.0D-9)) THEN
          GREEKMAT_TOTAL(1,N,1) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((GREEKMAT_TOTAL(1,N,1) < 0.0D0) .AND. &
                 (GREEKMAT_TOTAL(1,N,1) > -1.0D-9)) THEN
          GREEKMAT_TOTAL(1,N,1) = -1.0D-9
          SBUST(6) = .true.
        END IF
        !(3) Reconstruct the 1st-order phase func moment
        GREEKMAT_TOTAL(1,N,1) = 3.0D0*GREEKMAT_TOTAL(1,N,1)

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO

      ENDDO
!READ(*,*)

!mick mod 9/19/2017 - moved slant calculations from above (two places) to here
!                   - initialized all slant layer path thicknesses & solar transmissions here
!                   - added DO_SOLAR_SOURCES return (on hold right now)
!                   - copied updated code from LIDORT to calculate the slant layer path
!                     thicknesses & solar transmissions

      DELTAU_SLANT          = ZERO
      DELTAU_SLANT_UNSCALED = ZERO
      PARTAU_SLANT_UNSCALED = ZERO

      LEVELS_SOLARTRANS     = ZERO
      PARTIALS_SOLARTRANS   = ZERO
      SOLARBEAM_BOATRANS    = ZERO

!  If no solar terms, finish

      !IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Slant optical thickness values + scaling
!  ----------------------------------------

!  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT_UNSCALED(N,K,IB) = DELTAU_VERT_INPUT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

!rob fix 11/27/2014 - calculate SOLARBEAM_BOATRANS
!  Rob Fix, 7/18/17. level solar trans for all. Bug Fix 5/10/19 BOATRANS was set to DELS

      LEVELS_SOLARTRANS(0,1:NBEAMS) = ONE
      DO IB = 1, NBEAMS
         do N = 1, NLAYERS
           DELS = SUM(DELTAU_SLANT_UNSCALED(N,1:N,IB))
           IF ( DELS .le. MAX_TAU_SPATH ) LEVELS_SOLARTRANS(N,IB) = EXP(-DELS)
           if ( N==NLAYERS) SOLARBEAM_BOATRANS(IB) = LEVELS_SOLARTRANS(N,IB)
        enddo
      ENDDO

!  old code
!      DO IB = 1, NBEAMS
!         DELS = SUM(DELTAU_SLANT_UNSCALED(NLAYERS,1:NLAYERS,IB))
!         IF ( DELS .le. MAX_TAU_SPATH ) THEN
!            SOLARBEAM_BOATRANS(IB) = DELS
!         ENDIF
!      ENDDO

!  new Code for the PARTIALS_SOLARTRANS, PARTAU_SLANT_UNSCALED

        IF ( DO_PARTLAYERS ) THEN
          DO IB = 1, NBEAMS
            DO UT = 1, N_PARTLAYERS
              N   = PARTLAYERS_LAYERIDX(UT) ; N1 = N - 1
              XTD = DELTAU_VERT_INPUT(N) * PARTLAYERS_VALUES(UT)                ! Vertical OD in partial layer
              PARTAU_SLANT_UNSCALED(UT,N,IB) = XTD * PARTIAL_CHAPFACS(UT,N,IB)  ! Slant    OD in partial layer
              DO K = 1, N - 1
                PARTAU_SLANT_UNSCALED(UT,K,IB) = PARTIAL_CHAPFACS(UT,K,IB) * DELTAU_VERT_INPUT(K)  ! Slant ODs in other layers to TOA
              ENDDO
              DELS = SUM(PARTAU_SLANT_UNSCALED(UT,1:N,IB))
              IF ( DELS .le. MAX_TAU_SPATH ) PARTIALS_SOLARTRANS(UT,IB) = EXP(-DELS)
            ENDDO
          ENDDO
        ENDIF

!  Scale layer path thickness values, or not (just copy input)

      IF ( DO_DELTAM_SCALING ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            DO K = 1, N
              DELTAU_SLANT(N,K,IB) = DELTAU_SLANT_UNSCALED(N,K,IB) * FAC1(K)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            DO K = 1, N
              DELTAU_SLANT(N,K,IB) = DELTAU_SLANT_UNSCALED(N,K,IB)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish module

      RETURN
      END SUBROUTINE VLIDORT_DELTAMSCALE

!

      SUBROUTINE VLIDORT_SSALBINIT ( &
        NSTOKES, NLAYERS, NMOMENTS, MUELLER_INDEX, & ! Input numbers and indices
        OMEGA_TOTAL, GREEKMAT_TOTAL,               & ! Input optical
        OMEGA_GREEK )                                ! Output

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ

      IMPLICIT NONE

!  Input numbers and indices

      INTEGER, INTENT (IN) ::           NSTOKES 
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Input optical

      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GREEKMAT_TOTAL (0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  local variables

      INTEGER ::          N, L, O1, O2, OM
      DOUBLE PRECISION :: SUM

!  phase matrix-weighted OMEGA

      DO N = 1, NLAYERS
        DO L = 0, NMOMENTS
          DO O1 = 1, NSTOKES
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              SUM = OMEGA_TOTAL(N)*GREEKMAT_TOTAL(L,N,OM)
              OMEGA_GREEK(L,N,O1,O2)   = SUM
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_SSALBINIT

!

      SUBROUTINE VLIDORT_QSPREP ( &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_SPECIALIST_OPTION_3, & ! Input flags
        NLAYERS, NLAYERS_CUTOFF, NBEAMS, COS_SZANGLES, SUN_SZA_COSINES,    & ! Input numbers, SZA cosines
        DELTAU_VERT, TAUGRID, DELTAU_SLANT,                                & ! Input optical
        DO_REFLECTED_DIRECTBEAM, BEAM_CUTOFF, TRANS_SOLAR_BEAM,            & ! Output solar beam
        INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA )                          ! output solar beam

!  Pseudo-spherical approximation: Average-secant parameterization for solar beam atteunation

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_SZANGLES, MAXBEAMS, ZERO, ONE, MAX_TAU_SPATH

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_SPECIALIST_OPTION_3

!  Input numbers, SZA cosines

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NLAYERS_CUTOFF
      INTEGER, INTENT (IN) ::           NBEAMS
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )

!  Input Optical

      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  output solar beam control

      LOGICAL         , INTENT (OUT) :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      INTEGER         , INTENT (OUT) :: BEAM_CUTOFF ( MAXBEAMS )

!  output average-secant parameterization

      DOUBLE PRECISION, INTENT (OUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Local variables
!  ---------------

      INTEGER          :: N, K, IB
      DOUBLE PRECISION :: S_T_0, S_T_1, SEC0, TAU, TAU_SOLAR(MAXBEAMS)
      DOUBLE PRECISION :: TAUSLANT ( 0:MAXLAYERS, MAXBEAMS )

!  Specialist code
!  ---------------

!  Cutting off the solar source arbitrarily (only for Specialist Option 3)

      IF ( DO_SPECIALIST_OPTION_3 ) THEN
        DO IB = 1, NBEAMS
          BEAM_CUTOFF(IB) = NLAYERS_CUTOFF
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          BEAM_CUTOFF(IB) = NLAYERS
        ENDDO
      ENDIF

!  TOA  cosines

      DO IB = 1, NBEAMS
        LOCAL_CSZA(0,IB) = COS_SZANGLES(IB)
      ENDDO

!  plane-parallel case
!  -------------------

!  Bug 9/22/20. Local_csza must be set for all layers

      IF ( DO_PLANE_PARALLEL ) THEN

       DO IB = 1, NBEAMS
        SEC0 = ONE / COS_SZANGLES(IB)
        DO N = 1, NLAYERS
          TAUSLANT(N,IB) = TAUGRID(N) * SEC0
          IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              BEAM_CUTOFF(IB) = N
            ENDIF
            AVERAGE_SECANT(N,IB) = SEC0
            INITIAL_TRANS(N,IB)  = EXP ( - TAUGRID(N-1) * SEC0 )
!            LOCAL_CSZA(N,IB) = COS_SZANGLES(IB)     ! Line commented out, 9/22/20
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
!            LOCAL_CSZA(N,IB)     = ZERO             ! Line commented out, 9/22/20
          ENDIF
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)
        LOCAL_CSZA(1:NLAYERS,IB) = COS_SZANGLES(IB)  ! Line added here 9/22/20
       ENDDO

      ELSE

!  pseudo-spherical case
!  ---------------------

       DO IB = 1, NBEAMS

!  Get the total spherical attenuation from layer thickness sums

        TAUSLANT(0,IB) = ZERO
        DO N = 1, NLAYERS
          TAU = ZERO
          DO K = 1, N
            TAU = TAU + DELTAU_SLANT(N,K,IB)
          ENDDO
          TAUSLANT(N,IB) = TAU
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)

!  set up the average secant formulation

        S_T_1 = ZERO
        S_T_0 = ONE
!        BEAM_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              BEAM_CUTOFF(IB) = N
            ELSE
              S_T_1 = DEXP ( - TAUSLANT(N,IB) )
            ENDIF
            AVERAGE_SECANT(N,IB) = (TAUSLANT(N,IB)-TAUSLANT(N-1,IB)) / DELTAU_VERT(N)
            INITIAL_TRANS(N,IB)  = S_T_0
            S_T_0                = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
          ENDIF
        ENDDO

!  Set the Local solar zenith cosines
!  Distinguish between the refractive and non-refractive cases.
!  Bug 9/22/20. Local_csza must be set for all layers

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          DO N = 1, NLAYERS
             LOCAL_CSZA(N,IB) = SUN_SZA_COSINES(N,IB)
          ENDDO
        ELSE
          LOCAL_CSZA(1:NLAYERS,IB) = COS_SZANGLES(IB)
        ENDIF

!  (Old code, pre 9/22/20 bug)
!        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!          DO N = 1, NLAYERS
!            IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
!              LOCAL_CSZA(N,IB) = SUN_SZA_COSINES(N,IB)
!            ELSE
!              LOCAL_CSZA(N,IB) = ZERO
!            ENDIF
!          ENDDO
!        ELSE
!          DO N = 1, NLAYERS
!            IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
!              LOCAL_CSZA(N,IB) = COS_SZANGLES(IB)
!            ELSE
!              LOCAL_CSZA(N,IB) = ZERO
!            ENDIF
!          ENDDO
!        ENDIF

       ENDDO
      ENDIF

!  Set Direct Beam Flag and solar beam total attenuation to surface

      DO IB = 1, NBEAMS
        TRANS_SOLAR_BEAM(IB) = ZERO
        DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        IF ( .NOT.DO_SPECIALIST_OPTION_3 ) THEN
          IF ( TAU_SOLAR(IB) .LT. MAX_TAU_SPATH ) THEN
            TRANS_SOLAR_BEAM(IB) = EXP( - TAU_SOLAR(IB) )
            DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
          ENDIF
        ENDIF
      ENDDO

! debug
!      do n = 1, nlayers
!        write(*,*)AVERAGE_SECANT(N,1),  INITIAL_TRANS(N,1)
!      enddo
!      pause

!  finish

      RETURN
      END SUBROUTINE VLIDORT_QSPREP

!

      SUBROUTINE VLIDORT_PREPTRANS ( &
        DO_SOLUTION_SAVING, DO_USER_STREAMS,                         & ! Input flags
        DO_OBSERVATION_GEOMETRY, DO_TOA_CONTRIBS,                    & ! Input flags
        NSTREAMS, NLAYERS, N_PARTLAYERS, NBEAMS, N_USER_STREAMS,     & ! Input numbers 
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, PARTLAYERS_LAYERIDX, & ! Input level control
        QUAD_STREAMS, USER_SECANTS, DELTAU_VERT, PARTAU_VERT,        & ! Input streams and optical
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,                  & ! Input solar beam
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,              & ! Output transmittances (Discrete Ords.)
        T_DELT_MUBAR, T_UTUP_MUBAR, T_UTDN_MUBAR,                    & ! Output transmittances (Solar beams)
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Output transmittances (User-vza)
        CUMTRANS, ITRANS_USERM )                                       ! Output auxiliary       

!  Prepare transmittances and transmittance factors

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                 ZERO, ONE, MAX_TAU_QPATH, MAX_TAU_SPATH, MAX_TAU_UPATH

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLUTION_SAVING
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS

!  Input numbers

      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS

!  Input level control

      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

!  Input streams and optical

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT ( MAX_PARTLAYERS )

!  Input solar beam

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Output transmittances (Discrete Ords.)

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  Output transmittances (Solar beams)

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Output transmittances (User vzangles)

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Auxiliary output

      DOUBLE PRECISION, INTENT (OUT) :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER          :: N, UT, UM, IB, I, LUM
      DOUBLE PRECISION :: XT, SPHER, HELP

!  Local user index

      LUM = 1

!  Transmittance factors for discrete ordinate streams
!  ===================================================

!  New code by R. Spurr, RT Solutions, 12 April 2005
!    Required for the solution saving option
!  Partial layer code added 30 December 2005 by R. Spurr

!mick fix 7/23/2014 - initialized for packing

      T_DELT_DISORDS = ZERO
      T_DISORDS_UTDN = ZERO
      T_DISORDS_UTUP = ZERO

      ITRANS_USERM = ZERO
      T_DELT_USERM = ZERO
      T_UTDN_USERM = ZERO
      T_UTUP_USERM = ZERO

      CUMTRANS = ZERO

!  solution saving option
!  ----------------------

      IF ( DO_SOLUTION_SAVING ) THEN

!  whole layers

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            SPHER = DELTAU_VERT(N) / QUAD_STREAMS(I)
            IF ( SPHER .GT. MAX_TAU_QPATH ) THEN
              T_DELT_DISORDS(I,N) = ZERO
            ELSE
              T_DELT_DISORDS(I,N) = EXP ( - SPHER )
            ENDIF
          ENDDO
        ENDDO

!  Atmosphere Partial layers

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          XT = PARTAU_VERT(UT)
          DO I = 1, NSTREAMS
            HELP =  XT / QUAD_STREAMS(I)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DISORDS_UTDN(I,UT) = ZERO
            ELSE
              T_DISORDS_UTDN(I,UT) = EXP(-HELP)
            ENDIF
            HELP = ( DELTAU_VERT(N) - XT ) / QUAD_STREAMS(I)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DISORDS_UTUP(I,UT) = ZERO
            ELSE
              T_DISORDS_UTUP(I,UT) = EXP(-HELP)
            ENDIF
          ENDDO
        ENDDO

      ENDIF

!  Transmittance factors for average secant stream
!  ===============================================

!  start solar loop

      DO IB = 1, NBEAMS

!  Whole layer Transmittance factors
!  ---------------------------------

!  layer transmittance

       DO N = 1, NLAYERS
        IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
          T_DELT_MUBAR(N,IB) = ZERO
        ELSE
          SPHER = DELTAU_VERT(N) * AVERAGE_SECANT(N,IB)
          IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
            T_DELT_MUBAR(N,IB) = ZERO
          ELSE
            T_DELT_MUBAR(N,IB) = EXP ( - SPHER )
          ENDIF
        ENDIF
       ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  -----------------------------------------------------------------

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
           T_UTDN_MUBAR(UT,IB) = ZERO
           T_UTUP_MUBAR(UT,IB) = ZERO
         ELSE
           XT = PARTAU_VERT(UT)
           SPHER = XT * AVERAGE_SECANT(N,IB)
           IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
             T_UTDN_MUBAR(UT,IB) = ZERO
           ELSE
             T_UTDN_MUBAR(UT,IB) = EXP ( - SPHER )
           ENDIF
           SPHER = ( DELTAU_VERT(N) - XT ) * AVERAGE_SECANT(N,IB)
           IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
             T_UTUP_MUBAR(UT,IB) = ZERO
           ELSE
             T_UTUP_MUBAR(UT,IB) = EXP ( - SPHER )
           ENDIF
         ENDIF
        ENDDO

!  end solar beam loop

      ENDDO

!  Transmittances for User Streams
!  ===============================

!  return if not flagged

      IF ( .NOT. DO_USER_STREAMS ) RETURN

!  Initial transmittances divided by user streams
!  ----------------------------------------------

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
         DO N = 1, NLAYERS
          DO UM = 1, N_USER_STREAMS
           ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(UM)
          ENDDO
         ENDDO
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
         DO N = 1, NLAYERS
           ITRANS_USERM(N,LUM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(IB)
         ENDDO
        ENDDO
      ENDIF

!  Whole Layer transmittances

      DO N = 1, NLAYERS
       IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
        DO UM = 1, N_USER_STREAMS
         SPHER = DELTAU_VERT(N) * USER_SECANTS(UM)
         IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
          T_DELT_USERM(N,UM) = ZERO
         ELSE
          T_DELT_USERM(N,UM) = EXP ( - SPHER )
         ENDIF
        ENDDO
       ENDIF
      ENDDO

!  Cumulative tranmsittances (TOA contribution functions)
!     New section, 27 Janaury 2010

      if ( DO_TOA_CONTRIBS ) THEN
        DO UM = 1, N_USER_STREAMS
          CUMTRANS(1,UM) = ONE
          DO N = 1, NLAYERS - 1
            CUMTRANS(N+1,UM) = CUMTRANS(N,UM) * T_DELT_USERM(N,UM)
          ENDDO
        ENDDO
      endif

!  Partial Layer transmittances for off-grid optical depths

      DO UT = 1, N_PARTLAYERS
        N  = PARTLAYERS_LAYERIDX(UT)
        XT = PARTAU_VERT(UT)
        DO UM = 1, N_USER_STREAMS
          SPHER = XT * USER_SECANTS(UM)
          IF ( SPHER .GT. MAX_TAU_UPATH ) THEN
            T_UTDN_USERM(UT,UM) = ZERO
          ELSE
            T_UTDN_USERM(UT,UM) = EXP ( - SPHER )
          ENDIF
          SPHER = ( DELTAU_VERT(N) - XT ) * USER_SECANTS(UM)
          IF ( SPHER .GT. MAX_TAU_UPATH ) THEN
            T_UTUP_USERM(UT,UM) = ZERO
          ELSE
            T_UTUP_USERM(UT,UM) = EXP ( - SPHER )
          ENDIF
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_PREPTRANS

!

      SUBROUTINE VLIDORT_DIRECTRADIANCE ( DO_USER_STREAMS, DO_REFRACTIVE_GEOMETRY,  &
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING,          & ! input
            DO_SL_ISOTROPIC, DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION, & ! input
            NSTOKES, NSTREAMS, NBEAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK,         & ! input
            FOURIER, MUELLER_INDEX, FLUX_FACTOR, FLUXVEC, DELTA_FACTOR,             & ! input
            SZA_LOCAL_INPUT, COS_SZANGLES, TRANS_SOLAR_BEAM, DO_REFLEC_DIRBEAM,     & ! Input solar beam
            LAMBERTIAN_ALBEDO, BRDF_F_0, USER_BRDF_F_0,                             & ! input surface
            TRANS_ATMOS_FINAL, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,       & ! input sleave
            ATMOS_ATTN, DIRECT_BEAM, USER_DIRECT_BEAM, SL_QUADTERM, SL_USERTERM )     ! Output

!  Direct-beam reflected values.
!   --> Generalized to include the BRDF case as well, RT Solutions Inc., 26 July 2010
!   --. New Surface-Leaving arguments added, 17 May 2012

!  Streamlined for Version 2.8, 5 July 2016.
!    --> BRDF_F and USER_BRDF_F arguments removed.

!  4/9/19. revision to include DO_WATER_LEAVING flag
!  4/9/19. renamed, distinguish between reflected beam and surface-leaving radiance output
  
!  1/31/21. Version 2.8.3. 
!     -- BRDF and SLEAVE arrays are defined locally, each Fourier 
!     -- remove MAXMOMENTS from parameter inclusions
!     -- Use post-processing masks. Drop DO_OBSERVATION_GEOMETRY
!     -- TRANS_ATMOS_FINAL added argument (intent(inout)). 
!           -- This is established later for Fourier 0, but then used again here for Fourier > 0

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAXSTOKES, MAXSTREAMS, MAXBEAMS, MAX_SZANGLES, &
                                 MAX_USER_STREAMS, MAXSTOKES_SQ, ZERO, ONE, FOUR, PI4, DEG_TO_RAD

      IMPLICIT NONE

!  inputs
!  ------

!  Flags
!     -- 1/31/21. Version 2.8.3. Drop DO_OBSERVATION_GEOMETRY
!     -- 5/24/21. Version 2.8.3. Add the TF_iteration flag

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LEAVING
      
      LOGICAL, INTENT (IN) ::           DO_SL_ISOTROPIC
      LOGICAL, intent(in)  ::           DO_WATER_LEAVING
      LOGICAL, intent(in)  ::           DO_EXTERNAL_WLEAVE
      LOGICAL, intent(in)  ::           DO_TF_ITERATION

!  numbers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS

!   -- 1/31/21. Version 2.8.3. Post-processing masks

      INTEGER, intent(in)  ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  flux-factors, indices

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  DELTA_FACTOR

!   solar beam stuff

      DOUBLE PRECISION, INTENT (IN) ::  SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  TRANS_SOLAR_BEAM  ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::           DO_REFLEC_DIRBEAM ( MAXBEAMS )

!  surface reflectance inputs
!    -- 1/31/21. Version 2.8.3. BRDF and SLEAVE arrays are defined locally, each Fourier

      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F_0      ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )

!  New Surface-Leaving stuff 17 May 2012
!  Isotropic Surface leaving term (if flag set)
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams
!    -- 1/31/21. Version 2.8.3. SLEAVE arrays are defined locally, each Fourier

      DOUBLE PRECISION, INTENT (IN) ::  SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SLTERM_F_0       ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_SLTERM_F_0  ( MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  5/24/21. Version 2.8.3. Must include the Adjustable tranmsitted flux (intent(inout))

      DOUBLE PRECISION, intent(inout) :: TRANS_ATMOS_FINAL ( MAXBEAMS )

!  output arguments
!  ----------------

!  output, reflected direct-beams

      DOUBLE PRECISION, INTENT (OUT) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  SL terms (not for Waterleaving). 4/9/19

      DOUBLE PRECISION, intent(out) :: SL_QUADTERM ( MAXSTREAMS,       MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, intent(out) :: SL_USERTERM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: X0_FLUX, X0_BOA, ATTN, REFLEC, SUM, SL, HELP, LOCAL_TRANS ( MAXBEAMS )
      INTEGER          :: I, UI, O1, O2, IB, M, OM, LUI
      LOGICAL          :: DO_CALC_SLTERMS

!  Initialize
!  ----------

!  Safety first.
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

      DO IB = 1, NBEAMS
         DO I = 1, NSTREAMS
            DIRECT_BEAM(I,IB,1:NSTOKES) = ZERO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
           DO LUI = 1, N_PPSTREAMS
              UI = PPSTREAM_MASK(LUI,IB)
              USER_DIRECT_BEAM(UI,IB,1:NSTOKES) = ZERO
           ENDDO
        ENDIF
     ENDDO
     SL_QUADTERM = zero ; SL_USERTERM = zero

!  return if no surface or surface-leaving
!  (add surface-leaving 4/9/19)

      M = FOURIER
      IF ( .not. DO_INCLUDE_SURFACE .and. .not. DO_SURFACE_LEAVING ) RETURN

!  Calculation flag for surface-leaving

!  1/31/21. Version 2.8.3.. Condition relaxed to include all Fourier terms for water-leaving.
!    .. Land-surface Flourescence - Only for Fourier zero
!    -- Water-leaving, Only for Fourier > 0 (Adjusted), or all Fourier (Non-adjusted, external)

!    -- 1/31/21. Version 2.8.3. Generalize the Water-leaving to allow all Fourier m > 0 (non-isotropic)
!    -- 1/31/21. Version 2.8.3. For iterated non-isotropic and M > 0, need to multiply by TRANS_ATMOS_FINAL
!    -- 5/24/21. Version 2.8.3. Final determination

      DO_CALC_SLTERMS = .false. ; LOCAL_TRANS = one
      IF ( DO_SURFACE_LEAVING ) THEN
         IF ( DO_WATER_LEAVING ) THEN
            IF ( DO_EXTERNAL_WLEAVE ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  IF ( FOURIER .eq. 0 ) DO_CALC_SLTERMS = .true.
               ELSE
                  DO_CALC_SLTERMS = .true.
               ENDIF
            ELSE
               IF ( FOURIER.gt.0 ) then
                  DO_CALC_SLTERMS = .true.
                  IF ( DO_TF_ITERATION ) LOCAL_TRANS(1:NBEAMS) = TRANS_ATMOS_FINAL(1:NBEAMS)
               ENDIF
            ENDIF
         ELSE
           IF ( FOURIER.eq.0 )  DO_CALC_SLTERMS = .true.
         ENDIF
      ENDIF

!  here is the old condition
!      DO_CALC_SLTERMS = .false.
!      IF ( M.eq.0.and.DO_SURFACE_LEAVING ) then
!         IF ( .not.DO_WATER_LEAVING .or. (DO_WATER_LEAVING.and.DO_EXTERNAL_WLEAVE)) DO_CALC_SLTERMS = .true.
!      ENDIF

!  Attenuation of solar beam
!  -------------------------

!  New code to deal with refractive geometry case
!   R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
         IF ( DO_REFLEC_DIRBEAM(IB) ) THEN

!  Definition of BOA SZA cosine includes refractive geometry case

            IF ( DO_REFRACTIVE_GEOMETRY ) THEN
               X0_BOA = COS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
            ELSE
               X0_BOA = COS_SZANGLES(IB)
            ENDIF

!  There should be no flux factor here.
!    Bug fixed 18 November 2005. Earlier Italian Job!!
!       X0_FLUX        = FOUR * X0_BOA * FLUX_FACTOR / DELTA_FACTOR

            X0_FLUX        = FOUR * X0_BOA / DELTA_FACTOR
            X0_FLUX        = FLUX_FACTOR * X0_FLUX / PI4             ! New
            ATTN           = X0_FLUX * TRANS_SOLAR_BEAM(IB)
            ATMOS_ATTN(IB) = ATTN

!  Lambertian case
!  ---------------

!  Set value for Fourier = 0
!  Natural light only, so only the first Stokes component is nonzero
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask

            IF ( DO_LAMBERTIAN_SURFACE ) THEN
               REFLEC = ATTN * LAMBERTIAN_ALBEDO
               IF ( FOURIER .EQ. 0 ) THEN
                  DIRECT_BEAM(1:NSTREAMS,IB,1) = REFLEC
                  IF ( DO_USER_STREAMS ) THEN
                    DO LUI = 1, N_PPSTREAMS
                      UI = PPSTREAM_MASK(LUI,IB)
                      USER_DIRECT_BEAM(UI,IB,1) = REFLEC
                    ENDDO
                  ENDIF
               ENDIF
            ENDIF

!  Non-Lambertian case
!  -------------------

            IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN

!  Solar beam reflected into quad directions
!    -- 1/31/21. Version 2.8.3. BRDF arrays defined locally, remove M=Fourier index

               DO I = 1, NSTREAMS
                  DO O1 = 1, NSTOKES
                     SUM = ZERO
                     DO O2 = 1, NSTOKES
                        OM = MUELLER_INDEX(O1,O2) ; SUM = SUM + FLUXVEC(O2)*BRDF_F_0(OM,I,IB)
                     ENDDO
                     DIRECT_BEAM(I,IB,O1) = ATTN * SUM
                  ENDDO
               ENDDO

!  Solar beam reflected into User directions
!    -- 1/31/21. Version 2.8.3. BRDF arrays defined locally, remove M=Fourier index
!    -- 1/31/21. Version 2.8.3.  Use post-processing mask

               IF ( DO_USER_STREAMS ) THEN
                  DO LUI = 1, N_PPSTREAMS
                     UI = PPSTREAM_MASK(LUI,IB)
                     DO O1 = 1, NSTOKES
                        SUM = ZERO
                        DO O2 = 1, NSTOKES
                           OM = MUELLER_INDEX(O1,O2)
                           SUM = SUM + FLUXVEC(O2)*USER_BRDF_F_0(OM,UI,IB)
                        ENDDO
                        USER_DIRECT_BEAM(UI,IB,O1) = ATTN * SUM
                     ENDDO
                  ENDDO
               ENDIF

!  End Surface-type

            ENDIF

!  End reflected direct-beam clause

         ENDIF
      
!  New Surface-Leaving stuff 17 May 2012
!  -------------------------------------

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases
!  4/9/19. Not done here for  water-leaving contributions (done later)

!  1/31/21. Version 2.8.3. Some changes
!     -- SLEAVE Fourier arrays defined locally, FOURIER index dropped
!     -- Use post-processing mask
!     -- Calculation is done every Fourier for non-isotropic water-leaving.
!     -- for iterated non-isotropic, need to multiply by trans_atmos_final(Ibeam) for M > 0

         IF ( DO_CALC_SLTERMS ) THEN
            HELP = LOCAL_TRANS(IB) * FLUX_FACTOR / DELTA_FACTOR
            IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
               DO O1 = 1, NSTOKES
                  SL = SLTERM_ISOTROPIC(O1,IB) * HELP
                  SL_QUADTERM(1:NSTREAMS,IB,O1) = SL
                  IF ( DO_USER_STREAMS ) THEN
                     DO LUI = 1, N_PPSTREAMS
                        UI = PPSTREAM_MASK(LUI,IB) ; SL_USERTERM(UI,IB,O1) = SL
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO O1 = 1, NSTOKES
                  DO I = 1, NSTREAMS
                     SL_QUADTERM(I,IB,O1) = SLTERM_F_0(O1,I,IB) * HELP
                  ENDDO
                  IF ( DO_USER_STREAMS ) THEN
                     DO LUI = 1, N_PPSTREAMS
                        UI = PPSTREAM_MASK(LUI,IB)
                        SL_USERTERM(UI,IB,O1) = USER_SLTERM_F_0(O1,UI,IB) * HELP
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

!  end direct beam calculation

      ENDDO

!  finish

      RETURN
      END SUBROUTINE VLIDORT_DIRECTRADIANCE

!

      SUBROUTINE VLIDORT_PIMATRIX_SETUP ( &
        DO_REFRACTIVE_GEOMETRY, DO_USER_STREAMS, FOURIER,             & ! Input flags and Fourier
        NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS,           & ! Input numbers
        QUAD_STREAMS, USER_STREAMS, COS_SZANGLES, SUN_SZA_COSINES,    & ! Input streams and angles
        NMOMENTS, MUELLER_INDEX, DMAT,                                & ! auxiliary inputs
        PI_XQP, PI_XQM, PI_XUP, PI_XUM, PI_X0P,                       & ! Output Pi-Matrix stuff
        PI_XQM_POST, PI_XQM_PRE, PI_XQP_PRE, PI_XUM_POST, PI_XUP_PRE )  ! Output Pi-Matrix stuff

!  Notes
!  -----

!  This is equivalent to the Legendre setup modules in LIDORT

!---------old comments ---------------------------------
!  This needs modification for the case of Refractive atmosphere,
!  because the solar zenith angle is not constant.
!-------------------------------------------------------

!  PI-matrix setup, following recipe of Siewert (1982).
!  Single Fourier component only.

!  Tested against benchmark results in Vestrucci & Siewert (1984), Problem

!  original coding, September 2002, R. Spurr SAO
!  Multibeam SZA coding, July 2004, R. Spurr SAO
!  Coding for refractive geometry case, R. Spurr, RT Solutions, May 2005
!    -- Currently not enabled (2.7/2.8) --> SUN_SZA_COSINES, MUELLER_INDEX are not used

!  Streamlined for Version 2.8, 05 July 2016.
!     Currently not used.

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, MAXBEAMS, &
                                 MAXMOMENTS, MAX_SZANGLES, MAX_ALLSTRMS_P1, ZERO, HALF, ONE, TWO, THREE

      IMPLICIT NONE

!  Flags and Fourier


      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      INTEGER, INTENT (IN) ::           FOURIER

!  Numbers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS

      INTEGER, INTENT (IN) ::           NMOMENTS

!  Indices and Dmatrix from VLIDORT_DERIVE_INPUTS

      INTEGER         , INTENT (IN)  :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN)  :: DMAT          ( MAXSTOKES, MAXSTOKES )

!  Angles and streams

      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES    ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS    ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS    ( MAX_USER_STREAMS )

!  Output (generalized spherical functions)

      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM_POST ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP_PRE  ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Local matrices

      DOUBLE PRECISION :: XLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: YLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: ZLM ( 0:MAXMOMENTS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI  ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PHI ( 0:MAXMOMENTS, 0:MAXMOMENTS )
      DOUBLE PRECISION :: PINORM ( 0:MAXMOMENTS, 0:MAXMOMENTS )

!  Older save statements (pre version 2.7)

      DOUBLE PRECISION, SAVE :: DF_L ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_2LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LSQM4 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_RT_LP3XLM1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: ZHELP ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: RT_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: UPXSQ_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PLEG20 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: M2X_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: XUMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X_RT_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PIMM_KM ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X2LP1 ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PISIGN ( 0:MAXMOMENTS, 0:MAXMOMENTS )

!  local variables

      DOUBLE PRECISION :: FAC, XSQ, DF_LPM, DF_LP1MM, ZLM_VALUE
      DOUBLE PRECISION :: HX, HY, HZ, H1, H2, RCONST_20, RCONST_21
      DOUBLE PRECISION :: PL20, RL20, PL21, RL21, TL21, XM(MAX_ALLSTRMS_P1)
      INTEGER          :: M, I, I1, L, SK1, SK2, SK3, SK4, N
      INTEGER          :: NS, NSA, M1, IB, IP

!  Set integer M = Fourier number

      M = FOURIER
      NSA = 0

!  Local NSTOKES

!      LOCAL_NSTOKES = NSTOKES
!      IF ( NSTOKES .GT. 1 ) LOCAL_NSTOKES = 4

!  total number of angles

      IF ( DO_USER_STREAMS ) THEN
        NS = NSTREAMS + N_USER_STREAMS
      ELSE
        NS = NSTREAMS
      ENDIF

!  add solar streams, depends on the use of refractive geometry

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        NSA = NS + NBEAMS*NLAYERS
      ELSE
        NSA = NS + NBEAMS
      ENDIF

!  constants

      RCONST_20 = 0.25D0 * SQRT(6.0D0)
      RCONST_21 = 2.0D0 * RCONST_20

!  Coefficient matrix PINORM (for normalization) = [ (L-m)!/(L+m)! ] ^1/
!  ---------------------------------------------------------------------

!  .. first entry = 1 / (2m)!    ---- be careful with overflow

      PHI(M,M) = ONE
      FAC = ONE
      DO L = 1, 2*M
        FAC = FAC * DBLE(L)
      ENDDO
      PHI(M,M) = PHI(M,M) / FAC

!  .. Other entries by recurrence

      DO L = M + 1, NMOMENTS
        FAC = DBLE(L-M)/DBLE(L+M)
        PHI(L,M) = PHI(L-1,M)*FAC
      ENDDO

!  .. Square root

      DO L = M , NMOMENTS
        PINORM(L,M) = SQRT(PHI(L,M))
      ENDDO

!  Additional saved quantities (to commons in VLIDORT_PISETUP.VARS)
!  ----------------------------------------------------------------

!  Only for the first fundamental harmonic

      IF ( M .EQ. 0 ) THEN

!  Sign matrix

        DO M1 = 0, NMOMENTS
          DO L = M1, NMOMENTS
            IF (MOD((L-M1),2).EQ.0) THEN
              PISIGN(L,M1) = ONE
            ELSE
              PISIGN(L,M1) = -ONE
            ENDIF
          ENDDO
        ENDDO

!  floating point integer values

        DO L = 0, NMOMENTS
          DF_L(L)    = DBLE(L)
          DF_LP1(L)  = DBLE(L+1)
          DF_2LP1(L) = DBLE(2*L+1)
        ENDDO

        DF_LSQM4(2) = ZERO
        DO L = 3, NMOMENTS
          DF_LSQM4(L) = SQRT(DBLE(L*L-4))
        ENDDO

        DO L = 2, NMOMENTS
          DF_RT_LP3XLM1(L) = SQRT(DBLE((L+3)*(L-1)))
          ZHELP(L) = TWO*DF_2LP1(L)/DF_LP1(L)/DF_L(L)
        ENDDO

      ENDIF

!  local array of all streams (every Fourier component)

      DO I = 1, NSTREAMS
        XM(I) = QUAD_STREAMS(I)
      ENDDO
      IF ( DO_USER_STREAMS ) THEN
        DO I = 1, N_USER_STREAMS
          XM(I+NSTREAMS) = USER_STREAMS(I)
        ENDDO
      ENDIF

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        DO N = 1, NLAYERS
!          DO IB = 1, NBEAMS
!            IP = NS + NBEAMS * (N-1) + IB
!            XM(IP) = SUN_SZA_COSINES(N,IB)
!          ENDDO
!        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          IP = IB + NS
          XM(IP) = COS_SZANGLES(IB)
        ENDDO
      ENDIF

!  factors associated with stream values
!   Special case when XM(I) = ONE, factors are zeroed, and
!    those that give singularities are avoided later (see below)
!   R. Spurr and V. Natraj, 16 january 2006

      IF ( M .EQ. 0 )  THEN
        DO I = 1, NSA
          XSQ = XM(I) * XM(I)
          PLEG20(I)        = HALF*(THREE*XSQ-ONE)
          UMXSQ(I)         = ONE - XSQ
          XUMXSQ(I)        = UMXSQ(I) * XM(I)
          IF ( XM(I) .EQ. ONE ) THEN
           RT_UMXSQ(I)      = ZERO
           X_RT_UMXSQ(I)    = ZERO
           UPXSQ_D_UMXSQ(I) = ZERO
           M2X_D_UMXSQ(I)   = ZERO
          ELSE
           RT_UMXSQ(I)      = SQRT(UMXSQ(I))
           X_RT_UMXSQ(I)    = RT_UMXSQ(I) * XM(I)
           UPXSQ_D_UMXSQ(I) = (ONE + XSQ) / UMXSQ(I)
           M2X_D_UMXSQ(I)   = - TWO * XM(I) / UMXSQ(I)
          ENDIF
          DO L = 0, NMOMENTS
            X2LP1(L,I) = XM(I) * DF_2LP1(L)
          ENDDO
        ENDDO

!  D-matrix and Mueller indices. NOW DONE in VLIDORT_DERIVE_INPUTS
!        DO SK1 = 1, MAXSTOKES
!          DO SK2 = 1, MAXSTOKES
!            MUELLER_INDEX(SK1,SK2) = MAXSTOKES*(SK1-1) + SK2
!            DMAT(SK1,SK2) = ZERO
!          ENDDO
!          IF ( SK1.GT.2) THEN
!            DMAT(SK1,SK1) = -ONE
!          ELSE
!            DMAT(SK1,SK1) = ONE
!          ENDIF
!          MUELLER_DIAGONAL_INDEX(SK1) = MUELLER_INDEX(SK1,SK1)
!        ENDDO

      ENDIF

!  XYZ matrices
!  ------------

!  Inverse XLM diagonal matrices (Siewert (1982), Eq. 35a)
!  YLM diagonal matrices (Siewert (1982), Eq. 35b)
!  ZLM matrices. (Siewert (1982), Eq. 35c)

!  .. for the azimuth-independent harmonic

      IF ( M .EQ. 0 ) THEN

        DO L = 2, NMOMENTS
          XLM_DIAG(L,1) = DF_LP1(L)
          XLM_DIAG(L,2) = DF_RT_LP3XLM1(L)
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
          DO SK1 = 1, NSTOKES
            XLM_DIAG(L,SK1) = ONE/XLM_DIAG(L,SK1)
          ENDDO
        ENDDO

        DO L = 2, NMOMENTS
          YLM_DIAG(L,1) = DF_L(L)
          YLM_DIAG(L,2) = DF_LSQM4(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

!  .. for the other harmonics

      ELSE

        DO L = M, NMOMENTS
          DF_LP1MM = DBLE ( L + 1 - M )
          XLM_DIAG(L,1) = DF_LP1MM
          XLM_DIAG(L,1) = ONE/XLM_DIAG(L,1)
          IF ( L .EQ. 1 ) THEN
            XLM_DIAG(L,2) = ZERO
          ELSE
            XLM_DIAG(L,2) = DF_RT_LP3XLM1(L) * DF_LP1MM / DF_LP1(L)
            XLM_DIAG(L,2) = ONE / XLM_DIAG(L,2)
          ENDIF
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
        ENDDO

        DO SK1 = 1, NSTOKES
          YLM_DIAG(M,SK1) = ZERO
        ENDDO
        DO L = M + 1, NMOMENTS
          DF_LPM = DBLE ( L + M )
          YLM_DIAG(L,1) = DF_LPM
          YLM_DIAG(L,2) = DF_LPM * DF_LSQM4(L) / DF_L(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

        DO L = 2, NMOMENTS
          ZLM_VALUE = DBLE(M) * ZHELP(L)
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              ZLM(L,SK1,SK2) = ZERO
            ENDDO
          ENDDO
          ZLM(L,2,3) = ZLM_VALUE
          ZLM(L,3,2) = ZLM_VALUE
        ENDDO

      ENDIF

!  PI calculation
!  --------------

!  Initialise all the PI matrices for given harmonic

      DO I = 1, NSA
        DO L = M, NMOMENTS
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI(L,I,SK1,SK2) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!mick - alternate PI initialization
      PI = ZERO

!  M = 0 component. [Siewert (1982), Eqs (28) to (31)]

      IF ( M .EQ. 0 ) THEN

!  .. L = 0
        DO I = 1, NSA
          PI(0,I,1,1) = PI(0,I,1,1) + ONE
          PI(0,I,4,4) = PI(0,I,4,4) + ONE
        ENDDO
!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + XM(I)
          PI(1,I,4,4) = PI(1,I,4,4) + XM(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL20 = PLEG20(I)
          RL20 = RCONST_20*UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL20
          PI(2,I,2,2) = PI(2,I,2,2) + RL20
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                PI(L+1,I,SK1,SK2) = ( HX - HY ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M = 1 component. [Siewert (1982), Eqs (32) to (34), plus (35)]

      ELSE IF ( M .EQ. 1 ) THEN

!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + RT_UMXSQ(I)
          PI(1,I,4,4) = PI(1,I,4,4) + RT_UMXSQ(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL21 =   THREE     * X_RT_UMXSQ(I)
          RL21 = - RCONST_21 * X_RT_UMXSQ(I)
          TL21 = + RCONST_21 * RT_UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL21
          PI(2,I,2,2) = PI(2,I,2,2) + RL21
          PI(2,I,2,3) = PI(2,I,2,3) + TL21
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
          PI(2,I,3,2) = PI(2,I,2,3)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M > 1 components. [Siewert (1982), Eqs (36) to (38), plus (35)]
!  Limiting case of XM(I) = ONE requires special treatment to avoid NaN
!   R. Spurr and V. Natraj, 16 january 2006

      ELSE

!  .. L = M
        IF ( M .EQ. 2 ) THEN
          DO I = 1, NSA
            PIMM_11(I) =   THREE     * UMXSQ(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) =   RCONST_21
            ELSE
              PIMM_KM(I) =   RCONST_21 * UMXSQ(I)
            ENDIF
          ENDDO
        ELSE
          H1 = DF_2LP1(M-1)
          H2 = DF_LP1(M-1)*H1/DF_RT_LP3XLM1(M-1)
          DO I = 1, NSA
            PIMM_11(I) = H1 * RT_UMXSQ(I) * PIMM_11(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) = ZERO
            ELSE
              PIMM_KM(I) = H2 * RT_UMXSQ(I) * PIMM_KM(I)
            ENDIF
          ENDDO
        ENDIF

        DO I = 1, NSA
          PI(M,I,1,1) = PI(M,I,1,1) + PIMM_11(I)
          IF ( XM(I).EQ.ONE ) THEN
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * TWO
            PI(M,I,2,3) = PI(M,I,2,3) - PIMM_KM(I) * TWO
          ELSE
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * UPXSQ_D_UMXSQ(I)
            PI(M,I,2,3) = PI(M,I,2,3) + PIMM_KM(I) * M2X_D_UMXSQ(I)
          ENDIF
          PI(M,I,3,3) = PI(M,I,2,2)
          PI(M,I,4,4) = PI(M,I,1,1)
          PI(M,I,3,2) = PI(M,I,2,3)
        ENDDO
!  .. L = M + 1
        IF ( M .LT. NMOMENTS ) THEN
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(M,I)*PI(M,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(M,SK1,SK3)*PI(M,I,SK3,SK2)
                ENDDO
                PI(M+1,I,SK1,SK2) = ( HX + HZ ) * XLM_DIAG(M,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!  .. L > M + 1
        DO L = M + 1, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Normalized output.
!  ------------------

      DO L = M, NMOMENTS

!  .. at quadrature streams

        DO I = 1, NSTREAMS

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI_XQP(L,I,SK1,SK2) = PI(L,I,SK1,SK2) * PINORM(L,M)
            ENDDO
          ENDDO

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES

              H1 = ZERO
              H2 = ZERO
              DO SK3 = 1, NSTOKES
                H1 = H1 + DMAT(SK1,SK3) * PI_XQP(L,I,SK3,SK2)
                H2 = H2 + PI_XQP(L,I,SK1,SK3) * DMAT(SK3,SK2)
              ENDDO

              PI_XQP_PRE(L,I,SK1,SK2)  = H1
              PI_XQM_PRE(L,I,SK1,SK2)  = H1 * PISIGN(L,M)
              PI_XQM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

              H1 = ZERO
              DO SK3 = 1, NSTOKES
                H2 = ZERO
                DO SK4 = 1, NSTOKES
                  H2 = H2 + PI_XQP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                ENDDO
                H1 = H1 + DMAT(SK1,SK3)*H2
              ENDDO
              PI_XQM(L,I,SK1,SK2)  = H1 * PISIGN(L,M)

            ENDDO
          ENDDO

        ENDDO

!  .. at positive user_defined angles

        IF ( DO_USER_STREAMS ) THEN

          DO I = 1, N_USER_STREAMS
            I1 = I + NSTREAMS

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_XUP(L,I,SK1,SK2) = PI(L,I1,SK1,SK2) * PINORM(L,M)
              ENDDO
            ENDDO

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES

                H2 = ZERO
                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H1 = H1 + DMAT(SK1,SK3) * PI_XUP(L,I,SK3,SK2)
                  H2 = H2 + DMAT(SK3,SK2) * PI_XUP(L,I,SK1,SK3)
                ENDDO
                PI_XUP_PRE(L,I,SK1,SK2)  = H1
                PI_XUM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H2 = ZERO
                  DO SK4 = 1, NSTOKES
                    H2 = H2 + PI_XUP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                  ENDDO
                  H1 = H1 + DMAT(SK1,SK3)*H2
                ENDDO
                PI_XUM(L,I,SK1,SK2) = H1 * PISIGN(L,M)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

!  .. at solar zenith angles

!  depends on use of refractive geometry

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!          DO N = 1, NLAYERS
!            DO IB = 1, NBEAMS
!              IP = NS + NBEAMS * (N-1) + IB
!              DO SK1 = 1, NSTOKES
!                DO SK2 = 1, NSTOKES
!                  PI_X0P(L,IB,N,SK1,SK2) =
!     &                 PI(L,IP,SK1,SK2) * PINORM(L,M)
!                ENDDO
!              ENDDO
!            ENDDO
!          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            IP = IB + NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_X0P(L,IB,1,SK1,SK2) = &
                     PI(L,IP,SK1,SK2) * PINORM(L,M)
                DO N = 2, NLAYERS
                  PI_X0P(L,IB,N,SK1,SK2)= PI_X0P(L,IB,1,SK1,SK2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  end loop moments

      ENDDO

!  debug

!      write(77,*)'Fourier',M
!      DO L = M, NMOMENTS
!       write(77,*)'Moment ',L
!       DO SK1 = 1, 4
!        WRITE(77,'(I3,1p4e15.6)')SK1,(PI_XUP(L,4,SK1,SK2),SK2=1,4)
!       ENDDO
!      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_PIMATRIX_SETUP

!

      SUBROUTINE VLIDORT_PIMATRIX_SETUP_OMP ( &
        DO_REFRACTIVE_GEOMETRY, DO_USER_STREAMS, FOURIER,             & ! Input flags and Fourier
        NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS,           & ! Input numbers
        QUAD_STREAMS, USER_STREAMS, COS_SZANGLES, SUN_SZA_COSINES,    & ! Input streams and angles
        NMOMENTS, MUELLER_INDEX, DMAT,                                & ! auxiliary inputs
        PIMM_11, PIMM_KM, PI_XQP, PI_XQM, PI_XUP, PI_XUM, PI_X0P,     & ! Output Pi-Matrix stuff
        PI_XQM_POST, PI_XQM_PRE, PI_XQP_PRE, PI_XUM_POST, PI_XUP_PRE )  ! Output Pi-Matrix stuff

!  Notes
!  -----

!  This is equivalent to the Legendre setup modules in LIDORT

!---------old comments ---------------------------------
!  This needs modification for the case of Refractive atmosphere,
!  because the solar zenith angle is not constant.
!-------------------------------------------------------

!  PI-matrix setup, following recipe of Siewert (1982).
!  Single Fourier component only.

!  Tested against benchmark results in Vestrucci & Siewert (1984), Problem

!  original coding, September 2002, R. Spurr SAO
!  Multibeam SZA coding, July 2004, R. Spurr SAO

!  Coding for refractive geometry case, R. Spurr, RT Solutions, May 2005
!    -- Currently not enabled (2.7/2.8) --> SUN_SZA_COSINES, MUELLER_INDEX are not used

!  Adjusted for Version 2.7 (OpenMP). July 2014.
!  Streamlined for Version 2.8, 05 July 2016.

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, MAXBEAMS, &
                                 MAXMOMENTS, MAX_SZANGLES, MAX_ALLSTRMS_P1, ZERO, HALF, ONE, TWO, THREE

      IMPLICIT NONE

!  Flags and Fourier


      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      INTEGER, INTENT (IN) ::           FOURIER

!  Numbers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS

      INTEGER, INTENT (IN) ::           NMOMENTS

!  Indices and Dmatrix from VLIDORT_DERIVE_INPUTS

      INTEGER         , INTENT (IN)  :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN)  :: DMAT          ( MAXSTOKES, MAXSTOKES )

!  Angles and streams

      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES    ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS    ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS    ( MAX_USER_STREAMS )

!  InOut variables, necessary for the OpenMP thread-safety

      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_KM ( MAX_ALLSTRMS_P1 )

!  Output (generalized spherical functions)

      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM_POST ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP_PRE  ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Local matrices

      DOUBLE PRECISION :: XLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: YLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: ZLM ( 0:MAXMOMENTS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI  ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PHI ( 0:MAXMOMENTS, 0:MAXMOMENTS )
      DOUBLE PRECISION :: PINORM ( 0:MAXMOMENTS, 0:MAXMOMENTS )

      DOUBLE PRECISION, SAVE :: DF_L ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_2LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LSQM4 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_RT_LP3XLM1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: ZHELP ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: RT_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: UPXSQ_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PLEG20 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: M2X_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: XUMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X_RT_UMXSQ ( MAX_ALLSTRMS_P1 )
!mick fix 7/28/2014 - these two defined InOut for now to make VLIDORT threadsafe
      !DOUBLE PRECISION :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      !DOUBLE PRECISION :: PIMM_KM ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X2LP1 ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PISIGN ( 0:MAXMOMENTS, 0:MAXMOMENTS )

!  local variables

!     INTEGER          :: LOCAL_NSTOKES
      DOUBLE PRECISION :: FAC, XSQ, DF_LPM, DF_LP1MM, ZLM_VALUE
      DOUBLE PRECISION :: HX, HY, HZ, H1, H2, RCONST_20, RCONST_21
      DOUBLE PRECISION :: PL20, RL20, PL21, RL21, TL21, XM(MAX_ALLSTRMS_P1)
      INTEGER          :: M, I, I1, L, SK1, SK2, SK3, SK4, N
      INTEGER          :: NS, NSA, M1, IB, IP

!  Set integer M = Fourier number

      M = FOURIER
      NSA = 0

!  Local NSTOKES

!      LOCAL_NSTOKES = NSTOKES
!      IF ( NSTOKES .GT. 1 ) LOCAL_NSTOKES = 4

!  total number of angles

      IF ( DO_USER_STREAMS ) THEN
        NS = NSTREAMS + N_USER_STREAMS
      ELSE
        NS = NSTREAMS
      ENDIF

!  add solar streams, depends on the use of refractive geometry

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        NSA = NS + NBEAMS*NLAYERS
      ELSE
        NSA = NS + NBEAMS
      ENDIF

!  constants

      RCONST_20 = 0.25D0 * SQRT(6.0D0)
      RCONST_21 = 2.0D0 * RCONST_20

!  Coefficient matrix PINORM (for normalization) = [ (L-m)!/(L+m)! ] ^1/
!  ---------------------------------------------------------------------

!  .. first entry = 1 / (2m)!    ---- be careful with overflow

      PHI(M,M) = ONE
      FAC = ONE
      DO L = 1, 2*M
        FAC = FAC * DBLE(L)
      ENDDO
      PHI(M,M) = PHI(M,M) / FAC

!  .. Other entries by recurrence

      DO L = M + 1, NMOMENTS
        FAC = DBLE(L-M)/DBLE(L+M)
        PHI(L,M) = PHI(L-1,M)*FAC
      ENDDO

!  .. Square root

      DO L = M , NMOMENTS
        PINORM(L,M) = SQRT(PHI(L,M))
      ENDDO

!  Additional saved quantities
!  ---------------------------

!  Only for the first fundamental harmonic

!mick fix 7/28/2014 - do for each fourier right now to make VLIDORT threadsafe
      IF ( M .EQ. 0 ) THEN

!  Sign matrix

        DO M1 = 0, NMOMENTS
          DO L = M1, NMOMENTS
            IF (MOD((L-M1),2).EQ.0) THEN
              PISIGN(L,M1) = ONE
            ELSE
              PISIGN(L,M1) = -ONE
            ENDIF
          ENDDO
        ENDDO

!  floating point integer values

        DO L = 0, NMOMENTS
          DF_L(L)    = DBLE(L)
          DF_LP1(L)  = DBLE(L+1)
          DF_2LP1(L) = DBLE(2*L+1)
        ENDDO

        DF_LSQM4(2) = ZERO
        DO L = 3, NMOMENTS
          DF_LSQM4(L) = SQRT(DBLE(L*L-4))
        ENDDO

        DO L = 2, NMOMENTS
          DF_RT_LP3XLM1(L) = SQRT(DBLE((L+3)*(L-1)))
          ZHELP(L) = TWO*DF_2LP1(L)/DF_LP1(L)/DF_L(L)
        ENDDO

      ENDIF

!  local array of all streams (every Fourier component)

      DO I = 1, NSTREAMS
        XM(I) = QUAD_STREAMS(I)
      ENDDO
      IF ( DO_USER_STREAMS ) THEN
        DO I = 1, N_USER_STREAMS
          XM(I+NSTREAMS) = USER_STREAMS(I)
        ENDDO
      ENDIF

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        DO N = 1, NLAYERS
!          DO IB = 1, NBEAMS
!            IP = NS + NBEAMS * (N-1) + IB
!            XM(IP) = SUN_SZA_COSINES(N,IB)
!          ENDDO
!        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          IP = IB + NS
          XM(IP) = COS_SZANGLES(IB)
        ENDDO
      ENDIF

!  factors associated with stream values
!   Special case when XM(I) = ONE, factors are zeroed, and
!    those that give singularities are avoided later (see below)
!   R. Spurr and V. Natraj, 16 january 2006

!mick fix 7/28/2014 - do for each fourier right now to make VLIDORT threadsafe
      IF ( M .EQ. 0 )  THEN
        DO I = 1, NSA
          XSQ = XM(I) * XM(I)
          PLEG20(I)        = HALF*(THREE*XSQ-ONE)
          UMXSQ(I)         = ONE - XSQ
          XUMXSQ(I)        = UMXSQ(I) * XM(I)
          IF ( XM(I) .EQ. ONE ) THEN
           RT_UMXSQ(I)      = ZERO
           X_RT_UMXSQ(I)    = ZERO
           UPXSQ_D_UMXSQ(I) = ZERO
           M2X_D_UMXSQ(I)   = ZERO
          ELSE
           RT_UMXSQ(I)      = SQRT(UMXSQ(I))
           X_RT_UMXSQ(I)    = RT_UMXSQ(I) * XM(I)
           UPXSQ_D_UMXSQ(I) = (ONE + XSQ) / UMXSQ(I)
           M2X_D_UMXSQ(I)   = - TWO * XM(I) / UMXSQ(I)
          ENDIF
          DO L = 0, NMOMENTS
            X2LP1(L,I) = XM(I) * DF_2LP1(L)
          ENDDO
        ENDDO

!  D-matrix and Mueller indices. NOW DONE in VLIDORT_DERIVE_INPUTS
!        DO SK1 = 1, MAXSTOKES
!          DO SK2 = 1, MAXSTOKES
!            MUELLER_INDEX(SK1,SK2) = MAXSTOKES*(SK1-1) + SK2
!            DMAT(SK1,SK2) = ZERO
!          ENDDO
!          IF ( SK1.GT.2) THEN
!            DMAT(SK1,SK1) = -ONE
!          ELSE
!            DMAT(SK1,SK1) = ONE
!          ENDIF
!          MUELLER_DIAGONAL_INDEX(SK1) = MUELLER_INDEX(SK1,SK1)
!        ENDDO

      ENDIF

!  XYZ matrices
!  ------------

!  Inverse XLM diagonal matrices (Siewert (1982), Eq. 35a)
!  YLM diagonal matrices (Siewert (1982), Eq. 35b)
!  ZLM matrices. (Siewert (1982), Eq. 35c)

!  .. for the azimuth-independent harmonic

      IF ( M .EQ. 0 ) THEN

        DO L = 2, NMOMENTS
          XLM_DIAG(L,1) = DF_LP1(L)
          XLM_DIAG(L,2) = DF_RT_LP3XLM1(L)
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
          DO SK1 = 1, NSTOKES
            XLM_DIAG(L,SK1) = ONE/XLM_DIAG(L,SK1)
          ENDDO
        ENDDO

        DO L = 2, NMOMENTS
          YLM_DIAG(L,1) = DF_L(L)
          YLM_DIAG(L,2) = DF_LSQM4(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

!  .. for the other harmonics

      ELSE

        DO L = M, NMOMENTS
          DF_LP1MM = DBLE ( L + 1 - M )
          XLM_DIAG(L,1) = DF_LP1MM
          XLM_DIAG(L,1) = ONE/XLM_DIAG(L,1)
          IF ( L .EQ. 1 ) THEN
            XLM_DIAG(L,2) = ZERO
          ELSE
            XLM_DIAG(L,2) = DF_RT_LP3XLM1(L) * DF_LP1MM / DF_LP1(L)
            XLM_DIAG(L,2) = ONE / XLM_DIAG(L,2)
          ENDIF
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
        ENDDO

        DO SK1 = 1, NSTOKES
          YLM_DIAG(M,SK1) = ZERO
        ENDDO
        DO L = M + 1, NMOMENTS
          DF_LPM = DBLE ( L + M )
          YLM_DIAG(L,1) = DF_LPM
          YLM_DIAG(L,2) = DF_LPM * DF_LSQM4(L) / DF_L(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

        DO L = 2, NMOMENTS
          ZLM_VALUE = DBLE(M) * ZHELP(L)
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              ZLM(L,SK1,SK2) = ZERO
            ENDDO
          ENDDO
          ZLM(L,2,3) = ZLM_VALUE
          ZLM(L,3,2) = ZLM_VALUE
        ENDDO

      ENDIF

!  PI calculation
!  --------------

!  Initialise all the PI matrices for given harmonic

      DO I = 1, NSA
        DO L = M, NMOMENTS
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI(L,I,SK1,SK2) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!mick - alternate PI initialization
      PI = ZERO

!  M = 0 component. [Siewert (1982), Eqs (28) to (31)]

      IF ( M .EQ. 0 ) THEN

!  .. L = 0
        DO I = 1, NSA
          PI(0,I,1,1) = PI(0,I,1,1) + ONE
          PI(0,I,4,4) = PI(0,I,4,4) + ONE
        ENDDO
!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + XM(I)
          PI(1,I,4,4) = PI(1,I,4,4) + XM(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL20 = PLEG20(I)
          RL20 = RCONST_20*UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL20
          PI(2,I,2,2) = PI(2,I,2,2) + RL20
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                PI(L+1,I,SK1,SK2) = ( HX - HY ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M = 1 component. [Siewert (1982), Eqs (32) to (34), plus (35)]

      ELSE IF ( M .EQ. 1 ) THEN

!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + RT_UMXSQ(I)
          PI(1,I,4,4) = PI(1,I,4,4) + RT_UMXSQ(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL21 =   THREE     * X_RT_UMXSQ(I)
          RL21 = - RCONST_21 * X_RT_UMXSQ(I)
          TL21 = + RCONST_21 * RT_UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL21
          PI(2,I,2,2) = PI(2,I,2,2) + RL21
          PI(2,I,2,3) = PI(2,I,2,3) + TL21
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
          PI(2,I,3,2) = PI(2,I,2,3)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M > 1 components. [Siewert (1982), Eqs (36) to (38), plus (35)]
!  Limiting case of XM(I) = ONE requires special treatment to avoid NaN
!   R. Spurr and V. Natraj, 16 january 2006

      ELSE

!  .. L = M
        IF ( M .EQ. 2 ) THEN
          DO I = 1, NSA
            PIMM_11(I) =   THREE * UMXSQ(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) =   RCONST_21
            ELSE
              PIMM_KM(I) =   RCONST_21 * UMXSQ(I)
            ENDIF
          ENDDO
        ELSE
          H1 = DF_2LP1(M-1)
          H2 = DF_LP1(M-1)*H1/DF_RT_LP3XLM1(M-1)
          DO I = 1, NSA
            PIMM_11(I) = H1 * RT_UMXSQ(I) * PIMM_11(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) = ZERO
            ELSE
              PIMM_KM(I) = H2 * RT_UMXSQ(I) * PIMM_KM(I)
            ENDIF
          ENDDO
        ENDIF

        DO I = 1, NSA
          PI(M,I,1,1) = PI(M,I,1,1) + PIMM_11(I)
          IF ( XM(I).EQ.ONE ) THEN
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * TWO
            PI(M,I,2,3) = PI(M,I,2,3) - PIMM_KM(I) * TWO
          ELSE
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * UPXSQ_D_UMXSQ(I)
            PI(M,I,2,3) = PI(M,I,2,3) + PIMM_KM(I) * M2X_D_UMXSQ(I)
          ENDIF
          PI(M,I,3,3) = PI(M,I,2,2)
          PI(M,I,4,4) = PI(M,I,1,1)
          PI(M,I,3,2) = PI(M,I,2,3)
        ENDDO
!  .. L = M + 1
        IF ( M .LT. NMOMENTS ) THEN
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(M,I)*PI(M,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(M,SK1,SK3)*PI(M,I,SK3,SK2)
                ENDDO
                PI(M+1,I,SK1,SK2) = ( HX + HZ ) * XLM_DIAG(M,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!  .. L > M + 1
        DO L = M + 1, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Normalized output.
!  ------------------

      DO L = M, NMOMENTS

!  .. at quadrature streams

        DO I = 1, NSTREAMS

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI_XQP(L,I,SK1,SK2) = PI(L,I,SK1,SK2) * PINORM(L,M)
            ENDDO
          ENDDO

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES

              H1 = ZERO
              H2 = ZERO
              DO SK3 = 1, NSTOKES
                H1 = H1 + DMAT(SK1,SK3) * PI_XQP(L,I,SK3,SK2)
                H2 = H2 + PI_XQP(L,I,SK1,SK3) * DMAT(SK3,SK2)
              ENDDO

              PI_XQP_PRE(L,I,SK1,SK2)  = H1
              PI_XQM_PRE(L,I,SK1,SK2)  = H1 * PISIGN(L,M)
              PI_XQM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

              H1 = ZERO
              DO SK3 = 1, NSTOKES
                H2 = ZERO
                DO SK4 = 1, NSTOKES
                  H2 = H2 + PI_XQP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                ENDDO
                H1 = H1 + DMAT(SK1,SK3)*H2
              ENDDO
              PI_XQM(L,I,SK1,SK2)  = H1 * PISIGN(L,M)

            ENDDO
          ENDDO

        ENDDO

!  .. at positive user_defined angles

        IF ( DO_USER_STREAMS ) THEN

          DO I = 1, N_USER_STREAMS
            I1 = I + NSTREAMS

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_XUP(L,I,SK1,SK2) = PI(L,I1,SK1,SK2) * PINORM(L,M)
              ENDDO
            ENDDO

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES

                H2 = ZERO
                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H1 = H1 + DMAT(SK1,SK3) * PI_XUP(L,I,SK3,SK2)
                  H2 = H2 + DMAT(SK3,SK2) * PI_XUP(L,I,SK1,SK3)
                ENDDO
                PI_XUP_PRE(L,I,SK1,SK2)  = H1
                PI_XUM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H2 = ZERO
                  DO SK4 = 1, NSTOKES
                    H2 = H2 + PI_XUP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                  ENDDO
                  H1 = H1 + DMAT(SK1,SK3)*H2
                ENDDO
                PI_XUM(L,I,SK1,SK2) = H1 * PISIGN(L,M)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

!  .. at solar zenith angles

!  depends on use of refractive geometry

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!          DO N = 1, NLAYERS
!            DO IB = 1, NBEAMS
!              IP = NS + NBEAMS * (N-1) + IB
!              DO SK1 = 1, NSTOKES
!                DO SK2 = 1, NSTOKES
!                  PI_X0P(L,IB,N,SK1,SK2) =
!     &                 PI(L,IP,SK1,SK2) * PINORM(L,M)
!                ENDDO
!              ENDDO
!            ENDDO
!          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            IP = IB + NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_X0P(L,IB,1,SK1,SK2) = &
                     PI(L,IP,SK1,SK2) * PINORM(L,M)
                DO N = 2, NLAYERS
                  PI_X0P(L,IB,N,SK1,SK2)= PI_X0P(L,IB,1,SK1,SK2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  end loop moments

      ENDDO

!  debug

!      write(77,*)'Fourier',M
!      DO L = M, NMOMENTS
!       write(77,*)'Moment ',L
!       DO SK1 = 1, 4
!        WRITE(77,'(I3,1p4e15.6)')SK1,(PI_XUP(L,4,SK1,SK2),SK2=1,4)
!       ENDDO
!      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_PIMATRIX_SETUP_OMP

!

      SUBROUTINE EMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING,                                  & ! Input flags
        NLAYERS, NBEAMS, N_USER_STREAMS, TAYLOR_ORDER, N_PARTLAYERS, & ! Input numbers
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input Level control
        USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                      & ! Input optical, streams   
        BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,     & ! Input solar beam
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,      & ! Input User-stream Trans.
        EMULT_HOPRULE, SIGMA_M, SIGMA_P,                             & ! Output Multipliers
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                 ! Output Multipliers

!  Prepare multipliers for the Beam source terms

!  Version 2p7 notes
!     Rob Fix 05/10/13  - Introduce Taylor series routines (first version) HOPITAL_TOLERANCE replaced by TAYLOR_SMALL.
!     Rob Fix 02/19/14  - Small numbers analysis finalized using Taylor_Order parameter, as for LIDORT

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, ZERO, ONE, TAYLOR_SMALL

!  1/31/21. Version 2.8.3. Needed here

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_1

      IMPLICIT NONE

!  Inputs
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING

!  Numbers

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      INTEGER, INTENT (IN) ::           N_PARTLAYERS

!  Level control

      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

!  Streams, optical

      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT ( MAX_PARTLAYERS )

!  Solar beam

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  user-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  OUTPUT
!  ======

!  Diagnostics

      LOGICAL, INTENT (OUT) ::          EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Multipliers

      DOUBLE PRECISION, INTENT (OUT) :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER          :: N, UT, UM, IB
      DOUBLE PRECISION :: WDEL, WX, WUDEL, UDEL, EPS
      DOUBLE PRECISION :: DIFF, SB, SECMUM, SU, SD
      DOUBLE PRECISION :: UX_DN, UX_UP, WDEL_UXUP

!mick fix 7/23/2014 - EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN initialized for packing
!mick fix 1/5/2021  - SIGMA_M, SIGMA_P initialized for packing

      SIGMA_M  = ZERO
      SIGMA_P  = ZERO
      EMULT_UP = ZERO
      EMULT_DN = ZERO
      UT_EMULT_UP = ZERO
      UT_EMULT_DN = ZERO

!  Taylor-series flags for Downwelling EMULT
!  -----------------------------------------

      IF ( DO_DNWELLING ) THEN
       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_DN(N) ) THEN
         DO IB = 1, NBEAMS
          SB = AVERAGE_SECANT(N,IB)
          DO UM = 1, N_USER_STREAMS
            DIFF = DABS ( USER_SECANTS(UM) - SB )
            EMULT_HOPRULE(N,UM,IB) = ( DIFF .LT. TAYLOR_SMALL )
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDIF

!  sigma functions (all layers)
!  ----------------------------

      DO N = 1, NLAYERS
       DO IB = 1, NBEAMS
        SB = AVERAGE_SECANT(N,IB)
        DO UM = 1, N_USER_STREAMS
          SECMUM = USER_SECANTS(UM)
          SIGMA_P(N,UM,IB) = SB + SECMUM
          SIGMA_M(N,UM,IB) = SB - SECMUM
        ENDDO
       ENDDO
      ENDDO

!  upwelling External source function multipliers
!  ----------------------------------------------

      IF ( DO_UPWELLING ) THEN

!  whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_UP(UM,N,IB) = ZERO
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                WUDEL = WDEL * T_DELT_USERM(N,UM)
                SU = ( ONE - WUDEL ) / SIGMA_P(N,UM,IB)
                EMULT_UP(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SU
              ENDDO
            ENDIF
          ENDDO
         ENDIF
        ENDDO

!  Partial layer

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
            DO UM = 1, N_USER_STREAMS
              UT_EMULT_UP(UM,UT,IB) = ZERO
            ENDDO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            WDEL = T_DELT_MUBAR(N,IB)
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              WDEL_UXUP = UX_UP * WDEL
              SU = ( WX - WDEL_UXUP ) / SIGMA_P(N,UM,IB)
              UT_EMULT_UP(UM,UT,IB) = ITRANS_USERM(N,UM,IB) * SU
            ENDDO
          ENDIF
         ENDDO
        ENDDO

      ENDIF

!  downwelling External source function multipliers
!  ------------------------------------------------

      IF ( DO_DNWELLING ) THEN

!  whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_DN(UM,N,IB) = ZERO
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                UDEL = T_DELT_USERM(N,UM)
                IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                  EPS = SIGMA_M(N,UM,IB)
                  CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAU_VERT(N), WDEL, ONE, SD )
! 5/10/13         CALL TAYLOR_SERIES_1 ( -EPS, DELTAU_VERT(N), ONE, ONE, UDEL, SD )
! zero order       SD = DELTAU_VERT(N) * UDEL                          ! Old code
! first order      SD = SD*(ONE-HALF*DELTAU_VERT(N)*SIGMA_M(N,UM,IB))  ! Old code
                ELSE
                  SD = ( UDEL - WDEL ) / SIGMA_M(N,UM,IB)
                ENDIF
                EMULT_DN(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SD
              ENDDO
            ENDIF
          ENDDO
         ENDIF
        ENDDO

!  Partial layer

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
            DO UM = 1, N_USER_STREAMS
              UT_EMULT_DN(UM,UT,IB) = ZERO
            ENDDO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            DO UM = 1, N_USER_STREAMS
              UX_DN = T_UTDN_USERM(UT,UM)
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                EPS = SIGMA_M(N,UM,IB)
                CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAU_VERT(UT), WX, ONE, SD )
! 5/10/13          CALL TAYLOR_SERIES_1 ( -EPS, PARTAU_VERT(UT), ONE, ONE, UX_DN, SD )
! zero order       SD = PARTAU_VERT(UT) * UX_DN          ! Old code
              ELSE
                SD = ( UX_DN - WX ) / SIGMA_M(N,UM,IB)
              ENDIF
              UT_EMULT_DN(UM,UT,IB) = ITRANS_USERM(N,UM,IB) * SD
            ENDDO
          ENDIF
         ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE EMULT_MASTER

!

      SUBROUTINE EMULT_MASTER_OBSGEO ( &
        DO_UPWELLING, DO_DNWELLING,                                  & ! Input flags
        NLAYERS, NBEAMS, TAYLOR_ORDER, N_PARTLAYERS,                 & ! Input numbers
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input Level control
        USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                      & ! Input optical, streams   
        BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,     & ! Input solar beam
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,      & ! Input User-stream Trans.
        EMULT_HOPRULE, SIGMA_M, SIGMA_P,                             & ! Output Multipliers
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                 ! Output Multipliers

!  Prepare multipliers for the Beam source terms

!  Version 2p7 notes
!     Rob Fix 05/10/13  - Introduce Taylor series routines (first version) HOPITAL_TOLERANCE replaced by TAYLOR_SMALL.
!     Rob Fix 02/19/14  - Small numbers analysis finalized using Taylor_Order parameter, as for LIDORT

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, ZERO, ONE, TAYLOR_SMALL

!  1/31/21. Version 2.8.3. Needed here

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_1

      IMPLICIT NONE

!  Inputs
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING

!  Numbers


      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      INTEGER, INTENT (IN) ::           N_PARTLAYERS

!  Level control

      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

!  Streams, optical

      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT ( MAX_PARTLAYERS )

!  Solar beam

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  user-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  OUTPUT
!  ======

!  Diagnostics

      LOGICAL, INTENT (OUT) ::          EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Multipliers

      DOUBLE PRECISION, INTENT (OUT) :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER          :: N, UT, IB, LUM
      DOUBLE PRECISION :: WDEL, WX, WUDEL, UDEL, EPS
      DOUBLE PRECISION :: DIFF, SB, SECMUM, SU, SD
      DOUBLE PRECISION :: UX_DN, UX_UP, WDEL_UXUP

!  Local user index

      LUM = 1

!mick fix 7/23/2014 - EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN initialized for packing
!mick fix 1/5/2021  - SIGMA_M, SIGMA_P initialized for packing

      SIGMA_M  = ZERO
      SIGMA_P  = ZERO
      EMULT_UP = ZERO
      EMULT_DN = ZERO
      UT_EMULT_UP = ZERO
      UT_EMULT_DN = ZERO

!  Taylor-series flags for Downwelling EMULT
!  -----------------------------------------

      IF ( DO_DNWELLING ) THEN
       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_DN(N) ) THEN
         DO IB = 1, NBEAMS
           SB = AVERAGE_SECANT(N,IB)
           DIFF = DABS ( USER_SECANTS(IB) - SB )
           EMULT_HOPRULE(N,LUM,IB) = ( DIFF .LT. TAYLOR_SMALL )
         ENDDO
        ENDIF
       ENDDO
      ENDIF

!  sigma functions (all layers)
!  ----------------------------

      DO N = 1, NLAYERS
       DO IB = 1, NBEAMS
         SB = AVERAGE_SECANT(N,IB)
         SECMUM = USER_SECANTS(IB)
         SIGMA_P(N,LUM,IB) = SB + SECMUM
         SIGMA_M(N,LUM,IB) = SB - SECMUM
       ENDDO
      ENDDO

!  upwelling External source function multipliers
!  ----------------------------------------------

      IF ( DO_UPWELLING ) THEN

!  whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
              EMULT_UP(LUM,N,IB) = ZERO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              WUDEL = WDEL * T_DELT_USERM(N,IB)
              SU = ( ONE - WUDEL ) / SIGMA_P(N,LUM,IB)
              EMULT_UP(LUM,N,IB) = ITRANS_USERM(N,LUM,IB) * SU
            ENDIF
          ENDDO
         ENDIF
        ENDDO

!  Partial layer

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
            UT_EMULT_UP(LUM,UT,IB) = ZERO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            WDEL = T_DELT_MUBAR(N,IB)
            UX_UP = T_UTUP_USERM(UT,IB)
            WDEL_UXUP = UX_UP * WDEL
            SU = ( WX - WDEL_UXUP ) / SIGMA_P(N,LUM,IB)
            UT_EMULT_UP(LUM,UT,IB) = ITRANS_USERM(N,LUM,IB) * SU
          ENDIF
         ENDDO
        ENDDO

      ENDIF

!  downwelling External source function multipliers
!  ------------------------------------------------

!    .. Note use of L'Hopitals Rule
!       Retaining only the first order term

      IF ( DO_DNWELLING ) THEN

!  whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
              EMULT_DN(LUM,N,IB) = ZERO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              UDEL = T_DELT_USERM(N,IB)
              IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
                EPS = SIGMA_M(N,LUM,IB)
                CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAU_VERT(N), WDEL, ONE, SD )
! 5/10/13        CALL TAYLOR_SERIES_1 ( -EPS, DELTAU_VERT(N), ONE, ONE, UDEL, SD )
! zero order     SD = DELTAU_VERT(N) * UDEL                           ! Old code
! first order    SD = SD*(ONE-HALF*DELTAU_VERT(N)*SIGMA_M(N,LUM,IB))  ! Old code
              ELSE
                SD = ( UDEL - WDEL ) / SIGMA_M(N,LUM,IB)
              ENDIF
              EMULT_DN(LUM,N,IB) = ITRANS_USERM(N,LUM,IB) * SD
            ENDIF
          ENDDO
         ENDIF
        ENDDO

!  Partial layer

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
            UT_EMULT_DN(LUM,UT,IB) = ZERO
          ELSE
            WX    = T_UTDN_MUBAR(UT,IB)
            UX_DN = T_UTDN_USERM(UT,IB)
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              EPS = SIGMA_M(N,LUM,IB)
              CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAU_VERT(UT), WX, ONE, SD )
! 5/10/13     CALL TAYLOR_SERIES_1 ( -EPS, PARTAU_VERT(UT), ONE, ONE, UX_DN, SD )
! zero order   SD = PARTAU_VERT(UT) * UX_DN          ! Old code
            ELSE
              SD = ( UX_DN - WX ) / SIGMA_M(N,LUM,IB)
            ENDIF
            UT_EMULT_DN(LUM,UT,IB) = ITRANS_USERM(N,LUM,IB) * SD
          ENDIF
         ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE EMULT_MASTER_OBSGEO

      END MODULE vlidort_miscsetups_m

