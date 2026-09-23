
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
! #            VLIDORT_CHAPMAN                                  #
! #              LEVELS_GEOMETRY_PREPARE                        #
! #              PARTIALS_GEOMETRY_PREPARE                      #
! #                                                             #
! ###############################################################

      MODULE vlidort_geometry_m

      PRIVATE :: LEVELS_GEOMETRY_PREPARE, PARTIALS_GEOMETRY_PREPARE
      PUBLIC  :: VLIDORT_CHAPMAN

!  Version 2.8, The following geometry routines have been removed.
!                OUTGOING_SPHERGEOM_FINE_UP
!                OUTGOING_SPHERGEOM_FINE_DN
!                MULTI_OUTGOING_ADJUSTGEOM
!                OBSGEOM_OUTGOING_ADJUSTGEOM
!                LOSONLY_OUTGOING_ADJUSTGEOM  !  Version 2.6. Modified geometry routine for thermal emission at night

      CONTAINS

      SUBROUTINE VLIDORT_CHAPMAN ( &
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,            & !  Input
            NLAYERS, NBEAMS, FINEGRID, BEAM_SZAS,                 & !  Input
            N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, & !  Input
            EARTH_RADIUS, RFINDEX_PARAMETER,                      & !  Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,         & !  Input
            CHAPMAN_FACTORS, PARTIAL_CHAPFACS, SZA_LOCAL_INPUT,   & !  Output
            SUN_SZA_COSINES, FAIL, MESSAGE, TRACE )                 !  Output

!  The following options apply:

!   1. If the plane-parallel flag is on, no further inputs are required

!   2. If the plane-parallel flag is off, then Pseudo-spherical:
!       (a) Straight line geometry, must specify
!               Earth_radius, height grid
!       (b) Refractive geometry, must specify
!               Earth_radius, height grid
!               pressure grid, temperature grid

!  The logic will be checked before the module is called.

!  Newly programmed by R. Spurr, RT SOLUTIONS Inc. 5/5/05.

!    Based round a call to a pure geometry module which returns slant
!    path distances which was adapted for use in the Radiant model by
!    R. Spurr during an OCO L2 intensive April 24-29, 2005.

!  Major upgrade to include partial Chapman factors. 1/9/18 for Version 2.8
!    -- Cleaned up code
!    -- Introduced New subroutine for partials.

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, only : MAX_SZANGLES, MAXBEAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                 ZERO, ONE, DEG_TO_RAD

!  implicit none

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  flag for plane parallel case

      LOGICAL  , intent(in) :: DO_PLANE_PARALLEL

!  flag for refractive geometry

      LOGICAL  , intent(in) :: DO_REFRACTIVE_GEOMETRY

!  Number of layers

      INTEGER  , intent(in) :: NLAYERS

!  number of solar beams to be processed

      INTEGER  , intent(in)  :: NBEAMS

!  number of fine layers within coarse layers

      INTEGER  , intent(in)  :: FINEGRID(MAXLAYERS)

!  TOA solar zenith angles

      DOUBLE PRECISION, intent(in)  :: BEAM_SZAS ( MAXBEAMS )

!  New partial layer control, 1/9/18

      INTEGER  , intent(in)         :: N_PARTLAYERS
      INTEGER  , intent(in)         :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      DOUBLE PRECISION, intent(in)  :: PARTLAYERS_VALUES   (MAX_PARTLAYERS)

!  Earth radius (km)

      DOUBLE PRECISION, intent(in)  :: EARTH_RADIUS
        
!  Refractive index parametaer (Born-Wolf approximation)

      DOUBLE PRECISION, intent(in)  :: RFINDEX_PARAMETER

!  Coarse grids of heights, pressures and temperatures

      DOUBLE PRECISION, intent(in)  :: HEIGHT_GRID     (0:MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: PRESSURE_GRID   (0:MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: TEMPERATURE_GRID(0:MAXLAYERS)

!  Output arguments
!  ----------------

!  Chapman factors

      DOUBLE PRECISION, intent(inout) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Partial-layer Chapman Factors, 1/9/18

      DOUBLE PRECISION, intent(inout) :: PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  solar zenith angles at nadir

      DOUBLE PRECISION, intent(inout) :: SZA_LOCAL_INPUT(0:MAXLAYERS,MAXBEAMS)

!  average cosines (for the refractive geometry case)

      DOUBLE PRECISION, intent(inout) :: SUN_SZA_COSINES(MAXLAYERS,MAXBEAMS)

!  output status

      LOGICAL      , intent(out)  :: FAIL
      CHARACTER*(*), intent(out)  :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  number of iterations (refractive case only)
!      This is debug output

      INTEGER      :: ITERSAVE(MAXLAYERS)
      INTEGER      :: PITERSAVE(MAX_PARTLAYERS)

!  local height arrays

      DOUBLE PRECISION    :: RADII(0:MAXLAYERS)
      DOUBLE PRECISION    :: DELZ (MAXLAYERS)

!  other local variables

      INTEGER             :: IBEAM, N
      DOUBLE PRECISION    :: SUN0, TH_TOA, MU_TOA, SN_TOA
      CHARACTER*3         :: C3

!  set ups
!  -------

!  earth radii and heights differences

      RADII(0) = HEIGHT_GRID(0) + EARTH_RADIUS
      DO N = 1, NLAYERS
        RADII(N) = HEIGHT_GRID(N) + EARTH_RADIUS
        DELZ(N)  = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
      ENDDO

!  get spherical optical depths
!  ----------------------------

!  start beam loop

      DO IBEAM = 1, NBEAMS

!  BOA and nadir-TOA values

        SUN0    = BEAM_SZAS(IBEAM)
        TH_TOA  = SUN0 * DEG_TO_RAD
        MU_TOA = COS(TH_TOA)
        SN_TOA = SQRT ( ONE - MU_TOA * MU_TOA )

!  Levels geometry Chapman factor preparation

        CALL LEVELS_GEOMETRY_PREPARE &       
          ( IBEAM, NLAYERS, FINEGRID, SUN0, TH_TOA, MU_TOA, SN_TOA,                  & ! Input
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, RFINDEX_PARAMETER,            & ! Input
            EARTH_RADIUS, RADII, DELZ, HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, & ! Input
            CHAPMAN_FACTORS, SZA_LOCAL_INPUT, SUN_SZA_COSINES,                       & ! Output
            ITERSAVE, FAIL, MESSAGE )                                                  ! Output
        IF ( FAIL ) THEN
          C3 = '000' ; write(C3,'(I3)')IBEAM
          TRACE = 'LEVELS_GEOMETRY_PREPARE subroutine failed for Beam # '//C3
          RETURN
        ENDIF

!  New, 1/9/18 for Version 2.8, Partial-layer Chapman Factors
!  Add this line when you get the refractive geometry working
!            RFINDEX_PARAMETER, FINEGRID, HEIGHTS, PRESSURES, TEMPERATURES,                    & ! Input

        IF ( N_PARTLAYERS .gt. 0 ) then
          CALL PARTIALS_GEOMETRY_PREPARE &       
          ( IBEAM, NLAYERS, SUN0, TH_TOA, MU_TOA, SN_TOA,           & ! Input
            N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & !  Input
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, RADII, DELZ, & ! Input
            PARTIAL_CHAPFACS, PITERSAVE, FAIL, MESSAGE )              ! Output
          IF ( FAIL ) THEN
            C3 = '000' ; write(C3,'(I3)')IBEAM
            TRACE = 'PARTIALS_GEOMETRY_PREPARE subroutine failed for Beam # '//C3
            RETURN
          ENDIF
        ENDIF

!  end beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHAPMAN

!

      SUBROUTINE LEVELS_GEOMETRY_PREPARE & 
          ( IBEAM, NLAYERS, FINEGRID, SZA_GEOM_TRUE, TH_TOA, MU_TOA, SN_TOA,  & ! Input
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, RFINDEX_PARAMETER,     & ! Input
            REARTH, RADII, DELZ, HEIGHTS, PRESSURES, TEMPERATURES,            & ! Input
            CHAPMAN_FACTORS, SZA_LEVEL_OUTPUT, SUN_SZA_COSINES,               & ! Output
            ITERSAVE, FAIL, MESSAGE )                                           ! Output

!  Module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXBEAMS, MAXLAYERS, MAXFINELAYERS, ZERO, HALF, ONE, DEG_TO_RAD

!  Implicit none

      IMPLICIT NONE
        
!  Generate path CHAPMAN_FACTORS and SZA angles SZA_LEVEL_OUTPUT
!  for a curved ray-traced beam through a multilayer atmosphere.

!  Coarse layering is input to the module. Values of Z, P, T  are
!  given at the layer boundaries, with the first value (index 0) at TOA.

!  The refractive geometry is assumed to start at the TOA level.

!  We also require the earth radius and the refractive index parameter
!   (For the Born-Wolf approximation)

!  There is no refraction if the flag DO_REFRACTIVE_GEOMETRY is not set.
!  In this case we do not require pressure and temperature information.
!  The calculation will then be for geometric rays.

!  The plane parallel Flag and the refractive geometry flag should not
!  both be true - this should be checked outside.

!  In the refracting case, fine-gridding of pressure and temperature is
!  done internally, temperature is interpolated linearly with height,
!  and pressure log-linearly. The refraction uses Snell's law rule.
!  Finelayer gridding assumes equidistant heights within coarse layers
!  but the number of fine layers can be varied

!  Output is specified at coarse layer boundaries

!  Module is stand-alone.

!  Reprogrammed for the OCO L2 algorithm
!   R. Spurr, RT Solutions, Inc.   April 27, 2005

!  Intended use in VLIDORT and Radiant RT models.

!  Input arguments
!  ===============

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  number of coarse layers

      INTEGER  , intent(in)  :: NLAYERS

!  number of fine layers within coarse layers

      INTEGER  , intent(in)  :: FINEGRID(MAXLAYERS)

!  True solar zenith angle (degrees). Angle radians, cosine, sine

      DOUBLE PRECISION, intent(in)  :: SZA_GEOM_TRUE, TH_TOA, MU_TOA, SN_TOA

!  Earth radius (km)

      DOUBLE PRECISION, intent(in)  :: REARTH
        
!  local height arrays

      DOUBLE PRECISION, intent(in)  :: RADII(0:MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: DELZ (MAXLAYERS)

!  Refractive index parametaer (Born-Wolf approximation)

      DOUBLE PRECISION, intent(in)  :: RFINDEX_PARAMETER

!  flag for plane parallel case

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL
        
!  flag for refractive geometry

      LOGICAL  , intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Coarse grids of heights, pressures and temperatures

      DOUBLE PRECISION, intent(in)  :: HEIGHTS     (0:MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: PRESSURES   (0:MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: TEMPERATURES(0:MAXLAYERS)

!  Output arguments
!  ================

!  Path segments distances (km)

      DOUBLE PRECISION, intent(inout) :: CHAPMAN_FACTORS(MAXLAYERS,MAXLAYERS,MAXBEAMS)

!  Solar zenith angles at nadir

      DOUBLE PRECISION, intent(inout) :: SZA_LEVEL_OUTPUT(0:MAXLAYERS,MAXBEAMS)

!  Average cosines (for the refractive geometry case)

      DOUBLE PRECISION, intent(inout) :: SUN_SZA_COSINES(MAXLAYERS,MAXBEAMS)

!  Output status (FAIL/MESSAGE) and number of iterations (refractive case, DEBUG only)

      LOGICAL      , intent(out)  :: FAIL
      CHARACTER*(*), intent(out)  :: MESSAGE
      INTEGER      , intent(out)  :: ITERSAVE(MAXLAYERS)

!  Local variables
!  ===============
      
!  fine layer gridding for refraction

      DOUBLE PRECISION    :: ZRFINE(MAXLAYERS,0:MAXFINELAYERS)
      DOUBLE PRECISION    :: PRFINE(MAXLAYERS,0:MAXFINELAYERS)
      DOUBLE PRECISION    :: TRFINE(MAXLAYERS,0:MAXFINELAYERS)

!  help variables

      INTEGER             :: N, J, NRFINE, K, ITER, IB
      LOGICAL             :: LOOP
      DOUBLE PRECISION    :: MU_NEXT, Z1, Z0, Z, T1, T0, T, P1, P0, Q1, Q0, Q, FU, FL
      DOUBLE PRECISION    :: LAYER_DIST, MU_PREV, STH1, SINTH1, STH2, SINTH2, LOCAL_SUBTHICK, &
                             PHI, PHI_0, PHI_CUM, SINPHI, DELPHI, REFRAC, RATIO,              &
                             RE_LOWER, RE_UPPER, DIST, STH2D, SINTH2D, SNELL
      DOUBLE PRECISION    :: MU1, MU2

!  Standard temperature (K) and pressure (mbar).
!  Loschmidt's number (particles/cm2/km).

      DOUBLE PRECISION, PARAMETER  :: T_STANDARD = 273.16D0
      DOUBLE PRECISION, PARAMETER  :: P_STANDARD = 1013.25D0
      DOUBLE PRECISION, PARAMETER  :: STP_RATIO = T_STANDARD / P_STANDARD
      DOUBLE PRECISION, PARAMETER  :: RHO_STANDARD = 2.68675D+24

!  Some setup operations
!  =====================

!  initialise output

      IB = IBEAM
      SZA_LEVEL_OUTPUT(0:NLAYERS,IB) = zero
      DO N = 1, NLAYERS
        CHAPMAN_FACTORS(N,1:N,IB) = Zero
      ENDDO

      FAIL     = .FALSE.
      MESSAGE  = ' '
      ITERSAVE = 0

!  TOA value

      SZA_LEVEL_OUTPUT(0,IB) = SZA_GEOM_TRUE

!  initialize

      STH2D  = Zero
        
!  derive the fine values

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

        Z0 = HEIGHTS(0)
        P0 = PRESSURES(0)
        T0 = TEMPERATURES(0)
        Q0 = LOG(P0)
        DO N = 1, NLAYERS
          NRFINE = FINEGRID(N)
          LOCAL_SUBTHICK = DELZ(N) / DBLE(NRFINE)
          P1 = PRESSURES(N)
          Z1 = HEIGHTS(N)
          T1 = TEMPERATURES(N)
          Q1 = LOG(P1)
          ZRFINE(N,0) = Z0
          PRFINE(N,0) = P0
          TRFINE(N,0) = T0
          DO J = 1, NRFINE - 1
            Z  = Z0 - DBLE(J)*LOCAL_SUBTHICK
            FL = ( Z0 - Z ) / DELZ(N)
            FU = ONE - FL
            Q  = FL * Q1 + FU * Q0
            T  = FL * T0 + FU * T1
            PRFINE(N,J) = EXP ( Q )
            TRFINE(N,J) = T
            ZRFINE(N,J) = Z
          ENDDO
          PRFINE(N,NRFINE) = P1
          TRFINE(N,NRFINE) = T1
          ZRFINE(N,NRFINE) = Z1
          Z0 = Z1
          P0 = P1
          T0 = T1
          Q0 = Q1
        ENDDO
      ENDIF

!  plane-parallel case
!  ===================

      IF ( DO_PLANE_PARALLEL ) THEN
        DO N = 1, NLAYERS
          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE
          DO K = 1, N
            CHAPMAN_FACTORS(N,K,IB) = ONE / MU_TOA
          ENDDO
        ENDDO
        RETURN
      ENDIF
      
!  Refractive Geometry case
!  ========================

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

!  Rob fix 9/9/14. Sun zero case

        IF ( SN_TOA .eq. zero ) then
          DO N = 1, NLAYERS
            CHAPMAN_FACTORS(N,1:N,IB) = one
          ENDDO
          RETURN
        ENDIF

!  Non-zero case....(resume code)

!  Start value of SZA cosine

        MU_PREV = MU_TOA

!  start layer loop

        DO N = 1, NLAYERS

!  start values

          SINTH1 = SN_TOA * RADII(N) / RADII(0)
          STH1   = ASIN(SINTH1)
          PHI_0  = TH_TOA - STH1
          NRFINE = FINEGRID(N)

!  iteration loop

          ITER = 0
          LOOP = .TRUE.
          DO WHILE (LOOP.AND.ITER.LT.100)
            ITER = ITER + 1
            PHI_CUM = Zero
            RE_UPPER = ZRFINE(1,0) + REARTH
            RATIO  = PRFINE(1,0) * STP_RATIO / TRFINE(1,0)
            REFRAC = One + RFINDEX_PARAMETER * RATIO
            SNELL = REFRAC * RE_UPPER * SINTH1
            DO K = 1, N
              LAYER_DIST = Zero
              LOCAL_SUBTHICK = DELZ(K) / DBLE(NRFINE)
              DO J = 0, NRFINE - 1
                RATIO  = PRFINE(K,J) * STP_RATIO / TRFINE(K,J)
                REFRAC = One + RFINDEX_PARAMETER * RATIO
                RE_LOWER = RE_UPPER - LOCAL_SUBTHICK
                SINTH2 = SNELL/ (REFRAC * RE_UPPER )
                IF ( SINTH2.GT.One ) SINTH2 = One
                STH2 = ASIN(SINTH2)
                SINTH2D = RE_UPPER * SINTH2 / RE_LOWER
                IF ( SINTH2D .GT. One) THEN
                  MESSAGE = 'refraction yields angles > 90 some levels'
                  FAIL = .TRUE. ; RETURN
                ENDIF
                STH2D = ASIN(SINTH2D)
                PHI = STH2D - STH2
                SINPHI = SIN(PHI)
                PHI_CUM = PHI_CUM + PHI
                DIST = RE_UPPER * SINPHI / SINTH2D
                LAYER_DIST = LAYER_DIST +  DIST
                RE_UPPER = RE_LOWER
              ENDDO
              CHAPMAN_FACTORS(N,K,IB) = LAYER_DIST / DELZ(K)
            ENDDO

!  examine convergence

            DELPHI = PHI_0 - PHI_CUM
            LOOP = (ABS(DELPHI/PHI_CUM).GT.0.0001D0)

!  Fudge factors to speed up the iteration

            IF ( SZA_GEOM_TRUE .GT. 88.7D0 ) THEN
              STH1 = STH1 + 0.1D0 * DELPHI
              PHI_0 = TH_TOA - STH1
            ELSE IF ( SZA_GEOM_TRUE .LT. 80.0D0 ) THEN
              PHI_0 = PHI_CUM
              STH1 = TH_TOA - PHI_0
            ELSE
              STH1 = STH1 + 0.3D0 * DELPHI
              PHI_0 = TH_TOA - STH1
            ENDIF
            SINTH1 = SIN(STH1)

          ENDDO

!  failure

          IF ( LOOP ) THEN
            MESSAGE = 'refractive iteration not converged'
            FAIL = .TRUE. ; RETURN
          ENDIF

!  Update and save angle output

          MU_NEXT = COS(STH2D)
          MU_PREV = MU_NEXT
          SZA_LEVEL_OUTPUT(N,IB) = ACOS(MU_NEXT) / DEG_TO_RAD
          ITERSAVE(N) = ITER

        ENDDO

!mick fix 9/19/2017 - moved calculation of SUN_SZA_COSINES from VLIDORT_DERIVE_INPUT
!                     to here to resolve an I/O contradiction
!  Set average cosines in the refractive geometry case

        MU1 = COS(SZA_LEVEL_OUTPUT(0,IB)*DEG_TO_RAD)
        DO N = 1, NLAYERS
          MU2 = COS(SZA_LEVEL_OUTPUT(N,IB)*DEG_TO_RAD)
          SUN_SZA_COSINES(N,IB) = HALF * ( MU1 + MU2 )
          MU1 = MU2
        ENDDO

!  Straight line geometry
!  ======================

      ELSE
 
!  Rob fix 9/9/14. Sun zero case

        IF ( SN_TOA .eq. zero ) then
          DO N = 1, NLAYERS
            CHAPMAN_FACTORS(N,1:N,IB) = one
          ENDDO
          RETURN
        ENDIF

!  Non-zero case....(resume code)

        DO N = 1, NLAYERS

!  start values

          SINTH1 = SN_TOA * RADII(N) / RADII(0)
          STH1   = ASIN(SINTH1)
          RE_UPPER = RADII(0)

!  solar zenith angles are all the same = input value

          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE

! loop over layers K from 1 to layer N

          DO K = 1, N

!  sine-rule; PHI = earth-centered angle

            RE_LOWER = RE_UPPER - DELZ(K)
            SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
            STH2   = ASIN(SINTH2)
            PHI    = STH2 - STH1
            SINPHI = SIN(PHI)
            DIST = RE_UPPER * SINPHI / SINTH2
            CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)

!  re-set

            RE_UPPER = RE_LOWER
            SINTH1 = SINTH2
            STH1   = STH2

          ENDDO

!  finish main layer loop

        ENDDO

!  Finish

      ENDIF

!  end of routine

      RETURN
      END SUBROUTINE LEVELS_GEOMETRY_PREPARE

!

      SUBROUTINE PARTIALS_GEOMETRY_PREPARE & 
          ( IBEAM, NLAYERS, SZA_GEOM_TRUE, TH_TOA, MU_TOA, SN_TOA,  & ! Input
            N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & !  Input
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, RADII, DELZ, & ! Input
            PARTIAL_CHAPFACS, PITERSAVE, FAIL, MESSAGE )              ! Output

!  Add this line when you get the refractive geometry working
!            RFINDEX_PARAMETER, FINEGRID, HEIGHTS, PRESSURES, TEMPERATURES,                    & ! Input

!  Module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXLAYERS, MAXBEAMS, MAX_PARTLAYERS, ZERO, ONE, DEG_TO_RAD

!  Implicit none

      IMPLICIT NONE
        
!  Generate path PARTIAL_CHAPFACS
!  for a curved ray-traced beam through a multilayer atmosphere.

!  Coarse layering is input to the module. Values of Z, P, T  are
!  given at the layer boundaries, with the first value (index 0) at TOA.

!  The refractive geometry is assumed to start at the TOA level.
!  We also require the earth radius and the refractive index parameter
!   (For the Born-Wolf approximation)

!  There is no refraction if the flag DO_REFRACTIVE_GEOMETRY is not set.
!  In this case we do not require pressure and temperature information.
!  The calculation will then be for geometric rays.

!  The plane parallel Flag and the refractive geometry flag should not
!  both be true - this should be checked outside.

!  In the refracting case, fine-gridding of pressure and temperature is
!  done internally, temperature is interpolated linearly with height,
!  and pressure log-linearly. The refraction uses Snell's law rule.
!  Finelayer gridding assumes equidistant heights within coarse layers
!  but the number of fine layers can be varied

!  Output is specified at partial-layer boundaries
!    Module is stand-alone.

!  Input arguments
!  ===============

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  number of coarse layers

      INTEGER  , intent(in)  :: NLAYERS

!  New partial layer control, 1/9/18

      INTEGER  , intent(in)         :: N_PARTLAYERS
      INTEGER  , intent(in)         :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      DOUBLE PRECISION, intent(in)  :: PARTLAYERS_VALUES   (MAX_PARTLAYERS)

!  True solar zenith angle (degrees). Angle radians, cosine, sine

      DOUBLE PRECISION, intent(in)  :: SZA_GEOM_TRUE, TH_TOA, MU_TOA, SN_TOA

!  local height arrays

      DOUBLE PRECISION, intent(in)  :: RADII(0:MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: DELZ (MAXLAYERS)

!  flag for plane parallel case

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL
        
!  flag for refractive geometry

      LOGICAL  , intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Add these variables when ready

!  Refractive index parameter (Born-Wolf approximation)
!      DOUBLE PRECISION, intent(in)  :: RFINDEX_PARAMETER
!  number of fine layers within coarse layers
!      INTEGER  , intent(in)  :: FINEGRID(MAXLAYERS)
!  Coarse grids of heights, pressures and temperatures
!      DOUBLE PRECISION, intent(in)  :: HEIGHTS     (0:MAXLAYERS)
!      DOUBLE PRECISION, intent(in)  :: PRESSURES   (0:MAXLAYERS)
!      DOUBLE PRECISION, intent(in)  :: TEMPERATURES(0:MAXLAYERS)

!  Output arguments
!  ================

!  Partial-layer Chapman Factors, 1/9/18

      DOUBLE PRECISION   , intent(inout) :: PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  output status (FAIL/MESSAGE) and number of iterations (refractive case, DEBUG only)

      LOGICAL      , intent(out)  :: FAIL
      CHARACTER*(*), intent(out)  :: MESSAGE
      INTEGER      , intent(out)  :: PITERSAVE(MAX_PARTLAYERS)

!  Local variables
!  ===============
      
!  fine layer gridding for refraction
!      DOUBLE PRECISION    :: ZRFINE(MAXLAYERS,0:MAXFINELAYERS)
!      DOUBLE PRECISION    :: PRFINE(MAXLAYERS,0:MAXFINELAYERS)
!      DOUBLE PRECISION    :: TRFINE(MAXLAYERS,0:MAXFINELAYERS)

!  help variables

      INTEGER             :: N, K, IB, UT
      DOUBLE PRECISION    :: STH1, SINTH1, STH2, SINTH2, PHI, SINPHI,  RE_LOWER, RE_UPPER, DIST, STH2D, HP, RP

!      DOUBLE PRECISION    :: MU_PREV, PHI_0, PHI_CUM, DELPHI, SINTH2D, RATIO
!      INTEGER      :: J, NRFINE, ITER
!      LOGICAL      :: LOOP

!  Standard temperature (K) and pressure (mbar).
!  Loschmidt's number (particles/cm2/km).
!      DOUBLE PRECISION, PARAMETER  :: T_STANDARD = 273.16D0
!      DOUBLE PRECISION, PARAMETER  :: P_STANDARD = 1013.25D0
!      DOUBLE PRECISION, PARAMETER  :: STP_RATIO = T_STANDARD / P_STANDARD
!      DOUBLE PRECISION, PARAMETER  :: RHO_STANDARD = 2.68675D+24

!  Some setup operations
!  =====================

!  initialise output

      IB = IBEAM
      DO UT = 1, N_PARTLAYERS
        DO K = 1, NLAYERS
          PARTIAL_CHAPFACS(UT,K,IB) = Zero
        ENDDO
      ENDDO
      FAIL      = .FALSE.
      MESSAGE   = ' '
      PITERSAVE = 0

!  Temporary Fix - No refractive Geometry

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        MESSAGE = 'Partial Chapman Factors not enabled yet for refractive geometry'
        FAIL = .TRUE. ; RETURN
      ENDIF
  
!  initialize

      STH2D  = 0.0D0
        
!  plane-parallel case
!  ===================

      IF ( DO_PLANE_PARALLEL ) THEN
        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          DO K = 1, N
            PARTIAL_CHAPFACS(UT,K,IB) = ONE / MU_TOA
          ENDDO
        ENDDO
        RETURN
      ENDIF
      
!  Refractive Geometry case
!  ========================

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

!  P L A C E H O L D E R

!  Straight line geometry
!  ======================

      ELSE

 ! Sun zero case

        IF ( SN_TOA .eq. zero ) then
          DO UT = 1, N_PARTLAYERS
            N = PARTLAYERS_LAYERIDX(UT)
            PARTIAL_CHAPFACS(UT,1:N,IB) = ONE
          ENDDO
          RETURN
        ENDIF

!  Non-zero case....(resume code)

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          HP = DELZ(N) * PARTLAYERS_VALUES(UT)
          RP = RADII(N-1) - HP

!  start values

          SINTH1 = SN_TOA * RP / RADII(0)
          STH1   = ASIN(SINTH1)
          RE_UPPER = RADII(0)

! loop over layers K from 1 to layer N - 1
!  sine-rule; PHI = earth-centered angle

          DO K = 1, N-1
            RE_LOWER = RE_UPPER - DELZ(K)
            SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
            STH2   = ASIN(SINTH2)
            PHI    = STH2 - STH1
            SINPHI = SIN(PHI)
            DIST = RE_UPPER * SINPHI / SINTH2
            PARTIAL_CHAPFACS(UT,K,IB) = DIST / DELZ(K)
            RE_UPPER = RE_LOWER
            SINTH1 = SINTH2
            STH1   = STH2
          ENDDO

!  Partial layer

          RE_LOWER = RP
          SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
          STH2   = ASIN(SINTH2)
          PHI    = STH2 - STH1
          SINPHI = SIN(PHI)
          DIST = RE_UPPER * SINPHI / SINTH2
          PARTIAL_CHAPFACS(UT,N,IB) = DIST / HP

!  Check on the values, use sine-rule distancing.
!          SINTH1 = SN_TOA * RP / RADII(0)
!          STH1   = ASIN(SINTH1)
!          DIST = DOT_PRODUCT(PARTIAL_CHAPFACS(UT,1:N-1,IB),DELZ(1:N-1)) + HP*PARTIAL_CHAPFACS(UT,N,IB)
!          write(*,*)'Check ',UT,DIST,SIN(TH_TOA-STH1)*RADII(0)/SN_TOA

!  finish main partial layer loop

        ENDDO

!  Finish

      ENDIF

!  end of routine

      RETURN
      END SUBROUTINE PARTIALS_GEOMETRY_PREPARE

!  Finish

      END MODULE vlidort_geometry_m

