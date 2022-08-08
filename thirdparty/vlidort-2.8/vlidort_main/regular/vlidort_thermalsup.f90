
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
! #   setups                                                    #
! #          THERMAL_SETUP                                      #
! #                                                             #
! #   discrete ordinate particular integral                     #
! #          THERMAL_CLSOLUTION                                 #
! #                                                             #
! #   postprocessing source terms                               #
! #          THERMAL_STERMS_UP                                  #
! #          THERMAL_STERMS_DN                                  #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. No real changes, except..
!   -- Drop LOCAL_UM_START argument where appropriate. User-stream do-loops start with 1 now.

      MODULE vlidort_thermalsup_m

      PUBLIC :: THERMAL_SETUP,      &
                THERMAL_CLSOLUTION, &
                THERMAL_STERMS_UP,  &
                THERMAL_STERMS_DN

      CONTAINS

      SUBROUTINE THERMAL_SETUP ( &
        DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING,                 & ! Flags
        DO_PARTLAYERS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,      & ! Flags
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,     & ! Numbers basic
        THERMAL_BB_INPUT, USER_STREAMS,                              & ! thermal input, streams
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Level control
        OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Input transmittances
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                & ! output thermal setups
        T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )     ! output thermal direct solutions

!  SET-UP OF THERMAL EXPANSION COEFFICIENTS, ALWAYS DONE AFTER DELTA-M.

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_USER_STREAMS, ZERO, ONE

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING

      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS

!  thermal input and streams

      DOUBLE PRECISION, INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  level output control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )

!  Optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )

!  Transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  outputs
!  =======

!  Auxiliary setups

      DOUBLE PRECISION, INTENT (OUT) :: THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  thermal direct solutions

      DOUBLE PRECISION, INTENT (OUT) :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          N, S, UT, UM, NT
      DOUBLE PRECISION :: HELP ( MAXLAYERS )
      DOUBLE PRECISION :: XTAU, SUM, COSMUM
      DOUBLE PRECISION :: OMEGAS1 ( MAXLAYERS )

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_UP ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: T_MULT_DN ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )

!mick fix 7/23/2014 - initialized for packing
      T_DIRECT_UP = ZERO ; T_UT_DIRECT_UP = ZERO
      T_DIRECT_DN = ZERO ; T_UT_DIRECT_DN = ZERO

!  POWERS OF OPTICAL THICKNESS
!  ---------------------------

!  WHOLE LAYER

      DO N = 1, NLAYERS
        DELTAU_POWER(N,1) = ONE
        DO S = 2, N_THERMAL_COEFFS
          DELTAU_POWER(N,S) = DELTAU_VERT(N) * DELTAU_POWER(N,S-1)
        END DO
      END DO

!  PARTIAL LAYER

      IF ( DO_PARTLAYERS ) THEN
        DO UT = 1, N_PARTLAYERS
          XTAU = PARTAU_VERT(UT)
          XTAU_POWER(UT,1) = ONE
          DO S = 2, N_THERMAL_COEFFS
            XTAU_POWER(UT,S) = XTAU * XTAU_POWER(UT,S-1)
          END DO
        END DO
      ENDIF

!  INITIAL SET OF COEFFICIENTS

      DO N = 1, NLAYERS
        THERMCOEFFS(N,1) = THERMAL_BB_INPUT(N-1)
        HELP(N) = (THERMAL_BB_INPUT(N)-THERMAL_BB_INPUT(N-1)) / DELTAU_VERT(N)
      END DO

!  PIECEWISE CONTINUOUS FOR LINEAR REGIME

      IF ( N_THERMAL_COEFFS == 2 ) THEN
        DO N = 1, NLAYERS
          THERMCOEFFS(N,2) = HELP(N)
        END DO
      END IF

!  DERIVATIVE CONTINUITY FOR QUADRATIC REGIME
!    ( FIRST LAYER IS LINEAR; BETTER THAN USING A FREE DERIVATIVE AT TOA
!  IF ( N_THERMAL_COEFFS == 3 ) THEN
!     THERMCOEFFS(1,3) = ZERO
!     THERMCOEFFS(1,2) = HELP(1)
!     DO N = 1, N_COMP_LAYERS - 1
!        N1 = N + 1
!        THERMCOEFFS(N1,2) = THERMCOEFFS(N,2) + TWO * DELTAUS(N) * THERM
!        THERMCOEFFS(N1,3) = ( HELP(N1) - THERMCOEFFS(N1,2) ) / DELTAUS(
!     END DO
!  END IF

!  ALTERNATIVE SCHEME: BACKWARD LOOKING QUADRATICS
!    ( FIRST LAYER IS LINEAR; BETTER THAN USING A FREE DERIVATIVE AT TOA

!      IF ( N_THERMAL_COEFFS == 3 ) THEN
!       THERMCOEFFS(1,3) = ZERO
!       THERMCOEFFS(1,2) = HELP(1)
!       DO N = 1, NLAYERS - 1
!        N1 = N + 1
!        SUM = ( (THERMCOEFFS(N,1)-THERMCOEFFS(N1,1)) / DELTAU_VERT(N) ) + HELP(N1)
!       THERMCOEFFS(N1,3) = SUM / ( DELTAU_VERT(N1) + DELTAU_VERT(N) )
!        THERMCOEFFS(N1,2) = HELP(N1) - DELTAU_VERT(N1)*THERMCOEFFS(N1,3
!       END DO
!      END IF

!  DEBUG CHECK

!      XTAU = ZERO
!      DO N = 1, NLAYERS
!       SUM = DELTAU_VERT(N) / 10.0
!       DO S = 1, 10
!        XTAU = XTAU + SUM
!        WRITE(34,'(1P2E15.6)')  XTAU,THERMCOEFFS(N,1)+THERMCOEFFS(N,2)*S*SUM+THERMCOEFFS(N,3)*S*SUM*S*SUM
!       END DO
!      END DO

!  AUXILIARY QUANTITIES

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO N = 1, NLAYERS
          DO S = 1, N_THERMAL_COEFFS
            TCOM1(N,S) = THERMCOEFFS(N,S)
          END DO
        END DO
      ELSE
        DO N = 1, NLAYERS
          OMEGAS1(N) = ONE - OMEGA_TOTAL(N)
          DO S = 1, N_THERMAL_COEFFS
            TCOM1(N,S) = THERMCOEFFS(N,S) * OMEGAS1(N)
          END DO
        END DO
      ENDIF

! RETURN IF POST-PROCESSING NOT FLAGGED

      IF ( .NOT. DO_USER_STREAMS ) RETURN

!  ZERO DIRECT SOLUTIONS if working in MSMODE only, then return
!   (Zeroing is done at the beginning of the routine now)

      IF ( DO_MSMODE_THERMAL ) RETURN

!  SHORT HAND

      NT = N_THERMAL_COEFFS

!  UPWELLING DIRECT SOLUTION SOURCE TERMS
! ---------------------------------------

      IF ( DO_UPWELLING ) THEN

!  START USER STREAM LOOP

       DO UM = 1, N_USER_STREAMS
        COSMUM = USER_STREAMS(UM)

!  DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          T_MULT_UP(N,NT) = TCOM1(N,NT)
          DO S = NT - 1, 1, -1
           T_MULT_UP(N,S) = TCOM1(N,S) + S * COSMUM * T_MULT_UP(N,S+1)
          ENDDO
          SUM = T_MULT_UP(N,1)
          DO S = 2, NT
           SUM = SUM + T_MULT_UP(N,S) * DELTAU_POWER(N,S)
          ENDDO
          T_MULT_UP(N,0) = - SUM
          T_DIRECT_UP(UM,N) = T_MULT_UP(N,0) * T_DELT_USERM(N,UM) + T_MULT_UP(N,1)
         ENDIF
        ENDDO

!  DIRECT SOLUTION: PARTIAL LAYER SOURCE TERMS

        IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)
          T_UT_DIRECT_UP(UM,UT) = T_MULT_UP(N,0) * T_UTUP_USERM(UT,UM)
          SUM = ZERO
          DO S = 1, NT
           SUM = SUM + T_MULT_UP(N,S) * XTAU_POWER(UT,S)
          END DO
          T_UT_DIRECT_UP(UM,UT) = T_UT_DIRECT_UP(UM,UT) + SUM
         END DO
        ENDIF

!  END USER STREAM LOOP, AND UPWELLING

       ENDDO
      ENDIF

!  DOWNWELLING DIRECT SOLUTION SOURCE TERMS
! -----------------------------------------

      IF ( DO_DNWELLING ) THEN

!  START USER STREAM LOOP

       DO UM = 1, N_USER_STREAMS
        COSMUM = USER_STREAMS(UM)

!  DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          T_MULT_DN(N,NT) = TCOM1(N,NT)
          DO S = NT - 1, 1, -1
           T_MULT_DN(N,S) = TCOM1(N,S) - S * COSMUM * T_MULT_DN(N,S+1)
          END DO
          T_MULT_DN(N,0) = - T_MULT_DN(N,1)
          T_DIRECT_DN(UM,N) = T_MULT_DN(N,0) * T_DELT_USERM(N,UM)
          SUM = ZERO
          DO S = 1, NT
           SUM = SUM + T_MULT_DN(N,S) * DELTAU_POWER(N,S)
          END DO
          T_DIRECT_DN(UM,N) = T_DIRECT_DN(UM,N) + SUM
         END IF
        END DO

!  DIRECT SOLUTION: PARTIAL LAYER SOURCE TERMS

        IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)
          T_UT_DIRECT_DN(UM,UT) = T_MULT_DN(N,0) * T_UTDN_USERM(UT,UM)
          SUM = ZERO
          DO S = 1, NT
           SUM = SUM + T_MULT_DN(N,S) * XTAU_POWER(UT,S)
          END DO
          T_UT_DIRECT_DN(UM,UT) = T_UT_DIRECT_DN(UM,UT) + SUM
         END DO
        END IF

!  END USER STREAM LOOP, AND DOWNWELLING

       ENDDO
      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_SETUP

!

      SUBROUTINE THERMAL_CLSOLUTION ( &
        DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,              & ! input flags
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_THERMAL_TRANSONLY, & ! Input flags
        NSTREAMS, NLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,      & ! Input basic numbers
        NMOMENTS, NSTREAMS_2, N_PARTLAYERS,                       & ! Input numbers
        PARTLAYERS_LAYERIDX, STERM_LMASK_UP, STERM_LMASK_DN,      & ! Input level control
        QUAD_STREAMS, QUAD_HALFWTS, OMEGA_GREEK, SAB, DAB,        & ! Input optical and SAB/DAB
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,           & ! Input Discrete Ord. Trans.
        PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE,                       & ! Input PI matrices, 
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Input thermal setups
        T_WUPPER, T_WLOWER, UT_T_PARTIC,                          & ! Output thermal solutions
        U_TPOS1, U_TNEG1, U_TPOS2, U_TNEG2,                       & ! Output User thermal solutions
        STATUS, MESSAGE, TRACE )                                    ! Exception handling

!  CLASSICAL FUNCTION THERMAL PARTICULAR INTEGRAL, ALL LAYERS.
!  USES COEFFICIENT EXPANSION OF ATTENUATION.

!  1/31/21. Version 2.8.3. No real changes, except..
!   -- Drop LOCAL_UM_START argument where appropriate. User-stream do-loops start with 1 now.

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXMOMENTS,  &
                                 MAX_USER_STREAMS, MAX_THERMAL_COEFFS, MAXSTREAMS_2,            &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS,  ZERO, HALF, TWO

      USE LAPACK_TOOLS_m, Only : DGETRF, DGETRS

      IMPLICIT NONE

!  Flags 
!mick fix 3/30/2015 - added DO_MVOUT_ONLY to input
!mick fix 9/19/2017 - added DO_USER_STREAMS to input

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY 

!  removed 7/6/16 version 2.8
!      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT

      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  basic numbers

      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  other numbers

      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          N_PARTLAYERS

!  Level control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LMASK_DN ( MAXLAYERS )

!  Quadrature and Optical, SAB, DAB

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )

!  discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  PI matrices

      DOUBLE PRECISION, INTENT (IN) :: PI_XQP     ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUP     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  thermal setups

      DOUBLE PRECISION, INTENT (IN) :: THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  OUTPUT
!  ======

!  discrete ordinate thermal solutions

      DOUBLE PRECISION, INTENT (OUT) ::  T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  UT_T_PARTIC( MAXSTREAMS_2, MAX_PARTLAYERS )

!  user-stream thermal solutions

      DOUBLE PRECISION, INTENT (OUT) ::  U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )

!  exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

! LOCAL VARIABLES
! ---------------

      INTEGER ::           AA, AA1, I, I1, J, J1, L, M, O1
      INTEGER ::           S, N, NT, INFO, UM, UT
      DOUBLE PRECISION  :: SUM_M, SUM_P, TK, K1, SD, SU, A5
      DOUBLE PRECISION  :: TERM1, TERM2, SUM1, SUM2
      DOUBLE PRECISION  :: POS1, POS2, NEG1, NEG2
      CHARACTER (LEN=3) :: CI, C3

!  LOCAL MATRICES AND VECTORS

      DOUBLE PRECISION :: T_C_MINUS ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: T_C_PLUS  ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS )

      DOUBLE PRECISION :: HVEC1 ( MAXSTREAMS )
      DOUBLE PRECISION :: HVEC2 ( MAXSTREAMS )
      DOUBLE PRECISION :: JVEC1 ( MAXSTREAMS )
      DOUBLE PRECISION :: TVEC1 ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: TVEC2 ( MAXSTREAMS_2, MAXLAYERS )

      DOUBLE PRECISION :: T_HELP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: T_HELP2 ( 0:MAXMOMENTS )

      DOUBLE PRECISION :: TMAT   ( MAXSTREAMS, MAXSTREAMS )
      INTEGER          :: TPIVOT ( MAXSTREAMS )

! -------------------------------
!  ZERO THE BOUNDARY LAYER VALUES
! -------------------------------

      DO I = 1, NSTREAMS_2
       DO N = 1, NLAYERS
        T_WUPPER(I,N) = ZERO
        T_WLOWER(I,N) = ZERO
       ENDDO
      ENDDO

!  (1,1) COMPONENT ONLY

      O1 = 1

!  INITIALIZE EXCEPTION HANDLING

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  SHORTHAND

      NT = N_THERMAL_COEFFS

!  FOURIER COMPONENT = 0 (ALWAYS)

      M = 0

!  THERMAL TRANSMITTANCE ONLY, QUADRATURE SOLUTIONS
!  ================================================

      IF ( DO_THERMAL_TRANSONLY ) THEN

!  WHOLE LAYER SOLUTIONS

       DO N = 1, NLAYERS
        DO AA = 1, NSTREAMS
         AA1 = AA + NSTREAMS
         K1 = QUAD_STREAMS(AA)
         TK = T_DELT_DISORDS(AA,N)
         T_C_MINUS(AA,N,NT)  = K1 * THERMCOEFFS(N,NT)
         T_C_PLUS(AA,N,NT)   = K1 * THERMCOEFFS(N,NT)
         DO S = NT - 1, 1, -1
          T_C_MINUS(AA,N,S)= K1*(THERMCOEFFS(N,S)-S*T_C_MINUS(AA,N,S+1))
          T_C_PLUS(AA,N,S) = K1*(THERMCOEFFS(N,S)+S*T_C_PLUS (AA,N,S+1))
         END DO
         SUM_P = T_C_PLUS (AA,N,1)
         SUM_M = T_C_MINUS(AA,N,1)
         DO S = 2, NT
          SUM_M = SUM_M + T_C_MINUS(AA,N,S) * DELTAU_POWER(N,S)
          SUM_P = SUM_P + T_C_PLUS(AA,N,S)  * DELTAU_POWER(N,S)
         END DO
         T_C_MINUS(AA,N,0) = - T_C_MINUS(AA,N,1)
         T_C_PLUS(AA,N,0)  = - SUM_P
         T_WLOWER(AA,N)  = TK * T_C_MINUS(AA,N,0) + SUM_M
         T_WUPPER(AA1,N) = TK * T_C_PLUS(AA,N,0)  + T_C_PLUS(AA,N,1)
        END DO
       END DO

! removed Quad_output flag, Version 2.8, 7/6/16
!  OFFGRID: ONLY FOR QUADRATURE OR MEAN-VALUE OUTPUT
!       IF ( DO_QUAD_OUTPUT .OR. DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN

       IF ( DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( N_PARTLAYERS .GT. 0 ) THEN
         DO UT = 1, N_PARTLAYERS

!  REGULAR OFF-GRID SOLUTION

          N  = PARTLAYERS_LAYERIDX(UT)
          DO AA = 1, NSTREAMS
           AA1 = AA + NSTREAMS
           SD = T_C_MINUS(AA,N,0) * T_DISORDS_UTDN(AA,UT)
           SU = T_C_PLUS(AA,N,0)  * T_DISORDS_UTUP(AA,UT)
           DO S = 1, NT
            SD = SD + T_C_MINUS(AA,N,S) * XTAU_POWER(UT,S)
            SU = SU + T_C_PLUS(AA,N,S)  * XTAU_POWER(UT,S)
           END DO
           UT_T_PARTIC(AA,UT)  = SD
           UT_T_PARTIC(AA1,UT) = SU
          END DO

!  END OFFGRID OUTPUT

         ENDDO
        ENDIF
       END IF

!  RETURN THERMAL TRANSMITTANCE-ONLY

       RETURN

!  END THERMAL TRANSMITTANCE-ONLY CLAUSE

      ENDIF

! ---------------------------------------
! CLASSICAL SOLUTIONS FOR ALL LAYERS
! ---------------------------------------

      DO N = 1, NLAYERS

!  SOURCE CONSTANTS

        TERM1 = TWO * TCOM1(N,1)
        TERM2 = TWO * TCOM1(N,2)

!  H SOLUTIONS
!  -----------

!  SOLUTION MATRIX FOR THE REDUCED PROBLEM
!  ( MATRIX SHOULD BE SAVED IN THE LU DECOMPOSITION FORM)

        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            TMAT(I,J) = SAB(I,J,O1,O1,N) * QUAD_STREAMS(I)
          ENDDO
        ENDDO

!  L-U DECOMPOSITION OF THE SOLUTION MATRIX

        CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE = 'DGETRF CALL FOR H, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  H VECTOR_1 AND SOLUTION BY BACK-SUBSTITUTION

        DO I = 1, NSTREAMS
          HVEC1(I) = TERM1
        ENDDO

        CALL DGETRS  ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,HVEC1,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE = 'DGETRS CALL FOR H_1, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  H VECTOR_2 AND SOLUTION BY BACK-SUBSTITUTION

        DO I = 1, NSTREAMS
          HVEC2(I) = TERM2
        ENDDO
        CALL DGETRS ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT, HVEC2,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE = 'DGETRS CALL FOR H_2, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  J SOLUTION
!  ----------

!  SOLUTION MATRIX AND VECTOR FOR THE REDUCED PROBLEM
!  ( MATRIX SHOULD BE SAVED IN THE LU DECOMPOSITION FORM)

        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            TMAT(I,J) = - DAB(I,J,O1,O1,N)
          ENDDO
          JVEC1(I) = HVEC2(I)
        ENDDO

!  L-U DECOMPOSITION OF THE SOLUTION MATRIX

        CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE = 'DGETRF CALL FOR J, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  J VECTOR_1 SOLUTION BY BACK-SUBSTITUTION

        CALL DGETRS ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,JVEC1,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE = 'DGETRS CALL FOR J, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  SET SOLUTION
!  ============

!  EXPANSION

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TVEC1(I,N)  = HALF * (HVEC1(I) + JVEC1(I))
          TVEC1(I1,N) = HALF * (HVEC1(I) - JVEC1(I))
          TVEC2(I,N)  = HALF * HVEC2(I)
          TVEC2(I1,N) = TVEC2(I,N)
        ENDDO

!  VALUES AT THE LAYER BOUNDARIES

        DO I = 1, NSTREAMS_2
          T_WUPPER(I,N) = TVEC1(I,N)
          T_WLOWER(I,N) = TVEC1(I,N) + TVEC2(I,N) * DELTAU_POWER(N,2)
        ENDDO

!mick fix 9/19/2017 - added IF condition
        IF ( DO_USER_STREAMS ) THEN

!  USER SOLUTIONS - HELP VECTORS

          DO L = M, NMOMENTS
            SUM1 = ZERO
            SUM2 = ZERO
            DO  J = 1, NSTREAMS
              A5 = QUAD_HALFWTS(J)
              J1 = J + NSTREAMS
              POS1 = TVEC1(J1,N) * A5 * PI_XQP    (L,J,O1,O1)
              NEG1 = TVEC1(J,N)  * A5 * PI_XQM_PRE(L,J,O1,O1)
              SUM1 = SUM1 + POS1 + NEG1
              POS2 = TVEC2(J,N)  * A5 * PI_XQP    (L,J,O1,O1)
              NEG2 = TVEC2(J1,N) * A5 * PI_XQM_PRE(L,J,O1,O1)
              SUM2 = SUM2 + POS2 + NEG2
            ENDDO
            T_HELP1(L) = SUM1 * OMEGA_GREEK(L,N,O1,O1)
            T_HELP2(L) = SUM2 * OMEGA_GREEK(L,N,O1,O1)
!         if (n.eq.1) write(*,*)'Perturb',L,n,T_HELP1(L)
          ENDDO

!  UPWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

          IF ( DO_UPWELLING.AND.STERM_LMASK_UP(N) ) THEN
            DO UM = 1, N_USER_STREAMS
              POS1 = ZERO
              POS2 = ZERO
              DO L = M, NMOMENTS
                POS1 = POS1 + T_HELP1(L) * PI_XUP(L,UM,O1,O1)
                POS2 = POS2 + T_HELP2(L) * PI_XUP(L,UM,O1,O1)
              ENDDO
              U_TPOS1(UM,N) = POS1
              U_TPOS2(UM,N) = POS2
            ENDDO
          ENDIF

!  DNWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

          IF ( DO_DNWELLING.AND.STERM_LMASK_DN(N) ) THEN
            DO UM = 1, N_USER_STREAMS
              NEG1 = ZERO
              NEG2 = ZERO
              DO L = M, NMOMENTS
                NEG1 = NEG1 + T_HELP1(L)*PI_XUM(L,UM,O1,O1)
                NEG2 = NEG2 + T_HELP2(L)*PI_XUM(L,UM,O1,O1)
              ENDDO
              U_TNEG1(UM,N) = NEG1
              U_TNEG2(UM,N) = NEG2
!         if (n.eq.1) write(*,*)'Perturb',UM,n,U_TPOS1(UM,N),U_TNEG1(UM,N)
            ENDDO
          ENDIF

        ENDIF

!  END LAYER LOOP

      ENDDO

!  OFFGRID: ONLY FOR QUADRATURE OR MEAN-VALUE OUTPUT
!  =================================================

! removed Quad_output flag, Version 2.8, 7/6/16
!       IF ( DO_QUAD_OUTPUT .OR. DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN

      IF ( DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( N_PARTLAYERS .GT. 0 ) THEN
          DO UT = 1, N_PARTLAYERS
            N  = PARTLAYERS_LAYERIDX(UT)
            DO I = 1, NSTREAMS_2
              UT_T_PARTIC(I,UT) = TVEC1(I,N) + TVEC2(I,N) * XTAU_POWER(UT,2)
            ENDDO
          ENDDO
        ENDIF
      END IF

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_CLSOLUTION

!

      SUBROUTINE THERMAL_STERMS_UP ( &
        DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,    & ! Input flags
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input numbers
        N_ALLLAYERS_UP, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP,  & ! Input level control
        USER_STREAMS, T_DELT_USERM, T_UTUP_USERM,                 & ! Input User streams and transmittances
        DELTAU_POWER, XTAU_POWER, U_TPOS1, U_TPOS2,               & ! Input Thermal setups/solutions
        T_DIRECT_UP, T_UT_DIRECT_UP,                              & ! Input thermal direct solutions
        LAYER_TSUP_UP, LAYER_TSUP_UTUP )                            ! Output user RTE thermal

!  THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (UPWELLING)

!  1/31/21. Version 2.8.3. No real changes, except..
!   -- Drop LOCAL_UM_START argument where appropriate. User-stream do-loops start with 1 now.

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_THERMAL_COEFFS, ONE, PI4

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Level control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_UP
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )

!  User Streams/Transm.

      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Thermal setups and user solutions

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER ( MAXLAYERS,      MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )

!  Thermal direct solutions

      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_UP    ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  OUTPUT
!  ======

!  Postprocessed thermal solutions

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UP   ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  LOCAL VARIABLES
!  ----------------

      INTEGER          :: UM, N, UT, S, NT
      DOUBLE PRECISION :: SPAR, COSMUM, FAC, SUM

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_UP ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( GREEN'S FUNCTION SOLUTION )
!  --------------------------------------------------------------------

!  INITIAL MODULUS = 4.PI IF SOLAR SOURCES ARE INCLUDED

      FAC = ONE
      IF ( DO_SOLAR_SOURCES ) FAC = PI4

!  SHORTHAND

      NT = N_THERMAL_COEFFS

!  START USER ANGLE LOOP

      DO UM = 1, N_USER_STREAMS

!  LOCAL COSINE

       COSMUM = USER_STREAMS(UM)

!  DIRECT TERMS TO START

       DO N = N_ALLLAYERS_UP, NLAYERS
         LAYER_TSUP_UP(UM,N) = FAC * T_DIRECT_UP(UM,N)
       ENDDO
       IF ( DO_PARTLAYERS ) THEN
        DO UT = 1, N_PARTLAYERS
         LAYER_TSUP_UTUP(UM,UT) = FAC*T_UT_DIRECT_UP(UM,UT)
        ENDDO
       ENDIF

!  ONLY DO THE NEXT SECTION FOR SCATTERING SOLUTIONS
!   -- Version 2.8, 7/7/16 remove goto 
!       IF ( DO_THERMAL_TRANSONLY ) GO TO 678

       IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          T_MULT_UP(N,2) = U_TPOS2(UM,N)
          T_MULT_UP(N,1) = U_TPOS1(UM,N) + COSMUM * T_MULT_UP(N,2)
          SUM = T_MULT_UP(N,1)
          DO S = 2, NT
           SUM = SUM + T_MULT_UP(N,S) * DELTAU_POWER(N,S)
          ENDDO
          T_MULT_UP(N,0)   = - SUM
          SPAR = T_MULT_UP(N,0) * T_DELT_USERM(N,UM) + T_MULT_UP(N,1)
          LAYER_TSUP_UP(UM,N) = LAYER_TSUP_UP(UM,N) + SPAR * FAC
!     if ( n.eq.nlayers) write(*,*)'Perturb', um, n, layer_tsup_up(um,n)
         ENDIF
        ENDDO

!  PARTIAL LAYER SOURCES - ADD THERMAL DIFFUSE TERM

        IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
           N = PARTLAYERS_LAYERIDX(UT)
           SUM = T_MULT_UP(N,0) * T_UTUP_USERM(UT,UM)
           DO S = 1, NT
            SUM = SUM + T_MULT_UP(N,S) * XTAU_POWER(UT,S)
           END DO
           LAYER_TSUP_UTUP(UM,UT) = LAYER_TSUP_UTUP(UM,UT) + SUM*FAC
!if (ut.eq.1) write(*,*)'Perturb',ut,n,um,layer_tsup_utup(um,ut)
         ENDDO
        ENDIF

!  End scattering clause

       ENDIF

!  CONTINUATION POINT FOR AVOIDING SCATTERING CALCULATIONS
! 678   CONTINUE

!  END USER-STREAM LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_STERMS_UP

!

      SUBROUTINE THERMAL_STERMS_DN ( &
        DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,    & ! Input flags
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input numbers
        N_ALLLAYERS_DN, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_DN,  & ! Input level control
        USER_STREAMS, T_DELT_USERM, T_UTDN_USERM,                 & ! Input User streams and transmittances
        DELTAU_POWER, XTAU_POWER, U_TNEG1, U_TNEG2,               & ! Input Thermal setups/solutions
        T_DIRECT_DN, T_UT_DIRECT_DN,                              & ! Input thermal direct solutions
        LAYER_TSUP_DN, LAYER_TSUP_UTDN )                            ! Output user RTE thermal

!  THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (DOWNWELLING)

!  1/31/21. Version 2.8.3. No real changes, except..
!   -- Drop LOCAL_UM_START argument where appropriate. User-stream do-loops start with 1 now.

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_THERMAL_COEFFS, ONE, PI4

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Level control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_DN
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )

!  User Streams/Transm.

      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Thermal setups and user solutions

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER ( MAXLAYERS,      MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )

!  Thermal direct solutions

      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_DN    ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  OUTPUT
!  ======

!  Postprocessed thermal solutions

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_DN   ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  LOCAL VARIABLES
!  ----------------

      INTEGER          :: UM, N, UT, S, NT
      DOUBLE PRECISION :: SPAR, COSMUM, FAC, SUM

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_DN ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )

!  PARTICULAR SOLUTION LAYER SOURCE TERMS
!  --------------------------------------

!  INITIAL MODULUS = 4.PI IF SOLAR SOURCES ARE INCLUDED

      FAC = ONE
      IF ( DO_SOLAR_SOURCES ) FAC = PI4

!  SHORTHAND

      NT = N_THERMAL_COEFFS

!  START USER ANGLE LOOP

      DO UM = 1, N_USER_STREAMS

!  DIRECT TERMS TO START

       DO N = 1, N_ALLLAYERS_DN
         LAYER_TSUP_DN(UM,N) = FAC * T_DIRECT_DN(UM,N)
       ENDDO
       IF ( DO_PARTLAYERS ) THEN
        DO UT = 1, N_PARTLAYERS
         LAYER_TSUP_UTDN(UM,UT) = FAC*T_UT_DIRECT_DN(UM,UT)
        ENDDO
       ENDIF

!  ONLY DO THE NEXT SECTION FOR SCATTERING SOLUTIONS
!   -- Version 2.8, 7/7/16 remove goto 
!       IF ( DO_THERMAL_TRANSONLY ) GO TO 678

       IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  LOCAL COSINE

        COSMUM = USER_STREAMS(UM)

!  WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          T_MULT_DN(N,2) = U_TNEG2(UM,N)
          T_MULT_DN(N,1) = U_TNEG1(UM,N) - COSMUM * T_MULT_DN(N,2)
          T_MULT_DN(N,0) = - T_MULT_DN(N,1)
          SPAR = T_MULT_DN(N,0) * T_DELT_USERM(N,UM)
          SPAR = SPAR + T_MULT_DN(N,1)
          SPAR = SPAR + T_MULT_DN(N,2) * DELTAU_POWER(N,2)
          LAYER_TSUP_DN(UM,N) = LAYER_TSUP_DN(UM,N) + SPAR * FAC
!     if ( n.eq.nlayers) write(*,*)'Perturb', um, n, layer_tsup_dn(um,n)
         ENDIF
        ENDDO

!  PARTIAL LAYER SOURCE TERMS, ADD THERMAL DIFFUSE PART
 
       IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
           N  = PARTLAYERS_LAYERIDX(UT)
           SUM = T_MULT_DN(N,0) * T_UTDN_USERM(UT,UM)
           DO S = 1, NT
             SUM = SUM + T_MULT_DN(N,S) * XTAU_POWER(UT,S)
           END DO
           LAYER_TSUP_UTDN(UM,UT) = &
                 LAYER_TSUP_UTDN(UM,UT) + SUM*FAC
!if (ut.eq.1) write(*,*)'Perturb',ut,n,um,layer_tsup_utdn(um,ut)
         END DO
        END IF

!  End scattering clause

       ENDIF

!  CONTINUATION POINT FOR AVOIDING SCATTERING CALCULATIONS
! 678   CONTINUE

!  END USER-STREAM LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_STERMS_DN

!  End Module

      END MODULE vlidort_thermalsup_m

