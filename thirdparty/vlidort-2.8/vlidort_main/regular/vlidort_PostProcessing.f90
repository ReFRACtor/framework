
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

! ###########################################################
! #                                                         #
! # Subroutines in this Module                              #
! #                                                         #
! #   -- Source-term integration for post-processed field   #
! #                                                         #
! #          WHOLELAYER_STERM_UP                            #
! #          WHOLELAYER_STERM_DN                            #
! #          PARTLAYER_STERM_UP                             #
! #          PARTLAYER_STERM_DN                             #
! #                                                         #
! #   -- Post-processed discrete ordinate field             #
! #       Required for Mean-I and Flux output               #
! #                                                         #
! #          QUADINTENS_LEVEL_UP                            #
! #          QUADINTENS_LEVEL_DN                            #
! #          QUADINTENS_OFFGRID_UP                          #
! #          QUADINTENS_OFFGRID_DN                          #
! #                                                         #
! ###########################################################

!  1/31/21. Version 2.8.3. Green's Function changes
!    -- Need additional Taylor routines to be used, for Green's function post-processing
!    -- Whole/part Sourceterm routines completely rewritten       (Green's function included)
!    -- QuadIntens partial-layer subroutines completely rewritten (Green's function included)

!  1/31/21. Version 2.8.3. Other post-processing changes
!    -- Use post-processing mask to cut down on observational vs lattice/doublet geometry
!    -- This is a new module. 8 subroutines formerly in VLIDORT_INTENSITY.f90 now grouped here

      MODULE vlidort_PostProcessing_m

!  1/31/21. Version 2.8.3.  ==> Need Taylor series routines for Green's function post-processing

      USE vlidort_Taylor_m, only : TAYLOR_SERIES_1, TAYLOR_SERIES_2

!  everything public

      PUBLIC

      CONTAINS

!

      SUBROUTINE WHOLELAYER_STERM_UP ( &
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,         & ! Input flags
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SOURCETERM_FLAG,             & ! Input flags
        N, M, IBEAM, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,        & ! Input numbers
        DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,                 & ! Input Beam stuff
        ITRANS_USERM, T_DELT_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,  & ! Input Green's solution
        K_REAL, K_COMPLEX, LCON, MCON, UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2, & ! Input Homog solutions
        LAYER_TSUP_UP, UPAR_UP_1, UPAR_UP_2, SIGMA_P, EMULT_UP,                & ! Input Partic Integrals
        LAYERSOURCE, PMULT_UU, PMULT_UD )                                        ! output 

!  1/31/21. Version 2.8.3. Completely rewritten to include the Green's function treatment.
!    -- Add flag DO_CLASSICAL_SOLUTION, Outputs PMULT_UU, PMULT_UD (Green's function multipliers)
!    -- Additional inputs (lines 4-5) and SIGMA_P, EMULT_UP.
!    -- Add TAYLOR_SMALL to list of usable parameters
!    -- Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

      USE VLIDORT_PARS_m, only : MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXBEAMS, &
                                 MAXEVALUES, MAXSTRMSTKS, TAYLOR_SMALL, ZERO, ONE, PI4

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3. Greens function, add flag DO_CLASSICAL_SOLUTION

      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           SOURCETERM_FLAG

!  Numbers

      INTEGER, INTENT (IN) ::           N, M, IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES

!  1/31/21. Version 2.8.3. Introduce Post-processing masks. Remove OBSERVATION_GEOMETRY Flag

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Add TAYLOR Order 

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER

!  1/31/21. Version 2.8.3. Input optical depths

      DOUBLE PRECISION, INTENT (IN)   :: DELTAU_VERT  ( MAXLAYERS )

!  1/31/21. Version 2.8.3. Beam parameterization inputs needed for Green's
!    -- Initial/Layer transmittance factors for solar beams, and divided by user-cosines

      INTEGER         , INTENT (IN)   :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN)   :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined streams

      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Eigenvalue bookkeeping

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )

!  boundary value constants

      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M ( MAXEVALUES, MAXLAYERS )

!  User homogenous solutions and multipliers

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  User particular-integral solutions and multipliers/coefficients
!    ==> Thermal Layer source terms (direct + diffuse)
!    ==> Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  1/31/21. Version 2.8.3. Add SIGMA_P to multiplier input

      DOUBLE PRECISION, INTENT (IN) ::  EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_P  ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  output
!  ------

!  Layer source terms

      DOUBLE PRECISION, INTENT (INOUT) :: LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

!  1/31/21. Version 2.8.3.
!   -- Source function integrated Green function multipliers (whole layer)

      DOUBLE PRECISION, INTENT (INOUT) ::  PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (INOUT) ::  PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

      LOGICAL ::          L1
      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, M1, IB, LUM
      DOUBLE PRECISION :: H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION :: LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

!  1/31/21. Version 2.8.3. Additional variables for Green's function solution

      DOUBLE PRECISION :: SPAR(MAXSTOKES)
      DOUBLE PRECISION :: WDEL, ITRANS, ITRANSWDEL, SD, SU
      DOUBLE PRECISION :: EPS, YFAC, FAC1, FAC2, MULT

!  Not necessary now.
!      DOUBLE PRECISION, DIMENSION(MAXSTOKES)  ::  SGW =  (/ 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)
!      DOUBLE PRECISION, DIMENSION(MAXSTOKES)  ::  SGW =  (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)   THIS IS THE CORRECT ONE

!  Fourier number (debug only)

      M1 = M
      L1 = SOURCETERM_FLAG

!  Local user indices

      IB  = IBEAM

!   Very important to zero output term (bug solved 04/20/05)
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO O1 = 1, NSTOKES
            LAYERSOURCE(UM,O1) = ZERO
         ENDDO
      ENDDO

!      IF ( .not. L1 ) RETURN   ! @@@ Rob, wrong place.

!  Avoid this section if thermal transmittance only
!  remove GOTO 6789 label, Version 2.8. 7/5/16
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  @@@ Moved here

        IF ( .not. L1 ) RETURN

!  Offset

        KO1 = K_REAL(N) + 1

!  Whole layer source function Homogeneous only
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

        DO LUM = 1, N_PPSTREAMS
           UM = PPSTREAM_MASK(LUM,IB)
           DO O1 = 1, NSTOKES
              H_R = ZERO
              DO K = 1, K_REAL(N)
                 LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
                 MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
                 H_R = H_R + LUXR * HMULT_2(K,UM,N) + MUXR * HMULT_1(K,UM,N)
              ENDDO
              H_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                 K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                 LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - &
                          LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                 LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + &
                          LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                 MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - &
                          MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                 MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + &
                          MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                 H_CR = H_CR + LUX_CR * HMULT_2(K1,UM,N) - &
                               LUX_CI * HMULT_2(K2,UM,N) + &
                               MUX_CR * HMULT_1(K1,UM,N) - &
                               MUX_CI * HMULT_1(K2,UM,N)
              ENDDO
              LAYERSOURCE(UM,O1) = H_R + H_CR
           ENDDO
        ENDDO

!  End scattering

      ENDIF

!  Continuation point, removed Version 2.8 7/5/16
! 6789 continue

!  Add thermal emission term (direct and diffuse)
!     Modulus 4.pi if solar sources are included (taken care of earlier)
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         O1 = 1
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + LAYER_TSUP_UP(UM,N)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  CLASSICAL SOLUTION. Add particular integral contribution (solar term)
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( DO_CLASSICAL_SOLUTION ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
               SFOR2 =  EMULT_UP(LUM,N,IB) * UPAR_UP_2(UM,O1,N)
               LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
            ENDDO
         ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. GREENS FUNCTION SOLUTION
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Layer quantities

         WDEL       = T_DELT_MUBAR(N,IB)
         ITRANS     = INITIAL_TRANS(N,IB)
         ITRANSWDEL = - ITRANS * WDEL

!  Get the multipliers and store them.  Add contributions to the Greens function solution

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SPAR(1:NSTOKES) = ZERO
            DO K = 1, K_REAL(N)
               if ( ABS(GAMMA_M(K,N)) .lt. TAYLOR_SMALL ) THEN
                  EPS   = GAMMA_M(K,N)      ; FAC1 = ONE
                  YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL * T_DELT_USERM(N,UM)
                  CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAU_VERT(N), FAC1, FAC2, ONE, MULT )
                  SD = ITRANS_USERM (N,LUM,IB) * MULT 
               else
                  SD = ( ITRANS * HMULT_2(K,UM,N) - EMULT_UP(LUM,N,IB) ) / GAMMA_M(K,N)
               endif
               SU = ( ITRANSWDEL * HMULT_1(K,UM,N) + EMULT_UP(LUM,N,IB) ) / GAMMA_P(K,N)

!  Rob Change 1/51/21. Replace this code with better definitions of PMULT terms
!               PMULT_UD(K,UM,N) = SD * ATERM_SAVE(K,N)
!               PMULT_UU(K,UM,N) = SU * BTERM_SAVE(K,N)
!               do O1 = 1, NSTOKES
!                 T1 = UHOM_UPDN(UM,O1,K,N) * PMULT_UD(K,UM,N)
!                 T2 = SGW(O1) * UHOM_UPUP(UM,O1,K,N) * PMULT_UU(K,UM,N)
!                 SPAR(O1) = SPAR(O1) + T1 + T2
!               ENDDO

               PMULT_UD(K,UM,N) = SD
               PMULT_UU(K,UM,N) = SU
               do O1 = 1, NSTOKES
                 SPAR(O1) = SPAR(O1) + UHOM_UPDN(UM,O1,K,N) * PMULT_UD(K,UM,N) * ATERM_SAVE(K,N) &
                                     + UHOM_UPUP(UM,O1,K,N) * PMULT_UU(K,UM,N) * BTERM_SAVE(K,N)
               ENDDO

            ENDDO
            LAYERSOURCE(UM,1:NSTOKES) = LAYERSOURCE(UM,1:NSTOKES) + SPAR(1:NSTOKES)
         ENDDO

!  End Green's function calculation

      ENDIF

!  If operating in Ms-mode only, finished

      IF ( DO_MSMODE_VLIDORT ) RETURN

!  Full radiance mode, add single scatter part
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO O1 = 1, NSTOKES
            SFOR1 = UPAR_UP_1(UM,O1,N) * EMULT_UP(LUM,N,IB)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE WHOLELAYER_STERM_UP

!

      SUBROUTINE WHOLELAYER_STERM_DN ( &
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,         & ! Input flags
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SOURCETERM_FLAG,             & ! Input flags
        N, M, IBEAM, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,        & ! Input numbers
        DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,                 & ! Input Beam stuff
        ITRANS_USERM, T_DELT_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,  & ! Input Green's solution
        K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, HMULT_1, HMULT_2, & ! Input Homog solutions
        LAYER_TSUP_DN, UPAR_DN_1, UPAR_DN_2, SIGMA_M, EMULT_DN,                & ! Input Partic Integrals
        LAYERSOURCE, PMULT_DD, PMULT_DU )                                        ! output 

!  1/31/21. Version 2.8.3. Completely rewritten to include the Green's function treatment.
!    -- Add flag DO_CLASSICAL_SOLUTION, Outputs PMULT_DD, PMULT_DU (Green's function multipliers)
!    -- Additional inputs (lines 4-5) and SIGMA_M, EMULT_DN.
!    -- Add TAYLOR_SMALL to list of usable parameters
!    -- Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

      USE VLIDORT_PARS_m, only : MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXBEAMS, &
                                 MAXEVALUES, MAXSTRMSTKS, TAYLOR_SMALL, ZERO, ONE, PI4, TAYLOR_SMALL, SDU

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3. Greens function, add flag DO_CLASSICAL_SOLUTION

      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           SOURCETERM_FLAG

!  Numbers

      INTEGER, INTENT (IN) ::           N, M, IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES

!  1/31/21. Version 2.8.3. Introduce Post-processing masks. Remove OBSERVATION_GEOMETRY Flag

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Add TAYLOR Order 

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER

!  1/31/21. Version 2.8.3. Input optical depths

      DOUBLE PRECISION, INTENT (IN)   :: DELTAU_VERT  ( MAXLAYERS )

!  1/31/21. Version 2.8.3. Beam parameterization inputs needed for Green's
!    -- Initial/Layer transmittance factors for solar beams, and divided by user-cosines

      INTEGER         , INTENT (IN)   :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN)   :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined streams

      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Eigenvalue bookkeeping

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )

!  boundary value constants

      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M ( MAXEVALUES, MAXLAYERS )

!  User homogenous solutions and multipliers

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  User particular-integral solutions and multipliers/coefficients
!    ==> Thermal Layer source terms (direct + diffuse)
!    ==> Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  1/31/21. Version 2.8.3. Add SIGMA_M to multiplier input

      DOUBLE PRECISION, INTENT (IN) ::  EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_M  ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  output
!  ------

!  Layer source terms

      DOUBLE PRECISION, INTENT (INOUT) :: LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

!  1/31/21. Version 2.8.3.
!   -- Source function integrated Green function multipliers (whole layer)

      DOUBLE PRECISION, INTENT (INOUT) ::  PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (INOUT) ::  PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)

!  local variables

      LOGICAL ::          L1
      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, IB, LUM
      DOUBLE PRECISION :: H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION :: LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

!  1/31/21. Version 2.8.3. Additional variables for Green's function solution

      DOUBLE PRECISION :: SPAR(MAXSTOKES)
      DOUBLE PRECISION :: WDEL, ITRANS, ITRANSWDEL, SD, SU
      DOUBLE PRECISION :: EPS, YFAC, FAC1, FAC2, MULT

!  Local user indices

      IB = IBEAM
      L1 = SOURCETERM_FLAG

!  Zeroing is very important here
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO O1 = 1, NSTOKES
            LAYERSOURCE(UM,O1) = ZERO
         ENDDO
      ENDDO

!      IF ( .not. L1 ) RETURN   ! @@@ Rob, wrong place.

!  Avoid this section if thermal transmittance only
!  remove GOTO 6789 label, Version 2.8. 7/5/16
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  If no source term return
!  @@@ Moved here

        IF ( .not. L1 ) RETURN

!  Offset

        KO1 = K_REAL(N) + 1

!  Whole layer source function, homogeneous solution contribution
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

        DO LUM = 1, N_PPSTREAMS
           UM = PPSTREAM_MASK(LUM,IB)
           DO O1 = 1, NSTOKES
              H_R = ZERO
              DO K = 1, K_REAL(N)
                 LUXR = LCON(K,N) * UHOM_DNDN(UM,O1,K,N)
                 MUXR = MCON(K,N) * UHOM_DNUP(UM,O1,K,N)
                 H_R = H_R + LUXR * HMULT_1(K,UM,N) + MUXR * HMULT_2(K,UM,N)
              ENDDO
              H_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                 K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                 LUX_CR = LCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - &
                          LCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
                 LUX_CI = LCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + &
                          LCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
                 MUX_CR = MCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - &
                          MCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
                 MUX_CI = MCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + &
                          MCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
                 H_CR = H_CR + LUX_CR * HMULT_1(K1,UM,N) - &
                              LUX_CI * HMULT_1(K2,UM,N) + &
                              MUX_CR * HMULT_2(K1,UM,N) - &
                              MUX_CI * HMULT_2(K2,UM,N)
              ENDDO
              LAYERSOURCE(UM,O1) = H_R + H_CR
           ENDDO
        ENDDO

!  End scattering

      ENDIF

!  Continuation point
! 6789 continue

!  Add thermal emission term (direct and diffuse)
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         O1 = 1
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + LAYER_TSUP_DN(UM,N)*TM
         ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  CLASSICAL SOLUTION. Add particular integral contribution (solar term)
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( DO_CLASSICAL_SOLUTION ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
               SFOR2 =  EMULT_DN(LUM,N,IB) * UPAR_DN_2(UM,O1,N)
               LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
            ENDDO
         ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. GREENS FUNCTION SOLUTION
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Layer quantities

         WDEL       = T_DELT_MUBAR(N,IB)
         ITRANS     = INITIAL_TRANS(N,IB)
         ITRANSWDEL = - ITRANS * WDEL

!  Get the multipliers and store them. Add contributions to the Greens function solution

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SPAR(1:NSTOKES) = ZERO
            DO K = 1, K_REAL(N)
               if ( ABS(GAMMA_M(K,N)) .lt. TAYLOR_SMALL ) THEN
                  EPS   = GAMMA_M(K,N)      ; FAC1 = T_DELT_USERM(N,UM)
                  YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = WDEL
                  CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAU_VERT(N), FAC1, FAC2, ONE, MULT )
                  SD = ITRANS_USERM (N,LUM,IB) * MULT 
               else
                  SD = ( ITRANS * HMULT_1(K,UM,N) - EMULT_DN(LUM,N,IB) ) / GAMMA_M(K,N)
               endif
               SU = ( ITRANSWDEL * HMULT_2(K,UM,N) + EMULT_DN(LUM,N,IB) ) / GAMMA_P(K,N)

!  Rob Change 1/51/21. Replace this code with better definitions of PMULT terms
!               PMULT_DD(K,UM,N) = SD * ATERM_SAVE(K,N)
!               PMULT_DU(K,UM,N) = SU * BTERM_SAVE(K,N)
!               SPAR(1:NSTOKES) = SPAR(1:NSTOKES) + UHOM_DNDN(UM,1:NSTOKES,K,N) * PMULT_DD(K,UM,N) &
!                                                 + UHOM_DNUP(UM,1:NSTOKES,K,N) * PMULT_DU(K,UM,N)
               PMULT_DD(K,UM,N) = SD 
               PMULT_DU(K,UM,N) = SU
               SPAR(1:NSTOKES) = SPAR(1:NSTOKES) + UHOM_DNDN(UM,1:NSTOKES,K,N) * PMULT_DD(K,UM,N) * ATERM_SAVE(K,N) &
                                                 + UHOM_DNUP(UM,1:NSTOKES,K,N) * PMULT_DU(K,UM,N) * BTERM_SAVE(K,N)
            ENDDO
            LAYERSOURCE(UM,1:NSTOKES) = LAYERSOURCE(UM,1:NSTOKES) + SPAR(1:NSTOKES)
         ENDDO

!  End Green's function calculation

      ENDIF

!  If operating in MS-mode only, finished

      IF ( DO_MSMODE_VLIDORT ) RETURN

!  Full radiance mode, add single scatter part
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO O1 = 1, NSTOKES
            SFOR1 = UPAR_DN_1(UM,O1,N) * EMULT_DN(LUM,N,IB)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE WHOLELAYER_STERM_DN

!

      SUBROUTINE PARTLAYER_STERM_UP ( &
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,         & ! Input flags
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SOURCETERM_FLAG,             & ! Input flags
        N, UT, IBEAM, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,       & ! Input numbers
        DELTAU_VERT, PARTAU_VERT, T_UTUP_USERM, ITRANS_USERM,                  & ! Input delt
        BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                & ! Input Beam stuff
        ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, K_REAL, K_COMPLEX,           & ! input Green s function
        LCON, MCON, UHOM_UPDN, UHOM_UPUP, UT_HMULT_UU, UT_HMULT_UD,            & ! Input Homog solutions
        LAYER_TSUP_UTUP, UPAR_UP_1, UPAR_UP_2, UT_EMULT_UP, SIGMA_P,           & ! Input Partic Integrals
        LAYERSOURCE, UT_PMULT_UU, UT_PMULT_UD )                                  ! output 

!  1/31/21. Version 2.8.3. Completely rewritten to include the Green's function treatment.
!   -- Add flag DO_CLASSICAL_SOLUTION, Outputs UT_PMULT_UU, UT_PMULT_UD (Green's function multipliers)
!   -- Additional inputs (lines 4-5) and SIGMA_P.
!   -- Add TAYLOR_SMALL to list of usable parameters
!    -- Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

      USE VLIDORT_PARS_m, only : MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXBEAMS, &
                                 MAX_PARTLAYERS, MAXEVALUES, MAXSTRMSTKS, TAYLOR_SMALL, ZERO, ONE, PI4

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3. Greens function, add flag DO_CLASSICAL_SOLUTION

      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           SOURCETERM_FLAG

!  Numbers

      INTEGER, INTENT (IN) ::           N, UT, IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES

!  1/31/21. Version 2.8.3. Introduce Post-processing masks. Remove OBSERVATION_GEOMETRY Flag

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Add TAYLOR Order 

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER

!  1/31/21. Version 2.8.3.Input optical depths (Need partial now, 3/16/20GF )

      DOUBLE PRECISION, INTENT (IN)   :: DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN)   :: PARTAU_VERT  ( MAX_PARTLAYERS )

!   1/31/21. Version 2.8.3.Transmittance factors for user-defined streams, aaverage-secant streams

      DOUBLE PRECISION, INTENT (IN)   :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN)   :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Initial/Layer transmittance factors for solar beams, and divided by user-cosines

      INTEGER         , INTENT (IN)   :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN)   :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M ( MAXEVALUES, MAXLAYERS )

!  Eigenvalue bookkeeping

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )

!  boundary value constants

      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solution

      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  User homog and particular solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers. ==> 1/31/21. Version 2.8.3. Add Sigma_P

      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_P     ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Output
!  ------

!  source term

      DOUBLE PRECISION, INTENT (INOUT) :: LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

!  1/31/21. Version 2.8.3. Green's function multipliers, post processed

      DOUBLE PRECISION, INTENT (INOUT) ::  UT_PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (INOUT) ::  UT_PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

      LOGICAL ::          L1
      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, IB, LUM
      DOUBLE PRECISION :: H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION :: LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

!  1/31/21. Version 2.8.3. Additional variables for Green's function solution

      DOUBLE PRECISION :: SPAR(MAXSTOKES)
      DOUBLE PRECISION :: WDEL, ITRANS, ITRANSWDEL, SD, SU
      DOUBLE PRECISION :: EPS, YFAC, FAC1, FAC2, MULT1, MULT2

!  Local user indices

      IB = IBEAM
      L1 = SOURCETERM_FLAG

!  Zeroing is very important here
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      SFOR2 = ZERO
      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO O1 = 1, NSTOKES
            LAYERSOURCE(UM,O1) = ZERO
         ENDDO
      ENDDO

!      IF ( .not. L1 ) RETURN   ! @@@ Rob, wrong place.

!  Avoid this section if thermal transmittance only
!  remove GOTO 6789 label, Version 2.8. 7/5/16
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  If no source term return
!  @@@ Moved here

        IF ( .not. L1 ) RETURN

!  Offset

        KO1 = K_REAL(N) + 1

!  Whole layer source function, homogeneous contribution
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

        DO LUM = 1, N_PPSTREAMS
           UM = PPSTREAM_MASK(LUM,IB)
           DO O1 = 1, NSTOKES
              H_R = ZERO
              DO K = 1, K_REAL(N)
                 LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
                 MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
                 H_R = H_R + LUXR * UT_HMULT_UD(K,UM,UT) + MUXR * UT_HMULT_UU(K,UM,UT)
              ENDDO
              H_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                 K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                 LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - &
                          LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                 LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + &
                          LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                 MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - &
                          MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                 MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + &
                          MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                 H_CR = H_CR + LUX_CR * UT_HMULT_UD(K1,UM,UT) - &
                               LUX_CI * UT_HMULT_UD(K2,UM,UT) + &
                               MUX_CR * UT_HMULT_UU(K1,UM,UT) - &
                               MUX_CI * UT_HMULT_UU(K2,UM,UT)
              ENDDO
              LAYERSOURCE(UM,O1) = H_R + H_CR
           ENDDO
        ENDDO

!  End scattering

      ENDIF

!  Continuation point
! 6789 continue

!  Add thermal term (direct and diffuse)
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         O1 = 1
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + LAYER_TSUP_UTUP(UM,UT)*TM
         ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Particular integral contribution. Classical
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( DO_CLASSICAL_SOLUTION ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
              SFOR2 =  UT_EMULT_UP(LUM,UT,IB) * UPAR_UP_2(UM,O1,N)
              LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2
            ENDDO
         ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. Add the Green's function calculation
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Layer quantities

         WDEL       = T_DELT_MUBAR(N,IB)
         ITRANS     = INITIAL_TRANS(N,IB)
         ITRANSWDEL = - ITRANS * WDEL

!  Get the multipliers and store them. Add contributions to the Greens function solution

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SPAR(1:NSTOKES) = ZERO
            DO K = 1, K_REAL(N)
               if ( ABS(GAMMA_M(K,N)) .lt. TAYLOR_SMALL ) THEN
                  EPS   = GAMMA_M(K,N)     ; FAC1 = - T_UTDN_MUBAR(UT,IB) 
                  YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL * T_UTUP_USERM(UT,UM) 
                  CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTAU_VERT(UT), ZERO, FAC1, ONE, MULT1 )
                  CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAU_VERT(N),  ZERO, FAC2, ONE, MULT2 )
                  SD = ITRANS_USERM (N,LUM,IB) * ( MULT1 + MULT2 ) 
               else
                  SD = ( ITRANS * UT_HMULT_UD(K,UM,UT) - UT_EMULT_UP(LUM,UT,IB) ) / GAMMA_M(K,N)
               endif
               SU = ( ITRANSWDEL * UT_HMULT_UU(K,UM,UT) + UT_EMULT_UP(LUM,UT,IB) ) / GAMMA_P(K,N)

!  Rob Change 1/51/21. Replace this code with better definitions of PMULT terms
!               UT_PMULT_UD(K,UM,UT) = SD * ATERM_SAVE(K,N)
!               UT_PMULT_UU(K,UM,UT) = SU * BTERM_SAVE(K,N)
!               SPAR(1:NSTOKES) = SPAR(1:NSTOKES) + UHOM_UPDN(UM,1:NSTOKES,K,N) * UT_PMULT_UD(K,UM,UT) &
!                                                 + UHOM_UPUP(UM,1:NSTOKES,K,N) * UT_PMULT_UU(K,UM,UT)

               UT_PMULT_UD(K,UM,UT) = SD
               UT_PMULT_UU(K,UM,UT) = SU
               SPAR(1:NSTOKES) = SPAR(1:NSTOKES) + UHOM_UPDN(UM,1:NSTOKES,K,N) * UT_PMULT_UD(K,UM,UT) * ATERM_SAVE(K,N) &
                                                 + UHOM_UPUP(UM,1:NSTOKES,K,N) * UT_PMULT_UU(K,UM,UT) * BTERM_SAVE(K,N)

!  End eigenvalue loop

            ENDDO

!  Add Green's function result to the total

            LAYERSOURCE(UM,1:NSTOKES) = LAYERSOURCE(UM,1:NSTOKES) + SPAR(1:NSTOKES)

!  End user stream loop

         ENDDO

!  End Green's function calculation

      ENDIF

!  If NOT operating in MS-mode only, add single scatter part
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
               SFOR1 = UT_EMULT_UP(LUM,UT,IBEAM) * UPAR_UP_1(UM,O1,N)
               LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE PARTLAYER_STERM_UP

!

      SUBROUTINE PARTLAYER_STERM_DN ( &
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,         & ! Input flags
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SOURCETERM_FLAG,             & ! Input flags
        N, UT, IBEAM, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,       & ! Input numbers
        PARTAU_VERT, T_UTDN_USERM, ITRANS_USERM,                               & ! Input delt/trans
        BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                & ! Input Beam stuff
        ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, K_REAL, K_COMPLEX,           & ! input Greens function
        LCON, MCON, UHOM_DNDN, UHOM_DNUP, UT_HMULT_DD, UT_HMULT_DU,            & ! Input Homog solutions
        LAYER_TSUP_UTDN, UPAR_DN_1, UPAR_DN_2, UT_EMULT_DN, SIGMA_M,           & ! Input Partic Integrals
        LAYERSOURCE, UT_PMULT_DU, UT_PMULT_DD )                                  ! output 

!  1/31/21. Version 2.8.3. Completely rewritten to include the Green's function treatment.
!   -- Add flag DO_CLASSICAL_SOLUTION, Outputs UT_PMULT_DU, UT_PMULT_DD (Green's function multipliers)
!   -- Additional inputs (lines 4-5) and SIGMA_M.
!   -- Add TAYLOR_SMALL to list of usable parameters
!    -- Use post-processing mask system (N_PPSTREAMS, PPSTREAM_MASK).

      USE VLIDORT_PARS_m, only : MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXBEAMS, &
                                 MAX_PARTLAYERS, MAXEVALUES, MAXSTRMSTKS, TAYLOR_SMALL, ZERO, ONE, PI4, SDU

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3. Greens function, add flag DO_CLASSICAL_SOLUTION

      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           SOURCETERM_FLAG

!  Numbers

      INTEGER, INTENT (IN) ::           N, UT, IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES

!  1/31/21. Version 2.8.3. Introduce Post-processing masks. Remove OBSERVATION_GEOMETRY Flag

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Add TAYLOR Order 

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER

!  1/31/21. Version 2.8.3. Input optical depths

      DOUBLE PRECISION, INTENT (IN)   :: PARTAU_VERT  ( MAX_PARTLAYERS )

!   1/31/21. Version 2.8.3. Transmittance factors for user-defined streams, aaverage-secant streams

      DOUBLE PRECISION, INTENT (IN)   :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN)   :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Initial/Layer transmittance factors for solar beams, and divided by user-cosines

      INTEGER         , INTENT (IN)   :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN)   :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M ( MAXEVALUES, MAXLAYERS )

!  Eigenvalue bookkeeping

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )

!  boundary value constants

      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solution

      DOUBLE PRECISION, INTENT (IN) ::  LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  User homog and particular solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers. ==> 1/31/21. Version 2.8.3. Add Sigma_M

      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_M     ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Output
!  ------

!  source term

      DOUBLE PRECISION, INTENT (INOUT) :: LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

!  1/31/21. Version 2.8.3. Green's function multipliers, post processed

      DOUBLE PRECISION, INTENT (INOUT) ::  UT_PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (INOUT) ::  UT_PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

      LOGICAL ::          L1
      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, IB, LUM
      DOUBLE PRECISION :: H_R, H_CR, SFOR1, SFOR2, TM
      DOUBLE PRECISION :: LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI

!  1/31/21. Version 2.8.3. Additional variables for Green's function solution

      DOUBLE PRECISION :: SPAR(MAXSTOKES)
      DOUBLE PRECISION :: WDEL, ITRANS, ITRANSWDEL, SD, SU
      DOUBLE PRECISION :: EPS, YFAC, FAC1, FAC2, MULT

!  Local user indices

      IB  = IBEAM
      L1    = SOURCETERM_FLAG

!  Zeroing is very important here
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      SFOR2 = ZERO
      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO O1 = 1, NSTOKES
            LAYERSOURCE(UM,O1) = ZERO
         ENDDO
      ENDDO

!      IF ( .not. L1 ) RETURN   ! @@@ Rob, wrong place.

!  Avoid this section if thermal transmittance only
!  remove GOTO 6789 label, Version 2.8. 7/5/16
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  If no source term return
!  @@@ Moved here

        IF ( .not. L1 ) RETURN

!  Offset

        KO1 = K_REAL(N) + 1

!  Whole layer source function, homogeneous contribution
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

        DO LUM = 1, N_PPSTREAMS
           UM = PPSTREAM_MASK(LUM,IB)
           DO O1 = 1, NSTOKES
              H_R = ZERO
              DO K = 1, K_REAL(N)
                 LUXR = LCON(K,N) * UHOM_DNDN(UM,O1,K,N)
                 MUXR = MCON(K,N) * UHOM_DNUP(UM,O1,K,N)
                 H_R = H_R + LUXR * UT_HMULT_DD(K,UM,UT) + MUXR * UT_HMULT_DU(K,UM,UT)
              ENDDO
              H_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                 K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                 LUX_CR = LCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - &
                          LCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
                 LUX_CI = LCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + &
                          LCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
                 MUX_CR = MCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - &
                          MCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
                 MUX_CI = MCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + &
                          MCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
                 H_CR = H_CR + LUX_CR * UT_HMULT_DD(K1,UM,UT) - &
                               LUX_CI * UT_HMULT_DD(K2,UM,UT) + &
                               MUX_CR * UT_HMULT_DU(K1,UM,UT) - &
                               MUX_CI * UT_HMULT_DU(K2,UM,UT)
              ENDDO
              LAYERSOURCE(UM,O1) = H_R + H_CR
           ENDDO
        ENDDO

!  End scattering

      ENDIF

!  Continuation point
! 6789 continue

!  Add thermal term (direct and diffuse)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         O1 = 1
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + LAYER_TSUP_UTDN(UM,UT)*TM
         ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Particular integral contribution. Classical
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( DO_CLASSICAL_SOLUTION ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
              SFOR2 =  UT_EMULT_DN(LUM,UT,IB) * UPAR_DN_2(UM,O1,N)
              LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR2 
            ENDDO
         ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. Add the Green's function calculation
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Layer quantities

         WDEL       = T_DELT_MUBAR(N,IB)
         ITRANS     = INITIAL_TRANS(N,IB)
         ITRANSWDEL = - ITRANS * WDEL

!  Get the multipliers and store them. Add contributions to the Greens function solution
!    --bug found 3/1/21. Argument list for TAYLOR_SERIES_2 was wrong.

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SPAR(1:NSTOKES) = ZERO
            DO K = 1, K_REAL(N)
               if ( ABS(GAMMA_M(K,N)) .lt. TAYLOR_SMALL ) THEN
                  EPS   = GAMMA_M(K,N)      ; FAC1 = T_UTDN_USERM(UT,UM)  
                  YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = T_UTDN_MUBAR(UT,IB)
                  !CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTAU_VERT(UT), ZERO, FAC1, ONE, MULT )
                  CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTAU_VERT(UT), FAC1, FAC2, ONE, MULT ) 
                  SD = ITRANS_USERM (N,LUM,IB) * MULT
               else
                  SD = ( ITRANS * UT_HMULT_DD(K,UM,UT) - UT_EMULT_DN(LUM,UT,IB) ) / GAMMA_M(K,N)
               endif
               SU = ( ITRANSWDEL * UT_HMULT_DU(K,UM,UT) + UT_EMULT_DN(LUM,UT,IB) ) / GAMMA_P(K,N)

!  Rob Change 1/51/21. Replace this code with better definitions of PMULT terms
!               UT_PMULT_DD(K,UM,UT) = SD * ATERM_SAVE(K,N)
!               UT_PMULT_DU(K,UM,UT) = SU * BTERM_SAVE(K,N)
!               SPAR(1:NSTOKES) = SPAR(1:NSTOKES) + UHOM_DNDN(UM,1:NSTOKES,K,N) * UT_PMULT_DD(K,UM,UT) &
!                                                 + UHOM_DNUP(UM,1:NSTOKES,K,N) * UT_PMULT_DU(K,UM,UT)
               UT_PMULT_DD(K,UM,UT) = SD
               UT_PMULT_DU(K,UM,UT) = SU
               SPAR(1:NSTOKES) = SPAR(1:NSTOKES) + UHOM_DNDN(UM,1:NSTOKES,K,N) * UT_PMULT_DD(K,UM,UT) * ATERM_SAVE(K,N) &
                                                 + UHOM_DNUP(UM,1:NSTOKES,K,N) * UT_PMULT_DU(K,UM,UT) * BTERM_SAVE(K,N)

!  End eigenvalue loop

            ENDDO

!  Add Green's function result to the total

            LAYERSOURCE(UM,1:NSTOKES) = LAYERSOURCE(UM,1:NSTOKES) + SPAR(1:NSTOKES)

!  End user stream loop

         ENDDO

!  End Green's solution

      ENDIF

!  If NOT operating in MS-mode only, add single scatter part
!    -- 1/31/21. Version 2.8.3. Use post-processing masks.

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
               SFOR1 = UT_EMULT_DN(LUM,UT,IBEAM) * UPAR_DN_1(UM,O1,N)
               LAYERSOURCE(UM,O1) = LAYERSOURCE(UM,O1) + SFOR1
            ENDDO
          ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE PARTLAYER_STERM_DN

!

      SUBROUTINE QUADINTENS_LEVEL_UP ( &
          DO_THERMAL_TRANSONLY, DO_INCLUDE_BOAFLUX,                 & ! Input Flags
          NLEVEL, UTA, NSTOKES, NSTREAMS, NLAYERS,                  & ! Input numbers
          FLUX_MULTIPLIER, BOAFLUX, QUAD_STREAMS, T_DELT_DISORDS,   & ! Input Quad/Flux/Trans
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,    & ! Input Homog. solutions
          WLOWER, WUPPER, LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE, & ! Input thermal and PI solutions
          QSTOKES_F )                                                 ! OUTPUT

!  1/31/21. Version 2.8.3. Added Green's function capability; no change to the present subroutine

!   Version 2.8.1, Control for BOA illumination added, 3/23/19

      USE VLIDORT_PARS_m, Only : MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES, MAXLAYERS,  &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO

      IMPLICIT NONE

!  Flags and numbers

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_BOAFLUX ! New 3/23/19
      
      INTEGER, INTENT (IN) ::          NLEVEL, UTA
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Flux, quadrature and discrete ordinate trans.

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: BOAFLUX                       ! BOA Flux (new 3/23/19)
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  RTE solutions


      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, NL, I, I1, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPAR, SHOM, HOM1, HOM2, SHOM_R, FLUX
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI

!  For those optical depths at layer boundaries
!  --------------------------------------------

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we're
!  looking at the intensity at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling intensity
!  at the bottom of the atmosphere (treated separately).

      NL = NLEVEL
      N  = NL + 1

!  Lowest level contributions
!  ==========================

!  Lowest level, thermal transmittance only

      IF ( NL .EQ. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I,O1)
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix 2/10/2011 - zero-ize other elements of stokes vector
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For the lowest level, scattering solution
!  -----------------------------------------

      IF ( NL .EQ. NLAYERS ) THEN

!  Start main loops

        KO1 = K_REAL(NL) + 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NL)
              LXR = LCON(K,NL)*SOLA_XPOS(I1,O1,K,NL)
              MXR = MCON(K,NL)*SOLB_XNEG(I1,O1,K,NL)
              HOM1 = LXR * T_DELT_EIGEN(K,NL)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NL)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,NL) * SOLA_XPOS(I1,O1,K1,NL) - &
                        LCON(K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
              LXR_CI =  LCON(K1,NL) * SOLA_XPOS(I1,O1,K2,NL) + &
                        LCON(K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
              MXR_CR =  MCON(K1,NL) * SOLB_XNEG(I1,O1,K1,NL) - &
                        MCON(K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
              HOM1CR =  LXR_CR*T_DELT_EIGEN(K1,NL) &
                      - LXR_CI*T_DELT_EIGEN(K2,NL)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WLOWER(I1,O1,NL)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish streams/stokes loops

          ENDDO
        ENDDO

!  Add BOA flux if flagged
!    3/1/21. Version 2.8.3. Bug in the code, there was no I = 1, NSTREAMS loop.

        IF ( DO_INCLUDE_BOAFLUX ) THEN
          O1 = 1 ; FLUX = FLUX_MULTIPLIER * BOAFLUX
          DO I = 1, NSTREAMS
            QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX
          ENDDO
        ENDIF

!  End lowest level clause

      ENDIF

!  For other levels in the atmosphere
!  ==================================

!  For other levels, thermal transmittance only

      IF ( NL .NE. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I,O1)
          DO K = NLAYERS, N, -1
            TPROP = T_WUPPER(I1,K)/QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For other levels, scattering solution
!  -------------------------------------

      IF ( NL .NE. NLAYERS ) THEN

!  start main loops

        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR = LCON(K,N)*SOLA_XPOS(I1,O1,K,N)
              MXR = MCON(K,N)*SOLB_XNEG(I1,O1,K,N)
              HOM1 = LXR
              HOM2 = MXR * T_DELT_EIGEN(K,N)
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - &
                        LCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
              MXR_CR =  MCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - &
                        MCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
              MXR_CI =  MCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + &
                        MCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
              HOM1CR = LXR_CR
              HOM2CR =  MXR_CR * T_DELT_EIGEN(K1,N) &
                      - MXR_CI * T_DELT_EIGEN(K2,N)
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WUPPER(I1,O1,N)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish streams/stokes loops

          ENDDO
        ENDDO

!  Version 2.8.1, Add Transmittance of BOA flux. 3/23/19
       
        IF ( DO_INCLUDE_BOAFLUX ) THEN
           O1 = 1
           DO I = 1, NSTREAMS
              FLUX = BOAFLUX
              DO K = NLAYERS, N, -1
                 FLUX = FLUX * T_DELT_DISORDS(I,K)
              enddo
              QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * FLUX
            ENDDO
         ENDIF

!  End level clause
         
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADINTENS_LEVEL_UP

!

      SUBROUTINE QUADINTENS_LEVEL_DN ( &
          DO_THERMAL_TRANSONLY, DO_INCLUDE_TOAFLUX,              & ! Input flags
          NLEVEL, UTA, NSTOKES, NSTREAMS, FLUX_MULTIPLIER,       & ! Input numbers and Flux
          TOAFLUX, QUAD_STREAMS, T_DELT_DISORDS,                 & ! Input TOAFlux, stream transmittances
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, & ! Input Homog. solutions
          WLOWER, LCON, MCON, T_WLOWER,                          & ! Input thermal and PI solutions
          QSTOKES_F )                                              ! OUTPUT

!  1/31/21. Version 2.8.3. Added Green's function capability; no change to the present subroutine

!   Version 2.8.1, Control for TOA illumination added, 3/23/19
        
      USE VLIDORT_PARS_m, Only : MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES, MAXLAYERS,  &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO

      IMPLICIT NONE

!  Flags and numbers

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_TOAFLUX ! New 3/23/19
      
      INTEGER, INTENT (IN) ::          NLEVEL, UTA
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS

!  Flux multiplier, quadrature and discrete ordinate trans. TOA Flux (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) :: TOAFLUX
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  RTE solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )


!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP, FLUX
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR

!  For those optical depths at layer boundaries
!  --------------------------------------------

      N = NLEVEL

!  Downwelling diffuse radiation at TOA ( or N = 0 ) is zero
!  Version 2.8a, Add illumination at TOA with ISOTROPIC flux. 3/23/19

      IF ( NLEVEL .EQ. 0 ) THEN
         DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
               QSTOKES_F(UTA,I,O1) = ZERO
            ENDDO
            IF ( DO_INCLUDE_TOAFLUX ) THEN
               O1 = 1 ; FLUX = FLUX_MULTIPLIER * TOAFLUX
               QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX
            ENDIF
         ENDDO
         RETURN
      ENDIF

!  Other levels, Thermal transmittance-only solution

      IF ( NLEVEL .NE. 0 .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = ZERO
          DO K = 1, N
            TPROP = T_WLOWER(I,K)/QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  Other levels, scattering solution
!  ---------------------------------

      IF ( NLEVEL .NE. 0 ) THEN

!  start main loops

        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
              MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
              HOM1 = LXR * T_DELT_EIGEN(K,N)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

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

!  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WLOWER(I,O1,N)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish streams/stokes loops

          ENDDO
        ENDDO
       
!  Version 2.8.1, Add Transmittance of TOA flux. 3/23/19
       
        IF ( DO_INCLUDE_TOAFLUX ) THEN
          O1 = 1
          DO I = 1, NSTREAMS
            FLUX = TOAFLUX
            DO K = 1, N
              FLUX = FLUX * T_DELT_DISORDS(I,K)
            enddo
            QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * FLUX
          ENDDO
        ENDIF

!  End clause
        
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADINTENS_LEVEL_DN

!

      SUBROUTINE QUADINTENS_OFFGRID_UP ( &
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,            & ! Input flags
          DO_CLASSICAL_SOLUTION, DO_INCLUDE_BOAFLUX, N, UTA, UT, IBEAM,             & ! Input flags/indices
          NSTOKES, NSTREAMS, NSTKS_NSTRMS, NLAYERS, TAYLOR_ORDER, HAVE_MULT,        & ! Input numbers/Bookkeeping
          FLUX_MULTIPLIER, BOAFLUX, PARTAU_VERT, QUAD_STREAMS,                      & ! Input quad/Fluxes/delt
          BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                   & ! Input Beam for Greens
          T_DELT_DISORDS, T_DISORDS_UTUP, ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, & ! input Disords/Greens
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! Input Homog. solutions
          BVEC, LCON, MCON, T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,              & ! Input thermal and PI solutions
          QSTOKES_F, UT_GMULT_UP, UT_GMULT_DN  )                                      ! OUTPUT

!  1/31/21. Version 2.8.3. Added Green's function capability.
!    -- Controlled by flag DO_CLASSICAL_SOLUTION
!    -- Green;s function variables ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P now inputs
!    -- Require also Taylor Order, multiplier control (HAVE_MULT) and outputs UT_GMULT_UP, UT_GMULT_DN
!    -- Use BVEC input for the classical solution particular integral (Instead of WUPPER and T_WUPPER)

!   Version 2.8.1, Control for BOA illumination added, 3/23/19

      USE VLIDORT_PARS_m, Only : MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES, MAXLAYERS, MAXBEAMS, &
                                 MAX_PARTLAYERS, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  Flags and numbers

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_BOAFLUX ! New 3/23/19

!  1/31/21. Version 2.8.3. Greens function, add flag DO_CLASSICAL_SOLUTION

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Numbers. 

      INTEGER, INTENT (IN) ::          N, UTA, UT, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS

!  1/31/21. Version 2.8.3, add TAYLOR ORDER

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER

!  1/31/21. Version 2.8.3. Local flag for getting the multipliers
!    ==> if it comes in false, it goes out true !!!!!!!!!!!!!!

      LOGICAL, intent(inout)  ::       HAVE_MULT

!  Fluxes

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: BOAFLUX                       ! BOA Flux (new 3/23/19)

!  Quadratures

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Optical depths

      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3Initial/Layer transmittance factors for solar beams, and divided by user-cosines

      INTEGER         , INTENT (IN)   :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN)   :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  discrete ordinate trans.

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3. Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  RTE solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3. need BVEC for the classical solution articular integral.

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  BVP Constants

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  output
!  ------

!  Local Stokes vector

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  1/31/21. Version 2.8.3. Green functions multipliers for off-grid optical depths

      DOUBLE PRECISION, INTENT (INOUT) :: UT_GMULT_UP ( MAXEVALUES,MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UT_GMULT_DN ( MAXEVALUES,MAX_PARTLAYERS )

!  local variables
!  ---------------

      INTEGER ::          I, I1, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SHOM, HOM1, HOM2, SHOM_R, FLUX
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI
      DOUBLE PRECISION :: WX, ZW, CONST, EPS, SD, SU, SPAR

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I,O1)
          DO K = NLAYERS, N + 1, -1
            TPROP = T_WUPPER(I1,K) / QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          TPROP =  UT_T_PARTIC(I1,UT) / QUAD_STREAMS(I)
          THELP = THELP * T_DISORDS_UTUP(I,UT) + TPROP
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous solution

      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            LXR  = LCON(K,N)*SOLA_XPOS(I1,O1,K,N)
            MXR  = MCON(K,N)*SOLB_XNEG(I1,O1,K,N)
            HOM1 = LXR * T_UTDN_EIGEN(K,UT)
            HOM2 = MXR * T_UTUP_EIGEN(K,UT)
            SHOM_R = SHOM_R + HOM1 + HOM2
          ENDDO

!  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LXR_CR =  LCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - &
                      LCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
            LXR_CI =  LCON(K1,N) * SOLA_XPOS(I1,O1,K2,N) + &
                      LCON(K2,N) * SOLA_XPOS(I1,O1,K1,N)
            MXR_CR =  MCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - &
                      MCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
            MXR_CI =  MCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + &
                      MCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
            HOM1CR =    LXR_CR * T_UTDN_EIGEN(K1,UT) &
                      - LXR_CI * T_UTDN_EIGEN(K2,UT)
            HOM2CR =    MXR_CR * T_UTUP_EIGEN(K1,UT) &
                      - MXR_CI * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
          ENDDO

!  real part

          SHOM = SHOM_R + SHOM_CR
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * SHOM

!  Finish streams/stokes loops

        ENDDO
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = UT_T_PARTIC(I1, UT)
           QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  1/31/21. Version 2.8.3.  Green's function solution is newly implemented
!    -- No Beam solution Green's function if no source. Zero output for safety.

      IF ( N .GT. BEAM_CUTOFF(IBEAM) ) THEN
        UT_GMULT_DN = zero ; UT_GMULT_UP = zero ; RETURN
      ENDIF

!  1/31/21. Version 2.8.3. This section has been removed
!      IF ( DO_INCLUDE_THERMEMISS ) THEN
!        DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!            QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
!            SPAR = WUPPER(I1,O1,N)
!            IF ( O1.eq.1) SPAR = SPAR - T_WUPPER(I1,N)
!            SPAR = SPAR * T_UTDN_MUBAR(UT,IBEAM)
!            QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
!          ENDDO
!        ENDDO
!      ELSE
!        DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!            SPAR = WUPPER(I1,O1,N) * T_UTDN_MUBAR(UT,IBEAM)
!            QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
!          ENDDO
!        ENDDO
!      ENDIF

!  1/31/21. Version 2.8.3. Classical Solution
!    -- Replaces above commented-out code, Use BVEC directly.

      IF ( DO_CLASSICAL_SOLUTION ) THEN
         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
               SPAR = BVEC(I1,O1,N) * INITIAL_TRANS(N,IBEAM) * T_UTDN_MUBAR(UT,IBEAM)
               QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
            ENDDO
         ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. Green's function solution
!    ==> First Compute output multipliers UT_GMULT_DN/UT_GMULT_UP, then add to result

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN
         IF ( .not. HAVE_MULT ) THEN
            UT_GMULT_DN(:,UT) = ZERO ; UT_GMULT_UP(:,UT) = ZERO
            WX    = T_UTDN_MUBAR(UT,IBEAM)
            CONST = INITIAL_TRANS(N,IBEAM)
            DO K = 1, K_REAL(N)
               ZW = T_DELT_MUBAR(N,IBEAM) * T_UTUP_EIGEN(K,UT)
               SU =  ( WX - ZW ) / GAMMA_P(K,N)
               IF ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
                 EPS = GAMMA_M(K,N)
                 CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAU_VERT(UT), WX, ONE, SD )
               ELSE
                 SD =  ( T_UTDN_EIGEN(K,UT) - WX ) / GAMMA_M(K,N)
               ENDIF
               UT_GMULT_DN(K,UT) = UT_GMULT_DN(K,UT) + SD * ATERM_SAVE(K,N) * CONST
               UT_GMULT_UP(K,UT) = UT_GMULT_UP(K,UT) + SU * BTERM_SAVE(K,N) * CONST
            ENDDO
!mick debug - tempo turned off
            !HAVE_MULT = .true.
         ENDIF

!  Add the Green's function contributions

         DO O1 = 1, NSTOKES
            DO I = 1, NSTREAMS
               SPAR = DOT_PRODUCT ( UT_GMULT_DN(1:NSTKS_NSTRMS,UT),SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N) ) &
                    + DOT_PRODUCT ( UT_GMULT_UP(1:NSTKS_NSTRMS,UT),SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N) )
               IF ( O1.gt.2) SPAR = -SPAR
               QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
            ENDDO
         ENDDO
      ENDIF

!  Version 2.8a, Add Transmittance of BOA ISOTROPIC flux. 3/23/19
       
      IF ( DO_INCLUDE_BOAFLUX ) THEN
         O1 = 1
         DO I = 1, NSTREAMS
            FLUX = BOAFLUX
            DO K = NLAYERS, N, -1
               FLUX = FLUX * T_DELT_DISORDS(I,K)
            enddo
            FLUX = FLUX * T_DISORDS_UTUP(I,UT) 
            QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * FLUX
         ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADINTENS_OFFGRID_UP

!

      SUBROUTINE QUADINTENS_OFFGRID_DN ( &
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,            & ! Input flags
          DO_CLASSICAL_SOLUTION, DO_INCLUDE_TOAFLUX, N, UTA, UT, IBEAM,             & ! Input flags/indices
          NSTOKES, NSTREAMS, NSTKS_NSTRMS, TAYLOR_ORDER, HAVE_MULT,                 & ! Input numbers/Bookkeeping
          FLUX_MULTIPLIER, TOAFLUX, PARTAU_VERT, QUAD_STREAMS,                      & ! Input quad/Fluxes/delt
          BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                   & ! Input Beam for Greens
          T_DELT_DISORDS, T_DISORDS_UTDN, ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, & ! input Disords/Greens
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! Input Homog. solutions
          BVEC, LCON, MCON, T_WLOWER, UT_T_PARTIC,                                  & ! Input thermal and PI solutions
          QSTOKES_F, UT_GMULT_UP, UT_GMULT_DN  )                                      ! OUTPUT

!  1/31/21. Version 2.8.3. Added Green's function capability.
!    -- Controlled by flag DO_CLASSICAL_SOLUTION
!    -- Green;s function variables ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P now inputs
!    -- Require also Taylor Order, multiplier control (HAVE_MULT) and outputs UT_GMULT_UP, UT_GMULT_DN
!    -- Use BVEC input for the classical solution particular integral (Instead of WUPPER and T_WUPPER)

!   Version 2.8.1, Control for TOA illumination added, 3/23/19

      USE VLIDORT_PARS_m, Only : MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES, MAXLAYERS, MAXBEAMS, &
                                 MAX_PARTLAYERS, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, TAYLOR_SMALL, SDU

      IMPLICIT NONE

!  Flags

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_TOAFLUX ! New 3/23/19

!  1/31/21. Version 2.8.3. Greens function, add flag DO_CLASSICAL_SOLUTION

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Numbers.

      INTEGER, INTENT (IN) ::          N, UTA, UT, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS

!  1/31/21. Version 2.8.3. Add TAYLOR_ORDER

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER

!  1/31/21. Version 2.8.3. Local flag for getting the multipliers
!    ==> if it comes in false, it goes out true !!!!!!!!!!!!!!

      LOGICAL, intent(inout)  ::       HAVE_MULT

!  Fluxes. TOA Flux (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: TOAFLUX   ! new 3/23/19

!  quadrature.

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Optical depths

      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3.
!    -- Initial/Layer transmittance factors for solar beams, and divided by user-cosines

      INTEGER         , INTENT (IN)   :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN)   :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN)   :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Discrete ordinate trans.

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3. 
!    -- Saved quantities from the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  RTE solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3. Use BVEC for the classical solution (replaces use of WUPPER)

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  BVP constants of integration

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS )

!  output
!  ------

!  Local Stokes vector

      DOUBLE PRECISION, INTENT (INOUT) :: QSTOKES_F ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  1/31/21. Version 2.8.3. Green functions multipliers for off-grid optical depths

      DOUBLE PRECISION, INTENT (INOUT) :: UT_GMULT_UP ( MAXEVALUES,MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UT_GMULT_DN ( MAXEVALUES,MAX_PARTLAYERS )

!  local variables
!  ---------------

      INTEGER ::          I, I1, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP, FLUX
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI
      DOUBLE PRECISION :: WX, ZW, CONST, EPS, SD, SU, SPAR

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = ZERO
          DO K = 1, N - 1
            TPROP = T_WLOWER(I,K) / QUAD_STREAMS(I)
            THELP = THELP*T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          TPROP = UT_T_PARTIC(I,UT) / QUAD_STREAMS(I)
          THELP = THELP * T_DISORDS_UTDN(I,UT) + TPROP
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
!mick fix - 2/10/2011
          QSTOKES_F(UTA,I,2:MAXSTOKES) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous solution

      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
            MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
            HOM1 = LXR * T_UTDN_EIGEN(K,UT)
            HOM2 = MXR * T_UTUP_EIGEN(K,UT)
            SHOM_R = SHOM_R + HOM1 + HOM2
          ENDDO

!  complex homogeneous solutions

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
            MXR_CI =  MCON(K1,N) * SOLB_XNEG(I,O1,K2,N) + &
                      MCON(K2,N) * SOLB_XNEG(I,O1,K1,N)
            HOM1CR =    LXR_CR * T_UTDN_EIGEN(K1,UT) &
                      - LXR_CI * T_UTDN_EIGEN(K2,UT)
            HOM2CR =    MXR_CR * T_UTUP_EIGEN(K1,UT) &
                      - MXR_CI * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
          ENDDO

!  real part

          SHOM = SHOM_R + SHOM_CR
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * SHOM

!  Finish streams/stokes loops

        ENDDO
      ENDDO

!  Add the thermal particular solution  (if flagged)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          SPAR = UT_T_PARTIC(I,UT)
           QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  1/31/21. Version 2.8.3.  Green's function solution is newly implemented
!    -- No Beam solution Green's function if no source. Zero output for safety.

      IF ( N .GT. BEAM_CUTOFF(IBEAM) ) THEN
        UT_GMULT_DN = zero ; UT_GMULT_UP = zero ; RETURN
      ENDIF

!  1/31/21. Version 2.8.3. This section has been removed
!  Add the solar particular solution, in the presence of thermal sources or not.
!   -- Logic has been changed for Version 2.8a. 3/23/19      
!      if ( DO_SOLAR_SOURCES ) THEN
!         IF ( DO_INCLUDE_THERMEMISS) THEN
!            DO I = 1, NSTREAMS ; DO O1 = 1, NSTOKES
!               SPAR = WUPPER(I,O1,N) ; IF ( O1.eq.1) SPAR = SPAR - T_WUPPER(I,N)
!               QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR * T_UTDN_MUBAR(UT,IBEAM)
!            ENDDO ; ENDDO
!         ELSE
!            DO I = 1, NSTREAMS ; DO O1 = 1, NSTOKES
!               SPAR = WUPPER(I,O1,N) * T_UTDN_MUBAR(UT,IBEAM)
!               QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
!            ENDDO ; ENDDO
!         ENDIF
!      ENDIF
       
!  1/31/21. Version 2.8.3. Classical Solution
!    -- Replaces above commented-out code, Use BVEC directly.

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            SPAR = BVEC(I,O1,N) * INITIAL_TRANS(N,IBEAM) * T_UTDN_MUBAR(UT,IBEAM)
            QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
          ENDDO
        ENDDO
      ENDIF
  
!  1/31/21. Version 2.8.3.. Green's function solution
!    -- First Compute output multipliers UT_GMULT_DN/UT_GMULT_UP, then add to result

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN
         IF ( .not. HAVE_MULT ) THEN
            UT_GMULT_DN(:,UT) = ZERO ; UT_GMULT_UP(:,UT) = ZERO
            WX    = T_UTDN_MUBAR(UT,IBEAM)
            CONST = INITIAL_TRANS(N,IBEAM)
            DO K = 1, K_REAL(N)
               ZW = T_DELT_MUBAR(N,IBEAM) * T_UTUP_EIGEN(K,UT)
               SU =  ( WX - ZW ) / GAMMA_P(K,N)
               IF ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
                 EPS = GAMMA_M(K,N)
                 CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAU_VERT(UT), WX, ONE, SD )
               ELSE
                 SD =  ( T_UTDN_EIGEN(K,UT) - WX ) / GAMMA_M(K,N)
               ENDIF
               UT_GMULT_DN(K,UT) = SD * ATERM_SAVE(K,N) * CONST
               UT_GMULT_UP(K,UT) = SU * BTERM_SAVE(K,N) * CONST
            ENDDO
!mick debug - tempo turned off
            !HAVE_MULT = .true.
         ENDIF

!  Add the Green's function contributions

         DO O1 = 1, NSTOKES
            DO I = 1, NSTREAMS
               SPAR = DOT_PRODUCT ( UT_GMULT_UP(1:NSTKS_NSTRMS,UT),SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N) ) &
                    + DOT_PRODUCT ( UT_GMULT_DN(1:NSTKS_NSTRMS,UT),SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N) )
               IF ( O1.gt.2) SPAR = -SPAR
               QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
            ENDDO
         ENDDO
      ENDIF

!  Version 2.8a, Add Transmittance of TOA ISOTROPIC flux. 3/23/19
       
      IF ( DO_INCLUDE_TOAFLUX ) THEN
         O1 = 1
         DO I = 1, NSTREAMS
            FLUX = TOAFLUX
            DO K = 1, N - 1
               FLUX = FLUX * T_DELT_DISORDS(I,K)
            enddo
            FLUX = FLUX * T_DISORDS_UTDN(I,UT) 
            QSTOKES_F(UTA,I,O1) = QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * FLUX
         ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADINTENS_OFFGRID_DN

!  End module

      END MODULE vlidort_PostProcessing_m

