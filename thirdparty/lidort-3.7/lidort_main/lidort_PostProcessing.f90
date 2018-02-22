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
! #  This Version :   3.7 F90                               #
! #  Release Date :   June 2014                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! #       Surface-leaving, BRDF Albedo-scaling     (3.7)    # 
! #       Taylor series, BBF Jacobians, ThreadSafe (3.7)    #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

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

module lidort_PostProcessing

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  Taylor series routines

   USE lidort_Taylor_m, only : TAYLOR_SERIES_1, TAYLOR_SERIES_2

public

contains


SUBROUTINE WHOLELAYER_STERM_UP &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,       & ! Input
             IB, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,              & ! Input
             LAYER_PIS_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_P, DELTAU_VERT, & ! input
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,  & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WPOS1,          & ! Input
             LAYER_TSUP_UP, LCON, MCON, HMULT_1, HMULT_2, EMULT_UP,    & ! Input
             PMULT_UU, PMULT_UD, LAYER_SOURCE)                           ! Output

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAXBEAMS, &
                              ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Indices

      INTEGER  , intent(in)  :: N, IB

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  -------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined streams, aaverage-secant streams

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  RTE Solution inputs
!  -------------------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  solution multipliers (homogeneous, single-scatter)

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)

!  output
!  ------

!  Layer source terms

      REAL(fpk), intent(out) :: LAYER_SOURCE ( MAX_USER_STREAMS )
!      REAL(fpk), intent(out) :: MSCAT_LAYER_SOURCE ( MAX_USER_STREAMS )

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) ::  PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) ::  PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM, LUM
      REAL(fpk)  :: SPAR, SHOM, SFOR1, TM
      REAL(fpk)  :: WDEL, ITRANS, ITRANSWDEL, SD, SU
      REAL(fpk)  :: EPS, YFAC, FAC1, FAC2, MULT

!  No layer source term if no scattering in the layer
!   Very important to zero both output terms (bug solved 04/20/05)

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         LAYER_SOURCE(UM)       = ZERO
!         MSCAT_LAYER_SOURCE(UM) = ZERO
      ENDDO
      IF ( .NOT. SOURCETERM_FLAG ) RETURN

!  Homogeneous solutions
!  ---------------------

!  Only if scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
               LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
               MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
               SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N) + &
                             MCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N)
            ENDDO
            LAYER_SOURCE(UM) = SHOM
         ENDDO
      ENDIF

!  thermal emission term (direct and diffuse)
!  ------------------------------------------

      IF ( DO_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + LAYER_TSUP_UP(UM,N)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

!  Get the basic multipliers
!    mick/rob fix 9/12/13 - Taylor-series limiting case: PMULT_UD, AverageSecant(IB)--> Eigenvalue(AA)

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO AA = 1, NSTREAMS
            if ( ABS(GAMMA_M(AA,N)) .lt. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(AA,N)     ; FAC1 = ONE
               YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL * T_DELT_USERM(N,UM)
               CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAU_VERT(N), FAC1, FAC2, ONE, MULT )
               SD = ITRANS_USERM (N,LUM,IB) * MULT 
            else
               SD = ( ITRANS * HMULT_2(AA,UM,N) - EMULT_UP(LUM,N,IB) ) / GAMMA_M(AA,N)
            endif
            SU = ( ITRANSWDEL * HMULT_1(AA,UM,N) + EMULT_UP(LUM,N,IB) ) / GAMMA_P(AA,N)
            PMULT_UD(AA,UM,N) = SD * ATERM_SAVE(AA,N)
            PMULT_UU(AA,UM,N) = SU * BTERM_SAVE(AA,N)
         ENDDO
      ENDDO

!  Add contributions to the Green's function solution

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XPOS(UM,AA,N)*PMULT_UD(AA,UM,N) &
                        + U_XNEG(UM,AA,N)*PMULT_UU(AA,UM,N)
         ENDDO
         LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SPAR
      ENDDO

!  Options for adding the single scatter part of the beam solution
!  .. Full radiance mode, add single scatter part

      IF ( DO_MSMODE_LIDORT ) THEN
!        IF ( SAVE_LAYER_MSST ) THEN
!          DO UM = 1, N_USER_STREAMS
!            MSCAT_LAYER_SOURCE(UM) = LAYER_SOURCE(UM)
!          ENDDO
!        ENDIF
      ELSE
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SFOR1 = U_WPOS1(UM,N) * EMULT_UP(LUM,N,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SFOR1
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE WHOLELAYER_STERM_UP

!

SUBROUTINE WHOLELAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,       & ! Input
             IB, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,              & ! Input
             LAYER_PIS_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_M, DELTAU_VERT, & ! input
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,  & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WNEG1,          & ! Input
             LAYER_TSUP_DN, LCON, MCON, HMULT_1, HMULT_2, EMULT_DN,    & ! Input
             PMULT_DU, PMULT_DD, LAYER_SOURCE)                           ! Output

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAXBEAMS, &
                              ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Indices

      INTEGER  , intent(in)  :: N, IB

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  -------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_M      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined streams, aaverage-secant streams

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  RTE Solution inputs
!  -------------------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  solution multipliers (homogeneous, single-scatter)

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  output
!  ------

!  Layer source terms

      REAL(fpk), intent(out) :: LAYER_SOURCE ( MAX_USER_STREAMS )
!      REAL(fpk), intent(out) :: MSCAT_LAYER_SOURCE ( MAX_USER_STREAMS )

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) ::  PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) ::  PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM, LUM
      REAL(fpk)  :: SPAR, SHOM, SFOR1, TM
      REAL(fpk)  :: WDEL, ITRANS, ITRANSWDEL, SD, SU
      REAL(fpk)  :: EPS, YFAC, FAC1, FAC2, MULT

!  No layer source term if no scattering in the layer
!   Very important to zero both output terms (bug solved 04/20/05)

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         LAYER_SOURCE(UM)       = ZERO
      ENDDO
      IF ( .NOT. SOURCETERM_FLAG ) RETURN

!  Homogeneous solutions
!  ---------------------

!  Only if scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
              MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
              SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N) + &
                            MCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N)
            ENDDO
            LAYER_SOURCE(UM) = SHOM
         ENDDO
      ENDIF

!  thermal emission term (direct and diffuse)
!  ------------------------------------------

      IF ( DO_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + LAYER_TSUP_DN(UM,N)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

!  Get the basic multipliers and store them
!    mick/rob fix 9/12/13 - Taylor-series limiting case: PMULT_DD, AverageSecant(IB)--> Eigenvalue(AA)

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO AA = 1, NSTREAMS
            if ( ABS(GAMMA_M(AA,N)) .lt. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(AA,N)     ; FAC1 = T_DELT_USERM(N,UM)
               YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = WDEL
               CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAU_VERT(N), FAC1, FAC2, ONE, MULT )
               SD = ITRANS_USERM (N,LUM,IB) * MULT
            else
               SD = ( ITRANS * HMULT_1(AA,UM,N) - EMULT_DN(LUM,N,IB) ) / GAMMA_M(AA,N)
            endif
            SU = ( ITRANSWDEL * HMULT_2(AA,UM,N) + EMULT_DN(LUM,N,IB) ) / GAMMA_P(AA,N)
            PMULT_DD(AA,UM,N) = SD * ATERM_SAVE(AA,N)
            PMULT_DU(AA,UM,N) = SU * BTERM_SAVE(AA,N)
         ENDDO
      ENDDO

!  Add contributions to the Green's function solution

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XNEG(UM,AA,N)*PMULT_DD(AA,UM,N) &
                        + U_XPOS(UM,AA,N)*PMULT_DU(AA,UM,N)
         ENDDO
         LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SPAR
      ENDDO

!  Options
!  .. If operating in Ms-mode only, copy multiple scatter term
!  .. Full radiance mode, add single scatter part

      IF ( DO_MSMODE_LIDORT ) THEN
!        IF ( SAVE_LAYER_MSST ) THEN
!          DO UM = 1, N_USER_STREAMS
!            MSCAT_LAYER_SOURCE(UM) = LAYER_SOURCE(UM)
!          ENDDO
!        ENDIF
      ELSE
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SFOR1 = U_WNEG1(UM,N) * EMULT_DN(LUM,N,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SFOR1
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE WHOLELAYER_STERM_DN

!

SUBROUTINE PARTLAYER_STERM_UP &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                    & ! Input
             IB, UT, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,                       & ! Input
             LAYER_PIS_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_P, DELTAUS, PARTAUS,         & ! input
             INITIAL_TRANS, ITRANS_USERM, T_UTUP_USERM, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WPOS1, LAYER_TSUP_UTUP,      & ! Input
             LCON, MCON, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,                     & ! Input
             UT_PMULT_UU, UT_PMULT_UD, LAYER_SOURCE)                                  ! Output

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                              ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Indices

      INTEGER  , intent(in)  :: N, IB, UT

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  -------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAUS ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAUS ( MAX_PARTLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined streams, aaverage-secant streams

      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  RTE Solution Inputs
!  -------------------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_EMULT_UP(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  output
!  ------

!  source term

      REAL(fpk), intent(out)  :: LAYER_SOURCE(MAX_USER_STREAMS)

!  Green function multipliers, post processed

      REAL(fpk), intent(inout)  ::  UT_PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  ::  UT_PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM, LUM
      REAL(fpk)  :: SPAR, SHOM, SFOR, TM
      REAL(fpk)  :: WDEL, ITRANS, ITRANSWDEL, SD, SU
      REAL(fpk)  :: EPS, YFAC, FAC1, FAC2, MULT1, MULT2

!  No layer source term if no scattering in the layer

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         LAYER_SOURCE(UM)       = ZERO
      ENDDO
      if ( .NOT. SOURCETERM_FLAG ) RETURN

!  homogeneous solutions
!  ---------------------

!  Must be scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
              MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
              SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_UD(AA,UM,UT) + &
                            MCON_UXVEC(UM,AA)*UT_HMULT_UU(AA,UM,UT)
            ENDDO
            LAYER_SOURCE(UM) = SHOM
         ENDDO
      ENDIF

!  Add thermal term (direct and diffuse)
!  -------------------------------------

      IF ( DO_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM)+LAYER_TSUP_UTUP(UM,UT)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

!  Get the multipliers and store them

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO AA = 1, NSTREAMS
            if ( ABS(GAMMA_M(AA,N)) .lt. TAYLOR_SMALL ) THEN
              EPS   = GAMMA_M(AA,N)     ; FAC1 = - T_UTDN_MUBAR(UT,IB) 
              YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL * T_UTUP_USERM(UT,UM) 
              CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTAUS(UT), ZERO, FAC1, ONE, MULT1 )
              CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAUS(N),  ZERO, FAC2, ONE, MULT2 )
              SD = ITRANS_USERM (N,LUM,IB) * ( MULT1 + MULT2 )
            else
              SD = ( ITRANS * UT_HMULT_UD(AA,UM,UT) - UT_EMULT_UP(LUM,UT,IB) ) / GAMMA_M(AA,N)
            endif
            SU = ( ITRANSWDEL * UT_HMULT_UU(AA,UM,UT) + UT_EMULT_UP(LUM,UT,IB) ) / GAMMA_P(AA,N)
            UT_PMULT_UD(AA,UM,UT) = SD * ATERM_SAVE(AA,N)
            UT_PMULT_UU(AA,UM,UT) = SU * BTERM_SAVE(AA,N)
         ENDDO
      ENDDO

!  Add contributions to the Green's function solution

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XPOS(UM,AA,N)*UT_PMULT_UD(AA,UM,UT) &
                        + U_XNEG(UM,AA,N)*UT_PMULT_UU(AA,UM,UT)
         ENDDO
         LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SPAR
      ENDDO

!  If NOT operating in MS-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SFOR = U_WPOS1(UM,N) * UT_EMULT_UP(LUM,UT,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SFOR
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE PARTLAYER_STERM_UP

!

SUBROUTINE PARTLAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                    & ! Input
             IB, UT, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,                       & ! Input
             LAYER_PIS_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_M, PARTAUS,                  & ! input
             INITIAL_TRANS, ITRANS_USERM, T_UTDN_USERM, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WNEG1, LAYER_TSUP_UTDN,      & ! Input
             LCON, MCON, UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN,                     & ! Input
             UT_PMULT_DU, UT_PMULT_DD, LAYER_SOURCE)                                  ! Output

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                              ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Indices

      INTEGER  , intent(in)  :: N, IB, UT

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  -------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_M      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: PARTAUS ( MAX_PARTLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined streams, aaverage-secant streams

      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  RTE Solution Inputs
!  -------------------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  ::  ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  ::  BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  ::  U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  ::  U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  ::  U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  ::  LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Multipliers

      REAL(fpk), intent(in)  ::  UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  ::  UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  ::  UT_EMULT_DN(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  output
!  ------

!  source term

      REAL(fpk), intent(out)  :: LAYER_SOURCE(MAX_USER_STREAMS)

!  Green function multipliers, post processed solution

      REAL(fpk), intent(inout)  ::  UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  ::  UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM, LUM
      REAL(fpk)  :: SPAR, SHOM, SFOR, TM
      REAL(fpk)  :: WDEL, ITRANS, ITRANSWDEL, SD, SU
      REAL(fpk)  :: EPS, YFAC, FAC1, FAC2, MULT

!  No layer source term if no scattering in the layer

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         LAYER_SOURCE(UM) = ZERO
      ENDDO
      IF ( .not. SOURCETERM_FLAG ) RETURN

!  homogeneous solutions
!  ---------------------

!  Must be scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
              MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
              SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_DD(AA,UM,UT) + &
                            MCON_UXVEC(UM,AA)*UT_HMULT_DU(AA,UM,UT)
            ENDDO
            LAYER_SOURCE(UM) = SHOM
         ENDDO
      ENDIF

!  Add thermal term (direct and diffuse)
!  -------------------------------------

      IF ( DO_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM)+LAYER_TSUP_UTDN(UM,UT)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  Nothing further to do if no particular solution

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

!  Get the multipliers and store them

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO AA = 1, NSTREAMS
            if ( ABS(GAMMA_M(AA,N)) .lt. TAYLOR_SMALL ) THEN
              EPS   = GAMMA_M(AA,N)     ; FAC1 = T_UTDN_USERM(UT,UM) 
              YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = T_UTDN_MUBAR(UT,IB)
              CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTAUS(UT), FAC1, FAC2, ONE, MULT )
              SD = ITRANS_USERM (N,LUM,IB) * MULT
            else
              SD = ( ITRANS * UT_HMULT_DD(AA,UM,UT) - UT_EMULT_DN(LUM,UT,IB) ) / GAMMA_M(AA,N)
            endif
            SU = ( ITRANSWDEL * UT_HMULT_DU(AA,UM,UT) + UT_EMULT_DN(LUM,UT,IB) ) / GAMMA_P(AA,N)
            UT_PMULT_DD(AA,UM,UT) = SD * ATERM_SAVE(AA,N)
            UT_PMULT_DU(AA,UM,UT) = SU * BTERM_SAVE(AA,N)
         ENDDO
      ENDDO

!  Add contributions to the Greens function solution

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XNEG(UM,AA,N)*UT_PMULT_DD(AA,UM,UT) &
                        + U_XPOS(UM,AA,N)*UT_PMULT_DU(AA,UM,UT)
         ENDDO
         LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SPAR
      ENDDO

!  If NOT operating in MS-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SFOR = U_WNEG1(UM,N) * UT_EMULT_DN(LUM,UT,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SFOR
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE PARTLAYER_STERM_DN

!

SUBROUTINE QUADINTENS_LEVEL_UP                                        &
          ( DO_THERMAL_TRANSONLY, IPARTIC, UTA, NLEVEL, NSTREAMS,     & ! input
            NLAYERS, FLUX_MULTIPLIER, QUAD_STREAMS,                   & ! input
            T_DELT_DISORDS, BOA_THTONLY_SOURCE, T_WUPPER,             & ! input
            LCON_XVEC, MCON_XVEC,  WUPPER, WLOWER, T_DELT_EIGEN,      & ! input
            QUADINTENS )                                                ! output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAXLAYERS, MAXSTREAMS_2,   &
                              MAXSTREAMS,      MAXBEAMS,  MAX_DIRECTIONS, &
                              ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices

      INTEGER  , intent(in)  :: NLEVEL, UTA, IPARTIC

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)

!  Thermal transmittance-only source

      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Upper/Lower boundary

      REAL(fpk), intent(in)  :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
         (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  local variables
!  ---------------

      INTEGER    :: N, I, I1, AA, K
      REAL(fpk)  :: FM, SPAR, SHOM, HOM1, HOM2, THELP, QUAD

!  For those optical depths at layer boundaries
!  --------------------------------------------

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
!  looking at the intensity at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling intensity
!  at the bottom of the atmosphere (treated separately).

      N = NLEVEL + 1
      FM = FLUX_MULTIPLIER

!  homogeneous and particular solution contributions SHOM and SPAR

!  For the lowest level

      IF ( NLEVEL .EQ. NLAYERS ) THEN

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            THELP = BOA_THTONLY_SOURCE(I)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I1,AA,NLEVEL) * T_DELT_EIGEN(AA,NLEVEL)
              HOM2 = MCON_XVEC(I1,AA,NLEVEL)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            SPAR = WLOWER(I1,NLEVEL)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

!  For other levels in the atmosphere

      ELSE

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            QUAD = QUAD_STREAMS(I)
            THELP = BOA_THTONLY_SOURCE(I)
            DO K = NLAYERS, N, -1
              THELP = THELP*T_DELT_DISORDS(I,K) + T_WUPPER(I1,K)/QUAD
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I1,AA,N)
              HOM2 = MCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            SPAR = WUPPER(I1,N)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_LEVEL_UP

!

SUBROUTINE QUADINTENS_LEVEL_DN                                   &
          ( DO_THERMAL_TRANSONLY, IPARTIC, UTA, NLEVEL, NSTREAMS,& ! input
            FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,       & ! input
            T_WLOWER, LCON_XVEC, MCON_XVEC, WLOWER, T_DELT_EIGEN,& ! input
            QUADINTENS )                                           ! output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAXLAYERS, MAXSTREAMS_2,   &
                              MAXSTREAMS,      MAXBEAMS,  MAX_DIRECTIONS, &
                              ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices

      INTEGER  , intent(in)  :: NLEVEL, UTA, IPARTIC

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations 

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  local variables
!  ---------------

      INTEGER    :: N, I, AA, K
      REAL(fpk)  :: FM, SPAR, SHOM, HOM1, HOM2, QUAD, THELP

!  For those optical depths at layer boundaries
!  --------------------------------------------

      N = NLEVEL
      FM = FLUX_MULTIPLIER

!  Downwelling radiation at TOA ( or N = 0 ) is zero

      IF ( NLEVEL .EQ. 0 ) THEN

        DO I = 1, NSTREAMS
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = ZERO
        ENDDO

!  Other levels

      ELSE

!  Thermal transmittance solution, build from TOA downwards
!  Scattering solution, use the Discrete Ordinate solution

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            THELP = ZERO
            QUAD = QUAD_STREAMS(I)
            DO K = 1, N
              THELP = THELP*T_DELT_DISORDS(I,K) + T_WLOWER(I,K)/QUAD
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,DNIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            SPAR = WLOWER(I,N)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
              HOM2 = MCON_XVEC(I,AA,N)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,DNIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_LEVEL_DN

!

SUBROUTINE QUADINTENS_OFFGRID_UP &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! input
            IPARTIC, UTA, UT, N, NLAYERS, NSTREAMS, TAYLOR_ORDER,       & ! input
            FLUX_MULTIPLIER, INITIAL_TRANS, LAYER_PIS_CUTOFF,           & ! input
            QUAD_STREAMS, PARTAUS, T_UTUP_EIGEN, T_UTDN_EIGEN,          & ! input ! 5/6/13
            T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_DISORDS, T_DISORDS_UTUP, & ! input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                   & ! input
            T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,                  & ! input
            XPOS, LCON_XVEC, MCON_XVEC,                                 & ! input
            QUADINTENS, UT_GMULT_UP, UT_GMULT_DN )                        ! output

!Rob fix 5/6/13 - Add Partaus argument to QUADINTENS_OFFGRID_UP call

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAX_PARTLAYERS, MAXSTREAMS_2,   &
                              MAXLAYERS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, &
                              TAYLOR_SMALL, ZERO, ONE, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices

      INTEGER  , intent(in)  :: N, UTA, UT, IPARTIC

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!Rob Fix 5/6/13 - Add arguments
!  Input Optical depths required for Taylor-series limiting cases

      REAL(fpk), intent(in)  :: PARTAUS(MAX_PARTLAYERS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal transmittance-only source

      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed 3 outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(inout)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      INTEGER    :: I, I1, AA, K, IB
      REAL(fpk)  :: THELP, SPAR, PAR1, PAR2, SHOM, HOM1, HOM2, QUAD
      REAL(fpk)  :: WX, ZW, CONST, EPS, SD, SU

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          QUAD = QUAD_STREAMS(I)
          THELP = BOA_THTONLY_SOURCE(I)
          DO K = NLAYERS, N+1, -1
            THELP = THELP*T_DELT_DISORDS(I,K) + T_WUPPER(I1,K) / QUAD
          ENDDO
          THELP = THELP*T_DISORDS_UTUP(I,UT) + UT_T_PARTIC(I1,UT) / QUAD
          QUADINTENS(UTA,I,IPARTIC,UPIDX) = FLUX_MULTIPLIER*THELP
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          HOM1 = LCON_XVEC(I1,AA,N) * T_UTDN_EIGEN(AA,UT)
          HOM2 = MCON_XVEC(I1,AA,N) * T_UTUP_EIGEN(AA,UT)
          SHOM = SHOM + HOM1 + HOM2
        ENDDO
        QUADINTENS(UTA,I,IPARTIC,UPIDX) = FLUX_MULTIPLIER * SHOM
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = UT_T_PARTIC(I1, UT)
          QUADINTENS(UTA,I,IPARTIC,UPIDX) = &
            QUADINTENS(UTA,I,IPARTIC,UPIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  No Beam solution Green's function if no source. Zero output for safety.

      IF ( N .GT. LAYER_PIS_CUTOFF(IPARTIC) ) THEN
        UT_GMULT_DN = zero ; UT_GMULT_UP = zero ; RETURN
      ENDIF

!  Get the Multipliers for the partial solution.

      IB    = IPARTIC
      WX    = T_UTDN_MUBAR(UT,IB)
      CONST = INITIAL_TRANS(N,IB)
      DO AA = 1, NSTREAMS
        ZW    = T_DELT_MUBAR(N,IB) * T_UTUP_EIGEN(AA,UT)
        SU =  ( WX - ZW ) / GAMMA_P(AA,N)
        IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
           EPS = GAMMA_M(AA,N)
           CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAUS(UT), WX, ONE, SD )
        ELSE
           SD =  ( T_UTDN_EIGEN(AA,UT) - WX ) / GAMMA_M(AA,N)
        ENDIF
        UT_GMULT_DN(AA,UT) = SD * ATERM_SAVE(AA,N) * CONST
        UT_GMULT_UP(AA,UT) = SU * BTERM_SAVE(AA,N) * CONST
      ENDDO

!  Add the Green's function contributions

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          PAR1 = XPOS(I,AA,N)  * UT_GMULT_UP(AA,UT)
          PAR2 = XPOS(I1,AA,N) * UT_GMULT_DN(AA,UT)
          SPAR = SPAR + PAR1 + PAR2
        ENDDO
        QUADINTENS(UTA,I,IB,UPIDX) = QUADINTENS(UTA,I,IB,UPIDX) + FLUX_MULTIPLIER * SPAR
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_OFFGRID_UP

!

SUBROUTINE QUADINTENS_OFFGRID_DN &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! input
            IPARTIC, UTA, UT, N, NSTREAMS, TAYLOR_ORDER, DO_MULT,       & ! input
            FLUX_MULTIPLIER, INITIAL_TRANS, LAYER_PIS_CUTOFF,           & ! input
            QUAD_STREAMS, PARTAUS, T_UTUP_EIGEN, T_UTDN_EIGEN,          & ! input
            T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_DISORDS, T_DISORDS_UTDN, & ! input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                   & ! Input
            T_WLOWER, UT_T_PARTIC, XPOS, LCON_XVEC, MCON_XVEC,          & ! Input
            QUADINTENS, UT_GMULT_UP, UT_GMULT_DN )                        ! Output

!Rob fix 5/6/13 - Add Partaus argument to QUADINTENS_OFFGRID_DN call

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAX_PARTLAYERS, MAXSTREAMS_2,   &
                              MAXLAYERS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, &
                              TAYLOR_SMALL, ZERO, ONE, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices

      INTEGER  , intent(in)  :: N, UTA, UT, IPARTIC

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Local flag for getting the multipliers

      LOGICAL  , intent(in)  :: DO_MULT

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!Rob Fix 5/6/13 - Add arguments
!  Input Optical depths required for Taylor-series limiting cases

      REAL(fpk), intent(in)  :: PARTAUS(MAX_PARTLAYERS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed 3 outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Green functions multipliers for off-grid optical depths
!   These will only be generated as output if flagged

      REAL(fpk), intent(inout) :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      INTEGER    :: I, I1, AA, K, IB
      REAL(fpk)  :: THELP, SPAR, PAR1, PAR2, SHOM, HOM1, HOM2, QUAD
      REAL(fpk)  :: WX, ZW, CONST, EPS, SD, SU

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          QUAD = QUAD_STREAMS(I)
          THELP = ZERO
          DO K = 1, N-1
            THELP = THELP*T_DELT_DISORDS(I,K) + T_WLOWER(I,K) / QUAD
          ENDDO
          THELP = THELP*T_DISORDS_UTDN(I,UT) + UT_T_PARTIC(I,UT) / QUAD
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = FLUX_MULTIPLIER*THELP
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous

      DO I = 1, NSTREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          HOM1 = LCON_XVEC(I,AA,N) * T_UTDN_EIGEN(AA,UT)
          HOM2 = MCON_XVEC(I,AA,N) * T_UTUP_EIGEN(AA,UT)
          SHOM = SHOM + HOM1 + HOM2
        ENDDO
        QUADINTENS(UTA,I,IPARTIC,DNIDX) = FLUX_MULTIPLIER * SHOM
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          SPAR = UT_T_PARTIC(I,UT)
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = QUADINTENS(UTA,I,IPARTIC,DNIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  No Beam solution Green's function if no source.

      IF ( N .GT. LAYER_PIS_CUTOFF(IPARTIC) ) RETURN

!  Get the local Multipliers for the partial solution, only if not already obtained

      IF ( DO_MULT ) THEN
        IB    = IPARTIC
        WX    = T_UTDN_MUBAR(UT,IB)
        CONST = INITIAL_TRANS(N,IB)
        DO AA = 1, NSTREAMS
          ZW    = T_DELT_MUBAR(N,IB) * T_UTUP_EIGEN(AA,UT)
          SU =  ( WX - ZW ) / GAMMA_P(AA,N)
          IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
             EPS = GAMMA_M(AA,N)
             CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAUS(UT), WX, ONE, SD )
          ELSE
             SD =  ( T_UTDN_EIGEN(AA,UT) - WX ) / GAMMA_M(AA,N)
          ENDIF
          UT_GMULT_DN(AA,UT) = SD * ATERM_SAVE(AA,N) * CONST
          UT_GMULT_UP(AA,UT) = SU * BTERM_SAVE(AA,N) * CONST
        ENDDO
      ENDIF

!  Add the contributions to the Green's function solution

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          PAR1 = XPOS(I1,AA,N) * UT_GMULT_UP(AA,UT)
          PAR2 = XPOS(I,AA,N)  * UT_GMULT_DN(AA,UT)
          SPAR = SPAR + PAR1 + PAR2
        ENDDO
        QUADINTENS(UTA,I,IB,DNIDX) = QUADINTENS(UTA,I,IB,DNIDX) + FLUX_MULTIPLIER * SPAR
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_OFFGRID_DN

!  End Module

end module lidort_PostProcessing

