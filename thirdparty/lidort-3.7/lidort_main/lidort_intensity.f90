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
! #   -- Master subroutines for Radiances/Fluxs             #
! #                                                         #
! #          UPUSER_INTENSITY                               #
! #          DNUSER_INTENSITY                               #
! #          MIFLUX_INTENSITY                               #
! #                                                         #
! #   -- Master routines for post-processed TOA/BOA fields  #
! #                                                         #
! #          GET_TOASOURCE                                  #
! #          GET_BOASOURCE                                  #
! #                                                         #
! #   -- Master routines for azimuth Convergence            #
! #                                                         #
! #          LIDORT_CONVERGE                                #
! #          LIDORT_CONVERGE_OBSGEO                         #
! #                                                         #
! ###########################################################

module lidort_intensity

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  Parameter types

   USE LIDORT_PARS, only : fpk, ONE

!  Dependencies

   USE lidort_PostProcessing

!  Everything public

public 

contains

SUBROUTINE UPUSER_INTENSITY                                             &
           ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,       & ! Input
             DO_OBSERVATION_GEOMETRY, DO_LAYER_SCATTERING,              & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY, DO_INCLUDE_MVOUTPUT,  & ! Input
             FOURIER_COMPONENT, NSTREAMS, N_USER_STREAMS, NLAYERS,      & ! Input
             N_USER_LEVELS, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,    & ! Input
             UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                   & ! Input
             IPARTIC, TAYLOR_ORDER, FLUX_MULTIPLIER, QUAD_STREAMS,      & ! Input
             GAMMA_P, GAMMA_M, SIGMA_P,                                 & ! Input
             DELTAU_VERT, PARTAU_VERT,                                  & ! Input
             INITIAL_TRANS, ITRANS_USERM,                               & ! Input
             T_DELT_USERM, T_UTUP_USERM,                                & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR,                                & ! Input
             T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                  & ! Input
             T_DELT_DISORDS, T_DISORDS_UTUP,                            & ! Input
             T_WUPPER, UT_T_PARTIC, LAYER_TSUP_UP, LAYER_TSUP_UTUP,     & ! Input
             LAYER_PIS_CUTOFF, ATERM_SAVE, BTERM_SAVE,                  & ! Input
             XPOS, WUPPER, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,    & ! Input
             U_XPOS, U_XNEG, U_WPOS1, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
             UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                    & ! Input
             BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE,         & ! Input
             PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,              & ! Output
             FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                     & ! Output
             INTENSITY_F, QUADINTENS, CUMSOURCE_UP )                      ! Output

!  Upwelling post-processed Intensity Fourier component

!  dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
                              MAXSTREAMS_2, MAXMOMENTS, MAXBEAMS, MAX_USER_LEVELS,     &
                              MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_THERMEMISS

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Fourier component, beam index

      INTEGER  , intent(in)  :: FOURIER_COMPONENT
      INTEGER  , intent(in)  :: IPARTIC

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Flux multiplier

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  ------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT  ( MAX_PARTLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!   Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!   Transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )

!   Discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Upper/Lower boundary

      REAL(fpk), intent(in)  :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  BOA source terms

      REAL(fpk), intent(in)  :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout) :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  User-defined solutions

      REAL(fpk), intent(inout) :: INTENSITY_F &
        (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Cumulative source terms

      REAL(fpk), intent(inout) :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(inout) :: UT_PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_GMULT(MAX_PARTLAYERS)

      REAL(fpk), intent(inout) :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  help variables

      LOGICAL    :: SOURCETERM_FLAG
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER    :: UT, UTA, UM, LUM, M
      INTEGER    :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk)  :: LAYER_SOURCE ( MAX_USER_STREAMS )
!      REAL(fpk)  :: MSCAT_LAYER_SOURCE ( MAX_USER_STREAMS )
      REAL(fpk)  :: FINAL_SOURCE

!  Local post-processing control

      PPSTREAM_MASK(:,IPARTIC) = 0
      IF ( DO_OBSERVATION_GEOMETRY ) THEN
         N_PPSTREAMS = 1; PPSTREAM_MASK(1,IPARTIC) = IPARTIC
      else
         N_PPSTREAMS = N_USER_STREAMS
         do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IPARTIC) = UM
         enddo
      endif

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      M = FOURIER_COMPONENT
      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
           DO LUM = 1, N_PPSTREAMS
              INTENSITY_F(UTA,LUM,IPARTIC,UPIDX) = ZERO
           END DO
        ENDDO
      ENDIF

!  Initialize

      NUT = 0
      NC  = 0

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values

      IF ( DO_USER_STREAMS ) THEN
         NC = 0
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IPARTIC)
            CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM) + DIRECT_BOA_SOURCE(UM)
         ENDDO
      ENDIF

!  Recursion Loop in Source function integration
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART   = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            CALL WHOLELAYER_STERM_UP &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,         & ! Input
             IPARTIC, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,           & ! Input
             LAYER_PIS_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_P, DELTAU_VERT,   & ! input
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WPOS1,            & ! Input
             LAYER_TSUP_UP, LCON, MCON, HMULT_1, HMULT_2, EMULT_UP,      & ! Input
             PMULT_UU, PMULT_UD, LAYER_SOURCE)                             ! Output

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IPARTIC)
               CUMSOURCE_UP(UM,NC) = LAYER_SOURCE(UM) + &
                       T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
            ENDDO

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

!  Quadrature intensity calculation at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!Rob fix 5/6/13 - Add Partaus argument to QUADINTENS_OFFGRID_UP call

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            FLAGS_GMULT(UT) = .TRUE.
            CALL QUADINTENS_OFFGRID_UP &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! input
            IPARTIC, UTA, UT, N, NLAYERS, NSTREAMS, TAYLOR_ORDER,       & ! input
            FLUX_MULTIPLIER, INITIAL_TRANS, LAYER_PIS_CUTOFF,           & ! input
            QUAD_STREAMS, PARTAU_VERT, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! input ! 5/6/13
            T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_DISORDS, T_DISORDS_UTUP, & ! input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                   & ! input
            T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,                  & ! input
            XPOS, LCON_XVEC, MCON_XVEC,                                 & ! input
            QUADINTENS, UT_GMULT_UP, UT_GMULT_DN )                        ! output
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_UP &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                    & ! Input
             IPARTIC, UT, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,                  & ! Input
             LAYER_PIS_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_P, DELTAU_VERT, PARTAU_VERT, & ! input
             INITIAL_TRANS, ITRANS_USERM, T_UTUP_USERM, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WPOS1, LAYER_TSUP_UTUP,      & ! Input
             LCON, MCON, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,                     & ! Input
             UT_PMULT_UU, UT_PMULT_UD, LAYER_SOURCE)                                  ! Output

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IPARTIC)
               FINAL_SOURCE = LAYER_SOURCE(UM) + T_UTUP_USERM(UT,UM)*CUMSOURCE_UP(UM,NC)
               INTENSITY_F(UTA,LUM,IPARTIC,UPIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_LEVEL_UP                                 &
           ( DO_THERMAL_TRANSONLY, IPARTIC, UTA, NLEVEL, NSTREAMS,   & ! input
             NLAYERS, FLUX_MULTIPLIER, QUAD_STREAMS,                 & ! input
             T_DELT_DISORDS, BOA_THTONLY_SOURCE, T_WUPPER,           & ! input
             LCON_XVEC, MCON_XVEC,  WUPPER, WLOWER, T_DELT_EIGEN,    & ! input
             QUADINTENS )                                              ! output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IPARTIC)
                INTENSITY_F(UTA,LUM,IPARTIC,UPIDX) = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
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
END SUBROUTINE UPUSER_INTENSITY

!

SUBROUTINE DNUSER_INTENSITY                                            &
           ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,      & ! Input
             DO_OBSERVATION_GEOMETRY, DO_LAYER_SCATTERING,             & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY, DO_INCLUDE_MVOUTPUT, & ! Input
             FOURIER_COMPONENT, NSTREAMS, N_USER_STREAMS,              & ! Input
             N_USER_LEVELS, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,   & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                  & ! Input
             IPARTIC, TAYLOR_ORDER, FLUX_MULTIPLIER, QUAD_STREAMS,     & ! Input
             GAMMA_P, GAMMA_M, SIGMA_M,                                & ! Input
             DELTAU_VERT, PARTAU_VERT,                                 & ! Input
             INITIAL_TRANS, ITRANS_USERM,                              & ! Input
             T_DELT_USERM, T_UTDN_USERM,                               & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR,                               & ! Input
             T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                 & ! Input
             T_DELT_DISORDS, T_DISORDS_UTDN,                           & ! Input
             T_WLOWER, UT_T_PARTIC, LAYER_TSUP_DN, LAYER_TSUP_UTDN,    & ! Input
             LAYER_PIS_CUTOFF, ATERM_SAVE, BTERM_SAVE,                 & ! Input----
             XPOS, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,           & ! Input
             U_XPOS, U_XNEG, U_WNEG1, HMULT_1, HMULT_2, EMULT_DN,      & ! Input
             UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN, TOA_SOURCE,       & ! Input
             PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,             & ! Output
             FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                    & ! Output
             INTENSITY_F, QUADINTENS, CUMSOURCE_DN )                     ! Output

!  Downwelling post-processed Intensity Fourier component

!  Dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
                              MAXSTREAMS_2, MAXMOMENTS, MAXBEAMS, MAX_USER_LEVELS,     &
                              MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_THERMEMISS

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Fourier component, beam index

      INTEGER  , intent(in)  :: FOURIER_COMPONENT
      INTEGER  , intent(in)  :: IPARTIC

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  multipliers

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  ------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_M      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT  ( MAX_PARTLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!   Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!   Transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )

!   Discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_DN(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  TOA source terms

      REAL(fpk), intent(in)  :: TOA_SOURCE ( MAX_USER_STREAMS )

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  User-defined solutions

      REAL(fpk), intent(inout)  :: INTENSITY_F &
        (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Cumulative source terms

      REAL(fpk), intent(inout) :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(inout) :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_GMULT(MAX_PARTLAYERS)

      REAL(fpk), intent(inout) :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      LOGICAL    :: LOCAL_GMULT, SOURCETERM_FLAG
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER    :: UT, UTA, UM, NC, LUM, M
      INTEGER    :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk)  :: LAYER_SOURCE ( MAX_USER_STREAMS )
!      REAL(fpk)  :: MSCAT_LAYER_SOURCE ( MAX_USER_STREAMS )
      REAL(fpk)  :: FINAL_SOURCE

!  Local post-processing control

      PPSTREAM_MASK(:,IPARTIC) = 0
      IF ( DO_OBSERVATION_GEOMETRY ) THEN
         N_PPSTREAMS = 1; PPSTREAM_MASK(1,IPARTIC) = IPARTIC
      else
         N_PPSTREAMS = N_USER_STREAMS
         do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IPARTIC) = UM
         enddo
      endif

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      M = FOURIER_COMPONENT
      IF ( DO_USER_STREAMS ) THEN
         DO UTA = 1, N_USER_LEVELS
            DO LUM = 1, N_PPSTREAMS
               INTENSITY_F(UTA,LUM,IPARTIC,DNIDX) = ZERO
            ENDDO
         ENDDO
      ENDIF

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
         NC = 0
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IPARTIC)
            CUMSOURCE_DN(UM,NC) = TOA_SOURCE(UM)
         ENDDO
      ENDIF

!  Initialize

      NUT = 0
      NC  = 0

!  initialise cumulative source term loop

      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            NC = N
            SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

            CALL WHOLELAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,         & ! Input
             IPARTIC, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,           & ! Input
             LAYER_PIS_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_M, DELTAU_VERT,   & ! input
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WNEG1,            & ! Input
             LAYER_TSUP_DN, LCON, MCON, HMULT_1, HMULT_2, EMULT_DN,      & ! Input
             PMULT_DU, PMULT_DD, LAYER_SOURCE)                             ! Output

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IPARTIC)
               CUMSOURCE_DN(UM,NC) = LAYER_SOURCE(UM) + &
                       T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
            ENDDO

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!Rob fix 5/6/13 - Add Partaus argument to QUADINTENS_OFFGRID_DN call

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            LOCAL_GMULT = FLAGS_GMULT(UT)
            CALL QUADINTENS_OFFGRID_DN &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! input
            IPARTIC, UTA, UT, N, NSTREAMS, TAYLOR_ORDER, LOCAL_GMULT,   & ! input
            FLUX_MULTIPLIER, INITIAL_TRANS, LAYER_PIS_CUTOFF,           & ! input
            QUAD_STREAMS, PARTAU_VERT, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! input
            T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_DISORDS, T_DISORDS_UTDN, & ! input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                   & ! Input
            T_WLOWER, UT_T_PARTIC, XPOS, LCON_XVEC, MCON_XVEC,          & ! Input
            QUADINTENS, UT_GMULT_UP, UT_GMULT_DN )                        ! Output
            FLAGS_GMULT(UT) = .FALSE.
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                    & ! Input
             IPARTIC, UT, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,                  & ! Input
             LAYER_PIS_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_M, PARTAU_VERT,              & ! input
             INITIAL_TRANS, ITRANS_USERM, T_UTDN_USERM, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WNEG1, LAYER_TSUP_UTDN,      & ! Input
             LCON, MCON, UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN,                     & ! Input
             UT_PMULT_DU, UT_PMULT_DD, LAYER_SOURCE)                                  ! Output

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IPARTIC)
               FINAL_SOURCE = LAYER_SOURCE(UM) + T_UTDN_USERM(UT,UM)*CUMSOURCE_DN(UM,NC)
               INTENSITY_F(UTA,LUM,IPARTIC,DNIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_LEVEL_DN                               &
            ( DO_THERMAL_TRANSONLY, IPARTIC, UTA, NLEVEL, NSTREAMS,& ! input
              FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,       & ! input
              T_WLOWER, LCON_XVEC, MCON_XVEC, WLOWER, T_DELT_EIGEN,& ! input
              QUADINTENS )                                           ! output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IPARTIC)
                INTENSITY_F(UTA,LUM,IPARTIC,DNIDX) = FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
             ENDDO
          ENDIF

        ENDIF

!  debug

!        if ( uta .eq. n_out_usertaus ) then
!          do i = 1, nstreams
!             write(46,'(i4,f10.5,1p4e21.12)')i,xang(i), &
!               (QUADINTENS(UM,I,IPARTIC,DNIDX),UM=1,N_USER_LEVELS)
!          enddo
!        endif

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE DNUSER_INTENSITY

!

SUBROUTINE MIFLUX_INTENSITY                                         &
           ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTBEAM,     & ! Input
             IPARTIC, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,         & ! Input (remove Thread)
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,               & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,               & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS,                            & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,           & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,                & ! Input
             MEAN_INTENSITY, FLUX_INTEGRAL,                         & ! Output
             DNMEAN_DIRECT,  DNFLUX_DIRECT )                          ! Output

!  Flux and Actinic flux calculations
!  This routine has the thermal solution included.

!    Direct-beam contributions output separately
!       Beta-Coded 26 May 11, Added 24 August 2011 to Version 3.5.1

!  dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAXLAYERS, MAX_PARTLAYERS,  &
                              MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS,        &
                              ZERO, ONE, HALF, PI2, PI4, UPIDX, DNIDX

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Thread
!      INTEGER  , intent(in)  :: THREAD

!  Index

      INTEGER  , intent(in)  :: IPARTIC

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux factor

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  local solar zenith angle cosine

      REAL(fpk), intent(in)  :: LOCAL_CSZA ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Quadrature-defined solutions

      REAL(fpk), intent(in)  :: QUADINTENS &
            (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Output arguments
!  ----------------

!  Mean intensity (actinic flux)

      REAL(fpk), intent(inout) :: MEAN_INTENSITY (MAX_USER_LEVELS,MAXBEAMS,MAX_DIRECTIONS)

!  Flux

      REAL(fpk), intent(inout) :: FLUX_INTEGRAL  (MAX_USER_LEVELS,MAXBEAMS,MAX_DIRECTIONS)

!    Direct-beam contributions output separately, 26 May 11

      REAL(fpk), intent(inout) :: DNMEAN_DIRECT (MAX_USER_LEVELS,MAXBEAMS)
      REAL(fpk), intent(inout) :: DNFLUX_DIRECT (MAX_USER_LEVELS,MAXBEAMS)

!  local variables
!  ---------------

      INTEGER    :: I, UTA, UT, N
      REAL(fpk)  :: SMI, SFX, FMU0
      REAL(fpk)  :: DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI

!  Upwelling
!  ---------

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_USER_LEVELS
          SMI = ZERO
          SFX = ZERO
          DO I = 1, NSTREAMS
            SMI = SMI + QUAD_WEIGHTS(I)*QUADINTENS(UTA,I,IPARTIC,UPIDX)
            SFX = SFX + QUAD_STRMWTS(I)*QUADINTENS(UTA,I,IPARTIC,UPIDX)
          ENDDO
          MEAN_INTENSITY(UTA,IPARTIC,UPIDX) = SMI * HALF
          FLUX_INTEGRAL (UTA,IPARTIC,UPIDX) = SFX * PI2
        ENDDO
      ENDIF

!  Downwelling
!  -----------

      IF ( DO_DNWELLING ) THEN

!  Diffuse contribution

        DO UTA = 1, N_USER_LEVELS
          SMI = ZERO
          SFX = ZERO
          DO I = 1, NSTREAMS
            SMI = SMI + QUAD_WEIGHTS(I)*QUADINTENS(UTA,I,IPARTIC,DNIDX)
            SFX = SFX + QUAD_STRMWTS(I)*QUADINTENS(UTA,I,IPARTIC,DNIDX)
          ENDDO
          MEAN_INTENSITY(UTA,IPARTIC,DNIDX) = SMI * HALF
          FLUX_INTEGRAL (UTA,IPARTIC,DNIDX) = SFX * PI2
        ENDDO

!  nothing to do if no solar source direct beam contribution.

        IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) RETURN

!  add the direct beam contributions

!  loop over all the output optical depths

        DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Only contributions for layers above the PI cutoff
!    Direct-beam contributions output separately, 26 May 11

            IF ( N .LE. LAYER_PIS_CUTOFF(IPARTIC) ) THEN
              DIRECT_TRANS = INITIAL_TRANS(N,IPARTIC) * T_UTDN_MUBAR(UT,IPARTIC)
              DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
              FMU0 = LOCAL_CSZA(N,IPARTIC) * FLUX_FACTOR
              DIRECT_FLUX  = FMU0 * DIRECT_TRANS
              MEAN_INTENSITY(UTA,IPARTIC,DNIDX) = &
                   MEAN_INTENSITY(UTA,IPARTIC,DNIDX) + DIRECT_MEANI
              FLUX_INTEGRAL(UTA,IPARTIC,DNIDX)  = & 
                   FLUX_INTEGRAL(UTA,IPARTIC,DNIDX)  + DIRECT_FLUX
              DNMEAN_DIRECT(UTA,IPARTIC) = DIRECT_MEANI       ! Addition 5/26/11
              DNFLUX_DIRECT(UTA,IPARTIC) = DIRECT_FLUX        ! Addition 5/26/11
            ENDIF

!  For the on-grid values
!    Direct-beam contributions output separately, 26 May 11

          ELSE
            N = UTAU_LEVEL_MASK_DN(UTA)
            IF ( N .LE. LAYER_PIS_CUTOFF(IPARTIC) ) THEN
              IF ( N .EQ. 0 ) THEN
                DIRECT_TRANS = ONE
                FMU0 = LOCAL_CSZA(1,IPARTIC) * FLUX_FACTOR
              ELSE
                DIRECT_TRANS = INITIAL_TRANS(N,IPARTIC)*T_DELT_MUBAR(N,IPARTIC)
                FMU0 = LOCAL_CSZA(N,IPARTIC) * FLUX_FACTOR
              ENDIF
              DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
              DIRECT_FLUX  = FMU0 * DIRECT_TRANS
              MEAN_INTENSITY(UTA,IPARTIC,DNIDX) = &
                   MEAN_INTENSITY(UTA,IPARTIC,DNIDX) + DIRECT_MEANI
              FLUX_INTEGRAL(UTA,IPARTIC,DNIDX)  = &
                   FLUX_INTEGRAL(UTA,IPARTIC,DNIDX)  + DIRECT_FLUX
              DNMEAN_DIRECT(UTA,IPARTIC) = DIRECT_MEANI       ! Addition 5/26/11
              DNFLUX_DIRECT(UTA,IPARTIC) = DIRECT_FLUX        ! Addition 5/26/11
            ENDIF
          ENDIF

!  End loop over optical depth output values

        ENDDO

!  Finish Downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_INTENSITY

!

SUBROUTINE GET_TOASOURCE &
      ( DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, IBEAM, &
        TOA_SOURCE )

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, &
                              ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: IBEAM
      REAL(fpk), intent(out) :: TOA_SOURCE(MAX_USER_STREAMS)

!  local variables

      INTEGER         ::    UM

!  initialise TOA source function

      if ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = 1, N_USER_STREAMS
          TOA_SOURCE(UM) = ZERO
        ENDDO
      else
        TOA_SOURCE(IBEAM) = ZERO
      endif

!  Finish

      RETURN
END SUBROUTINE GET_TOASOURCE

!
!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE

SUBROUTINE GET_BOASOURCE                                                &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_THERMAL,    & ! Input @@@@
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM, & ! Input
            DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS, DO_OBSERVATION_GEOMETRY, & ! Input
            NSTREAMS, NLAYERS, N_USER_STREAMS, IPARTIC, FOURIER,        & ! Input
            SURFACE_FACTOR, QUAD_STRMWTS, QUAD_WEIGHTS, T_DELT_DISORDS, & ! Input
            LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_WLOWER,       & ! Input
            ALBEDO, BRDF_F, USER_BRDF_F, USER_DIRECT_BEAM,              & ! Input
            SURFBB, EMISSIVITY, USER_EMISSIVITY,                        & ! Input
            BOA_SOURCE, DIRECT_BOA_SOURCE,                              & ! Output
            BOA_THTONLY_SOURCE, IDOWNSURF )                               ! Output

!  Bottom of the atmosphere source term
!  This routine has the thermal solution included.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS,   MAX_USER_STREAMS, MAXSTREAMS, &
                              MAXSTREAMS_2, MAXLAYERS,        MAXBEAMS,   &
                              ZERO, ONE

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE
      LOGICAL  , intent(in)  :: DO_MSMODE_THERMAL
!  @@@@@@@@@@@ End Robfix 13 January 2012.

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  surface multiplier, albedo and Fourier/beam indices

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO
      INTEGER  , intent(in)  :: FOURIER
      INTEGER  , intent(in)  :: IPARTIC

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )

!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F &
        ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F &
        ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )

!  Discrete ordinate tranmsittances (thermal transmittance only)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  General thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Direct beam solutions

      REAL(fpk), intent(in)  :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Subroutine output arguments
!  ---------------------------

!  BOA source terms

      REAL(fpk), intent(out) :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk), intent(out) :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  Thermal transmittance-only source

      REAL(fpk), intent(out) :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(fpk), intent(out) :: IDOWNSURF(MAXSTREAMS)

!  local variables
!  ---------------

      LOGICAL      :: DO_QTHTONLY
      INTEGER      :: M, N, J, I, UM, AA, K, LUM
      REAL(fpk)    :: PAR, HOM, REFLEC, KMULT, THELP(MAXSTREAMS)

!  Local indices

      M   = FOURIER
      N   = NLAYERS
      LUM = 1

!  initialise boa source function

      IF ( DO_USER_STREAMS ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = 1, N_USER_STREAMS
            BOA_SOURCE(UM)        = ZERO
            DIRECT_BOA_SOURCE(UM) = ZERO
          ENDDO
        ELSE
          BOA_SOURCE(IPARTIC)        = ZERO
          DIRECT_BOA_SOURCE(IPARTIC) = ZERO
        ENDIF
      ENDIF

!  Special flag, thermal tranmsittance-only, Mean-value output

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .AND.DO_INCLUDE_MVOUTPUT )
      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          BOA_THTONLY_SOURCE(I) = ZERO
        ENDDO
      ENDIF

!  Can return if Fourier > 0 and Lambertian or thermal-only case

      IF ( FOURIER.GT.0.and..not.DO_BRDF_SURFACE ) RETURN

!  First calculate reflectance integrand
!  -------------------------------------

!  Always need this for the Surface weighting functions

!  Thermal transmittance-only solution, build from TOA downwards
!     --> Develop reflectance integrand  a(j).I(-j)

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          THELP(I) = ZERO
          DO K = 1, NLAYERS
            THELP(I) = THELP(I)*T_DELT_DISORDS(I,K) + T_WLOWER(I,K)
          ENDDO
          IDOWNSURF(I) = QUAD_WEIGHTS(I) * THELP(I)
        ENDDO
      ENDIF

!  Full solution with scattering:
!      Downward intensity at computational angles (beam/homog)
!     --> Develop reflectance integrand  a(j).x(j).I(-j)

      IF ( .not.DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          PAR = WLOWER(I,N)
          HOM = ZERO
          DO AA = 1, NSTREAMS
            HOM = HOM + LCON_XVEC(I,AA,N)*T_DELT_EIGEN(AA,N) + &
                        MCON_XVEC(I,AA,N)
          ENDDO
          IDOWNSURF(I) = QUAD_STRMWTS(I) * ( PAR + HOM )
        ENDDO
      ENDIF

!  No reflectance from surface if no surface!
!    ( This includes the Dark Surface case)

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  ###### Lambertian reflectance (same for all user-streams)

      IF ( .not. DO_BRDF_SURFACE ) THEN
        KMULT = SURFACE_FACTOR * ALBEDO
        IF ( FOURIER .EQ. 0 ) THEN
          REFLEC = ZERO
          DO J = 1, NSTREAMS
            REFLEC = REFLEC + IDOWNSURF(J)
          ENDDO
          REFLEC = KMULT * REFLEC
          IF ( DO_USER_STREAMS ) THEN
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = 1, N_USER_STREAMS
                BOA_SOURCE(UM) = REFLEC
              ENDDO
            ELSE
              BOA_SOURCE(IPARTIC) = REFLEC
            ENDIF
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = REFLEC
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!  ###### BRDF reflectance

      IF ( DO_BRDF_SURFACE ) THEN
        IF ( DO_USER_STREAMS ) THEN
          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = 1, N_USER_STREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + IDOWNSURF(J)* USER_BRDF_F(M,UM,J)
              ENDDO
              BOA_SOURCE(UM) = REFLEC * SURFACE_FACTOR
            ENDDO
          ELSE
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + IDOWNSURF(J)* USER_BRDF_F(M,IPARTIC,J)
            ENDDO
            BOA_SOURCE(IPARTIC) = REFLEC * SURFACE_FACTOR
          ENDIF
        ENDIF
        IF ( DO_QTHTONLY ) THEN
          DO I = 1, NSTREAMS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + IDOWNSURF(J) * BRDF_F(M,I,J)
            ENDDO
            BOA_THTONLY_SOURCE(I) = KMULT * REFLEC
          ENDDO
        ENDIF
      ENDIF

!  Add direct beam if flagged

      IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = 1, N_USER_STREAMS
            DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM,IPARTIC)
          ENDDO
        ELSE
          DIRECT_BOA_SOURCE(IPARTIC) = USER_DIRECT_BEAM(LUM,IPARTIC)
        ENDIF
      ENDIF

!  Add surface emission term if flagged

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Use DO_MSMODE_THERMAL flag to control Direct Surface emission

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          IF ( DO_USER_STREAMS.and..not.DO_MSMODE_THERMAL ) THEN ! @@ 1/13/12
            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = 1, N_USER_STREAMS
                BOA_SOURCE(UM) = BOA_SOURCE(UM)+SURFBB*USER_EMISSIVITY(UM)
              ENDDO
            ELSE
              BOA_SOURCE(IPARTIC) = BOA_SOURCE(IPARTIC)+SURFBB*USER_EMISSIVITY(IPARTIC)
            ENDIF
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) + &
                   SURFBB * EMISSIVITY(I)
            ENDDO
          ENDIF
        ELSE
          REFLEC = SURFBB * ( ONE - ALBEDO )
          IF ( DO_USER_STREAMS.and..not.DO_MSMODE_THERMAL ) THEN ! @@ 1/13/12
            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = 1, N_USER_STREAMS
                BOA_SOURCE(UM) = BOA_SOURCE(UM) + REFLEC
              ENDDO
            ELSE
              BOA_SOURCE(IPARTIC) = BOA_SOURCE(IPARTIC) + REFLEC
            ENDIF
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) + REFLEC
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!  debug

!      um = 97
!      if (   do_fdtest ) um = 98
!      if ( fourier_component.eq.0.and.ibeam.eq.1) then
!         write(um,'(1p4e17.9)')BOA_SOURCE(1)+DIRECT_BOA_SOURCE(1)
!      ENDIF

!  Finish

      RETURN
END SUBROUTINE GET_BOASOURCE

!

SUBROUTINE LIDORT_CONVERGE                                             &
      ( DO_UPWELLING, DO_NO_AZIMUTH, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, & ! Input
        DO_SS_EXTERNAL,                                                & ! Input New 15 March 2012
        DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                & ! Input
        DO_DOUBLE_CONVTEST, N_CONVTESTS, LIDORT_ACCURACY,              & ! Input
        N_USER_STREAMS, N_USER_LEVELS, N_USER_RELAZMS,                 & ! Input
        NSTREAMS, IBEAM, FOURIER_COMPONENT,                            & ! Input (Remove thread)
        UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,        & ! Input
        AZMFAC, INTENSITY_F, INTENSITY_SS, INTENSITY_DB,               & ! Input
        INTENSITY, FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )            ! Output

!  convergence testing on the Radiance intensity

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAX_DIRECTIONS, MAX_GEOMETRIES,  &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS,     &
                              ZERO, UPIDX, DNIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  Local flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH
      LOGICAL  , intent(in)  :: DO_RAYLEIGH_ONLY
      LOGICAL  , intent(in)  :: DO_ALL_FOURIER

      LOGICAL  , intent(in)  :: DO_SSFULL
      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING
!  New 15 March 2012
      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL

!  Convergence control

      LOGICAL  , intent(in)  :: DO_DOUBLE_CONVTEST
      INTEGER  , intent(in)  :: N_CONVTESTS
      REAL(fpk), intent(in)  :: LIDORT_ACCURACY

!  Fourier component and beam. Thread removed

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_RELAZMS

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Bookkeeping: Offsets for geometry indexing

      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local number of azimuths and azimuth factors

      INTEGER  , intent(in)  :: LOCAL_N_USERAZM
      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(fpk), intent(in)  :: INTENSITY_F &
          (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter solutions

      REAL(fpk), intent(in)  :: INTENSITY_SS &
          (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Direct-beam results

      REAL(fpk), intent(in)  :: INTENSITY_DB &
          (MAX_USER_LEVELS,MAX_GEOMETRIES)

!  modified/output variables
!  -------------------------

!  Intensity

      REAL(fpk), intent(inout) :: INTENSITY &
            (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Number of saved Fourier components

      INTEGER  , intent(inout) :: FOURIER_SAVED ( MAXBEAMS )

!  Modified output for testing convergence

      LOGICAL  , intent(inout) :: LOCAL_ITERATION
      INTEGER  , intent(inout) :: TESTCONV

!  local variables
!  ---------------

!  local variables

      INTEGER      :: COUNT, COUNT_A
      INTEGER      :: I, IDIR, UT, UA, W, V
      REAL(fpk)    :: TNEW, ACCUR, TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .not. DO_SSFULL ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = UMOFF(IBEAM,I) + UA
                  INTENSITY(UT,V,W) = INTENSITY_F(UT,I,IBEAM,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
         !write(*,*)INTENSITY(5,V,W)
         FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
         DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = UMOFF(IBEAM,I) + UA
                  INTENSITY(UT,V,W) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_USER_LEVELS
             DO I = 1, N_USER_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                 V = UMOFF(IBEAM,I) + UA
                INTENSITY(UT,V,W) = &
                   INTENSITY(UT,V,W) + INTENSITY_SS(UT,V,W)
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDIF

!  Add the Direct bounce to the upwelling

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        IF ( DO_UPWELLING ) THEN
         !IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                INTENSITY(UT,V,UPIDX) = &
                 INTENSITY(UT,V,UPIDX) + INTENSITY_DB(UT,V)
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!  For each azimuth, add Fourier component

          DO UA = 1, LOCAL_N_USERAZM

!     - for direction, user optical depth, out stream

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = UMOFF(IBEAM,I) + UA
                  TOLD = INTENSITY(UT,V,W)
                  TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F(UT,I,IBEAM,W)
                  INTENSITY(UT,V,W) = TOLD + TAZM
                ENDDO
              ENDDO
            ENDDO

          ENDDO

!  Examine convergence on intensity only 
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AND
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO UA = 1, N_USER_RELAZMS
            COUNT_A = 0
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = UMOFF(IBEAM,I) + UA
                  TOLD = INTENSITY(UT,V,W)
                  TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F(UT,I,IBEAM,W)
                  TNEW = TOLD + TAZM
                  IF ( TAZM .NE. ZERO ) THEN
                    ACCUR     = DABS(TAZM/TNEW)
                    IF ( ACCUR .LT. LIDORT_ACCURACY ) THEN
                      COUNT   = COUNT + 1
                      COUNT_A = COUNT_A + 1
                    ENDIF
                  ELSE
                    COUNT   = COUNT + 1
                    COUNT_A = COUNT_A + 1
                  ENDIF
                  INTENSITY(UT,V,W)     = TNEW
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
                LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER_COMPONENT .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
          ENDIF
        ENDIF

!  For all Fourier, keep saving the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_CONVERGE

SUBROUTINE LIDORT_CONVERGE_OBSGEO                                      &
      ( DO_UPWELLING, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER,                & ! Input
        DO_SS_EXTERNAL, DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,& ! Input
        DO_DOUBLE_CONVTEST, N_CONVTESTS, LIDORT_ACCURACY,              & ! Input
        N_USER_LEVELS, NSTREAMS, IBEAM, FOURIER_COMPONENT,             & ! Input (Remove Thread)
        N_DIRECTIONS, WHICH_DIRECTIONS,                                & ! Input
        AZMFAC, INTENSITY_F, INTENSITY_SS, INTENSITY_DB,               & ! Input
        INTENSITY, FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )            ! Output

!  convergence testing on the Radiance intensity

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAX_DIRECTIONS, MAX_GEOMETRIES,  &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS,     &
                              ZERO, UPIDX, DNIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  Local flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_RAYLEIGH_ONLY
      LOGICAL  , intent(in)  :: DO_ALL_FOURIER

      LOGICAL  , intent(in)  :: DO_SSFULL
      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING
!  New 15 March 2012
      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL

!  Convergence control

      LOGICAL  , intent(in)  :: DO_DOUBLE_CONVTEST
      INTEGER  , intent(in)  :: N_CONVTESTS
      REAL(fpk), intent(in)  :: LIDORT_ACCURACY

!  Fourier component and beam. Thread removed

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  azimuth factors

      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(fpk), intent(in)  :: INTENSITY_F &
          (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter solutions

      REAL(fpk), intent(in)  :: INTENSITY_SS &
          (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Direct-beam results

      REAL(fpk), intent(in)  :: INTENSITY_DB &
          (MAX_USER_LEVELS,MAX_GEOMETRIES)

!  modified/output variables
!  -------------------------

!  Intensity

      REAL(fpk), intent(inout) :: INTENSITY &
            (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Number of saved Fourier components

      INTEGER  , intent(inout) :: FOURIER_SAVED ( MAXBEAMS )

!  Modified output for testing convergence

      LOGICAL  , intent(inout) :: LOCAL_ITERATION
      INTEGER  , intent(inout) :: TESTCONV

!  local variables
!  ---------------

!  local variables

      INTEGER      :: COUNT, IDIR, UT, W, LUM, LUA
      REAL(fpk)    :: TNEW, ACCUR, TOLD, TAZM

!  Local user indices

      LUM = 1
      LUA = 1

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .not. DO_SSFULL ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              INTENSITY(UT,IBEAM,W) = INTENSITY_F(UT,LUM,IBEAM,W)
            ENDDO
          ENDDO
        ELSE
         FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
         DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              INTENSITY(UT,IBEAM,W) = ZERO
            ENDDO
          ENDDO
        ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_USER_LEVELS
             INTENSITY(UT,IBEAM,W) = &
                 INTENSITY(UT,IBEAM,W) + INTENSITY_SS(UT,IBEAM,W)
           ENDDO
         ENDDO
       ENDIF

!  Add the Direct bounce to the upwelling

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        IF ( DO_UPWELLING ) THEN
         !IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN
          DO UT = 1, N_USER_LEVELS
             INTENSITY(UT,IBEAM,UPIDX) = &
                 INTENSITY(UT,IBEAM,UPIDX) + INTENSITY_DB(UT,IBEAM)
          ENDDO
         ENDIF
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!     - for direction, user optical depth

           DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                 TOLD = INTENSITY(UT,IBEAM,W)
                 TAZM = AZMFAC(LUM,IBEAM,LUA)*INTENSITY_F(UT,LUM,IBEAM,W)
                 INTENSITY(UT,IBEAM,W) = TOLD + TAZM
              ENDDO
           ENDDO

!  Examine convergence on intensity only 
!  -------------------------------------

!  convergence test applied to ALL directions user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
                TOLD = INTENSITY(UT,IBEAM,W)
                TAZM = AZMFAC(LUM,IBEAM,LUA)*INTENSITY_F(UT,LUM,IBEAM,W)
                TNEW = TOLD + TAZM
                IF ( TAZM .NE. ZERO ) THEN
                   ACCUR = DABS(TAZM/TNEW)
                   IF ( ACCUR .LT. LIDORT_ACCURACY ) THEN
                      COUNT = COUNT + 1
                   ENDIF
                ELSE
                   COUNT = COUNT + 1
                ENDIF
                INTENSITY(UT,IBEAM,W) = TNEW
             ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
              LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER_COMPONENT .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
          ENDIF
        ENDIF

!  For all Fourier, keep saving the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_CONVERGE_OBSGEO

!  End

end module lidort_intensity
