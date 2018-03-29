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

module lidort_lp_wfatmos

!  Parameter types

   USE LIDORT_PARS, only : fpk, TAYLOR_SMALL

!  Dependencies

   USE lidort_lp_PostProcessing

!private
public

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #   -- Master routines for Column Jacobians              #
! #                                                        #
! #            UPUSER_PROFILEWF                            #
! #            DNUSER_PROFILEWF                            #
! #            MIFLUX_PROFILEWF                            #
! #                                                        #
! #   -- Master routines for post-processed TOA/BOA fields #
! #                                                        #
! #            GET_LP_TOASOURCE                            #
! #            GET_LP_BOASOURCE                            #
! #                                                        #
! #   -- Master routines for azimuth Convergence           #
! #                                                        #
! #            LIDORT_LP_CONVERGE                          #
! #            LIDORT_LP_CONVERGE_OBSGEO                   #
! #                                                        #
! ##########################################################

contains

SUBROUTINE UPUSER_PROFILEWF &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES,                        & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,              & ! Input
              DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT,                    & ! Input
              DO_OBSERVATION_GEOMETRY, FLUX_MULTIPLIER,                 & ! Input
              NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,         & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! Input
              UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                  & ! Input
              FOURIER_COMPONENT, IBEAM, VARIATION_INDEX, K_PARAMETERS,  & ! Input
              TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,    & ! Input
              QUAD_STREAMS, USER_SECANTS, DO_LAYER_SCATTERING,          & ! Input
              LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,          & ! Input
              T_DELT_MUBAR,   T_UTDN_MUBAR,                      & ! Input
              T_DELT_USERM,   T_UTUP_USERM,                      & ! Input
              T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN,      & ! Input
              T_DELT_DISORDS, T_DISORDS_UTUP,                    & ! Input
              GAMMA_M, GAMMA_P, SIGMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON,            & ! Input
              T_WUPPER, BOA_THTONLY_SOURCE,                      & ! Input
              U_XPOS, U_XNEG, U_WPOS, HMULT_1, HMULT_2, EMULT_UP,        & ! Input
              UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                    & ! Input
              PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,              & ! Input
              UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_UP,                    & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                & ! Input
              LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                  & ! Input
              L_T_DELT_USERM,   L_T_UTUP_USERM,                   & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTUP,                 & ! Input
              L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Input
              L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,                 & ! Input
              L_XPOS, L_XNEG, L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC,  & ! Input
              L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                         & ! Input
              L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,              & ! Input
              L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP, L_UT_T_PARTIC, L_T_WUPPER, & ! Input
              L_BOA_THTONLY_SOURCE, LP_BOA_MSSOURCE, LP_BOA_DBSOURCE,        & ! Input
              FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,                & ! Output
              PROFILEWF_F, QUADPROFILEWF )                                     ! Output

!  Upwelling post-processed Profile Jacobian Fourier component

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAXBEAMS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS,       &
                              MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Observational geometry control

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Flux multiplier = F/4.pi

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

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

!  Input Fourier number and beam index
!  surface factor (2 for m = 0, 1 otherwise). Not required.

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, IBEAM

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS
      INTEGER  , intent(in)  :: VARIATION_INDEX

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Regular Inputs
!  --------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: GAMMA_M ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_P ( MAXSTREAMS, MAXLAYERS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: BTERM_SAVE ( MAXSTREAMS, MAXLAYERS )

!  Initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature streams, only required for Transonly

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  Solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Cumulative source terms

      REAL(fpk), intent(in)  :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(in)  :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS,MAXBEAMS, MAX_ATMOSWFS  )

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, solar beam

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LP_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LP_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LP_UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized direct thermal solution (upwelling)

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UP    ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTUP  ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized BOA source terms

      REAL(fpk), intent(in)  :: LP_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Stuff for thermal transonly computations (Quadratures)
!  ------------------------------------------------------

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS   ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTUP   ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WUPPER   ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized partial-layer Thermal solution

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Special case thermal transmittance - BOA source term + linearization

      REAL(fpk), intent(in)  ::   BOA_THTONLY_SOURCE (MAXSTREAMS)
      REAL(fpk), intent(in)  :: L_BOA_THTONLY_SOURCE (MAXSTREAMS,MAX_ATMOSWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Profile weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADPROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS,  MAX_USER_LEVELS, &
           MAXSTREAMS,   MAXBEAMS,   MAX_DIRECTIONS )

!  Profile weighting functions at user angles

      REAL(fpk), intent(inout) :: PROFILEWF_F &
         ( MAX_ATMOSWFS,     MAXLAYERS, MAX_USER_LEVELS, &
           MAX_USER_STREAMS, MAXBEAMS,  MAX_DIRECTIONS)

!  Linearized Green's function multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_LP_GMULT(MAX_PARTLAYERS)

!  Linearized Green functions multipliers for off-grid optical depths
!   Will only be calculated as output, if the flag has been set

      REAL(fpk), intent(inout) :: LP_UT_GMULT_UP &
         (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LP_UT_GMULT_DN &
         (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL    :: SOURCETERM_FLAG
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER    :: UTA, UM, Q, NC, UT, IB, K, LUM, M
      INTEGER    :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

      REAL(fpk)  :: L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_FINAL_SOURCE

!  Local post-processing control

      PPSTREAM_MASK(:,IBEAM) = 0
      IF ( DO_OBSERVATION_GEOMETRY ) THEN
         N_PPSTREAMS = 1; PPSTREAM_MASK(1,IBEAM) = IBEAM
      else
         N_PPSTREAMS = N_USER_STREAMS
         do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IBEAM) = UM
         enddo
      endif

!  Indices

      IB = IBEAM
      M  = FOURIER_COMPONENT
      K   = VARIATION_INDEX

!  Zero all Fourier components - New rule, better for safety

!      IF ( DO_USER_STREAMS ) THEN       ! @@@ removed for safety
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, K_PARAMETERS
            DO LUM = 1, N_PPSTREAMS
              PROFILEWF_F(Q,K,UTA,LUM,IB,UPIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
!      ENDIF                             ! @@@ removed for safety

!  Initialize post-processing recursion
!  ====================================

!  start the recursion

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
               L_CUMUL_SOURCE(UM,Q) = LP_BOA_MSSOURCE(UM,Q) + LP_BOA_DBSOURCE(UM,Q)
            ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

            SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            CALL LP_WHOLELAYER_STERM_UP &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, & ! input
              DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,            & ! input, FLags + Order
              IB, N, K, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,  & ! input, Numbers
              DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM,        & ! input, Optical + User_streams
              LAYER_PIS_CUTOFF,  INITIAL_TRANS,    T_DELT_MUBAR,             & ! input, Beam stuff
              LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,          & ! input, Beam stuff
              U_XPOS, U_XNEG, U_WPOS, LCON, MCON,                            & ! input, RTE Solutions
              L_KEIGEN, L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,           & ! input, Linearized RTE solutios
              GAMMA_M, GAMMA_P, SIGMA_P, L_LAYER_TSUP_UP,                    & ! input, solutions
              ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,            & ! input, Greens Function Saved values
              HMULT_1, HMULT_2, EMULT_UP, PMULT_UU, PMULT_UD,                & ! input, Multipliers
              L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                             & ! input, Linearized HMULT/EMULT
              L_LAYER_SOURCE )                                                 ! output

!  Cumulative sourceterm

            IF ( N.EQ.K ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, K_PARAMETERS
                     L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)          + &
                           T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,Q)   + &
                         L_T_DELT_USERM(N,UM,Q) *   CUMSOURCE_UP(UM,NC-1)
                  ENDDO
               ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, K_PARAMETERS
                     L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  + &
                          T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
                  ENDDO
               ENDDO
            ENDIF

!  End layer recursion

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

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            FLAGS_LP_GMULT(UT) = .TRUE.
            CALL QUADPROFILEWF_OFFGRID_UP  &
           ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             TAYLOR_ORDER, NSTREAMS, NLAYERS, IB, UTA, UT, N, K, K_PARAMETERS,              & ! Input
             QUAD_STREAMS, FLUX_MULTIPLIER, PARTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,   & ! Input
             INITIAL_TRANS, AVERAGE_SECANT,  T_DELT_MUBAR, T_UTDN_MUBAR,                    & ! Input
             T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_DISORDS, T_DISORDS_UTUP,                    & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, UT_GMULT_UP, UT_GMULT_DN,            & ! Input
             XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,          & ! Input
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,         & ! Input
             L_KEIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_T_DELT_DISORDS, L_T_DISORDS_UTUP,  & ! Input
             L_ATERM_SAVE, L_BTERM_SAVE, L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,              & ! Input
             L_UT_T_PARTIC, L_T_WUPPER, L_BOA_THTONLY_SOURCE,                               & ! Input
             QUADPROFILEWF, LP_UT_GMULT_UP, LP_UT_GMULT_DN )                                  ! Output
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

            CALL LP_PARTLAYER_STERM_UP  &
           ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,       & ! input, Flags
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                  & ! input, FLags + Order
             IB, UT, N, K, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,    & ! input, Numbers
             PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTUP_USERM, & ! input, Optical
             LAYER_PIS_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,         & ! input, Beam stuff
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,                & ! input, Beam stuff
             U_XPOS, U_XNEG, U_WPOS, LCON, MCON,                                  & ! input, RTE Solutions
             L_KEIGEN, L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,                 & ! input, Linearized RTE solutios
             GAMMA_M, GAMMA_P, SIGMA_P, L_LAYER_TSUP_UTUP,                        & ! input, solutions
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,                  & ! input, Greens Function values
             UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, UT_PMULT_UU, UT_PMULT_UD,     & ! input
             L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,                        & ! input
             L_LAYER_SOURCE )                                                       ! output

!  Cumulative and final

            IF ( N.EQ.K ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)             + &
                          T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q) + &
                        L_T_UTUP_USERM(UT,UM,Q)*   CUMSOURCE_UP(UM,NC)
                    PROFILEWF_F(Q,K,UTA,LUM,IB,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
               ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q) + T_UTUP_USERM(UT,UM) * L_CUMUL_SOURCE(UM,Q)
                    PROFILEWF_F(Q,K,UTA,LUM,IB,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADPROFILEWF_LEVEL_UP                              &
            ( DO_THERMAL_TRANSONLY, FLUX_MULTIPLIER,                 & ! Input
              NLAYERS, NSTREAMS, IB, UTA, NLEVEL, K, K_PARAMETERS,   & ! Input
              QUAD_STREAMS, L_XPOS, L_XNEG, L_WLOWER, L_WUPPER,      & ! Input
              LCON, LCON_XVEC, NCON_XVEC,   T_DELT_EIGEN,            & ! Input
              MCON, MCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN,            & ! Input
              T_DELT_DISORDS, L_T_DELT_DISORDS, T_WUPPER, L_T_WUPPER,& ! Input
              BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,              & ! Input
              QUADPROFILEWF )                                          ! Output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, K_PARAMETERS
                   PROFILEWF_F(Q,K,UTA,LUM,IB,UPIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
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
END SUBROUTINE UPUSER_PROFILEWF

!

SUBROUTINE DNUSER_PROFILEWF &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES,                         & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,               & ! Input
              DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT,                     & ! Input
              DO_OBSERVATION_GEOMETRY, FLUX_MULTIPLIER,                  & ! Input
              NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                   & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
              UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                   & ! Input
              FOURIER_COMPONENT, IBEAM, VARIATION_INDEX, K_PARAMETERS,   & ! Input
              TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,     & ! Input
              QUAD_STREAMS, USER_SECANTS, DO_LAYER_SCATTERING,           & ! Input
              LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,           & ! Input
              T_DELT_MUBAR,   T_UTDN_MUBAR,                      & ! Input
              T_DELT_USERM,   T_UTDN_USERM,                      & ! Input
              T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN,      & ! Input
              T_DELT_DISORDS, T_DISORDS_UTDN,                    & ! Input
              GAMMA_M, GAMMA_P, SIGMA_M, ATERM_SAVE, BTERM_SAVE, & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WLOWER,          & ! Input
              U_XPOS, U_XNEG, U_WNEG, HMULT_1, HMULT_2, EMULT_DN,        & ! Input
              UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN,                    & ! Input
              PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,              & ! Input
              UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_DN,                    & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                & ! Input
              LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                  & ! Input
              L_T_DELT_USERM,   L_T_UTDN_USERM,                   & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTDN,                 & ! Input
              L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Input
              L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,                 & ! Input
              L_XPOS, L_XNEG, L_WLOWER, NCON_XVEC, PCON_XVEC,            & ! Input
              L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                         & ! Input
              L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,              & ! Input
              L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN, L_UT_T_PARTIC,         & ! Input
              L_T_WLOWER, LP_TOA_SOURCE,                                 & ! Input
              FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,            & ! Output
              PROFILEWF_F, QUADPROFILEWF )                                 ! Output

!  Downwelling post-processed Profile Jacobian Fourier component

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAXBEAMS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS,       &
                              MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Observational geometry control

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Flux multiplier = F/4.pi

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Input Fourier number and beam index

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, IBEAM

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS
      INTEGER  , intent(in)  :: VARIATION_INDEX

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Regular Inputs
!  --------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: GAMMA_M ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_P ( MAXSTREAMS, MAXLAYERS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: BTERM_SAVE ( MAXSTREAMS, MAXLAYERS )

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature streams, only required for Transonly

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Cumulative source terms

      REAL(fpk), intent(in)  :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(in)  :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized transmittances, solar beam

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LP_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LP_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LP_UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Thermal solutions for the Trans-only special case

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_DN &
         (MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTDN &
         (MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized TOA source terms

      REAL(fpk), intent(in)  :: LP_TOA_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Stuff for thermal transonly computations (Quadratures)
!  ------------------------------------------------------

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTDN (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN &
         (MAXSTREAMS, MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: T_WLOWER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Particular solutions

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC &
         (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Profile weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADPROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS,  MAX_USER_LEVELS, &
           MAXSTREAMS,   MAXBEAMS,   MAX_DIRECTIONS )

!  Profile weighting functions at user angles

      REAL(fpk), intent(inout) :: PROFILEWF_F &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
           MAX_USER_STREAMS, MAXBEAMS,  MAX_DIRECTIONS)

!  Linearized Green functions multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_LP_GMULT(MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL    :: LOCAL_LP_GMULT
      LOGICAL    :: SOURCETERM_FLAG
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER    :: UTA, UM, Q, NC, UT, IB, K, LUM, M
      INTEGER    :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

      REAL(fpk)  :: L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_FINAL_SOURCE

!  Local post-processing control

      PPSTREAM_MASK(:,IBEAM) = 0
      IF ( DO_OBSERVATION_GEOMETRY ) THEN
         N_PPSTREAMS = 1; PPSTREAM_MASK(1,IBEAM) = IBEAM
      else
         N_PPSTREAMS = N_USER_STREAMS
         do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IBEAM) = UM
         enddo
      endif

!  Initialise local indices

      IB = IBEAM
      M  = FOURIER_COMPONENT
      K   = VARIATION_INDEX

!  Zero all Fourier component output

!      IF ( DO_USER_STREAMS ) THEN       @@@ Removed for safety
      DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
            DO LUM = 1, N_PPSTREAMS
               PROFILEWF_F(Q,K,UTA,LUM,IB,DNIDX) = ZERO
            ENDDO
         ENDDO
      ENDDO
!      ENDIF                             @@@ Removed for safety

!  Initialize post-processing recursion
!  ====================================

! Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            DO Q = 1, K_PARAMETERS
              L_CUMUL_SOURCE(UM,Q) = LP_TOA_SOURCE(UM,Q)
            ENDDO
         ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
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

!  Layer source term

            SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            CALL LP_WHOLELAYER_STERM_DN  &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, & ! input
              DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,            & ! input, FLags + Order
              IB, N, K, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,  & ! input, Numbers
              DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM,        & ! input, Optical + User_streams
              LAYER_PIS_CUTOFF, AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR, & ! input, Beam stuff
              LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,          & ! input, Beam stuff
              U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                            & ! input, RTE Solutions
              L_KEIGEN, L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,           & ! input, Linearized RTE solutios
              GAMMA_M, GAMMA_P, SIGMA_M, L_LAYER_TSUP_DN,                    & ! input, solutions
              ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,            & ! input, Greens Function Saved values
              HMULT_1, HMULT_2, EMULT_DN, PMULT_DU, PMULT_DD,                & ! input, Multipliers
              L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                             & ! input, Linearized HMULT/EMULT
              L_LAYER_SOURCE )                                                 ! output

!  Cumulative source

            IF ( N.EQ.K ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO Q = 1, K_PARAMETERS
                    L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  + &
                         T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q) + &
                        L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
                  ENDDO
                ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO Q = 1, K_PARAMETERS
                    L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + &
                         T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
                  ENDDO
                ENDDO
            ENDIF

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

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            LOCAL_LP_GMULT = FLAGS_LP_GMULT(UT)
            CALL QUADPROFILEWF_OFFGRID_DN                                   &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
              TAYLOR_ORDER, NSTREAMS, IB, UTA, UT, N, K, K_PARAMETERS,                       & ! Input
              QUAD_STREAMS, FLUX_MULTIPLIER, PARTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,   & ! Input
              LOCAL_LP_GMULT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,     & ! Input
              T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_DISORDS, T_DISORDS_UTDN,                    & ! Input
              GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, UT_GMULT_UP, UT_GMULT_DN,            & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WLOWER,                              & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,         & ! Input
              L_KEIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_T_DELT_DISORDS, L_T_DISORDS_UTDN,  & ! Input
              L_ATERM_SAVE, L_BTERM_SAVE, L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,              & ! Input
              L_UT_T_PARTIC, L_T_WLOWER,                                                     & ! Input
              QUADPROFILEWF, LP_UT_GMULT_UP, LP_UT_GMULT_DN )                                  ! Output
            FLAGS_LP_GMULT(UT) = .FALSE.
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LP_PARTLAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! input, Flags
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,               & ! input, FLags + Order
             IB, UT, N, K, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK, & ! input, Numbers
             PARTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTDN_USERM,           & ! input, Optical
             LAYER_PIS_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! input, Beam stuff
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,             & ! input, Beam stuff
             U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                               & ! input, RTE Solutions
             L_KEIGEN, L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,              & ! input, Linearized RTE solutios
             GAMMA_M, GAMMA_P, SIGMA_M, L_LAYER_TSUP_UTDN,                     & ! input, solutions
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,               & ! input, Greens Function Saved values
             UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, UT_PMULT_DU, UT_PMULT_DD,  & ! Input, Multipliers
             L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,                     & ! Input, Multipliers
             L_LAYER_SOURCE )                                                    ! output

!  Cumulative and final

            IF ( N.EQ.K ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)            + &
                        T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,Q) + &
                      L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,NC)
                    PROFILEWF_F(Q,K,UTA,LUM,IB,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
               ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q) + &
                        T_UTDN_USERM(UT,UM) * L_CUMUL_SOURCE(UM,Q)
                    PROFILEWF_F(Q,K,UTA,LUM,IB,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
               ENDDO
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADPROFILEWF_LEVEL_DN                  &
            ( DO_THERMAL_TRANSONLY, FLUX_MULTIPLIER,     & ! Input
              NSTREAMS, IB, UTA, NLEVEL, K, K_PARAMETERS,& ! Input
              QUAD_STREAMS, L_XPOS, L_XNEG, L_WLOWER,    & ! Input
              LCON, LCON_XVEC, NCON_XVEC, T_DELT_EIGEN,  & ! Input
              MCON, PCON_XVEC, L_T_DELT_EIGEN,           & ! Input
              T_DELT_DISORDS, L_T_DELT_DISORDS, T_WLOWER,& ! Input
              L_T_WLOWER,                                & ! Input
              QUADPROFILEWF )                              ! Output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IBEAM)
                DO Q = 1, K_PARAMETERS
                   PROFILEWF_F(Q,K,UTA,LUM,IB,DNIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
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
END SUBROUTINE DNUSER_PROFILEWF

!

SUBROUTINE MIFLUX_PROFILEWF &
           ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTBEAM,         & ! Input
             IB, K, K_PARAMETERS, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR, & ! Input (Remove thread)
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                   & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS,                                & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,               & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADPROFILEWF,                 & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,        & ! Input
             MINT_PROFILEWF, FLUX_PROFILEWF )                            ! Output

!  Profile Jacobians for the hemispherically integrated fields

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                              MAXSTREAMS,   MAXBEAMS,  MAX_PARTLAYERS,  &
                              MAX_DIRECTIONS, ZERO, HALF, PI2, PI4, UPIDX, DNIDX

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

      INTEGER  , intent(in)  :: IB

!  linearization control

      INTEGER  , intent(in)  :: K, K_PARAMETERS

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

!  Linearized initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized transmittances

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Quadrature-defined weighting functions

      REAL(fpk), intent(in)  :: QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS )

!  Output arguments
!  ----------------

!mick fix 6/29/11 - changed these two from "out" to "inout"

!  Mean intensity (actinic flux)

      REAL(fpk), intent(inout) :: MINT_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                   MAXBEAMS,     MAX_DIRECTIONS )

!  Flux

      REAL(fpk), intent(inout) :: FLUX_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                   MAXBEAMS,     MAX_DIRECTIONS )

!  local variables
!  ----------------

      INTEGER    :: I, UTA, UT, Q, N
      REAL(fpk)  :: SM, SF, FMU0
      REAL(fpk)  :: L_TRANS, L_DIRECT_FLUX, L_DIRECT_MEANI

!  mean intensity and flux
!  -----------------------

!  Upwelling

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
          SM = ZERO
          SF = ZERO
          DO I = 1, NSTREAMS
           SM = SM + QUAD_WEIGHTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX)
           SF = SF + QUAD_STRMWTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX)
          ENDDO
          MINT_PROFILEWF(Q,K,UTA,IB,UPIDX) = SM * HALF
          FLUX_PROFILEWF(Q,K,UTA,IB,UPIDX) = SF * PI2
         ENDDO
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN

!  Diffuse term contribution

        DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
          SM = ZERO
          SF = ZERO
          DO I = 1, NSTREAMS
           SM = SM + QUAD_WEIGHTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX)
           SF = SF + QUAD_STRMWTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX)
          ENDDO
          MINT_PROFILEWF(Q,K,UTA,IB,DNIDX) = SM * HALF
          FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX) = SF * PI2
         ENDDO
        ENDDO

!  nothing to do if no solar sources

        IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) RETURN

!  For the downward direction, add the direct beam contributions

        DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Only contributions for layers above the PI cutoff
!    L_INITIAL_TRANS is a logarithmic derivative

            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              FMU0 = LOCAL_CSZA(N,IB) * FLUX_FACTOR
              DO Q = 1, K_PARAMETERS
                L_TRANS = LP_T_UTDN_MUBAR(UT,K,IB,Q) + LP_INITIAL_TRANS(N,K,IB,Q) * T_UTDN_MUBAR(UT,IB)
                L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
                L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                L_DIRECT_FLUX  = FMU0 * L_TRANS
                MINT_PROFILEWF(Q,K,UTA,IB,DNIDX) = MINT_PROFILEWF(Q,K,UTA,IB,DNIDX) + L_DIRECT_MEANI
                FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX) = FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX) + L_DIRECT_FLUX
              ENDDO
            ENDIF

!  For the on-grid balues

          ELSE
            N = UTAU_LEVEL_MASK_DN(UTA)
            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              IF ( N.GT.0 ) THEN
                FMU0 = LOCAL_CSZA(N,IB) * FLUX_FACTOR
                DO Q = 1, K_PARAMETERS
                  L_TRANS = LP_T_DELT_MUBAR(N,K,IB,Q) + LP_INITIAL_TRANS(N,K,IB,Q) * T_DELT_MUBAR(N,IB)
                  L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
                  L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                  L_DIRECT_FLUX  = FMU0 * L_TRANS
                  MINT_PROFILEWF(Q,K,UTA,IB,DNIDX) = MINT_PROFILEWF(Q,K,UTA,IB,DNIDX) + L_DIRECT_MEANI
                  FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX) = FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX) + L_DIRECT_FLUX
                ENDDO
              ENDIF
            ENDIF
          ENDIF

!  Close loops

        ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_PROFILEWF

!

SUBROUTINE GET_LP_TOASOURCE &
      ( DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, IBEAM, K_PARAMETERS, &
        LP_TOA_SOURCE )

!  Linearized Top of the atmosphere source term

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Observational geometry control

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Control integer

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Beam

      INTEGER  , intent(in)  :: IBEAM

!  Linearization

      INTEGER  , intent(in)  :: K_PARAMETERS

!  Output

      REAL(fpk), intent(out) :: LP_TOA_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       :: UM, Q

!  initialise TOA source function
!  ------------------------------

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            LP_TOA_SOURCE(UM,Q) = ZERO
          ENDDO
        ENDDO
      ELSE
        DO Q = 1, K_PARAMETERS
          LP_TOA_SOURCE(IBEAM,Q) = ZERO
        ENDDO
      ENDIF

!  Finish

END SUBROUTINE GET_LP_TOASOURCE

!

SUBROUTINE GET_LP_BOASOURCE                                           &
      ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,& ! Input
        DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,                         & ! Input
        DO_OBSERVATION_GEOMETRY,                                      & ! Input
        NLAYERS, NSTREAMS, N_USER_STREAMS, FOURIER_COMPONENT,         & ! Input
        IBEAM, K, K_PARAMETERS, QUAD_STRMWTS, QUAD_WEIGHTS,           & ! Input
        SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                  & ! Input
        USER_DIRECT_BEAM, DELTAU_SLANT, L_DELTAU_VERT,                & ! Input
        LCON, MCON, LCON_XVEC, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input
        T_WLOWER, L_XPOS, L_XNEG, L_T_WLOWER, L_WLOWER,               & ! Input
        NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_T_DELT_DISORDS,       & ! Input
        L_BOA_MSSOURCE, L_BOA_DBSOURCE, L_BOA_THTONLY_SOURCE )          ! Output

!  Linearized Bottom of the atmosphere source term

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_2, &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  control integers

      INTEGER  , intent(in)  :: NLAYERS, NSTREAMS, N_USER_STREAMS

!  surface multiplier, albedo and Fourier/beam indices

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO
      INTEGER  , intent(in)  :: FOURIER_COMPONENT
      INTEGER  , intent(in)  :: IBEAM

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )

!    incident quadrature streams, reflected quadrature streams

      REAL(FPK), intent(in)  :: BRDF_F &
         ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!    incident quadrature streams, reflected user streams

      REAL(FPK), intent(in)  :: USER_BRDF_F &
         ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )

!  linearization control

      INTEGER  , intent(in)  :: K, K_PARAMETERS

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Direct beam solutions

      REAL(fpk), intent(in)  :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Linearized deltaQUAD_STRMWTSu

      REAL(fpk), intent(in)  :: L_DELTAU_VERT(MAX_ATMOSWFS, MAXLAYERS )

!  Slant optical thickness values

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Discrete ordinate transmittances

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: T_WLOWER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER &
         (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: L_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL    :: DO_QTHTONLY
      INTEGER    :: M, N, J, I, UM, AA, Q, IB, NN, LUM
      REAL(fpk)  :: DOWN (MAXSTREAMS)
      REAL(fpk)  :: L_DOWN (MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: REFLEC, L_BEAM, FAC, KMULT, SPAR
      REAL(fpk)  :: SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

!  Starting section
!  ----------------

!  Local indices

      M   = FOURIER_COMPONENT
      N   = NLAYERS
      IB  = IBEAM
      LUM = 1

!  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .and. DO_INCLUDE_MVOUTPUT )

!  initialise linearized BOA source functions

      IF ( DO_USER_STREAMS ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_BOA_MSSOURCE(UM,Q) = ZERO
              L_BOA_DBSOURCE(UM,Q) = ZERO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, K_PARAMETERS
            L_BOA_MSSOURCE(IB,Q) = ZERO
            L_BOA_DBSOURCE(IB,Q) = ZERO
          ENDDO
        ENDIF
      ENDIF

!  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            L_BOA_THTONLY_SOURCE(I,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Exit if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  1. Thermal Transmittance only
!     %%%%%%%%%%%%%%%%%%%%%%%%%%

!  Linearization of downwelling quadrature field at surface
!   ---Thermal transmittance solution, build from TOA downwards

      IF ( DO_THERMAL_TRANSONLY ) THEN

!  Initialise

        DO I = 1, NSTREAMS
          DOWN(I) = ZERO
          DO Q = 1, K_PARAMETERS
            L_DOWN(I,Q) = ZERO
          ENDDO
        ENDDO

!  Build

        DO NN = 1, NLAYERS
          IF ( K.EQ.NN ) THEN
            DO I = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                L_DOWN(I,Q) = L_DOWN(I,Q) *   T_DELT_DISORDS(I,NN) &
                              + DOWN(I)   * L_T_DELT_DISORDS(I,NN,Q) &
                              + L_T_WLOWER(I,NN,Q)
              ENDDO
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                L_DOWN(I,Q) = L_DOWN(I,Q) * T_DELT_DISORDS(I,NN) 
              ENDDO
            ENDDO
          ENDIF
          DO I = 1, NSTREAMS
            DOWN(I) = DOWN(I)*T_DELT_DISORDS(I,NN) + T_WLOWER(I,NN)
          ENDDO
        ENDDO

!  2. Scattering solutions
!     %%%%%%%%%%%%%%%%%%%%

!  Linearization of downwelling quadrature field at surface
!    Set reflectance integrand  a(j).x(j).L_DOWN(-j)

      ELSE

!  Two cases:
!  If  K = N, this is also the layer that is varying --> Extras!
!  If  N > K with variations in layer K above N
!    Add Particular integral linearization

        IF ( K.EQ.N ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
                HOM4 = PCON_XVEC(I,AA,N,Q)
                HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = ZERO
              IF ( DO_SOLAR_SOURCES .OR. DO_INCLUDE_THERMEMISS ) THEN
                SPAR = L_WLOWER(I,N,Q)
              ENDIF
              L_DOWN(I,Q) = SPAR + SHOM
            ENDDO
          ENDDO
        ELSE IF (K.LT.N.AND.K.NE.0) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = PCON_XVEC(I,AA,N,Q)
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              SPAR = ZERO
              IF ( DO_SOLAR_SOURCES ) THEN
                SPAR = L_WLOWER(I,N,Q)
              ENDIF
              L_DOWN(I,Q) = SPAR + SHOM
             ENDDO
          ENDDO
        ENDIF

      ENDIF

!    Set reflectance integrand  a(j).x(j).L_DOWN(-j)  Scattering solutions
!    Set reflectance integrand       a(j).L_DOWN(-j)  Thermal tranmsittance

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS
            L_DOWN(I,Q) = L_DOWN(I,Q) * QUAD_WEIGHTS(I)
          ENDDO
        ENDDO
      ELSE
        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS
          L_DOWN(I,Q) = L_DOWN(I,Q) * QUAD_STRMWTS(I)
          ENDDO
        ENDDO
      ENDIF

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  .. integrate reflectance, same for all user-streams in Lambertian case

      IF ( .not. DO_BRDF_SURFACE ) THEN

        KMULT = SURFACE_FACTOR * ALBEDO
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO Q = 1, K_PARAMETERS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + L_DOWN(J,Q)
            ENDDO
            REFLEC = KMULT * REFLEC
            IF ( DO_USER_STREAMS ) THEN
              IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
                DO UM = 1, N_USER_STREAMS
                  L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
                ENDDO
              ELSE
                L_BOA_MSSOURCE(IB,Q) = L_BOA_MSSOURCE(IB,Q) + REFLEC
              ENDIF
            ENDIF
            IF ( DO_QTHTONLY ) THEN
              DO I = 1, NSTREAMS
               L_BOA_THTONLY_SOURCE(I,Q) = &
                 L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  .. integrate reflectance, BRDF case

      ELSE IF ( DO_BRDF_SURFACE ) THEN

        DO Q = 1, K_PARAMETERS
          IF ( DO_USER_STREAMS ) THEN
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = 1, N_USER_STREAMS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  REFLEC = REFLEC + L_DOWN(J,Q) * USER_BRDF_F(M,UM,J)
                ENDDO
                REFLEC = SURFACE_FACTOR * REFLEC
                L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
              ENDDO
            ELSE
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + L_DOWN(J,Q) * USER_BRDF_F(M,IB,J)
              ENDDO
              REFLEC = SURFACE_FACTOR * REFLEC
              L_BOA_MSSOURCE(IB,Q) = L_BOA_MSSOURCE(IB,Q) + REFLEC
            ENDIF
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + L_DOWN(J,Q) * BRDF_F(M,I,J)
              ENDDO
              REFLEC = SURFACE_FACTOR* REFLEC
              L_BOA_THTONLY_SOURCE(I,Q) = &
                   L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
            ENDDO
          ENDIF
        ENDDO

      ENDIF

!  Add direct beam if flagged.
!   This is different from the column WF case

      IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = 1, N_USER_STREAMS
            FAC = - USER_DIRECT_BEAM(UM,IB) * DELTAU_SLANT(N,K,IB)
            DO Q = 1, K_PARAMETERS
              L_BEAM = L_DELTAU_VERT(Q,K) * FAC
              L_BOA_DBSOURCE(UM,Q) = L_BEAM
            ENDDO
          ENDDO
        ELSE
          FAC = - USER_DIRECT_BEAM(LUM,IB) * DELTAU_SLANT(N,K,IB)
          DO Q = 1, K_PARAMETERS
            L_BEAM = L_DELTAU_VERT(Q,K) * FAC
            L_BOA_DBSOURCE(IB,Q) = L_BEAM
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
END SUBROUTINE GET_LP_BOASOURCE

!

SUBROUTINE LIDORT_LP_CONVERGE                         &
           ( DO_UPWELLING, DO_SS_EXTERNAL, DO_SSFULL, & ! Input
             DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,     & ! Input
             DO_NO_AZIMUTH, AZMFAC, LOCAL_N_USERAZM,  & ! Input
             N_USER_STREAMS, N_USER_LEVELS, NLAYERS,  & ! Input
             IBEAM, FOURIER_COMPONENT,                & ! Input (remove thread)
             UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS,   & ! Input
             LAYER_VARY_FLAG, LAYER_VARY_NUMBER,      & ! Input
             PROFILEWF_F, PROFILEWF_SS, PROFILEWF_DB, & ! Input
             PROFILEWF )                                ! Output

!  Just upgrades the weighting function Fourier cosine series

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, MAX_GEOMETRIES,   &
                              MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  Local flags

      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH
      LOGICAL  , intent(in)  :: DO_UPWELLING
!  New 15 March 2012
      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL
      LOGICAL  , intent(in)  :: DO_SSFULL
      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING

!  FOurier component and beam

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Bookkeeping: Offsets for geometry indexing

      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local number of azimuths and azimuth factors

      INTEGER  , intent(in)  :: LOCAL_N_USERAZM
      REAL(fpk), intent(in)  :: AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Linearization control

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Fourier-component Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                              MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_SS ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Direct-bounce Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES )

!  output
!  ------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: PROFILEWF ( MAX_ATMOSWFS,   MAXLAYERS,  MAX_USER_LEVELS, &
                                              MAX_GEOMETRIES, MAX_DIRECTIONS )

!  local variables

      INTEGER       :: I, IDIR, UT, UA, Q, W, V, N

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Full single scatter calculation is initialized to zero here (Version 2.3)

!  Bulk/column atmospheri! weighting functions (Version 3.3)
!  ---------------------------------------------------------

!  This section newly written, 26 September 2006, installed 14 May 2007.

!  Diffuse field at all output angles

        IF ( .NOT. DO_SSFULL ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    PROFILEWF(Q,N,UT,V,W) = PROFILEWF_F(Q,N,UT,I,IBEAM,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ELSE
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    PROFILEWF(Q,N,UT,V,W) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
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
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    PROFILEWF(Q,N,UT,V,W) = &
                      PROFILEWF(Q,N,UT,V,W) + PROFILEWF_SS(Q,N,UT,V,W)
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  Add the Direct bounce to the upwelling

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        IF ( DO_UPWELLING ) THEN
         !IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN
          DO N = 1, NLAYERS
           IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                   PROFILEWF(Q,N,UT,V,UPIDX) = &
                     PROFILEWF(Q,N,UT,V,UPIDX) + PROFILEWF_DB(Q,N,UT,V)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF
        ENDIF

!  If no_azimuth, then exit

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Add next Fourier component to output

        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO UA = 1, LOCAL_N_USERAZM
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
               V = UMOFF(IBEAM,I) + UA
               PROFILEWF(Q,N,UT,V,W) = PROFILEWF(Q,N,UT,V,W) + &
                     PROFILEWF_F(Q,N,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LP_CONVERGE

!

SUBROUTINE LIDORT_LP_CONVERGE_OBSGEO                  &
           ( DO_UPWELLING, DO_SS_EXTERNAL, DO_SSFULL, & ! Input
             DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,     & ! Input
             DO_NO_AZIMUTH, AZMFAC,                   & ! Input (argument list corrected)
             N_USER_LEVELS, NLAYERS, IBEAM, FOURIER,  & ! Input (Remove thread)
             N_DIRECTIONS, WHICH_DIRECTIONS,          & ! Input (argument list corrected)
             LAYER_VARY_FLAG, LAYER_VARY_NUMBER,      & ! Input
             PROFILEWF_F, PROFILEWF_SS, PROFILEWF_DB, & ! Input
             PROFILEWF )                                ! Output

!  Just upgrades the weighting function Fourier cosine series

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, MAX_GEOMETRIES,   &
                              MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  Local flags

      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH
      LOGICAL  , intent(in)  :: DO_UPWELLING
!  New 15 March 2012
      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL
      LOGICAL  , intent(in)  :: DO_SSFULL
      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING

!  FOurier component and beam

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Local number of azimuths and azimuth factors

      REAL(fpk), intent(in)  :: AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Linearization control

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Fourier-component Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                              MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_SS ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Direct-bounce Profile weighting functions at user angles

      REAL(fpk), intent(in)  :: PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES )

!  output
!  ------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: PROFILEWF ( MAX_ATMOSWFS,   MAXLAYERS,  MAX_USER_LEVELS, &
                                              MAX_GEOMETRIES, MAX_DIRECTIONS )

!  local variables

      INTEGER       :: IDIR, UT, Q, W, N, LUM, LUA

!  Local user indices

      LUM = 1
      LUA = 1

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Full single scatter calculation is initialized to zero here (Version 2.3)

!  Bulk/column atmospheri! weighting functions (Version 3.3)
!  ---------------------------------------------------------

!  This section newly written, 26 September 2006, installed 14 May 2007.

!  Diffuse field at all output angles

        IF ( .NOT. DO_SSFULL ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                PROFILEWF(Q,N,UT,IBEAM,W) = PROFILEWF_F(Q,N,UT,LUM,IBEAM,W)
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ELSE
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                PROFILEWF(Q,N,UT,IBEAM,W) = ZERO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
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
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                PROFILEWF(Q,N,UT,IBEAM,W) = &
                  PROFILEWF(Q,N,UT,IBEAM,W) + PROFILEWF_SS(Q,N,UT,IBEAM,W)
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  Add the Direct bounce to the upwelling

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        IF ( DO_UPWELLING ) THEN
         !IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN
          DO N = 1, NLAYERS
           IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              DO UT = 1, N_USER_LEVELS
                PROFILEWF(Q,N,UT,IBEAM,UPIDX) = &
                  PROFILEWF(Q,N,UT,IBEAM,UPIDX) + PROFILEWF_DB(Q,N,UT,IBEAM)
              ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF
        ENDIF

!  If no_azimuth, then exit

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Add next Fourier component to output

        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              PROFILEWF(Q,N,UT,IBEAM,W) = PROFILEWF(Q,N,UT,IBEAM,W) + &
                    PROFILEWF_F(Q,N,UT,LUM,IBEAM,W)*AZMFAC(LUM,IBEAM,LUA)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LP_CONVERGE_OBSGEO

!  End

end module lidort_lp_wfatmos
