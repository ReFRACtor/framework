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

MODULE LIDORT_L_INPUTS

!  Parameter types

   USE LIDORT_pars
   USE lidort_aux, only : GFINDPAR, FINDPAR_ERROR, LEN_STRING, &
                          GAULEG, HPSORT, INDEXX

!private
public

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_L_INPUT_MASTER (master), calling:         #
! #              LIDORT_L_INIT_INPUTS                           #
! #                                                             #
! #    These routines are called by the Main LIDORT modules     #
! #                                                             #
! #            LIDORT_L_CHECK_INPUT_DIMS                        #
! #                                                             #
! ###############################################################

!  Version 3.7, Internal threading removed

contains

SUBROUTINE LIDORT_L_INPUT_MASTER ( &
        FILNAM,            & ! INPUT
        LIDORT_FixIn,      & ! OUTPUTS
        LIDORT_ModIn,      & ! OUTPUTS
        LIDORT_LinFixIn,   & ! OUTPUTS
        LIDORT_LinModIn,   & ! OUTPUTS
        LIDORT_InputStatus ) ! OUTPUTS

!  Parameter types
!      USE LIDORT_pars, only: fpk, MAXLAYERS, MAXBEAMS,  MAX_MESSAGES, &
!                             MAX_USER_RELAZMS, MAX_USER_STREAMS, &
!                             MAX_USER_LEVELS, MAX_USER_OBSGEOMS, &
!                             LIDORT_SUCCESS, LIDORT_SERIOUS, LIDORT_INUNIT

      USE LIDORT_Inputs_def
      USE LIDORT_LinInputs_def
      USE LIDORT_Outputs_def
      USE LIDORT_Inputs

!  Implicit none

      implicit none

!  Input argument (input filename)
!  -------------------------------

      CHARACTER (LEN=*), intent(in) :: FILNAM

!  Outputs

      TYPE(LIDORT_Fixed_Inputs), INTENT (OUT)             :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT (OUT)          :: LIDORT_ModIn
      TYPE(LIDORT_Fixed_LinInputs), INTENT (OUT)          :: LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs), INTENT (OUT)       :: LIDORT_LinModIn
      TYPE(LIDORT_Input_Exception_Handling), INTENT (OUT) :: LIDORT_InputStatus

!  LIDORT local variables
!  ++++++++++++++++++++++

!  1. CONTROL FLAGS
!  ================

!  Basic top-level control
!   -- Thermal to be reintroduced later 2010

      LOGICAL ::    DO_SOLAR_SOURCES
      LOGICAL ::    DO_THERMAL_EMISSION

!  directional control

      LOGICAL ::    DO_UPWELLING
      LOGICAL ::    DO_DNWELLING

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL ::    DO_FULLRAD_MODE

!  stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL ::    DO_USER_STREAMS

!  single scatter corrections. (Includes direct beam reflection)
!    - NADIR    : Plane-parallel line-of-sight
!    - OUTGOING : Line-of-sight in curved atmosphere

      LOGICAL ::    DO_SSCORR_NADIR         ! May be re-set after Checking
      LOGICAL ::    DO_SSCORR_OUTGOING      ! May be re-set after Checking

!  New 15 March 2012
      LOGICAL ::    DO_SS_EXTERNAL

!  Flag for Full-up single scatter calculation. (No MS field)
!    One of the above two SSCORR flags must be set

      LOGICAL ::    DO_SSFULL

!  Flag for performing a complete separate delta-M truncation on the
!  single scatter corrrection  calculations. **** Use with CAUTION.

      LOGICAL ::    DO_SSCORR_TRUNCATION

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL ::    DO_PLANE_PARALLEL

!  Flag for use of BRDF surface
!    - If not set, default to Lambertian surface

      LOGICAL ::    DO_BRDF_SURFACE

!  Surface emission flag

      LOGICAL ::    DO_SURFACE_EMISSION

!  Transmittance only for thermal mode.

      LOGICAL ::    DO_THERMAL_TRANSONLY

!  mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL ::    DO_ADDITIONAL_MVOUT     ! May be re-set after Checking

!  mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL ::    DO_MVOUT_ONLY           ! May be re-set after Checking

!  particular solution control
!    Removed in stripped down version, March 2008
!      LOGICAL ::    DO_CLASSICAL_SOLUTION

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set. 

      LOGICAL ::    DO_CHAPMAN_FUNCTION     ! May be re-set after Checking

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set. 

      LOGICAL ::    DO_REFRACTIVE_GEOMETRY  ! May be re-set after Checking

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL ::    DO_DELTAM_SCALING       ! May be re-set after Checking

!  double convergence test flag

      LOGICAL ::    DO_DOUBLE_CONVTEST      ! May be re-set after Checking

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL ::    DO_SOLUTION_SAVING      ! May be re-set after Checking
      LOGICAL ::    DO_BVP_TELESCOPING      ! May be re-set after Checking

!  scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!
!    - Isotropic only, if set, phase function is 1

      LOGICAL ::    DO_RAYLEIGH_ONLY        ! May be re-set after Checking
      LOGICAL ::    DO_ISOTROPIC_ONLY       ! May be re-set after Checking

!  Debug and testing flags
!   - Normally should not be set

      LOGICAL ::    DO_NO_AZIMUTH           ! May be re-set after Checking
      LOGICAL ::    DO_ALL_FOURIER

!  New 17 May 2012
      LOGICAL ::    DO_SURFACE_LEAVING
      LOGICAL ::    DO_SL_ISOTROPIC

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      LOGICAL ::    DO_OBSERVATION_GEOMETRY

!  2. CONTROL INTEGERS
!  ===================

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7
      
      INTEGER :: TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  number of computational layers

      INTEGER :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER :: NFINELAYERS

!  number of Legendre phase function expansion moments

      INTEGER :: NMOMENTS_INPUT          ! May be re-set after Checking

!  number of solar beams to be processed

      INTEGER :: NBEAMS

!  Number of user-defined relative azimuths

      INTEGER :: N_USER_RELAZMS

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER :: N_USER_STREAMS

!  Number of User-defined vertical levels for  output

      INTEGER :: N_USER_LEVELS

!  Number of thermal coefficients (2 should be the default)

      INTEGER :: N_THERMAL_COEFFS

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      INTEGER :: N_USER_OBSGEOMS

!  3. CONTROL NUMBERS
!  ==================

!  Flux factor ( should be 1 or pi ). Same for all beams.

      REAL(fpk) :: FLUX_FACTOR

!  accuracy for convergence of Fourier series

      REAL(fpk) :: LIDORT_ACCURACY

!  Zenith tolerance (nearness of output zenith cosine to 1.0 )
!    removed 02 June 2010
!      REAL(fpk) :: ZENITH_TOLERANCE

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(fpk) :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(fpk) :: RFINDEX_PARAMETER

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      REAL(fpk) :: GEOMETRY_SPECHEIGHT

!  BOA solar zenith angles (degrees)

      REAL(fpk) :: BEAM_SZAS ( MAXBEAMS )

!  user-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      REAL(fpk) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  User-defined viewing zenith angles input (degrees) 

      REAL(fpk) :: USER_ANGLES ( MAX_USER_STREAMS )

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      REAL(fpk) :: USER_LEVELS ( MAX_USER_LEVELS )

!  User-defined Observation Geometry angle input
!   New variable, 25 OCtober 2012, for Observational Geometry input

      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: USER_OBSGEOMS

!  Linearized inputs

      LOGICAL ::    DO_COLUMN_LINEARIZATION
      LOGICAL ::    DO_PROFILE_LINEARIZATION
      LOGICAL ::    DO_SURFACE_LINEARIZATION
      LOGICAL ::    DO_SLEAVE_WFS

!  BlackBody Jacobian Flags, Introduced March 18th 2014

      LOGICAL  ::   DO_ATMOS_LBBF
      LOGICAL  ::   DO_SURFACE_LBBF

!  Derived inputs

      LOGICAL ::    LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER ::    LAYER_VARY_NUMBER( MAXLAYERS )

      INTEGER ::    N_TOTALCOLUMN_WFS
      INTEGER ::    N_SURFACE_WFS
      INTEGER ::    N_SLEAVE_WFS

!  Exception handling
!  ------------------

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER       :: STATUS
      INTEGER       :: NMESSAGES
      CHARACTER*120 :: MESSAGES( 0:MAX_MESSAGES )
      CHARACTER*120 :: ACTIONS ( 0:MAX_MESSAGES )

!  ============================================================
!  ============================================================

!  SPECIAL NOTE on variable GEOMETRY_SPECHEIGHT

!    This is only required when the outgoing sphericity correction is
!    in operation. Otherwise, the regular pseudo-spherical correction
!    (wiht or without an exact single-scatter correction) uses the same
!    set of angles all the up the nadir from the bottom of atmosphere.

!     This height is normally set equal to the height at the lowest
!     level of the atmosphere: GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
!     In this case, no adjustment to the geometrical inputs is needed
!     for the outgoing sphericity correction.

!     If there is a situation GEOMETRY_SPECHEIGHT < HEIGHT_GRID(NLAYERS),
!     then an adjustment to the geometrical inputs is needed for the
!     outgoing sphericity correction. This adjustment is internal and
!     the model output will still be given at the geometrical angles
!     as specified by the user, even though these angles may not be the
!     ones at which the calculations were done. This situation will occur
!     when we are given a BOA geometry but we want to make a calculation
!     for a reflecting surface (such as a cloud-top) which is above the
!     BOA level. In this case, GEOMETRY_SPECHEIGHT = 0.0, and the lowest
!     height HEIGHT_GRID(NLAYERS) = cloud-top.

!     This height cannot be greater than HEIGHT_GRID(NLAYERS). If this is
!     the case, this height will be set equal to HEIGHT_GRID(NLAYERS), and
!     the calculation will go through without the adjustment. A warning
!     about this incorrect input choice will be sent to LOGFILE.

!  ============================================================
!  ============================================================

!  local variables
!  ---------------

      INTEGER    :: STATUS_SUB, FILUNIT
!      INTEGER    :: LEN_STRING
!      EXTERNAL      LEN_STRING

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  initialize variables
!    New: Observation Geometry variables, 25 October 2012

      CALL LIDORT_INIT_INPUTS &
          ( DO_SOLAR_SOURCES,        DO_THERMAL_EMISSION,              & ! output
            DO_UPWELLING,            DO_DNWELLING,                     & ! output
            DO_FULLRAD_MODE,         DO_USER_STREAMS,                  & ! output
            DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,               & ! output
            DO_SS_EXTERNAL,          DO_OBSERVATION_GEOMETRY,          & ! output (line modified 10/25/12)
            DO_SSFULL,               DO_SSCORR_TRUNCATION,             & ! output
            DO_PLANE_PARALLEL,       DO_THERMAL_TRANSONLY,             & ! output
            DO_BRDF_SURFACE,         DO_SURFACE_EMISSION,              & ! output
            DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                    & ! output
            DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,           & ! output
            DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,               & ! output
            DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,                & ! output
            DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,               & ! output
            DO_ALL_FOURIER,          DO_NO_AZIMUTH,                    & ! output
            DO_SURFACE_LEAVING,      DO_SL_ISOTROPIC,                  & ! output
            TAYLOR_ORDER, NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,           & ! output (modified 10/10/13)
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_OBSGEOMS, N_USER_LEVELS, & ! output (modified 10/25/12)
            N_THERMAL_COEFFS, FLUX_FACTOR, LIDORT_ACCURACY,                         & ! input/output
            EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,                   & ! input/output
            BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_OBSGEOMS, USER_LEVELS )        ! output (modified 10/25/12)

!  initialize linearization variables

      CALL LIDORT_L_INIT_INPUTS &
          ( DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION, DO_ATMOS_LBBF, & ! output
            DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS, DO_SURFACE_LBBF,         & ! output
            N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,            & ! output
            N_SURFACE_WFS, N_SLEAVE_WFS )                                       ! output

!  Open file

      FILUNIT = LIDORT_INUNIT
      OPEN(LIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  read standard inputs

      CALL LIDORT_READ_INPUTS &
          ( DO_SOLAR_SOURCES,        DO_THERMAL_EMISSION,          & ! input/output
            DO_UPWELLING,            DO_DNWELLING,                 & ! input/output
            DO_FULLRAD_MODE,         DO_USER_STREAMS,              & ! input/output
            DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,           & ! input/output
            DO_SS_EXTERNAL,          DO_OBSERVATION_GEOMETRY,      & ! input/output (modified 10/25/12)
            DO_SSFULL,               DO_SSCORR_TRUNCATION,         & ! input/output
            DO_PLANE_PARALLEL,       DO_THERMAL_TRANSONLY,         & ! input/output
            DO_BRDF_SURFACE,         DO_SURFACE_EMISSION,          & ! input/output
            DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                & ! input/output
            DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,       & ! input/output
            DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,           & ! input/output
            DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,            & ! input/output
            DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,           & ! input/output
            DO_ALL_FOURIER,          DO_NO_AZIMUTH,                & ! input/output
            DO_SURFACE_LEAVING,      DO_SL_ISOTROPIC,              & ! input/output
            TAYLOR_ORDER, NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,           & ! input/output (modified 10/10/13)
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_OBSGEOMS, N_USER_LEVELS, & ! output (modified 10/25/12)
            N_THERMAL_COEFFS, FLUX_FACTOR, LIDORT_ACCURACY,                         & ! input/output
            EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,                   & ! input/output
            BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_OBSGEOMS, USER_LEVELS,       & ! output (modified 10/25/12)
            STATUS_SUB, NMESSAGES, MESSAGES, ACTIONS )                                ! Output

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        LIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
        LIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
        LIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
        LIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS
        CLOSE(FILUNIT)
        RETURN
      ENDIF

!  normal execution: Copy all variables and return
!  -----------------------------------------------

!  First close file

      CLOSE(FILUNIT)

!  THE COPYING TASK
!  ================

!  Fixed Boolean inputs

      LIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE      = DO_FULLRAD_MODE
      LIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = DO_SSCORR_TRUNCATION

!  New 15 march 2012
      LIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL       = DO_SS_EXTERNAL

      LIDORT_FixIn%Bool%TS_DO_SSFULL            = DO_SSFULL

      LIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION  = DO_THERMAL_EMISSION
      LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION  = DO_SURFACE_EMISSION

      LIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL    = DO_PLANE_PARALLEL
      LIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE      = DO_BRDF_SURFACE

      LIDORT_FixIn%Bool%TS_DO_UPWELLING         = DO_UPWELLING
      LIDORT_FixIn%Bool%TS_DO_DNWELLING         = DO_DNWELLING

!  New 17 May 2012
      LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING   = DO_SURFACE_LEAVING
      LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC      = DO_SL_ISOTROPIC

!  Fixed control inputs

!  New, 10/10/13
      LIDORT_FixIn%Cont%TS_TAYLOR_ORDER     = TAYLOR_ORDER

      LIDORT_FixIn%Cont%TS_NSTREAMS         = NSTREAMS
      LIDORT_FixIn%Cont%TS_NLAYERS          = NLAYERS
      LIDORT_FixIn%Cont%TS_NFINELAYERS      = NFINELAYERS
      LIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = N_THERMAL_COEFFS
      LIDORT_FixIn%Cont%TS_LIDORT_ACCURACY  = LIDORT_ACCURACY

!  Fixed Beam inputs

      LIDORT_FixIn%SunRays%TS_FLUX_FACTOR = FLUX_FACTOR

!  Fixed User value inputs

      LIDORT_FixIn%UserVal%TS_N_USER_LEVELS = N_USER_LEVELS

!  Moved to the Modified structure, 10/25/12
!      LIDORT_FixIn%UserVal%TS_N_USER_STREAMS = N_USER_STREAMS

!  Fixed Chapman function inputs

      LIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = RFINDEX_PARAMETER

!  Modified Boolean inputs

      LIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR        = DO_SSCORR_NADIR
      LIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING     = DO_SSCORR_OUTGOING

      LIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION

      LIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      LIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY      = DO_ISOTROPIC_ONLY
      LIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH          = DO_NO_AZIMUTH
      LIDORT_ModIn%MBool%TS_DO_ALL_FOURIER         = DO_ALL_FOURIER

      LIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING

      LIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      LIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING

      LIDORT_ModIn%MBool%TS_DO_USER_STREAMS        = DO_USER_STREAMS

      LIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = DO_MVOUT_ONLY

      LIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = DO_THERMAL_TRANSONLY

!  New 25 October 2012
      LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = DO_OBSERVATION_GEOMETRY

!  Modified control inputs

      LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT = NMOMENTS_INPUT

!  Modified Beam inputs

      LIDORT_ModIn%MSunRays%TS_NBEAMS      = NBEAMS
      LIDORT_ModIn%MSunRays%TS_BEAM_SZAS   = BEAM_SZAS

!  Modified User value inputs
!   N_USER_STREAMS is moved from the Fixed input. 10/25/12
!   Observational Goemetry inputs introduced 10/25/12

      LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS      = N_USER_RELAZMS
      LIDORT_ModIn%MUserVal%TS_USER_RELAZMS        = USER_RELAZMS
      LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS      = N_USER_STREAMS
      LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT   = USER_ANGLES
      LIDORT_ModIn%MUserVal%TS_USER_LEVELS         = USER_LEVELS
      LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

      LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS     = N_USER_OBSGEOMS
      LIDORT_ModIn%MUserVal%TS_USER_OBSGEOM_INPUT  = USER_OBSGEOMS

!  Modified Chapman function inputs

      LIDORT_ModIn%MChapman%TS_EARTH_RADIUS = EARTH_RADIUS

!  Fixed Linearized Control inputs

      LIDORT_LinFixIn%Cont%TS_DO_COLUMN_LINEARIZATION  = &
                              DO_COLUMN_LINEARIZATION
      LIDORT_LinFixIn%Cont%TS_DO_PROFILE_LINEARIZATION = &
                              DO_PROFILE_LINEARIZATION
      LIDORT_LinFixIn%Cont%TS_DO_SURFACE_LINEARIZATION = &
                              DO_SURFACE_LINEARIZATION
      LIDORT_LinFixIn%Cont%TS_DO_SLEAVE_WFS     = &
                              DO_SLEAVE_WFS

      LIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS = &
                              N_TOTALCOLUMN_WFS

      LIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG   = &
                              LAYER_VARY_FLAG
      LIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER = &
                              LAYER_VARY_NUMBER

      LIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS     = &
                              N_SURFACE_WFS
      LIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS      = &
                              N_SLEAVE_WFS

!  BlackBody Jacobian Flags, Introduced March 18th 2014

      LIDORT_LinFixIn%Cont%TS_DO_ATMOS_LBBF     = DO_ATMOS_LBBF
      LIDORT_LinFixIn%Cont%TS_DO_SURFACE_LBBF   = DO_SURFACE_LBBF

!  Modified Linearized Control inputs (none at present - 1 dummy variable)

      LIDORT_LinModIn%Dummy = 0

!  Exception handling

      LIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
      LIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
      LIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
      LIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS

!  Return

      RETURN

!  Open file error
!  ===============

!  continuation point

300   CONTINUE

!  Final status

      STATUS = LIDORT_SERIOUS

!  Final message and action

      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for ' // &
                             FILNAM(1:LEN_STRING(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right File!!'

!  Close file

      CLOSE(FILUNIT)

!  Only copy the exception handling

      LIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
      LIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
      LIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
      LIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS

!  Return

      RETURN

!  Finish

END SUBROUTINE LIDORT_L_INPUT_MASTER

SUBROUTINE LIDORT_L_INIT_INPUTS &
          ( DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION, DO_ATMOS_LBBF, & ! output
            DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS, DO_SURFACE_LBBF,         & ! output
            N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, & ! output
            N_SURFACE_WFS, N_SLEAVE_WFS )                            ! output

!  Initializes linearized control inputs for LIDORT
!  ------------------------------------------------

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXLAYERS

      IMPLICIT NONE

!  Output Arguments
!  ================

      LOGICAL, intent(out) ::   DO_COLUMN_LINEARIZATION
      LOGICAL, intent(out) ::   DO_PROFILE_LINEARIZATION
      LOGICAL, intent(out) ::   DO_SURFACE_LINEARIZATION
      LOGICAL, intent(out) ::   DO_SLEAVE_WFS

!  BlackBody Jacobian Flags, Introduced March 18th 2014

      LOGICAL, intent(out) ::   DO_ATMOS_LBBF
      LOGICAL, intent(out) ::   DO_SURFACE_LBBF

      INTEGER, intent(out) ::   N_TOTALCOLUMN_WFS
      LOGICAL, intent(out) ::   LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER, intent(out) ::   LAYER_VARY_NUMBER(MAXLAYERS)

      INTEGER, intent(out) ::   N_SURFACE_WFS
      INTEGER, intent(out) ::   N_SLEAVE_WFS

!  Initialize linearized control variables

      DO_COLUMN_LINEARIZATION  = .FALSE.
      DO_PROFILE_LINEARIZATION = .FALSE.
      DO_SURFACE_LINEARIZATION = .FALSE.
      DO_SLEAVE_WFS            = .FALSE.
      DO_ATMOS_LBBF            = .FALSE.
      DO_SURFACE_LBBF          = .FALSE.

      N_TOTALCOLUMN_WFS = 0

      LAYER_VARY_FLAG   = .FALSE.
      LAYER_VARY_NUMBER = 0

      N_SURFACE_WFS  = 0
      N_SLEAVE_WFS   = 0

! finish

      RETURN
END SUBROUTINE LIDORT_L_INIT_INPUTS

!

SUBROUTINE LIDORT_L_CHECK_INPUT_DIMS &
      ( LIDORT_LinFixIn, LIDORT_LinModIn, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_ATMOSWFS, MAX_SURFACEWFS, &
                              MAX_SLEAVEWFS, MAX_MESSAGES, &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

      USE LIDORT_LinInputs_def

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  LIDORT input structures

      TYPE(LIDORT_Fixed_LinInputs), INTENT (IN) :: &
        LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs), INTENT (IN) :: &
        LIDORT_LinModIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS
      NM = NMESSAGES

!  Check LIDORT linearized input dimensions
!    against maximum dimensions
!  =========================================

      IF ( LIDORT_LinFixIn%Cont%TS_DO_COLUMN_LINEARIZATION ) THEN
        IF ( LIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for column WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_TOTALCOLUMN_WFS or increase MAX_ATMOSWFS in LIDORT.PARS'
         STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

!      IF ( LIDORT_LinFixIn%Cont%TS_DO_PROFILE_LINEARIZATION ) THEN
!        IF ( LIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
!         NM = NM + 1
!         MESSAGES(NM) = &
!             'Bad error: Insuffient dimensioning for profile WFs'
!         ACTIONS(NM)  = &
!             'Action: Decrease N_TOTALPROFILE_WFS or increase MAX_ATMOSWFS in LIDORT.PARS'
!         STATUS = LIDORT_SERIOUS
!        ENDIF
!      ENDIF

      IF ( LIDORT_LinFixIn%Cont%TS_DO_SURFACE_LINEARIZATION ) THEN
        IF ( LIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACE_WFS in LIDORT.PARS'
         STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( LIDORT_LinFixIn%Cont%TS_DO_SLEAVE_WFS ) THEN
        IF ( LIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS in LIDORT.PARS'
         STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

END SUBROUTINE LIDORT_L_CHECK_INPUT_DIMS

END MODULE LIDORT_L_INPUTS

