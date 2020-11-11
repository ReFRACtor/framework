! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #        --           -            -        -        -    #
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_LPS_MASTER (top-level master)             #
! #            LIDORT_LPS_FOURIER_MASTER                        #
! #                                                             #
! ###############################################################

!  Notes for Version 3.5
!  ---------------------

!  Construction August 26-31, 2010
!  R. Spurr, RT SOLUTIONS Inc., 9 Channing Street, Cambridge MA 02138

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/03/14  - Validated against 3p6

!  Additional notes, Version 3.7
!  -----------------------------

!  Internal threading removed, 02 May 2014

MODULE LIDORT_LPS_masters

!  Module for calling LIDORT, radiances and Profile Jacobians

!  Parameter types

   USE LIDORT_pars

!  I/O Structures for LIDORT
!  =========================

   USE LIDORT_IO_DEFS
   USE LIDORT_LIN_IO_DEFS

!  LIDORT module dependencies

   USE lidort_aux, only            : LEN_STRING

   USE lidort_inputs,         only : LIDORT_CHECK_INPUT_DIMS, LIDORT_CHECK_INPUT, &
                                     LIDORT_CHECK_INPUT_OPTICAL, LIDORT_DERIVE_INPUT
   USE lidort_l_inputs,       only : LIDORT_L_CHECK_INPUT_DIMS

   USE lidort_geometry,       only : OBSGEOM_OUTGOING_ADJUSTGEOM, MULTI_OUTGOING_ADJUSTGEOM, &
                                     LOSONLY_OUTGOING_ADJUSTGEOM, LIDORT_CHAPMAN

   USE lidort_corrections   , only : LIDORT_SSCORR_NADIR, LIDORT_DBCORRECTION

   USE lidort_lp_corrections, only : LIDORT_LP_SSCORR_NADIR, LIDORT_LP_SSCORR_OUTGOING, &
                                     LIDORT_LAP_DBCORRECTION
   USE lidort_ls_corrections, only : LIDORT_LS_DBCORRECTION

   USE lidort_miscsetups           ! Need all these

   USE lidort_la_miscsetups        ! Need all these
   USE lidort_lp_miscsetups        ! Need 2 public routines

   USE lidort_solutions            ! Need all these
   USE lidort_lpc_solutions        ! Need all these

   USE lidort_bvproblem,      only : BVP_MATRIXSETUP_MASTER,    BVP_SOLUTION_MASTER, &
                                     BVPTEL_MATRIXSETUP_MASTER, BVPTEL_SOLUTION_MASTER

   USE lidort_lp_bvproblem,   only : LP_BVP_SOLUTION_MASTER, LP_BVPTEL_SOLUTION_MASTER

   USE lidort_intensity            ! Need all these

   USE lidort_lp_wfatmos           ! Need all these

   USE lidort_ls_wfsurface,   only : SURFACEWF_MASTER, LIDORT_LS_CONVERGE, LIDORT_LS_CONVERGE_OBSGEO

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
   USE lidort_ls_wfsleave
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   USE lidort_thermalsup           ! Need all these

   USE lidort_l_thermalsup         ! Need all these

   USE lidort_lbbf_jacobians_m     ! Need this, New 18 March 2014

   USE lidort_writemodules         ! Need all these

   USE lidort_l_writemodules       ! Need all these

IMPLICIT NONE

!  everything public, except the Fourier Master
!PRIVATE :: LIDORT_LPS_FOURIER_MASTER

PUBLIC

CONTAINS

SUBROUTINE LIDORT_LPS_master ( &
        LIDORT_FixIn,    & ! INPUTS
        LIDORT_ModIn,    & ! INPUTS (possibly modified)
        LIDORT_Sup,      & ! INPUTS/OUTPUTS
        LIDORT_Out,      & ! OUTPUTS
        LIDORT_LinFixIn, & ! INPUTS
        LIDORT_LinModIn, & ! INPUTS (possibly modified)
        LIDORT_LinSup,   & ! INPUTS/OUTPUTS
        LIDORT_LinOut )    ! OUTPUTS

!  Parameter types

      USE LIDORT_pars

!  Implicit none

      implicit none

!  LIDORT input structures

      TYPE(LIDORT_Fixed_Inputs)   , INTENT (IN)    :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT (INOUT) :: LIDORT_ModIn

!  LIDORT supplements structure

      TYPE(LIDORT_Sup_InOut), INTENT (INOUT)       :: LIDORT_Sup

!  LIDORT output structure

      !TYPE(LIDORT_Outputs), INTENT (INOUT)         :: LIDORT_Out
      TYPE(LIDORT_Outputs), INTENT (OUT)           :: LIDORT_Out

!  LIDORT linearized input structures

      TYPE(LIDORT_Fixed_LinInputs)   , INTENT (IN)    :: LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs), INTENT (INOUT) :: LIDORT_LinModIn

!  LIDORT linearized supplements structure

      TYPE(LIDORT_LinSup_InOut), INTENT (INOUT)       :: LIDORT_LinSup

!  LIDORT linearized output structure

      !TYPE(LIDORT_LinOutputs), INTENT (INOUT)         :: LIDORT_LinOut
      TYPE(LIDORT_LinOutputs), INTENT (OUT)         :: LIDORT_LinOut

!  LIDORT local variables
!  ++++++++++++++++++++++

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

!  Boolean
!  =======

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL ::    DO_FULLRAD_MODE

!  Flag for external single scatter calculation. New 15 March 2012

      LOGICAL ::    DO_SS_EXTERNAL

!  Flag for Full-up single scatter calculation. (No MS field)
!    One of the above two SSCORR flags must be set

      LOGICAL ::    DO_SSFULL

!  Basic top-level solar beam control

      LOGICAL ::    DO_SOLAR_SOURCES

!  Surface and thermal emission flags

      LOGICAL ::    DO_THERMAL_EMISSION
      LOGICAL ::    DO_SURFACE_EMISSION

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL ::    DO_PLANE_PARALLEL

!  Flag for use of BRDF surface
!    - If not set, default to Lambertian surface

      LOGICAL ::    DO_BRDF_SURFACE

!  Directional control

      LOGICAL ::    DO_UPWELLING
      LOGICAL ::    DO_DNWELLING

!  Surface leaving flags - New 17 May 2012

      LOGICAL ::    DO_SURFACE_LEAVING
      LOGICAL ::    DO_SL_ISOTROPIC

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      LOGICAL ::    DO_OBSERVATION_GEOMETRY

!  Control
!  =======

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7
      
      INTEGER :: TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  Number of computational layers

      INTEGER :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER :: NFINELAYERS

!  Number of thermal coefficients (2 should be the default)

      INTEGER :: N_THERMAL_COEFFS

!  Accuracy for convergence of Fourier series

      REAL(fpk) :: LIDORT_ACCURACY

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      INTEGER :: N_USER_OBSGEOMS

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 15 May 2015, Version 3.7, following 2stream 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(fpk) :: TCUTOFF

!  Sunrays
!  =======

!  Flux factor ( should be 1 or pi ). Same for all beams.

      REAL(fpk) :: FLUX_FACTOR

!  Number of solar beams to be processed

      INTEGER :: NBEAMS

!  BOA solar zenith angles (degrees)

      REAL(fpk) :: BEAM_SZAS ( MAXBEAMS )

!  UserValues
!  ==========

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER :: N_USER_STREAMS

!  User-defined viewing zenith angles input (degrees) 

      REAL(fpk) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Number of user-defined relative azimuths

      INTEGER :: N_USER_RELAZMS

!  User-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      REAL(fpk) :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  Number of User-defined vertical levels for  output

      INTEGER :: N_USER_LEVELS

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      REAL(fpk) :: USER_LEVELS   (MAX_USER_LEVELS)

!  User-defined Observation Geometry angle input
!   New variable, 25 OCtober 2012, for Observational Geometry input

      REAL(fpk) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3)

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      REAL(fpk) :: GEOMETRY_SPECHEIGHT

!  Chapman
!  =======

!  Multilayer Height inputs
!   (only required for the Chapman function calculations)

      REAL(fpk) :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Multilayer atmospheric inputs (P andT)
!   (only required for the Chapman function calculations, refractive geometry)

      REAL(fpk) :: PRESSURE_GRID    ( 0:MAXLAYERS )
      REAL(fpk) :: TEMPERATURE_GRID ( 0:MAXLAYERS )

!  Number of fine layer gradations 
!    (only for Chapman function calculations with refractive geometry)

      INTEGER   :: FINEGRID ( MAXLAYERS )

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(fpk) :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(fpk) :: RFINDEX_PARAMETER

!  Optical
!  =======

!  Multilayer optical property (bulk) inputs

      REAL(fpk) :: DELTAU_VERT_INPUT  ( MAXLAYERS )
      REAL(fpk) :: OMEGA_TOTAL_INPUT  ( MAXLAYERS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk) :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )

!  Thermal Black Body functions

      REAL(fpk) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Lambertian Surface control

      REAL(fpk) :: LAMBERTIAN_ALBEDO

!  Surface Black body inputs

      REAL(fpk) :: SURFACE_BB_INPUT

!  Write
!  =====

!      LOGICAL ::    DO_WRITE_INPUT
!      LOGICAL ::    DO_WRITE_SCENARIO
!      LOGICAL ::    DO_WRITE_FOURIER
!      LOGICAL ::    DO_WRITE_RESULTS

!  --------------------------
!  Standard Inputs - Variable
!  --------------------------

!  Boolean
!  =======

!  Single scatter corrections. (Includes direct beam reflection)
!    - NADIR    : Plane-parallel line-of-sight
!    - OUTGOING : Line-of-sight in curved atmosphere

      LOGICAL ::    DO_SSCORR_NADIR         ! May be re-set after Checking
      LOGICAL ::    DO_SSCORR_OUTGOING      ! May be re-set after Checking

!  Flag for performing a complete separate delta-M truncation on the
!  Single scatter corrrection  calculations. **** Use with CAUTION.

      LOGICAL ::    DO_SSCORR_TRUNCATION         ! May be re-set after Checking

!  Double convergence test flag

      LOGICAL ::    DO_DOUBLE_CONVTEST      ! May be re-set after Checking

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set.

      LOGICAL ::    DO_REFRACTIVE_GEOMETRY  ! May be re-set after Checking

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set.

      LOGICAL ::    DO_CHAPMAN_FUNCTION     ! May be re-set after Checking

!  Scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!
!    - Isotropic only, if set, phase function is 1

      LOGICAL ::    DO_RAYLEIGH_ONLY        ! May be re-set after Checking
      LOGICAL ::    DO_ISOTROPIC_ONLY       ! May be re-set after Checking

!  Debug and testing flags
!   - Normally should not be set

      LOGICAL ::    DO_NO_AZIMUTH           ! May be re-set after Checking
      LOGICAL ::    DO_ALL_FOURIER

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL ::    DO_DELTAM_SCALING       ! May be re-set after Checking

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL ::    DO_SOLUTION_SAVING      ! May be re-set after Checking
      LOGICAL ::    DO_BVP_TELESCOPING      ! May be re-set after Checking

!  Stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL ::    DO_USER_STREAMS

!  Mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL ::    DO_ADDITIONAL_MVOUT     ! May be re-set after Checking

!  Mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL ::    DO_MVOUT_ONLY           ! May be re-set after Checking

!  Flag for use of BRDF surface
!  Transmittance only for thermal mode.

      LOGICAL ::    DO_THERMAL_TRANSONLY

!  Control
!  =======

!  Number of Legendre phase function expansion moments

      INTEGER ::    NMOMENTS_INPUT

!  Chapman
!  =======

!  Chapman factors (from pseudo-spherical geometry)

      REAL(fpk) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  -----------------------
!  Standard Supplement I/O
!  -----------------------

!  BRDF Inputs
!  ===========

!  Exact (direct bounce) BRDF

      REAL(fpk) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of BRDF, in the following order

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: BRDF_F_0      ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: BRDF_F        ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

      REAL(fpk) :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk) :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Emissivity

      REAL(fpk) :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk) :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  SS and DB I/O
!  =============

!  Single scatter intensity

      REAL(fpk) :: INTENSITY_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Direct bounce intensity

      REAL(fpk) :: INTENSITY_DB (MAX_USER_LEVELS,MAX_GEOMETRIES)

!  Surface-Leaving Inputs, 17 May 12
!  =================================

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Exact Surface-Leaving term

      REAL(fpk) ::  SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk) ::  SLTERM_F_0      ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) ::  USER_SLTERM_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  ----------------
!  Standard Outputs
!  ----------------

!  Main
!  ====

!  Intensity Results at all angles and optical depths

      REAL(fpk) :: INTENSITY ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Results for mean-value output

      REAL(fpk) :: MEAN_INTENSITY (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
      REAL(fpk) :: FLUX_INTEGRAL  (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)

!  Direct-beam contributions output separately, 26 May 11, 24 August 2011

      REAL(fpk) :: DNMEAN_DIRECT (MAX_USER_LEVELS, MAXBEAMS)
      REAL(fpk) :: DNFLUX_DIRECT (MAX_USER_LEVELS, MAXBEAMS)

!  Ancillary Output
!  ================

!  Fourier numbers used

      INTEGER      :: FOURIER_SAVED ( MAXBEAMS ) 

!  Number of geometries (bookkeeping output)

      INTEGER      :: N_GEOMETRIES

!  Exception handling
!  ==================

!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER              :: STATUS_INPUTCHECK
      INTEGER              :: NCHECKMESSAGES
      CHARACTER(Len=120)   :: CHECKMESSAGES(0:MAX_MESSAGES)
      CHARACTER(Len=120)   :: ACTIONS (0:MAX_MESSAGES)

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER              :: STATUS_CALCULATION
      CHARACTER(Len=120)   :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  -------------------------
!  Linearized Inputs - Fixed
!  -------------------------

!  Control
!  =======

!  Control linearization

      LOGICAL :: DO_COLUMN_LINEARIZATION
      LOGICAL :: DO_PROFILE_LINEARIZATION
      LOGICAL :: DO_SURFACE_LINEARIZATION

!  Control for atmospheric linearizations, layer by layer

      LOGICAL :: LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER :: LAYER_VARY_NUMBER(MAXLAYERS)

!  Total number of column Jacobians

      INTEGER :: N_TOTALCOLUMN_WFS

!  Total number of profile Jacobians. Local only, not a proxy. Version 3.7

      INTEGER :: N_TOTALPROFILE_WFS

!  Total number of Surface Jacobians

      INTEGER :: N_SURFACE_WFS

!  Control for  Blackbody Jacobians, New 18 March 2014

      LOGICAL :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  Optical
!  =======

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk) :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS)
      REAL(fpk) :: L_DELTAU_VERT_INPUT(MAX_ATMOSWFS,MAXLAYERS)
      REAL(fpk) :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS)

!  -------------------------
!  Linearized Supplement I/O
!  -------------------------

!  BRDF Inputs
!  ===========

!  Linearized Exact (direct bounce) BRDF

      REAL(fpk) :: LS_EXACTDB_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized Fourier components of BRDF, in the following order

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: LS_BRDF_F        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

      REAL(fpk) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Emissivity

      REAL(fpk) :: LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTREAMS )
      REAL(fpk) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  SS and DB I/O
!  =============

!  Single scatter column weighting functions at user angles

      REAL(fpk) :: COLUMNWF_SS ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Direct Bounce column weighting functions at user angles

      REAL(fpk) :: COLUMNWF_DB ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Single scatter profile weighting functions at user angles

      REAL(fpk) :: PROFILEWF_SS ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Direct Bounce profile weighting functions at user angles

      REAL(fpk) :: PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Direct Bounce surface weighting functions at user angles

      REAL(fpk) :: SURFACEWF_DB  ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Surface-Leaving Inputs, 22 Aug 12
!  =================================

!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
      LOGICAL          ::   DO_SLEAVE_WFS
      INTEGER          ::   N_SLEAVE_WFS

      REAL(fpk) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(fpk) :: LSSL_SLTERM_USERANGLES ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      REAL(fpk) :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: LSSL_USER_SLTERM_F_0   ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  ------------------
!  Linearized Outputs
!  ------------------

!  Profile Jacobians
!  =================

!  Profile weighting functions

      REAL(fpk) :: PROFILEWF &
        ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Mean-intensity and flux weighting functions

      REAL(fpk) :: MINT_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk) :: FLUX_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

!  Surface Jacobians
!  =================

!  Surface weighting functions

      REAL(fpk) :: SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Mean-intensity and flux weighting functions

      REAL(fpk) :: MINT_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk) :: FLUX_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

!  BLACKBODY Jacobians, New 18 March 2014
!  ======================================

!  Postprocessed Jacobians.

      REAL(fpk) :: ABBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS)
      REAL(fpk) :: SBBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS)

!  Flux Jacobians.

      REAL(fpk) :: ABBWFS_FLUXES ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS)
      REAL(fpk) :: SBBWFS_FLUXES ( MAX_USER_LEVELS, 2, MAX_DIRECTIONS)

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

! ######################################################################
!                       Local Arguments
! ######################################################################

!  Local bookkeeping
!  ----------------

!               Intent(In) To the Fourier  routine
!               Intent(In) To the Converge routine

!  Mode of operation

      LOGICAL         :: DO_MSMODE_LIDORT

! @@@@@@@@@@@@@ Robfix, 13 January 2012.
!               Add MSMODE_THERMAL flag
      logical         :: DO_MSMODE_THERMAL
! @@@@@@@@@@@@@ End Robfix, 13 January 2012.

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER         :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER         :: NSTREAMS_2
      INTEGER         :: NTOTAL
      INTEGER         :: N_SUBDIAG, N_SUPDIAG

!  Quadrature weights and abscissae, and product

      REAL(fpk)       :: QUAD_STREAMS (MAXSTREAMS)
      REAL(fpk)       :: QUAD_WEIGHTS (MAXSTREAMS)
      REAL(fpk)       :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      REAL(fpk)       :: USER_ANGLES  (MAX_USER_STREAMS)
      REAL(fpk)       :: USER_STREAMS  (MAX_USER_STREAMS)
      REAL(fpk)       :: USER_SECANTS  (MAX_USER_STREAMS)

!  Output optical depth masks and indices

      LOGICAL         :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER         :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER         :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER         :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

!  Off-grid optical depths (values, masks, indices)

      INTEGER         :: N_PARTLAYERS
      INTEGER         :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)
      REAL(fpk)       :: PARTLAYERS_VALUES       (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL         :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL         :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Indexing numbers

      INTEGER         :: N_VIEWING

!  Offsets for geometry indexing

      INTEGER         :: IBOFF(MAXBEAMS)
      INTEGER         :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(fpk)       :: SUNLAYER_COSINES(MAXLAYERS,MAXBEAMS)

!  Local solar zenith angles Cosines (regular case)

      REAL(fpk)       :: BEAM_COSINES(MAXBEAMS)

!  Solar beam flags (always internal)

      LOGICAL         :: DO_MULTIBEAM (MAXBEAMS,0:MAXFOURIER)

!  Number of directions (1 or 2) and directional array

      INTEGER         :: N_DIRECTIONS
      INTEGER         :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Number of convergence tests

      INTEGER         :: N_CONVTESTS

!  Adjusted geometries. New, 2007.
!  -------------------------------

!               Intent(Out) from the Adjust-geometry routine
!               Intent(In)  to   the Correction      routine

      REAL(fpk)       :: USER_ANGLES_ADJUST (MAX_USER_STREAMS)
      REAL(fpk)       :: BEAM_SZAS_ADJUST   (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)
      REAL(fpk)       :: USER_RELAZMS_ADJUST(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Arrays for setups and Corrections
!  ---------------------------------

!               Intent(In) To the Fourier routine

!  Local flags for the solution saving option

      INTEGER         :: LAYER_MAXMOMENTS (MAXLAYERS)

!  Initial transmittances * (secants)

      REAL(fpk)       :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Saved arrays for truncation factor and Delta-M scaling

      REAL(fpk)       :: TRUNC_FACTOR(MAXLAYERS)
      REAL(fpk)       :: FAC1(MAXLAYERS)

!  Derived Slant optical thickness inputs

      REAL(fpk)       :: TAUSLANT    ( 0:MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: TAUGRID     ( 0:MAXLAYERS )
      REAL(fpk)       :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: DELTAU_SLANT_UNSCALED( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Scaled SSAs and phase function moments

      REAL(fpk)       :: OMEGA_TOTAL    ( MAXLAYERS )
      REAL(fpk)       :: PHASMOMS_TOTAL ( 0:MAXMOMENTS, MAXLAYERS )

!  Linearized Input optical depths after delta-M scaling

      REAL(fpk)       :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk)       :: L_PHASMOMS_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS )

!  Derived Slant optical thickness inputs

      REAL(fpk)       :: L_DELTAU_SLANT (MAX_ATMOSWFS,MAXLAYERS,MAXLAYERS,MAXBEAMS)

!  Linearized truncation factor

      REAL(fpk)       :: L_TRUNC_FACTOR(MAX_ATMOSWFS,MAXLAYERS)

!  L'Hopital's rule logical variables

      LOGICAL         :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Coefficient functions for user-defined angles

      REAL(fpk)       :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk)       :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Fourier component output
!  ------------------------

!               Intent(Out) from the Fourier routine
!               Intent(in)  to   the Converge routine

!  Fourier comonents User-defined solutions

      REAL(fpk)       :: INTENSITY_F &
         (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Fourier-component Profile weighting functions at user angles

      REAL(fpk)       :: PROFILEWF_F &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_USER_STREAMS,MAXBEAMS, MAX_DIRECTIONS)

!  Fourier-component surface weighting functions at user angles

      REAL(fpk)       :: SURFACEWF_F ( MAX_SURFACEWFS,  MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Help arrays from the SS/DB correction routines
!  ==============================================

!  Saved Legendre polynomials. Rob Fix 11/18/14. Removed from here.
!      REAL(fpk)       :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
!      REAL(fpk)       :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
!  Saved TMS (Nakajima-Tanaka) factor.Rob Fix 11/18/14. Removed from here.
!      REAL(fpk)       :: TMS(MAXLAYERS)
!  Local truncation factors for additional DELTAM scaling .Rob Fix 11/18/14. Removed from here.
!      REAL(fpk)       :: SSFDEL ( MAXLAYERS )
!  Exact Phase function calculations. Rob Fix 11/18/14. Removed from here.
!      REAL(fpk)       :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
!      REAL(fpk)       :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)

!  Cumulative single scatter source terms

      REAL(fpk)       :: SS_CUMSOURCE_UP(MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(fpk)       :: SS_CUMSOURCE_DN(MAX_GEOMETRIES,0:MAXLAYERS)

!  Atmospheric attenuation before reflection

      REAL(fpk)       :: ATTN_DB_SAVE(MAX_GEOMETRIES)

!  Exact direct beam source terms

      REAL(fpk)       :: EXACTDB_SOURCE(MAX_GEOMETRIES)

!  Cumulative direct bounce source terms

      REAL(fpk)       :: DB_CUMSOURCE(MAX_GEOMETRIES,0:MAXLAYERS)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(fpk)       :: BOA_ATTN(MAX_GEOMETRIES)

!  Outgoing sphericity stuff
!  Whole and part-layer LOS transmittance factors

      REAL(fpk)       :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(fpk)       :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Outgoing sphericity stuff
!    Linearized output. Output required for later on.

      REAL(fpk)       :: L_UP_LOSTRANS   (MAXLAYERS,     MAX_ATMOSWFS,MAX_GEOMETRIES)
      REAL(fpk)       :: L_UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(fpk)       :: L_BOA_ATTN(MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Arrays required at the Top level
!  ================================

!               Intent(In) To the Fourier routine

!  Input optical properties after delta-M scaling

      REAL(fpk)       :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk)       :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk)       :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )

!  Linearized input optical properties after delta-M scaling

      LOGICAL         :: DO_PHFUNC_VARIATION ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk)       :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk)       :: L_OMEGA_MOMS   ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )

!  Local input solar zenith angles by levels
!  ( Only required for refractive geometry attenuation of the solar beam)
!  These will be set internally if the refraction flag is set.

      REAL(fpk)       :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER         :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk)       :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )

!  Solar beam attenuation

      REAL(fpk)       :: SOLAR_BEAM_OPDEP ( MAXBEAMS )

!  Rob fix 11/27/14. Proxy for new output.
!     NOTE. SOLARBEAM_BOATRANS is computed with Unscaled optical depths, SOLAR_BEAM_OPDEP with scaled ODs.
!           [ They are the same if there is no Deltam-scaling]

      REAL(fpk)       :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk)       :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk)       :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Reflectance flags

      LOGICAL         :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local flags for the solution saving option

      LOGICAL         :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Local flags,  BVP telescoping enhancement

      LOGICAL         :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Masking for regular case. Required again for linearization

      INTEGER         :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER         :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Telescoping initial flag (modified argument), Layer bookkeeping
!  Number of telescoped layers, active layers,  Size of BVP matrix 

      LOGICAL         :: DO_BVTEL_INITIAL
      INTEGER         :: NLAYERS_TEL
      INTEGER         :: ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER         :: N_BVTELMATRIX_SIZE

!  Set up for band matrix compression

      INTEGER         :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)
      INTEGER         :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Transmittance Setups
!  --------------------

!               Intent(In) To the Fourier routine

!  Discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk)       :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage

      REAL(fpk)       :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk)       :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      REAL(fpk)       :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage 

      REAL(fpk)       :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk)       :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk)       :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!               Intent(In) To the Fourier routine

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(fpk)       :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk)       :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Partial layer multipliers

      REAL(fpk)       :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk)       :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  LINEARIZED Arrays required at the Top level
!  ===========================================

!  Linearized Transmittance Setups
!  -------------------------------

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk)       :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk)       :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized transmittances, solar beam

      REAL(fpk)       :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Beam multipliers
!  ---------------------------

!  Linearized whole layer multipliers

      REAL(fpk)       :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(fpk)       :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized help arrays for Green's function

      REAL(fpk)       :: HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Thermal Setup outputs
!  ---------------------

!  Optical depth powers

      REAL(fpk)       :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk)       :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Thermal coefficients

      REAL(fpk)       :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk)       :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk)       :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )

      REAL(fpk)       :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized optical depth powers

      REAL(fpk)       :: L_DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Linearized Thermal coefficients

      REAL(fpk)       :: L_THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Linearized Tranmsittance solutions

      REAL(fpk)       :: L_T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      REAL(fpk)       :: L_T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Help variables
!  ==============

!  Local flags

      LOGICAL         :: LOCAL_DO_NO_AZIMUTH
      LOGICAL         :: SAVE_DO_NO_AZIMUTH
      LOGICAL         :: LOCAL_ITERATION
      LOGICAL         :: ADJUST_SURFACE
      LOGICAL         :: DO_PROCESS_FOURIER

!  Local fourier index

      INTEGER         :: FOURIER_COMPONENT
      INTEGER         :: N_FOURIER_COMPONENTS

!  Local integers

      INTEGER         :: UA, UM, IB, TESTCONV, L, NSOURCES, LUM, LUA, NMOMSINP
      INTEGER         :: LOCAL_N_USERAZM, STATUS_SUB

!  Flux multiplier

      REAL(fpk)       :: SS_FLUX_MULTIPLIER

!  Total number of surface Jacobians (both surface and surface-leaving)

      INTEGER         :: N_TOTALSURFACE_WFS

!  Modified eradius

      REAL(fpk)       :: MODIFIED_ERADIUS

!  Fourier cosine arguments

      REAL(fpk)       :: AZM_ARGUMENT, DFC
      REAL(fpk)       :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Local convergence control


      INTEGER         :: IBEAM_COUNT, IBEAM
      LOGICAL         :: BEAM_ITERATION ( MAXBEAMS )
      INTEGER         :: BEAM_TESTCONV  ( MAXBEAMS )

!  Local linearization variables

      INTEGER         :: LAYER_TO_VARY, N_PARAMETERS

!  Local error handling

      LOGICAL         :: FAIL

!  Test variables

!      LOGICAL         :: DO_FDTEST=.FALSE.

      LOGICAL         :: DO_DEBUG_INPUT=.FALSE.
!      LOGICAL         :: DO_DEBUG_INPUT=.TRUE.

!  For debug
!      INTEGER          :: I

!  Initialize some variables
!  -------------------------

!mick debug
!write(*,*) 'inside LIDORT_LPS_master'

!  Main status

      STATUS_CALCULATION = LIDORT_SUCCESS
      STATUS_INPUTCHECK  = LIDORT_SUCCESS

!mick fix 6/29/11 - initialize "Input checks"
!  Input checks

      NCHECKMESSAGES = 0
      CHECKMESSAGES  = ' '
      ACTIONS        = ' '

!  Model calculation

      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  Local user indices

      LUM = 1
      LUA = 1

!  Check input dimensions
!  ----------------------

!  Regular

      CALL LIDORT_CHECK_INPUT_DIMS &
      ( LIDORT_FixIn, LIDORT_ModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  Linearized

      CALL LIDORT_L_CHECK_INPUT_DIMS &
      ( LIDORT_LinFixIn, LIDORT_LinModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  Fixed Boolean inputs

      DO_FULLRAD_MODE      = LIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE
      DO_SSCORR_TRUNCATION = LIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION

      DO_SS_EXTERNAL       = LIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL
      DO_SSFULL            = LIDORT_FixIn%Bool%TS_DO_SSFULL

      DO_THERMAL_EMISSION  = LIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION
      DO_SURFACE_EMISSION  = LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION

      DO_PLANE_PARALLEL    = LIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL
      DO_BRDF_SURFACE      = LIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE

      DO_UPWELLING         = LIDORT_FixIn%Bool%TS_DO_UPWELLING
      DO_DNWELLING         = LIDORT_FixIn%Bool%TS_DO_DNWELLING

!  New 17 May 2012

      DO_SURFACE_LEAVING   = LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC      = LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  Fixed control inputs. Taylor Ordering control, added 10/10/13
!     TCUTOFF introduced 15 May 2015

      TAYLOR_ORDER     = LIDORT_FixIn%Cont%TS_TAYLOR_ORDER
      NSTREAMS         = LIDORT_FixIn%Cont%TS_NSTREAMS
      NLAYERS          = LIDORT_FixIn%Cont%TS_NLAYERS
      NFINELAYERS      = LIDORT_FixIn%Cont%TS_NFINELAYERS
      N_THERMAL_COEFFS = LIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS
      LIDORT_ACCURACY  = LIDORT_FixIn%Cont%TS_LIDORT_ACCURACY
      TCUTOFF          = LIDORT_FixIn%Cont%TS_THERMAL_CUTOFF

!  Fixed beam inputs

      FLUX_FACTOR         = LIDORT_FixIn%Sunrays%TS_FLUX_FACTOR

!  Fixed user value inputs
!     Observation-Geometry input control. New 25 October 2012, R. Spurr, RT SOLUTIONS Inc.
!     Moved to the Modified structure, 10/25/12
!     N_USER_STREAMS      = LIDORT_FixIn%UserVal%TS_N_USER_STREAMS

      N_USER_LEVELS       = LIDORT_FixIn%UserVal%TS_N_USER_LEVELS

!  Fixed Chapman function inputs

      HEIGHT_GRID(0:NLAYERS)      = LIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:NLAYERS)
      PRESSURE_GRID(0:NLAYERS)    = LIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:NLAYERS)
      TEMPERATURE_GRID(0:NLAYERS) = LIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(0:NLAYERS)
      FINEGRID(1:NLAYERS)         = LIDORT_FixIn%Chapman%TS_FINEGRID(1:NLAYERS)

      RFINDEX_PARAMETER   = LIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER

!  Fixed optical inputs

      NMOMSINP = LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT
      DELTAU_VERT_INPUT(1:NLAYERS)               = LIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(1:NLAYERS)
      PHASMOMS_TOTAL_INPUT(0:NMOMSINP,1:NLAYERS) = LIDORT_FixIn%Optical%TS_PHASMOMS_TOTAL_INPUT(0:NMOMSINP,1:NLAYERS)
      THERMAL_BB_INPUT(0:NLAYERS)                = LIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(0:NLAYERS)

      LAMBERTIAN_ALBEDO = LIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO
      SURFACE_BB_INPUT  = LIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT

!  Modified Boolean inputs

      DO_SSCORR_NADIR        = LIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR
      DO_SSCORR_OUTGOING     = LIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING

      DO_DOUBLE_CONVTEST     = LIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST

      DO_SOLAR_SOURCES       = LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES

      DO_REFRACTIVE_GEOMETRY = LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      DO_CHAPMAN_FUNCTION    = LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION

      DO_RAYLEIGH_ONLY       = LIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
      DO_ISOTROPIC_ONLY      = LIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY
      DO_NO_AZIMUTH          = LIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH
      DO_ALL_FOURIER         = LIDORT_ModIn%MBool%TS_DO_ALL_FOURIER

      DO_DELTAM_SCALING      = LIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING

      DO_SOLUTION_SAVING     = LIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING
      DO_BVP_TELESCOPING     = LIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING

      DO_USER_STREAMS        = LIDORT_ModIn%MBool%TS_DO_USER_STREAMS

      DO_ADDITIONAL_MVOUT    = LIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      DO_MVOUT_ONLY          = LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

      DO_THERMAL_TRANSONLY   = LIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      DO_OBSERVATION_GEOMETRY   = LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY

!  Modified control inputs

      NMOMENTS_INPUT = LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT

!  Modified beam inputs

      NBEAMS              = LIDORT_ModIn%MSunrays%TS_NBEAMS
      BEAM_SZAS(1:NBEAMS) = LIDORT_ModIn%MSunrays%TS_BEAM_SZAS(1:NBEAMS)

!  Modified User value inputs
!    * N_USER_STREAMS is moved from the Fixed input. 10/25/12
!    * Observational Geometry inputs introduced 10/25/12

      N_USER_RELAZMS = LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS(1:N_USER_RELAZMS)      = &
        LIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS)

      N_USER_STREAMS = LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS
      USER_ANGLES_INPUT(1:N_USER_STREAMS) = &
        LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT(1:N_USER_STREAMS)

      USER_LEVELS(1:N_USER_LEVELS)        = &
        LIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)

      N_USER_OBSGEOMS  = LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS

      USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)  = &
         LIDORT_ModIn%MUserVal%TS_USER_OBSGEOM_INPUT(1:N_USER_OBSGEOMS,1)
      USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)  = &
         LIDORT_ModIn%MUserVal%TS_USER_OBSGEOM_INPUT(1:N_USER_OBSGEOMS,2)
      USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)  = &
         LIDORT_ModIn%MUserVal%TS_USER_OBSGEOM_INPUT(1:N_USER_OBSGEOMS,3)

      GEOMETRY_SPECHEIGHT = LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT

!  Modified Chapman function inputs

      !CHAPMAN_FACTORS     = LIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS
      EARTH_RADIUS        = LIDORT_ModIn%MChapman%TS_EARTH_RADIUS

!  Modified optical inputs

      OMEGA_TOTAL_INPUT(1:NLAYERS) = LIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS)

!  BRDF Supplement Inputs

      EXACTDB_BRDFUNC(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
        LIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC&
          (1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)

      BRDF_F_0(0:MAXMOMENTS,1:NSTREAMS,1:NBEAMS)              = &
        LIDORT_Sup%BRDF%TS_BRDF_F_0(0:MAXMOMENTS,1:NSTREAMS,1:NBEAMS)
      BRDF_F(0:MAXMOMENTS,1:NSTREAMS,1:NSTREAMS)              = &
        LIDORT_Sup%BRDF%TS_BRDF_F(0:MAXMOMENTS,1:NSTREAMS,1:NSTREAMS)

      USER_BRDF_F_0(0:MAXMOMENTS,1:N_USER_STREAMS,1:NBEAMS)   = &
        LIDORT_Sup%BRDF%TS_USER_BRDF_F_0&
          (0:MAXMOMENTS,1:N_USER_STREAMS,1:NBEAMS)
      USER_BRDF_F(0:MAXMOMENTS,1:N_USER_STREAMS,1:NSTREAMS)   = &
        LIDORT_Sup%BRDF%TS_USER_BRDF_F&
          (0:MAXMOMENTS,1:N_USER_STREAMS,1:NSTREAMS)

      EMISSIVITY(1:NSTREAMS)            = LIDORT_Sup%BRDF%TS_EMISSIVITY(1:NSTREAMS)
      USER_EMISSIVITY(1:N_USER_STREAMS) = LIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1:N_USER_STREAMS)

!  Surface-Leaving Supplement inputs
!    (This code introduced 17 May 2012)

      SLTERM_ISOTROPIC(1:NBEAMS)                                    = &
        LIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NBEAMS)
      SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
        LIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)

      SLTERM_F_0(0:MAXMOMENTS,1:NSTREAMS,1:NBEAMS)            = &
        LIDORT_Sup%SLEAVE%TS_SLTERM_F_0(0:MAXMOMENTS,1:NSTREAMS,1:NBEAMS)
      USER_SLTERM_F_0(0:MAXMOMENTS,1:N_USER_STREAMS,1:NBEAMS) = &
        LIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(0:MAXMOMENTS,1:N_USER_STREAMS,1:NBEAMS)

!  New 12 March 2012
!   IF SS results already available copy them !
!  SS Supplement Inputs

      IF ( DO_SS_EXTERNAL ) THEN
         INTENSITY_SS = LIDORT_Sup%SS%TS_INTENSITY_SS
         INTENSITY_DB = LIDORT_Sup%SS%TS_INTENSITY_DB
      ENDIF

!  Linearized control inputs

      DO_COLUMN_LINEARIZATION      = &
        LIDORT_LinFixIn%Cont%TS_DO_COLUMN_LINEARIZATION
      DO_PROFILE_LINEARIZATION     = &
        LIDORT_LinFixIn%Cont%TS_DO_PROFILE_LINEARIZATION
      DO_SURFACE_LINEARIZATION     = &
        LIDORT_LinFixIn%Cont%TS_DO_SURFACE_LINEARIZATION

      LAYER_VARY_FLAG(1:NLAYERS)   = &
        LIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)
      LAYER_VARY_NUMBER(1:NLAYERS) = &
        LIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS)

!  Local value of number of total profile WFs

      N_TOTALPROFILE_WFS = 0
      DO L = 1, NLAYERS
        N_TOTALPROFILE_WFS  = MAX(N_TOTALPROFILE_WFS,LAYER_VARY_NUMBER(L))
      ENDDO
!      N_TOTALPROFILE_WFS = LIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS

      N_TOTALCOLUMN_WFS  = LIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS
      N_SURFACE_WFS      = LIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS

!  Linearized optical inputs

      L_DELTAU_VERT_INPUT(1:N_TOTALPROFILE_WFS,1:NLAYERS) = &
        LIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1:N_TOTALPROFILE_WFS,1:NLAYERS)
      L_OMEGA_TOTAL_INPUT(1:N_TOTALPROFILE_WFS,1:NLAYERS) = &
        LIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1:N_TOTALPROFILE_WFS,1:NLAYERS)
      L_PHASMOMS_TOTAL_INPUT(1:N_TOTALPROFILE_WFS,0:NMOMENTS_INPUT,1:NLAYERS) = &
        LIDORT_LinFixIn%Optical%TS_L_PHASMOMS_TOTAL_INPUT(1:N_TOTALPROFILE_WFS,0:NMOMENTS_INPUT,1:NLAYERS)

!  Linearized BRDF Inputs

      LS_EXACTDB_BRDFUNC(1:N_SURFACE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
        LIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC(1:N_SURFACE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)

      LS_BRDF_F_0(1:N_SURFACE_WFS,0:MAXMOMENTS,1:NSTREAMS,1:NBEAMS) = &
        LIDORT_LinSup%BRDF%TS_LS_BRDF_F_0(1:N_SURFACE_WFS,0:MAXMOMENTS,1:NSTREAMS,1:NBEAMS)

      LS_BRDF_F(1:N_SURFACE_WFS,0:MAXMOMENTS,1:NSTREAMS,1:NSTREAMS) = &
        LIDORT_LinSup%BRDF%TS_LS_BRDF_F(1:N_SURFACE_WFS,0:MAXMOMENTS,1:NSTREAMS,1:NSTREAMS)

      LS_USER_BRDF_F_0(1:N_SURFACE_WFS,0:MAXMOMENTS,1:N_USER_STREAMS,1:NBEAMS) = &
        LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0(1:N_SURFACE_WFS,0:MAXMOMENTS,1:N_USER_STREAMS,1:NBEAMS)

      LS_USER_BRDF_F(1:N_SURFACE_WFS,0:MAXMOMENTS,1:N_USER_STREAMS,1:NSTREAMS) = &
        LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F(1:N_SURFACE_WFS,0:MAXMOMENTS,1:N_USER_STREAMS,1:NSTREAMS)

      LS_EMISSIVITY(1:N_SURFACE_WFS,1:NSTREAMS)            = &
        LIDORT_LinSup%BRDF%TS_LS_EMISSIVITY(1:N_SURFACE_WFS,1:NSTREAMS)

      LS_USER_EMISSIVITY(1:N_SURFACE_WFS,1:N_USER_STREAMS) = &
        LIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY(1:N_SURFACE_WFS,1:N_USER_STREAMS)

!  New 12 March 2012
!   IF SS results already available copy them !
!  Linearized SS Inputs

      IF ( DO_SS_EXTERNAL ) THEN
         !SS atmosphere weighting functions

         COLUMNWF_SS   = LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS
         COLUMNWF_DB   = LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB

         PROFILEWF_SS  = LIDORT_LinSup%SS%Atmos%TS_PROFILEWF_SS
         PROFILEWF_DB  = LIDORT_LinSup%SS%Atmos%TS_PROFILEWF_DB

         !SS surface weighting functions

         !SURFACEWF_SS  = LIDORT_LinSup%SS%Surf%TS_SURFACEWF_SS
         SURFACEWF_DB  = LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@

      DO_SLEAVE_WFS  = LIDORT_LinFixIn%Cont%TS_DO_SLEAVE_WFS
      N_SLEAVE_WFS   = LIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS

      LSSL_SLTERM_ISOTROPIC  = LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC
      LSSL_SLTERM_USERANGLES = LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES
      LSSL_SLTERM_F_0        = LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0
      LSSL_USER_SLTERM_F_0   = LIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0

      N_TOTALSURFACE_WFS = N_SURFACE_WFS + N_SLEAVE_WFS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  BlackBody Jacobian Flags, New, 18 March 2014

      DO_ATMOS_LBBF   = LIDORT_LinFixIn%Cont%TS_DO_ATMOS_LBBF
      DO_SURFACE_LBBF = LIDORT_LinFixIn%Cont%TS_DO_SURFACE_LBBF

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  LIDORT input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL LIDORT_DEBUG_INPUT_MASTER()
        IF (DO_COLUMN_LINEARIZATION  .OR. &
            DO_PROFILE_LINEARIZATION .OR. &
            DO_SURFACE_LINEARIZATION) &
          CALL LIDORT_DEBUG_LIN_INPUT_MASTER()
      END IF

!mick fix 6/29/11 - initialize output array entries
!  Main outputs. Direct_DN output also initialized (8/24/11)

      INTENSITY(1:N_USER_LEVELS,:,:) = ZERO
      MEAN_INTENSITY(1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO
      FLUX_INTEGRAL(1:N_USER_LEVELS,1:NBEAMS,:)   = ZERO

      DNMEAN_DIRECT(1:N_USER_LEVELS,1:NBEAMS) = ZERO
      DNFLUX_DIRECT(1:N_USER_LEVELS,1:NBEAMS) = ZERO

!  New 15 March 2012

      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         INTENSITY_SS = ZERO
         INTENSITY_DB = ZERO
      ENDIF

      FOURIER_SAVED(1:NBEAMS) = 0
      N_GEOMETRIES            = 0

      PROFILEWF(1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,:,:) = ZERO
      MINT_PROFILEWF(1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO
      FLUX_PROFILEWF(1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO

      SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,:) = ZERO
      MINT_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO
      FLUX_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO

!  New 15 March 2012

      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         !SS atmosphere profile weighting functions
         COLUMNWF_SS  = ZERO
         COLUMNWF_DB  = ZERO
         PROFILEWF_SS = ZERO
         PROFILEWF_DB = ZERO
         !SS surface weighting functions
         !SURFACEWF_SS = ZERO
         SURFACEWF_DB = ZERO
      ENDIF

!  New 18 March 2014. BLACKBODY Linearization
!  ------------------------------------------

      ABBWFS_JACOBIANS(1:N_USER_LEVELS,:,:,:)  = ZERO
      ABBWFS_FLUXES   (1:N_USER_LEVELS,:,:,:)  = ZERO
      SBBWFS_JACOBIANS(1:N_USER_LEVELS,:,:)    = ZERO
      SBBWFS_FLUXES   (1:N_USER_LEVELS,:,:)    = ZERO

!  Single scatter correction: flux multiplier
!    Now always F / 4pi

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

!  %% Major revision of I/O output list, 25 October 2012.
!     ---- Observational Geometry control, New, 25 October 2012
!     ---- Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, DO_USER_STREAMS, DO_NO_AZIMUTH

!  %% DO_FULLRAD_MODE argument added (Second line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged

!  %% Revision, 10/10/13. Include TAYLOR_ORDER parameter, revise argument list

      CALL LIDORT_CHECK_INPUT                                                 &
      ( DO_UPWELLING, DO_DNWELLING, DO_THERMAL_EMISSION,                      & ! Fixed inputs
        DO_FULLRAD_MODE, DO_SSFULL, DO_SS_EXTERNAL, DO_PLANE_PARALLEL,        & ! Fixed inputs 
        DO_SOLAR_SOURCES, DO_OBSERVATION_GEOMETRY, DO_THERMAL_TRANSONLY,      & ! Input, possibly modified
        DO_USER_STREAMS, DO_NO_AZIMUTH, DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,  & ! Input, possibly modified
        DO_DELTAM_SCALING, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,               & ! Input, possibly modified
        DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION, DO_MVOUT_ONLY,           & ! Input, possibly modified
        DO_ADDITIONAL_MVOUT, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,          & ! Input, possibly modified
        TAYLOR_ORDER, NSTREAMS, N_USER_LEVELS, NLAYERS, HEIGHT_GRID,          & ! Fixed Inputs
        N_USER_OBSGEOMS, USER_OBSGEOMS, NBEAMS, BEAM_SZAS,                    & ! Input, possibly modified
        N_USER_STREAMS, N_USER_RELAZMS, USER_ANGLES_INPUT, USER_RELAZMS,      & ! Input, possibly modified
        NMOMENTS_INPUT, USER_LEVELS, EARTH_RADIUS, GEOMETRY_SPECHEIGHT,       & ! Input, possibly modified
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )                    ! Output

!  Exception handling

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ELSE IF ( STATUS_SUB .EQ. LIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = LIDORT_WARNING
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!mick fix 11/8/2013 - these outputs set at the end
        !LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        !LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        !LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
      ENDIF

!  Check input optical values (IOPs, albedo)

      CALL LIDORT_CHECK_INPUT_OPTICAL                          &
         ( NLAYERS, LAMBERTIAN_ALBEDO, DELTAU_VERT_INPUT,      & ! Input
           OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,            & ! Input
           STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS)   ! output

!  Exception handling

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  if there's no azimuth dependence, just do one value in azimuth loop

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
      ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

!  Number of sources

      IF ( DO_SOLAR_SOURCES ) THEN
        NSOURCES = NBEAMS
      ELSE
        NSOURCES = 1
      ENDIF

!  save some offsets for indexing geometries
!   This section revised for the Observational Geometry option

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        N_VIEWING    = N_USER_STREAMS * LOCAL_N_USERAZM
        N_GEOMETRIES = NSOURCES * N_VIEWING
        DO IBEAM = 1, NBEAMS
          IBOFF(IBEAM) = N_VIEWING * ( IBEAM - 1 )
          DO UM = 1, N_USER_STREAMS
            UMOFF(IBEAM,UM) = IBOFF(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
          END DO
        END DO
      ELSE
        N_VIEWING    = N_USER_OBSGEOMS
        N_GEOMETRIES = N_USER_OBSGEOMS
        IBOFF = 0 ; UMOFF = 0
      ENDIF

!  Geometry adjustment
!  -------------------

!  Adjust surface condition

      ADJUST_SURFACE = .FALSE.
      IF ( DO_SSCORR_OUTGOING ) THEN
        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
         ADJUST_SURFACE = .TRUE.
        ENDIF
      ENDIF

!  Perform adjustment

      MODIFIED_ERADIUS = EARTH_RADIUS + GEOMETRY_SPECHEIGHT
!mick hold - 9/26/2012
!      IF ( DO_SOLAR_SOURCES ) THEN

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          CALL MULTI_OUTGOING_ADJUSTGEOM                              &
         ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS,              & ! Input
           N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                & ! Input
           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,    & ! Input
           USER_ANGLES_INPUT,  BEAM_SZAS, USER_RELAZMS,               & ! Input
           USER_ANGLES_ADJUST, BEAM_SZAS_ADJUST, USER_RELAZMS_ADJUST, & ! Output
           FAIL, MESSAGE, TRACE_1 )                                     ! Output
      ELSE
          CALL OBSGEOM_OUTGOING_ADJUSTGEOM                            &
         ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS,              & ! Input
           N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                & ! Input
           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,    & ! Input
           USER_ANGLES_INPUT,  BEAM_SZAS, USER_RELAZMS,               & ! Input
           USER_ANGLES_ADJUST, BEAM_SZAS_ADJUST, USER_RELAZMS_ADJUST, & ! Output
           FAIL, MESSAGE, TRACE_1 )                                     ! Output
      ENDIF

!      ELSE
!        CALL LOSONLY_OUTGOING_ADJUSTGEOM                           &
!         ( MAX_USER_STREAMS, N_USER_STREAMS,                       & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE, & ! Input
!           USER_ANGLES_INPUT,                                      & ! Input
!           USER_ANGLES_ADJUST,                                     & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                  ! Output
!      ENDIF

!  Exception handling

      if ( fail ) then
        TRACE_2 = ' Failure in multi_outgoing_adjustgeom'
        TRACE_3 = ' ** Called in LIDORT_LPS_MASTER'
        STATUS_CALCULATION = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
        LIDORT_Out%Status%TS_MESSAGE = MESSAGE
        LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
        LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
        LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
        return
      endif

!  Chapman function calculation
!  ----------------------------

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN

          CALL LIDORT_CHAPMAN                             &
          ( DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,    & ! Input
            NLAYERS, NBEAMS, FINEGRID, BEAM_SZAS,         & ! Input
            EARTH_RADIUS, RFINDEX_PARAMETER,              & ! Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, & ! Input
            CHAPMAN_FACTORS, SZA_LOCAL_INPUT,             & ! output
            FAIL, MESSAGE, TRACE_1 )                        ! output

          IF (FAIL) THEN
            TRACE_2 = 'Failure in LIDORT_CHAPMAN'
            TRACE_3 = ' ** Called in LIDORT_LPS_MASTER'
            STATUS_CALCULATION = LIDORT_SERIOUS
            LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            LIDORT_Out%Status%TS_MESSAGE = MESSAGE
            LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF

        ENDIF
      ENDIF

!  Get derived inputs
!  ==================

!  Miscellaneous and layer input.
!   ( 6/21/10)  Important point: Must use ADJUSTED VZangles as input here.
!   ( 11/19/14) DISAGREE WITH THIS POINT------------------!!!!!!!!!!!!!!!!!

      CALL LIDORT_DERIVE_INPUT                                       &
      ( DO_FULLRAD_MODE, DO_USER_STREAMS,                            & ! Input
        DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY, DO_OBSERVATION_GEOMETRY,& ! Input
        DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,              & ! Input
        DO_UPWELLING, DO_DNWELLING, DO_REFRACTIVE_GEOMETRY,          & ! Input
        DO_THERMAL_TRANSONLY,                                        & ! Input
        DO_DOUBLE_CONVTEST, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,  & ! In/Out
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,                          & ! In/Out
        OMEGA_TOTAL_INPUT,                                           & ! In/Out
        NLAYERS, NSTREAMS, NBEAMS, NMOMENTS_INPUT,                   & ! Input
        N_USER_STREAMS, N_USER_RELAZMS, USER_ANGLES_INPUT,           & ! Input   FORMERLY USER_ANGLES_ADJUST
        N_USER_LEVELS,  USER_LEVELS, BEAM_SZAS, SZA_LOCAL_INPUT,     & ! Input
        DO_MSMODE_LIDORT, NMOMENTS, BEAM_COSINES, SUNLAYER_COSINES,  & ! Output
        N_CONVTESTS, N_DIRECTIONS, WHICH_DIRECTIONS,                 & ! Output
        NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG, N_PARTLAYERS,      & ! Output
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                    & ! Output
        USER_ANGLES, USER_STREAMS, USER_SECANTS,                     & ! Output
        PARTLAYERS_OUTFLAG,PARTLAYERS_OUTINDEX,PARTLAYERS_LAYERIDX,  & ! Output
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_VALUES,   & ! Output
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN )                       ! Output

! @@@@@@@@@@@@@ Robfix, 13 January 2012.
!               Set MSMODE_THERMAL flag
      DO_MSMODE_THERMAL = (.not.DO_FULLRAD_MODE) .and. &
             ( DO_SURFACE_EMISSION .and.DO_THERMAL_EMISSION )
! @@@@@@@@@@@@@ End Robfix, 13 January 2012.

!  #################
!  Set up operations
!  #################

!   MISCSETUPS (6 subroutines)  :
!   -----------------------------

!       Performance Setup,
!       Delta-M scaling,
!       average-secant formulation,
!       transmittance setup
!       Beam solution multipliers
!       Thermal setups (if flagged)

!  Performance set-up

      CALL LIDORT_PERFORMANCE_SETUP                                        &
        ( DO_SOLAR_SOURCES, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,           & ! input
          DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,                          & ! input
          DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                             & ! input
          NLAYERS, NMOMENTS, NMOMENTS_INPUT, PHASMOMS_TOTAL_INPUT,         & ! input
          LAYER_MAXMOMENTS, DO_LAYER_SCATTERING, BVP_REGULAR_FLAG,         & ! Output
          STATUS_SUB, MESSAGE, TRACE_1)                                      ! Output

!  Exception handling on this module
!   Even though this is a warning, must exit with Serious condition

      IF ( STATUS_SUB .EQ. LIDORT_WARNING ) THEN
        TRACE_2 = 'Failure in LIDORT_PERFORMANCE_SETUP'
        TRACE_3 = ' ** Called in LIDORT_LPS_MASTER'
        STATUS_CALCULATION = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
        LIDORT_Out%Status%TS_MESSAGE = MESSAGE
        LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
        LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
        LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
        RETURN
      ENDIF 

!  Delta-m scaling of input quantities

      CALL LIDORT_DELTAMSCALE                                        &
       ( DO_SOLAR_SOURCES, DO_DELTAM_SCALING,                        & ! input
         NLAYERS, N_PARTLAYERS, NMOMENTS, NBEAMS, N_USER_LEVELS,     & ! input
         PARTLAYERS_OUTFLAG, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, & ! input
         CHAPMAN_FACTORS, DELTAU_VERT_INPUT,                         & ! input
         OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                    & ! input
         DELTAU_VERT, PARTAU_VERT, TAUGRID,                          & ! Output
         OMEGA_TOTAL, OMEGA_MOMS, PHASMOMS_TOTAL,                    & ! Output
         FAC1, TRUNC_FACTOR,                                         & ! Output
         DELTAU_SLANT, DELTAU_SLANT_UNSCALED, SOLARBEAM_BOATRANS )     ! Output      ! Rob fix 11/27/14 added argument                       

!  Prepare quasi-spherical attenuation

      CALL LIDORT_QSPREP                                                &
       ( DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,   & ! input
         NLAYERS, NBEAMS, BEAM_COSINES, SUNLAYER_COSINES,               & ! input
         TAUGRID, DELTAU_VERT, DELTAU_SLANT,                            & ! input
         INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,   & ! Output
         TAUSLANT, SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM )            ! Output

!  Transmittances and Transmittance factors

      CALL LIDORT_PREPTRANS                                            &
        ( DO_SOLAR_SOURCES, DO_SOLUTION_SAVING, DO_USER_STREAMS,       & ! input
          DO_OBSERVATION_GEOMETRY,                                     & ! input
          NSTREAMS, N_USER_STREAMS, NBEAMS, NLAYERS, N_PARTLAYERS,     & ! input
          QUAD_STREAMS, DELTAU_VERT, PARTLAYERS_LAYERIDX, PARTAU_VERT, & ! input
          USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,        & ! input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,             & ! input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,              & ! Output
          T_DELT_MUBAR,   T_UTUP_MUBAR,   T_UTDN_MUBAR,                & ! Output
          T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,                & ! Output
          ITRANS_USERM  )                                                ! Output

!  Linearization of Delta-M scaled inputs

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_LA_DELTAMSCALE                                &
       ( DO_DELTAM_SCALING, NLAYERS, NMOMENTS, NBEAMS,            & ! Input
         CHAPMAN_FACTORS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,     & ! Input
         DELTAU_VERT_INPUT,    L_DELTAU_VERT_INPUT,               & ! Input
         OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,               & ! Input
         PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,            & ! Input
         OMEGA_MOMS, FAC1, TRUNC_FACTOR,                          & ! Input
         L_DELTAU_VERT, L_DELTAU_SLANT,                           & ! Output
         L_OMEGA_TOTAL, L_OMEGA_MOMS, L_PHASMOMS_TOTAL,           & ! Output
         L_TRUNC_FACTOR, DO_PHFUNC_VARIATION )                      ! Output
      ENDIF

!  Linearization of transmittances

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_LA_PREPTRANS                                 &
        ( DO_SOLUTION_SAVING, DO_USER_STREAMS,                   & ! Input
          NLAYERS, NSTREAMS, N_USER_STREAMS,                     & ! Input
          N_PARTLAYERS, PARTLAYERS_LAYERIDX,                     & ! Input
          QUAD_STREAMS, USER_SECANTS,                            & ! Input
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                    & ! Input
          DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,               & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,        & ! Input
          T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,          & ! Input
          L_T_DELT_DISORDS, L_T_DISORDS_UTUP,  L_T_DISORDS_UTDN, & ! Output
          L_T_DELT_USERM,   L_T_UTUP_USERM,    L_T_UTDN_USERM )    ! Output
      ENDIF

!  Linearization of pseudo-spherical setup

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_LP_PREPTRANS                                  &
        ( DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS,       & ! Input
          PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,                 & ! Input
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,  & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,        & ! Input
          T_DELT_MUBAR, T_UTDN_MUBAR,                             & ! Input
          LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR, LP_INITIAL_TRANS,    & ! Output
          LP_AVERAGE_SECANT, HELP_AQ, HELP_BQ )                     ! Output
      ENDIF

!   EMULT_MASTER  : Beam source function multipliers. Not required for the
!                  Full SS calculation in outgoing mode
!  Rob Fix, 10/10/13 - Introduce Taylor Order parameter

      IF ( DO_SOLAR_SOURCES  ) THEN
        IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

          CALL LIDORT_EMULT_MASTER  &
       ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,     & ! Input
         TAYLOR_ORDER, N_USER_STREAMS, NBEAMS, NLAYERS,           & ! Input
         N_PARTLAYERS, PARTLAYERS_LAYERIDX,                       & ! Input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! input
         DELTAU_VERT, PARTAU_VERT, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! input
         T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,            & ! input
         ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF,          & ! input
         SIGMA_M, SIGMA_P, EMULT_HOPRULE,                         & ! Output
         EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )             ! Output

          IF ( DO_PROFILE_LINEARIZATION ) THEN
            CALL LIDORT_LP_EMULT_MASTER  &
         ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,     & ! Input
           DO_PLANE_PARALLEL, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,   & ! Inputs
           TAYLOR_ORDER, N_USER_STREAMS, NBEAMS, NLAYERS,           & ! Input
           N_PARTLAYERS, PARTLAYERS_LAYERIDX,                       & ! Input
           USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! input
           DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,                 & ! input
           T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                & ! input
           ITRANS_USERM, LAYER_PIS_CUTOFF,  T_DELT_MUBAR,           & ! input
           SIGMA_M, SIGMA_P, EMULT_HOPRULE,                         & ! input
           EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,            & ! input
           LP_AVERAGE_SECANT, LP_INITIAL_TRANS,                     & ! input
           LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                       & ! input
           L_T_DELT_USERM,  L_T_UTDN_USERM, L_T_UTUP_USERM,         & ! input
           LP_EMULT_UP, LP_UT_EMULT_UP,                             & ! output
           LP_EMULT_DN, LP_UT_EMULT_DN )                              ! output
          ENDIF

        ENDIF
      ENDIF

!      write(*,*)'EMULT',EMULT_DN(3,15,3), UT_EMULT_DN(3,2,3)
!      write(*,*)'EMULT',EMULT_DN(4,15,3), UT_EMULT_DN(4,2,3)
!      write(*,*)'EMULT',EMULT_DN(5,15,3), UT_EMULT_DN(5,2,3)

!  Thermal setups
!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL argument to THERMAL_SETUP_PLUS
!  @@@@@@@@@@@ Robfix 19 March 2012. Do the Regular setups if no COLUMN linearization
!  Thermal Cutoff argument added 5/15/15

      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( DO_PROFILE_LINEARIZATION ) THEN
          CALL THERMAL_SETUP_PLUS                                   &
          ( DO_UPWELLING, DO_DNWELLING, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL, &
            DO_USER_STREAMS, DO_PROFILE_LINEARIZATION,              & ! input
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                 & ! input
            NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, TCUTOFF,    & ! input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_USER_STREAMS,     & ! input
            N_THERMAL_COEFFS, USER_STREAMS, THERMAL_BB_INPUT,       & ! input
            DELTAU_VERT, L_DELTAU_VERT, PARTAU_VERT, OMEGA_TOTAL,   & ! input
            L_OMEGA_TOTAL, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,& ! input
            L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,         & ! input
            deltau_power, xtau_power, thermcoeffs,                  & ! output
            L_deltau_power, L_xtau_power, L_thermcoeffs,            & ! output
            t_direct_up,    L_t_direct_up,                          & ! output
            t_direct_dn,    L_t_direct_dn,                          & ! output
            t_ut_direct_up, L_t_ut_direct_up,                       & ! output
            t_ut_direct_dn, L_t_ut_direct_dn )                        ! output
        ELSE
          CALL THERMAL_SETUP                                           &
          ( DO_UPWELLING, DO_DNWELLING, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL, &
            DO_USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,             & ! input
            NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,                          & ! input
            N_USER_STREAMS, N_THERMAL_COEFFS, TCUTOFF,                           & ! input
            USER_STREAMS, DELTAU_VERT, PARTAU_VERT, OMEGA_TOTAL,                 & ! input
            THERMAL_BB_INPUT, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,          & ! input
            DELTAU_POWER, XTAU_POWER, THERMCOEFFS,                               & ! output
            T_DIRECT_UP, T_UT_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_DN )             ! output
        ENDIF
      ENDIF

!  #####################
!  Correction operations
!  #####################

!  Single scatter corrections (pre-calculation)
!  ============================================

!     Must be done after MISCSETUPS and EMULT_MASTER, as we need
!      multipliers and transmittance factors for SUN and LOS paths.
!      Code added 6 May 2005. Replaces call in Master routine.
!      Version 3.1. Added call to the new outgoing sphericity correction

      IF ( DO_USER_STREAMS ) THEN

!  regular nadir-scattering SS correction.
!  Rob Fix 11/18/14. Revision of code
!     OUTPUT Removal and localization of SS_LEG variables, to reduce Stacksize.
!     OUTPUT Removal of TMS, SSFDEL, EXACTSCAT_UP, EXACTSCAT_DN.
!     LP_SSCORR is now a "Plus" routine (Intensity and Jacobians)

        IF ( DO_SSCORR_NADIR .AND. .NOT. DO_SS_EXTERNAL ) THEN

!  Intensity only or Intensity+Profile Linearization

          IF ( DO_PROFILE_LINEARIZATION ) THEN
            CALL LIDORT_LP_SSCORR_NADIR                                          &
              ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,             & ! Input
                DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, DO_REFRACTIVE_GEOMETRY, & ! Input
                SS_FLUX_MULTIPLIER, NLAYERS, NMOMENTS_INPUT, NMOMENTS, NBEAMS,   & ! Input
                N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                   & ! Input
                LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                              & ! Input
                BEAM_SZAS, SUNLAYER_COSINES, USER_STREAMS, USER_RELAZMS,         & ! Input
                N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,     & ! Input
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,    & ! Input
                UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, DO_PHFUNC_VARIATION,     & ! Input
                TRUNC_FACTOR, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,           & ! Input
                L_OMEGA_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                     & ! Input
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                    & ! Input
                T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                        & ! Input
                LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN,        & ! Input
                L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,                  & ! Input
                SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, INTENSITY_SS, PROFILEWF_SS )     ! Output
          ELSE
            CALL LIDORT_SSCORR_NADIR                                             &
              ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,             & ! Input
                DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, DO_REFRACTIVE_GEOMETRY, & ! Input
                SS_FLUX_MULTIPLIER, NLAYERS, NMOMENTS_INPUT, NBEAMS,             & ! Input
                N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                   & ! Input
                BEAM_SZAS, SUNLAYER_COSINES, USER_STREAMS, USER_RELAZMS,         & ! Input
                N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,     & ! Input
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,    & ! Input
                UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                          & ! Input
                TRUNC_FACTOR, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,           & ! Input
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                    & ! Input
                T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                        & ! Input
                SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, INTENSITY_SS )                   ! Output
          ENDIF
        ENDIF

!  New outgoing sphericity correction
!  ----------------------------------

        IF ( DO_SSCORR_OUTGOING .AND. .NOT. DO_SS_EXTERNAL ) THEN

          CALL LIDORT_LP_SSCORR_OUTGOING                                            &
              ( DO_UPWELLING, DO_DNWELLING, DO_PROFILE_LINEARIZATION,               & ! Input
                DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, DO_OBSERVATION_GEOMETRY,   & ! Input
                SS_FLUX_MULTIPLIER, NLAYERS, NFINELAYERS, NMOMENTS_INPUT, NMOMENTS, & ! Input
                NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,              & ! Input
                LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                                 & ! Input
                BEAM_SZAS_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,          & ! Input
                N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,        & ! Input
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,       & ! Input
                N_PARTLAYERS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,               & ! Input
                HEIGHT_GRID, EARTH_RADIUS, DELTAU_VERT, PARTAU_VERT,                & ! Input
                TRUNC_FACTOR, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,              & ! Input
                DO_PHFUNC_VARIATION, L_DELTAU_VERT,                                 & ! Input
                L_OMEGA_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                        & ! Input
                  UP_LOSTRANS,   UP_LOSTRANS_UT,   BOA_ATTN,                        & ! Output
                L_UP_LOSTRANS, L_UP_LOSTRANS_UT, L_BOA_ATTN,                        & ! Output
                INTENSITY_SS, PROFILEWF_SS,                                         & ! Output
                STATUS_SUB, MESSAGE, TRACE_1  )                                       ! Output

          IF ( STATUS_SUB .ne. LIDORT_SUCCESS ) THEN
            TRACE_2 = 'Error from LP_LIDORT_SSCORR_OUTGOING'
            TRACE_3 = ' ** Called in LIDORT_LPS_MASTER'
            STATUS_CALCULATION = LIDORT_SERIOUS
            LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            LIDORT_Out%Status%TS_MESSAGE = MESSAGE
            LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF 

        ENDIF

!  Exact (Direct bounce) surface term if single scatter is being done
!     renamed after introduction of BRDF option, 23 March 2010

        IF ( ( DO_SSCORR_NADIR .or. DO_SSCORR_OUTGOING ) &
              .AND. .NOT. DO_SS_EXTERNAL ) THEN

! Intensity term

          CALL LIDORT_DBCORRECTION                                       &
        ( DO_SSCORR_OUTGOING, DO_BRDF_SURFACE, DO_OBSERVATION_GEOMETRY,  & ! Input
          DO_REFRACTIVE_GEOMETRY, DO_UPWELLING,                          & ! Input
          DO_REFLECTED_DIRECTBEAM, SS_FLUX_MULTIPLIER, NLAYERS, NBEAMS,  & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                 & ! Input
          BEAM_SZAS, BEAM_SZAS_ADJUST, SZA_LOCAL_INPUT,                  & ! Input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP,                       & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
          LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC,                            & ! Input
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,                           & ! Input
          SLTERM_ISOTROPIC, SLTERM_USERANGLES,                           & ! Input
          SOLAR_BEAM_OPDEP, BOA_ATTN,                                    & ! Input
          T_DELT_USERM, UP_LOSTRANS, T_UTUP_USERM, UP_LOSTRANS_UT,       & ! Input
          INTENSITY_DB, ATTN_DB_SAVE, EXACTDB_SOURCE, DB_CUMSOURCE )       ! Output

!  Profile Linearization

          IF ( DO_PROFILE_LINEARIZATION ) THEN

            CALL LIDORT_LAP_DBCORRECTION                                    &
        ( DO_SSCORR_OUTGOING, DO_BRDF_SURFACE,                              & ! Input
          DO_REFRACTIVE_GEOMETRY, DO_UPWELLING,                             & ! Input
          DO_REFLECTED_DIRECTBEAM, DO_OBSERVATION_GEOMETRY,                 & ! Input
          SS_FLUX_MULTIPLIER, NLAYERS, NBEAMS,                              & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                    & ! Input
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                               & ! Input
          BEAM_SZAS, BEAM_SZAS_ADJUST, SZA_LOCAL_INPUT,                     & ! Input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP,                          & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,     & ! Input
          SOLAR_BEAM_OPDEP, DELTAU_SLANT, L_DELTAU_VERT, L_BOA_ATTN,        & ! Input
          T_DELT_USERM, T_UTUP_USERM, L_T_DELT_USERM, L_T_UTUP_USERM,       & ! Input
          UP_LOSTRANS, UP_LOSTRANS_UT, L_UP_LOSTRANS, L_UP_LOSTRANS_UT,     & ! Input
          LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, DB_CUMSOURCE,                 & ! Input
          PROFILEWF_DB )                                                      ! Output

          ENDIF

!  Surface linearization

          IF ( DO_SURFACE_LINEARIZATION ) THEN

            CALL LIDORT_LS_DBCORRECTION                                 &
        ( DO_SSCORR_OUTGOING, DO_UPWELLING, DO_BRDF_SURFACE,            & ! Input
          DO_REFLECTED_DIRECTBEAM, DO_OBSERVATION_GEOMETRY,             & ! Input
          NLAYERS, NBEAMS, N_SURFACE_WFS,                               & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                & ! Input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP, SS_FLUX_MULTIPLIER,  & ! Input
          LS_EXACTDB_BRDFUNC, ATTN_DB_SAVE,                             & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
          T_DELT_USERM, UP_LOSTRANS, T_UTUP_USERM, UP_LOSTRANS_UT,      & ! Input
          SURFACEWF_DB )                                                  ! Output

          ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@ Addition of SLEAVE WFs to Corrected Directbeam @@@@@@@@@@
!    R. Spurr, 22 August 2012
!
          IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then
            CALL LIDORT_LSSL_DBCORRECTION ( &
              DO_SSCORR_OUTGOING, DO_UPWELLING, DO_SL_ISOTROPIC,       & ! Input
              DO_OBSERVATION_GEOMETRY,                                 & ! Input
              NLAYERS, NBEAMS, N_SLEAVE_WFS, N_SURFACE_WFS,            & ! Input
              N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,           & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 & ! Input
              PARTLAYERS_LAYERIDX, UTAU_LEVEL_MASK_UP,                 & ! Input
              N_GEOMETRIES, UMOFF, DO_REFLECTED_DIRECTBEAM,            & ! Input
              SS_FLUX_MULTIPLIER,                                      & ! Input
              LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES,           & ! Input
              T_DELT_USERM, T_UTUP_USERM, UP_LOSTRANS, UP_LOSTRANS_UT, & ! Input
              SURFACEWF_DB )                                             ! Output
          ENDIF

!@@@@@@@@@@ END Addition of SLEAVE Corrected Directbeam @@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        ENDIF

!  End User streams clause

      ENDIF

!  ####################
!   MAIN FOURIER LOOP
!  ####################

!  Initialise Fourier loop
!  =======================

!  Set Number of Fourier terms (NMOMENTS = Maximum).
!    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

!  Local flags

      IF ( .NOT. DO_SOLAR_SOURCES  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_ISOTROPIC_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_MVOUT_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_SSFULL  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
        DO_PROCESS_FOURIER  = .false.
      ELSE
        DO_PROCESS_FOURIER  = .true.
      ENDIF

!  set Fourier number (1 or 2 for Rayleigh only, otherwise 2*Nstreams-1))
!    Local do no azimuth skips convergence. 4/9/20 Change number of Fouriers

      IF ( LOCAL_DO_NO_AZIMUTH ) THEN
        N_FOURIER_COMPONENTS = 0
      ELSE
        IF ( DO_RAYLEIGH_ONLY  ) THEN
          N_FOURIER_COMPONENTS = min(2,NMOMENTS)
        ELSE
          N_FOURIER_COMPONENTS = NMOMENTS
        ENDIF
      ENDIF

!  re-set no-azimuth flag

      DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

!  Fourier loop
!  ============

!  initialize

      LOCAL_ITERATION   = .TRUE.
      FOURIER_COMPONENT = -1
      TESTCONV          = 0

!  set up solar beam flags

      DO IBEAM = 1, NBEAMS
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO L = 0, MAXFOURIER
          DO_MULTIBEAM   ( IBEAM, L ) = .TRUE.
        ENDDO
      ENDDO

!  start loop
!  ----------

      DO WHILE ( LOCAL_ITERATION .AND. FOURIER_COMPONENT.LT.N_FOURIER_COMPONENTS )

!  Fourier counter

        FOURIER_COMPONENT = FOURIER_COMPONENT + 1

!  Local start of user-defined streams
!  Now set = 1. Fudging in earlier versions caused Havoc.
!  Here is the older version fudge
!        IF ( FOURIER_COMPONENT. EQ. 0 ) THEN
!          LOCAL_UM_START = 1
!        ELSE
!          LOCAL_UM_START = N_OUT_STREAMS - N_CONV_STREAMS + 1
!        ENDIF
!        LOCAL_UM_START = 1

!  azimuth cosine factor, using adjust geometries.

        IF ( FOURIER_COMPONENT .GT. 0 ) THEN
          DFC = DBLE(FOURIER_COMPONENT)
          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UA = 1, LOCAL_N_USERAZM
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  AZM_ARGUMENT = USER_RELAZMS_ADJUST(UM,IB,UA) * DFC
                  AZMFAC(UM,IB,UA) = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO IB = 1, NBEAMS
!              AZM_ARGUMENT = USER_RELAZMS_ADJUST(1,IB,1) * DFC
              AZM_ARGUMENT = USER_RELAZMS(IB) * DFC
              AZMFAC(LUM,IB,LUA) = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
            ENDDO
          ENDIF
        ENDIF

!  Main call to Lidort Fourier module
!  ----------------------------------

!  Only if flagged

        IF ( DO_PROCESS_FOURIER ) THEN

!mick debug
!write(*,*)'entering LIDORT_LPS_FOURIER_MASTER'

        !write(*,*)' ..calculating fourier component',FOURIER_COMPONENT

!  @@@@@@@@@@@ Robfix 13 January 2012, Added argument 'DO_MSMODE_THERMAL'
!  @@@@@@@@@@@ Robfix 10 October 2013 - Add TAYLOR_ORDER argument (Line 1)
!  Thermal Cutoff argument added 5/15/15

          CALL LIDORT_LPS_FOURIER_MASTER &
          ( FOURIER_COMPONENT, TAYLOR_ORDER, DO_OBSERVATION_GEOMETRY,                                                      & !Input
            DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_BRDF_SURFACE, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,       & !Input
            DO_THERMAL_TRANSONLY, DO_PLANE_PARALLEL, DO_USER_STREAMS, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,              & !Input
            DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_REFRACTIVE_GEOMETRY, DO_MSMODE_LIDORT, DO_MULTIBEAM, DO_MSMODE_THERMAL, & !Input
            NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, NBEAMS, NMOMENTS, NSTREAMS_2, NTOTAL,                        & !Input
            N_SUPDIAG, N_SUBDIAG, N_PARTLAYERS, N_THERMAL_COEFFS, N_DIRECTIONS, WHICH_DIRECTIONS, TCUTOFF,                 & !Input
            FLUX_FACTOR, BEAM_COSINES, SUNLAYER_COSINES, LOCAL_CSZA, LAYER_PIS_CUTOFF,                                     & !Input
            QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, USER_STREAMS, USER_SECANTS, SIGMA_P, SIGMA_M,                        & !Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                                                  & !Input
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                                & !Input
            DELTAU_VERT, PARTAU_VERT, OMEGA_MOMS, DELTAU_SLANT, SOLAR_BEAM_OPDEP, L_DELTAU_VERT, L_OMEGA_MOMS,             & !Input
            DELTAU_POWER, XTAU_POWER, THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_THERMCOEFFS,                            & !Input
            T_DIRECT_UP, T_UT_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_DN,                                                      & !Input
            L_T_DIRECT_UP, L_T_UT_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_DN,                                              & !Input
            INITIAL_TRANS, ITRANS_USERM, AVERAGE_SECANT, LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                              & !Input
            T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN, L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_DISORDS_UTDN,          & !Input
            T_DELT_MUBAR, T_UTDN_MUBAR, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, HELP_AQ, HELP_BQ,                                & !Input
            T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,                      & !Input
            EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN,        & !Input
            DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                        & !Input
            DO_ATMOS_LBBF, DO_SURFACE_LBBF, & ! NEW LINE 18 MARCH 2014
            LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,                                               & !Input
            N_SURFACE_WFS, LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,                                       & !Input
            SURFACE_BB_INPUT,EMISSIVITY,USER_EMISSIVITY,LS_EMISSIVITY,LS_USER_EMISSIVITY,                                  & !Input
            DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_USERANGLES, SLTERM_F_0, USER_SLTERM_F_0,         & !Input
            DO_SLEAVE_WFS, N_SLEAVE_WFS,LSSL_SLTERM_ISOTROPIC,LSSL_SLTERM_USERANGLES,LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0, & !Input
            DO_REFLECTED_DIRECTBEAM, DO_LAYER_SCATTERING, BVP_REGULAR_FLAG, LCONMASK, MCONMASK, BMAT_ROWMASK,              & !InOut
            DO_BVTEL_INITIAL, NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE, BTELMAT_ROWMASK,                             & !InOut
            INTENSITY_F, MEAN_INTENSITY, FLUX_INTEGRAL, DNMEAN_DIRECT, DNFLUX_DIRECT,                                      & !Output
            PROFILEWF_F, MINT_PROFILEWF, FLUX_PROFILEWF, SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,                      & !Output
            ABBWFS_JACOBIANS, ABBWFS_FLUXES, SBBWFS_JACOBIANS, SBBWFS_FLUXES, & ! NEW LINE 18 MARCH 2014
            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                                                          !Output

!mick debug
!write(*,*)'leaving LIDORT_LPS_FOURIER_MASTER'

!  error handling

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            TRACE_3 = ' ** Called in LIDORT_LPS_MASTER'
            STATUS_CALCULATION = LIDORT_SERIOUS
            LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            LIDORT_Out%Status%TS_MESSAGE = MESSAGE
            LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF

!  End fourier processing

        ENDIF

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!   -- only done for beams which are still not converged
!      This is controlled by flag DO_MULTIBEAM

!   -- new criterion, SS is added for Fourier = 0, as this means that
!      higher-order terms will be relatively smaller, which implies
!      faster convergence in some circumstances (generally not often).

        IBEAM_COUNT = 0
        DO IBEAM = 1, NBEAMS
          IF ( DO_MULTIBEAM ( IBEAM, FOURIER_COMPONENT ) ) THEN

!  Lattice calculation

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN

!  Convergence and radiance summation

              CALL LIDORT_CONVERGE                                       &
              ( DO_UPWELLING, DO_NO_AZIMUTH,                             & ! Input
                DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_SS_EXTERNAL,        & ! Input
                DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,          & ! Input
                DO_DOUBLE_CONVTEST, N_CONVTESTS, LIDORT_ACCURACY,        & ! Input
                N_USER_STREAMS, N_USER_LEVELS, N_USER_RELAZMS,           & ! Input
                NSTREAMS, IBEAM, FOURIER_COMPONENT,                      & ! Input
                UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,  & ! Input
                AZMFAC, INTENSITY_F, INTENSITY_SS, INTENSITY_DB,         & ! Input
                INTENSITY, FOURIER_SAVED,                                & ! output
                BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )              ! output

!  Fourier summation of Profile Jacobians

              IF ( DO_PROFILE_LINEARIZATION ) THEN
                CALL LIDORT_LP_CONVERGE                    &
                ( DO_UPWELLING, DO_SS_EXTERNAL, DO_SSFULL, & ! Input
                  DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,     & ! Input
                  DO_NO_AZIMUTH, AZMFAC, LOCAL_N_USERAZM,  & ! Input
                  N_USER_STREAMS, N_USER_LEVELS, NLAYERS,  & ! Input
                  IBEAM, FOURIER_COMPONENT,                & ! Input
                  UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS,   & ! Input
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER,      & ! Input
                  PROFILEWF_F, PROFILEWF_SS, PROFILEWF_DB, & ! Input
                  PROFILEWF )                                ! Output
              ENDIF

!  Fourier summation of Surface Jacobians

              IF ( DO_SURFACE_LINEARIZATION ) THEN
                CALL LIDORT_LS_CONVERGE                                      &
                ( DO_UPWELLING, DO_BRDF_SURFACE, DO_SS_EXTERNAL,             & ! Input
                  DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,            & ! Input
                  N_USER_STREAMS, N_USER_LEVELS, N_SURFACE_WFS, N_SLEAVE_WFS,& ! Input
                  LOCAL_N_USERAZM, AZMFAC, IBEAM, FOURIER_COMPONENT,         & ! Input
                  UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS,                     & ! Input
                  SURFACEWF_F,  SURFACEWF_DB,                                & ! Input
                  SURFACEWF )                                                  ! Output
              ENDIF

!  Observational geometry calculation

            ELSE

!  Convergence and radiance summation

              CALL LIDORT_CONVERGE_OBSGEO                                         &
               ( DO_UPWELLING, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_SS_EXTERNAL,  & ! Input
                 DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                  & ! Input
                 DO_DOUBLE_CONVTEST, N_CONVTESTS, LIDORT_ACCURACY,                & ! Input
                 N_USER_LEVELS, NSTREAMS, IBEAM, FOURIER_COMPONENT,               & ! Input
                 N_DIRECTIONS, WHICH_DIRECTIONS, AZMFAC,                          & ! Input
                 INTENSITY_F, INTENSITY_SS, INTENSITY_DB,                         & ! Input
                 INTENSITY, FOURIER_SAVED,                                        & ! output
                 BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )                      ! output

!  Fourier summation of Profile Jacobians

              IF ( DO_PROFILE_LINEARIZATION ) THEN
                CALL LIDORT_LP_CONVERGE_OBSGEO             &
                ( DO_UPWELLING, DO_SS_EXTERNAL, DO_SSFULL, & ! Input
                  DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,     & ! Input
                  DO_NO_AZIMUTH, AZMFAC,                   & ! Input
                  N_USER_LEVELS, NLAYERS,                  & ! Input
                  IBEAM, FOURIER_COMPONENT,                & ! Input
                  N_DIRECTIONS, WHICH_DIRECTIONS,          & ! Input
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER,      & ! Input
                  PROFILEWF_F, PROFILEWF_SS, PROFILEWF_DB, & ! Input
                  PROFILEWF )                                ! Output
              ENDIF

!  Fourier summation of Surface Jacobians

              IF ( DO_SURFACE_LINEARIZATION ) THEN
                CALL LIDORT_LS_CONVERGE_OBSGEO                               &
                ( DO_UPWELLING, DO_BRDF_SURFACE, DO_SS_EXTERNAL,             & ! Input
                  DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,            & ! Input
                  N_USER_LEVELS, N_SURFACE_WFS, N_SLEAVE_WFS,                & ! Input
                  AZMFAC, IBEAM, FOURIER_COMPONENT,                          & ! Input
                  N_DIRECTIONS, WHICH_DIRECTIONS,                            & ! Input
                  SURFACEWF_F,  SURFACEWF_DB,                                & ! Input
                  SURFACEWF )                                                  ! Output
              ENDIF

            ENDIF

!  Check number of beams already converged

            IF ( BEAM_ITERATION(IBEAM) ) THEN
              IBEAM_COUNT = IBEAM_COUNT + 1
            ELSE
              DO L = FOURIER_COMPONENT+1,MAXFOURIER
                DO_MULTIBEAM (IBEAM,L) = .FALSE.
              ENDDO
            ENDIF

!  end beam count loop

          ENDIF
        END DO

!  If all beams have converged, stop iteration

        IF ( IBEAM_COUNT .EQ. 0 ) LOCAL_ITERATION = .FALSE.

!  Fourier output
!  --------------

!  Open file if Fourier = 0
!  Write Standard Fourier output
!  Close file if iteration has finished

!  New comment:
!    If the SS correction is set, Fourier=0 will include SS field

!        IF ( DO_WRITE_FOURIER ) THEN
!          FUNIT = LIDORT_FUNIT
!          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
!            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
!          ENDIF
!          CALL LIDORT_WRITEFOURIER ( FUNIT, FOURIER_COMPONENT )
!          IF ( DO_PROFILE_LINEARIZATION ) THEN
!            CALL LIDORT_LP_WRITEFOURIER ( DO_INCLUDE_SURFACE, FUNIT, FOURIER_COMPONENT )
!          ENDIF
!          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
!        ENDIF

!  end iteration loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  OUTPUT COPYING TASK, IN/OUT variables
!  =====================================

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

!  Modified control inputs

      LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT = NMOMENTS_INPUT

!  Modified beam inputs

      LIDORT_ModIn%MSunRays%TS_NBEAMS              = NBEAMS
      LIDORT_ModIn%MSunRays%TS_BEAM_SZAS(1:NBEAMS) = BEAM_SZAS(1:NBEAMS)

!  Modified user value inputs

      LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS      = N_USER_RELAZMS
      LIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS) = &
        USER_RELAZMS(1:N_USER_RELAZMS)

      LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS      = N_USER_STREAMS
      !LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT   = USER_ANGLES

      LIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS) = &
        USER_LEVELS(1:N_USER_LEVELS)

      LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

!  New Code, 26 October 2012
      LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = DO_OBSERVATION_GEOMETRY
      LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS      = N_USER_OBSGEOMS

      LIDORT_ModIn%MUserVal%TS_USER_OBSGEOM_INPUT(1:N_USER_OBSGEOMS,:) = &
        USER_OBSGEOMS(1:N_USER_OBSGEOMS,:)

!  Modified Chapman function inputs

      !LIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS = CHAPMAN_FACTORS
      LIDORT_ModIn%MChapman%TS_EARTH_RADIUS    = EARTH_RADIUS

!  Modified optical inputs

      LIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS) = &
        OMEGA_TOTAL_INPUT(1:NLAYERS)

!  OUTPUT COPYING TASK, pure OUT variables
!  =======================================

!mick debug
!write(*,*) 'copying output to type structures'

!mick fix 6/29/11 - pass output array entries
!   Rob : Direct-beam contributions output separately, 26 May 11, 24 August 2011

!  Radiance and mean

      LIDORT_Out%Main%TS_INTENSITY(1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
        INTENSITY(1:N_USER_LEVELS,1:N_GEOMETRIES,:)
      LIDORT_Out%Main%TS_MEAN_INTENSITY(1:N_USER_LEVELS,1:NBEAMS,:)  = &
        MEAN_INTENSITY(1:N_USER_LEVELS,1:NBEAMS,:)
      LIDORT_Out%Main%TS_FLUX_INTEGRAL(1:N_USER_LEVELS,1:NBEAMS,:)   = &
        FLUX_INTEGRAL(1:N_USER_LEVELS,1:NBEAMS,:)

      LIDORT_Out%Main%TS_DNMEAN_DIRECT(1:N_USER_LEVELS,1:NBEAMS) = &
        DNMEAN_DIRECT(1:N_USER_LEVELS,1:NBEAMS)  ! Added 8/24/11
      LIDORT_Out%Main%TS_DNFLUX_DIRECT(1:N_USER_LEVELS,1:NBEAMS) = &
        DNFLUX_DIRECT(1:N_USER_LEVELS,1:NBEAMS)  ! Added 8/24/11

!  new 12 March 2012
!   IF SS results already available, no need to copy them !

      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
           INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,:)
         LIDORT_Sup%SS%TS_INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES) = &
           INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)
      ENDIF

!  Bookkeeping

      LIDORT_Out%Main%TS_FOURIER_SAVED(1:NBEAMS) = FOURIER_SAVED(1:NBEAMS)
      LIDORT_Out%Main%TS_N_GEOMETRIES            = N_GEOMETRIES

!  Profile Jacobians

      LIDORT_LinOut%Atmos%TS_PROFILEWF&
        (1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
          PROFILEWF(1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)
      LIDORT_LinOut%Atmos%TS_MINT_PROFILEWF&
        (1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:NBEAMS,:)       = &
          MINT_PROFILEWF(1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:NBEAMS,:)
      LIDORT_LinOut%Atmos%TS_FLUX_PROFILEWF&
        (1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:NBEAMS,:)       = &
          FLUX_PROFILEWF(1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:NBEAMS,:)

!  Surface Jacobians

      LIDORT_LinOut%Surf%TS_SURFACEWF&
        (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
          SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)
      LIDORT_LinOut%Surf%TS_MINT_SURFACEWF&
        (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)       = &
          MINT_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)
      LIDORT_LinOut%Surf%TS_FLUX_SURFACEWF&
        (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)       = &
          FLUX_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)

!  new 12 March 2012
!   IF SS results already available, no need to copy them !

      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         !SS atmosphere profile weighting functions

         LIDORT_LinSup%SS%Atmos%TS_PROFILEWF_SS&
           (1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
           PROFILEWF_SS(1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)
         LIDORT_LinSup%SS%Atmos%TS_PROFILEWF_DB&
           (1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:N_GEOMETRIES) = &
           PROFILEWF_DB(1:MAX_ATMOSWFS,1:NLAYERS,1:N_USER_LEVELS,1:N_GEOMETRIES)

         !SS surface weighting functions

         !LIDORT_LinSup%SS%Surf%TS_SURFACEWF_SS&
         !  (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
         !  SURFACEWF_SS(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)
         LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB&
           (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)   = &
           SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)
      ENDIF

!  BLACKBODY JACOBIANS, New 18 March 2014. 

      LIDORT_LinOut%Atmos%TS_ABBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_STREAMS,0:NLAYERS,1:2) = &
          ABBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_STREAMS,0:NLAYERS,1:2)
      LIDORT_LinOut%Atmos%TS_ABBWFS_FLUXES(1:N_USER_LEVELS,1:2,0:NLAYERS,1:2) = &
          ABBWFS_FLUXES(1:N_USER_LEVELS,1:2,0:NLAYERS,1:2)

      LIDORT_LinOut%Surf%TS_SBBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_STREAMS,1:2) = &
          SBBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_STREAMS,1:2)
      LIDORT_LinOut%Surf%TS_SBBWFS_FLUXES(1:N_USER_LEVELS,1:2,1:2) = &
          SBBWFS_FLUXES(1:N_USER_LEVELS,1:2,1:2)

!  Solar Beam Transmittance to BOA
!  rob fix 11/27/2014, for diagnostic use only

      LIDORT_Out%Main%TS_SOLARBEAM_BOATRANS (1:NBEAMS) = SOLARBEAM_BOATRANS (1:NBEAMS)

!  Exception handling

      LIDORT_Out%Status%TS_STATUS_INPUTCHECK  = STATUS_INPUTCHECK
      LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION

      LIDORT_Out%Status%TS_NCHECKMESSAGES                  = &
        NCHECKMESSAGES
      LIDORT_Out%Status%TS_CHECKMESSAGES(0:NCHECKMESSAGES) = &
        CHECKMESSAGES(0:NCHECKMESSAGES)
      LIDORT_Out%Status%TS_ACTIONS(0:NCHECKMESSAGES)       = &
        ACTIONS(0:NCHECKMESSAGES)

      LIDORT_Out%Status%TS_MESSAGE = MESSAGE
      LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
      LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
      LIDORT_Out%Status%TS_TRACE_3 = TRACE_3

!  END OUTPUT COPYING TASK
!  =======================

!mick debug
!write(*,*) 'leaving LIDORT_LPS_master'

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE LIDORT_DEBUG_INPUT_MASTER()

      CALL LIDORT_WRITE_STD_INPUT ( &
        DO_FULLRAD_MODE,DO_SSCORR_TRUNCATION,DO_SS_EXTERNAL,&
        DO_SSFULL,DO_THERMAL_EMISSION,DO_SURFACE_EMISSION,&
        DO_PLANE_PARALLEL,DO_BRDF_SURFACE,DO_UPWELLING,DO_DNWELLING,&
        DO_SURFACE_LEAVING,DO_SL_ISOTROPIC,&
        NSTREAMS,NLAYERS,NFINELAYERS,N_THERMAL_COEFFS,TCUTOFF,&
        LIDORT_ACCURACY,FLUX_FACTOR,N_USER_LEVELS,&
        HEIGHT_GRID,PRESSURE_GRID,TEMPERATURE_GRID,&
        FINEGRID,RFINDEX_PARAMETER,&
        DELTAU_VERT_INPUT,PHASMOMS_TOTAL_INPUT,&
        THERMAL_BB_INPUT,LAMBERTIAN_ALBEDO,SURFACE_BB_INPUT,&
        DO_SSCORR_NADIR,DO_SSCORR_OUTGOING,DO_DOUBLE_CONVTEST,&
        DO_SOLAR_SOURCES,DO_REFRACTIVE_GEOMETRY,DO_CHAPMAN_FUNCTION,&
        DO_RAYLEIGH_ONLY,DO_ISOTROPIC_ONLY,DO_NO_AZIMUTH,DO_ALL_FOURIER,&
        DO_DELTAM_SCALING,DO_SOLUTION_SAVING,DO_BVP_TELESCOPING,&
        DO_USER_STREAMS,DO_ADDITIONAL_MVOUT,&
        DO_MVOUT_ONLY,DO_THERMAL_TRANSONLY,DO_OBSERVATION_GEOMETRY,&
        NMOMENTS_INPUT,NBEAMS,BEAM_SZAS,N_USER_RELAZMS,USER_RELAZMS,&
        N_USER_STREAMS,USER_ANGLES_INPUT,USER_LEVELS,&
        GEOMETRY_SPECHEIGHT,N_USER_OBSGEOMS,USER_OBSGEOMS,&
        EARTH_RADIUS,OMEGA_TOTAL_INPUT )

      IF (DO_BRDF_SURFACE) THEN
        CALL LIDORT_WRITE_SUP_BRDF_INPUT ( &
          NSTREAMS,NBEAMS,N_USER_STREAMS,N_USER_RELAZMS,&
          EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,&
          EMISSIVITY,USER_EMISSIVITY)
      END IF

      IF (DO_SS_EXTERNAL) THEN
        CALL LIDORT_WRITE_SUP_SS_INPUT ( &
          N_USER_LEVELS,&
          INTENSITY_SS,INTENSITY_DB)
      END IF

      IF (DO_SURFACE_LEAVING) THEN
        CALL LIDORT_WRITE_SUP_SLEAVE_INPUT ( &
          NSTREAMS,NBEAMS,N_USER_STREAMS,N_USER_RELAZMS,&
          SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)
      END IF

      END SUBROUTINE LIDORT_DEBUG_INPUT_MASTER

      SUBROUTINE LIDORT_DEBUG_LIN_INPUT_MASTER()

      CALL LIDORT_WRITE_LIN_INPUT ( &
        NLAYERS,NMOMENTS_INPUT,&
        DO_COLUMN_LINEARIZATION,DO_PROFILE_LINEARIZATION,&
        DO_SURFACE_LINEARIZATION,DO_SLEAVE_WFS,&
        LAYER_VARY_FLAG,LAYER_VARY_NUMBER,&
        N_TOTALCOLUMN_WFS,N_SURFACE_WFS,N_SLEAVE_WFS,&
        L_DELTAU_VERT_INPUT,L_OMEGA_TOTAL_INPUT,&
        L_PHASMOMS_TOTAL_INPUT)

      IF (DO_BRDF_SURFACE .AND. DO_SURFACE_LINEARIZATION) THEN
        CALL LIDORT_WRITE_LIN_SUP_BRDF_INPUT ( &
          NSTREAMS,NBEAMS,N_USER_STREAMS,N_USER_RELAZMS,N_SURFACE_WFS,&
          LS_EXACTDB_BRDFUNC,LS_BRDF_F_0,LS_BRDF_F,LS_USER_BRDF_F_0,LS_USER_BRDF_F,&
          LS_EMISSIVITY,LS_USER_EMISSIVITY)
      END IF

      IF (DO_SS_EXTERNAL) THEN
        CALL LIDORT_WRITE_LPS_SUP_SS_INPUT ( &
          DO_PROFILE_LINEARIZATION,DO_SURFACE_LINEARIZATION,&
          NLAYERS,N_USER_LEVELS,N_TOTALPROFILE_WFS,N_SURFACE_WFS,&
          PROFILEWF_SS,PROFILEWF_DB,SURFACEWF_DB)
      END IF

      IF (DO_SURFACE_LEAVING .AND. DO_SURFACE_LINEARIZATION &
          .AND. DO_SLEAVE_WFS) THEN
        CALL LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT ( &
          NSTREAMS,NBEAMS,N_USER_STREAMS,N_USER_RELAZMS,N_SLEAVE_WFS,&
          LSSL_SLTERM_ISOTROPIC,LSSL_SLTERM_USERANGLES,&
          LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0)
      END IF

      END SUBROUTINE LIDORT_DEBUG_LIN_INPUT_MASTER

   END SUBROUTINE LIDORT_LPS_master

!  @@@@@@@@@@@ Robfix 13 January 2012 - Add argument 'DO_MSMODE_THERMAL' 
!  @@@@@@@@@@@ Robfix 10 October 2013 - Add TAYLOR_ORDER argument (Line 1)
!  Thermal Cutoff argument added 5/15/15

SUBROUTINE LIDORT_LPS_FOURIER_MASTER&
      ( FOURIER_COMPONENT, TAYLOR_ORDER, DO_OBSERVATION_GEOMETRY,                                                      & !Input
        DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_BRDF_SURFACE, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,       & !Input
        DO_THERMAL_TRANSONLY, DO_PLANE_PARALLEL, DO_USER_STREAMS, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,              & !Input
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_REFRACTIVE_GEOMETRY, DO_MSMODE_LIDORT, DO_MULTIBEAM, DO_MSMODE_THERMAL, & !Input
        NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, NBEAMS, NMOMENTS, NSTREAMS_2, NTOTAL,                        & !Input
        N_SUPDIAG, N_SUBDIAG, N_PARTLAYERS, N_THERMAL_COEFFS, N_DIRECTIONS, WHICH_DIRECTIONS, TCUTOFF,                 & !Input
        FLUX_FACTOR, BEAM_COSINES, SUNLAYER_COSINES, LOCAL_CSZA, LAYER_PIS_CUTOFF,                                     & !Input
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, USER_STREAMS, USER_SECANTS, SIGMA_P, SIGMA_M,                        & !Input
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                                                  & !Input
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                                & !Input
        DELTAU_VERT, PARTAU_VERT, OMEGA_MOMS, DELTAU_SLANT, SOLAR_BEAM_OPDEP, L_DELTAU_VERT, L_OMEGA_MOMS,             & !Input
        DELTAU_POWER, XTAU_POWER, THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_THERMCOEFFS,                            & !Input
        T_DIRECT_UP, T_UT_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_DN,                                                      & !Input
        L_T_DIRECT_UP, L_T_UT_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_DN,                                              & !Input
        INITIAL_TRANS, ITRANS_USERM, AVERAGE_SECANT, LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                              & !Input
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN, L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_DISORDS_UTDN,          & !Input
        T_DELT_MUBAR, T_UTDN_MUBAR, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, HELP_AQ, HELP_BQ,                                & !Input
        T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,                      & !Input
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN,        & !Input
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                        & !Input
        DO_ATMOS_LBBF, DO_SURFACE_LBBF, & ! NEW LINE 18 MARCH 2014
        LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,                                               & !Input
        N_SURFACE_WFS, LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,                                       & !Input
        SURFACE_BB_INPUT,EMISSIVITY,USER_EMISSIVITY,LS_EMISSIVITY,LS_USER_EMISSIVITY,                                  & !Input
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_USERANGLES, SLTERM_F_0, USER_SLTERM_F_0,         & !Input
        DO_SLEAVE_WFS, N_SLEAVE_WFS,LSSL_SLTERM_ISOTROPIC,LSSL_SLTERM_USERANGLES,LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0, & !Input
        DO_REFLECTED_DIRECTBEAM, DO_LAYER_SCATTERING, BVP_REGULAR_FLAG, LCONMASK, MCONMASK, BMAT_ROWMASK,              & !InOut
        DO_BVTEL_INITIAL, NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE, BTELMAT_ROWMASK,                             & !InOut
        INTENSITY_F, MEAN_INTENSITY, FLUX_INTEGRAL, DNMEAN_DIRECT, DNFLUX_DIRECT,                                      & !Output
        PROFILEWF_F, MINT_PROFILEWF, FLUX_PROFILEWF, SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,                      & !Output
         ABBWFS_JACOBIANS, ABBWFS_FLUXES, SBBWFS_JACOBIANS, SBBWFS_FLUXES, & ! NEW LINE 18 MARCH 2014
       STATUS, MESSAGE, TRACE_1, TRACE_2 )                                                                              !Output

!  Complete Fourier component calculation for the Profile-Jacobian Code.

!  Parameter types

      USE LIDORT_pars

!  Implicit none

      implicit none

!  Input
!  -----

!  Observational Geometry flag

      LOGICAL,intent(in)   :: DO_OBSERVATION_GEOMETRY

!  Basic top-level control

      LOGICAL, intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL, intent(in)  :: DO_THERMAL_EMISSION

!  Surface emission flag

      LOGICAL, intent(in)  :: DO_SURFACE_EMISSION

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  :: DO_THERMAL_TRANSONLY

!  Beam particular solution pseudo-spherical options

      LOGICAL, intent(in)  :: DO_PLANE_PARALLEL

!  Stream angle flag

      LOGICAL, intent(in)  :: DO_USER_STREAMS

!  Beam particular solution pseudo-spherical options

      LOGICAL, intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Surface control (New, 23 March 2010)

      LOGICAL, intent(in)  :: DO_BRDF_SURFACE

!  Directional control

      LOGICAL, intent(in)  :: DO_UPWELLING
      LOGICAL, intent(in)  :: DO_DNWELLING

! @@@@@@@@@@@@@ Robfix, 13 January 2012.
!               Add MSMODE_THERMAL flag to calling statement
      logical, intent(in)  :: DO_MSMODE_THERMAL
! @@@@@@@@@@@@@ End Robfix, 13 January 2012.

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7.
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER, intent(in)  :: NSTREAMS

!  Number of computational layers

      INTEGER, intent(in)  :: NLAYERS

!  Number of solar beams to be processed

      INTEGER, intent(in)  :: NBEAMS

!  Local input solar zenith angles Cosines

      REAL(fpk), intent(in)  :: SUNLAYER_COSINES(MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: BEAM_COSINES(MAXBEAMS)

!  User-defined zenith angle input 

      INTEGER, intent(in)  :: N_USER_STREAMS

!  User-defined vertical level output
!    New system. IF input = 0.1, this means in layer 1, but only 0.1 down

      INTEGER, intent(in)  :: N_USER_LEVELS

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, intent(in)  :: LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER, intent(in)  :: LAYER_VARY_NUMBER(MAXLAYERS)

!  Total number of Surface Jacobians

      INTEGER, intent(in)  :: N_SURFACE_WFS

!  Control for  Blackbody Jacobians, New 18 March 2014

      LOGICAL  , intent(in)  :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  Number of thermal coefficients

      INTEGER, intent(in)  :: N_THERMAL_COEFFS

!  Number of directions (1 or 2) and directional array

      INTEGER, intent(in)  :: N_DIRECTIONS
      INTEGER, intent(in)  :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Thermal Cutoff parameter. Should be set to 1.0e-08.
!  This is actually a minimum value of layer optical thickness
!     Added, 5/14/15 for use with Thermal Solutions

      REAL(fpk),intent(in) :: TCUTOFF

!  Control linearization

      LOGICAL, intent(in)  :: DO_PROFILE_LINEARIZATION
      LOGICAL, intent(in)  :: DO_SURFACE_LINEARIZATION

!  Performance enhancement

      LOGICAL, intent(inout)  :: DO_SOLUTION_SAVING
      LOGICAL, intent(inout)  :: DO_BVP_TELESCOPING

!  Mean value control

      LOGICAL, intent(in)  :: DO_ADDITIONAL_MVOUT
      LOGICAL, intent(in)  :: DO_MVOUT_ONLY

!  Bookkeeping arguments
!  ---------------------

!  Mode of operation

      LOGICAL, intent(in)  :: DO_MSMODE_LIDORT

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER, intent(in)  :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER, intent(in)  :: NSTREAMS_2
      INTEGER, intent(in)  :: NTOTAL
      INTEGER, intent(in)  :: N_SUBDIAG, N_SUPDIAG

!  Number of partial layers

      INTEGER  , intent(in) :: N_PARTLAYERS

!  Quadrature weights and abscissae, and product

      REAL(fpk), intent(in)   :: QUAD_STREAMS (MAXSTREAMS)
      REAL(fpk), intent(in)   :: QUAD_WEIGHTS (MAXSTREAMS)
      REAL(fpk), intent(in)   :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      REAL(fpk), intent(in)   :: USER_STREAMS (MAX_USER_STREAMS)
      REAL(fpk), intent(in)   :: USER_SECANTS (MAX_USER_STREAMS)

!  Offgrid output optical depth masks and indices

      LOGICAL, intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL, intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL, intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Solar beam flags (always internal)

      LOGICAL, intent(in)  :: DO_MULTIBEAM (MAXBEAMS,0:MAXFOURIER)

!  Input Fourier component number

      INTEGER, intent(in)  :: FOURIER_COMPONENT

!  Absoute flux factor

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Surface Blackbody

      REAL(fpk), intent(in)  :: SURFACE_BB_INPUT

!  Lambertian Surface control

      REAL(fpk), intent(in)  :: LAMBERTIAN_ALBEDO

!  Fourier components of BRDF, in the following order
!    ( New code, 23 March 2010 )

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: BRDF_F_0      ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: BRDF_F        ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

      REAL(fpk), intent(in)  :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Fourier components of BRDF, in the following order

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in) :: LS_BRDF_F        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

      REAL(fpk), intent(in) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Emissivity

      REAL(fpk), intent(in)  :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Linearized Emissivity

      REAL(fpk), intent(in)  :: LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTREAMS )
      REAL(fpk), intent(in)  :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  New Surface-Leaving stuff 17 May 2012

      LOGICAL, intent(in) ::    DO_SURFACE_LEAVING
      LOGICAL, intent(in) ::    DO_SL_ISOTROPIC

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(in) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Exact Surface-Leaving term

      REAL(fpk), intent(in) ::  SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), intent(in) ::  SLTERM_F_0      ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in) ::  USER_SLTERM_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
      LOGICAL, intent(in)          ::   DO_SLEAVE_WFS
      INTEGER, intent(in)          ::   N_SLEAVE_WFS

      REAL(fpk), intent(in) :: LSSL_SLTERM_ISOTROPIC   ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(fpk), intent(in) :: LSSL_SLTERM_USERANGLES  ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      REAL(fpk), intent(in) :: LSSL_SLTERM_F_0         ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in) :: LSSL_USER_SLTERM_F_0    ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  Arrays required at the Top level
!  ================================

!                    Intent(In) to the Fourier routine

!  Input optical depths after delta-M scaling

      REAL(fpk), intent(in) :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk), intent(in) :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk), intent(in) :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )
      REAL(fpk), intent(in) :: DELTAU_SLANT   ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), intent(in) :: L_OMEGA_MOMS  ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in) :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(in) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in) :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Solar beam attenuation

      REAL(fpk), intent(in)   :: SOLAR_BEAM_OPDEP ( MAXBEAMS )

!  Reflectance flags

      LOGICAL, intent(InOut)  :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local flags for the solution saving option

      LOGICAL, intent(InOut)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Local flags,  BVP telescoping enhancement

      LOGICAL, intent(InOut)  :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Masking for regular case. Required again for linearization

      INTEGER, intent(inout)  :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER, intent(inout)  :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Telescoping initial flag (modified argument), Layer bookkeeping
!  Number of telescoped layers, active layers,  Size of BVP matrix 

      LOGICAL, intent(inout)  :: DO_BVTEL_INITIAL
      INTEGER, intent(inout)  :: NLAYERS_TEL
      INTEGER, intent(inout)  :: ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, intent(inout)  :: N_BVTELMATRIX_SIZE

!  Set up for band matrix compression

      INTEGER, intent(inout)  :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)
      INTEGER, intent(inout)  :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Transmittance Setups
!  --------------------

!                    Intent(In) to the Fourier routine

!  Discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!                 Intent(In) to Fourier Routine

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Coefficient functions for user-defined angles. Other constants
!Rob fix 5/6/13 - added these 3 arrays to input list

      REAL(fpk), intent(in)  :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Partial layer multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Linearized Transmittance Setups
!  -------------------------------

!                    Intent(In) to the Fourier routine

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized transmittances, solar beam

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Beam multipliers
!  ---------------------------

!                    Intent(In) to the Fourier routine

!  Linearized whole layer multipliers

      REAL(fpk), intent(in)  :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(fpk), intent(in)  :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized help arrays for Green's function

      REAL(fpk), intent(in)  :: HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Thermal Setup outputs (Input here)
!  ---------------------

!  Optical depth powers

      REAL(fpk), intent(in)  :: deltau_power (maxlayers,     max_thermal_coeffs)
      REAL(fpk), intent(in)  :: xtau_power (max_partlayers,max_thermal_coeffs)

!  Thermal coefficients

      REAL(fpk), intent(in)  :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk), intent(in)  :: t_direct_up ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: t_direct_dn ( MAX_USER_STREAMS, MAXLAYERS )

      REAL(fpk), intent(in)  :: t_ut_direct_up (MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: t_ut_direct_dn (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized optical depth powers

      REAL(fpk), intent(in)  :: L_deltau_power (maxlayers,     max_thermal_coeffs,max_atmoswfs)
      REAL(fpk), intent(in)  :: L_xtau_power (max_partlayers,max_thermal_coeffs,max_atmoswfs)

!  Linearized Thermal coefficients

      REAL(fpk), intent(in)  :: L_THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS,max_atmoswfs)

!  Linearized Tranmsittance solutions

      REAL(fpk), intent(in)  :: L_t_direct_up ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_t_direct_dn ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      REAL(fpk), intent(in)  :: L_t_ut_direct_up ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_t_ut_direct_dn ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Outputs
!  =======

!  Fourier component solutions

      REAL(fpk), intent(out) :: INTENSITY_F (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!mick fix 6/29/11 - changed these two from "out" to "inout"
!Rob  fix 8/24/11 - make sure Direct DN outputs are "inout"

!  Results for mean-value output

      REAL(fpk), intent(inout) :: MEAN_INTENSITY (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
      REAL(fpk), intent(inout) :: FLUX_INTEGRAL  (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)

!  Direct-beam contributions output separately, 26 May 11, 24 August 2011

      REAL(fpk), intent(inout)  :: DNMEAN_DIRECT (MAX_USER_LEVELS, MAXBEAMS)
      REAL(fpk), intent(inout)  :: DNFLUX_DIRECT (MAX_USER_LEVELS, MAXBEAMS)

!  Fourier-component Profile weighting functions at user angles

      REAL(fpk), intent(out) :: PROFILEWF_F &
         ( MAX_ATMOSWFS,     MAXLAYERS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS,  MAX_DIRECTIONS)

!  Mean-intensity and Flux weighting functions

      REAL(fpk), intent(inout) :: MINT_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), intent(inout) :: FLUX_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

!  Surface weighting functions at user angles

      REAL(fpk), intent(out) :: SURFACEWF_F &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Flux and mean-intensity Surface weighting functions

      REAL(fpk), intent(inout) :: MINT_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
      REAL(fpk), intent(inout) :: FLUX_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)

!  BLACKBODY Postprocessed Jacobians and Fluxes, New, 18 March 2014

      REAL(fpk), INTENT(INOUT) :: ABBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS)
      REAL(fpk), INTENT(INOUT) :: ABBWFS_FLUXES    ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS)

      REAL(fpk), INTENT(INOUT) :: SBBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS)
      REAL(fpk), INTENT(INOUT) :: SBBWFS_FLUXES    ( MAX_USER_LEVELS, 2, MAX_DIRECTIONS)

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER, intent(out)         :: STATUS
      CHARACTER*(*), intent(out)   :: MESSAGE, TRACE_1, TRACE_2

!  Local Arrays for argument passing
!  =================================

!  Atmospheric attenuation

      REAL(fpk)    :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solutions

      REAL(fpk)    :: DIRECT_BEAM      ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk)    :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Multiplier arrays
!  -----------------

!  Coefficient functions for user-defined angles

      REAL(fpk)    :: ZETA_P(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: ZETA_M(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk)    :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk)    :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk)    :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Legendre Setups
!  ---------------

!  At quadrature angles

      REAL(fpk)    :: LEG_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: LEG_M(MAXSTREAMS,0:MAXMOMENTS)

!  At beam angles. LEG0_M holds stored quantities.

      REAL(fpk)    :: LEG0_P(0:MAXMOMENTS)
      REAL(fpk)    :: LEG0_M(0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Legendre polynomial products

      REAL(fpk)    :: PLMI_PLMJ_P(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: PLMI_PLMJ_M(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: PLMI_X0_P(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)
      REAL(fpk)    :: PLMI_X0_M(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Polynomial-weight products

      REAL(fpk)    :: WT_LEGP(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: WT_LEGM(MAXSTREAMS,0:MAXMOMENTS)

!  Legendre functions on User defined polar angles

      REAL(fpk)    :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: U_LEG_M(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  Local matrices for eigenvalue computation

      REAL(fpk)    :: SAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: DAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: EIGENMAT_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: EIGENVEC_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: DIFVEC_SAVE  (MAXSTREAMS,MAXSTREAMS)

!  (Positive) Eigenvalues

      REAL(fpk)    :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk)    :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(fpk)    :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Saved help variables

      REAL(fpk)    :: U_HELP_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: U_HELP_M(MAXSTREAMS,0:MAXMOMENTS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk)    :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Help arrays for reflected solutions

      REAL(fpk)    :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: H_XNEG(MAXSTREAMS,MAXSTREAMS)

!  Boundary Value Problem
!  ----------------------

!  Single Matrix, Band-matrices

      REAL(fpk)    :: SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk)    :: BANDMAT2   (MAXBANDTOTAL,MAXTOTAL)
      REAL(fpk)    :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER      :: IPIVOT     (MAXTOTAL)
      INTEGER      :: SIPIVOT    (MAXSTREAMS_2)
      INTEGER      :: IPIVOTTEL  (MAXTOTAL)

!  Particular integrals
!  --------------------

!  General beam solutions at the boundaries

      REAL(fpk)    :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)    :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Help array for reflected solutions

      REAL(fpk)    :: H_WLOWER(MAXSTREAMS)

!  Green's function particular integral arrays

      REAL(fpk)    :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)
      REAL(fpk)    :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: DMI(MAXSTREAMS), DPI(MAXSTREAMS)
      REAL(fpk)    :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: BGP(MAXSTREAMS,MAXLAYERS)

!  Layer ! and D functions
!  Green function Multipliers for solution

      REAL(fpk)    :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Saved help variables

      REAL(fpk)    :: W_HELP(0:MAXMOMENTS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk)    :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk)    :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (part layer)

      REAL(fpk)    :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      LOGICAL      :: FLAGS_GMULT(MAX_PARTLAYERS)
      REAL(fpk)    :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Output from Boundary Value Problem
!  ----------------------------------

!  Column vectors for solving BCs

      REAL(fpk)    :: COL2    (MAXTOTAL,MAXBEAMS)
      REAL(fpk)    :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(fpk)    :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Solution constants of integration, and related quantities

      REAL(fpk)    :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Post-processing variables
!  -------------------------

!  BOA source terms

      REAL(fpk)    :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk)    :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(fpk)    :: IDOWNSURF(MAXSTREAMS)

!  TOA source term

      REAL(fpk)    :: TOA_SOURCE(MAX_USER_STREAMS)

!  Quadrature-defined solutions

      REAL(fpk)  :: QUADINTENS (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Cumulative source terms

      REAL(fpk)  :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(fpk)  :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Local Linearized arrays for argument passing, each Fourier or M = 0
!  ===================================================================

!  Linearized (Positive) Eigenvalues

      REAL(fpk)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvector solutions

      REAL(fpk)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized norms for the Green function solution

      REAL(fpk)  :: L_NORM_SAVED(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk)  :: L_HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk)  :: LP_U_WPOS1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: LP_U_WNEG1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk)  :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

      REAL(fpk)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized TOA and BOA source terms

      REAL(fpk)  :: LP_TOA_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: LP_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: LP_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Profile weighting functions at quadrature angles

      REAL(fpk)  :: QUADPROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Linearized Green's function multipliers for off-grid optical depths

      LOGICAL    :: FLAGS_LP_GMULT(MAX_PARTLAYERS)
      REAL(fpk)  :: LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Additional variables for the thermal solutions
!  ----------------------------------------------

!  Solutions to the Thermal RT equations

      REAL(fpk)       :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)       :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)       :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Linearized Solutions to the Thermal RT equations

      REAL(fpk)       :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_UT_T_PARTIC &
         (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Saved quantities for the Green function solution

      REAL(fpk)       :: TTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: T_C_MINUS(MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk)       :: T_C_PLUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk)       :: L_TTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_T_C_MINUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_T_C_PLUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Layer source terms (direct + diffuse). Upwelling

      REAL(fpk)       :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Layer source terms (direct + diffuse). Upwelling

      REAL(fpk)       :: L_LAYER_TSUP_UP   ( MAX_USER_STREAMS, MAXLAYERS,      MAX_ATMOSWFS)
      REAL(fpk)       :: L_LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  Layer source terms (direct + diffuse). Downwelling

      REAL(fpk)       :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Layer source terms (direct + diffuse). Downwelling

      REAL(fpk)       :: L_LAYER_TSUP_DN   ( MAX_USER_STREAMS, MAXLAYERS,      MAX_ATMOSWFS)
      REAL(fpk)       :: L_LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  Thermal transmittance-only source for BOA. + Linearization

      REAL(fpk)       :: BOA_THTONLY_SOURCE ( MAXSTREAMS )
      REAL(fpk)       :: L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@

!FROM LIDORT_LSDL_DBSETUPS
      REAL(fpk)       :: LSSL_DIRECT_BEAM      ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk)       :: LSSL_USER_DIRECT_BEAM ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Local help variables
!  --------------------

!   @@@@@@@@ Rob Fix 2/11/11, define Omega thermal

      REAL(fpk)       :: L_OMEGA_THERMAL ( MAX_ATMOSWFS, MAXLAYERS )

!  Local direct beam reflectance

      LOGICAL         :: DO_LOCALBEAM ( MAXBEAMS )

!  Indices

      INTEGER         :: LAYER, IBEAM, I, N, Q

!  local inclusion flags

      LOGICAL         :: DO_INCLUDE_MVOUTPUT
      LOGICAL         :: DO_INCLUDE_DIRECTBEAM
      LOGICAL         :: DO_INCLUDE_SURFACE
      LOGICAL         :: DO_INCLUDE_SURFACEWF

      LOGICAL         :: DO_INCLUDE_SURFEMISS
      LOGICAL         :: DO_INCLUDE_THERMEMISS

!  Flux multiplier and Fourier component numbers

      REAL(fpk)       :: FLUX_MULTIPLIER
      REAL(fpk)       :: DELTA_FACTOR
      REAL(fpk)       :: SURFACE_FACTOR, ALBEDO

!  Local linearization control

      LOGICAL         :: DOVARY
      LOGICAL         :: DO_SIMULATION_ONLY
      INTEGER         :: N_PARAMETERS, LAYER_TO_VARY, N_LAYER_WFS

!  error tracing

      INTEGER         :: STATUS_SUB
      CHARACTER*2     :: CF

!  progress

      LOGICAL, PARAMETER :: LPS_CHK = .FALSE.
      !LOGICAL, PARAMETER :: LPS_CHK = .TRUE.

!  debug
      INTEGER :: UT

!  ##############
!  initialization
!  ##############

!  module status and message initialization

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  Set local flags
!  ---------------

!  Local simulation only
!   Definition extended, 18 March 2014

      DO_SIMULATION_ONLY = ( .NOT. DO_PROFILE_LINEARIZATION .AND. &
                             .NOT. DO_SURFACE_LINEARIZATION .and. &
                   .not.do_ATMOS_LBBF .and. .not. do_surface_LBBF )

!  Albedo (Lambertian case). Dark surface = Lambertian with albedo 0

      ALBEDO = LAMBERTIAN_ALBEDO

!  inclusion of thermal surface emission term. only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)
!    should be true for every component if the BRDF case

      DO_INCLUDE_SURFACEWF = .TRUE.
      DO_INCLUDE_SURFACE   = .TRUE.
      IF ( .not. DO_BRDF_SURFACE ) THEN
        IF ( FOURIER_COMPONENT .NE. 0 ) THEN
          DO_INCLUDE_SURFACE   = .FALSE.
          DO_INCLUDE_SURFACEWF = .FALSE.
!mick fix 1/24/12 - lower portion of IF block commented out.
!                   "DO_INCLUDE_SURFACE = .TRUE." needed when doing
!                   analytic surfacewf even when ALBEDO = ZERO
        !ELSE
        !  IF ( ALBEDO .EQ. ZERO ) THEN
        !    DO_INCLUDE_SURFACE = .FALSE.
        !  ENDIF
        ENDIF
      ENDIF

!  Direct beam flag (only if above albedo flag has been set)

      DO IBEAM = 1, NBEAMS
        DO_LOCALBEAM(IBEAM) = .FALSE.
      ENDDO
      IF ( DO_SOLAR_SOURCES .and. DO_INCLUDE_SURFACE ) THEN
        DO IBEAM = 1, NBEAMS
          DO_LOCALBEAM(IBEAM) = DO_REFLECTED_DIRECTBEAM(IBEAM)
        ENDDO
      ENDIF

!  surface reflectance factors

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

!  Flux multipliers
!   = 1 / 4.pi with beam sources,  = 1 for Thermal alone.

      FLUX_MULTIPLIER   = DELTA_FACTOR

!  inclusion of mean value output

      DO_INCLUDE_MVOUTPUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_MVOUTPUT = .TRUE.
        ENDIF
      ENDIF

!  Initialise BVP telescoping. Important.

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        DO_BVTEL_INITIAL = DO_BVP_TELESCOPING
      ENDIF

!  @@@@@@@@ Rob Fix 2/11/11, define Omega thermal
!  mick fix 8/13/11 - added if block wrapper
!  mick fix 1/19/12 - refined if block computation and moved it below
!                     defining of DO_INCLUDE_THERMEMISS
!  set omega thermal

!      IF (DO_PROFILE_LINEARIZATION) THEN
!        DO N  = 1, MAXLAYERS
!           DO Q = 1, MAX_ATMOSWFS
!              L_OMEGA_THERMAL(Q,N) = L_OMEGA_MOMS(Q,N,0)/OMEGA_MOMS(N,0)
!           ENDDO
!        ENDDO
!      ENDIF

      IF ( DO_PROFILE_LINEARIZATION .AND. DO_INCLUDE_THERMEMISS ) THEN
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_OMEGA_THERMAL(Q,N) = L_OMEGA_MOMS(Q,N,0)/OMEGA_MOMS(N,0)
            ENDDO
          ENDIF
        ENDDO
      END IF

!  Reflected Direct beam attenuation

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL LIDORT_DIRECTBEAM                           &
        ( DO_OBSERVATION_GEOMETRY,                       & ! Input
          DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,           & ! input
          DO_REFRACTIVE_GEOMETRY, DO_USER_STREAMS,       & ! input
          NSTREAMS, N_USER_STREAMS, NBEAMS, NLAYERS,     & ! input
          FOURIER_COMPONENT, DELTA_FACTOR, FLUX_FACTOR,  & ! input
          BEAM_COSINES, SUNLAYER_COSINES,                & ! input
          ALBEDO, BRDF_F_0, USER_BRDF_F_0,               & ! input
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,           & ! input
          SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0, & ! input
          SOLAR_BEAM_OPDEP, DO_LOCALBEAM,                & ! input
          ATMOS_ATTN, DIRECT_BEAM, USER_DIRECT_BEAM )      ! Output

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@ Addition of SLEAVE weighting function setups @@@@@@@@@@@@
!    R. Spurr, 22 August 2012

        IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then
          CALL LIDORT_LSSL_DBSETUPS ( &
            DO_OBSERVATION_GEOMETRY, DO_SL_ISOTROPIC, DO_USER_STREAMS,    &
            DO_REFLECTED_DIRECTBEAM, FOURIER_COMPONENT,                   &
            NSTREAMS, NBEAMS, N_USER_STREAMS, N_SLEAVE_WFS,               &
            FLUX_FACTOR, DELTA_FACTOR,                                    &
            LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0, &
            LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM )
        ENDIF

!@@@@@@@@@@ END Addition of SLEAVE weighting function setups @@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  End solar sources only clause for corrections

      ENDIF

!  Get Legendre polynomials for this Fourier component
!   --Not requried for the thermal transmittance case

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL LIDORT_LEGENDRE_SETUP                                       &
           ( DO_REFRACTIVE_GEOMETRY, FOURIER_COMPONENT,                  & ! Input
             NSTREAMS, NBEAMS, NMOMENTS, NLAYERS,                        & ! Input
             BEAM_COSINES, SUNLAYER_COSINES, QUAD_STREAMS, QUAD_WEIGHTS, & ! Input
             PLMI_PLMJ_P, PLMI_PLMJ_M, PLMI_X0_P, PLMI_X0_M,             & ! Output
             LEG_P, LEG_M, LEG0_P, LEG0_M, WT_LEGP, WT_LEGM )              ! Output

        IF ( DO_USER_STREAMS ) THEN
          CALL LIDORT_USERLEGENDRE_SETUP         &
             ( N_USER_STREAMS, NMOMENTS,         & ! Input
               USER_STREAMS, FOURIER_COMPONENT,  & ! Input
               U_LEG_P, U_LEG_M )                  ! Output
        ENDIF
      ENDIF

!  ########################################
!  RT differential equation Eigensolutions
!  ########################################

      IF (LPS_CHK) write(*,*) 'at hom solution'

!  Avoid this section if no Scattering

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Start layer loop

       DO LAYER = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer

        CALL LIDORT_HOM_SOLUTION                                  &
          ( DO_SOLUTION_SAVING, NSTREAMS, NMOMENTS,               & ! Input
            LAYER, FOURIER_COMPONENT, DO_LAYER_SCATTERING,        & ! Input
            OMEGA_MOMS, QUAD_STREAMS, QUAD_WEIGHTS,               & ! Input
            PLMI_PLMJ_P, PLMI_PLMJ_M,                             & ! Input
            SAB, DAB, EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE,  & ! Output
            KEIGEN, XPOS, XNEG,                                   & ! Output
            STATUS_SUB, MESSAGE, TRACE_1 )                          ! Output

!  .. error tracing

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 =  'Called in LIDORT_LPS_FOURIER_MASTER, Fourier component '//CF
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Get Post-processing ("user") solutions for this layer

        IF  ( STERM_LAYERMASK_UP(LAYER) .OR. &
              STERM_LAYERMASK_DN(LAYER) ) THEN
          IF ( DO_USER_STREAMS ) THEN
            CALL LIDORT_HOM_USERSOLUTION                        &
          ( NSTREAMS, N_USER_STREAMS, NMOMENTS,                 & ! Input
            LAYER, FOURIER_COMPONENT, DO_LAYER_SCATTERING,      & ! Input
            XPOS, XNEG, OMEGA_MOMS, WT_LEGP, WT_LEGM, U_LEG_P,  & ! Input
            U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )                  ! Output
          ENDIF
        ENDIF

!  Linearized solutions
!  --------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Parameter control

          DOVARY        = LAYER_VARY_FLAG(LAYER)
          N_PARAMETERS  = LAYER_VARY_NUMBER(LAYER)

!  Linearized discrete ordinate solutions

          CALL LIDORT_L_HOM_SOLUTION                      &
         ( DO_SOLUTION_SAVING, NSTREAMS, NMOMENTS,        & ! Input
           LAYER, FOURIER_COMPONENT, DO_LAYER_SCATTERING, & ! Input
           DOVARY, N_PARAMETERS, L_OMEGA_MOMS,            & ! Input
           QUAD_STREAMS, QUAD_WEIGHTS,                    & ! Input
           PLMI_PLMJ_P, PLMI_PLMJ_M, KEIGEN, SAB, DAB,    & ! Input
           EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE,     & ! Input
           L_KEIGEN, L_XPOS, L_XNEG,                      & ! Input
           STATUS_SUB, MESSAGE, TRACE_1 )                   ! Output

!  .. error tracing

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER_COMPONENT
            TRACE_2 =  'Called in LIDORT_LPS_FOURIER_MASTER, Fourier component '//CF
            STATUS  = LIDORT_SERIOUS
            RETURN
          ENDIF

!  Get Linearizations of ("user") solutions for this layer

          CALL LIDORT_L_HOM_USERSOLUTION                    & 
         ( NSTREAMS, N_USER_STREAMS, NMOMENTS, LAYER,       & ! Input
           FOURIER_COMPONENT, DOVARY, N_PARAMETERS,         & ! Input
           STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,          & ! Input
           DO_LAYER_SCATTERING, OMEGA_MOMS, L_OMEGA_MOMS,   & ! Input
           L_XPOS, L_XNEG, WT_LEGP, WT_LEGM,                & ! Input
           U_LEG_P, U_HELP_P, U_HELP_M,                     & ! Input
           L_U_XPOS, L_U_XNEG  )                              ! Output

!  End linearization clause

        ENDIF

!  end layer loop

      ENDDO

!  prepare eigenstream tranmsittances, and linearizations

      CALL LIDORT_HOM_EIGENTRANS                                         &
        ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,          & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
          FOURIER_COMPONENT, DO_LAYER_SCATTERING,                        & ! Input
          DELTAU_VERT, PARTAU_VERT, KEIGEN,                              & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                & ! Input
          T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN )                   ! Output

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_L_HOM_EIGENTRANS                                  &
        ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,                       & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,               & ! Input
          FOURIER_COMPONENT, DO_LAYER_SCATTERING, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, & ! Input
          DELTAU_VERT,   PARTAU_VERT,  KEIGEN, L_DELTAU_VERT, L_KEIGEN,               & ! Input
          T_DELT_EIGEN,     T_UTUP_EIGEN,     T_UTDN_EIGEN,                           & ! Input
          L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_DISORDS_UTDN,                       & ! Input
          L_T_DELT_EIGEN,   L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN )                          ! Output
      ENDIF

!  Prepare solution norms, and linearizations

      CALL LIDORT_HOM_NORMS                          &
           ( NSTREAMS, NLAYERS, QUAD_STRMWTS, XPOS,  & ! Input
             NORM_SAVED )                              ! Output

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_L_HOM_NORMS                               &
          ( NSTREAMS, NLAYERS, QUAD_STRMWTS,                  & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, XPOS, L_XPOS, & ! Input
            L_NORM_SAVED )                                      ! Output
      ENDIF

!  Prepare homogeneous solution multipliers, and linearizations

!Rob Fix 5/6/13   - Now calculating ZETA_M and ZETA_P as output
!Rob Fix 5/6/13   - Introduce limiting case scenarios (Taylor series)
!Rob Fix 10/10/13 - Introduce Taylor order parameter, finalize Taylor expansions for Version 3.7

      CALL HMULT_MASTER                                                   &
        ( FOURIER_COMPONENT, DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER,    & ! Input
          NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,               & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! Input
          USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,           & ! Input
          KEIGEN, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_USERM, & ! Input
          T_UTUP_USERM, T_UTDN_USERM, DELTAU_VERT, PARTAU_VERT,           & ! Input
          ZETA_M, ZETA_P, HMULT_1, HMULT_2,                               & ! Output
          UT_HMULT_UU, UT_HMULT_UD,                                       & ! Output
          UT_HMULT_DU, UT_HMULT_DD )                                        ! Output

!     if ( fourier_component.eq.0 ) then
!        i = 7 ; layer = 15 ; ut = 2
!        write(*,*)'HMULT_1',HMULT_1(I,3,layer),HMULT_1(I,4,layer),HMULT_1(I,5,layer)
!        write(*,*)'UT_HMULT_UU',UT_HMULT_UU(I,3,ut),UT_HMULT_UU(I,4,ut),UT_HMULT_UU(I,5,ut)
!        write(*,*)'UT_HMULT_DD',UT_HMULT_DD(I,3,ut),UT_HMULT_DD(I,4,ut),UT_HMULT_DD(I,5,ut)
!        write(*,*)'ZETA_M ',ZETA_M(I,3,layer), ZETA_M(I,4,layer), ZETA_M(I,5,layer)
!     endif

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL L_HMULT_MASTER                                             &
       ( DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER,                      & ! Input
         NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,              & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,          & ! Input
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, ZETA_M, ZETA_P, L_KEIGEN,  & ! Input
         DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,                       & ! Input
         T_DELT_EIGEN, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,        & ! Input
         HMULT_1, UT_HMULT_UU, UT_HMULT_UD,                             & ! Input
         HMULT_2, UT_HMULT_DU, UT_HMULT_DD,                             & ! Input
         L_T_DELT_EIGEN,   L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,            & ! Input
         L_T_DELT_USERM,   L_T_UTUP_USERM,   L_T_UTDN_USERM,            & ! Input
         L_HMULT_1, L_UT_HMULT_UU, L_UT_HMULT_UD,                       & ! Output
         L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD )                        ! Output
      ENDIF

!  ############################################
!   boundary value problem - MATRIX PREPARATION
!  ############################################

!      write(*,*)'bvpmatrix',BVP_REGULAR_FLAG(FOURIER_COMPONENT)

      IF (LPS_CHK) write(*,*) 'at bvp prep'

!  standard case using compression of band matrices, etc..

      IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

        CALL BVP_MATRIXSETUP_MASTER                                            &
        ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER_COMPONENT,                  & ! Inputs
          NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,             & ! Inputs
          QUAD_STRMWTS, SURFACE_FACTOR, ALBEDO, BRDF_F, XPOS, XNEG, T_DELT_EIGEN,  & ! Inputs
          H_XPOS, H_XNEG, LCONMASK, MCONMASK,                                      & ! Output
          BMAT_ROWMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                          & ! Output
          STATUS_SUB, MESSAGE, TRACE_1 )                                             ! Output

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 = 'Error from BVP_MATRIXSETUP_MASTER, '//   &
               'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
          STATUS = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Telescoped case

      ELSE

        CALL BVPTEL_MATRIXSETUP_MASTER                                &
           ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS,                   & ! Input
             NSTREAMS_2, N_SUPDIAG, N_SUBDIAG, FOURIER_COMPONENT,     & ! Input
             DO_LAYER_SCATTERING, XPOS, XNEG, T_DELT_EIGEN,           & ! Input
             DO_BVTEL_INITIAL, NLAYERS_TEL,                           & ! Output
             ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,                       & ! Output
             BTELMAT_ROWMASK, BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, & ! Output
             STATUS_SUB, MESSAGE, TRACE_1 )                             ! Output

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 = 'Error from BVPTEL_MATRIXSETUP_MASTER, '// &
           'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
          STATUS = LIDORT_SERIOUS
          RETURN
        ENDIF

       ENDIF

!  End scattering solution clause

      ENDIF

!  ################
!  Thermal Solution
!  ################

!    All solutions will be scaled up by factor 4.pi if solar beams
!       are also included in the solution. Linearizations if flagged
!     Added TCUTOFF, 5/14/15 for use with Thermal Solutions

      IF ( DO_INCLUDE_THERMEMISS ) THEN

        IF (LPS_CHK) write(*,*) 'at thermal prep'

!  1. Find the Green's function solution (also thermal transmittance only)
!     I     omega_moms(1,0), L_omega_moms(1,1,0),     ! replaced line

        IF ( DO_PROFILE_LINEARIZATION ) THEN
          CALL THERMAL_GFSOLUTION_PLUS                           &
          ( DO_UPWELLING, DO_DNWELLING, DO_ADDITIONAL_MVOUT,     & ! Input
            DO_THERMAL_TRANSONLY, DO_PROFILE_LINEARIZATION,      & ! Input
            NSTREAMS, N_THERMAL_COEFFS, NLAYERS, N_PARTLAYERS,   & ! Input
            PARTLAYERS_LAYERIDX, QUAD_STREAMS, QUAD_WEIGHTS,     & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                  & ! Input
            OMEGA_MOMS(1,0), L_OMEGA_THERMAL,                    & ! Input
            DELTAU_POWER, XTAU_POWER, THERMCOEFFS, TCUTOFF,      & ! Input
            L_DELTAU_POWER, L_XTAU_POWER, L_THERMCOEFFS,         & ! Input
            T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,      & ! Input
            T_DELT_EIGEN,   T_UTDN_EIGEN,   T_UTUP_EIGEN,        & ! Input
            KEIGEN, XPOS, NORM_SAVED,                            & ! Input
            L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,& ! Input
            L_T_DELT_EIGEN,   L_T_UTDN_EIGEN,   L_T_UTUP_EIGEN,  & ! Input
            L_KEIGEN, L_XPOS, L_NORM_SAVED,                      & ! Input
            T_C_PLUS, T_C_MINUS, TTERM_SAVE,                     & ! Output
            UT_T_PARTIC, T_WUPPER, T_WLOWER,                     & ! Output
            L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,               & ! Output
            L_UT_T_PARTIC, L_T_WUPPER, L_T_WLOWER )                ! Output
        ELSE
          CALL THERMAL_GFSOLUTION                                &
          ( DO_UPWELLING, DO_DNWELLING, DO_ADDITIONAL_MVOUT,     & ! Input
            DO_THERMAL_TRANSONLY, NSTREAMS, N_THERMAL_COEFFS,    & ! Input
            NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,          & ! Input
            QUAD_STREAMS, QUAD_WEIGHTS, OMEGA_MOMS(1,0),         & ! Input
            DELTAU_POWER, XTAU_POWER, THERMCOEFFS, TCUTOFF,      & ! Input
            T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,      & ! Input
            T_DELT_EIGEN,   T_UTDN_EIGEN,   T_UTUP_EIGEN,        & ! Input
            KEIGEN, XPOS, NORM_SAVED,                            & ! Input
            T_C_PLUS, T_C_MINUS, TTERM_SAVE,                     & ! Output
            UT_T_PARTIC, T_WUPPER, T_WLOWER )                      ! Output
        ENDIF

!  2. Compute thermal layer source terms. Upwelling
!     Added TCUTOFF, 5/14/15 for use with Thermal Solutions

        IF ( DO_UPWELLING ) THEN
          IF ( DO_PROFILE_LINEARIZATION ) THEN
            CALL THERMAL_STERMS_UP_PLUS                              &
            ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY,                & ! Input
              DO_PROFILE_LINEARIZATION, STERM_LAYERMASK_UP,          & ! Input
              NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,            & ! Input
              NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS, TCUTOFF,   & ! Input
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, USER_STREAMS,      & ! Input
              DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,& ! Input
                T_DELT_USERM,   T_UTUP_USERM,   U_XPOS,   U_XNEG,    & ! Input
              L_T_DELT_USERM, L_T_UTUP_USERM, L_U_XPOS, L_U_XNEG,    & ! Input
                T_C_PLUS,   T_C_MINUS,   TTERM_SAVE,                 & ! Input
              L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,                 & ! Input
                T_DIRECT_UP,   T_UT_DIRECT_UP,                       & ! Input
              L_T_DIRECT_UP, L_T_UT_DIRECT_UP,                       & ! Input
                HMULT_1,   HMULT_2,   UT_HMULT_UU,   UT_HMULT_UD,    & ! Input
              L_HMULT_1, L_HMULT_2, L_UT_HMULT_UU, L_UT_HMULT_UD,    & ! Input
                LAYER_TSUP_UP,   LAYER_TSUP_UTUP,                    & ! Output
              L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP )                     ! Output
          ELSE
            CALL THERMAL_STERMS_UP                                         &
            ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, STERM_LAYERMASK_UP,  & ! Input
              NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,                  & ! Input
              NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS, TCUTOFF,         & ! Input
              USER_STREAMS, DELTAU_POWER, XTAU_POWER,                      & ! Input
              T_DELT_USERM, T_UTUP_USERM, U_XPOS, U_XNEG,                  & ! Input
              T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_DIRECT_UP, T_UT_DIRECT_UP,& ! Input
              HMULT_1, HMULT_2, UT_HMULT_UU,  UT_HMULT_UD,                 & ! Input
              LAYER_TSUP_UP, LAYER_TSUP_UTUP )                               ! Output
           ENDIF
         ENDIF

!  3. Compute thermal layer source terms. Downwelling
!     Added TCUTOFF, 5/14/15 for use with Thermal Solutions

        IF ( DO_DNWELLING ) THEN
          IF ( DO_PROFILE_LINEARIZATION ) THEN
            CALL THERMAL_STERMS_DN_PLUS                              &
            ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY,                & ! Input
              DO_PROFILE_LINEARIZATION, STERM_LAYERMASK_DN,          & ! Input
              NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,            & ! Input
              NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS, TCUTOFF,   & ! Input
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, USER_STREAMS,      & ! Input
              DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,& ! Input
                T_DELT_USERM,   T_UTDN_USERM,   U_XPOS,   U_XNEG,    & ! Input
              L_T_DELT_USERM, L_T_UTDN_USERM, L_U_XPOS, L_U_XNEG,    & ! Input
                T_C_PLUS,   T_C_MINUS,   TTERM_SAVE,                 & ! Input
              L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,                 & ! Input
                T_DIRECT_DN,   T_UT_DIRECT_DN,                       & ! Input
              L_T_DIRECT_DN, L_T_UT_DIRECT_DN,                       & ! Input
                HMULT_1,   HMULT_2,   UT_HMULT_DU,   UT_HMULT_DD,    & ! Input
              L_HMULT_1, L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD,    & ! Input
                LAYER_TSUP_DN,   LAYER_TSUP_UTDN,                    & ! Output
              L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN )                     ! Output
          ELSE
            CALL THERMAL_STERMS_DN                                         &
            ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, STERM_LAYERMASK_DN,  & ! Input
              NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,                  & ! Input
              NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS, TCUTOFF,         & ! Input
              USER_STREAMS, DELTAU_POWER, XTAU_POWER,                      & ! Input
              T_DELT_USERM, T_UTDN_USERM, U_XPOS, U_XNEG,                  & ! Input
              T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_DIRECT_DN, T_UT_DIRECT_DN,& ! Input
              HMULT_1, HMULT_2, UT_HMULT_DU,  UT_HMULT_DD,                 & ! Input
              LAYER_TSUP_DN, LAYER_TSUP_UTDN )                               ! Output
          ENDIF
        ENDIF

      ENDIF

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Skip this thermal-only section if there are solar sources

      IF ( .not. DO_SOLAR_SOURCES ) THEN

        IF (LPS_CHK) write(*,*) 'at thermal rads'

!  Only one solution, local direct_beam flag NOT set

        IBEAM = 1
        DO_INCLUDE_DIRECTBEAM = .FALSE.

!  Avoid the thermal scattering solutions if "transmittance-only"

        IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Thermal-only. Find the BVP solution and intensity field
!  set the BVP PI solution at the lower/upper boundaries

          DO LAYER = 1, NLAYERS
            DO I = 1, NSTREAMS_2
              WUPPER(I,LAYER) = T_WUPPER(I,LAYER)
              WLOWER(I,LAYER) = T_WLOWER(I,LAYER)
            ENDDO
          ENDDO

!  Solve the boundary value problem

          IF (LPS_CHK) write(*,*) 'at thermal bvp'

!mick fix 6/29/11 - removed COL2 & SCOL2 from call
!mick fix 4/17/2014 - reinserted COL2 & SCOL2
          CALL BVP_SOLUTION_MASTER                                    &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,& ! Input
            DO_INCLUDE_DIRECTBEAM, NSTREAMS, NLAYERS, NSTREAMS_2,     & ! Input
            NTOTAL, N_SUPDIAG, N_SUBDIAG, FOURIER_COMPONENT, IBEAM,   & ! Input
            QUAD_STRMWTS, SURFACE_FACTOR, ALBEDO, BRDF_F,             & ! Input
            SURFACE_BB_INPUT, EMISSIVITY,                             & ! Input
            XPOS, XNEG, WUPPER, WLOWER, DIRECT_BEAM,                  & ! Input
            BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                         & ! Input
            COL2, SCOL2,                                              & ! Input
            H_WLOWER, LCON, MCON, LCON_XVEC, MCON_XVEC,               & ! Output
            STATUS_SUB, MESSAGE, TRACE_1 )                              ! Output

!  Error handling

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER_COMPONENT
            TRACE_2 = 'Error return from BVP_SOLUTION_MASTER, '// &
                      'Called in LIDORT_FOURIER_MASTER, Fourier # '//CF
            STATUS = LIDORT_SERIOUS
            RETURN
          ENDIF

!  Finish calculating thermal scattering solution

        ENDIF

        IF (LPS_CHK) write(*,*) 'at thermal post proc'

!  Post-processing - upwelling thermal-only field

        IF ( DO_UPWELLING ) THEN

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE

          CALL GET_BOASOURCE                                            &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_THERMAL,    & ! Input @@@@
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM, & ! Input
            DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS, DO_OBSERVATION_GEOMETRY, & ! Input
            NSTREAMS, NLAYERS, N_USER_STREAMS, IBEAM, FOURIER_COMPONENT,& ! Input
            SURFACE_FACTOR, QUAD_STRMWTS, QUAD_WEIGHTS, T_DELT_DISORDS, & ! Input
            LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_WLOWER,       & ! Input
            ALBEDO, BRDF_F, USER_BRDF_F, USER_DIRECT_BEAM,              & ! Input
            SURFACE_BB_INPUT, EMISSIVITY, USER_EMISSIVITY,              & ! Input
            BOA_SOURCE, DIRECT_BOA_SOURCE,                              & ! Output
            BOA_THTONLY_SOURCE, IDOWNSURF )                               ! Output

!Rob fix 5/6/13 - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced

          CALL UPUSER_INTENSITY                                    &
          ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,   & ! Input
            DO_OBSERVATION_GEOMETRY, DO_LAYER_SCATTERING,          & ! Input
            DO_THERMAL_EMISSION,                                   & ! Input
            DO_THERMAL_TRANSONLY, DO_INCLUDE_MVOUTPUT,             & ! Input
            FOURIER_COMPONENT, NSTREAMS, N_USER_STREAMS, NLAYERS,  & ! Input
            N_USER_LEVELS, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,& ! Input
            UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,               & ! Input
            IBEAM, TAYLOR_ORDER, FLUX_MULTIPLIER, QUAD_STREAMS,    & ! Input
            GAMMA_P, GAMMA_M, SIGMA_P,                             & ! Input
            DELTAU_VERT, PARTAU_VERT,                              & ! Input
            INITIAL_TRANS, ITRANS_USERM,                           & ! Input
            T_DELT_USERM, T_UTUP_USERM,                            & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR,                            & ! Input
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,              & ! Input
            T_DELT_DISORDS, T_DISORDS_UTUP,                        & ! Input
            T_WUPPER, UT_T_PARTIC, LAYER_TSUP_UP, LAYER_TSUP_UTUP, & ! Input
            LAYER_PIS_CUTOFF, ATERM_SAVE, BTERM_SAVE,              & ! Input
            XPOS, WUPPER, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,& ! Input
            U_XPOS, U_XNEG, U_WPOS1, HMULT_1, HMULT_2, EMULT_UP,   & ! Input
            UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                & ! Input
            BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE,     & ! Input
            PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,          & ! Output
            FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                 & ! Output
            INTENSITY_F, QUADINTENS, CUMSOURCE_UP )                  ! Output

        ENDIF

!  Post-processing - Downwelling thermal-only field

        IF ( DO_DNWELLING ) THEN

          CALL GET_TOASOURCE ( DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, IBEAM, TOA_SOURCE )

!Rob fix 5/6/13   - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced

          CALL DNUSER_INTENSITY                                       &
          ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,      & ! Input
            DO_OBSERVATION_GEOMETRY, DO_LAYER_SCATTERING,             & ! Input
            DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY,                & ! Input
            DO_INCLUDE_MVOUTPUT,                                      & ! Input
            FOURIER_COMPONENT, NSTREAMS, N_USER_STREAMS,              & ! Input
            N_USER_LEVELS, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,   & ! Input
            UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                  & ! Input
            IBEAM, TAYLOR_ORDER, FLUX_MULTIPLIER, QUAD_STREAMS,       & ! Input
            GAMMA_P, GAMMA_M, SIGMA_M,                                & ! Input
            DELTAU_VERT, PARTAU_VERT,                                 & ! Input
            INITIAL_TRANS, ITRANS_USERM,                              & ! Input
            T_DELT_USERM, T_UTDN_USERM,                               & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR,                               & ! Input
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                 & ! Input
            T_DELT_DISORDS, T_DISORDS_UTDN,                           & ! Input
            T_WLOWER, UT_T_PARTIC, LAYER_TSUP_DN, LAYER_TSUP_UTDN,    & ! Input
            LAYER_PIS_CUTOFF, ATERM_SAVE, BTERM_SAVE,                 & ! Input
            XPOS, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,           & ! Input
            U_XPOS, U_XNEG, U_WNEG1, HMULT_1, HMULT_2, EMULT_DN,      & ! Input
            UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN, TOA_SOURCE,       & ! Input
            PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,             & ! Output
            FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                    & ! Output
            INTENSITY_F, QUADINTENS, CUMSOURCE_DN )                     ! Output

        ENDIF

!  mean value output for thermal-only field
!    Direct-beam contributions output separately, 26 May 11, 24 August 2011

        IF ( DO_INCLUDE_MVOUTPUT ) THEN

          IF (LPS_CHK) write(*,*) 'at thermal mean vals'

          CALL MIFLUX_INTENSITY                                 &
          ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTBEAM,  & ! Input
           IBEAM, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,         & ! Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,            & ! Input
            UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,            & ! Input
            QUAD_WEIGHTS, QUAD_STRMWTS,                         & ! Input
            INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,        & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,             & ! Input
            MEAN_INTENSITY, FLUX_INTEGRAL,                      & ! Output
            DNMEAN_DIRECT,  DNFLUX_DIRECT )                       ! Output

        ENDIF

!  Thermal only. Avoid weighting functions all together
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        IF ( .not. DO_SIMULATION_ONLY ) THEN

         IF (LPS_CHK) write(*,*) 'at thermal jacs'

!  Thermal Only: Atmospheric Profile weighting functions
!  =====================================================

         IF ( DO_PROFILE_LINEARIZATION ) THEN

          IF (LPS_CHK) write(*,*) 'at thermal lp bvp'

!  Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations

          DO LAYER_TO_VARY = 1, NLAYERS
           IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
            N_LAYER_WFS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

!  Avoidance of BVP problem for transmittance only

            IF ( .not. DO_THERMAL_TRANSONLY ) THEN

            IF (LPS_CHK) write(*,*) 'at thermal lp post proc'

!  Solve the linearized BVP (No telescoped solution)
!  Rob Fix, 01/06/14, Reorganized calling arguments, including TAYLOR_ORDER

            CALL LP_BVP_SOLUTION_MASTER &
            ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_LOCALBEAM(IBEAM),     & ! Input
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_PLANE_PARALLEL,   & ! Input
              DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER_COMPONENT, IBEAM,  & ! Input, Flags and order
              NLAYERS, NSTREAMS, NSTREAMS_2, NTOTAL,                        & ! Input, Numbers
              N_SUBDIAG, N_SUPDIAG, LAYER_TO_VARY, N_LAYER_WFS,             & ! Input, Numbers
              DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
              SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,    & ! Input, Surface Stuff
              INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
              BANDMAT2, SMAT2, IPIVOT, SIPIVOT, LCONMASK, MCONMASK,         & ! Input, BVP Bandmat
              T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
              GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
              L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous  
              L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens/Thermal
              NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,         & ! Output, Linearized Constants + PI
              STATUS_SUB, MESSAGE, TRACE_1 )                                  ! Output

             IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER_COMPONENT
               TRACE_2 = 'Error return from LP_BVP_SOLUTION_MASTER, '// &
                         'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS
               RETURN
             ENDIF

!  Continuation point for avoiding scattering solution

            ENDIF

!  Post-processing for the Upwelling profile weighting functions
!  -------------------------------------------------------------

            IF ( DO_UPWELLING ) THEN

!  Profile-Linearized BOA source term

             CALL GET_LP_BOASOURCE                                           &
             ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input
               DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,& ! Input
               DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,                         & ! Input
               DO_OBSERVATION_GEOMETRY,                                      & ! Input
               NLAYERS, NSTREAMS, N_USER_STREAMS, FOURIER_COMPONENT,         & ! Input
               IBEAM, LAYER_TO_VARY, N_LAYER_WFS, QUAD_STRMWTS, QUAD_WEIGHTS,& ! Input
               SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                  & ! Input
               USER_DIRECT_BEAM, DELTAU_SLANT, L_DELTAU_VERT,                & ! Input
               LCON, MCON, LCON_XVEC, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input
               T_WLOWER, L_XPOS, L_XNEG, L_T_WLOWER, L_WLOWER,               & ! Input
               NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_T_DELT_DISORDS,       & ! Input
               LP_BOA_MSSOURCE, LP_BOA_DBSOURCE, L_BOA_THTONLY_SOURCE )        ! Output

!  Linearized upwelling PROFILE weighting functions

             CALL UPUSER_PROFILEWF &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES,                        & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,              & ! Input
              DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT,                    & ! Input
              DO_OBSERVATION_GEOMETRY, FLUX_MULTIPLIER,                 & ! Input
              NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,         & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! Input
              UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                  & ! Input
              FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS,     & ! Input
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
              U_XPOS, U_XNEG, U_WPOS1, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
              UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                    & ! Input
              PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,              & ! Input
              UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_UP,                    & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                & ! Input
              LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                  & ! Input
              L_T_DELT_USERM,   L_T_UTUP_USERM,                   & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTUP,                 & ! Input
              L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Input
              L_U_XPOS, L_U_XNEG, LP_U_WPOS1, NCON, PCON,                & ! Input
              L_XPOS, L_XNEG, L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC,  & ! Input
              L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                         & ! Input
              L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,              & ! Input
              L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP, L_UT_T_PARTIC, L_T_WUPPER, & ! Input
              L_BOA_THTONLY_SOURCE, LP_BOA_MSSOURCE, LP_BOA_DBSOURCE,        & ! Input
              FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,                & ! Output
              PROFILEWF_F, QUADPROFILEWF )                                     ! Output

            ENDIF

!  Post-processing for PROFILE weighting functions (Downwelling)
!  ------------------------------------------------------------

            IF ( DO_DNWELLING ) THEN

!  TOA source term linearized

             CALL GET_LP_TOASOURCE &
               ( DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, IBEAM, N_LAYER_WFS, &
                 LP_TOA_SOURCE )

!  Downwelling weighting functions

             CALL DNUSER_PROFILEWF                                       &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES,                         & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,               & ! Input
              DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT,                     & ! Input
              DO_OBSERVATION_GEOMETRY, FLUX_MULTIPLIER,                  & ! Input
              NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                   & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
              UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                   & ! Input
              FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS,      & ! Input
              TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,     & ! Input
              QUAD_STREAMS, USER_SECANTS, DO_LAYER_SCATTERING,           & ! Input
              LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,           & ! Input
              T_DELT_MUBAR,   T_UTDN_MUBAR,                      & ! Input
              T_DELT_USERM,   T_UTDN_USERM,                      & ! Input
              T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN,      & ! Input
              T_DELT_DISORDS, T_DISORDS_UTDN,                    & ! Input
              GAMMA_M, GAMMA_P, SIGMA_M, ATERM_SAVE, BTERM_SAVE, & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WLOWER,          & ! Input
              U_XPOS, U_XNEG, U_WNEG1, HMULT_1, HMULT_2, EMULT_DN,       & ! Input
              UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN,                    & ! Input
              PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,              & ! Input
              UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_DN,                    & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                & ! Input
              LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                  & ! Input
              L_T_DELT_USERM,   L_T_UTDN_USERM,                   & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTDN,                 & ! Input
              L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Input
              L_U_XPOS, L_U_XNEG, LP_U_WNEG1, NCON, PCON,                & ! Input
              L_XPOS, L_XNEG, L_WLOWER, NCON_XVEC, PCON_XVEC,            & ! Input
              L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                         & ! Input
              L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,              & ! Input
              L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN, L_UT_T_PARTIC,         & ! Input
              L_T_WLOWER, LP_TOA_SOURCE,                                 & ! Input
              FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,            & ! Output
              PROFILEWF_F, QUADPROFILEWF )                                 ! Output

            ENDIF

!  Post-processing for Profile weighting functions (Mean/Flux)
!  ----------------------------------------------------------

            IF ( DO_INCLUDE_MVOUTPUT ) THEN

              IF (LPS_CHK) write(*,*) 'at thermal lp mean val jacs'

              CALL MIFLUX_PROFILEWF                                &
              ( DO_UPWELLING, DO_DNWELLING, DO_LOCALBEAM(IBEAM),   & ! Input
                IBEAM, LAYER_TO_VARY, N_LAYER_WFS,                 & ! Input
                NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,              & ! Input
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,           & ! Input
                UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,           & ! Input
                QUAD_WEIGHTS, QUAD_STRMWTS,                        & ! Input
                INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,       & ! Input
                T_DELT_MUBAR, T_UTDN_MUBAR, QUADPROFILEWF,         & ! Input
                LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,& ! Input
                MINT_PROFILEWF, FLUX_PROFILEWF )                     ! Output
            ENDIF

!  Finish loop over layers with variation

           ENDIF
          ENDDO

!  End atmospheric profile weighting functions

         ENDIF

!  Surface Property weighting functions (regardless of Dark surface condition)
!  ====================================

         IF ( DO_SURFACE_LINEARIZATION .AND. DO_INCLUDE_SURFACEWF .AND. &
              (N_SURFACE_WFS > 0) ) THEN

            IF (LPS_CHK) write(*,*) 'at thermal surf jacs'

!  Surface property WF master

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add argument to SURFACEWF_MASTER

            CALL SURFACEWF_MASTER                                           &
            ( DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE,                  & ! Input
              DO_USER_STREAMS, DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY,  & ! Input
              DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT, & ! Input
              DO_MSMODE_THERMAL, DO_OBSERVATION_GEOMETRY,                   & ! Input
              NLAYERS, NTOTAL, N_SUPDIAG, N_SUBDIAG, N_USER_LEVELS,         & ! Input
              NSTREAMS, NSTREAMS_2, N_USER_STREAMS, N_SURFACE_WFS,          & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
              UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                       & ! Input
              FOURIER_COMPONENT, IBEAM, FLUX_MULTIPLIER,                    & ! Input
              SURFACE_BB_INPUT, LS_EMISSIVITY, LS_USER_EMISSIVITY,          & ! Input
              SURFACE_FACTOR, ALBEDO, LS_BRDF_F, LS_BRDF_F_0,               & ! Input
              USER_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0,                & ! Input
              QUAD_WEIGHTS, QUAD_STRMWTS, ATMOS_ATTN, IDOWNSURF,            & ! Input
              LCONMASK, MCONMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,         & ! Input
              LCON, MCON, XPOS, XNEG, H_XPOS, H_XNEG, H_WLOWER,             & ! Input
              T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                     & ! Input
              T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                     & ! Input
              T_DELT_DISORDS, T_DISORDS_UTUP, U_XPOS, U_XNEG, HMULT_1,      & ! Input
              HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,  & ! Input
              SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,                  & ! Output
              STATUS_SUB, MESSAGE, TRACE_1 )                                  ! Output

!  Error handling

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER_COMPONENT
              TRACE_2 = 'Error return from SURFACEWF_MASTER, '// &
                        'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
              STATUS = LIDORT_SERIOUS
              RETURN
            ENDIF

!  end of surface weighting functions

         ENDIF

!  New, 18 March 2014. Linearization for BLACKBODY

        IF ( ( DO_ATMOS_LBBF .or. DO_SURFACE_LBBF ) .and. Fourier_component.eq.0 ) THEN

           IF (LPS_CHK) write(*,*) 'at thermal lbbf jacs'

           if ( n_partlayers .gt. 0 ) then
              call lidort_lbbf_jacobians_wpartials &
               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,       & ! Inputs 3/25
                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,               & ! Inputs
                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,    & ! input
                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                        & ! input
                 NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,           & ! input
                 NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,                   & ! input
                 N_PARTLAYERS, PARTLAYERS_LAYERIDX,                          & ! Input 3/21
                 PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                    & ! Input 3/21
                 UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                     & ! Input
                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,       & ! Input
                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                   & ! input
                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                & ! input
                 EMISSIVITY, USER_EMISSIVITY,                                & ! input
                 FLUX_MULTIPLIER, OMEGA_MOMS(:,0), DELTAU_VERT, PARTAU_VERT, & ! Input 3/21
                 T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,             & ! inputs 3/25
                 KEIGEN, TTERM_SAVE, XPOS, XNEG,                             & ! inputs 3/21
                 T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                   & ! inputs 3/21
                 BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                           & ! Input
                 T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                   & ! Inputs 3/21
                 U_XPOS, U_XNEG, HMULT_1, HMULT_2,                           & ! inputs
                 UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,         & ! Inputs 3/21
                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                            & ! Output
                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                            & ! Output
                 STATUS_SUB, MESSAGE, TRACE_1 )                                ! Output
           else
              call lidort_lbbf_jacobians_whole &
               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,    & ! Inputs 3/25
                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,            & ! Inputs
                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, & ! input
                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                     & ! input
                 NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,        & ! input
                 NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,                & ! input
                 UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                  & ! Input
                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! Input
                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                & ! input
                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,             & ! input
                 EMISSIVITY, USER_EMISSIVITY,                             & ! inputs
                 FLUX_MULTIPLIER, OMEGA_MOMS(:,0), DELTAU_VERT,           & ! Input
                 T_DELT_DISORDS, T_DELT_EIGEN, KEIGEN, XPOS, XNEG,        & ! inputs
                 TTERM_SAVE, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Input
                 T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2,          & ! inputs
                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                         & ! Output
                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                         & ! Output
                 STATUS_SUB, MESSAGE, TRACE_1 )                             ! Output
           endif

!  Error handling

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER_COMPONENT
              TRACE_2 = 'Error return from lidort_lbbf_jacobians, thermal only'// &
                        'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
              STATUS = LIDORT_SERIOUS
              RETURN
           ENDIF

!  ENd of BB Jacobians

        ENDIF

!  End clause for simulation only

        ENDIF

!  All done. Finish Thermal only solution, so exit

        RETURN

!  End clause for thermal-only solution

      ENDIF

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

      IF (LPS_CHK) write(*,*) 'at complete rads'

!  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

!  Only calculate if still not converged

        IF ( DO_MULTIBEAM(IBEAM,FOURIER_COMPONENT) ) THEN

!  Solar beam Particular solutions (Green's function)
!  --------------------------------------------------

!  start layer loop

          DO LAYER = 1, NLAYERS

!  Control

            IF ( DO_PROFILE_LINEARIZATION ) THEN
              DOVARY        = LAYER_VARY_FLAG(LAYER)
              N_PARAMETERS  = LAYER_VARY_NUMBER(LAYER)
            ENDIF

!  Green's function solution

            CALL LIDORT_GBEAM_SOLUTION ( TAYLOR_ORDER,                  & ! Input (10/10/13, argument added)
            NSTREAMS, NSTREAMS_2, NMOMENTS, LAYER, FOURIER_COMPONENT,   & ! Input
            FLUX_FACTOR, IBEAM, DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,  & ! Input
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, DELTAU_VERT,   & ! Input
            QUAD_WEIGHTS, OMEGA_MOMS, KEIGEN, XPOS, T_DELT_EIGEN,       & ! Input
            NORM_SAVED, PLMI_X0_P, PLMI_X0_M,                           & ! Input
            GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE,         & ! Output
            CFUNC, DFUNC, AGM, BGP, GFUNC_UP, GFUNC_DN,                 & ! Output
            WUPPER, WLOWER )                                              ! Output

!  Linearized Green's function solution

            IF ( DO_PROFILE_LINEARIZATION ) THEN
              CALL LIDORT_L_GBEAM_SOLUTION                   &
         ( NSTREAMS, NMOMENTS, LAYER, FOURIER_COMPONENT,     & ! Input
           FLUX_FACTOR, IBEAM, DOVARY, N_PARAMETERS,         & ! Input
           DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,            & ! Input
           QUAD_WEIGHTS, L_OMEGA_MOMS, PLMI_X0_P, PLMI_X0_M, & ! Input
           NORM_SAVED, XPOS, L_NORM_SAVED, L_XPOS,           & ! Input
           DMI, DPI, ATERM_SAVE, BTERM_SAVE,                 & ! Input
           L_ATERM_SAVE, L_BTERM_SAVE )                        ! Output
            ENDIF

!  Post-processing solutions

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR. &
                  STERM_LAYERMASK_DN(LAYER) ) THEN

              IF ( DO_USER_STREAMS ) THEN

                CALL LIDORT_GBEAM_USERSOLUTION                          &
                ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,  & ! Input
                  N_USER_STREAMS, NMOMENTS,                             & ! Input
                  LAYER, FOURIER_COMPONENT, IBEAM, FLUX_FACTOR,         & ! Input
                  DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,                & ! Input
                  OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,                 & ! Input
                  U_WPOS1, U_WNEG1, W_HELP )                              ! Output

                IF ( DO_PROFILE_LINEARIZATION ) THEN
                  CALL LIDORT_L_GBEAM_USERSOLUTION             &
                  ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,& ! Input
                    N_USER_STREAMS, NMOMENTS, LAYER, FOURIER_COMPONENT, & ! Input
                    IBEAM, FLUX_FACTOR, DOVARY, N_PARAMETERS,           & ! Input
                    DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,              & ! Input
                    L_OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,             & ! Input
                    LP_U_WPOS1, LP_U_WNEG1 )                              ! Output
                ENDIF

!  End post-processing

              END IF
            END IF

!  end layer loop

          END DO

!  Add thermal solutions if flagged
!    Beware the 4.PI modulus on the thermal contribution (already scaled up)

          IF ( DO_INCLUDE_THERMEMISS ) THEN
           DO LAYER = 1, NLAYERS
            DO I = 1, NSTREAMS_2
             WUPPER(I,LAYER) = WUPPER(I,LAYER) + T_WUPPER(I,LAYER)
             WLOWER(I,LAYER) = WLOWER(I,LAYER) + T_WLOWER(I,LAYER)
            ENDDO
           ENDDO
          ENDIF

!  Adding the linearized thermal solutions is done later.....

!  Solve boundary value problem
!  ----------------------------

          IF (LPS_CHK) write(*,*) 'at complete bvp'

!  Get the Regular BVP solution

          IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

!mick fix 6/29/11 - removed COL2 & SCOL2 from call
!mick fix 4/17/2014 - reinserted COL2 & SCOL2
           CALL BVP_SOLUTION_MASTER                                    &
           ( DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,& ! Input
             DO_LOCALBEAM(IBEAM), NSTREAMS, NLAYERS, NSTREAMS_2,       & ! Input
             NTOTAL, N_SUPDIAG, N_SUBDIAG, FOURIER_COMPONENT, IBEAM,   & ! Input
             QUAD_STRMWTS, SURFACE_FACTOR, ALBEDO, BRDF_F,             & ! Input
             SURFACE_BB_INPUT, EMISSIVITY,                             & ! Input
             XPOS, XNEG, WUPPER, WLOWER, DIRECT_BEAM,                  & ! Input
             BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                         & ! Input
             COL2, SCOL2,                                              & ! Input
             H_WLOWER, LCON, MCON, LCON_XVEC, MCON_XVEC,               & ! Output
             STATUS_SUB, MESSAGE, TRACE_1 )                              ! Output

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
             write(CF,'(I2)')FOURIER_COMPONENT
             TRACE_2 = 'Error return from BVP_SOLUTION_MASTER, '// &
                'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
             STATUS = LIDORT_SERIOUS
             RETURN
           ENDIF 

!  Get the telescoped boundary value result

          ELSE

!mick fix 6/29/11 - removed COLTEL2 & SCOL2 from call
!mick fix 4/17/2014 - reinserted COLTEL2 & SCOL2
           CALL BVPTEL_SOLUTION_MASTER                            &
             ( IBEAM, NSTREAMS, NLAYERS,                          & ! Input
               NSTREAMS_2, N_SUBDIAG, N_SUPDIAG, WUPPER, WLOWER,  & ! Input
               XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input
               NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,    & ! Input
               BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,            & ! Input
               COLTEL2, SCOL2,                                    & ! Input
               LCON, MCON, LCON_XVEC,  MCON_XVEC,                 & ! Output
               STATUS_SUB, MESSAGE, TRACE_1 )                       ! Output

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
             write(CF,'(I2)')FOURIER_COMPONENT
             TRACE_2 = 'Error return from BVPTEL_SOLUTION_MASTER, '// &
                'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
             STATUS = LIDORT_SERIOUS
             RETURN
           ENDIF

          ENDIF

!  Radiance Field Post Processing
!  ------------------------------

          IF (LPS_CHK) write(*,*) 'at complete post proc'

!  Direct beam inclusion flag:
!   This now has the DBCORRECTION option: if the DBCORRECTION Flag
!   is set, then we will be doing exact calculations of the reflected
!   directbeam, so we do not need to include it in the Post-processing.
!   However, the direct beam will need to be included in the basic RT
!   solution (the BVP), and this is controlled separately by the
!   DO_REFLECTED_DIRECTBEAM(IBEAM) flags.
!     R. Spurr, RT Solutions, Inc., 19 August 2005.

!          DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. &
!             (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION))

          DO_INCLUDE_DIRECTBEAM = &
            ( DO_UPWELLING .AND. DO_LOCALBEAM(IBEAM) ) &
                     .AND. .NOT. DO_MSMODE_LIDORT

!  upwelling

!mick fix 6/29/11 - initialize FLAGS_GMULT
          FLAGS_GMULT = .FALSE.

          IF ( DO_UPWELLING ) THEN

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE

            CALL GET_BOASOURCE                                            &
            ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_THERMAL,    & ! Input @@@@
              DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM, & ! Input
              DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS, DO_OBSERVATION_GEOMETRY, & ! Input
              NSTREAMS, NLAYERS, N_USER_STREAMS, IBEAM, FOURIER_COMPONENT,& ! Input
              SURFACE_FACTOR, QUAD_STRMWTS, QUAD_WEIGHTS, T_DELT_DISORDS, & ! Input
              LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_WLOWER,       & ! Input
              ALBEDO, BRDF_F, USER_BRDF_F, USER_DIRECT_BEAM,              & ! Input
              SURFACE_BB_INPUT, EMISSIVITY, USER_EMISSIVITY,              & ! Input
              BOA_SOURCE, DIRECT_BOA_SOURCE,                              & ! Output
              BOA_THTONLY_SOURCE, IDOWNSURF )                               ! Output

!Rob fix 5/6/13   - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced

            CALL UPUSER_INTENSITY                                    &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,   & ! Input
              DO_OBSERVATION_GEOMETRY, DO_LAYER_SCATTERING,          & ! Input
              DO_INCLUDE_THERMEMISS,                                 & ! Input
              DO_THERMAL_TRANSONLY, DO_INCLUDE_MVOUTPUT,             & ! Input
              FOURIER_COMPONENT, NSTREAMS, N_USER_STREAMS, NLAYERS,  & ! Input
              N_USER_LEVELS, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,& ! Input
              UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,               & ! Input
              IBEAM, TAYLOR_ORDER, FLUX_MULTIPLIER, QUAD_STREAMS,    & ! Input
              GAMMA_P, GAMMA_M, SIGMA_P,                             & ! Input
              DELTAU_VERT, PARTAU_VERT,                              & ! Input
              INITIAL_TRANS, ITRANS_USERM,                           & ! Input
              T_DELT_USERM, T_UTUP_USERM,                            & ! Input
              T_DELT_MUBAR, T_UTDN_MUBAR,                            & ! Input
              T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,              & ! Input
              T_DELT_DISORDS, T_DISORDS_UTUP,                        & ! Input
              T_WUPPER, UT_T_PARTIC, LAYER_TSUP_UP, LAYER_TSUP_UTUP, & ! Input
              LAYER_PIS_CUTOFF, ATERM_SAVE, BTERM_SAVE,              & ! Input
              XPOS, WUPPER, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,& ! Input
              U_XPOS, U_XNEG, U_WPOS1, HMULT_1, HMULT_2, EMULT_UP,   & ! Input
              UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                & ! Input
              BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE,     & ! Input
              PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,          & ! Output
              FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                 & ! Output
              INTENSITY_F, QUADINTENS, CUMSOURCE_UP )                  ! Output

          ENDIF

!  Downwelling

          IF ( DO_DNWELLING ) THEN

            CALL GET_TOASOURCE ( DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, IBEAM, TOA_SOURCE )

!Rob fix 5/6/13   - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced

            CALL DNUSER_INTENSITY &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,      & ! Input
              DO_OBSERVATION_GEOMETRY, DO_LAYER_SCATTERING,             & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,              & ! Input
              DO_INCLUDE_MVOUTPUT,                                      & ! Input
              FOURIER_COMPONENT, NSTREAMS, N_USER_STREAMS,              & ! Input
              N_USER_LEVELS, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,   & ! Input 
              UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                  & ! Input
              IBEAM, TAYLOR_ORDER, FLUX_MULTIPLIER, QUAD_STREAMS,       & ! Input
              GAMMA_P, GAMMA_M, SIGMA_M,                                & ! Input
              DELTAU_VERT, PARTAU_VERT,                                 & ! Input
              INITIAL_TRANS, ITRANS_USERM,                              & ! Input
              T_DELT_USERM, T_UTDN_USERM,                               & ! Input
              T_DELT_MUBAR, T_UTDN_MUBAR,                               & ! Input
              T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                 & ! Input
              T_DELT_DISORDS, T_DISORDS_UTDN,                           & ! Input
              T_WLOWER, UT_T_PARTIC, LAYER_TSUP_DN, LAYER_TSUP_UTDN,    & ! Input
              LAYER_PIS_CUTOFF, ATERM_SAVE, BTERM_SAVE,                 & ! Input
              XPOS, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,           & ! Input
              U_XPOS, U_XNEG, U_WNEG1, HMULT_1, HMULT_2, EMULT_DN,      & ! Input
              UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN, TOA_SOURCE,       & ! Input
              PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,             & ! Output
              FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                    & ! Output
              INTENSITY_F, QUADINTENS, CUMSOURCE_DN )                     ! Output

          ENDIF

!     if ( fourier_component.eq.0 ) then
!         I = 7 ; layer = 15 ; ut = 2
!         write(*,*)'Greens',ibeam,average_secant(15,ibeam),PMULT_UD(I,3,15),UT_PMULT_UD(I,3,ut),&
!                                                           PMULT_DD(I,3,15),UT_PMULT_DD(I,3,ut) 
!     endif

!  mean value output
!    Direct-beam contributions output separately, 26 May 11, 24 August 2011

          IF ( DO_INCLUDE_MVOUTPUT ) THEN

            IF (LPS_CHK) write(*,*) 'at complete mean vals'

            CALL MIFLUX_INTENSITY                                 &
            ( DO_UPWELLING, DO_DNWELLING, DO_LOCALBEAM(IBEAM),    & ! input
              IBEAM, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,        & ! input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,            & ! input
              UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,            & ! input
              QUAD_WEIGHTS, QUAD_STRMWTS,                         & ! input
              INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,        & ! input
              T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,             & ! input
              MEAN_INTENSITY, FLUX_INTEGRAL,                      & ! Output
              DNMEAN_DIRECT,  DNFLUX_DIRECT )                       ! Output

          ENDIF

!  Finish this Beam solution, if only intensity is required

          IF ( .not. DO_SIMULATION_ONLY ) THEN

          IF (LPS_CHK) write(*,*) 'at complete jacs'

!  Profile weighting functions
!  ===========================

          IF ( DO_PROFILE_LINEARIZATION ) THEN

           IF (LPS_CHK) write(*,*) 'at complete lp bvp'

!  Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations

           DO LAYER_TO_VARY = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
             N_LAYER_WFS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

!  Regular BVP linearization

             IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

!  Get the linearized BVP solution for this beam component
!  Rob Fix, 01/06/14, Reorganized calling arguments, including TAYLOR_ORDER

              CALL LP_BVP_SOLUTION_MASTER                                   &
            ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_LOCALBEAM(IBEAM),     & ! Input
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_PLANE_PARALLEL,   & ! Input
              DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER_COMPONENT, IBEAM,  & ! Input, Flags and order
              NLAYERS, NSTREAMS, NSTREAMS_2, NTOTAL,                        & ! Input, Numbers
              N_SUBDIAG, N_SUPDIAG, LAYER_TO_VARY, N_LAYER_WFS,             & ! Input, Numbers
              DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
              SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,    & ! Input, Surface Stuff
              INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
              BANDMAT2, SMAT2, IPIVOT, SIPIVOT, LCONMASK, MCONMASK,         & ! Input, BVP Bandmat
              T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
              GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
              L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous  
              L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens/Thermal
              NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,         & ! Output, Linearized Constants + PI
              STATUS_SUB, MESSAGE, TRACE_1 )                                  ! Output

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from LP_BVP_SOLUTION_MASTER, '// &
                   'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
                STATUS = LIDORT_SERIOUS
                RETURN
              ENDIF

!  telescoped BVP linearization

             ELSE

!  Rob Fix, 01/06/14, Reorganized calling arguments, including TAYLOR_ORDER

              CALL LP_BVPTEL_SOLUTION_MASTER &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, FOURIER_COMPONENT, & ! Input, Flags, Order
             IBEAM, NLAYERS, NSTREAMS, NSTREAMS_2, LAYER_TO_VARY, N_LAYER_WFS,        & ! Input, Numbers
             N_SUPDIAG, N_SUBDIAG, ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,    & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,                  & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Input, Beam Quantities
             BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,                        & ! Input, BVProblem
             T_DELT_DISORDS, L_T_DELT_DISORDS,  LCON_XVEC, MCON_XVEC,       & ! Input, Extinction, nonactive
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                          & ! Input, Homogeneous + constants
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,  & ! Input, Greens Function
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                      & ! Input, Linearized Homogeneous
             L_ATERM_SAVE, L_BTERM_SAVE,                                    & ! Input, Linearized Greens
             L_WUPPER, L_WLOWER, NCON, PCON, NCON_XVEC, PCON_XVEC,          & ! Output
             STATUS_SUB, MESSAGE, TRACE_1 )                                   ! Output

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from LP_BVPTEL_SOLUTION_MASTER, '// &
                   'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
                STATUS = LIDORT_SERIOUS
                RETURN
              ENDIF

             ENDIF

!  Post-processing for Profile weighting functions (upwelling)
!  -----------------------------------------------------------

             IF (LPS_CHK) write(*,*) 'at complete lp post proc'

!mick fix 6/29/11 - initialize FLAGS_LP_GMULT
             FLAGS_LP_GMULT = .TRUE.

             IF ( DO_UPWELLING ) THEN

!  Linearized BOA source term

              CALL GET_LP_BOASOURCE                                           &
              ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,& ! Input
                DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,                         & ! Input
                DO_OBSERVATION_GEOMETRY,                                      & ! Input
                NLAYERS, NSTREAMS, N_USER_STREAMS, FOURIER_COMPONENT,         & ! Input
                IBEAM, LAYER_TO_VARY, N_LAYER_WFS, QUAD_STRMWTS, QUAD_WEIGHTS,& ! Input
                SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                  & ! Input
                USER_DIRECT_BEAM, DELTAU_SLANT, L_DELTAU_VERT,                & ! Input
                LCON, MCON, LCON_XVEC, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input
                T_WLOWER, L_XPOS, L_XNEG, L_T_WLOWER, L_WLOWER,               & ! Input
                NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_T_DELT_DISORDS,       & ! Input
                LP_BOA_MSSOURCE, LP_BOA_DBSOURCE, L_BOA_THTONLY_SOURCE )        ! Output

!  Linearized upwelling weighting functions

              CALL UPUSER_PROFILEWF &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES,                        & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,              & ! Input
              DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT,                    & ! Input
              DO_OBSERVATION_GEOMETRY, FLUX_MULTIPLIER,                 & ! Input
              NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,         & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! Input
              UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                  & ! Input
              FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS,     & ! Input
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
              U_XPOS, U_XNEG, U_WPOS1, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
              UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                    & ! Input
              PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,              & ! Input
              UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_UP,                    & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                & ! Input
              LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                  & ! Input
              L_T_DELT_USERM,   L_T_UTUP_USERM,                   & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTUP,                 & ! Input
              L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Input
              L_U_XPOS, L_U_XNEG, LP_U_WPOS1, NCON, PCON,                & ! Input
              L_XPOS, L_XNEG, L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC,  & ! Input
              L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                         & ! Input
              L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,              & ! Input
              L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP, L_UT_T_PARTIC, L_T_WUPPER, & ! Input
              L_BOA_THTONLY_SOURCE, LP_BOA_MSSOURCE, LP_BOA_DBSOURCE,        & ! Input
              FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,                & ! Output
              PROFILEWF_F, QUADPROFILEWF )                                     ! Output

             ENDIF

!  Post-processing for profile weighting functions (Downwelling)
!  ------------------------------------------------------------

             IF ( DO_DNWELLING ) THEN

!  TOA source term linearized

              CALL GET_LP_TOASOURCE &
                ( DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, IBEAM, N_LAYER_WFS, &
                  LP_TOA_SOURCE )

!  Downwelling weighting functions

              CALL DNUSER_PROFILEWF &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES,                         & ! Input
              DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,               & ! Input
              DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT,                     & ! Input
              DO_OBSERVATION_GEOMETRY, FLUX_MULTIPLIER,                  & ! Input
              NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                   & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
              UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                   & ! Input
              FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS,      & ! Input
              TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,     & ! Input
              QUAD_STREAMS, USER_SECANTS, DO_LAYER_SCATTERING,           & ! Input
              LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,           & ! Input
              T_DELT_MUBAR,   T_UTDN_MUBAR,                      & ! Input
              T_DELT_USERM,   T_UTDN_USERM,                      & ! Input
              T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN,      & ! Input
              T_DELT_DISORDS, T_DISORDS_UTDN,                    & ! Input
              GAMMA_M, GAMMA_P, SIGMA_M, ATERM_SAVE, BTERM_SAVE, & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WLOWER,          & ! Input
              U_XPOS, U_XNEG, U_WNEG1, HMULT_1, HMULT_2, EMULT_DN,       & ! Input
              UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN,                    & ! Input
              PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,              & ! Input
              UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_DN,                    & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                & ! Input
              LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR,                  & ! Input
              L_T_DELT_USERM,   L_T_UTDN_USERM,                   & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTDN,                 & ! Input
              L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Input
              L_U_XPOS, L_U_XNEG, LP_U_WNEG1, NCON, PCON,                & ! Input
              L_XPOS, L_XNEG, L_WLOWER, NCON_XVEC, PCON_XVEC,            & ! Input
              L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                         & ! Input
              L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,              & ! Input
              L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN, L_UT_T_PARTIC,         & ! Input
              L_T_WLOWER, LP_TOA_SOURCE,                                 & ! Input
              FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,            & ! Output
              PROFILEWF_F, QUADPROFILEWF )                                 ! Output

             ENDIF

!  Post-processing for Profile weighting functions (Mean/Flux)
!  -----------------------------------------------------------

             IF ( DO_INCLUDE_MVOUTPUT ) THEN

               IF (LPS_CHK) write(*,*) 'at complete lp mean val jacs'

               CALL MIFLUX_PROFILEWF                                      &
               ( DO_UPWELLING, DO_DNWELLING, DO_LOCALBEAM(IBEAM),         & ! Input
                 IBEAM, LAYER_TO_VARY, N_LAYER_WFS,                       & ! Input
                 NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,                    & ! Input
                 PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 & ! Input
                 UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                 & ! Input
                 QUAD_WEIGHTS, QUAD_STRMWTS,                              & ! Input
                 INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,             & ! Input
                 T_DELT_MUBAR, T_UTDN_MUBAR, QUADPROFILEWF,               & ! Input
                 LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,      & ! Input
                 MINT_PROFILEWF, FLUX_PROFILEWF )                           ! Output

!mick debug
!write(*,*) 'leaving MIFLUX_PROFILEWF'

             ENDIF

!  Finish loop over layers with variation

            ENDIF
           ENDDO

!  End atmospheric profile weighting functions

          ENDIF

!   Surface Property Jacobians (only if surface included)
!  ===========================

          IF ( DO_SURFACE_LINEARIZATION .AND. DO_INCLUDE_SURFACEWF ) THEN

           IF (LPS_CHK) write(*,*) 'at complete surf jacs'

           IF ( N_SURFACE_WFS > 0 ) THEN

!  Surface property WF master

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add argument to SURFACEWF_MASTER

             CALL SURFACEWF_MASTER                                           &
             ( DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE,                  & ! Input
               DO_USER_STREAMS, DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY,  & ! Input
               DO_LOCALBEAM(IBEAM), DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT,   & ! Input
               DO_MSMODE_THERMAL, DO_OBSERVATION_GEOMETRY,                   & ! Input
               NLAYERS, NTOTAL, N_SUPDIAG, N_SUBDIAG, N_USER_LEVELS,         & ! Input
               NSTREAMS, NSTREAMS_2, N_USER_STREAMS, N_SURFACE_WFS,          & ! Input
               PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
               UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                       & ! Input
               FOURIER_COMPONENT, IBEAM, FLUX_MULTIPLIER,                    & ! Input
               SURFACE_BB_INPUT, LS_EMISSIVITY, LS_USER_EMISSIVITY,          & ! Input
               SURFACE_FACTOR, ALBEDO, LS_BRDF_F, LS_BRDF_F_0,               & ! Input
               USER_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0,                & ! Input
               QUAD_WEIGHTS, QUAD_STRMWTS, ATMOS_ATTN, IDOWNSURF,            & ! Input
               LCONMASK, MCONMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,         & ! Input
               LCON, MCON, XPOS, XNEG, H_XPOS, H_XNEG, H_WLOWER,             & ! Input
               T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                     & ! Input
               T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                     & ! Input
               T_DELT_DISORDS, T_DISORDS_UTUP, U_XPOS, U_XNEG, HMULT_1,      & ! Input
               HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,  & ! Input
               SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,                  & ! Output
               STATUS_SUB, MESSAGE, TRACE_1 )                                  ! Output

!  Error handling

             IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER_COMPONENT
               TRACE_2 = 'Error return from SURFACEWF_MASTER, '// &
                  'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS
               RETURN
             ENDIF

!  End surface-property linearizations

           ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@ Addition of SLEAVE weighting function @@@@@@@@@@@@@@@@@@@
!    R. Spurr, 22 August 2012
!
           IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then

             IF (LPS_CHK) write(*,*) 'at complete sleave jacs'

!  Direct-beam Flag (repeated here, just so we know!)

             DO_INCLUDE_DIRECTBEAM = &
               ( DO_UPWELLING .AND. DO_LOCALBEAM(IBEAM) ) &
                 .AND. .NOT.DO_MSMODE_LIDORT

!  Rob Fix @@@ 11 Sep 12, Add New Line to the Argument list

              CALL LIDORT_LSSL_WFS ( &
                DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT,                  &
                DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING,                 &
                DO_SL_ISOTROPIC, DO_OBSERVATION_GEOMETRY,                    &
                N_SLEAVE_WFS, N_SURFACE_WFS, NSTREAMS, NLAYERS,              &
                N_USER_STREAMS, FOURIER_COMPONENT, IBEAM,                    &
                N_USER_LEVELS, N_DIRECTIONS, WHICH_DIRECTIONS,               &
                NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTREAMS_2,                    &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                     &
                UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
                DO_BRDF_SURFACE, ALBEDO, USER_BRDF_F,                        &
                SURFACE_FACTOR, FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS, &
                LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM,            &
                XPOS, XNEG,                                         &
                T_DELT_EIGEN, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,     &
                T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,           &
                T_UTUP_EIGEN, T_UTDN_EIGEN,  HMULT_1, HMULT_2,      &
                U_XPOS, U_XNEG,                                     &
                UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD, &
                SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,        &
                STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

             IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER_COMPONENT
               TRACE_2 = 'Error return from LIDORT_LSSL_WFS, '// &
                'Beam Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS
               RETURN
             ENDIF

!  End surface-leaving linearizations

           ENDIF

!@@@@@@@@@@ END Addition of SLEAVE weighting function @@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  End of surface weighting functions

          ENDIF

!  New, 18 March 2014. Linearization for BLACKBODY. Only if IB = 1, FOURIER = 0

           IF ( ( DO_ATMOS_LBBF .or. DO_SURFACE_LBBF ) .and. IBEAM.eq.1 .and. Fourier_component.eq.0 )  THEN

             IF (LPS_CHK) write(*,*) 'at complete lbbf jacs'

              if ( n_partlayers .gt. 0 ) then
                 call lidort_lbbf_jacobians_wpartials &
                  ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,       & ! Inputs 3/25
                    DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,               & ! Inputs
                    DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,    & ! input
                    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                        & ! input
                    NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,           & ! input
                    NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,                   & ! input
                    N_PARTLAYERS, PARTLAYERS_LAYERIDX,                          & ! Input 3/21
                    PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                    & ! Input 3/21
                    UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                     & ! Input
                    USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,       & ! Input
                    QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                   & ! input
                    SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                & ! input
                    EMISSIVITY, USER_EMISSIVITY,                                & ! input
                    FLUX_MULTIPLIER, OMEGA_MOMS(:,0), DELTAU_VERT, PARTAU_VERT, & ! Input 3/21
                    T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,             & ! inputs 3/25
                    KEIGEN, TTERM_SAVE, XPOS, XNEG,                             & ! inputs 3/21
                    T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                   & ! inputs 3/21
                    BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                           & ! Input
                    T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                   & ! Inputs 3/21
                    U_XPOS, U_XNEG, HMULT_1, HMULT_2,                           & ! inputs
                    UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,         & ! Inputs 3/21
                    ABBWFS_JACOBIANS, ABBWFS_FLUXES,                            & ! Output
                    SBBWFS_JACOBIANS, SBBWFS_FLUXES,                            & ! Output
                    STATUS_SUB, MESSAGE, TRACE_1 )                                ! Output
              else
                 call lidort_lbbf_jacobians_whole &
                  ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,    & ! Inputs 3/25
                    DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,            & ! Inputs
                    DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, & ! input
                    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                     & ! input
                    NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,        & ! input
                    NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,                & ! input
                    UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                  & ! Input
                    USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! Input
                    QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                & ! input
                    SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,             & ! input
                    EMISSIVITY, USER_EMISSIVITY,                             & ! inputs
                    FLUX_MULTIPLIER, OMEGA_MOMS(:,0), DELTAU_VERT,           & ! Input
                    T_DELT_DISORDS, T_DELT_EIGEN, KEIGEN, XPOS, XNEG,        & ! inputs
                    TTERM_SAVE, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Input
                    T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2,          & ! inputs
                    ABBWFS_JACOBIANS, ABBWFS_FLUXES,                         & ! Output
                    SBBWFS_JACOBIANS, SBBWFS_FLUXES,                         & ! Output
                    STATUS_SUB, MESSAGE, TRACE_1 )                             ! Output
              endif

!  Error handling

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                 write(CF,'(I2)')FOURIER_COMPONENT
                 TRACE_2 = 'Error return from lidort_lbbf_jacobians, With Solar Term '// &
                           'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
                 STATUS = LIDORT_SERIOUS
                 RETURN
               ENDIF

!  ENd of BB Jacobians

            ENDIF

!  End complete weighting function clause (NON simulation-only!)

          ENDIF

!  End loop over beam solutions

        END IF
      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE LIDORT_LPS_FOURIER_MASTER

end MODULE LIDORT_LPS_masters
