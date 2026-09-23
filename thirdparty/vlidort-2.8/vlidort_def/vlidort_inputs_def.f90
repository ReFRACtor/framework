
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

      MODULE VLIDORT_Inputs_def_m

!  This module contains the following VLIDORT input structures,
!  with intents :

!        VLIDORT_Fixed_Boolean    nested in VLIDORT_Fixed_Inputs
!        VLIDORT_Fixed_Control    nested in VLIDORT_Fixed_Inputs
!        VLIDORT_Fixed_Sunrays    nested in VLIDORT_Fixed_Inputs
!     VLIDORT_Fixed_UserValues    nested in VLIDORT_Fixed_Inputs
!        VLIDORT_Fixed_Chapman    nested in VLIDORT_Fixed_Inputs
!        VLIDORT_Fixed_Optical    nested in VLIDORT_Fixed_Inputs
!          VLIDORT_Fixed_Write    nested in VLIDORT_Fixed_Inputs
!         VLIDORT_Fixed_Inputs    Intent(In)

!     VLIDORT_Modified_Boolean    nested in VLIDORT_Modified_Inputs
!     VLIDORT_Modified_Control    nested in VLIDORT_Modified_Inputs
!     VLIDORT_Modified_Sunrays    nested in VLIDORT_Modified_Inputs
!  VLIDORT_Modified_UserValues    nested in VLIDORT_Modified_Inputs
!     VLIDORT_Modified_Chapman    nested in VLIDORT_Modified_Inputs
!     VLIDORT_Modified_Optical    nested in VLIDORT_Modified_Inputs
!      VLIDORT_Modified_Inputs    Intent(InOut)

      USE VLIDORT_PARS_m, Only : fpk, MAXMOMENTS_INPUT, MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_GEOMETRIES, &
                                 MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,    &
                                 MAX_USER_LEVELS, MAXSTOKES_SQ

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Boolean

!  Full Radiance calculation, SS + MS fields
!    if False, then VLIDORT only does a multiple scatter calculation.

      LOGICAL     :: TS_DO_FULLRAD_MODE

!  Surface and thermal emission flags

      LOGICAL     :: TS_DO_THERMAL_EMISSION
      LOGICAL     :: TS_DO_SURFACE_EMISSION

!  Beam particular solution pseudo-spherical options

      LOGICAL     :: TS_DO_PLANE_PARALLEL

!  Surface control (New, 23 March 2010)

!      LOGICAL     :: TS_DO_BRDF_SURFACE

!  Directional control

      LOGICAL     :: TS_DO_UPWELLING
      LOGICAL     :: TS_DO_DNWELLING

!  Stream angle flag. removed for Version 2.8, 7/6/16
!      LOGICAL     :: TS_DO_QUAD_OUTPUT

!  Contributions (RT Solutions Use Only)

      LOGICAL     :: TS_DO_TOA_CONTRIBS

!  Surface

      LOGICAL     :: TS_DO_LAMBERTIAN_SURFACE

!  Special Options (RT Solutions Use Only)

      LOGICAL     :: TS_DO_SPECIALIST_OPTION_1
      LOGICAL     :: TS_DO_SPECIALIST_OPTION_2
      LOGICAL     :: TS_DO_SPECIALIST_OPTION_3

!  Surface leaving Control. New 17 May 2012

      LOGICAL     :: TS_DO_SURFACE_LEAVING
      LOGICAL     :: TS_DO_SL_ISOTROPIC

!   Water leaving flag added 28 October 2015 (2.7a)
!   Fluorescence  flag added 31 January 2016  (2.8)

      LOGICAL     :: TS_DO_WATER_LEAVING
      LOGICAL     :: TS_DO_FLUORESCENCE

!   Water-leaving transmittance iteration flag added 7/6/16 (2.8)

      LOGICAL     :: TS_DO_TF_ITERATION 

!  Water-leaving output flag. 3/18/19 for Version 2.8.1
      
      LOGICAL     :: TS_DO_WLADJUSTED_OUTPUT 

!  TOA/BOA Illumination flags. 3/23/19 for Version 2.8.1

      LOGICAL     :: TS_DO_TOA_ILLUMINATION
      LOGICAL     :: TS_DO_BOA_ILLUMINATION

!   4/26/19 Added control for the media problem. Version 3.8a
!     Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!       1 = Isotropic illumination from Top, 2 = Isotropic illumination from BOA

      LOGICAL     :: TS_DO_ALBTRN_MEDIA(2)

!  Flag for the Planetary problem calculation.
!    -- This is independent of the previous flags !!!!!

      LOGICAL     :: TS_DO_PLANETARY_PROBLEM

!  1/31/21. Version 2.8.3. Flag for calculating MSSTs output
!    -- Additional output of the multiple-scattering source terms (MSSTs)
!    -- This must be set by hand, it is not a configuration file read

      LOGICAL     :: TS_DO_MSSTS

!  1/31/21. Version 2.8.3. Flag for using an NSTOKES = 2 calculation for Fourier 0
!    -- This must be set by hand, it is not a configuration file read

      LOGICAL     :: TS_DO_FOURIER0_NSTOKES2

      END TYPE VLIDORT_Fixed_Boolean

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Control

!  Taylor ordering parameter. Should be set to 2 or 3
!     Added, 2/19/14 for Taylor-series expansions

      INTEGER   :: TS_TAYLOR_ORDER

!  Number of Stokes parameters

      INTEGER   :: TS_NSTOKES

!  Number of discrete ordinate streams

      INTEGER   :: TS_NSTREAMS

!  Number of computational layers

      INTEGER   :: TS_NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing spherical correction algorithm)

      INTEGER   :: TS_NFINELAYERS

!  Number of thermal coefficients (2 should be the default)

      INTEGER   :: TS_N_THERMAL_COEFFS

!  Accuracy for convergence of Fourier series

      REAL(fpk) :: TS_VLIDORT_ACCURACY

!  1/31/21. Version 2.8.3. ASYMTX Tolerance variable
!    -- This must be set by hand, it is not a configuration file read

      REAL(fpk) :: TS_ASYMTX_TOLERANCE
      
!  Special Options (RT Solutions Use Only)

      INTEGER   :: TS_NLAYERS_NOMS
      INTEGER   :: TS_NLAYERS_CUTOFF

!  Water-leaving: Control for iterative calculation of transmittanaces
!    Variables added for Version 2.8, 7/6/16

      INTEGER   :: TS_TF_MAXITER
      REAL(fpk) :: TS_TF_CRITERION

!  TOA/BOA Isotropic Illumination. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized
      
      REAL(fpk) :: TS_TOA_ILLUMINATION
      REAL(fpk) :: TS_BOA_ILLUMINATION

      END TYPE VLIDORT_Fixed_Control

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Sunrays

!  Flux factor ( should be 1 or pi ). Same for all solar beams.

      REAL(fpk)  :: TS_FLUX_FACTOR

      END TYPE VLIDORT_Fixed_Sunrays

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_UserValues

!  User-defined zenith angle input
!    Now moved to Modified User Values, 25 October 2012
!      INTEGER                                 :: TS_N_USER_STREAMS

!  User-defined vertical level output

      INTEGER                                 :: TS_N_USER_LEVELS

      END TYPE VLIDORT_Fixed_UserValues

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Chapman

!  Multilayer Height inputs in [km]
!   Required for the Chapman function calculations

      REAL(fpk), dimension (0:MAXLAYERS) :: TS_HEIGHT_GRID

!  Multilayer atmospheric inputs. Pressures in [mb], temperatures in [K]
!   Required for the Chapman function calculations, refractive geometry

      REAL(fpk), dimension (0:MAXLAYERS) :: TS_PRESSURE_GRID
      REAL(fpk), dimension (0:MAXLAYERS) :: TS_TEMPERATURE_GRID

!  Number of fine layer gradations
!   Required for the Chapman function calculations, refractive geometry

      INTEGER,   dimension (MAXLAYERS)   :: TS_FINEGRID

!  Refractive index parameter
!    (only for Chapman function calculations with refractive index)

      REAL(fpk) :: TS_RFINDEX_PARAMETER

      END TYPE VLIDORT_Fixed_Chapman

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Optical

!  Multilayer optical property (bulk)

      REAL(fpk), dimension ( MAXLAYERS ) :: TS_DELTAU_VERT_INPUT

!  Phase function Legendre-polynomial expansion coefficients
!    Include all that you require for exact single scatter calculations

      REAL(fpk), dimension ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ ) :: TS_GREEKMAT_TOTAL_INPUT

!  F-matrix input as an alternative to the use of expansion coefficients in
!   the single-scatter (SSCORR and FO) codes.
!  Introduced for Version 2.8 by R. Spurr, 02/08/16
!mick mod 9/19/2017 - changed names from "FMATRIX_INPUT_UP/DN" to "FMATRIX_UP/DN"
!                     for consistency with the companion linearized input variables

      REAL(fpk), dimension ( MAXLAYERS, MAX_GEOMETRIES,  6 ) :: TS_FMATRIX_UP
      REAL(fpk), dimension ( MAXLAYERS, MAX_GEOMETRIES,  6 ) :: TS_FMATRIX_DN

!  Lambertian Surface control

      REAL(fpk) :: TS_LAMBERTIAN_ALBEDO

!  Thermal Black Body functions

      REAL(fpk), dimension ( 0:MAXLAYERS ) :: TS_THERMAL_BB_INPUT

!  Surface Black body inputs

      REAL(fpk) :: TS_SURFACE_BB_INPUT

!  Special LTE variables (RT Solutions Use Only).
!   This has been superseded in Version 2.7. No longer required
!      REAL(fpk), dimension ( 2, MAXLAYERS ) :: TS_LTE_DELTAU_VERT_INPUT
!      REAL(fpk), dimension ( 0:MAXLAYERS )  :: TS_LTE_THERMAL_BB_INPUT

!  Rob Fix 3/18/15. Add Wavelength (Microns) as a Diagnostic
!    -- This can now be read from configuration file...(Version 2.8)
!    -- This should always be set, even if not used
!    -- This MUST BE SET when using VBRDF and/or VSLEAVE supplements with
!         wavelength-dependent inputs

      REAL(fpk) :: TS_ATMOS_WAVELENGTH

      END TYPE VLIDORT_Fixed_Optical

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Write

!  Debug control

      LOGICAL ::            TS_DO_DEBUG_WRITE

!  Input file

      LOGICAL ::            TS_DO_WRITE_INPUT
      CHARACTER (LEN=60) :: TS_INPUT_WRITE_FILENAME

!  Scene file

      LOGICAL ::            TS_DO_WRITE_SCENARIO
      CHARACTER (LEN=60) :: TS_SCENARIO_WRITE_FILENAME

!  Fourier component results file

      LOGICAL ::            TS_DO_WRITE_FOURIER
      CHARACTER (LEN=60) :: TS_FOURIER_WRITE_FILENAME

!  Results file

      LOGICAL ::            TS_DO_WRITE_RESULTS
      CHARACTER (LEN=60) :: TS_RESULTS_WRITE_FILENAME

      END TYPE VLIDORT_Fixed_Write

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Inputs

      TYPE(VLIDORT_Fixed_Boolean)    :: Bool
      TYPE(VLIDORT_Fixed_Control)    :: Cont
      TYPE(VLIDORT_Fixed_Sunrays)    :: Sunrays
      TYPE(VLIDORT_Fixed_UserValues) :: UserVal
      TYPE(VLIDORT_Fixed_Chapman)    :: Chapman
      TYPE(VLIDORT_Fixed_Optical)    :: Optical
      TYPE(VLIDORT_Fixed_Write)      :: Write

      END TYPE VLIDORT_Fixed_Inputs

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Boolean

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 2.8. 9/18/16, 3/1/17, 9/19/17.
!mick mod 9/19/2017 - reordered FO variables in a manner similar to LIDORT3.8

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      LOGICAL     :: TS_DO_FOCORR

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first.

      LOGICAL     :: TS_DO_FOCORR_EXTERNAL

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.
!mick mod 9/19/2017 - now defined in VLIDORT internally

      !LOGICAL     :: TS_DO_FOCORR_ALONE

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first. DO_FOCORR_EXTERNAL should be off.

      LOGICAL     :: TS_DO_FOCORR_NADIR
      LOGICAL     :: TS_DO_FOCORR_OUTGOING

!  Additional SSCORR flags
!  -----------------------

!  Flag for performing SSCORR truncation. Single-scattering
!     - Kept here, but disabled for Version 2.8.

      LOGICAL     :: TS_DO_SSCORR_TRUNCATION

!  Flag for using the F-matrix in the Single-scatter calculations (instead of Greekmat Coefficients)
!     - Introduced for Version 2.8, 7/7/16.  R. Spurr
!     - DO_FOCORR must be set first.
!     - Does not apply to thermal case..

      LOGICAL     :: TS_DO_SSCORR_USEFMAT

!  Other variables
!  ---------------

 !  Additional Control for Externalized Water-leaving inputs.
!     Introduced 3/18/19 for Version 2.8.1
     
      LOGICAL     :: TS_DO_EXTERNAL_WLEAVE

!  Double convergence test flag

      LOGICAL     :: TS_DO_DOUBLE_CONVTEST

!  Basic top-level solar beam control

      LOGICAL     :: TS_DO_SOLAR_SOURCES

!  1/31/21. Version 2.8.3. Add the classical solution flag
!     ==> Could be modified, as the Greens function solution applies only for SOLAR SOURCES alone.

      LOGICAL     :: TS_DO_CLASSICAL_SOLUTION

!  Pseudo-spherical input control

      LOGICAL     :: TS_DO_REFRACTIVE_GEOMETRY
      LOGICAL     :: TS_DO_CHAPMAN_FUNCTION

!  Scatterers and phase function control
!      Isotropic, no_azimuth and All-Fourier flags disabled.

      LOGICAL     :: TS_DO_RAYLEIGH_ONLY
!      LOGICAL    :: TS_DO_ISOTROPIC_ONLY
!      LOGICAL    :: TS_DO_NO_AZIMUTH
!      LOGICAL    :: TS_DO_ALL_FOURIER

!  Delta-M scaling flag

      LOGICAL     :: TS_DO_DELTAM_SCALING

!  2 new flags in Version 3.0

      LOGICAL     :: TS_DO_SOLUTION_SAVING
      LOGICAL     :: TS_DO_BVP_TELESCOPING

!  Stream angle flags. Flag renamed, Version 2.5

!      LOGICAL    :: TS_DO_USER_STREAMS
      LOGICAL     :: TS_DO_USER_VZANGLES

!  Mean value control

      LOGICAL     :: TS_DO_ADDITIONAL_MVOUT
      LOGICAL     :: TS_DO_MVOUT_ONLY

!  Transmittance only for thermal mode.

      LOGICAL     :: TS_DO_THERMAL_TRANSONLY

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      LOGICAL     :: TS_DO_OBSERVATION_GEOMETRY

!  1/31/21. Version 2.8.3. Doublet-Geometry input control.

      LOGICAL     :: TS_DO_DOUBLET_GEOMETRY

      END TYPE VLIDORT_Modified_Boolean

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Control

!  Number of Legendre phase function expansion moments
!       May be re-set after Checking

      INTEGER   :: TS_NGREEK_MOMENTS_INPUT

      END TYPE VLIDORT_Modified_Control

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Sunrays

!  Number of solar beams to be processed

!      INTEGER    :: TS_NBEAMS
      INTEGER    :: TS_N_SZANGLES

!  Bottom-of-atmosphere solar zenith angles, DEGREES

!      REAL(fpk), dimension (MAXBEAMS) :: TS_BEAM_SZAS
      REAL(fpk), dimension (MAX_SZANGLES) :: TS_SZANGLES

      END TYPE VLIDORT_Modified_Sunrays

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_UserValues

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER                                 :: TS_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: TS_USER_RELAZMS

!  User-defined zenith angle input. Quantities renamed, Version 2.5

!      INTEGER                                  :: TS_N_USER_STREAMS
!      REAL(fpk), dimension (MAX_USER_STREAMS)  :: TS_USER_ANGLES_INPUT
      INTEGER                                  :: TS_N_USER_VZANGLES
      REAL(fpk), dimension (MAX_USER_VZANGLES) :: TS_USER_VZANGLES_INPUT

!  User-defined vertical level output. From Top-of-atmosphere. Example for 4 outputs:
!     USER_LEVELS(1) = 0.0           --> Top-of-atmosphere
!     USER_LEVELS(2) = 1.1           --> One tenth of the way down into Layer 2
!     USER_LEVELS(3) = 3.5           --> One half  of the way down into Layer 4
!     USER_LEVELS(4) = dble(NLAYERS) --> Bottom of atmosphere

      REAL(fpk), dimension (MAX_USER_LEVELS)  :: TS_USER_LEVELS

!  Geometry specification height

      REAL(fpk)                               :: TS_GEOMETRY_SPECHEIGHT

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      INTEGER                                 :: TS_N_USER_OBSGEOMS

!  User-defined Observation Geometry angle input
!     New variable, 25 October 2012, for Observational Geometry input

      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: TS_USER_OBSGEOMS_INPUT

!  Doublet-Geometry input control and angles
!    -- 1/31/21. Version 2.8.3. (2/16/21). Newly introduced

      INTEGER                                    :: TS_N_USER_DOUBLETS
      REAL(fpk), dimension (MAX_USER_VZANGLES,2) :: TS_USER_DOUBLETS

      END TYPE VLIDORT_Modified_UserValues

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Chapman

!  Output from Chapman function calculations - can also be input (not used)

      !REAL(fpk), dimension (MAXLAYERS,MAXLAYERS,MAXBEAMS) :: &
      !  TS_CHAPMAN_FACTORS

!  Earth radius in [km] for Chapman function calculation of TAUTHICK_INPUT

      REAL(fpk) :: TS_EARTH_RADIUS

      END TYPE VLIDORT_Modified_Chapman

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Optical

!  Multilayer optical property (bulk)

      REAL(fpk), dimension ( MAXLAYERS ) :: TS_OMEGA_TOTAL_INPUT

      END TYPE VLIDORT_Modified_Optical

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Inputs

      TYPE(VLIDORT_Modified_Boolean)    :: MBool
      TYPE(VLIDORT_Modified_Control)    :: MCont
      TYPE(VLIDORT_Modified_Sunrays)    :: MSunrays
      TYPE(VLIDORT_Modified_UserValues) :: MUserVal
      TYPE(VLIDORT_Modified_Chapman)    :: MChapman
      TYPE(VLIDORT_Modified_Optical)    :: MOptical

      END TYPE VLIDORT_Modified_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_Fixed_Boolean, &
                VLIDORT_Fixed_Control, &
                VLIDORT_Fixed_Sunrays, &
                VLIDORT_Fixed_UserValues, &
                VLIDORT_Fixed_Chapman, &
                VLIDORT_Fixed_Optical, &
                VLIDORT_Fixed_Write, &
                VLIDORT_Fixed_Inputs, &
                VLIDORT_Modified_Boolean, &
                VLIDORT_Modified_Control, &
                VLIDORT_Modified_Sunrays, &
                VLIDORT_Modified_UserValues, &
                VLIDORT_Modified_Chapman, &
                VLIDORT_Modified_Optical, &
                VLIDORT_Modified_Inputs

      END MODULE VLIDORT_Inputs_def_m

