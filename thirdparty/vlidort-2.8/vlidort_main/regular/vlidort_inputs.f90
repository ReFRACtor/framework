
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
! #            VLIDORT_INPUT_MASTER (master), calls             #
! #             VLIDORT_INIT_INPUTS                             #
! #             VLIDORT_READ_INPUTS                             #
! #                                                             #
! #            VLIDORT_Sup_Init                                 #
! #            VLIDORT_BRDF_Sup_Init                            #
! #            VLIDORT_SLEAVE_Sup_Init                          #
! #            VLIDORT_SS_Sup_Init                              #
! #                                                             #
! #    These routines are called by the Main VLIDORT module     #
! #                                                             #
! #            VLIDORT_CHECK_INPUT_DIMS                         #
! #            VLIDORT_CHECK_INPUT                              #
! #            VLIDORT_CHECK_INPUT_OPTICAL                      #
! #            VLIDORT_DERIVE_INPUT                             #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. SEVERAL INNOVATIONS and CHANGES
!    -- Converge routines have been moved to their own module (vlidort_converge.f90)
!    -- Need additional Taylor routines to be used, for Green's function post-processing
!    -- MSST option is now included, preserves output for sphericity correction
!    -- DO_DOUBLET_GEOMETRY post-processing operation is in force.
!    -- Whole/part Sourceterm routines completely rewritten (Green's function included)
!    -- QuadIntens subroutines completely rewritten (Green's function included)
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier

!  1/31/21. Version 2.8.3. This Module
!    -- Initialize, read/copy 3 new Boolean flags (DOUBLET_GEOMETRY, MSSTS, CLASSICAL_SOLUTION)
!    -- Initialize and read Outgoing sphericity criticality flag and attenuation criterion
!    -- Apply appropriate checks on these variables (e.g. nstokes=1/3 for Green's function)

      MODULE vlidort_inputs_m

      !PRIVATE
      PUBLIC

      CONTAINS

      SUBROUTINE VLIDORT_INPUT_MASTER ( &
        FILNAM,             & ! INPUT
        VLIDORT_FixIn,      & ! OUTPUTS
        VLIDORT_ModIn,      & ! OUTPUTS
        VLIDORT_InputStatus ) ! OUTPUTS

!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination
        
!  Parameter types

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, &
                                 MAX_USER_LEVELS, MAX_USER_OBSGEOMS, MAXLAYERS, MAX_MESSAGES, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_INUNIT

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_Outputs_def_m

      IMPLICIT NONE

!  Inputs
!  ------

      CHARACTER (LEN=*), intent(in) :: FILNAM

!  Outputs
!  -------

      TYPE(VLIDORT_Fixed_Inputs)            , INTENT (OUT) :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs)         , INTENT (OUT) :: VLIDORT_ModIn
      TYPE(VLIDORT_Input_Exception_Handling), INTENT (OUT) :: VLIDORT_InputStatus

!  Local variables
!  ===============

      INTEGER       :: FILUNIT 
      INTEGER       :: ISTAT, STATUS_SUB 

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER       :: STATUS
      INTEGER       :: NMESSAGES
      CHARACTER*120 :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*120 :: ACTIONS (0:MAX_MESSAGES)

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Initialize variables
!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added (1/31/16)
!     Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.
!     Version 2.8: FO/SS flags upgraded. 9/19/16, 3/1/17. Argument list rearranged.
!                  VLIDORT input type structures now used for argument passing - 9/19/2017
!                  All VLIDORT fixed and modified inputs now initialized - 9/19/2017
!     Version 2.8.3 (1/31/21). Initialize 3 new Boolean flags

      CALL VLIDORT_INIT_INPUTS ( &
         VLIDORT_FixIn, VLIDORT_ModIn ) !Outputs

!  Open VLIDORT config file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,IOSTAT=ISTAT,STATUS='OLD')

!  Open file error

      IF (ISTAT .GT. 0) THEN

!  Define the exception

         STATUS = VLIDORT_SERIOUS
         NMESSAGES = NMESSAGES + 1
         MESSAGES(NMESSAGES) = 'Openfile failure for ' // Trim(FILNAM)
         ACTIONS(NMESSAGES)  = 'Find the Right File!!'

!  Copy to exception handling

         VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
         VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
         VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
         VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS

         RETURN
      ENDIF

!  Read standard inputs

      CALL VLIDORT_READ_INPUTS ( &
         VLIDORT_FixIn, VLIDORT_ModIn,& !InOut
         STATUS_SUB,                  & !Output
         NMESSAGES, MESSAGES, ACTIONS ) !InOut

!  Read file error

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) STATUS = VLIDORT_SERIOUS

!  Copy to exception handling

      VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
      VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
      VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
      VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS
      CLOSE(FILUNIT)

!  Return

      RETURN

!  Finish

      END SUBROUTINE VLIDORT_input_master

!

      SUBROUTINE VLIDORT_INIT_INPUTS &
      ( VLIDORT_FixIn, VLIDORT_ModIn ) !Outputs

!  Initialize variables

!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added 1/31/16).
!     Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.
!     Version 2.8: FO/SS flags upgraded. 9/19/16, 3/1/17. Argument list rearranged.
!                  VLIDORT input type structures now used for argument passing - 9/19/2017
!                  All VLIDORT fixed and modified inputs now initialized - 9/19/2017

!  Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!    -- Introduce flag for including TOA isotropic illumination
!    -- Also include value of this illumination
   
!  1/31/21. Version 2.8.3. This Subroutine.
!    -- Initialize 3 new Boolean flags (DOUBLET_GEOMETRY, MSSTS, CLASSICAL_SOLUTION)

!  Initialises all inputs for VLIDORT
!  ---------------------------------

!  Module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
                                 MAX_USER_LEVELS, MAXBEAMS, ZERO
      USE VLIDORT_Inputs_def_m

!  Implicit none

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT standard input structures

      TYPE(VLIDORT_Fixed_Inputs), INTENT(OUT)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(OUT) :: VLIDORT_ModIn

!  Initialize inputs
!  =================

!  VLIDORT Fixed Boolean
!  ---------------------

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE        = .FALSE.

      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION    = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION    = .FALSE.

      VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL      = .FALSE.

!  Surface control (New, 23 March 2010).  removed for Version 2.8, 7/6/16
      !VLIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE        = .FALSE.

      VLIDORT_FixIn%Bool%TS_DO_UPWELLING           = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_DNWELLING           = .FALSE.

!  Stream angle flag. removed for Version 2.8, 7/6/16

      !VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT         = .FALSE.

!  Contributions (RT Solutions Use Only)

      VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS        = .FALSE.

      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = .FALSE.

!  1/31/21. Version 2.8.3. (First introduced 11/20/19).
!    -- Initialize the DO_MSSTS flag. This is not a configuration file-read.

      VLIDORT_FixIn%Bool%TS_DO_MSSTS               = .FALSE.

!  1/31/21. Version 2.8.3. (RTS 2/16/21).
!    -- Initialize the DO_FOURIER0_NSTOKES2 flag. This is not a configuration file-read.

      VLIDORT_FixIn%Bool%TS_DO_FOURIER0_NSTOKES2   = .FALSE.

!  Special Options (RT Solutions Use Only). These are not configuration file-reads.

      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1 = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2 = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3 = .FALSE.

!  Surface leaving Control. New 17 May 2012

      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING     = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC        = .FALSE.

!  Water leaving flag added 28 October 2015 (2.7a)
!  Fluorescence  flag added 31 January 2016  (2.8)

      VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING       = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_FLUORESCENCE        = .FALSE.

!  Water-leaving transmittance iteration flag added 7/6/16 (2.8)

      VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION        = .FALSE.

!  Water-leaving output flag. 3/18/19 for Version 2.8.1
      
      VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT   = .FALSE.

!  Planetary problem, and ALBTRN_MEDIA
!   4/26/19 Added control for the media problem. Version 2.8.1

      VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA         = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM         = .FALSE.

!  TOA/BOA Illumination flags. 3/23/19 for Version 2.8.1

      VLIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION         = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION         = .FALSE.

!  VLIDORT Modified Boolean
!  ------------------------

!  Flags for FO and SS Corr. Newly re-arranged, Version 3.8, 3/3/17

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

!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined internally

      VLIDORT_ModIn%MBool%TS_DO_FOCORR               = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL      = .FALSE.
      !VLIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE         = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR         = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING      = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION    = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT       = .FALSE.

!  1/31/21. Version 2.8.3.. Doublet-Geometry input control, initialize

      VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY     = .FALSE.

!  Additional Control for Externalized water-leaving inputs. 3/18/19 for Version 2.8.1

      VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE      = .FALSE.

!  1/31/21. Version 2.8.3. Initialize new variable for RT Solution method
!            TRUE = Classical solution of PI, FALSE = Green's function solution

      VLIDORT_ModIn%MBool%TS_DO_CLASSICAL_SOLUTION   = .FALSE.

!  Other flags

      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST      = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES        = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY  = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION     = .FALSE.

!     Isotropic, no_azimuth and All-Fourier flags disabled.

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY        = .FALSE.
      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY       = .FALSE.
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH           = .FALSE.
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER          = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING       = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING      = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING      = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES        = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT     = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY           = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY    = .FALSE.

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = .FALSE.

!  VLIDORT Fixed Control
!  ---------------------

!  Taylor ordering parameter. Should be set to 2 or 3
!     Added, 2/19/14 for Taylor-series expansions

      VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER     = 0

      VLIDORT_FixIn%Cont%TS_NSTOKES          = 0
      VLIDORT_FixIn%Cont%TS_NSTREAMS         = 0
      VLIDORT_FixIn%Cont%TS_NLAYERS          = 0
      VLIDORT_FixIn%Cont%TS_NFINELAYERS      = 0
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = 0

      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY = ZERO

!  Special Options (RT Solutions Use Only)

      VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS     = 0
      VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF   = 0

!  Water-leaving: Control for iterative calculation of transmittanaces
!    Variables added for Version 2.8, 7/6/16

      VLIDORT_FixIn%Cont%TS_TF_MAXITER       = 0
      VLIDORT_FixIn%Cont%TS_TF_CRITERION     = ZERO

!  TOA/BOA Illumination. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized
      
      VLIDORT_FixIn%Cont%TS_TOA_ILLUMINATION       = ZERO
      VLIDORT_FixIn%Cont%TS_BOA_ILLUMINATION       = ZERO

!  VLIDORT Modified Control
!  ------------------------

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 0

!  VLIDORT Fixed Sunrays
!  ---------------------

      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR   = ZERO

!  VLIDORT Modified Sunrays
!  ------------------------

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES   = 0
      VLIDORT_ModIn%MSunrays%TS_SZANGLES     = ZERO

!  VLIDORT Fixed UserValues
!  ------------------------

      VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = 0

!  VLIDORT Modified UserValues
!  ---------------------------

      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS       = 0
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS         = ZERO

      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES      = 0
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT  = ZERO

!  User-defined vertical level output. From Top-of-atmosphere. Example for 4 outputs:
!     USER_LEVELS(1) = 0.0           --> Top-of-atmosphere
!     USER_LEVELS(2) = 1.1           --> One tenth of the way down into Layer 2
!     USER_LEVELS(3) = 3.5           --> One half  of the way down into Layer 4
!     USER_LEVELS(4) = dble(NLAYERS) --> Bottom of atmosphere

      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS          = ZERO

      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT  = ZERO

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS      = 0
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT  = ZERO

!  VLIDORT Fixed Chapman
!  ---------------------

      VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID       = ZERO

      VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID     = ZERO
      VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID  = ZERO

      VLIDORT_FixIn%Chapman%TS_FINEGRID          = 0

      VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = ZERO

!  VLIDORT Modified Chapman
!  ------------------------

!  Output from Chapman function calculations - can also be input (not used)

      !VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS = ZERO

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS     = ZERO

!  VLIDORT Fixed Optical
!  ---------------------

      VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT     = ZERO

      VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT  = ZERO

!  F-matrix input as an alternative to the use of expansion coefficients in
!   the single-scatter (SSCORR and FO) codes.
!  Introduced for Version 2.8 by R. Spurr, 02/08/16
!mick mod 9/19/2017 - changed names from "FMATRIX_INPUT_UP/DN" to "FMATRIX_UP/DN"
!                     for consistency with the companion linearized input variables

      VLIDORT_FixIn%Optical%TS_FMATRIX_UP            = ZERO
      VLIDORT_FixIn%Optical%TS_FMATRIX_DN            = ZERO

      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO     = ZERO

      VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT      = ZERO
      VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT      = ZERO

!  Special LTE variables (RT Solutions Use Only).
!   This has been superseded in Version 2.7. No longer required

      !VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT = ZERO
      !VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT  = ZERO

!  Rob Fix 3/18/15. Add Wavelength (Microns) as a Diagnostic

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH      = ZERO

!  VLIDORT Modified Optical
!  ------------------------

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT    = ZERO

!  VLIDORT Fixed Write
!  -------------------

      VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE          = .FALSE.

      VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT          = .FALSE.
      VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME    = ' '

      VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO       = .FALSE.
      VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME = ' '

      VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER        = .FALSE.
      VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME  = ' '

      VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS        = .FALSE.
      VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME  = ' '

! Finish

      RETURN
      END SUBROUTINE VLIDORT_INIT_INPUTS

!

      SUBROUTINE VLIDORT_Sup_Init ( VLIDORT_Sup )

      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT supplement input structure

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT) :: VLIDORT_Sup

!  Initialize VLIDORT supplement inputs
!  ====================================

      CALL VLIDORT_BRDF_Sup_Init   ( VLIDORT_Sup )
      CALL VLIDORT_SLEAVE_Sup_Init ( VLIDORT_Sup )
      CALL VLIDORT_SS_Sup_Init     ( VLIDORT_Sup )

!  Finish

      END SUBROUTINE VLIDORT_Sup_Init

!

      SUBROUTINE VLIDORT_BRDF_Sup_Init ( VLIDORT_Sup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT supplement input structure

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT) :: VLIDORT_Sup

!  Initialize VLIDORT brdf supplement inputs
!  =========================================

      VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = ZERO
      VLIDORT_Sup%BRDF%TS_BRDF_F_0        = ZERO
      VLIDORT_Sup%BRDF%TS_BRDF_F          = ZERO
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = ZERO
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F     = ZERO

      VLIDORT_Sup%BRDF%TS_EMISSIVITY      = ZERO
      VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY = ZERO

!  Finish

      END SUBROUTINE VLIDORT_BRDF_Sup_Init

!

      SUBROUTINE VLIDORT_SLEAVE_Sup_Init ( VLIDORT_Sup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT supplement input structure

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT) :: VLIDORT_Sup

!  Initialize VLIDORT sleave supplement inputs
!  ===========================================

      VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC  = ZERO
      VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES = ZERO
      VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0        = ZERO
      VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0   = ZERO

!  Finish

      END SUBROUTINE VLIDORT_SLEAVE_Sup_Init

!

      SUBROUTINE VLIDORT_SS_Sup_Init ( VLIDORT_Sup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT supplement input structure

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT) :: VLIDORT_Sup

!  Initialize VLIDORT single-scatter supplement inputs
!  ===================================================

      VLIDORT_Sup%SS%TS_STOKES_SS = ZERO
      VLIDORT_Sup%SS%TS_STOKES_DB = ZERO

!  Finish

      END SUBROUTINE VLIDORT_SS_Sup_Init

!

      SUBROUTINE VLIDORT_READ_INPUTS ( &
      VLIDORT_FixIn, VLIDORT_ModIn,    & !InOut
      STATUS,                          & !Output
      NMESSAGES, MESSAGES, ACTIONS )     !InOut

!  Read standard inputs

!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added (1/31/16)
!     Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.
!     Version 2.8: FO/SS flags upgraded. 9/19/16, 3/1/17. Argument list rearranged.

!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination

!  1/31/21. Version 2.8.3. 
!    -- Control for reading DO_CLASSICAL_SOLUTION, DO_DOUBLET_GEOMETRY
!    -- Local arguments declared with extenxive commentary and reorganized

!  Read all control inputs for VLIDORT
!  -----------------------------------

      USE VLIDORT_PARS_m, Only : MAXMOMENTS_INPUT, MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAX_SZANGLES,   &
                                 MAX_MESSAGES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_USER_LEVELS, &
                                 MAX_USER_OBSGEOMS, MAX_THERMAL_COEFFS, MAXFINELAYERS,               &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_INUNIT, ONE

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_AUX_m , Only : GFINDPAR, FINDPAR_ERROR

      IMPLICIT NONE

!  InOut
!  -----

      TYPE(VLIDORT_Fixed_Inputs), INTENT(INOUT)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(INOUT) :: VLIDORT_ModIn

!  Outputs
!  -------

!  Exception handling
!     Message Length should be at least 120 Characters

      INTEGER, INTENT(OUT) ::                STATUS
      INTEGER, INTENT(INOUT) ::              NMESSAGES
      CHARACTER (LEN=*), INTENT(INOUT) ::    MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT(INOUT) ::    ACTIONS  ( 0:MAX_MESSAGES )

!  Local variables
!  ---------------

!  1. BOOLEANS
!  ===========

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL ::   DO_FULLRAD_MODE

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 2.8, 3/1/17

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      LOGICAL   :: DO_FOCORR

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first.

      LOGICAL   :: DO_FOCORR_EXTERNAL

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.
!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined in VLIDORT internally

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first.

      LOGICAL   :: DO_FOCORR_NADIR
      LOGICAL   :: DO_FOCORR_OUTGOING

!  1/31/21. Version 2.8.3.
!    -- Outgoing sphericity - criticality flag and attenuation criterion

      LOGICAL           :: DO_FOCORR_DOCRIT
      DOUBLE PRECISION  :: DO_FOCORR_ACRIT

!  Additional SSCORR flags
!  -----------------------

!  Flag for performing SSCORR truncation. Single-scattering
!     - Kept here, but disabled for Version 2.8.
!   1/31/21. Version 2.8.3. Flag removed.
!      LOGICAL   :: DO_SSCORR_TRUNCATION

!  Flag for using Phase Matrices in Single-scatter calculations (instead of GSF Expansion Coefficients)
!     - Introduced for Version 2.8, 3/3/17.  R. Spurr
!     - DO_FOCORR must be set first.
!     - Does not apply to thermal case.

      LOGICAL   :: DO_SSCORR_USEFMAT

!  Local SS variables (Version 2.7 and earlier)
!      LOGICAL   :: DO_SSCORR_NADIR
!      LOGICAL   :: DO_SSCORR_OUTGOING
!      LOGICAL   :: DO_SSCORR_TRUNCATION
!      LOGICAL   :: DO_SS_EXTERNAL
!      LOGICAL   :: DO_SSFULL 

!  Other flags
!  -----------

!  Basic top-level control

      LOGICAL :: DO_SOLAR_SOURCES
      LOGICAL :: DO_THERMAL_EMISSION

!  directional control

      LOGICAL :: DO_UPWELLING
      LOGICAL :: DO_DNWELLING

!  stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL :: DO_USER_VZANGLES

!  Observation-Geometry input control. 10/25/12

      LOGICAL :: DO_OBSERVATION_GEOMETRY

!  1/31/21. Version 2.8.3. Add DO_DOUBLET_GEOMETRY flag

      LOGICAL :: DO_DOUBLET_GEOMETRY

!  1/31/21. Version 2.8.3. New variable for RT Solution method
!            TRUE = Classical solution of PI, FALSE = Green's function solution

      LOGICAL :: DO_CLASSICAL_SOLUTION

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL :: DO_PLANE_PARALLEL

!  Transmittance only for thermal mode.

      LOGICAL :: DO_THERMAL_TRANSONLY

!  Flag for use of Lambertian surface
!    - If not set, default to BRDF surface

      LOGICAL :: DO_LAMBERTIAN_SURFACE

!  Surface emission flag

      LOGICAL :: DO_SURFACE_EMISSION

!  mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL :: DO_ADDITIONAL_MVOUT

!  mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL :: DO_MVOUT_ONLY

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set. 

      LOGICAL :: DO_CHAPMAN_FUNCTION

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set. 

      LOGICAL :: DO_REFRACTIVE_GEOMETRY

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL :: DO_DELTAM_SCALING

!  double convergence test flag

      LOGICAL :: DO_DOUBLE_CONVTEST

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL :: DO_SOLUTION_SAVING
      LOGICAL :: DO_BVP_TELESCOPING

!  scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!

      LOGICAL :: DO_RAYLEIGH_ONLY

!  Surface leaving control
!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.8: DO_WATER_LEAVING (added 10/28/15) and DO_FLUORESCENCE (added 1/31/16) flags.

      LOGICAL :: DO_SURFACE_LEAVING
      LOGICAL :: DO_SL_ISOTROPIC

!  Version 2.8: DO_WATER_LEAVING (added 10/28/15) and DO_FLUORESCENCE (added 1/31/16) flags.

      LOGICAL :: DO_WATER_LEAVING
      LOGICAL :: DO_FLUORESCENCE

!  Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.

      LOGICAL  ::  DO_TF_ITERATION
      INTEGER  ::  TF_MAXITER
      DOUBLE PRECISION  ::   TF_CRITERION
      
!  Water-leaving output flag. 3/18/19 for Version 2.8.1

      LOGICAL  ::  DO_WLADJUSTED_OUTPUT

!  4/29/19. Planetary problem and media properties flags,  Version 2.8.1

      LOGICAL  ::  DO_ALBTRN_MEDIA(2)
      LOGICAL  ::  DO_PLANETARY_PROBLEM

!  TOA/BOA illumination control, New Version 2.8.1, 3/23/19

      LOGICAL  ::  DO_TOA_ILLUMINATION
      LOGICAL  ::  DO_BOA_ILLUMINATION

      DOUBLE PRECISION :: TOA_ILLUMINATION
      DOUBLE PRECISION :: BOA_ILLUMINATION

!  Additional Control for Externalized water-leaving inputs. 3/18/19 for Version 2.8.1

      LOGICAL  ::  DO_EXTERNAL_WLEAVE

!  Following removed from Version 2.8, 7/6/16
!      LOGICAL ::             DO_QUAD_OUTPUT


! debug output control

      LOGICAL ::             DO_DEBUG_WRITE
      LOGICAL ::             DO_WRITE_INPUT
      LOGICAL ::             DO_WRITE_SCENARIO
      LOGICAL ::             DO_WRITE_FOURIER
      LOGICAL ::             DO_WRITE_RESULTS

      CHARACTER (LEN=60) ::  INPUT_WRITE_FILENAME
      CHARACTER (LEN=60) ::  SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60) ::  FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60) ::  RESULTS_WRITE_FILENAME

!  2. CONTROL INTEGERS
!  ===================

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7
      
      INTEGER :: TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  number of computational layers

      INTEGER :: NLAYERS

!  number of Stokes vector elements

      INTEGER :: NSTOKES

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER :: NFINELAYERS

!  number of GSF F-matrix expansion moments (External to the model)

      INTEGER :: NGREEK_MOMENTS_INPUT

!  number of solar beams to be processed

      INTEGER :: N_SZANGLES

!  Number of user-defined relative azimuths

      INTEGER :: N_USER_RELAZMS

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER :: N_USER_VZANGLES

!  Number of User-defined vertical levels for  output

      INTEGER :: N_USER_LEVELS

!  Number of thermal coefficients (2 should be the default)

      INTEGER :: N_THERMAL_COEFFS

!  Number of observational geometries. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      INTEGER :: N_USER_OBSGEOMS

!  Number of doublet geometries
!mick mod 1/5/2021 - declared N_USER_DOUBLETS

      INTEGER :: N_USER_DOUBLETS

!  3. CONTROL NUMBERS
!  ==================

!  Flux factor ( should be 1 or pi ). Same for all beams.

      DOUBLE PRECISION :: FLUX_FACTOR

!  accuracy for convergence of Fourier series

      DOUBLE PRECISION :: VLIDORT_ACCURACY

!  Zenith tolerance (nearness of output zenith cosine to 1.0 )
!    removed 02 June 2010
!      DOUBLE PRECISION :: ZENITH_TOLERANCE

!  Atmospheric wavelength, new Version 2.8, bookkeeping

      DOUBLE PRECISION :: ATMOS_WAVELENGTH

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      DOUBLE PRECISION :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      DOUBLE PRECISION :: RFINDEX_PARAMETER

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      DOUBLE PRECISION :: GEOMETRY_SPECHEIGHT

!  Lambertian Albedo

      DOUBLE PRECISION :: LAMBERTIAN_ALBEDO

!  BOA solar zenith angles (degrees)

      DOUBLE PRECISION :: SZANGLES ( MAX_SZANGLES )

!  user-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      DOUBLE PRECISION :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined viewing zenith angles input (degrees) 

      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      DOUBLE PRECISION :: USER_LEVELS  (MAX_USER_LEVELS)

!  User-defined Observation Geometry angle input
!   New variable, 25 OCtober 2012, for Observational Geometry input

      DOUBLE PRECISION :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  User-defined doublet geometry input
!mick mod 1/5/2021 - declared USER_DOUBLETS

      DOUBLE PRECISION :: USER_DOUBLETS (MAX_USER_VZANGLES,2)

!  Other variables
!  ===============

      CHARACTER (LEN=9), PARAMETER :: PREFIX = 'VLIDORT -'
      LOGICAL            :: ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER            :: I, FILUNIT, NM

!  Initialize status

      STATUS = VLIDORT_SUCCESS
      ERROR  = .FALSE.
      NM     = 0

!  These are already initialized in calling routine
!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
!      ACTIONS(0)      = 'No Action required for this Task'

!  File unit

      FILUNIT = VLIDORT_INUNIT

!  1. READ ALL CONTROL VARIABLES (BOOLEAN INPUTS)
!  ==============================================

!  Operation modes
!  ---------------

!  Full Stokes vector calculation, SS + MS fields
!    if False, then VLIDORT only does a multiple scatter calculation.

      PAR_STR = 'Do full Stokes vector calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FULLRAD_MODE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 2.8. 9/18/16, 3/1/17.

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      PAR_STR = 'Do First-Order (FO) correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All other options depend on this flag.

      IF ( DO_FOCORR ) THEN

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first. New 15 March 2012.

        PAR_STR = 'Do external First-Order correction?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_EXTERNAL
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.
!mick mod 9/19/2017 - Turned off.  Now defined in VLIDORT internally

        !PAR_STR = 'Do First-Order correction alone?'
        !IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_ALONE
        !CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Control for executing the FO code Version 1.5.
!    -  Only if the external field is turned off.

        IF ( .not. DO_FOCORR_EXTERNAL ) THEN

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first. DO_FOCORR_EXTERNAL should be off.

          PAR_STR = 'Do nadir First-Order correction?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_NADIR
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do outgoing First-Order correction?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_OUTGOING
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Flag for using the F-matrix in the SS calculations (instead of Greekmat Coefficients)
!     - Introduced for Version 2.8, 7/7/16.  R. Spurr

          PAR_STR = 'Do Fmatrix usage in single scatter correction?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SSCORR_USEFMAT
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Additional truncation scaling for single-scatter corrections
!   --- Disabled for Version 2.8, do we require it again ???
!        PAR_STR = 'Do truncation scaling on single scatter corrections?'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SSCORR_TRUNCATION
!        CALL FINDPAR_ERROR  ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End First-Order options

        ENDIF

      ENDIF

!  Direct beam correction (BRDF options only)
!      PAR_STR = 'Do direct beam correction?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DBCORRECTION
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Multiple scatter source function output control. Removed. 30 March 2007
!      PAR_STR = 'Output multiple scatter layer source functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) SAVE_LAYER_MSST
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  4/29/19 for Version 2.8.1. AlbTrn Media control
      
      PAR_STR = 'Do Media properties calculation with TOA illumination?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ALBTRN_MEDIA(1)
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      PAR_STR = 'Do Media properties calculation with BOA illumination?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ALBTRN_MEDIA(2)
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!   4/29/19 for Version 2.8.1. Planetary problem input

      PAR_STR = 'Do planetary problem calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_PLANETARY_PROBLEM
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Solar beam control
!  ------------------

!  Basic control

      PAR_STR = 'Use solar sources?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SOLAR_SOURCES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Other control options for the solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

!  Pseudo-spherical control

        PAR_STR = 'Do plane-parallel treatment of direct beam?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_PLANE_PARALLEL
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Internal Chapman function calculation

        PAR_STR = 'Do internal Chapman function calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_CHAPMAN_FUNCTION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Refractive atmosphere

        PAR_STR = 'Do refractive geometry?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_REFRACTIVE_GEOMETRY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  1/31/21/ Version 2.8.3. Read ew variable for RT Solution method
!            TRUE = Classical solution of PI, FALSE = Green's function solution

        PAR_STR = 'Do classical particular-integral solar-beam solution?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_CLASSICAL_SOLUTION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End control

      ENDIF

!  Direct-beam control
!      PAR_STR = 'Include direct beam?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DIRECT_BEAM
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Thermal Control
!  ---------------

!  Atmospheric thermal emission, Basic control

      PAR_STR = 'Do thermal emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_THERMAL_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_THERMAL_EMISSION ) THEN

!  Thermal sources, transmittance only

        PAR_STR = 'Do thermal emission, transmittance only?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_THERMAL_TRANSONLY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of coefficients (includes a dimensioning check)

        PAR_STR = 'Number of thermal coefficients'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_THERMAL_COEFFS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
             'Entry under "Number of thermal coefficients" >'// &
             ' allowed Maximum dimension'
            ACTIONS(NM)  = &
             'Re-set input value or increase MAX_THERMAL_COEFFS dimension '// &
             'in VLIDORT_PARS'
          STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
        ENDIF
      ENDIF

!  Surface emission control

      PAR_STR = 'Do surface emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Surface control
!  ---------------

!  Lambertian surface

      PAR_STR = 'Do Lambertian surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LAMBERTIAN_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!    START Surface Leaving control (New 17 May 2012)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

!  Basic control

      PAR_STR = 'Do surface-leaving term?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SURFACE_LEAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_SURFACE_LEAVING ) THEN
         
!  Additional Control for Externalized water-leaving inputs. Introduced 3/18/19 for Version 2.8.1
!    -- This is a general flag, you can still adjust the water-leaving (if flagged)         
         
         PAR_STR = 'Do external water-leaving production?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_EXTERNAL_WLEAVE
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Isotropic control

         PAR_STR = 'Do isotropic surface-leaving term?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SL_ISOTROPIC
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Water-leaving control. New 28 October 2015

         PAR_STR = 'Do Water-leaving option?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WATER_LEAVING
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  This section added 7/6/16 for Version 2.8.
!   The 3 TF inputs will control the water-leaving Transmittance calculation 
!      Also included now is the Water-leaving output flag. 3/18/19 for Version 2.8.1

         IF ( DO_WATER_LEAVING ) THEN
           PAR_STR = 'Do iterative calculation of Water-leaving transmittance?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_TF_ITERATION
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
 
           PAR_STR = 'Flag for output of transmittance-adjusted water-leaving radiances'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WLADJUSTED_OUTPUT
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

           IF ( DO_TF_ITERATION ) THEN

             PAR_STR = 'Maximum number of iterations in calculation of Water-leaving transmittance'
             IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TF_MAXITER
             CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

             PAR_STR = 'Convergence criterion for iterative calculation of Water-leaving transmittance'
             IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TF_CRITERION
             CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

           ENDIF
         ENDIF

!  Fluorescence control. New 31 January 2016
!   ( This is mutually exclusive from Water-leaving; condition will be checked )

         PAR_STR = 'Do Fluorescence option?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!     END Surface Leaving control
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!  TOA/BOA Illumination. 3/23/19 for Version 2.8.1  NEW SECTION
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

      PAR_STR = 'Do TOA Illumination for Airglow?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_TOA_ILLUMINATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_TOA_ILLUMINATION ) THEN
         PAR_STR = ' TOA Illumination Flux (sun-normalized)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TOA_ILLUMINATION
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

      PAR_STR = 'Do BOA Illumination for nighttime?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_BOA_ILLUMINATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_BOA_ILLUMINATION ) THEN
         PAR_STR = ' BOA Illumination Flux (sun-normalized)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) BOA_ILLUMINATION
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Performance control
!  -------------------

!  Delta-M scaling. Now set regardless of source
!    Surely a mistake (Versions 2.7 and earlier) --> Should only be set for solar beam sources

!      IF ( DO_SOLAR_SOURCES ) THEN
      PAR_STR = 'Do delta-M scaling?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DELTAM_SCALING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      ENDIF

!  Double convergence test

      PAR_STR = 'Do double convergence test?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DOUBLE_CONVTEST
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Solution saving mode.
!    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Do solution saving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SOLUTION_SAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Boundary value problem (BVP) telescope mode.
!    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Do boundary-value telescoping?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_BVP_TELESCOPING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  1/31/21. Version 2.8.3. Expand this section on post-processing options
!    -- Originally added 10/25/12 for observational Geometry, 4/25/20 for Doublet Geometry
!    -- Options turned off if no solar sources.
!mick mod 1/5/2021 - modifed this section to put asking for obsgeo and doublet inputs on an equal footing
!                    (similar to VBRDF and VSLEAVE supplements)

      IF ( DO_SOLAR_SOURCES ) THEN
         PAR_STR = 'Do Observation Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_OBSERVATION_GEOMETRY
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

         PAR_STR = 'Do Doublet Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DOUBLET_GEOMETRY
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

         IF ( DO_OBSERVATION_GEOMETRY .and. DO_DOUBLET_GEOMETRY) THEN
           NM = NM + 1
           MESSAGES(NM) = 'Not allowed to have both Observation and Doublet Geometry options'
           ACTIONS(NM)  = 'Re-set input to one or the other of these flags'
           STATUS       = VLIDORT_SERIOUS
           NMESSAGES    = NM
           RETURN
         ENDIF
      ELSE
         DO_OBSERVATION_GEOMETRY = .FALSE.
         DO_DOUBLET_GEOMETRY     = .FALSE.
      ENDIF

!  User-defined output control
!  ---------------------------

!  Directional output control

      PAR_STR = 'Do upwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_UPWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do downwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DNWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Stream angle and optical depth output control
!  ---------------------------------------------

!  User-defined viewing zenith angle

!  Removed.
!      PAR_STR = 'Include quadrature angles in output?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR) READ (FILUNIT,*,ERR=998) DO_QUAD_OUTPUT
!      CALL FINDPAR_ERROR( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_USER_VZANGLES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Old code, version 2.3 and earlier............
!      PAR_STR = 'User-defined optical depths?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_USER_TAUS
!      CALL FINDPAR_ERROR( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Layer boundary optical depths?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_LBOUND_TAUS
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Mean-value output control
!  -------------------------

      PAR_STR = 'Do mean-value output additionally?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ADDITIONAL_MVOUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do only mean-value output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_MVOUT_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Numerical control (azimuth series)
!  ----------------------------------

!  Scatterers and phase function control

      PAR_STR='Do Rayleigh atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_RAYLEIGH_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Removed isotropic-only option, 17 January 2006.
!      PAR_STR='Isotropic atmosphere only?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ISOTROPIC_ONLY
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  No azimuth dependence (TRUE means Fourier m = 0 only )
!    Removed. Now this is a bookkeeping variable.
!    17 January 2006.
!      PAR_STR = 'No azimuth dependence in the solution?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All possible Fourier components (2N-1). Debug only
!      PAR_STR = 'Compute all Fourier components?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ALL_FOURIER
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Write control
!  -------------

!  Output write flags

      PAR_STR = 'Do debug write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DEBUG_WRITE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do input control write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WRITE_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do input scenario write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WRITE_SCENARIO
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do Fourier component output write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WRITE_FOURIER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'DO results write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WRITE_RESULTS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Output filenames

      PAR_STR = 'filename for input write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,'(a)',ERR=998) INPUT_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'filename for scenario write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,'(a)',ERR=998) SCENARIO_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'filename for Fourier output write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,'(a)',ERR=998) FOURIER_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'filename for main output'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,'(a)',ERR=998) RESULTS_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2. READ ALL THE CONTROL NUMBERS (INTEGER INPUTS)
!  ================================================

!  Streams/layers/finelayers/moments (INTEGER input)
!  -------------------------------------------------

!  Taylor order parameter added 2/19/14, Version 2p7

      PAR_STR = 'Number of small-number terms in Taylor series expansions'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TAYLOR_ORDER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of Stokes parameters

      PAR_STR = 'Number of Stokes vector components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NSTOKES .GT. MAXSTOKES ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under' &
              //' "Number of Stokes parameters" >'//' allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value to 4 or less'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Number of computational streams

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
               'Entry under "Number of half-space streams" >'// ' allowed Maximum dimension'
        ACTIONS(NM)  = &
             'Re-set input value or increase MAXSTREAMS dimension '// 'in VLIDORT_PARS'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Number of atmospheric layers

      PAR_STR = 'Number of atmospheric layers'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NLAYERS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
               'Entry under "Number of atmospheric layers" >'// ' allowed Maximum dimension'
        ACTIONS(NM)  = &
             'Re-set input value or increase MAXLAYERS dimension '// 'in VLIDORT_PARS'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Change for Version 2.8. 3/1/17
!mick fix 9/19/2017 - now added DO_FOCORR & DO_FOCORR_EXTERNAL if conditions
!                     to DO_FOCORR_OUTGOING if condition

      !IF ( DO_SSCORR_OUTGOING .OR. DO_SSFULL ) THEN
      IF ( DO_FOCORR ) THEN
        IF ( .NOT.DO_FOCORR_EXTERNAL ) THEN
          IF ( DO_FOCORR_OUTGOING ) THEN

!  Number of fine layers

            PAR_STR = 'Number of fine layers (outgoing sphericity option only)'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NFINELAYERS
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

            IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Entry under "Number of fine layers..." >'//' allowed Maximum dimension'
              ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension in VLIDORT_PARS'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

!  Input Geometry specification height

            PAR_STR = 'Input geometry specification height (km)'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
               READ (FILUNIT,*,ERR=998) GEOMETRY_SPECHEIGHT
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          ENDIF
        ENDIF
      ENDIF

!  Number of scattering matrix expansion coefficients

      PAR_STR = 'Number of scattering matrix expansion coefficients'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NGREEK_MOMENTS_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of input Scattering Matrix expansion coefficients" >'//' allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXMOMENTS_INPUT dimension '//'in VLIDORT_PARS'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Geometry inputs
!  ---------------

!  1/31/21. Version 2.8.3. Add whole new section for the DOUBLET GEOMETRY settings.

!  Observational Geometry control.  New, 25 October 2012

      IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  Number of Observational Geometry inputs

        PAR_STR = 'Number of Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_USER_OBSGEOMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of Observation Geometry inputs" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          RETURN
        ENDIF

!  Observational Geometry inputs

        PAR_STR = 'Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_OBSGEOMS
             READ (FILUNIT,*,ERR=998) USER_OBSGEOMS(I,1), USER_OBSGEOMS(I,2), USER_OBSGEOMS(I,3)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS and DO_USER_VZANGLES

        N_SZANGLES       = N_USER_OBSGEOMS
        N_USER_VZANGLES  = N_USER_OBSGEOMS
        N_USER_RELAZMS   = N_USER_OBSGEOMS
        DO_USER_VZANGLES = .TRUE.

!   Automatic setting of SZANGLES, USER_VZANGLES, and USER_RELAZMS

        SZANGLES     (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
        USER_VZANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        USER_RELAZMS (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)

!  1/31/21. Version 2.8.3. Add whole new section for the DOUBLET GEOMETRY settings.

      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN

!  1. Number of Solar zenith angles
!     ---- check not exceeding dimensioned number

        PAR_STR = 'Number of solar zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_SZANGLES
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_SZANGLES .GT. MAX_SZANGLES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of solar zenith angles" >'//' allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_SZANGLES dimension '//'in VLIDORT_PARS'
          STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
        ENDIF

!  BOA solar zenith angle inputs

        PAR_STR = 'Solar zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_SZANGLES
            READ (FILUNIT,*,ERR=998) SZANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!mick mod 1/5/2021 - replaced reading of lattice VZAs and RAAs sections here with these
!                    doublet sections (similar to VBRDF and VSLEAVE supplements) 

!  Number of Doublet Geometries

        PAR_STR = 'Number of Doublet Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_DOUBLETS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

        IF ( N_USER_DOUBLETS .GT. MAX_USER_VZANGLES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of doublet geometry inputs > Maximum dimension MAX_USER_VZANGLES'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_VZANGLES dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          RETURN
        ENDIF

!  Read Doublet Geometry values

        PAR_STR = 'Doublet Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_DOUBLETS
             READ (FILUNIT,*,ERR=998)USER_DOUBLETS(I,1), USER_DOUBLETS(I,2)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of N_USER_STREAMS, N_USER_RELAZMS

         N_USER_VZANGLES  = N_USER_DOUBLETS
         N_USER_RELAZMS   = N_USER_DOUBLETS
         DO_USER_VZANGLES = .TRUE.

!  Automatic setting of user angles

         USER_VZANGLES(1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,1)
         USER_RELAZMS (1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,2)

!  Lattice Geometry control
!  1/31/21. Version 2.8.3. Add DO_DOUBLET_GEOMETRY flag to this option

      ELSE IF ( .not.DO_OBSERVATION_GEOMETRY .and. .not.DO_DOUBLET_GEOMETRY ) THEN

!  1. Number of Solar zenith angles
!     ---- check not exceeding dimensioned number

        PAR_STR = 'Number of solar zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_SZANGLES
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_SZANGLES .GT. MAX_SZANGLES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of solar zenith angles" >'//' allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_SZANGLES dimension '//'in VLIDORT_PARS'
          STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
        ENDIF

!  BOA solar zenith angle inputs

        PAR_STR = 'Solar zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_SZANGLES
            READ (FILUNIT,*,ERR=998) SZANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2. User defined viewing zenith angles (should be positive)
!     ---- check not exceeding dimensioned number

        IF ( DO_USER_VZANGLES ) THEN
          PAR_STR = 'Number of user-defined viewing zenith angles'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_VZANGLES
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          IF ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES ) THEN
            NM = NM + 1
            MESSAGES(NM) = &
               'Entry under "Number of .....viewing zenith angles" >'// &
               ' allowed Maximum dimension'
            ACTIONS(NM)  = &
               'Re-set input value or increase MAX_USER_VZANGLES dimension '// &
               'in VLIDORT_PARS'
            STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
          ENDIF

          PAR_STR = 'User-defined viewing zenith angles (degrees)'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            DO I = 1, N_USER_VZANGLES
              READ (FILUNIT,*,ERR=998) USER_VZANGLES(I)
            ENDDO
          ENDIF
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  3. Number of azimuths
!     ---- check not exceeding dimensioned number

        PAR_STR = 'Number of user-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of ..... azimuth angles" >'//' allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension '//'in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS ; NMESSAGES    = NM ; RETURN
        ENDIF

!  Azimuth angles

        PAR_STR = 'User-defined relative azimuth angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End lattice geometry clause

      ENDIF

!  Continuation point for skipping lattice geometry input-reads. 
!5665 CONTINUE   !!! Goto statement removed, version 2.8, replaced by If-block.

!  4. Number of output levels
!     ---- check not exceeding dimensioned number

!  This is designed to ensure that output is not related to optical
!  depth (which is wavelength-dependent); we use a height-based system.

!  User defined boundary (whole layer) and off-boundary (partial layer)
!  output choices are specified as follows.
!       USER_LEVELS(1) = 0.0    Top of the first layer
!       USER_LEVELS(2) = 1.0    Bottom of the first layer
!       USER_LEVELS(3) = 17.49  Output is in Layer 18, at a distance of
!                               0.49 of the way down from top (in height

      PAR_STR = 'Number of user-defined output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) N_USER_LEVELS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
             'Entry under "Number of ..... output levels" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_LEVELS dimension in VLIDORT_PARS'
        STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
      ENDIF

!  Vertical output levels

      PAR_STR = 'User-defined output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_LEVELS
          READ (FILUNIT,*,ERR=998) USER_LEVELS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  3. READ ALL THE FLOATING-POINT INPUTS
!  =====================================

!  Flux constant. Should be set to 1 if no solar sources.
!   Formerly a Book-keeping variable set to 1 in "derive inputs"
!   Now (July 2009) allowed to vary because of thermal emission
!   Must take physical values if using solar + thermal.

      IF ( DO_SOLAR_SOURCES ) THEN
        PAR_STR = 'Solar flux constant'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) FLUX_FACTOR
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ELSE
        FLUX_FACTOR = ONE
      ENDIF

!  Note the following possibility (not yet allowed for)
!      PAR_STR = 'TOA flux vector'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) (FLUXVEC(I),I=1,MAXSTOKES)
!      CALL FINDPAR_ERROR( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Accuracy criterion

      PAR_STR = 'Fourier series convergence'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) VLIDORT_ACCURACY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Zenith tolerance level. Now removed. 17 January 2006.

!  Atmospheric Wavelength [Microns].
!  Configuration file input added for Version 2.8, 1/31/16.
!  Just a diagnostic, only used for comparison with BRDF/SLEAVE wavelengths

      PAR_STR = 'Atmospheric Wavelength [Microns]'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) ATMOS_WAVELENGTH
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Earth radius and RF parameter
!  ---- only for Chapman function calculation
!mick mod 9/19/2017 - added DO_SOLAR_SOURCES if condition

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN
          PAR_STR = 'Earth radius (km)'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
               READ (FILUNIT,*,ERR=998) EARTH_RADIUS
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
            PAR_STR = 'Refractive index parameter'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
                 READ (FILUNIT,*,ERR=998) RFINDEX_PARAMETER
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          ENDIF
        ENDIF
      ENDIF

!  Lambertian surface

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        PAR_STR = 'Lambertian albedo'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) LAMBERTIAN_ALBEDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Define VLIDORT std inputs
!  =========================

!mick mod 9/19/2017 - initializing of all std VLIDORT type structure input variables is done
!                     in subroutine VLIDORT_INIT_INPUTS; only those read in from the VLIDORT 
!                     config file are possibly modified here.  IF conditions are applied where
!                     appropriate.
      
!  TOA/BOA Illumination. 3/23/19 for Version 2.8.1

!  VLIDORT Fixed Boolean Inputs

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE       = DO_FULLRAD_MODE
      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION   = DO_THERMAL_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION   = DO_SURFACE_EMISSION
      IF ( DO_SOLAR_SOURCES ) &
         VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL  = DO_PLANE_PARALLEL 
      VLIDORT_FixIn%Bool%TS_DO_UPWELLING          = DO_UPWELLING
      VLIDORT_FixIn%Bool%TS_DO_DNWELLING          = DO_DNWELLING
      !VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT        = DO_QUAD_OUTPUT
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING    = DO_SURFACE_LEAVING
      
      IF ( DO_SURFACE_LEAVING ) THEN
         VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC    = DO_SL_ISOTROPIC
         VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   = DO_WATER_LEAVING
         IF ( DO_WATER_LEAVING ) THEN
            VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT = DO_WLADJUSTED_OUTPUT   ! New 3/18/19 Version 2.8.1
            VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION = DO_TF_ITERATION
            IF ( DO_TF_ITERATION ) THEN
               VLIDORT_FixIn%Cont%TS_TF_MAXITER   = TF_MAXITER
               VLIDORT_FixIn%Cont%TS_TF_CRITERION = TF_CRITERION
            ENDIF
         ENDIF
         VLIDORT_FixIn%Bool%TS_DO_FLUORESCENCE    = DO_FLUORESCENCE
      ENDIF

!  TOA/BOA Illumination. New code for Version 2.8.1, 3/18/19 
      
      VLIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION   = DO_TOA_ILLUMINATION
      if ( DO_TOA_ILLUMINATION ) THEN
         VLIDORT_FixIn%Cont%TS_TOA_ILLUMINATION   = TOA_ILLUMINATION
      endif
      VLIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION   = DO_BOA_ILLUMINATION
      if ( DO_BOA_ILLUMINATION ) THEN
         VLIDORT_FixIn%Cont%TS_BOA_ILLUMINATION   = BOA_ILLUMINATION
      endif

!  Flag for the Planetary problem calculation. Version 2.8.1, 4/29/19 

      VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM   = DO_PLANETARY_PROBLEM

!  Flags for the media properties calculation. Version 2.8.1, 4/29/19 

      VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA   = DO_ALBTRN_MEDIA

!  VLIDORT Modified Boolean Inputs
!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined internally
!mick fix 6/4/2019 - added "DO_SURFACE_LEAVING" IF block for DO_EXTERNAL_WLEAVE

      VLIDORT_ModIn%MBool%TS_DO_FOCORR                   = DO_FOCORR
      IF ( DO_FOCORR ) THEN
         VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL       = DO_FOCORR_EXTERNAL
         !VLIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE          = DO_FOCORR_ALONE
         IF ( .NOT.DO_FOCORR_EXTERNAL ) THEN
            VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR       = DO_FOCORR_NADIR 
            VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING    = DO_FOCORR_OUTGOING
            IF ( DO_FOCORR_OUTGOING ) THEN
              VLIDORT_FixIn%Cont%TS_NFINELAYERS             = NFINELAYERS
              VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT
            ENDIF
            VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT     = DO_SSCORR_USEFMAT
            !VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION  = DO_SSCORR_TRUNCATION
         ENDIF
      ENDIF
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST          = DO_DOUBLE_CONVTEST
      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES            = DO_SOLAR_SOURCES

!  1/31/21. Version 2.8.3.
!    -- Copy new variable DO_CLASSICAL_SOLUTION for RT Solution method
!    -- Copy new variable DO_DOUBLET_GEOMETRY
!    -- Only if solar sources is turned on
!mick mod 1/5/2021 - added IF section to copy doublet inputs

      IF ( DO_SOLAR_SOURCES ) THEN
         VLIDORT_ModIn%MBool%TS_DO_CLASSICAL_SOLUTION    = DO_CLASSICAL_SOLUTION

         VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION      = DO_CHAPMAN_FUNCTION
         IF ( DO_CHAPMAN_FUNCTION ) &
            VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS       = EARTH_RADIUS

         VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY   = DO_REFRACTIVE_GEOMETRY
         IF ( DO_CHAPMAN_FUNCTION .AND. DO_REFRACTIVE_GEOMETRY ) &
            VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = RFINDEX_PARAMETER

         VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY  = DO_OBSERVATION_GEOMETRY
         VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY      = DO_DOUBLET_GEOMETRY
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS    = N_USER_OBSGEOMS
            VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3) = &
                                      USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3)
         ENDIF
         IF ( DO_DOUBLET_GEOMETRY ) THEN
            VLIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS    = N_USER_DOUBLETS
            VLIDORT_ModIn%MUserVal%TS_USER_DOUBLETS(1:N_USER_DOUBLETS,1:2) = &
                                      USER_DOUBLETS(1:N_USER_DOUBLETS,1:2)
         ENDIF
      ENDIF

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY            = DO_RAYLEIGH_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY           = DO_ISOTROPIC_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH               = DO_NO_AZIMUTH
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER              = DO_ALL_FOURIER
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING           = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING          = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING          = DO_BVP_TELESCOPING
      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES            = DO_USER_VZANGLES
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT         = DO_ADDITIONAL_MVOUT
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY               = DO_MVOUT_ONLY
      IF ( DO_THERMAL_EMISSION ) &
         VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY     = DO_THERMAL_TRANSONLY

      ! New 4/22/19 Version 2.8.1.
      IF ( DO_SURFACE_LEAVING ) THEN
         VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE       = DO_EXTERNAL_WLEAVE
      ELSE
         VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE       = .false.  
      ENDIF

!  VLIDORT Fixed Control Inputs

      VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER    = TAYLOR_ORDER
      VLIDORT_FixIn%Cont%TS_NSTOKES         = NSTOKES
      VLIDORT_FixIn%Cont%TS_NSTREAMS        = NSTREAMS
      VLIDORT_FixIn%Cont%TS_NLAYERS         = NLAYERS
      IF ( DO_THERMAL_EMISSION ) &
         VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = N_THERMAL_COEFFS
      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY = VLIDORT_ACCURACY

!  VLIDORT Modified Control Inputs

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = NGREEK_MOMENTS_INPUT

!  VLIDORT Fixed Sunrays Inputs:

      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR  = FLUX_FACTOR

!  VLIDORT Modified Sunrays Inputs

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES = N_SZANGLES
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:N_SZANGLES) = SZANGLES(1:N_SZANGLES)

!  VLIDORT Fixed UserValues Inputs

      VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = N_USER_LEVELS

!  VLIDORT Modified UserValues Inputs

      !IF ( .NOT. DO_NO_AZIMUTH ) THEN
         VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS = N_USER_RELAZMS
         VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS) = &
                                   USER_RELAZMS(1:N_USER_RELAZMS) 
      !ENDIF
      IF ( DO_USER_VZANGLES ) THEN
         VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = N_USER_VZANGLES
         VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES) = &
                                   USER_VZANGLES(1:N_USER_VZANGLES)
      ENDIF
      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS) = &
                                USER_LEVELS(1:N_USER_LEVELS)

!  VLIDORT Fixed Optical Inputs

      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = LAMBERTIAN_ALBEDO
      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH  = ATMOS_WAVELENGTH

!  VLIDORT Fixed Write Inputs

      VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE          = DO_DEBUG_WRITE
      VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT          = DO_WRITE_INPUT
      VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME    = INPUT_WRITE_FILENAME
      VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO       = DO_WRITE_SCENARIO
      VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME = SCENARIO_WRITE_FILENAME
      VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER        = DO_WRITE_FOURIER
      VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME  = FOURIER_WRITE_FILENAME
      VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS        = DO_WRITE_RESULTS
      VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME  = RESULTS_WRITE_FILENAME

!  Normal return

!mick fix
      NMESSAGES = NM

      RETURN

!  Line read error - abort immediately

998   CONTINUE
      NM = NM + 1
      STATUS       = VLIDORT_SERIOUS
      MESSAGES(NM) = 'Read failure for entry below String: ' //Trim(Adjustl(PAR_STR))
      ACTIONS(NM)  = 'Re-set value: Entry wrongly formatted in Input file'
      NMESSAGES    = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_READ_INPUTS

!

      SUBROUTINE VLIDORT_CHECK_INPUT_DIMS &
      ( VLIDORT_FixIn, VLIDORT_ModIn, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  module, dimensions and numbers
!   Rob Fix 3/18/15. Added MAX_GEOMETRIES to this list..

      USE VLIDORT_pars_m, only : MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                                 MAX_SZANGLES, MAXSTREAMS, MAXLAYERS, MAX_GEOMETRIES,  &
                                 MAXMOMENTS_INPUT, MAX_THERMAL_COEFFS, MAXSTOKES,      &
                                 MAX_USER_OBSGEOMS, MAXFINELAYERS, MAX_MESSAGES,       &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE VLIDORT_Inputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)   , INTENT (IN) :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (IN) :: VLIDORT_ModIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS
      NM = NMESSAGES

!  Check VLIDORT input dimensions against maximum dimensions
!  =========================================================

!  1a. Basic dimensions - always checked

      IF ( VLIDORT_FixIn%Cont%TS_NSTOKES .GT. MAXSTOKES ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number Stokes parameters NSTOKES > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value to 4 or less'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_FixIn%Cont%TS_NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams NSTREAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_FixIn%Cont%TS_NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of layers NLAYERS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of user vertical output levels N_USER_LEVELS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_LEVELS dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Expansion coefficients NGREEK_MOMENTS_INPUT > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXMOMENTS_INPUT dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

!  1b. Basic dimensions - conditionally checked

      IF ( VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING  ) THEN
        IF ( VLIDORT_FixIn%Cont%TS_NFINELAYERS .GT. MAXFINELAYERS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of fine layers NFINELAYERS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
        IF ( VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of thermal coefficients" > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_THERMCOEFFS dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  2a. Geometry dimensions - always checked

      IF ( VLIDORT_ModIn%MSunrays%TS_N_SZANGLES .GT. MAX_SZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles N_SZANGLES > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_SZANGLES dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of relative azimuths N_USER_RELAZMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

!  2b. Geometry dimensions - conditionally checked

      IF ( VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES ) THEN
        IF ( VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES .GT. MAX_USER_VZANGLES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of user streams N_USER_VZANGLES > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_VZANGLES dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN
        IF ( VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'ObsGeo Case: Number of observation geometries N_USER_OBSGEOMS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Rob Fix 3/18/15. Geometry dimension to be checked for the Lattice case
!  1/31/21. Version 2.8.3.  Add DO_DOUBLET_GEOMETRY flag to this clause.

      IF ( .NOT.VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY .AND. &
           .NOT.VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY ) THEN
         IF ( MAX_GEOMETRIES .NE. MAX_SZANGLES*MAX_USER_VZANGLES*MAX_USER_RELAZMS ) then
          NM = NM + 1
          MESSAGES(NM) = 'Lattice Case: MAX_GEOMETRIES not equal to MAX_SZANGLES*MAX_USER_VZANGLES*MAX_USER_RELAZMS'
          ACTIONS(NM)  = 'Re-set MAX_GEOMETRIES dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  1/31/21. Version 2.8.3. Geometry dimension to be checked for the Doublet case
!mick mod 1/5/2021 - modified doublet case IF condition

      IF ( VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY ) THEN
        !IF ( MAX_GEOMETRIES.ne.MAX_SZANGLES*MAX_USER_VZANGLES ) then
        IF ( VLIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS .GT. MAX_USER_VZANGLES ) then
          NM = NM + 1
          MESSAGES(NM) = 'Doublet Case: Number of doublet geometries N_USER_DOUBLETS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_VZANGLES dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHECK_INPUT_DIMS

!
      
      SUBROUTINE VLIDORT_CHECK_INPUT ( &
           DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                     & ! Input
           DO_THERMAL_TRANSONLY, DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_FLUORESCENCE,      & ! Input
           DO_WATER_LEAVING, DO_TF_ITERATION, DO_WLADJUSTED_OUTPUT, DO_EXTERNAL_WLEAVE,           & ! Boolean 
           DO_FULLRAD_MODE, DO_RAYLEIGH_ONLY, DO_SSCORR_USEFMAT, DO_DIRECT_BEAM,                  & ! Boolean
           DO_USER_VZANGLES, DO_DELTAM_SCALING, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY,     & ! Boolean
           DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_ALONE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,   & ! Boolean
           DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION, & ! Boolean
           DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,            & ! Boolean
           DO_TOA_CONTRIBS, DO_ALL_FOURIER, DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, DO_MSSTS, & ! Boolean
           TAYLOR_ORDER, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, TF_MAXITER,                   & ! Integer
           N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS, N_USER_OBSGEOMS,      & ! Integer
           N_OUT_STREAMS, NLAYERS_NOMS, NLAYERS_CUTOFF, NGREEK_MOMENTS_INPUT, & ! Integer
           SZANGLES, USER_VZANGLES, USER_RELAZMS, USER_OBSGEOMS, USER_LEVELS, & ! Floating point
           OUT_ANGLES, TF_CRITERION, EARTH_RADIUS, GEOMETRY_SPECHEIGHT,       & ! Floating point
           HEIGHT_GRID, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,              & ! Floating point
           STATUS, NMESSAGES, MESSAGES, ACTIONS )                               ! Exception handling

!  1/31/21. Version 2.8.3. 
!   -- Add DO_CLASSICAL_SOLUTION to list arguments (Line 7). Also add NSTOKES (line 10)
!   -- Add DO_MSSTS flag to this list. Line 9, final Boolean. Arguments rearranged
!   -- Add DOUBLET_GEOMETRY to this list
!   -- Green's function only for solar sources with NSTOKES = 1 or 3. Add NSTOKES input.
!   -- Several restrictions on use of MSSTS option, need to be checked

      USE VLIDORT_PARS_m, Only : MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ, MAX_SZANGLES, MAX_USER_RELAZMS, &
                                 MAX_USER_VZANGLES, MAX_USER_LEVELS, MAX_USER_OBSGEOMS, MAX_USER_STREAMS,   &
                                 MAX_ALLSTRMS_P1, MAX_TAYLOR_TERMS, MAX_MESSAGES,                           &
                                 VLIDORT_WARNING, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, ONE, OMEGA_SMALLNUM

!  Revised Call 3/17/17

      USE VLIDORT_AUX_m , Only : RQSORT_IDX

      IMPLICIT NONE

!  %% DO_FULLRAD_MODE argument added (First line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged
      
!  Order of Taylor series (including terms up to EPS^n). Introduced 2/19/14 for Version 2.7
!  %% Revision, 3/3/17.   Complete overhaul of FO and SSCORR flags (Version 2.8)

!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 3/18/19 for Version 2.8.1. 
!  Add Flag DO_EXTERNAL_WLEAVE   (Water-leaving output). 4/22/19 for Version 2.8.1. 

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  %% 05 March 2013. Added FullRad mode flag

      LOGICAL  , intent(in)  :: DO_FULLRAD_MODE

!  Plane-parallel solution

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Thermal inputs

      LOGICAL  , intent(in)  :: DO_THERMAL_EMISSION

!  surface

      LOGICAL  , intent(in)  :: DO_LAMBERTIAN_SURFACE
      LOGICAL  , INTENT (IN) :: DO_SURFACE_LEAVING

!  Added 1/31/16 and 7/7/16 for Version 2.8
!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 3/18/19 for Version 2.8.1. 

      LOGICAL  , INTENT (IN) :: DO_FLUORESCENCE
      LOGICAL  , INTENT (IN) :: DO_WATER_LEAVING
      LOGICAL  , INTENT (IN) :: DO_TF_ITERATION
      LOGICAL  , INTENT (IN) :: DO_WLADJUSTED_OUTPUT
      LOGICAL  , INTENT (IN) :: DO_EXTERNAL_WLEAVE

!  1/31/21. Version 2.8.3.  MSST flag declaration

      LOGICAL, INTENT (IN) ::             DO_MSSTS

!  Specialist flags

      LOGICAL, INTENT (IN) ::             DO_SPECIALIST_OPTION_2
      LOGICAL, INTENT (IN) ::             DO_SPECIALIST_OPTION_3
      LOGICAL, INTENT (IN) ::             DO_TOA_CONTRIBS

!  Modified input flags
!  --------------------

!  Basic top-level control

      LOGICAL, intent(inout)  ::          DO_SOLAR_SOURCES    ! May be re-set with a Warning

!  FO and SSCORR Booleans

      LOGICAL, INTENT (INOUT) ::          DO_FOCORR
      LOGICAL, INTENT (INOUT) ::          DO_FOCORR_EXTERNAL
      LOGICAL, INTENT (INOUT) ::          DO_FOCORR_ALONE
      LOGICAL, INTENT (INOUT) ::          DO_FOCORR_NADIR
      LOGICAL, INTENT (INOUT) ::          DO_FOCORR_OUTGOING
      LOGICAL, INTENT (INOUT) ::          DO_SSCORR_USEFMAT

!  1/31/21. Version 2.8.3.
!    -- Add CLASSICAL_SOLUTION flag, declared here

      LOGICAL, INTENT (INOUT) ::          DO_CLASSICAL_SOLUTION

!  Beam particular solution pseudo-spherical options

      LOGICAL, INTENT (INOUT) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (INOUT) ::          DO_CHAPMAN_FUNCTION

!  RT Mode flags
!    1/31/21. Version 2.8.3.  DO_DOUBLET_GEOMETRY flag declaration

      LOGICAL, INTENT (IN) ::             DO_DOUBLET_GEOMETRY
      LOGICAL, INTENT (IN) ::             DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (INOUT) ::          DO_MVOUT_ONLY
      LOGICAL, INTENT (INOUT) ::          DO_USER_VZANGLES
      LOGICAL, INTENT (INOUT) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (INOUT) ::          DO_THERMAL_TRANSONLY

!  Performance flags

      LOGICAL, INTENT (INOUT) ::          DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (INOUT) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (INOUT) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT (INOUT) ::          DO_BVP_TELESCOPING

!  INTEGERS. (TF_MAXITER added 7/7/16, Version 2.8)
!    1/31/21. Version 2.8.3.  NSTOKES declaration added

      INTEGER, INTENT (IN) ::             TAYLOR_ORDER
      INTEGER, INTENT (IN) ::             NSTOKES
      INTEGER, INTENT (IN) ::             NSTREAMS
      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             TF_MAXITER

      INTEGER, INTENT (IN) ::             NLAYERS_NOMS
      INTEGER, INTENT (IN) ::             NLAYERS_CUTOFF
      INTEGER, INTENT (INOUT) ::          NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (INOUT) ::          N_SZANGLES
      INTEGER, INTENT (INOUT) ::          N_USER_RELAZMS
      INTEGER, INTENT (INOUT) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::             N_USER_OBSGEOMS
      INTEGER, INTENT (IN) ::             N_USER_LEVELS

!  FLOATING POINT. (TF_CRITERION added 7/7/16, Version 2.8)

      DOUBLE PRECISION, INTENT (IN) ::    OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT (INOUT) :: SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN)    :: TF_CRITERION
      DOUBLE PRECISION, INTENT (INOUT) :: EARTH_RADIUS
      DOUBLE PRECISION, INTENT (IN) ::    HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: GEOMETRY_SPECHEIGHT

      DOUBLE PRECISION, INTENT (IN) ::    USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )
      DOUBLE PRECISION, INTENT (INOUT) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT (INOUT) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (INOUT) :: USER_LEVELS ( MAX_USER_LEVELS )

!  PURE OUTPUT VARIABLES

      LOGICAL, INTENT (OUT) ::            DO_ALL_FOURIER
      LOGICAL, INTENT (OUT) ::            DO_DIRECT_BEAM

      INTEGER, INTENT (OUT) ::            N_OUT_STREAMS
      DOUBLE PRECISION, INTENT (OUT) ::   OUT_ANGLES ( MAX_USER_STREAMS )

!  Following flags removed from Version 2.8, 7/6/16
!      LOGICAL, INTENT (IN) ::             DO_QUAD_OUTPUT
!      LOGICAL, INTENT (OUT) ::            DO_CLASSICAL_SOLUTION

!  Exception handling.

      INTEGER, INTENT (OUT) ::              STATUS
      INTEGER, INTENT (INOUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGES( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT (INOUT) ::  ACTIONS ( 0:MAX_MESSAGES )

!  Local variables

      INTEGER ::           I, L, N, UTA, NSTART, NALLSTREAMS, NM
      INTEGER ::           INDEX_ANGLES ( MAX_USER_STREAMS )
      DOUBLE PRECISION  :: XT, ALL_ANGLES ( MAX_USER_STREAMS )
      CHARACTER (LEN=3) :: C3
      CHARACTER (LEN=2) :: C2
      LOGICAL ::           LOOP

!  Initialize output status

      STATUS = VLIDORT_SUCCESS

!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of VLIDORT Basic Input'
!      ACTIONS(0)      = 'No Action required for this Task'

      NM = NMESSAGES

!  Automatic input
!    Flux factor set to unity. (21 December 2005)

      DO_ALL_FOURIER        = .FALSE.
!      DO_CLASSICAL_SOLUTION = .TRUE.  ! Removed for Version 2.8
      DO_DIRECT_BEAM        = .TRUE.

!  Check top level Solar/Thermal options, set warnings
!  ---------------------------------------------------

!  Check thermal or Solar sources present

      IF ( .NOT.DO_SOLAR_SOURCES.AND..NOT.DO_THERMAL_EMISSION ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: DO_SOLAR_SOURCES,DO_THERMAL_EMISSION both not set'
        ACTIONS(NM)  = 'Abort: must set DO_SOLAR_SOURCES and/or DO_THERMAL_EMISSION'
        STATUS = VLIDORT_SERIOUS
      ENDIF

!  Switch off several flags with thermal-only option
!    Set default, regardless of whether solar sources are on. Revised, 2.8
!mick fix 9/19/2017 - turned off

      !IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
      ! IF ( DO_FOCORR .and. DO_FOCORR_NADIR ) THEN
      !   NM = NM + 1
      !   MESSAGES(NM) = 'Switch off FO Nadir correction, not needed for thermal-only'
      !   ACTIONS(NM)  = 'Warning: FOCORR_NADIR correction flag turned off internally'
      !   STATUS = VLIDORT_WARNING
      !   DO_FOCORR_NADIR = .FALSE.
      !  ENDIF
      !ENDIF

!  Set number of solar angles N_SZANGLES to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
        IF ( N_SZANGLES .NE. 1 ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: N_SZANGLES not set to 1 for thermal-only'
         ACTIONS(NM) = 'Warning: N_SZANGLES is set to 1 internally'
         STATUS = VLIDORT_WARNING
         N_SZANGLES = 1
        ENDIF
      ENDIF

!  Set number of azimuths N_USER_RELAZMS to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
        IF ( N_USER_RELAZMS .NE. 1 ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: N_USER_RELAZMS not set to 1 for thermal-only'
         ACTIONS(NM)  = 'Warning: N_USER_RELAZMS is set to 1 internally'
         STATUS = VLIDORT_WARNING
         N_USER_RELAZMS = 1
        ENDIF
      ENDIF

!  1/31/21. Version 2.8.3.  New check for validity of the Green's function choice
!    -- Cannot be used with thermal emission ==> Warning flag
!    -- Cannot be used if NSTOKES = 4        ==> Fatal flag

      IF ( .not. DO_CLASSICAL_SOLUTION ) then
        IF ( DO_THERMAL_EMISSION ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Greens function solution cannot be used for thermal emission'
          ACTIONS(NM)  = 'Warning  : Set DO_CLASSICAL_SOLUTION to TRUE internally'
          STATUS = VLIDORT_WARNING
          DO_CLASSICAL_SOLUTION = .true.
        ENDIF
      ENDIF

      IF ( .not. DO_CLASSICAL_SOLUTION ) then
        IF ( DO_SOLAR_SOURCES .and. NSTOKES.eq.4 ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Greens function solution not valid (YET) for NSTOKES = 4'
          ACTIONS(NM)  = 'Abort    : Set DO_CLASSICAL_SOLUTION to TRUE'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  New code for the Observational Geometry
!  ---------------------------------------

      IF ( DO_OBSERVATION_GEOMETRY ) THEN

        IF ( .NOT.DO_SOLAR_SOURCES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: DO_SOLAR_SOURCES not set for Observation Geometry option'
          ACTIONS(NM)  = 'Abort: must set DO_SOLAR_SOURCE'
          STATUS = VLIDORT_SERIOUS
        ENDIF

!%%% No reason you cannot have Observation Geometry with Cross-Over
!%%% Clause commented out, 05 March 2013
!        IF ( DO_THERMAL_EMISSION ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = &
!            'Bad input: DO_THERMAL_EMISSION should not be set for Observation Geometry option'
!          ACTIONS(NM)  = &
!            'Abort: must turn off DO_THERMAL EMISSION'
!          STATUS = VLIDORT_SERIOUS
!        ENDIF
!%%% End Clause commented out, 05 March 2013

!  Observational Geometry control
!   New, 25 October 2012
!     ---- Automatic setting of N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS,
!          and DO_USER_VZANGLES

        N_SZANGLES       = N_USER_OBSGEOMS
        N_USER_VZANGLES  = N_USER_OBSGEOMS
        N_USER_RELAZMS   = N_USER_OBSGEOMS
        DO_USER_VZANGLES = .true.

!     ---- Automatic setting of SZANGLES, USER_VZANGLES, and USER_RELAZMS

        SZANGLES     (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
        USER_VZANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        USER_RELAZMS (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)

      ENDIF

!  Check Taylor series parameter is not out of range.
!   This code was added 2/19/14 for Version 2.7

      IF ( (TAYLOR_ORDER .GT. MAX_TAYLOR_TERMS-2) .OR. (TAYLOR_ORDER .LT. 0) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Taylor series Order parameter out of range'
        ACTIONS(NM)  = 'Re-set input value, should be in range 0-4'
        STATUS       = VLIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  Check inputs (both file-read and derived)
!  -----------------------------------------

!  Check Chapman function options

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( .NOT. DO_CHAPMAN_FUNCTION ) THEN
        IF ( DO_PLANE_PARALLEL ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Chapman Function not set, plane parallel'
          ACTIONS(NM) = 'Warning: Chapman function set internally'
          STATUS = VLIDORT_WARNING
          DO_CHAPMAN_FUNCTION = .TRUE.
        ELSE
          NM = NM + 1
          MESSAGES(NM)   = 'Chapman Function not set, pseudo-spherical'
          ACTIONS(NM) = 'Have you set the CHAPMAN_FACTORS values?'
          STATUS = VLIDORT_SERIOUS
        ENDIF
       ENDIF
      ENDIF

!  Check sphericity corrections for Solar....Cannot both be turned on
!    --------- New code 31 January 2007
!mick fix 9/19/2017 - turned off "DO_SOLAR_SOURCES" IF condition
!                   - included additional DO_FOCORR_NADIR/DO_FOCORR_OUTGOING check 

! 6/17/20. Bug. Must allow for plane parallel FO capability

      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_FOCORR ) THEN

          IF ( DO_FOCORR_NADIR .and. DO_FOCORR_OUTGOING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Cannot have both FO corrections on for Solar single-scatter'
            ACTIONS(NM)  = 'Turn off DO_FOCORR_NADIR and/or DO_FOCORR_OUTGOING'
            STATUS = VLIDORT_SERIOUS
          ENDIF

! 6/17/20 Comment out this clause
!          IF (      ( DO_FOCORR_NADIR .AND. DO_FOCORR_OUTGOING ) .OR. &
!               .NOT.( DO_FOCORR_NADIR .OR.  DO_FOCORR_OUTGOING )         )  THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Cannot have both FO correction types on or off when FO correction is in use'
!            ACTIONS(NM)  = 'Set one of DO_FOCORR_NADIR or DO_FOCORR_OUTGOING'
!            STATUS = VLIDORT_SERIOUS
!          ENDIF

        ENDIF
      !ENDIF

!  New for Version 2.8, Check surface leaving inputs
!  -------------------------------------------------

      IF ( DO_SURFACE_LEAVING) THEN

!  Check compatibility. 1/31/16

        IF ( DO_WATER_LEAVING .and. DO_FLUORESCENCE ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Cannot have both Water-leaving and Fluorescence turned on'
          ACTIONS(NM)    = 'Turn off DO_WATER_LEAVING or DO_FLUORESCENCE'
          STATUS = VLIDORT_SERIOUS
        ENDIF

        IF ( .not.DO_WATER_LEAVING .and. .not.DO_FLUORESCENCE ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Surface-leaving is ON, but both Water-leaving and Fluorescence are OFF!'
          ACTIONS(NM)    = 'Turn on either DO_WATER_LEAVING or DO_FLUORESCENCE'
          STATUS = VLIDORT_SERIOUS
        ENDIF

!  Check Use of external water-leaving source. 4/22/19 for Version 2.8.1
        
        IF ( DO_EXTERNAL_WLEAVE .and..not.DO_WATER_LEAVING  ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'External Water-leaving source flag is ON, but main Water-leaving Flag is OFF!'
          ACTIONS(NM)    = 'Turn on DO_WATER_LEAVING if you want to use DO_EXTERNAL_WLEAVE'
          STATUS = VLIDORT_SERIOUS
        ENDIF
        IF ( DO_EXTERNAL_WLEAVE .and. DO_WATER_LEAVING .and. DO_TF_ITERATION  ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Both External Water-leaving source flag and TF_ITERATION flag are on'
          ACTIONS(NM)    = 'Turn off DO_TF_ITERATION if you want to use DO_EXTERNAL_WLEAVE'
          STATUS = VLIDORT_SERIOUS
        ENDIF

!  Check Water-leaving output flag. 3/18/19 for Version 2.8.1
       
        IF ( DO_WLADJUSTED_OUTPUT .and..not.DO_WATER_LEAVING  ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Water-leaving output is ON, but main Water-leaving Flag is OFF!'
          ACTIONS(NM)    = 'Turn on DO_WATER_LEAVING if you want to use DO_WLADJUSTED_OUTPUT'
          STATUS = VLIDORT_SERIOUS
        ENDIF

!  Check TF conditions. 7/8/16

        IF ( DO_WATER_LEAVING .and. DO_TF_ITERATION ) THEN
          IF ( TF_MAXITER .gt. 10 ) then
            NM = NM + 1
            MESSAGES(NM)   = 'Water-leaving Transmittance calculation, cannot have more than 10 iterations'
            ACTIONS(NM)    = 'Set input TF_MAXITER to 10 or less'
            STATUS = VLIDORT_SERIOUS
          ENDIF
          IF ( TF_CRITERION .gt. 0.01d0 .or. TF_CRITERION.lt. 1.0d-06 ) then
            NM = NM + 1
            MESSAGES(NM)   = 'Water-leaving Transmittance calculation, criterion not in range [0.0000001 to 0.01]'
            ACTIONS(NM)    = 'Set input TF_CRITERION in range [0.0000001,0.01]'
            STATUS = VLIDORT_SERIOUS
          ENDIF
        ENDIF

      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  05 March 2013  Start %%%%%%%%%%%%%%%%%

!   New check to make sure output is MS-only when DO_FULLRAD_MODE is off
!     -- Do not want Internal FO Corrections in this case.
!     -- Turn off Internal FO Corrections, with Warning.

!  Old code, pre 2.8
!      IF ( .not.DO_SS_EXTERNAL .and. .not.DO_FULLRAD_MODE ) then
!        IF ( DO_SSCORR_NADIR .or. DO_SSCORR_OUTGOING ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Internal FO calculation must be Off when MS-only is desired'
!          ACTIONS(NM)  = 'DO_SSCORR_NADIR/DO_SSCORR_OUTGOING flags turned off Internally'
!          STATUS = VLIDORT_WARNING
!          DO_SSCORR_NADIR    = .false.
!          DO_SSCORR_OUTGOING = .false.
!        ENDIF
!      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  05 March 2013 Finish %%%%%%%%%%%%%%%%

!   R. Spurr, RT SOLUTIONS Inc.  Revised for Version 2.8
!mick mod 9/19/2017 - turned off now to allow tightening of some VLIDORT input
!                     definitions related to type of rad solution VLIDORT returns in
!                     different cases

      !IF ( DO_FOCORR .and. .not.DO_FULLRAD_MODE ) then
      !   NM = NM + 1
      !   MESSAGES(NM) = 'Internal FO calculation must be off when MS-only is desired'
      !   ACTIONS(NM)  = 'DO_FOCORR/DO_FOCORR_NADIR/DO_FOCORR_OUTGOING flags turned off Internally'
      !   STATUS = VLIDORT_WARNING
      !   DO_FOCORR          = .false.
      !   DO_FOCORR_NADIR    = .false.
      !   DO_FOCORR_OUTGOING = .false.
      !ENDIF

!  New 15 March 2012
!    Turn off FOCORR flags if the external FO calculation applies

      IF ( DO_FOCORR_EXTERNAL ) then
        IF ( DO_FOCORR_NADIR ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'External FO calculation: Cannot have Nadir single scatter correction'
          ACTIONS(NM)  = 'Turn off DO_FOCORR_NADIR flag'
          STATUS = VLIDORT_SERIOUS
        ENDIF
        IF ( DO_FOCORR_OUTGOING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'External FO calculation: Cannot have Outgoing single scatter correction'
          ACTIONS(NM)  = 'Turn off DO_FOCORR_OUTGOING flag'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  New 02 Jul 2013. Check removed, 19 March 2015, with advent of Lattice option
!  Check dedicated FO code / ObsGeo mode compatibility

!      IF (DO_FO_CALC) THEN
!        IF ( .NOT. DO_OBSERVATION_GEOMETRY ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Internal SS calculation: Must set Observation Geometry option'
!          ACTIONS(NM)  = 'Turn on DO_OBSERVATION_GEOMETRY flag'
!          STATUS = VLIDORT_SERIOUS
!        ENDIF
!      ENDIF

!  Check beam mode operation

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_PLANE_PARALLEL ) THEN
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: plane-parallel and refractive flags both set'
         ACTIONS(NM)  = 'Warning: turn off Refraction internally'
         STATUS       = VLIDORT_WARNING
         DO_REFRACTIVE_GEOMETRY = .FALSE.
        ENDIF
       ENDIF
      ENDIF

!  Check consistency of mean value input control
!  ---------------------------------------------

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Cannot have both mean-value flags set'
          ACTIONS(NM)  = 'Warning: disable DO_MVOUT_ONLY flag internally'
          STATUS = VLIDORT_WARNING
          DO_MVOUT_ONLY = .FALSE.
        ENDIF
      ENDIF

!  Removed DO_NO_AZIMUTH check. 17 January 2006

      IF ( .NOT.DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          IF ( DO_USER_VZANGLES ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Mean-value option needs quadratures only'
            ACTIONS(NM)  = 'Warning: DO_USER_VZANGLES flag disabled internally'
            STATUS = VLIDORT_WARNING
            DO_USER_VZANGLES = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  1/31/21. Version 2.8.3.  Use of MSST flags. Section added check consistency
!     1. TOA Upwelling or Downwelling only (not both), Observational Geometry
!     2. Full-radiance mode, FOCORR Outgoing, no external FO

!   ---- Specialist option. Originally Removed 30 March 2007.
!   ---- Reinstated and expanded check, under a different name (DO_MSSTS) 11/20/19
!   ---- Check revised to include possibility of Downwelling 12/18/19

!  here is the old removed code from 2007...............................
!      IF ( SAVE_LAYER_MSST ) THEN
!        IF ( .NOT. DO_USER_VZANGLES ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Bad input: MSCAT. source term - USER_VZANGLES flag not set'
!          ACTIONS(NM)  = 'Check DO_USER_VZANGLES and SAVE_LAYER_MSST flags'
!          STATUS = VLIDORT_SERIOUS
!        ENDIF
!      ENDIF

!  here is the new code (11/20/19, revised 12/18/19)

      IF ( DO_MSSTS ) THEN

!  Either Upwelling or downwelling (not both), and with FullRad Mode

        IF ( DO_UPWELLING .and. DO_DNWELLING ) THEN
          NM = NM + 1 ;STATUS = VLIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: Upwelling and downwelling flags both set for the MSSTS option'
          ACTIONS(NM)  = 'Check DO_UPWELLING and DO_DNWELLING flags, turn off one of them!!'
        ELSE IF ( ( DO_UPWELLING .or. DO_DNWELLING ) .and..not.DO_FULLRAD_MODE ) THEN
          NM = NM + 1 ;STATUS = VLIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: FullRad_mode flag MUST BE set for the MSSTS option'
          ACTIONS(NM)  = 'Check DO_FULLRAD_MODE flag'
        ENDIF

!  Only 1 output level, either TOA (upwelling) or BOA (downwelling)

        IF ( N_USER_LEVELS.ne.1 ) THEN
          NM = NM + 1 ;STATUS = VLIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: single TOA or BOA output only for the MSSTS option'
          ACTIONS(NM)  = 'Check N_USER_LEVELS, and re-set it to 1'
        ELSE
          IF ( DO_UPWELLING .and. (USER_LEVELS(1).ne.zero) ) THEN
            NM = NM + 1 ;STATUS = VLIDORT_SERIOUS
            MESSAGES(NM) = 'Bad input: single Upwelling TOA output only for the MSSTS option'
            ACTIONS(NM)  = 'Check USER_LEVELS(1) input, must be set to 0.0'
          ELSE IF ( DO_DNWELLING .and. (USER_LEVELS(1).ne.DBLE(NLAYERS)) ) THEN
            NM = NM + 1 ;STATUS = VLIDORT_SERIOUS
            MESSAGES(NM) = 'Bad input: single Downwelling BOA output only for the MSSTS option'
            ACTIONS(NM)  = 'Check USER_LEVELS(1) input, must be set to REAL(NLAYERS,fpk)'
          ENDIF
        ENDIF

!  Check on use of MSSTS, must have DO_OBSERVATION_GEOMETRY and FOCORR_OUTGOING flags set

        IF ( .NOT.DO_OBSERVATION_GEOMETRY ) THEN
          NM = NM + 1 ;STATUS = VLIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: Observational geometry not set correctly for the MSSTS option'
          ACTIONS(NM)  = 'Check DO_OBSERVATION_GEOMETRY flag'
        ENDIF

        IF ( .NOT. DO_FOCORR .or. (DO_FOCORR.and..not.DO_FOCORR_OUTGOING) ) THEN
          NM = NM + 1 ; STATUS = VLIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: FOCORR OUTGOING flags must be set for MSSTS option'
          ACTIONS(NM)  = 'Reset DO_FOCORR and/or DO_FOCORR_OUTGOING flags'
        ENDIF

      ENDIF

!  1/31/21. Version 2.8.3.  End of section checking MSST usage

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Check consistency of BVP_TELESCOPING and SOLUTION_SAVING flags
!  ---Warning. Set solution-saving internally

      IF (DO_BVP_TELESCOPING .AND. .NOT.DO_SOLUTION_SAVING) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: BVP telescoping -> solution saving must be set'
        ACTIONS(NM)  = 'Warning:  Solution saving was set internally'
        STATUS = VLIDORT_WARNING
        DO_SOLUTION_SAVING = .TRUE.
      ENDIF

!  Check consistency of Rayleigh-only and Isotropic-only cases
!   Removed 17 January 2006, Isotropic option removed.

!      IF ( DO_RAYLEIGH_ONLY .AND. DO_ISOTROPIC_ONLY ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'Bad input: Isotropic_only & Rayleigh-only flags both set'
!        ACTIONS(NM)  = 'Check DO_RAYLEIGH_ONLY and DO_ISOTROPIC_ONLY flags'
!        STATUS = VLIDORT_SERIOUS
!      ENDIF

!  -----------Note the following in the scalar code -------------------

!  No Delta-M scaling with Rayleigh only
!   ---Warning. Turn off delta-M scaling.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: No delta-M scaling with Rayleigh-only'
          ACTIONS(NM)  = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS = VLIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  No Delta-M scaling with Isotropic only
!   ---Warning. Turn off delta-M scaling.
!   Removed 17 January 2006, Isotropic option removed.

!      IF ( DO_ISOTROPIC_ONLY ) THEN
!        IF ( DO_DELTAM_SCALING ) THEN
!         NM = NM + 1
!          MESSAGES(NM) = 'Bad input: No delta-M scaling with Isotropic-only'
!          ACTIONS(NM)  = 'Warning: DO_DELTAM_SCALING turned off internally'
!          STATUS = VLIDORT_WARNING
!          DO_DELTAM_SCALING = .FALSE.
!        ENDIF
!      ENDIF

!  Check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: no directional input is set'
        ACTIONS(NM)  = 'Check DO_UPWELLING & DO_DNWELLING: one must be set!'
        STATUS = VLIDORT_SERIOUS
      ENDIF

!  Check number of input expansion coefficient moments (non-Rayleigh)
!  ------------------------------------------------------------------

!  Checks for general scattering case
!    Isotropic part Removed 17 January 2006, Isotropic option removed.

      IF ( .NOT.DO_RAYLEIGH_ONLY .AND. .NOT. DO_FOCORR_ALONE ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          IF ( NGREEK_MOMENTS_INPUT.LT.2*NSTREAMS ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Fewer than 2N expansion moments with delta-M'
            ACTIONS(NM)  = 'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N internally'
            STATUS = VLIDORT_WARNING
            NGREEK_MOMENTS_INPUT = 2*NSTREAMS
          ENDIF
        ELSE
!          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN   !  Changed, Version 2.8
          IF ( DO_FOCORR ) THEN
            IF ( NGREEK_MOMENTS_INPUT.LT.2*NSTREAMS-1 ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Fewer than 2N-1 expansion moments without delta-M'
              ACTIONS(NM)  = 'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N-1 internally'
              STATUS = VLIDORT_WARNING
              NGREEK_MOMENTS_INPUT = 2*NSTREAMS - 1
            ENDIF
          ENDIF
        ENDIF

      ELSE

!  Checks for Rayleigh-only option
!   All warnings.

        IF ( DO_RAYLEIGH_ONLY ) THEN

          IF ( NGREEK_MOMENTS_INPUT.NE.2 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Rayleigh-only, expansion momemts NOT = 2'
            ACTIONS(NM)  = 'Warning: Set NGREEK_MOMENTS_INPUT = 2 internally'
            STATUS       = VLIDORT_WARNING
            NGREEK_MOMENTS_INPUT = 2
          ENDIF

          IF ( DO_BVP_TELESCOPING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Bvp telescoping not possible, Rayleigh only'
            ACTIONS(NM)  = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS       = VLIDORT_WARNING
            DO_BVP_TELESCOPING = .FALSE.
          ENDIF

          IF ( DO_SOLUTION_SAVING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Solution saving not possible, Rayleigh only'
            ACTIONS(NM)  = 'Warning: Turn off SOLUTION_SAVING internally'
            STATUS       = VLIDORT_WARNING
            DO_SOLUTION_SAVING = .FALSE.
          ENDIF

        ENDIF

!  Checks for Isotropic only option. All removed, 17 January 2006.

      ENDIF

!  Reset solution saving and BVP telescoping flags
!  Do not need the isotropic-only options

!      DO_SOLUTION_SAVING =  ( DO_SOLUTION_SAVING .AND. &
!          ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )
!      DO_BVP_TELESCOPING =  ( DO_BVP_TELESCOPING .AND. &
!          ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )

      DO_SOLUTION_SAVING =  ( DO_SOLUTION_SAVING .AND. .NOT.DO_RAYLEIGH_ONLY )
      DO_BVP_TELESCOPING =  ( DO_BVP_TELESCOPING .AND. .NOT.DO_RAYLEIGH_ONLY )

!  BVP telescoping doesn't work with non-Lambertian surfaces
!   Condition relaxed, Version 2.8, now possible.
!      IF (  DO_BVP_TELESCOPING ) THEN
!        IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'BVP telescoping disabled for non-Lambertian'
!          ACTIONS(NM)  = 'Turn off DO_BVP_TELESCOPING flag, done internally'
!          STATUS = VLIDORT_WARNING
!          DO_BVP_TELESCOPING = .FALSE.
!        ENDIF
!      ENDIF

!  Check azimuth-only conditions
!  -----------------------------

!  Check no-Azimuth flag. Now set internally
!    ---WARNING. Do-no-Azimuth Flag turned on
!      IF ( .NOT.DO_NO_AZIMUTH ) THEN
!        IF ( DO_USER_VZANGLES. AND. N_USER_VZANGLES.EQ.1 ) THEN
!          IF ( USER_ANGLES_INPUT(1) .EQ. ZERO ) THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Bad input: zenith-sky output requires no azimuth'
!            ACTIONS(NM)  = 'Warning: DO_NO_AZIMUTH flag set true internally'
!            STATUS = VLIDORT_WARNING
!          ENDIF
!        ENDIF
!      ENDIF

!  Check: OLD single scattering correction and Do Rayleigh
!    ---WARNING. SS Flag turned off
!  Check only required for the diffuse field calculations (Version 2.3)

!  @@@ Rob Fix 5/20/13. NO turn-off DO_SSCORR_NADIR for Rayleigh-BRDF situation
!     if DO_RAYLEIGH_ONLY.and..NOT.DO_LAMBERTIAN_SURFACE and.not.DO_SS_EXTERNAL,
!     then one of the SSCORR flags should be set internally
!  Old Code -------
!      IF ( DO_SSCORR_NADIR ) THEN
!        IF ( DO_RAYLEIGH_ONLY .AND. .NOT. DO_SSFULL ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Bad input: No SS correction for Rayleigh only'
!          ACTIONS(NM)  = 'Warning: DO_SSCORR_NADIR turned off internally'
!          STATUS = VLIDORT_WARNING
!          DO_SSCORR_NADIR = .FALSE.
!        ENDIF
!      ENDIF
!  New code --------  ! Rob Fix 8/18/15
!      IF ( DO_RAYLEIGH_ONLY .AND. .NOT.DO_LAMBERTIAN_SURFACE ) THEN
!        IF ( .NOT.DO_SSCORR_OUTGOING .OR. .NOT.DO_SSCORR_NADIR ) THEN
!          IF ( .NOT.DO_SS_EXTERNAL ) THEN
!            if ( .not.DO_SSCORR_OUTGOING ) then
!              NM = NM + 1
!              MESSAGES(NM) = 'Bad input: Rayleigh-only+BRDF: set an SSCORR flag!'
!              ACTIONS(NM)  = 'Warning: DO_SSCORR_NADIR turned on internally'
!              STATUS = VLIDORT_WARNING
!              DO_SSCORR_NADIR = .true.   ! Rob Fix 8/18/15
!            endif
!          ENDIF
!        ENDIF
!      ENDIF
!  @@@ Rob Fix 5/20/13. End of Fix ===================================

!  Version 2.8, Fix 9/19/16. 3/1/17.

      IF ( DO_RAYLEIGH_ONLY .AND. .NOT.DO_LAMBERTIAN_SURFACE ) THEN
        IF ( .NOT.DO_FOCORR .AND. .NOT.DO_FOCORR_EXTERNAL ) THEN
           NM = NM + 1
           MESSAGES(NM) = 'Bad input: Rayleigh-only+BRDF: Need to set FOCORR flags!'
           ACTIONS(NM)  = 'Warning: DO_FOCORR and DO_FOCORR_NADIR turned on internally'
           STATUS = VLIDORT_WARNING
           DO_FOCORR       = .true. 
           DO_FOCORR_NADIR = .true. 
        ENDIF
      ENDIF

!  Just the single scatter, enabled 25 September 2007. Name changed, Version 2.8
!   Single scatter corrections must be turned on
!   No longer required, Version 2.8
!      IF ( DO_SSCORR_ALONE ) THEN
!        IF ( .NOT.DO_SSCORR_NADIR .AND. .NOT.DO_SSCORR_OUTGOING ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Bad input: Full SS, must have one SSCORR flag set'
!          ACTIONS(NM)  = 'Full SS: default to use outgoing SS correction'
!          STATUS = VLIDORT_WARNING
!          DO_SSCORR_NADIR    = .FALSE.
!          DO_SSCORR_OUTGOING = .TRUE.
!        ENDIF
!      ENDIF

!  Just the FO correction alone, enabled 25 September 2007. Name changed, Version 2.8
!   Diffuse-field Delta-M scaling must be turned off

      IF ( DO_FOCORR_ALONE ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: First-Order ALONE --> diffuse-field delta-M on'
          ACTIONS(NM)  = 'First-Order ALONE: internal default to deltam_scaling = false'
          STATUS = VLIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  Check thermal inputs
!  --------------------

!  If thermal transmittance only, check thermal flag

      IF ( DO_SOLAR_SOURCES .or. .NOT. DO_THERMAL_EMISSION ) THEN
       IF ( DO_THERMAL_TRANSONLY ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: No thermal, must turn off transonly flag'
         ACTIONS(NM)  = 'Warning: DO_THERMAL_TRANSONLY turned off internally'
         STATUS = VLIDORT_WARNING
         DO_THERMAL_TRANSONLY = .FALSE.
       ENDIF
      ENDIF

!  Switch off a bunch of flags
!    Solution saving is required though!

      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( DO_THERMAL_TRANSONLY ) THEN
          !DO_RAYLEIGH_ONLY  = .FALSE.
          DO_DELTAM_SCALING  = .FALSE.
          DO_DELTAM_SCALING  = .FALSE.
          DO_SOLUTION_SAVING = .TRUE.
        ENDIF
      ENDIF

!  No solar sources for thermal transmittance

      IF ( DO_THERMAL_TRANSONLY ) THEN
        IF ( DO_SOLAR_SOURCES ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: thermal tranmsittance, must turn off solar'
         ACTIONS(NM)  = 'Warning: DO_SOLAR_SOURCES turned off internally'
         STATUS = VLIDORT_WARNING
         DO_SOLAR_SOURCES = .FALSE.
       ENDIF
      ENDIF

!  Check viewing geometry input
!  ----------------------------

!  Check Earth radius (Chapman function only)
!    ---WARNING. Default value of 6371.0 will be set

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          IF ( EARTH_RADIUS.LT.6320.0D0 .OR. &
               EARTH_RADIUS.GT.6420.0D0 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Earth radius outside of [6320-6420]'
            ACTIONS(NM)  = 'Warning: default value of 6371.0 was set'
            STATUS = VLIDORT_WARNING
            EARTH_RADIUS = 6371.0D0
          ENDIF
        ENDIF
      ENDIF

!  Check dimensioning on Legendre numbers (refractive geometry only)

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        NALLSTREAMS = N_SZANGLES*NLAYERS + NSTREAMS + N_USER_VZANGLES
        IF ( NALLSTREAMS .GT. MAX_ALLSTRMS_P1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Dimensioning error for refractive beam angles'
          ACTIONS(NM)  = 'Increase dimension MAX_ALLSTRMS_P1'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Check GEOMETRY_SPECHEIGHT (only for outgoing sphericity correction)
!    GEOMETRY_SPECHEIGHT cannot be greater than HEIGHT_GRID(NLAYERS)

      IF ( DO_FOCORR .and. DO_FOCORR_OUTGOING ) THEN
        IF ( GEOMETRY_SPECHEIGHT .GT. HEIGHT_GRID(NLAYERS) ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'GEOMETRY_SPECHEIGHT must be =< Input BOA-HEIGHT '
          ACTIONS(NM)  = 'Warning: Internal Re-set of GEOMETRY_SPECHEIGHT '
          STATUS = VLIDORT_WARNING
          GEOMETRY_SPECHEIGHT  = HEIGHT_GRID(NLAYERS)
        ENDIF
      ENDIF

!  Check solar zenith angle input

!mick hold - 9/26/2012
      !IF (DO_SOLAR_SOURCES ) THEN
        DO I = 1, N_SZANGLES
          IF ( SZANGLES(I) .LT. ZERO .OR. &
               SZANGLES(I) .GE. 90.0D0 ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: out-of-range solar angle, no. '//C2
            ACTIONS(NM)  = 'Look at SZANGLES input, should be < 90 & > 0'
            STATUS = VLIDORT_SERIOUS
          ENDIF
        ENDDO
      !ELSE IF ( .NOT.DO_SOLAR_SOURCES .AND. DO_THERMAL_EMISSION ) THEN
      !  SZANGLES(1:N_SZANGLES) = ZERO
      !ENDIF

!  Check relative azimuths

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I) .GT. 360.0D0   .OR. &
             USER_RELAZMS(I) .LT. ZERO ) THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: out-of-range azimuth angle, no. '//C2
          ACTIONS(NM)  = 'Look at azimuth angle input, should be in [0,360]'
          LOOP = .FALSE.
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Limits on user-defined options

!      IF ( .NOT. DO_USER_VZANGLES .AND..NOT.DO_QUAD_OUTPUT ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'Bad input: No angular stream output is specified'
!        ACTIONS(NM)  = 'Check DO_USER_VZANGLES and DO_QUAD_OUTPUT flags'
!        STATUS = VLIDORT_SERIOUS
!      ENDIF

!  Check user-defined stream angles (should always be [0,90])

      IF ( DO_USER_VZANGLES ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.N_USER_VZANGLES)
          I = I + 1
          IF ( USER_VZANGLES(I) .GT. 90.0D0   .OR. &
               USER_VZANGLES(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: out-of-range user stream, no. '//C2
            ACTIONS(NM)  = 'Look at user viewing zenith angle input'
            LOOP         = .FALSE.
            STATUS       = VLIDORT_SERIOUS
          ENDIF
        ENDDO
      ENDIF

!  Re-order the input angles
!  -------------------------

!  This section has been truncated. 28 March 2007

!mick fix - moved "N_OUT_STREAMS = N_USER_VZANGLES" outside if block
!           & initialized OUT_ANGLES to ZERO in case DO_USER_VZANGLES = .FALSE.
      N_OUT_STREAMS  = N_USER_VZANGLES
      OUT_ANGLES     = ZERO
      IF ( DO_USER_VZANGLES ) THEN
        IF ( N_OUT_STREAMS .EQ. 1 ) THEN
          OUT_ANGLES(1) =  USER_VZANGLES(1)
        ELSE
!%%%%%%%%%%%%%%%% ROB clause 05 March 2013 %%%%%%%%%%%%%%%%%
!%%%%%% Re-ordering only for Lattice computations %%%%%%%%%%

!  revised Call, 3/17/17, with RQSORT_IDX

           IF ( .NOT.DO_OBSERVATION_GEOMETRY ) THEN
              DO I = 1, N_USER_VZANGLES
                ALL_ANGLES(I)   = USER_VZANGLES(I)
                INDEX_ANGLES(I) = I                 ! 3/17/17  Need to initialize
              ENDDO
              CALL RQSORT_IDX(N_OUT_STREAMS,USER_VZANGLES,N_OUT_STREAMS,INDEX_ANGLES)
              DO I = 1, N_OUT_STREAMS
                OUT_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
                USER_VZANGLES(I) = OUT_ANGLES(I)
              ENDDO
           ELSE
              DO I = 1, N_OUT_STREAMS
                OUT_ANGLES(I) = USER_VZANGLES(I)
              ENDDO
           ENDIF

!%%%%%%%%%%%%%%%%%%%%%%% End ROB clause 05 March 2013 %%%%%%%%%%%%%
        ENDIF
      ENDIF

!  Check height grid input (Chapman function only)

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.NLAYERS)
          I = I + 1
          IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Height-grid not monotonically decreasing; Layer '//C2
            ACTIONS(NM)  = 'Look at Height-grid input'
            LOOP = .FALSE.
            STATUS = VLIDORT_SERIOUS
          ENDIF
        ENDDO
      ENDIF

!  Check vertical outputs
!  ----------------------

!  Check vertical output levels (should always be within atmosphere!)

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_LEVELS)
        I = I + 1
        IF ( USER_LEVELS(I) .GT. DBLE(NLAYERS) .OR. &
             USER_LEVELS(I) .LT. ZERO )  THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Out of range for level choice # '//C2
          ACTIONS(NM)  = 'Re-set level output '
          LOOP = .FALSE.
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check repetition of vertical output choices

      UTA = 0
      LOOP = .TRUE.
      DO WHILE ( LOOP .AND. UTA .LT. N_USER_LEVELS )
        UTA = UTA + 1
        XT = USER_LEVELS(UTA)
        NSTART = 0
        DO N = 1, N_USER_LEVELS
          IF ( XT .EQ. USER_LEVELS(N)) NSTART = NSTART + 1
        ENDDO
        IF ( NSTART .NE. 1 ) THEN
          LOOP = .FALSE.
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: repetition of vertical output choice'
          ACTIONS(NM)  = 'Re-set level output '
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check geophysical scattering inputs
!  -----------------------------------

!  Check single scatter albedos

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
       DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
          NM = NM + 1
          WRITE(C3,'(I3)')L
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS = VLIDORT_SERIOUS
        ELSE IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
          NM = NM + 1
          WRITE(C3,'(I3)')L
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS = VLIDORT_SERIOUS
        ENDIF
       ENDDO
      ENDIF

!  Solar beam, cannot be too small
!mick mod 9/19/2017 - turned this check off (one above encompasses this one)

      !IF ( DO_SOLAR_SOURCES ) THEN
      !  DO L = 1, NLAYERS
      !   IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
      !    WRITE(C3,'(I3)')L
      !    NM = NM + 1
      !    MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
      !    ACTIONS(NM)  = 'Check SS-albedo input'
      !    STATUS = VLIDORT_SERIOUS
      !   ENDIF
      !  ENDDO
      !ENDIF

!  Check first phase function moments
!    Version 2.8, Upgrade, 3/1/17. if not using the Greekmat coefficients in FOCORR

      IF ( ( .NOT.DO_THERMAL_TRANSONLY ) .OR. .NOT.( DO_FOCORR.AND.DO_SSCORR_USEFMAT ) ) THEN
       DO L = 1, NLAYERS
        IF ( GREEKMAT_TOTAL_INPUT(0,L,1).NE.ONE ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'First phase moment (GREEK_11) not 1 for layer '//C3
          ACTIONS(NM)  = 'Check First phase function moment'
          STATUS = VLIDORT_SERIOUS
        ENDIF
       ENDDO
      ENDIF

!  Additional check on non-negativity of (1,1) moments.
!    Furnished 27 September 2007, as a result of input-checking.
!      DO N = 1, NLAYERS
!        DO L = 1, NGREEK_MOMENTS_INPUT
!          IF ( GREEKMAT_TOTAL_INPUT(L,N,1).LT.ZERO ) THEN
!            WRITE(C3,'(I3)')N
!            NM = NM + 1
!            MESSAGES(NM) = 'Some moments (GREEK_11) are NEGATIVE for layer '//C3
!            ACTIONS(NM)  = 'Check Greek moments input - some bad values!'
!            STATUS = VLIDORT_SERIOUS
!          ENDIF
!        ENDDO
!      ENDDO

!  Specialist options, Check the given input

      IF ( DO_SPECIALIST_OPTION_2 ) THEN
        IF ( NLAYERS_NOMS .GT. NLAYERS-1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input specialist option 2: NLAYERS_NOMS'
          ACTIONS(NM)  = 'Check NLAYERS_NOMS must be less than NLAYERS'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_SPECIALIST_OPTION_3 ) THEN
        IF ( NLAYERS_CUTOFF .lt. 2 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input specialist option 3: NLAYERS_CUTOFF'
          ACTIONS(NM)  = 'Check NLAYERS_CUTOFF must be greater than 1'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  TOA contributions
!    TOA Upwelling must be set, if you are using this flag

      IF ( DO_TOA_CONTRIBS ) THEN
        IF ( .not. DO_UPWELLING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input TOA contributions, Upwelling Not set'
          ACTIONS(NM)  = 'Must set the DO_UPWELLING flag'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_TOA_CONTRIBS ) THEN
        N = 0
        DO UTA = 1, N_USER_LEVELS
          if ( USER_LEVELS(UTA) .ne. 0.0d0 ) N = N + 1
        ENDDO
        IF ( N .EQ. N_USER_LEVELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input TOA contributions, Level TOA not set'
          ACTIONS(NM)  = 'Must set one level output to be TOA (0.0)'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  1/31/21/ Version 2.8.3. New section, Checking on Doublet Geometry

      IF ( DO_DOUBLET_GEOMETRY ) THEN
        IF ( DO_TOA_CONTRIBS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input Doublet Geometry: TOA contributions not allowed '
          ACTIONS(NM)  = 'Turn off TOA_CONTRIBS flag'
          STATUS = VLIDORT_SERIOUS
        ENDIF
        IF ( DO_THERMAL_EMISSION ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input Doublet Geometry: No thermal emission with Doublet Geometry'
          ACTIONS(NM)  = 'Turn off DO_THERMAL_EMISSION flags'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!mick fix - define NMESSAGES
!  Number of messages

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHECK_INPUT

!

      SUBROUTINE VLIDORT_CHECK_INPUT_OPTICAL &
          ( NLAYERS, LAMBERTIAN_ALBEDO, DELTAU_VERT_INPUT,  & ! Input
            OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,        & ! Input
            STATUS, NMESSAGES, MESSAGES, ACTIONS )            ! output

!  Check the threaded optical property inputs

!  Module, dimensions and numbers

      USE VLIDORT_pars_m, only : fpk, MAXLAYERS, MAXMOMENTS_INPUT, MAX_MESSAGES, &
                                 MAXSTOKES_SQ, VLIDORT_SUCCESS, VLIDORT_SERIOUS, &
                                 ZERO, ONE, OMEGA_SMALLNUM

      IMPLICIT NONE

!  Module input
!  ------------

      INTEGER  , intent(in)  :: NLAYERS

!  Multilayer optical property (bulk) inputs

      REAL(fpk), intent(in)  :: OMEGA_TOTAL_INPUT  ( MAXLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_VERT_INPUT  ( MAXLAYERS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk), intent(in)  :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Lambertian Surface control

      REAL(fpk), intent(in)  :: LAMBERTIAN_ALBEDO

!  Module output
!  -------------

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER          :: L, NM
      CHARACTER*3      :: C3

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of VLIDORT Optical Input'
!      ACTIONS(0)      = 'No Action required for this Task'

      NM = NMESSAGES

!  Check Thread-dependent optical property inputs
!  ----------------------------------------------

!  Make sure the Lambertian surface is in range

      IF ( LAMBERTIAN_ALBEDO.LT.ZERO .OR. LAMBERTIAN_ALBEDO.GT.ONE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Lambertian albedo not in range [0,1]'
        ACTIONS(NM)  = 'Check albedo input'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

!  Check optical thickness values

      DO L = 1, NLAYERS
        IF ( DELTAU_VERT_INPUT(L).LE.ZERO ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: optical thickness <= 0, layer '//C3
          ACTIONS(NM)  = 'Check optical thickness input'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check single scatter albedos

      DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Solar beam, cannot be too small

      DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check first phase function moments

      DO L = 1, NLAYERS
        IF ( GREEKMAT_TOTAL_INPUT(0,L,1).NE.ONE ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'First Greek matrix moment not 1 for layer '//C3
          ACTIONS(NM)  = 'Check First Greek matrix moment input'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Number of messages

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHECK_INPUT_OPTICAL

!

      SUBROUTINE VLIDORT_DERIVE_INPUT ( &
        DO_FULLRAD_MODE, DO_UPWELLING, DO_DNWELLING, DO_FOCORR,                  & ! Input Boolean
        DO_FOCORR_ALONE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_SSCORR_USEFMAT, & ! Input Boolean
        DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_RAYLEIGH_ONLY,              & ! Input Boolean
        DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY,           & ! Input Boolean
        DO_THERMAL_TRANSONLY, DO_SPECIALIST_OPTION_2, DO_SOLUTION_SAVING,      & ! Input Boolean
        DO_BVP_TELESCOPING, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_NO_AZIMUTH, & ! Input Boolean
        DO_DOUBLE_CONVTEST, DO_ALL_FOURIER, DO_DIRECT_BEAM, DO_DBCORRECTION,   & ! Input Boolean
        NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SZANGLES,  NLAYERS_NOMS,  & ! Input Integer
        NGREEK_MOMENTS_INPUT, N_USER_RELAZMS, N_USER_VZANGLES,  N_OUT_STREAMS, & ! Input Integer
        SZANGLES, USER_VZANGLES, USER_LEVELS,                                  & ! Input Floating point
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,            & ! Input Floating point
        COS_SZANGLES, SIN_SZANGLES, QUAD_STREAMS,                              & ! OUTPUT
        QUAD_WEIGHTS,QUAD_STRMWTS, QUAD_HALFWTS, QUAD_SINES, QUAD_ANGLES,      & ! OUTPUT
        DO_MSMODE_VLIDORT, NMOMENTS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, & ! OUTPUT
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, NBEAMS, NPARTICSOLS, NSTOKES_SQ,         & ! OUTPUT
        FLUXVEC, DMAT, MUELLER_INDEX, GREEKMAT_INDEX, DO_REAL_EIGENSOLVER,     & ! OUTPUT
        BVP_REGULAR_FLAG, LAYER_MAXMOMENTS, DO_LAYER_SCATTERING,               & ! OUTPUT
        N_CONVTESTS, N_CONV_STREAMS, N_DIRECTIONS, WHICH_DIRECTIONS,           & ! OUTPUT
        USER_STREAMS, USER_SINES, USER_SECANTS, PARTLAYERS_OUTFLAG,            & ! OUTPUT
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,           & ! OUTPUT
        DO_PARTLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & ! OUTPUT
        N_LAYERSOURCE_UP, N_LAYERSOURCE_DN, N_ALLLAYERS_UP, N_ALLLAYERS_DN,    & ! OUTPUT
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, TAUGRID_INPUT,          & ! OUTPUT
        STATUS_SUB, MESSAGE )

!  1/31/21. Version 2.8.3. Introduce DO_DOUBLET_GEOMETRY input
!   -- re-ordered first 4 lines of input
!   -- DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY are not needed.

      USE VLIDORT_PARS_m, Only : MAXMOMENTS_INPUT, MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXLAYERS,      &
                                 MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_LEVELS, MAX_USER_STREAMS,  &
                                 MAX_PARTLAYERS, MAX_DIRECTIONS, MAX_PSOLS, MAXSTOKES_SQ,             &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_WARNING,                   &
                                 ZERO, ONE, HALF, THREE, DEG_TO_RAD, UPIDX, DNIDX

!  Revised Call, 3/17/17
!mick fix 9/19/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input
!                   - removed SZA_LOCAL_INPUT & SUN_SZA_COSINES from I/O

      USE VLIDORT_AUX_m, Only : GETQUAD2, RSSORT

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::             DO_FULLRAD_MODE

      LOGICAL, INTENT (IN) ::             DO_FOCORR
      LOGICAL, INTENT (IN) ::             DO_FOCORR_ALONE
!mick fix 9/19/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input
      LOGICAL, INTENT (IN) ::             DO_FOCORR_NADIR
      LOGICAL, INTENT (IN) ::             DO_FOCORR_OUTGOING
      LOGICAL, INTENT (IN) ::             DO_SSCORR_USEFMAT

      LOGICAL, INTENT (IN) ::             DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::             DO_REFRACTIVE_GEOMETRY

      LOGICAL, INTENT (IN) ::             DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::             DO_UPWELLING
      LOGICAL, INTENT (IN) ::             DO_DNWELLING

!  1/31/21. Version 2.8.3. Introduce DO_DOUBLET_GEOMETRY input

      LOGICAL, INTENT (IN) ::             DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::             DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::             DO_DOUBLET_GEOMETRY
      LOGICAL, INTENT (IN) ::             DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::             DO_SPECIALIST_OPTION_2

      LOGICAL, INTENT (INOUT) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT (INOUT) ::          DO_BVP_TELESCOPING
      LOGICAL, INTENT (INOUT) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (INOUT) ::          DO_MVOUT_ONLY

      LOGICAL, INTENT (INOUT) ::          DO_DOUBLE_CONVTEST
      LOGICAL, INTENT (INOUT) ::          DO_ALL_FOURIER
      LOGICAL, INTENT (INOUT) ::          DO_DIRECT_BEAM
      LOGICAL, INTENT (INOUT) ::          DO_DBCORRECTION
      LOGICAL, INTENT (INOUT) ::          DO_NO_AZIMUTH

!  Input integers

      INTEGER, INTENT (IN) ::             NSTOKES
      INTEGER, INTENT (IN) ::             NSTREAMS
      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (IN) ::             N_SZANGLES

      INTEGER, INTENT (IN) ::             N_USER_RELAZMS
      INTEGER, INTENT (IN) ::             N_USER_VZANGLES
      INTEGER, INTENT (IN) ::             N_USER_LEVELS

      INTEGER, INTENT (IN) ::             NLAYERS_NOMS
      INTEGER, INTENT (IN) ::             N_OUT_STREAMS

!  Input real variables. Adjusted angles removed, version 2.8

      DOUBLE PRECISION, INTENT (IN) ::    SZANGLES ( MAX_SZANGLES )
!      DOUBLE PRECISION, INTENT (IN) ::    SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::    USER_VZANGLES ( MAX_USER_VZANGLES )
!      DOUBLE PRECISION, INTENT (IN) ::    USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )

      DOUBLE PRECISION, INTENT (IN) ::    DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT (INOUT) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: USER_LEVELS ( MAX_USER_LEVELS )

!  OUTPUT

      DOUBLE PRECISION, INTENT (OUT) ::   FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) ::   COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (OUT) ::   SIN_SZANGLES ( MAX_SZANGLES )
!      DOUBLE PRECISION, INTENT (OUT) ::   SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )

      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_SINES   ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   QUAD_ANGLES  ( MAXSTREAMS )

      LOGICAL, INTENT (OUT) ::            DO_MSMODE_VLIDORT
      INTEGER, INTENT (OUT) ::            NMOMENTS
      INTEGER, INTENT (OUT) ::            NSTREAMS_2
      INTEGER, INTENT (OUT) ::            NTOTAL
      INTEGER, INTENT (OUT) ::            N_SUBDIAG
      INTEGER, INTENT (OUT) ::            N_SUPDIAG
      INTEGER, INTENT (OUT) ::            NSTOKES_SQ
      INTEGER, INTENT (OUT) ::            NSTKS_NSTRMS
      INTEGER, INTENT (OUT) ::            NSTKS_NSTRMS_2
      INTEGER, INTENT (OUT) ::            NBEAMS
      INTEGER, INTENT (OUT) ::            NPARTICSOLS

      INTEGER, INTENT (OUT) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) ::   DMAT ( MAXSTOKES, MAXSTOKES )

      INTEGER, INTENT (OUT) ::            GREEKMAT_INDEX ( 6 )
      LOGICAL, INTENT (OUT) ::            DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (OUT) ::            BVP_REGULAR_FLAG ( 0:MAXMOMENTS )
      INTEGER, INTENT (OUT) ::            LAYER_MAXMOMENTS ( MAXLAYERS )
      LOGICAL, INTENT (OUT) ::            DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

      INTEGER, INTENT (OUT) ::            N_CONVTESTS
      INTEGER, INTENT (OUT) ::            N_CONV_STREAMS
      INTEGER, INTENT (OUT) ::            N_DIRECTIONS

      INTEGER, INTENT (OUT) ::            WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (OUT) ::   USER_STREAMS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   USER_SINES    ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) ::   USER_SECANTS  ( MAX_USER_STREAMS )

      LOGICAL, INTENT (OUT) ::            PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (OUT) ::            PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (OUT) ::            UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (OUT) ::            UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (OUT) ::            DO_PARTLAYERS
      INTEGER, INTENT (OUT) ::            N_PARTLAYERS
      INTEGER, INTENT (OUT) ::            PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::   PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

      INTEGER, INTENT (OUT) ::            N_LAYERSOURCE_UP
      INTEGER, INTENT (OUT) ::            N_LAYERSOURCE_DN
      INTEGER, INTENT (OUT) ::            N_ALLLAYERS_UP
      INTEGER, INTENT (OUT) ::            N_ALLLAYERS_DN
      LOGICAL, INTENT (OUT) ::            STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (OUT) ::            STERM_LAYERMASK_DN ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) ::   DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) ::   TAUGRID_INPUT ( 0:MAXLAYERS )

      INTEGER, INTENT (OUT) ::            STATUS_SUB
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE

!  Local variables
!  ---------------

      DOUBLE PRECISION :: DT, RT
      INTEGER ::          I, UT, N, UTA, NSTART, M, NAP, NS
      INTEGER ::          O1, O2, O1_S, MS2, P, NA, Q, QC, L, EC
      LOGICAL ::          LOOP, LOCAL_NADIR_ONLY
      INTEGER ::          CONV_START_STREAMS

      DOUBLE PRECISION :: DMAT_PSOLS ( MAXSTOKES, MAXSTOKES, MAX_PSOLS )
      DOUBLE PRECISION :: DMAT_PSOLS_FLUXVEC ( MAXSTOKES,MAX_PSOLS )

!  Set additional numbers (derived input)
!  ======================

!  Set status

      STATUS_SUB = VLIDORT_SUCCESS
      MESSAGE    = ' '

!  Automatic input

      DO_ALL_FOURIER        = .FALSE.
!      DO_CLASSICAL_SOLUTION = .TRUE.   ! Removed from the Code, Version 2.8, 7/6/16
      DO_DIRECT_BEAM        = .TRUE.

!  FOCORR_ALONE flag cancels some other flags. Not the DB correction !!!!
!    (Revision, version 2.8, new name for FOCORR_ALONE, used to be DO_SSFULL)

      IF ( DO_FOCORR_ALONE ) THEN
        DO_DOUBLE_CONVTEST  = .FALSE.
        DO_SOLUTION_SAVING  = .FALSE.
        DO_BVP_TELESCOPING  = .FALSE.
        DO_ADDITIONAL_MVOUT = .FALSE.
        DO_MVOUT_ONLY       = .FALSE.
      ENDIF

!  Set DB correction flag (this section 11 October 2010)

      DO_DBCORRECTION = .false.
      IF ( DO_FOCORR ) DO_DBCORRECTION = .true.

!  Mode of operation
!   SS outgoing sphericity option, added 31 January 2007
!   SS full calculation option added 25 September 2007.
!  Version 2.8. MSMODE VLIDORT is never True now.
!mick mod 9/19/2017 - modified IF structure to accomodate some modification of some
!                     VLIDORT input definitions related to the type of rad solution VLIDORT
!                     returns in different cases.

      DO_MSMODE_VLIDORT = .FALSE.

!      IF ( DO_FULLRAD_MODE ) THEN
!!        IF ( .NOT. DO_FOCORR_ALONE ) THEN
!!          IF ( DO_FOCORR_NADIR .OR. DO_FOCORR_OUTGOING ) THEN
!!            DO_MSMODE_VLIDORT = .TRUE.
!!          ENDIF
!!        ENDIF
!      ELSE
!        IF ( .NOT. DO_FOCORR_ALONE ) THEN
!          DO_MSMODE_VLIDORT = .TRUE.
!        ENDIF
!      ENDIF

      IF ( DO_FULLRAD_MODE ) THEN
        IF ( DO_FOCORR ) THEN
          !Case: MS trunc + FO corr
          DO_MSMODE_VLIDORT = .TRUE.
        ENDIF
      ELSE
        IF ( .NOT.DO_FOCORR ) THEN
          !Case: MS trunc only
          DO_MSMODE_VLIDORT = .TRUE.
        ENDIF
      ENDIF

!  Directional indices

      IF ( DO_UPWELLING .AND. DO_DNWELLING ) THEN
        N_DIRECTIONS = 2
        WHICH_DIRECTIONS(1) = UPIDX
        WHICH_DIRECTIONS(2) = DNIDX
      ELSE
        N_DIRECTIONS = 1
        WHICH_DIRECTIONS(2) = 0
        IF ( DO_UPWELLING ) THEN
          WHICH_DIRECTIONS(1) = UPIDX
        ELSE IF ( DO_DNWELLING) THEN
          WHICH_DIRECTIONS(1) = DNIDX
        ENDIF
      ENDIF

!  New section. Setting the DO_NO_AZIMUTH flag
!    Rt Solutions. 17 January 2006. R. Spurr and V. Natraj.

!  DO_NO_AZIMUTH should be set internally only when:
!     (a) NSTOKES = 1 and nadir view only, and no SS correction
!     (b) DO_MVOUT_ONLY flag is true.
!     (c) SSFULL calculation flag is set

      LOCAL_NADIR_ONLY = .FALSE.
!mick fix 3/30/2015 - added outer if block
      IF ( DO_USER_STREAMS ) THEN
        IF ( N_USER_VZANGLES.EQ.1 .AND. USER_VZANGLES(1).EQ.ZERO ) THEN
          LOCAL_NADIR_ONLY = .TRUE.
        ENDIF
      ELSE
        LOCAL_NADIR_ONLY = .TRUE.
      END IF

!  this condition simplified for Version 2.8

      DO_NO_AZIMUTH = .FALSE.
      IF ( ( NSTOKES.EQ.1 .AND. LOCAL_NADIR_ONLY .AND. .NOT.DO_FOCORR ) .OR. DO_MVOUT_ONLY ) THEN
        DO_NO_AZIMUTH = .TRUE.
      ENDIF
!      IF ( ( NSTOKES.EQ.1 .AND. LOCAL_NADIR_ONLY .AND. &
!             (.NOT.DO_FOCORR_ALONE .OR. .NOT.DO_FOCORR_NADIR .OR. .NOT.DO_FOCORR_OUTGOING) ) &
!           .OR. DO_MVOUT_ONLY  ) THEN
!        DO_NO_AZIMUTH = .TRUE.
!      ENDIF

!  Flux vector set unity input.

      FLUXVEC(1) = ONE
      FLUXVEC(2) = ZERO
      FLUXVEC(3) = ZERO
      FLUXVEC(4) = ZERO

!  Convert Surface emission input.............. NOT USED....
!  Input values should be Watts/ sq m
!      IF ( DO_SURFACE_EMISSION ) THEN
!        FP_SURFBB = PI4 * SURFBB
!      ENDIF

!  Number of moments. Isotropic option removed 17 January 2006.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = 2
      ENDIF
      IF ( .NOT.DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = MIN ( 2 * NSTREAMS - 1, NGREEK_MOMENTS_INPUT )
      ENDIF

!  Total quadratures (up and down)

      NSTREAMS_2 = 2*NSTREAMS

!  Additional quantities (Stokes)

      NSTOKES_SQ     = NSTOKES * NSTOKES
      NSTKS_NSTRMS   = NSTOKES * NSTREAMS
      NSTKS_NSTRMS_2 = NSTOKES * NSTREAMS_2

!  Mueller index

      DO O1 = 1, MAXSTOKES
        O1_S = MAXSTOKES*(O1 - 1)
        DO O2 = 1, MAXSTOKES
          MUELLER_INDEX(O1,O2) = O1_S + O2
        ENDDO
      ENDDO

!  Greek matrix index

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

!  Check number of particular solution modes
!   Current default is NPARTICSOLS = 1

      IF ( FLUXVEC(3).EQ.ZERO.AND.FLUXVEC(4).EQ.ZERO ) THEN
        NPARTICSOLS = 1
      ELSE
        NPARTICSOLS = 2
      ENDIF

!  D matrices, and multiply by Fluxvector

      DO O1 = 1, MAXSTOKES
        DFLUX(O1) = ZERO
        DO O2 = 1, MAXSTOKES
          DMAT(O1,O2) = ZERO
          DO P = 1, NPARTICSOLS
            DMAT_PSOLS(O1,O2,P) = ZERO
          ENDDO
        ENDDO
      ENDDO
      MS2 = MAXSTOKES/2
      DO O1 = 1, MS2
        O2 = O1 + MS2
        DMAT(O1,O1) = ONE
        DMAT(O2,O2) = -ONE
        DMAT_PSOLS(O1,O1,1) = ONE
        DMAT_PSOLS_FLUXVEC(O1,1) = FLUXVEC(O1)
        DFLUX(O1) = FLUXVEC(O1)
      ENDDO
      IF ( NPARTICSOLS .EQ. 2 ) THEN
        DO O1 = 1, MS2
          O2 = O1 + MS2
          DMAT_PSOLS(O2,O2,2) = ONE
          DMAT_PSOLS_FLUXVEC(O2,1) = FLUXVEC(O2)
          DFLUX(O2) = FLUXVEC(O2)
        ENDDO
      ENDIF

!  Set Quadrature abscissae and weights
!    - REVISED for Version 2.8, 3/17/17

      IF ( NSTREAMS .EQ. 1 ) THEN
        QUAD_STREAMS(1) = SQRT ( ONE / THREE )
        QUAD_WEIGHTS(1) = ONE
      ELSE
        CALL GETQUAD2 ( ZERO, ONE, NSTREAMS, QUAD_STREAMS, QUAD_WEIGHTS )
      ENDIF

!  Set auxiliary quantities

      DO I = 1, NSTREAMS
        QUAD_STRMWTS(I) = QUAD_STREAMS(I)*QUAD_WEIGHTS(I)
        QUAD_HALFWTS(I) = HALF * QUAD_WEIGHTS(I)
        QUAD_ANGLES(I)  = ACOS(QUAD_STREAMS(I))/DEG_TO_RAD
        QUAD_SINES(I)   = SQRT(ONE-QUAD_STREAMS(I)*QUAD_STREAMS(I))
      ENDDO

!  Size of boundary value problem matrices and vectors

      NTOTAL = NLAYERS*NSTKS_NSTRMS_2

!  Number of sub and super diagonals in band matrix (boundary value problem)

      IF ( NLAYERS .EQ. 1 ) THEN
        N_SUBDIAG = 2*NSTKS_NSTRMS - 1
        N_SUPDIAG = 2*NSTKS_NSTRMS - 1
      ELSE
        N_SUBDIAG = 3*NSTKS_NSTRMS - 1
        N_SUPDIAG = 3*NSTKS_NSTRMS - 1
      ENDIF

!  Solar zenith angle cosines/sines

      NBEAMS = N_SZANGLES
      DO I = 1, N_SZANGLES
        COS_SZANGLES(I) = COS(SZANGLES(I) * DEG_TO_RAD )
        SIN_SZANGLES(I) = SQRT(ONE-COS_SZANGLES(I)*COS_SZANGLES(I))
      ENDDO

!mick fix 9/19/2017 - moved calculation of SUN_SZA_COSINES from here
!                     to LEVELS_GEOMETRY_PREPARE to resolve an I/O contradiction
!  Set average cosines in the refractive geometry case

!      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        DO I = 1, N_SZANGLES
!          MU1 = COS(SZA_LOCAL_INPUT(0,I)*DEG_TO_RAD)
!          DO N = 1, NLAYERS
!            MU2 = COS(SZA_LOCAL_INPUT(N,I)*DEG_TO_RAD)
!            SUN_SZA_COSINES(N,I) = HALF * ( MU1 + MU2 )
!            MU1 = MU2
!          ENDDO
!        ENDDO
!      ENDIF

!  Set performance flags
!  ---------------------

!  New section for Version 2.0 and higher.

!  Set the layer scattering flags, set the BVP telescoping flags
!    Specialist Option 2: Set solution saving. BVP Telescoping. Remove 4545 point
!                         Set the layer-scattering flags to be false for all N up to NLAYERS_N

!  Initialise (M = Fourier index) to normal mode.

      DO M = 0, NMOMENTS
        BVP_REGULAR_FLAG(M) = .TRUE.
        DO N = 1, NLAYERS
          DO_LAYER_SCATTERING(M,N) = .TRUE.
        ENDDO
      ENDDO

!  Specialist Option 2: Set solution saving. BVP Telescoping. Remove 4545 point
!   Set the layer-scattering flags to be false for all N up to NLAYERS_N

      IF ( DO_SPECIALIST_OPTION_2 ) THEN
        DO_SOLUTION_SAVING = .TRUE.
        DO_BVP_TELESCOPING = .TRUE.
        DO M = 0, NMOMENTS
          BVP_REGULAR_FLAG(M) = .FALSE.
          DO N = 1, NLAYERS_NOMS
            DO_LAYER_SCATTERING(M,N) = .FALSE.
          ENDDO
          DO N = NLAYERS_NOMS + 1, NLAYERS
            DO_LAYER_SCATTERING(M,N) = .TRUE.
          ENDDO
        ENDDO
      ENDIF

!  Special clause for transmittance-only calculation
!   New flag, 31 July 2007. Move on to other bookkeeping

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO_SOLUTION_SAVING = .TRUE.
        DO N = 1, NLAYERS
          OMEGA_TOTAL_INPUT(N) = ZERO
          DO_LAYER_SCATTERING(0,N) = .FALSE.
        ENDDO
!        GO TO 5557          ! removed Version 2.8
      ENDIF 

!  For M > 2 terms, if solution saving flag is set...
!  .... examine Greek matrix entries (1,1) - no scattering if
!      they are all zero for L > M - 1.
!  [Equivalent to examining phase function moments]
!  Addition of Solution_saving flag, 22 November 2009.
!    Bug spotted by V. Natraj, 20 November 2009

!  No, should be set always, otherwise linearizations fail

      IF ( DO_SOLAR_SOURCES ) THEN
!       IF ( DO_SOLUTION_SAVING ) THEN
        DO M = 3, NMOMENTS
          QC = NMOMENTS - M + 1
          DO N = 1, NLAYERS
            Q = 0
            DO L = M, NMOMENTS
              IF(GREEKMAT_TOTAL_INPUT(L,N,1).EQ.ZERO)Q=Q+1
            ENDDO
            DO_LAYER_SCATTERING(M,N) = (Q.LT.QC)
          ENDDO
        ENDDO
!       ENDIF
      ENDIF

!  BVP telescoping (only if do_solution_saving is set)

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SOLUTION_SAVING ) THEN
        IF ( DO_BVP_TELESCOPING ) THEN
          DO M = 3, NMOMENTS
            Q = 0
            DO N = 1, NLAYERS
              IF (.NOT.DO_LAYER_SCATTERING(M,N))Q = Q + 1
            ENDDO
            IF ( Q.GT.1) BVP_REGULAR_FLAG(M) = .FALSE.
          ENDDO
        ENDIF
       ENDIF
      ENDIF

!    Set of Telescoped layers must be contiguous

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SOLUTION_SAVING ) THEN
        IF ( DO_BVP_TELESCOPING ) THEN
         DO M = 3, NMOMENTS
          NS = 0
          NAP = 0
          DO N = 1, NLAYERS
           IF ( DO_LAYER_SCATTERING(M,N) ) THEN
            NS = NS + 1
            NA = N
            IF ( NS.GT.1)  THEN
              IF ( NA.NE.NAP+1 ) THEN
                STATUS_SUB = VLIDORT_WARNING
                GO TO 4564
              ENDIF
            ENDIF
            NAP = NA
           ENDIF
          ENDDO
         ENDDO
        ENDIF

!  Collect warning and re-set default option

 4564   continue
        if ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
         MESSAGE = 'Telescoped layers not contiguous: turn off option'
         DO_BVP_TELESCOPING = .FALSE.
         DO M = 3, NMOMENTS
           BVP_REGULAR_FLAG(M) = .TRUE.
         ENDDO
        ENDIF

!  End clause

       ENDIF
      ENDIF

!  Continuation point for avoiding the default, if using Specialist # 2
!    -- Removed Version 2.8, not needed
! 4545 CONTINUE

!  Save calculation time for single scattering.
!   Remove isotropic-only option. 17 January 2006.
!   Bug. 22 November 2006. L must be < Ngreek_moments_input in DO WHILE
!   Bug 31 January 2007. LAYER_MAXMOMENTS should be = nmoments_input
!   Version 2.8, simplified

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_FOCORR .and..not. DO_SSCORR_USEFMAT ) THEN           ! 2.8
!       IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN         ! 2.7 and earlier
        IF ( DO_RAYLEIGH_ONLY ) THEN
          DO N = 1, NLAYERS
            LAYER_MAXMOMENTS(N) = 2
          ENDDO
        ELSE
         DO N = 1, NLAYERS
            L = 2
            LOOP = .TRUE.
            DO WHILE (LOOP.AND.L.LT.NGREEK_MOMENTS_INPUT)
              L = L + 1
              LOOP= (GREEKMAT_TOTAL_INPUT(L,N,1).NE.ZERO)
            ENDDO
!            LAYER_MAXMOMENTS(N) = L - 1
            LAYER_MAXMOMENTS(N) = L
          ENDDO
        ENDIF
       ENDIF
      ENDIF

!  Set the Eigensolver Mask. New section 21 December 2005.
!   1. Always set the Real Eigensolver for Rayleigh-only or Scalar-only
!   2. Otherwise, examine the occurrence of the "epsilon" Greek constant

      IF ( DO_RAYLEIGH_ONLY .OR. NSTOKES.EQ.1 ) THEN
       DO M = 0, NMOMENTS
        DO N = 1, NLAYERS
         DO_REAL_EIGENSOLVER(M,N) = .TRUE.
        ENDDO
       ENDDO
      ELSE
       DO M = 0, NMOMENTS
        DO N = 1, NLAYERS
         IF ( DO_LAYER_SCATTERING(M,N) ) THEN
           EC = 0
           DO L = M, NMOMENTS
!mick fix 1/22/2013 - added NSTOKES IF condition
            !IF (GREEKMAT_TOTAL_INPUT(L,N,12).NE.ZERO) EC = EC + 1
            IF (NSTOKES.EQ.4) THEN
              IF (GREEKMAT_TOTAL_INPUT(L,N,12).NE.ZERO ) EC = EC + 1
            ENDIF
           ENDDO
           DO_REAL_EIGENSOLVER(M,N) = (EC.EQ.ZERO)
!mick fix 2/17/11 - added ELSE
         ELSE
           DO_REAL_EIGENSOLVER(M,N) = .TRUE.
         ENDIF
        ENDDO
       ENDDO
      ENDIF

!  Debug

!      DO M = 0, NMOMENTS
!       write(*,'(I4,L2,5x,7L2)')M,BVP_REGULAR_FLAG(M),
!     &             (DO_LAYER_SCATTERING(M,N),N=1,7)
!        ENDDO
!      PAUSE

!  Continuation point.
!    -- 1/31/21. Version 2.8.3. Removed

! 5557 continue

!  User stream cosines and secants
!   Adjusted values removed, Version 2.8

      IF ( DO_USER_STREAMS ) THEN
        DO I = 1, N_USER_VZANGLES
          USER_STREAMS(I) = COS(DEG_TO_RAD*USER_VZANGLES(I))
!          USER_STREAMS(I) = COS(DEG_TO_RAD*USER_VZANGLES_ADJUST(I))
          USER_SECANTS(I) = ONE / USER_STREAMS(I)
          USER_SINES(I)   = SQRT(ONE-USER_STREAMS(I)*USER_STREAMS(I))
        ENDDO
      ENDIF

!  Number of tests to be applied for convergence
!  Number of streams for convergence
!    - Zenith tolerance no longer applies for the vector code
!    - Set CONV_START_STREAMS = 1. Do not examine for "Zenith_tolerance"

      CONV_START_STREAMS = 1
      N_CONV_STREAMS = N_OUT_STREAMS - CONV_START_STREAMS + 1

!  1/31/21. Version 2.8.3. Add Doublet geometry option

      IF ( DO_OBSERVATION_GEOMETRY ) THEN
        N_CONVTESTS = N_DIRECTIONS * N_USER_LEVELS
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        N_CONVTESTS = N_CONV_STREAMS * N_DIRECTIONS
        N_CONVTESTS = N_CONVTESTS * N_USER_LEVELS
      ELSE
        N_CONVTESTS = N_USER_RELAZMS * N_CONV_STREAMS * N_DIRECTIONS
        N_CONVTESTS = N_CONVTESTS * N_USER_LEVELS
      ENDIF

!  Sort out User vertical level outputs
!  ------------------------------------

!  Sort in ascending order. 
!    - Revised Call, 3/17/17

      IF ( N_USER_LEVELS .GT. 1 ) THEN
        CALL RSSORT(N_USER_LEVELS,USER_LEVELS)
      ENDIF

!  Mark all output levels not equal to layer boundary values

      NSTART = 0
      UT = 0
      DO UTA = 1, N_USER_LEVELS
        DT = USER_LEVELS(UTA)
        RT = DT - DBLE(INT(DT))
        N = INT(DT) + 1
        IF ( RT.GT.ZERO) THEN
          UT = UT + 1
          PARTLAYERS_OUTFLAG(UTA)  = .TRUE.
          PARTLAYERS_OUTINDEX(UTA) = UT
          PARTLAYERS_LAYERIDX(UT)  = N
          UTAU_LEVEL_MASK_UP(UTA) = N
          UTAU_LEVEL_MASK_DN(UTA) = N - 1
          PARTLAYERS_VALUES(UT)    = RT
        ELSE
          PARTLAYERS_OUTFLAG(UTA)  = .FALSE.
          PARTLAYERS_OUTINDEX(UTA) =   0
          UTAU_LEVEL_MASK_UP(UTA) = N - 1
          UTAU_LEVEL_MASK_DN(UTA) = N - 1
        ENDIF
      ENDDO
      N_PARTLAYERS = UT
      DO_PARTLAYERS = ( N_PARTLAYERS .NE. 0 )

!  Set masking and number of layer source terms
!  --------------------------------------------

!   .. for upwelling

!mick fix 6/25/2012 - initialize outside
      DO N = 1, NLAYERS
        STERM_LAYERMASK_UP(N) = .FALSE.
      ENDDO
      IF ( DO_UPWELLING ) THEN
        !DO N = 1, NLAYERS
        !  STERM_LAYERMASK_UP(N) = .FALSE.
        !ENDDO
        UTA = 1
        UT  = 1
        IF ( .NOT. PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_LAYERSOURCE_UP = UTAU_LEVEL_MASK_UP(UTA) + 1
          N_ALLLAYERS_UP   = N_LAYERSOURCE_UP
        ELSE
          N_LAYERSOURCE_UP = PARTLAYERS_LAYERIDX(UT) + 1
          N_ALLLAYERS_UP   = N_LAYERSOURCE_UP - 1
        ENDIF
        DO N = NLAYERS, N_ALLLAYERS_UP, -1
          STERM_LAYERMASK_UP(N) = .TRUE.
        ENDDO
      ENDIF

!   .. for downwelling

!mick fix 6/25/2012 - initialize outside
      DO N = 1, NLAYERS
        STERM_LAYERMASK_DN(N) = .FALSE.
      ENDDO
      IF ( DO_DNWELLING ) THEN
        !DO N = 1, NLAYERS
        !  STERM_LAYERMASK_DN(N) = .FALSE.
        !ENDDO
        UTA = N_USER_LEVELS
        UT  = N_PARTLAYERS
        IF ( .NOT. PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_LAYERSOURCE_DN = UTAU_LEVEL_MASK_DN(UTA)
          N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
        ELSE
          N_LAYERSOURCE_DN = PARTLAYERS_LAYERIDX(UT)
          N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
        ENDIF
        DO N = 1, N_ALLLAYERS_DN
          STERM_LAYERMASK_DN(N) = .TRUE.
        ENDDO
      ENDIF

!  Define TauGrid_Input
!  --------------------

      TAUGRID_INPUT(0) = ZERO
      DO N = 1, NLAYERS
        TAUGRID_INPUT(N) = TAUGRID_INPUT(N-1) + DELTAU_VERT_INPUT(N)
      END DO

!  PUT THE OFFSETS HERE

!     11 December 2020. PLACEHOLDER

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_DERIVE_INPUT

      END MODULE vlidort_inputs_m

