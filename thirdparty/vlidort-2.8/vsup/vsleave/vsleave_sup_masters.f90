
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
! #            VSLEAVE_INPUTMASTER                              #
! #            VSLEAVE_MAINMASTER (master)                      #
! #                                                             #
! ###############################################################

!  1/31/21, Version 2.8.3.
!  Water-Leaving implementation has been re-written, several new things.

!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output
!    ==> Inputs : azimuth information (azms/nazms), number of do_Fourier_output-term azimuth quadratures
!    ==> Inputs : Doublet geometry option and angles
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_fOQ_BSF, calculates the terms needed for Fourier output  
!    ==> removed the Ta tayleigh stuff

      MODULE vsleave_sup_masters_m

      PRIVATE
      PUBLIC :: VSLEAVE_INPUTMASTER,&
                VSLEAVE_MAINMASTER

      CONTAINS

      SUBROUTINE VSLEAVE_INPUTMASTER ( &
        FILNAM, VSLEAVE_Sup_In, &
        VSLEAVE_Sup_InputStatus )

!  Input routine for VSLEAVE program

!  Version 2.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 December 2012. 
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  1/31/21, Version 2.8.3.
!  Water-Leaving implementation has been re-written, several new things.

!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output)
!    ==> Inputs : azimuth information (azms/nazms), number of do_Fourier_output-term azimuth quadratures
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_fOQ_BSF, calculates the terms needed for Fourier output  
!    ==> removed the Ta tayleigh stuff
!    ==> Explicit treatment of Doublet-geometry output

      USE VLIDORT_pars_m, Only : MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS, &
                                 MAXSTREAMS, MAXBEAMS, MAX_MESSAGES, FPK, ZERO, ONE, &
                                 VLIDORT_SUCCESS, VLIDORT_WARNING, VLIDORT_SERIOUS, VLIDORT_INUNIT

      USE VSLEAVE_FINDPAR_m

      USE vsleave_sup_inputs_def_m
      USE vsleave_sup_outputs_def_m

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN) :: FILNAM

      TYPE(VSLEAVE_Sup_inputs), INTENT(OUT) :: VSLEAVE_Sup_In

      TYPE(VSLEAVE_Input_Exception_Handling), INTENT(OUT) :: VSLEAVE_Sup_InputStatus

!  Local variables
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Rough Surface flag for Water-leaving. New 10/05/15 Rob Fix

      LOGICAL :: DO_ROUGHSURFACE

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Fluorescence flag

      LOGICAL :: DO_FLUORESCENCE

!  Solar sources flag

      LOGICAL :: DO_SOLAR_SOURCES

!  Path to VSLEAVE_DATA
!mick fix 9/19/2017 - added VSLEAVE_DATAPATH

      CHARACTER (LEN=200) :: VSLEAVE_DATAPATH

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Geometry and control
!  --------------------

!  Geometry flags 
!  1/31/21, Version 2.8.3. Doublet geometry option added

      LOGICAL   :: DO_USER_OBSGEOMS
      LOGICAL   :: DO_DOUBLET_GEOMETRY

!  Number of Stokes components

      INTEGER :: NSTOKES

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  Number of solar beams to be processed

      INTEGER   :: NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk) :: BEAM_SZAS (MAXBEAMS)

!  User-defined relative azimuths

      INTEGER   :: N_USER_RELAZMS
      REAL(fpk) :: USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input 

      INTEGER   :: N_USER_STREAMS
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  Observational geometry inputs

      INTEGER   :: N_USER_OBSGEOMS
      REAL(fpk) :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Local Doublet Geometry control and angles
!    -- 1/31/21. Version 2.8.3. Newly added

      INTEGER    :: N_USER_DOUBLETS
      REAL(fpk)  :: USER_DOUBLETS (MAX_USER_STREAMS,2)

!  Water-leaving basic variables
!  -----------------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  1/31/21. Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    ==> First introduced, 21 November 2020.
!    ==> Add Azimuth-dependent local flag DO_AZIMUTHDEP

      LOGICAL   :: DO_AZIMUTHDEP

!  1/31/21, Version 2.8.3. Fourier dependence in the water-leaving radiance
!    ==> Add Fourier output local flag DO_FOURIER_OUTPUT

      LOGICAL   :: DO_FOURIER_OUTPUT

!  Rough surface variables only (Now under separate control)
!  ---------------------------------------------------------

!  Changed for Version 2.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Removed, Version 2.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 2.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude 

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: FL_Amplitude755

!  Flag for using Data Gaussian parameters

      LOGICAL   :: FL_DO_DataGaussian

!  Input (non-data) Gaussians

      REAL(fpk) :: FL_InputGAUSSIANS(3,2)

!  Exception handling
!  ------------------

!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  local variables
!  ===============

      CHARACTER (LEN=12), PARAMETER :: PREFIX = 'VSLEAVESUP -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            I, FILUNIT, NM

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=TRIM(FILNAM),ERR=300,STATUS='OLD')

!  Initialize inputs
!  =================

!mick mod 9/19/2017 - reorganized initializations to mirror input type structure
!mick fix 9/19/2017 - added VSLEAVE_DATAPATH

!  Control flags

      DO_SLEAVING         = .FALSE.
      DO_ISOTROPIC        = .FALSE.
      DO_EXACT            = .FALSE.  !@@ New line
      DO_EXACTONLY        = .FALSE.
      DO_FLUORESCENCE     = .FALSE.

      DO_USER_STREAMS     = .FALSE.
      DO_SOLAR_SOURCES    = .FALSE.  !@@ New line
      DO_USER_OBSGEOMS    = .FALSE.
      DO_DOUBLET_GEOMETRY = .FALSE.

!  Data path

      VSLEAVE_DATAPATH  = ' '

!  Integer control

      NSTOKES  = 0
      NSTREAMS = 0

!  Conventional geometry

      NBEAMS          = 0
      N_USER_STREAMS  = 0
      N_USER_RELAZMS  = 0
      BEAM_SZAS       = ZERO
      USER_ANGLES     = ZERO
      USER_RELAZMS    = ZERO

!  Observational Geometry

      N_USER_OBSGEOMS = 0
      USER_OBSGEOMS   = ZERO

!  1/31/21, Version 2.8.3. Zero Doublet geometry flag, number and geometries

      N_USER_DOUBLETS = 0
      USER_DOUBLETS   = ZERO

!  Water-leaving variables

      SALINITY   = ZERO
      CHLORCONC  = ZERO
      WAVELENGTH = ZERO

      WINDSPEED  = ZERO
      WINDDIR    = ZERO

!  1/31/21, Version 2.8.3. Azimuth and Fourier dependence in the water-leaving radiance
!    ==> First introduced, 21 November 2020. Initialize flag DO_AZIMUTHDEP
!    ==> First introduced, 21 February 2021. Initialize flag DO_FOURIER_OUTPUT

      DO_ROUGHSURFACE   = .FALSE.
      DO_AZIMUTHDEP     = .FALSE.
      DO_FOURIER_OUTPUT = .FALSE.

      DO_GlintShadow   = .FALSE.
      DO_FoamOption    = .FALSE.
      DO_FacetIsotropy = .FALSE.

!  Fluorescence variables

      FL_WAVELENGTH      = ZERO
      FL_LATITUDE        = ZERO
      FL_LONGITUDE       = ZERO
      FL_EPOCH           = 0
      FL_Amplitude755    = ZERO
      FL_DO_DataGaussian = .FALSE.
      FL_InputGAUSSIANS  = ZERO

!  Geometry and Input Control
!  ==========================

!  !@@ Solar sources is True, always

      DO_SOLAR_SOURCES = .TRUE.

!  user-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  number of Stokes components

      PAR_STR = 'Number of Stokes vector components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR (ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  1/31/21. Version 2.8.3. GEOMETRY SECTION REORGANIZED, SAME AS IN BRDF SUPPLEMENT

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Observational/Doublet Geometry Choice
!  =====================================

!  Only if you are doing solar source, and only if you want user angle output

      IF ( DO_SOLAR_SOURCES .and. DO_USER_STREAMS ) THEN

!  !@@ New flag, Observational Geometry

         PAR_STR = 'Do Observation Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_USER_OBSGEOMS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  1/31/21, Version 2.8.3. Add Doublet geometry option

         PAR_STR = 'Do Doublet Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_DOUBLET_GEOMETRY
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  1/31/21, Version 2.8.3. Safety check on doublet geometry and observational geometry inputs

         IF ( DO_USER_OBSGEOMS .and. DO_DOUBLET_GEOMETRY) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Not allowed to have both Observation and Doublet Geometry options'
            ACTIONS(NM)  = 'Re-set input to one or the other of these flags'
            STATUS       = VLIDORT_SERIOUS
            NMESSAGES    = NM
            GO TO 764
         ENDIF

!  End Obsgeom/Doublet control flags

      ENDIF

!  Observational Geometry control
!  ==============================

!  1/31/21. Version 2.8.3. Section re-written 
!mick mod 1/5/2021 - added defining of DO_USER_STREAMS in section below

      IF ( DO_USER_OBSGEOMS ) THEN

!  Number of Observation Geometry inputs

        PAR_STR = 'Number of Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_OBSGEOMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of ObsGeometry inputs > Maximum dimension'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_OBSGEOMS dimension'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF

!  Observational Geometry inputs

        PAR_STR = 'Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_OBSGEOMS
             READ (FILUNIT,*,ERR=998) USER_OBSGEOMS(I,1), USER_OBSGEOMS(I,2), USER_OBSGEOMS(I,3)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

         NBEAMS          = N_USER_OBSGEOMS
         N_USER_STREAMS  = N_USER_OBSGEOMS
         N_USER_RELAZMS  = N_USER_OBSGEOMS
         DO_USER_STREAMS = .TRUE.

!  Automatic setting of BEAM_SZAS, USER_VZANGLES, and USER_RELAZMS

         BEAM_SZAS    (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
         USER_ANGLES  (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         USER_RELAZMS (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)

!  Finish - go to control point

         GO TO 5667

!  End Observational geometry clause

      ENDIF

!  Solar beam input (Lattice or Doublet)
!  =====================================

!  Number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  TOA solar zenith angle inputs

      PAR_STR = 'Solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Doublet Geometry input
!  ======================

!  1/31/21. Version 2.8.3. Section completely new

      IF ( DO_DOUBLET_GEOMETRY ) THEN

!  Number of Doublet Geometries

        PAR_STR = 'Number of Doublet Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_DOUBLETS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_DOUBLETS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of doublet geometry inputs > Maximum dimension MAX_USER_STREAMS'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_STREAMS dimension'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF

!  Doublet Geometry values

        PAR_STR = 'Doublet Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_DOUBLETS
             READ (FILUNIT,*,ERR=998)USER_DOUBLETS(I,1), USER_DOUBLETS(I,2)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of N_USER_STREAMS, N_USER_RELAZMS

        N_USER_STREAMS  = N_USER_DOUBLETS
        N_USER_RELAZMS  = N_USER_DOUBLETS

!  Automatic setting of user angles

        USER_ANGLES  (1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,1)
        USER_RELAZMS (1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,2)

!  Finish - go to control point

        GO TO 5667

!  End Doublet geometry clause

      ENDIF

!  Lattice Geometry controls
!  =========================
 
!  Make explicit the control here

      IF ( DO_USER_STREAMS .and. .not.DO_USER_OBSGEOMS .and. .not.DO_DOUBLET_GEOMETRY ) then

!  Number of User defined viewing zenith angles

        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  User-defined viewing zenith angles (should be positive)

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of azimuth angles

        PAR_STR = 'Number of user-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          NM = NM + 1
          MESSAGES(NM) =  'Number of relative azimuth angles > maximum dimension'
          ACTIONS(NM)  =  'Re-set input value or increase MAX_USER_RELAZMS dimension'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF

!  User defined Azimuth angles

        PAR_STR = 'User-defined relative azimuth angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End lattice input geometry clause

      ENDIF

!  Continuation point for Skipping the Lattice or Doublet inputs

5667  continue

!  Surface stuff
!  =============

!  General SLEAVING input
!  ----------------------

!  Basic flag

      PAR_STR = 'Do surface-leaving Contributions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SLEAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Skip if not doing surface leaving
!   @@ Rob Fix 04 August 2014

      if ( .not.DO_SLEAVING ) GOTO 652

!  Isotropic flag

      PAR_STR = 'Do Isotropic surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_ISOTROPIC
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  !@@ Overall-Exact flag

      PAR_STR = 'Do Overall-Exact surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
         READ (FILUNIT,*,ERR=998)DO_EXACT
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Exact only flag. Only if above is set (!@@)

      IF ( DO_EXACT ) THEN
        PAR_STR = 'Do Exact-only (no Fourier-term contributions)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_EXACTONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Basic source

      PAR_STR = 'Do surface-leaving Fluorescence?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Continuation point

652   continue

!  Inputs for Water-leaving (Non-Fluorescence case)
!  ------------------------------------------------

      IF ( DO_SLEAVING.and..not.DO_FLUORESCENCE ) THEN

!  Rob Fix 10/05/15, Water-leaving Rough-Surface Control
!    Update 1/5/16. Only for Non-isotropic water-leaving

!    -- 1/31/21. Version 2.8.3. Condition relaxed, now possible to have Rough surface with Isotropic

        PAR_STR = 'Do rough-surface water-leaving?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_ROUGHSURFACE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  1/31/21. Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Read control flag DO_AZIMUTHDEP, if non-isotropic

        if ( .not. DO_ISOTROPIC ) then
           PAR_STR = 'Do Azimuth-dependent water-leaving output?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
                READ (FILUNIT,*,ERR=998) DO_AZIMUTHDEP
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        endif

!  1/31/21, Version 2.8.3. Fourier dependence in the water-leaving radiance
!    ==> Read control flag DO_FOURIER_OUTPUT, if non-isotropic

        if ( .not. DO_ISOTROPIC ) then
           PAR_STR = 'Do Fourier dependence in water-leaving output?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
                READ (FILUNIT,*,ERR=998) DO_FOURIER_OUTPUT
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        endif

!  Basic variables: salinity, chlorophyll concentration, wavelength
!   The flat-surface WL terms can now be non-isotropic ! (10/05/15)

        PAR_STR = 'Ocean water salinity [ppt]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) SALINITY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Chlorophyll concentration in [mg/M]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) CHLORCONC
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Wavelength in [Microns]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Rob Fix, 1/5/16. Wind-speed and Whitecap (foam) option still good for all Water-leaving
!    (no longer part of the Rough surface non-isotropic case)

        PAR_STR = 'Windspeed in [m/s]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WINDSPEED
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Do whitecap (foam) calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) Do_FoamOption
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  New for Version 2.7. Glint, Control Flags for Facet Isotropy and Shadowing
!  Rob Fix 10/05/15, Water-leaving Rough-Surface Non-Isotropic Control, Upgraded 1/5/16
!     GlintShadow, FacetIsotropy (flags) and wind direction (Latter needs checking)
!   -- 1/31/21. Version 2.8.3.  Rough Surface also possible with Isotropic, Remove this line....
!          .... if ( do_roughsurface .and. .not. do_isotropic ) then

        if ( do_roughsurface ) then

          PAR_STR = 'Do glint calculation with facet isotropy?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) Do_FacetIsotropy
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do glint calculation with shadowing?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) DO_GlintShadow
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check Added 9/27/14. Multiple SZAS not allowed with Facet ANISOTROPY
!      LOGIC changed 12/18/15 to read wind directions.....

          IF ( .not. Do_FacetIsotropy ) then
            if ( NBEAMS.gt.1 ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Facet Anisotropy (wind direction) not Allowed for NBEAMS > 1'
              ACTIONS(NM)  = 'Either re-set NBEAMS = 1 or re-set Do_FacetIsotropy = .true.'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            else
              PAR_STR = 'Wind directions (degrees) relative to sun positions'
              IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
                DO I = 1, NBEAMS
                  READ (FILUNIT,*,ERR=998) WINDDIR(I)
                ENDDO
              ENDIF
              CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
            endif
          ENDIF

!  End non-isotropic clause, rough surface clause

        endif
      
!  Inputs for Fluorescence Case
!  ----------------------------

      ELSE IF ( DO_SLEAVING.and.DO_FLUORESCENCE ) THEN

!  Temporary Check. MUST BE ISOTROPIC FOR NOW

        IF ( .not. DO_ISOTROPIC ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'DO_ISOTROPIC was set to .FALSE. in fluorescence case'
          ACTIONS(NM)  = 'Tempo! Set DO_ISOTROPIC to .TRUE. if doing fluorescence'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Use of Data Gaussians (New, 8 August 2012)
!    IF NOT SET, YOU MUST USE YOUR OWN PARAMETERS

        PAR_STR = 'Do Data Gaussians in Fluorescence?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_DO_DataGaussian
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Amplitude for FS755 (Nominally, this is one)

        PAR_STR = 'Amplitude for Fluorescence model at 755 nm'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_Amplitude755
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Lat/Long, day-of-year, wavelength

        PAR_STR = 'Latitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LATITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Longitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LONGITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Epoch for Fluorescence model'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_EPOCH(1:6)
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Wavelength for Fluorescence model in [nm]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy general control variables
!mick fix 9/19/2017 - added VSLEAVE_DATAPATH

      VSLEAVE_Sup_In%SL_DO_SLEAVING      = DO_SLEAVING
      VSLEAVE_Sup_In%SL_DO_ISOTROPIC     = DO_ISOTROPIC
      VSLEAVE_Sup_In%SL_DO_ROUGHSURFACE  = DO_ROUGHSURFACE
      VSLEAVE_Sup_In%SL_DO_EXACT         = DO_EXACT         !@@
      VSLEAVE_Sup_In%SL_DO_EXACTONLY     = DO_EXACTONLY
      VSLEAVE_Sup_In%SL_DO_FLUORESCENCE  = DO_FLUORESCENCE
      VSLEAVE_Sup_In%SL_DO_SOLAR_SOURCES = DO_SOLAR_SOURCES   !@@
      VSLEAVE_Sup_In%SL_DO_USER_STREAMS  = DO_USER_STREAMS

      VSLEAVE_Sup_In%SL_VSLEAVE_DATAPATH = VSLEAVE_DATAPATH

!  Copy General Geometry results

      VSLEAVE_Sup_In%SL_NSTOKES           = NSTOKES
      VSLEAVE_Sup_In%SL_NSTREAMS          = NSTREAMS

      VSLEAVE_Sup_In%SL_NBEAMS            = NBEAMS
      VSLEAVE_Sup_In%SL_BEAM_SZAS         = BEAM_SZAS
      VSLEAVE_Sup_In%SL_N_USER_RELAZMS    = N_USER_RELAZMS
      VSLEAVE_Sup_In%SL_USER_RELAZMS      = USER_RELAZMS
      VSLEAVE_Sup_In%SL_N_USER_STREAMS    = N_USER_STREAMS
      VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT = USER_ANGLES

!  1/31/21, Version 2.8.3. Copy Observational geometry information

      VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS  = DO_USER_OBSGEOMS
      VSLEAVE_Sup_In%SL_N_USER_OBSGEOMS   = N_USER_OBSGEOMS
      VSLEAVE_Sup_In%SL_USER_OBSGEOMS     = USER_OBSGEOMS

!  1/31/21, Version 2.8.3. Copy Doublet geometry information

      VSLEAVE_Sup_In%SL_DO_DOUBLET_GEOMETRY = DO_DOUBLET_GEOMETRY
      VSLEAVE_Sup_In%SL_N_USER_DOUBLETS     = N_USER_DOUBLETS
      VSLEAVE_Sup_In%SL_USER_DOUBLETS       = USER_DOUBLETS

!  Copy Water-leaving inputs
!  -------------------------

!  Original, basic set

      VSLEAVE_Sup_In%SL_SALINITY         = SALINITY
      VSLEAVE_Sup_In%SL_CHLORCONC        = CHLORCONC
      VSLEAVE_Sup_In%SL_WAVELENGTH       = WAVELENGTH

!  Version 2.7 changes for Rough surface input

!      VSLEAVE_Sup_In%SL_NSTREAMS_AZQUAD  = NSTREAMS_AZQUAD
      VSLEAVE_Sup_In%SL_WINDSPEED        = WINDSPEED
      VSLEAVE_Sup_In%SL_WINDDIR          = WINDDIR

      VSLEAVE_Sup_In%SL_DO_GlintShadow   = DO_GlintShadow
      VSLEAVE_Sup_In%SL_DO_FoamOption    = DO_FoamOption
      VSLEAVE_Sup_In%SL_DO_FacetIsotropy = DO_FacetIsotropy

!  1/31/21. Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Set control flag DO_AZIMUTHDEP

      VSLEAVE_Sup_In%SL_AZIMUTHDEP = DO_AZIMUTHDEP

!  1/31/21, Version 2.8.3. Fourier dependence in the water-leaving radiance
!    ==> Set control flag DO_FOURIER_OUTPUT

      VSLEAVE_Sup_In%SL_DO_FOURIER_OUTPUT = DO_FOURIER_OUTPUT

!  Copy Fluorescence inputs
!  ------------------------

      VSLEAVE_Sup_In%SL_FL_LATITUDE        = FL_LATITUDE
      VSLEAVE_Sup_In%SL_FL_LONGITUDE       = FL_LONGITUDE
      VSLEAVE_Sup_In%SL_FL_WAVELENGTH      = FL_WAVELENGTH
      VSLEAVE_Sup_In%SL_FL_EPOCH           = FL_EPOCH
      VSLEAVE_Sup_In%SL_FL_Amplitude755    = FL_Amplitude755
      VSLEAVE_Sup_In%SL_FL_DO_DataGaussian = FL_DO_DataGaussian
      VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS  = FL_InputGAUSSIANS

!  Exception handling

      VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      VSLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  line read error - abort immediately

998   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)
      GO TO 764

!  Final error copying

764   CONTINUE

      VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      VSLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VSLEAVE_INPUTMASTER

!

      SUBROUTINE VSLEAVE_MAINMASTER ( &
        VSLEAVE_Sup_In,         & ! Inputs
        VSLEAVE_Sup_Out,        & ! Outputs
        VSLEAVE_Sup_OutputStatus) ! Outputs

!  Prepares the Surface Leaving necessary for VLIDORT.

!  Version 2.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 December 2012. 
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  1/31/21, Version 2.8.3.
!  Water-Leaving implementation has been re-written, several new things.

!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output
!    ==> Inputs : azimuth information (azms/nazms), number of do_Fourier_output-term azimuth quadratures
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_fOQ_BSF, calculates the terms needed for Fourier output  
!    ==> removed the Ta tayleigh stuff
!    ==> Explicit treatment of doublet goemetry

      USE VLIDORT_pars_m, Only : MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS, &
                                 MAXSTREAMS, MAXBEAMS, MAX_MESSAGES, MAXMOMENTS, FPK, ZERO, ONE, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE vsleave_sup_inputs_def_m
      USE vsleave_sup_outputs_def_m

      USE vsleave_sup_aux_m        , only : GETQUAD2
      USE vsleave_sup_routines_2_m , only : WaterLeaving_2EE,          &
                                            get_fluorescence_755,      &
                                            solar_spec_irradiance

      IMPLICIT NONE

!  Input structure
!  ---------------

      TYPE(VSLEAVE_Sup_Inputs), INTENT(IN)   :: VSLEAVE_Sup_In

!  Output structure
!  ----------------

      TYPE(VSLEAVE_Sup_Outputs), INTENT(OUT) :: VSLEAVE_Sup_Out
!mick mod 9/19/2017 - added output exception handling for Version 2.8   
      TYPE(VSLEAVE_Output_Exception_Handling), INTENT(OUT) :: VSLEAVE_Sup_OutputStatus

!  VLIDORT local variables
!  +++++++++++++++++++++++

!  Input arguments
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Rough Surface flag for Water-leaving. New 10/05/15 Rob Fix

      LOGICAL :: DO_ROUGHSURFACE

!  Fluorescence flag

      LOGICAL :: DO_FLUORESCENCE

!  Solar sources flag

      LOGICAL :: DO_SOLAR_SOURCES

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Geometry and control
!  --------------------

!  Number of Stokes components

      INTEGER :: NSTOKES

!  Number of discrete ordinate streams, quadrature
!   Version 2.7, added the qudrature arrays

      INTEGER   :: NSTREAMS
      REAL(fpk) :: STREAMS (MAXSTREAMS)
      REAL(fpk) :: WEIGHTS (MAXSTREAMS)

!  Local angle control (general)

      INTEGER   :: NBEAMS
      INTEGER   :: N_USER_STREAMS
      INTEGER   :: N_USER_RELAZMS
      REAL(fpk) :: BEAM_SZAS   (MAXBEAMS)
      REAL(fpk) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  Local Observational Geometry control and angles

      LOGICAL    :: DO_USER_OBSGEOMS
      INTEGER    :: N_USER_OBSGEOMS
      REAL(fpk)  :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Local Doublet Geometry control and angles
!    -- 1/31/21. Version 2.8.3. Newly added

      LOGICAL    :: DO_DOUBLET_GEOMETRY
!      INTEGER    :: N_USER_DOUBLETS
!      REAL(fpk)  :: USER_DOUBLETS (MAX_USER_STREAMS,2)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  1/31/21. Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Add control flag do_Azimuth_Output

      LOGICAL   :: do_Azimuth_Output

!  1/31/21, Version 2.8.3. Fourier dependence in the water-leaving radiance
!    ==> Add control flag do_Fourier_Output

      LOGICAL   :: DO_FOURIER_OUTPUT

!  Changed for Version 2.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Removed, Version 2.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 2.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: FL_Amplitude755

!  Flag for using Data Gaussians

      LOGICAL   :: FL_DO_DataGaussian

!  Output variables in structures (Commented out)
!  ------------------------------

!  Exact Surface-Leaving term

!   REAL(fpk), dimension ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS ) :: SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

!   REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )       :: SLTERM_F_0
!   REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS ) :: USER_SLTERM_F_0

!  Other local variables
!  =====================

!  General
!  -------

!  Exception handling
!     Message Length should be at least 120 Characters
!mick mod 9/19/2017 - exception handling upgraded for Version 2.8

      LOGICAL :: FAIL
      CHARACTER (LEN=200) :: MESSAGE

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )

!  Set VSLEAVE DATA PATH

      CHARACTER (LEN=200) :: VSLEAVE_DATAPATH

!  Water-leaving model
!  -------------------

!  1/31/21. Version 2.8.3. No more Ta stuff

!  Files for Data

      CHARACTER (LEN=200) :: FoQFile
!      CHARACTER (LEN=200) :: TaRayFile

!  (Version 2.6 code). Water-leaving model
!      REAL :: WAV,CHL,RW,SAL,A,REFR,REFI,N12,RWB,TDS,TDV

!  Approximate Tranmsittance Flag. Moved here, 05 October 15 (RJDS)
!    -  If set, WaterLeaving will return the Gordon/Wang simple approximation
!    -  If not set, Get rayleigh-atmosphere transmittance from data base
!      logical, parameter :: do_Approximate_Ta = .true.
!      LOGICAL, PARAMETER :: do_Approximate_Ta = .false.

!  Isotropic value. Fast calculation

      REAL(fpk)    :: WLeaving_ISO ( MAXBEAMS )

!  Input solar, output stream angles. 1/31/21. Version 2.8.3. ALL FOURIERS!

      REAL(fpk)    :: WLeaving_SDF ( 0:MAXMOMENTS, MAXBEAMS, MAXSTREAMS )

!  input solar, output view angles. 1/31/21. Version 2.8.3. ALL FOURIERS!

      REAL(fpk)    :: WLeaving_SVF ( 0:MAXMOMENTS, MAXBEAMS, MAX_USER_STREAMS )

!  1/31/21. Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Add Azimuth-dependent term WLeaving_SVA, controlled by flag do_Azimuth_Output
!    ==> MS contributions are still azimuth-independent (for now)

      REAL(fpk)    :: WLeaving_SVA ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Atmospheric Transmittance
!      REAL(fpk)    :: TaSav ( MAXBEAMS, 4 )

!  1/31/21. Version 2.8.3. Local dimensioning for Fouriers
!    - Must be even numbers, one twice the other

     INTEGER, parameter :: Maxazmquads = 100
     INTEGER, parameter :: Maxaqhalf   = 50

!  Fluorescence model
!  ------------------

!  Files for Data

      CHARACTER*200 :: Fluofile, Solfile1, Solfile2

!  Fluorescence Isotropic Surface leaving term

      REAL(fpk) :: Fluor_ISOTROPIC

!  Fluorescence Gaussian parameters
!     Parameters of the fluorescence Gaussian spectral shape model.
!           Gaussian    A (Wm−2 μm−1 sr−1) Lambda(nm) Sigma(nm)
!              1           1.445           736.8        21.2
!              2           0.868           685.2        9.55

      REAL(FPK) :: FL_DataGAUSSIANS(3,2), FL_GAUSSIANS(3,2)
      data FL_DataGAUSSIANS(1,1) / 1.445d0 /
      data FL_DataGAUSSIANS(2,1) / 736.8d0 /
      data FL_DataGAUSSIANS(3,1) / 21.2d0  /
      data FL_DataGAUSSIANS(1,2) / 0.868d0 /
      data FL_DataGAUSSIANS(2,2) / 685.2d0 /
      data FL_DataGAUSSIANS(3,2) / 9.55d0  /

!  Solar spectral radiance model wavelength

      REAL(FPK) :: ssr_wvl

!  Help variables

      INTEGER   :: I, IB, K, UM, LUM, LUA, NMOMENTS, nazmquads, localm
      REAL(FPK) :: Fs755(MAXBEAMS), FL_SunSpec, FsSum
      REAL(FPK) :: ampli, lamda, sigma, arg, gauss, var, ff

!  Start Code
!  ==========

!  1/31/21. Version 2.8.3. Set the Obsgeom/Doublet indices here

      LUA = 1 ; LUM = 1

!  Initialize Exception handling
!  -----------------------------

!mick mod 9/19/2017 - initialize exception handling as part of upgrade for Version 2.8

      STATUS = VLIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Execution of VLIDORT VSLEAVE Sup Master'

!  Copy from input structure
!  -------------------------

!  Set Vsleave_DataPath. Input setting, 12/28/15

!      Vsleave_DataPath = 'VSLEAVE_DATA/'
      Vsleave_DataPath = Trim(VSLEAVE_Sup_In%SL_VSLEAVE_DATAPATH)

!  Copy Top-level general Control inputs

      DO_USER_STREAMS = VSLEAVE_Sup_In%SL_DO_USER_STREAMS
      DO_SLEAVING     = VSLEAVE_Sup_In%SL_DO_SLEAVING

      DO_EXACT        = VSLEAVE_Sup_In%SL_DO_EXACT          !@@
      DO_EXACTONLY    = VSLEAVE_Sup_In%SL_DO_EXACTONLY
      DO_FLUORESCENCE = VSLEAVE_Sup_In%SL_DO_FLUORESCENCE
      DO_ISOTROPIC    = VSLEAVE_Sup_In%SL_DO_ISOTROPIC
      DO_ROUGHSURFACE = VSLEAVE_Sup_In%SL_DO_ROUGHSURFACE   ! New, 10/05/15 Rob Fix
      DO_SOLAR_SOURCES = VSLEAVE_Sup_In%SL_DO_SOLAR_SOURCES

!  Set number of stokes elements and streams
!   Stream Quadrature new for Version 2.7. Revised Call, 3/17/17

      NSTOKES  = VSLEAVE_Sup_In%SL_NSTOKES
      NSTREAMS = VSLEAVE_Sup_In%SL_NSTREAMS
      CALL GETQUAD2 ( ZERO, ONE, NSTREAMS, STREAMS(1:nstreams), WEIGHTS(1:nstreams) )
      NMOMENTS = 2 * NSTREAMS - 1

!  Geometry proxy settings. 1/31/21. Version 2.8.3. Rewritten
!  -----------------------

!   !@@ Either set from User Observational Geometry
!          Or Copy from Usual lattice input

      DO_USER_OBSGEOMS    = VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS
      DO_DOUBLET_GEOMETRY = VSLEAVE_Sup_In%SL_DO_DOUBLET_GEOMETRY

!  1/31/21, Version 2.8.3. Add settings for the Doublet geometry option

      IF ( DO_USER_OBSGEOMS ) THEN
        N_USER_OBSGEOMS = VSLEAVE_Sup_In%SL_N_USER_OBSGEOMS
        USER_OBSGEOMS   = VSLEAVE_Sup_In%SL_USER_OBSGEOMS
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS          = N_USER_OBSGEOMS
          N_USER_STREAMS  = N_USER_OBSGEOMS
          N_USER_RELAZMS  = N_USER_OBSGEOMS
          BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
          USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
          USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = N_USER_OBSGEOMS
          USER_ANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        ENDIF
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS            = VSLEAVE_Sup_In%SL_NBEAMS
          BEAM_SZAS         = VSLEAVE_Sup_In%SL_BEAM_SZAS
          N_USER_STREAMS    = VSLEAVE_Sup_In%SL_N_USER_DOUBLETS
          N_USER_RELAZMS    = N_USER_STREAMS
          USER_ANGLES (1:N_USER_STREAMS) = VSLEAVE_Sup_In%SL_USER_DOUBLETS(1:N_USER_STREAMS,1)
          USER_RELAZMS(1:N_USER_RELAZMS) = VSLEAVE_Sup_In%SL_USER_DOUBLETS(1:N_USER_STREAMS,2)
        ELSE
!  NOT ALLOWED
!          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
!          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
!          N_USER_STREAMS = VBRDF_Sup_In%BS_N_USER_STREAMS
!          USER_ANGLES = VBRDF_Sup_In%BS_USER_ANGLES_INPUT
        ENDIF
      ELSE
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS         = VSLEAVE_Sup_In%SL_NBEAMS
          BEAM_SZAS      = VSLEAVE_Sup_In%SL_BEAM_SZAS
          N_USER_RELAZMS = VSLEAVE_Sup_In%SL_N_USER_RELAZMS
          USER_RELAZMS   = VSLEAVE_Sup_In%SL_USER_RELAZMS
          N_USER_STREAMS = VSLEAVE_Sup_In%SL_N_USER_STREAMS
          USER_ANGLES    = VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = VSLEAVE_Sup_In%SL_N_USER_STREAMS
          USER_ANGLES    = VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT
        ENDIF
      ENDIF

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      SALINITY        = VSLEAVE_Sup_In%SL_SALINITY
      CHLORCONC       = VSLEAVE_Sup_In%SL_CHLORCONC
      WAVELENGTH      = VSLEAVE_Sup_In%SL_WAVELENGTH

!  1/31/21, Version 2.8.3. Fourier dependence in the water-leaving radiance
!    ==> Copy control flag do_Fourier_Output

      DO_FOURIER_OUTPUT = VSLEAVE_Sup_In%SL_DO_FOURIER_OUTPUT

!  1/31/21, Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Copy control flag do_Azimuth_Output

      do_Azimuth_Output = VSLEAVE_Sup_In%SL_AZIMUTHDEP

!  Version 2.7 changes

!      NSTREAMS_AZQUAD = VSLEAVE_Sup_In%SL_NSTREAMS_AZQUAD
      WINDSPEED       = VSLEAVE_Sup_In%SL_WINDSPEED
      WINDDIR         = VSLEAVE_Sup_In%SL_WINDDIR

      DO_GlintShadow   = VSLEAVE_Sup_In%SL_DO_GlintShadow
      DO_FoamOption    = VSLEAVE_Sup_In%SL_DO_FoamOption
      DO_FacetIsotropy = VSLEAVE_Sup_In%SL_DO_FacetIsotropy

!  Copy Fluorescence inputs
!  ------------------------

!  Main variables

      FL_Wavelength   = VSLEAVE_Sup_In%SL_FL_Wavelength
      FL_Latitude     = VSLEAVE_Sup_In%SL_FL_Latitude
      FL_Longitude    = VSLEAVE_Sup_In%SL_FL_Longitude
      FL_Epoch        = VSLEAVE_Sup_In%SL_FL_Epoch
      FL_Amplitude755 = VSLEAVE_Sup_In%SL_FL_Amplitude755
      FL_DO_DataGaussian = VSLEAVE_Sup_In%SL_FL_DO_DataGaussian

!mick fix 8/31/2012 - added outer if block
      if (DO_FLUORESCENCE) then
        if ( FL_DO_DataGaussian ) then
           FL_GAUSSIANS(1:3,1) = FL_DataGAUSSIANS(1:3,1)
           FL_GAUSSIANS(1:3,2) = FL_DataGAUSSIANS(1:3,2)
        else
           FL_GAUSSIANS(1:3,1) = VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,1)
           FL_GAUSSIANS(1:3,2) = VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,2)
        endif
      endif

!  Main code
!  ---------

!  Zero the output

      VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC  = ZERO
      VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES = ZERO
      VSLEAVE_Sup_Out%SL_SLTERM_F_0        = ZERO
      VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0   = ZERO
      VSLEAVE_Sup_Out%SL_TRANS_ATMOS       = ZERO  !   New, 10/5/15

!  Return if not called
!   ROb Fix, 04 August 2014

      IF ( .not. DO_SLEAVING ) RETURN

!  Fluorescence
!  ============

      IF ( DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) Stop 'Non-isotropic not allowed yet if doing fluorescence'

!  F_755 data file. Rob Fix, use directed path. 10/5/15

!        Fluofile = 'vlidort_v_test/data/fluor_data_2009_fortran.dat'
        Fluofile = Trim(Vsleave_DataPath)//'/fluor_data_2009_fortran.dat'

!  Solar Files. Rob Fix, use directed path. 10/5/15

        Solfile1 = Trim(Vsleave_DataPath)//'/wehrli85.dat'
        Solfile2 = Trim(Vsleave_DataPath)//'/ref_solar_irradiance_whi-2008_ver2.dat'

!  Get solar spectral radiance, in (W m−2 μm−1 sr−1), to normalize data
!      Rob Fix, use directed path. 10/5/15

        !FL_SunSpec = 1.0d0  ! Temporary
        ssr_wvl = FL_Wavelength*1.0d-3 !convert from nm to um
        FL_SunSpec = solar_spec_irradiance( ssr_wvl, Solfile1, Solfile2 )

!  factor: After  some fiddling, this is 1.0 (July 30th, 2012)
!    4pi factor is required in DB correction ter,
!         FF = PI4
        FF = 1.0d0

!  For each Solar zenith angle

        DO IB = 1, NBEAMS

 !  Get the F_755 data from the subroutine

          CALL get_fluorescence_755 &
   ( FL_Latitude, FL_Longitude, FL_Epoch, BEAM_SZAS(IB), FluoFile, Fs755(IB) )

!  Apply Gaussians

          FsSum = zero
          do k = 1, 2
            ampli = FL_Gaussians(1,k)
            lamda = FL_Gaussians(2,k)
            sigma = FL_Gaussians(3,k)
            var = 0.5d0/sigma/sigma
            arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
            Gauss = zero
            if ( arg.lt.88.0d0 ) gauss = ampli * dexp ( - arg )
            FsSum = FsSum + Gauss
          enddo

!  Assign output Fluorescence (Apply Amplitude)
!  multiply by Fs755, and normalize to solar spectrum
!   FF is the fudge factor

          Fluor_ISOTROPIC = FF * FsSum * Fs755(IB) / FL_SunSpec
          VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,IB) = FL_Amplitude755 * Fluor_ISOTROPIC

!  End Beam loop

        ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( .not. DO_FLUORESCENCE ) THEN

!  Rob Fix, use directed path. 10/5/15. 1/31/21, Version 2.8.3. Skip Ta stuff

         FoQFile   = Trim(Vsleave_DataPath)//'/values0.dat'
!         TaRayFile = Trim(Vsleave_DataPath)//'/Rayleigh_TransAtmos_Table_270900_14Szas.all'

!  1/31/21, Version 2.8.3.  Local setting of quadrature

         nazmquads = Maxazmquads

!  Call the routine for Version 2.7. Updated, 01-05 October 2015.
!        01 October, new flag for Roughsurface
!        05 October, Data file names, TaSav output, Approximate_Ta flag

!  1/31/21, Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    ==> First introduced, 21 November 2020.
!    ==> Add Azimuth-dependent term WLeaving_SVA, controlled by flag do_Azimuth_Output
!    ==> MS contributions are still azimuth-independent (for now)

         Call WaterLeaving_2EE &
          ( Maxbeams, Max_User_Streams, Max_User_Relazms, Maxstreams,            & ! Dimensions
            Maxmoments, Maxazmquads, Maxaqhalf,                                  & ! Dimensions
            FoQFile, do_Isotropic, do_Azimuth_Output, do_Fourier_Output,         & ! Flags and file
            Do_RoughSurface, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,    & ! Glitter flags
            Wavelength, Salinity, ChlorConc, Windspeed, WindDir,                 & ! Physical inputs
            nbeams, n_user_streams, n_user_relazms, nstreams, nazmquads,         & ! Geometry
            beam_szas, user_angles, user_relazms, streams,                       & ! geometry
            WLeaving_ISO, WLeaving_SDF, WLeaving_SVF, WLeaving_SVA, fail, message )       

!  old code
!         Call WaterLeaving_2E &
!         ( Maxbeams, Max_User_Streams, Max_User_Relazms, Maxstreams,               &
!           FoQFile, TaRayFile, do_Approximate_Ta, do_Isotropic, do_Azimuth_Output, &
!           Do_RoughSurface, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,       &
!           Wavelength, Salinity, ChlorConc, Windspeed, WindDir,                    &
!           nbeams, n_user_streams, n_user_relazms, nstreams, beam_szas, user_angles, user_relazms, streams, &
!           WLeaving_ISO, WLeaving_SD, WLeaving_SV, WLeaving_SVA, TaSav, fail, message )

         if ( fail ) then
!mick mod 9/19/2017 - exception handling upgraded for Version 2.8
           !write(*,'(A)')Trim(message) ; Stop 'WaterLeaving_2 failed'
           STATUS = VLIDORT_SERIOUS
           NMESSAGES = NMESSAGES + 1
           MESSAGES(NMESSAGES) = 'Fatal - WaterLeaving_2EE failed in VSLEAVE_MAINMASTER'
           NMESSAGES = NMESSAGES + 1
           MESSAGES(NMESSAGES) = Trim(message)
           GO TO 899
         endif

!  Copy to Type structure arrays
!    If Isotropic, VLIDORT takes care of the rest

!  1/31/21. Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Add Azimuth-dependent term WLeaving_SVA, controlled by flag do_Azimuth_Output
!    ==> Add All Fourier terms, if flagged ( do_Fourier_output = .true., Non-Isotropic, use nmoments)
!    ==> Introduce doublet and observational geometry output settings

         VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,1:NBEAMS) = WLeaving_ISO(1:NBEAMS)
         localm = nmoments ; if ( do_Isotropic.or..not.do_Fourier_output ) localm = 0
         do ib = 1, nbeams
            if ( do_user_obsgeoms ) then
              VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1,lum,lua,ib)     = WLeaving_SVA(ib,ib,ib)
              VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0:localm,1,ib,ib) = WLeaving_SVF(0:localm,ib,ib)
            else if ( do_doublet_geometry ) then
              do um = 1, n_user_streams
                VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1,um,lua,ib)      = WLeaving_SVA(ib,um,um)
                VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0:localm,1,um,ib) = WLeaving_SVF(0:localm,ib,um)
              enddo
            else
              do um = 1, n_user_streams
                VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1,um,1:n_user_relazms,ib) = WLeaving_SVA(ib,um,1:n_user_relazms)
                VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0:localm,1,um,ib)         = WLeaving_SVF(0:localm,ib,um)
              enddo
            endif
            do i = 1, nstreams
              VSLEAVE_Sup_Out%SL_SLTERM_F_0(0:nmoments,1,i,ib) = WLeaving_SDF(0:nmoments,ib,i)
            enddo
         enddo

!  Atmospheric Transmittance (Diagnostic output). 1/31/21. Version 2.8.3.Skip Ta stuff, set to 1
!mick fix 9/19/2017 - reduced SL_TRANS_ATMOS to one dimension

         VSLEAVE_Sup_Out%SL_TRANS_ATMOS(1:NBEAMS) = one
!         VSLEAVE_Sup_Out%SL_TRANS_ATMOS(1:NBEAMS) = TaSav(1:NBEAMS,1)

!  Here is the compressed Version 2.6 code ---------------------------------------
!    INDWAT/MORCASIWAT calls . use single precision routine
!        SAL = REAL(SALINITY) ; WAV = REAL(WAVELENGTH)
!        CALL INDWAT(WAV,SAL,refr,refi)
!        CHL = REAL(CHLORCONC) ; WAV = REAL(WAVELENGTH)
!        CALL MORCASIWAT(WAV,CHL,RW,.false.)
!     Perfect Transmittance. Add change in solid angle through surface
!     that accounts for 1/(n12*n12) decrease in directional reflectivity
!        a   = 0.485 ; tds = 1.0 ; tdv = 1.0
!        n12 = refr*refr + refi*refi  ; n12 = sqrt(n12)
!        Rwb=(1.0/(n12*n12))*tds*tdv*Rw/(1-a*Rw)
!        VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,1:NBEAMS)= DBLE(Rwb)

!  end water leaving

      endif

!mick mod 9/19/2017 - added these two exception handling sections for Version 2.8

!  Continuation point for error finish

899   continue

!  Write Exception handling to output structure

      VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT   = STATUS
      VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES = NMESSAGES
      VSLEAVE_Sup_OutputStatus%SL_OUTPUTMESSAGES  = MESSAGES

!  Finish

      RETURN
      END SUBROUTINE VSLEAVE_MAINMASTER

      END MODULE vsleave_sup_masters_m
