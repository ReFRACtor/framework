
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
! #            VBRDF_VSLEAVE_INPUT_CHECK                        #
! #            VBRDF_VSLEAVE_INPUT_CHECK_ERROR                  #
! #                                                             #
! ###############################################################

!  Upgrade Version 2.7. Notes by R. Spurr 5 May 2014
!  -------------------------------------------------

!    ** A new VBRDF_VSLEAVE_INPUT_CHECKER routine has been added, 
!       to ensure consistency between VSleave and VBRDF inputs, when the 
!       "New CM" ocean-treatment is in force

!  Upgrade Version 2.8
!  -------------------

!  Mick Mod 9/19/17.
!  VLIDORT standard VBRDF & VSLEAVE sup accessories separated into VBRDF, VSLEAVE,
!    and JOINT modules.
!  Renamed subroutine VBRDF_VSLEAVE_INPUT_CHECKER to VBRDF_VSLEAVE_INPUT_CHECK
!  Added subroutine VBRDF_VSLEAVE_INPUT_CHECK_ERROR:
!    ** Subroutine to perform error handling for VBRDF_VSLEAVE_INPUT_CHECK

      MODULE vlidort_joint_sup_accessories_m

      PRIVATE
      PUBLIC :: VBRDF_VSLEAVE_INPUT_CHECK, &
                VBRDF_VSLEAVE_INPUT_CHECK_ERROR

      CONTAINS

      SUBROUTINE VBRDF_VSLEAVE_INPUT_CHECK ( &
        VSLEAVE_Sup_In,             & ! Inputs
        VBRDF_Sup_In,               & ! Inputs
        VBRDF_VSLEAVECheck_Status )   ! Outputs

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_MESSAGES, ONE, &
                                 SMALLNUM, VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE VSLEAVE_Sup_Inputs_def_m
      USE VBRDF_Sup_Inputs_def_m
      USE VLIDORT_Outputs_def_m

      IMPLICIT NONE

      TYPE(VSLEAVE_Sup_Inputs), INTENT(IN)           :: VSLEAVE_Sup_In
      TYPE(VBRDF_Sup_Inputs), INTENT(IN)             :: VBRDF_Sup_In

      TYPE(VLIDORT_Exception_Handling), INTENT(OUT)  :: VBRDF_VSLEAVECheck_Status

!  ---------------
!  Local variables
!  ---------------

!  VSLEAVE supplement inputs
!  -------------------------

!  Surface-leaving control flags #1 (used for gatekeeping checks)

      LOGICAL :: SL_DO_ISOTROPIC
      LOGICAL :: SL_DO_FLUORESCENCE

!  Number of solar zenith angles

      INTEGER :: SL_NBEAMS

!  Surface-leaving control flags #2
!   Rob Fix, 3/7/17 added Rough Surface, Version 2.8

      LOGICAL :: SL_DO_GlintShadow
      LOGICAL :: SL_DO_FoamOption
      LOGICAL :: SL_DO_FacetIsotropy
      LOGICAL :: SL_DO_RoughSurface

!  Salinity

      DOUBLE PRECISION :: SL_SALINITY

!  Wind-speed and directions

      DOUBLE PRECISION :: SL_WINDSPEED
      DOUBLE PRECISION :: SL_WINDDIR ( MAX_SZANGLES )

!  Rob Fix 3/18/15. VSLEAVE Wavelength for the NewCM Glint Model [Micron]
!  Rob Fix 3/17/17. SLEAVE Wavelength for Fluorescence [nm]. (Version 2.8.3. No longer required)

      DOUBLE PRECISION :: SL_WAVELENGTH
!      DOUBLE PRECISION :: SL_FL_Wavelength

!  VBRDF supplement inputs
!  -----------------------

!  Number of solar zenith angles

      INTEGER :: BS_NBEAMS

!  Surface-leaving control flags

      LOGICAL :: BS_DO_GlintShadow
      LOGICAL :: BS_DO_FoamOption
      LOGICAL :: BS_DO_FacetIsotropy

!  Salinity

      DOUBLE PRECISION :: BS_SALINITY

!  Wind-speed and directions

      DOUBLE PRECISION :: BS_WINDSPEED
      DOUBLE PRECISION :: BS_WINDDIR ( MAX_SZANGLES )

!  Rob Fix 3/18/15. VBRDF Wavelength for the NewCM Glint Models [Microns]
!  update for NewGCM, 12/28/15

      LOGICAL          :: BS_DO_NewCMGLINT
      LOGICAL          :: BS_DO_NewGCMGLINT
      DOUBLE PRECISION :: BS_WAVELENGTH

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  VSLEAVE
!  -------

!  VSLEAVE Control inputs #1 (used for gatekeeping checks)

      SL_DO_ISOTROPIC     = VSLEAVE_Sup_In%SL_DO_ISOTROPIC
      SL_DO_FLUORESCENCE  = VSLEAVE_Sup_In%SL_DO_FLUORESCENCE

!  VSLEAVE Geometry inputs

      SL_NBEAMS           = VSLEAVE_Sup_In%SL_NBEAMS

!  VSLEAVE Control inputs #2

      SL_DO_GlintShadow   = VSLEAVE_Sup_In%SL_DO_GlintShadow
      SL_DO_FoamOption    = VSLEAVE_Sup_In%SL_DO_FoamOption
      SL_DO_FacetIsotropy = VSLEAVE_Sup_In%SL_DO_FacetIsotropy
      SL_DO_RoughSurface  = VSLEAVE_Sup_In%SL_DO_RoughSurface

!  VSLEAVE Other inputs

      SL_SALINITY         = VSLEAVE_Sup_In%SL_SALINITY
      SL_WINDSPEED        = VSLEAVE_Sup_In%SL_WINDSPEED
      SL_WINDDIR          = VSLEAVE_Sup_In%SL_WINDDIR

! Rob Fix 3/18/15
      SL_WAVELENGTH       = VSLEAVE_Sup_In%SL_WAVELENGTH

!  VBRDF
!  -----

! Rob Fix 3/18/15. Update for NewGCM, 12/28/15
      BS_DO_NewCMGLINT    = VBRDF_Sup_In%BS_DO_NewCMGLINT
      BS_DO_NewGCMGLINT   = VBRDF_Sup_In%BS_DO_NewGCMGLINT
      BS_WAVELENGTH       = VBRDF_Sup_In%BS_WAVELENGTH

!  VBRDF Geometry inputs

      BS_NBEAMS           = VBRDF_Sup_In%BS_NBEAMS

!  VBRDF Control inputs. New CM options

      BS_DO_GlintShadow   = VBRDF_Sup_In%BS_DO_GlintShadow
      BS_DO_FoamOption    = VBRDF_Sup_In%BS_DO_FoamOption
      BS_DO_FacetIsotropy = VBRDF_Sup_In%BS_DO_FacetIsotropy

!  VBRDF Other inputs. New CM options

      BS_SALINITY         = VBRDF_Sup_In%BS_SALINITY
      BS_WINDSPEED        = VBRDF_Sup_In%BS_WINDSPEED
      BS_WINDDIR          = VBRDF_Sup_In%BS_WINDDIR

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = VLIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of VSLEAVE/VBRDF compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Check on Beams
!  ==============

!  1/5/16. Very important

      IF ( BS_NBEAMS .ne. SL_NBEAMS ) then
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar beams in VBRDF NOT SAME as number in VSLEAVE'
        ACTIONS(NM)  = 'Make them the same!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS ;  GOTO 500
      ENDIF

!  NewCMGLINT/VSLEAVE gatekeeper checks
!  ===================================

!  (if any fails, don't look at further inputs unless directed, Rob Fix 1/5/16)

      IF ( BS_DO_NewCMGLINT .or. BS_DO_NewGCMGLINT ) then

!      IF ( SL_DO_ISOTROPIC ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'New Cox-Munk glint VBRDF is active and surface-leaving DO_ISOTROPIC is active'
!        ACTIONS(NM)  = 'Deactivate surface-leaving isotropy!'
!        STATUS_INPUTCHECK = VLIDORT_SERIOUS
!      ENDIF

!  Fluorescence not allowed, exit immediately

        IF ( SL_DO_FLUORESCENCE ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'New Cox-Munk glint VBRDF is active and surface-leaving DO_FLUORESCENCE is active'
          ACTIONS(NM)  = 'Deactivate surface-leaving fluorescence!'
          STATUS_INPUTCHECK = VLIDORT_SERIOUS ; GOTO 500
        ENDIF

!  Salinity, wavelength, Foam Option and Windspeed must agree, not matter what
!  Rob Fix 3/18/15, updated, 12/28/15  for NewGCM, 1/5/16 for better consistency
!  Rob Fix 3/7/17. Check wavelengths using adjacency (SMALLNUM = 1.0e-09). both Wavelengths in [Microns]

        IF ( SL_DO_FoamOption .neqv. BS_DO_FoamOption ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Surface-leaving foam option flag does not agree'
          ACTIONS(NM)  = 'Check flag compatibility!'
          STATUS_INPUTCHECK = VLIDORT_SERIOUS
        ENDIF

        IF ( SL_SALINITY .ne. BS_SALINITY) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Surface-leaving salinity does not agree'
          ACTIONS(NM)  = 'Check SL_SALINITY and BS_SALINITY input'
          STATUS_INPUTCHECK = VLIDORT_SERIOUS
        ENDIF

        if ( ABS ( (BS_WAVELENGTH/SL_WAVELENGTH) - ONE ) .gt. SMALLNUM ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'NewCM/NewGCM Glint: Wavelength does not agree with Water-leaving wavelength'
          ACTIONS(NM)  = 'Check input values of BS_WAVELENGTH and SL_WAVELENGTH'
          STATUS_INPUTCHECK = VLIDORT_SERIOUS
        endif

        IF ( SL_WINDSPEED .ne. BS_WINDSPEED) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Surface-leaving wind-speed does not agree'
          ACTIONS(NM)  = 'Check SL_WINDSPEED and BS_WINDSPEED input'
          STATUS_INPUTCHECK = VLIDORT_SERIOUS
        ENDIF

!  If any mismatches after these four, exit immediately.

        IF ( STATUS_INPUTCHECK == VLIDORT_SERIOUS ) GOTO 500

!  Remaining consistency checks, for the Water-Leaving Rough Surface (non-isotropic) option
!   Wind direction check only if there is NOT Facet Isotropy

        IF ( SL_DO_RoughSurface ) then

          IF ( SL_DO_FacetIsotropy .neqv. BS_DO_FacetIsotropy ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Surface-leaving facet isotropy flag does not agree'
            ACTIONS(NM)  = 'Check flag compatibility!'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          ENDIF

          IF ( SL_DO_GlintShadow .neqv. BS_DO_GlintShadow ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Surface-leaving glint shadow flag does not agree'
            ACTIONS(NM)  = 'Check flag compatibility!'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          ENDIF

          IF ( .not.SL_DO_FacetIsotropy .and. .not.BS_DO_FacetIsotropy ) THEN
            DO I = 1, BS_NBEAMS
              if ( SL_WINDDIR(I) .ne. BS_WINDDIR(I) ) THEN
                write(C2,'(I2)')I ; NM = NM + 1
                MESSAGES(NM) = 'Wind direction angle does not agree, # '//C2
                ACTIONS(NM)  = 'Check SL_WINDDIR and BS_WINDDIR input'
                STATUS_INPUTCHECK = VLIDORT_SERIOUS
              endif
            ENDDO
          ENDIF

        ENDIF

!  If any mismatches after these three, exit immediately.

        IF ( STATUS_INPUTCHECK == VLIDORT_SERIOUS ) GOTO 500

!  End of New CM and surface-leaving consistency check

      ENDIF

500   CONTINUE

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      VBRDF_VSLEAVECheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      VBRDF_VSLEAVECheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      VBRDF_VSLEAVECheck_Status%TS_CHECKMESSAGES     = MESSAGES
      VBRDF_VSLEAVECheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VBRDF_VSLEAVE_INPUT_CHECK

!

      SUBROUTINE VBRDF_VSLEAVE_INPUT_CHECK_ERROR ( ERRORFILE, VBRDF_VSLEAVECheck_Status )

!  Module, dimensions and numbers

      USE VBRDF_Sup_aux_m, ONLY : VBRDF_ERRUNIT  !(note: could have used VSLEAVE_ERRUNIT instead)
      USE VLIDORT_Outputs_def_m, ONLY : VLIDORT_Exception_Handling

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  VBRDF/VSLEAVE input check status

      TYPE(VLIDORT_Exception_Handling), intent(in) :: VBRDF_VSLEAVECheck_Status

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = VBRDF_ERRUNIT

!  Write VBRDF/VSLEAVE input compatibility errors to VBRDF error file
!  (Note: could write them to either the VBRDF or VSLEAVE error file here) 

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*)' FATAL: VBRDFSup and VSLEAVESup inputs are incompatible'
      WRITE(W,*)'  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)')'    ** Number of messages = ',VBRDF_VSLEAVECheck_Status%TS_NCHECKMESSAGES
      DO N = 1, VBRDF_VSLEAVECheck_Status%TS_NCHECKMESSAGES
        WRITE(W,'(A,I3,A,A)')'Message # ',N,': ',&
          ADJUSTL(TRIM(VBRDF_VSLEAVECheck_Status%TS_CHECKMESSAGES(N)))
        WRITE(W,'(A,I3,A,A)')'Action  # ',N,': ',&
          ADJUSTL(TRIM(VBRDF_VSLEAVECheck_Status%TS_ACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Checking fail: Look at file ' // TRIM(ERRORFILE)
      STOP ' Program stoppage from call to VBRDF_VSLEAVE_INPUT_CHECK_ERROR'

      END SUBROUTINE VBRDF_VSLEAVE_INPUT_CHECK_ERROR

!  Finish Module

      END MODULE vlidort_joint_sup_accessories_m
