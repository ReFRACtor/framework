
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
! #            VLIDORT_VBRDF_INPUT_CHECK                        #
! #            VLIDORT_VBRDF_INPUT_CHECK_ERROR                  #
! #            SET_VLIDORT_VBRDF_INPUTS                         #
! #                                                             #
! ###############################################################

!  Upgrade Version 2.7. Notes by R. Spurr 5 May 2014
!  -------------------------------------------------

!    ** Although the VBRDF supplement has had major upgrades
!       for Version 2.7, the "Input_Check" subroutine listed here
!       has not changed - the geometrical inputs for the supplement
!       are as they have always been!

!  Upgrade Version 2.8
!  -------------------

!  Rob Fix 3/7/17. Wavelength checks introduced.

!  Mick Mod 9/19/17.
!  VLIDORT standard VBRDF & VSLEAVE sup accessories separated into VBRDF, VSLEAVE,
!    and JOINT modules.
!  Renamed subroutine VLIDORT_VBRDF_INPUT_CHECKER to VLIDORT_VBRDF_INPUT_CHECK.
!  Added subroutine VLIDORT_VBRDF_INPUT_CHECK_ERROR:
!    ** Subroutine to perform error handling for VLIDORT_VBRDF_INPUT_CHECK
!  Added subroutine SET_VLIDORT_VBRDF_INPUTS:
!    ** Subroutine to define the main VLIDORT VBRDF inputs by defining
!       them using the corresponding VLIDORT VBRDF supplement outputs 

      MODULE vlidort_vbrdf_sup_accessories_m

      PRIVATE
      PUBLIC :: VLIDORT_VBRDF_INPUT_CHECK, &
                VLIDORT_VBRDF_INPUT_CHECK_ERROR, &
                SET_VLIDORT_VBRDF_INPUTS

      CONTAINS

      SUBROUTINE VLIDORT_VBRDF_INPUT_CHECK ( &
        VBRDF_Sup_In,             & ! Inputs
        VLIDORT_FixIn,            & ! Inputs
        VLIDORT_ModIn,            & ! Inputs
        VLIDORT_VBRDFCheck_Status ) ! Outputs

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_MESSAGES, &
                                 ONE, SMALLNUM, VLIDORT_SUCCESS, VLIDORT_SERIOUS


!  Other

      USE VBRDF_Sup_Inputs_def_m
      USE VLIDORT_Inputs_def_m
      USE VLIDORT_Outputs_def_m

      IMPLICIT NONE

      TYPE(VBRDF_Sup_Inputs), INTENT(IN)             :: VBRDF_Sup_In

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)        :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (IN)     :: VLIDORT_ModIn

      TYPE(VLIDORT_Exception_Handling), INTENT(OUT)  :: VLIDORT_VBRDFCheck_Status

!  ---------------
!  Local variables
!  ---------------

!  VBRDF supplement inputs
!  -----------------------

!  User stream, BRDF surface and surface emission flags

      LOGICAL ::             BS_DO_USER_STREAMS
      LOGICAL ::             BS_DO_BRDF_SURFACE
      LOGICAL ::             BS_DO_SURFACE_EMISSION

!  Number of Stokes vector components

      INTEGER ::             BS_NSTOKES

!  Number of discrete ordinate streams

      INTEGER ::             BS_NSTREAMS

!  BOA solar zenith angles

      INTEGER   ::           BS_NBEAMS
      DOUBLE PRECISION ::    BS_BEAM_SZAS (MAX_SZANGLES)

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER   ::           BS_N_USER_RELAZMS
      DOUBLE PRECISION ::    BS_USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input

      INTEGER   ::           BS_N_USER_STREAMS
      DOUBLE PRECISION ::    BS_USER_ANGLES_INPUT (MAX_USER_VZANGLES)

!  Rob Fix 3/18/15. VBRDF Wavelength for the NewCM Glint Models [Microns]
!   Updated for New GCM, 12/28/15

      LOGICAL   ::           BS_DO_NewCMGLINT
      LOGICAL          ::    BS_DO_NewGCMGLINT
      DOUBLE PRECISION ::    BS_WAVELENGTH

!  VLIDORT Main inputs
!  -------------------

!  VLIDORT_Fixed_Boolean

      LOGICAL ::             DO_SURFACE_EMISSION
      LOGICAL ::             DO_LAMBERTIAN_SURFACE

!  VLIDORT_Fixed_Control

      INTEGER ::             NSTREAMS
      INTEGER ::             NSTOKES

!  VLIDORT_Modified_Sunrays

      INTEGER ::             N_SZANGLES
      DOUBLE PRECISION ::    SZANGLES ( MAX_SZANGLES )

!  VLIDORT_Modified_UserValues

      INTEGER ::             N_USER_RELAZMS
      DOUBLE PRECISION ::    USER_RELAZMS ( MAX_USER_RELAZMS )
      INTEGER ::             N_USER_VZANGLES
      DOUBLE PRECISION ::    USER_VZANGLES ( MAX_USER_VZANGLES )

!  VLIDORT_Modified_Boolean

      LOGICAL ::             DO_USER_VZANGLES
      LOGICAL ::             DO_MVOUT_ONLY

!  Rob Fix 3/18/15. VLIDORT optical Wavelength[Microns]

      DOUBLE PRECISION ::    ATMOS_WAVELENGTH

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

!  VBRDF Inputs
!  ------------

!  VBRDF Control inputs

      BS_DO_USER_STREAMS     = VBRDF_Sup_In%BS_DO_USER_STREAMS
      BS_DO_BRDF_SURFACE     = VBRDF_Sup_In%BS_DO_BRDF_SURFACE
      BS_DO_SURFACE_EMISSION = VBRDF_Sup_In%BS_DO_SURFACE_EMISSION

!  VBRDF Geometry inputs

      BS_NSTOKES             = VBRDF_Sup_In%BS_NSTOKES
      BS_NSTREAMS            = VBRDF_Sup_In%BS_NSTREAMS
      BS_NBEAMS              = VBRDF_Sup_In%BS_NBEAMS
      BS_BEAM_SZAS           = VBRDF_Sup_In%BS_BEAM_SZAS
      BS_N_USER_RELAZMS      = VBRDF_Sup_In%BS_N_USER_RELAZMS
      BS_USER_RELAZMS        = VBRDF_Sup_In%BS_USER_RELAZMS
      BS_N_USER_STREAMS      = VBRDF_Sup_In%BS_N_USER_STREAMS
      BS_USER_ANGLES_INPUT   = VBRDF_Sup_In%BS_USER_ANGLES_INPUT

! Rob Fix 3/18/15. Update for New GCM, 12/28/15

      BS_DO_NewCMGLINT       = VBRDF_Sup_In%BS_DO_NewCMGLINT
      BS_DO_NewGCMGLINT      = VBRDF_Sup_In%BS_DO_NewGCMGLINT
      BS_WAVELENGTH          = VBRDF_Sup_In%BS_WAVELENGTH

!  VLIDORT Main inputs
!  -------------------

!  VLIDORT Fixed Boolean inputs

      DO_SURFACE_EMISSION   = VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION
      DO_LAMBERTIAN_SURFACE = VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE

!  VLIDORT Fixed Control inputs

      NSTOKES    = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTREAMS   = VLIDORT_FixIn%Cont%TS_NSTREAMS

!  VLIDORT Modified Sunrays inputs

      N_SZANGLES = VLIDORT_ModIn%MSunRays%TS_N_SZANGLES
      SZANGLES   = VLIDORT_ModIn%MSunRays%TS_SZANGLES

!  VLIDORT Modified User Value inputs

      N_USER_RELAZMS  = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS    = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS

      N_USER_VZANGLES = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      USER_VZANGLES   = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT

!  VLIDORT Modified Boolean inputs

      DO_USER_VZANGLES = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_MVOUT_ONLY    = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

!  VLIDORT wavelength. Version 2.8, added 3/18/15

      ATMOS_WAVELENGTH = VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = VLIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of BRDF/MAIN compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks

      IF ( BS_DO_BRDF_SURFACE .neqv. (.not.DO_LAMBERTIAN_SURFACE) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'BRDF surface not set for VLIDORT Main'
        ACTIONS(NM)  = 'DO_LAMBERTIAN_SURFACE should be False!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_SURFACE_EMISSION .neqv. DO_SURFACE_EMISSION ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface emission flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_USER_STREAMS .neqv. DO_USER_VZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams/VZangles flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_USER_STREAMS .neqv. (.not.DO_MVOUT_ONLY) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams and DO_MVOUT_ONLY not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_NSTREAMS .ne. NSTREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of discrete ordinates does not agree'
        ACTIONS(NM)  = 'Check NSTREAMS input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_NSTOKES .ne. NSTOKES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Stokes components does not agree'
        ACTIONS(NM)  = 'Check NSTOKES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

!  Angles

      IF ( BS_NBEAMS .ne. N_SZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Solar beams does not agree'
        ACTIONS(NM)  = 'Check BS_NBEAMS and N_SZANGLES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_SZANGLES
          if ( BS_BEAM_SZAS(I) .ne. SZANGLES(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Solar beam angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_BEAM_SZAS and SZANGLES input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( BS_N_USER_STREAMS .ne. N_USER_VZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check N_USER_STREAMS and N_USER_VZANGLES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_VZANGLES
          if ( BS_USER_ANGLES_INPUT(I) .ne. USER_VZANGLES(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'View zenith angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_USER_ANGLES & USER_VZANGLES input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( BS_N_USER_RELAZMS .ne. N_USER_RELAZMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing azimuth angles does not agree'
        ACTIONS(NM)  = 'Check BS_N_USER_RELAZMS & N_USER_RELAZMS input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_RELAZMS
          if ( BS_USER_RELAZMS(I) .ne.USER_RELAZMS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Azimuth angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_USER_RELAZMS & USER_RELAZMS input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

!  Rob Fix 3/18/15. Check wavelengths, NewCM/NewGCM Glint only, Update to NewGCM, 12/28/15
!  Rob Fix 3/7/17. Check wavelengths using adjacency (SMALLNUM = 1.0e-09)

      IF ( BS_DO_NewCMGLINT .or. BS_DO_NewGCMGLINT ) THEN
         if ( ABS((BS_WAVELENGTH/ATMOS_WAVELENGTH)-ONE) .gt. SMALLNUM ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'NewCM/NewGCM Glint: Wavelength does not agree with atmospheric optical-property wavelength'
            ACTIONS(NM)  = 'Check input values of BS_WAVELENGTH and ATMOS_WAVELENGTH'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
         endif
      endif

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      VLIDORT_VBRDFCheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      VLIDORT_VBRDFCheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      VLIDORT_VBRDFCheck_Status%TS_CHECKMESSAGES     = MESSAGES
      VLIDORT_VBRDFCheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_VBRDF_INPUT_CHECK

!

      SUBROUTINE VLIDORT_VBRDF_INPUT_CHECK_ERROR ( ERRORFILE, VLIDORT_VBRDFCheck_Status )

!  Module, dimensions and numbers

      USE VBRDF_Sup_aux_m, ONLY : VBRDF_ERRUNIT
      USE VLIDORT_Outputs_def_m, ONLY : VLIDORT_Exception_Handling

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  VLIDORT/VBRDF input check status

      TYPE(VLIDORT_Exception_Handling), intent(in) :: VLIDORT_VBRDFCheck_Status

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = VBRDF_ERRUNIT

!  Write VLIDORT/VBRDF input compatibility errors to VBRDF error file

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*)' FATAL: VLIDORT and VBRDFSup inputs are incompatible'
      WRITE(W,*)'  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)')'    ** Number of messages = ',VLIDORT_VBRDFCheck_Status%TS_NCHECKMESSAGES
      DO N = 1, VLIDORT_VBRDFCheck_Status%TS_NCHECKMESSAGES
        WRITE(W,'(A,I3,A,A)')'Message # ',N,': ',&
          ADJUSTL(TRIM(VLIDORT_VBRDFCheck_Status%TS_CHECKMESSAGES(N)))
        WRITE(W,'(A,I3,A,A)')'Action  # ',N,': ',&
          ADJUSTL(TRIM(VLIDORT_VBRDFCheck_Status%TS_ACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Checking fail: Look at file ' // TRIM(ERRORFILE)
      STOP ' Program stoppage from call to VLIDORT_VBRDF_INPUT_CHECK_ERROR'

      END SUBROUTINE VLIDORT_VBRDF_INPUT_CHECK_ERROR

!

      SUBROUTINE SET_VLIDORT_VBRDF_INPUTS ( &
        VBRDF_Sup_Out, VLIDORT_FixIn, VLIDORT_ModIn, & !Inputs
        VLIDORT_Sup )                                  !Outputs

!  This subroutine defines the main VLIDORT VBRDF inputs using the corresponding
!  VLIDORT VBRDF supplement outputs (std only)

!  Use Modules

      USE VBRDF_Sup_Outputs_def_m

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m

      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Inputs

      TYPE(VBRDF_Sup_Outputs), INTENT(IN)       :: VBRDF_Sup_Out
      TYPE(VLIDORT_Fixed_Inputs), INTENT(IN)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(IN) :: VLIDORT_ModIn

!  Outputs

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT)    :: VLIDORT_Sup

!  Error output

      !LOGICAL, INTENT(OUT)                     :: FAIL
      !INTEGER, INTENT(INOUT)                   :: N_MESSAGES
      !CHARACTER (LEN=*), INTENT(INOUT)         :: MESSAGES ( MAX_MESSAGES )

!  ---------------
!  Local variables
!  ---------------

      INTEGER :: NSTOKES, N_SZANGLES, NUSERS, NAZIMS, NDISOS, NMOMS
      LOGICAL :: DO_USER_STREAMS, DO_SURFACE_EMISSION

!  Start program

!  Check some inputs (none at present)

      !FAIL = .FALSE.

!  Define some local variables

      DO_USER_STREAMS     = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_SURFACE_EMISSION = VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION

      NSTOKES    = VLIDORT_FixIn%Cont%TS_NSTOKES
      N_SZANGLES = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      NDISOS     = VLIDORT_FixIn%Cont%TS_NSTREAMS
      NMOMS      = 2*NDISOS - 1

      IF ( DO_USER_STREAMS ) THEN
        NUSERS = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
        NAZIMS = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      ENDIF

!  Set VLIDORT standard VBRDF inputs

      VLIDORT_Sup%BRDF%TS_BRDF_F_0(0:NMOMS,:,1:NDISOS,1:N_SZANGLES)  = &
        VBRDF_Sup_Out%BS_BRDF_F_0 (0:NMOMS,:,1:NDISOS,1:N_SZANGLES)
      VLIDORT_Sup%BRDF%TS_BRDF_F  (0:NMOMS,:,1:NDISOS,1:NDISOS) = &
        VBRDF_Sup_Out%BS_BRDF_F   (0:NMOMS,:,1:NDISOS,1:NDISOS)
      IF ( DO_SURFACE_EMISSION ) THEN
        VLIDORT_Sup%BRDF%TS_EMISSIVITY(1:NSTOKES,1:NDISOS) = &
          VBRDF_Sup_Out%BS_EMISSIVITY (1:NSTOKES,1:NDISOS)
      ENDIF

      IF ( DO_USER_STREAMS ) THEN
        VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC(:,1:NUSERS,1:NAZIMS,1:N_SZANGLES) = &
          VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC (:,1:NUSERS,1:NAZIMS,1:N_SZANGLES)
        VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0  (0:NMOMS,:,1:NUSERS,1:N_SZANGLES)  = &
          VBRDF_Sup_Out%BS_USER_BRDF_F_0   (0:NMOMS,:,1:NUSERS,1:N_SZANGLES)
        VLIDORT_Sup%BRDF%TS_USER_BRDF_F    (0:NMOMS,:,1:NUSERS,1:NDISOS) = &
          VBRDF_Sup_Out%BS_USER_BRDF_F     (0:NMOMS,:,1:NUSERS,1:NDISOS)
        IF ( DO_SURFACE_EMISSION ) THEN
          VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1:NSTOKES,1:NUSERS) = &
            VBRDF_Sup_Out%BS_USER_EMISSIVITY (1:NSTOKES,1:NUSERS)
        ENDIF
      ENDIF

      END SUBROUTINE SET_VLIDORT_VBRDF_INPUTS

!  Finish Module

      END MODULE vlidort_vbrdf_sup_accessories_m
