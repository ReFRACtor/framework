
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
! #            VLIDORT_L_INPUT_MASTER (master), calling:        #
! #             (4 routines from vlidort_inputs)                #
! #            VLIDORT_L_INIT_INPUTS                            #
! #            VLIDORT_LinSup_Init                              #
! #            VLIDORT_BRDF_LinSup_Init                         #
! #            VLIDORT_SLEAVE_LinSup_Init                       #
! #            VLIDORT_SS_LinSup_Init                           #
! #                                                             #
! #            VLIDORT_L_READ_INPUTS                            #
! #                                                             #
! #    These routines called by Main VLIDORT_LPS/LCS module     #
! #                                                             #
! #            VLIDORT_L_CHECK_INPUT_DIMS                       #
! #            VLIDORT_L_CHECK_INPUT                            #
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

!  1/31/21. Version 2.8.3. This Module, no actual code changes.

      MODULE vlidort_l_inputs_m

      !PRIVATE
      PUBLIC

      CONTAINS

      SUBROUTINE VLIDORT_L_INPUT_MASTER ( &
        FILNAM,             & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_LinFixIn,   & ! Outputs
        VLIDORT_LinModIn,   & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_USER_LEVELS, &
                                 MAX_USER_OBSGEOMS, MAX_ATMOSWFS, MAXLAYERS, MAX_MESSAGES, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_INUNIT

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_LinInputs_def_m
      USE VLIDORT_Outputs_def_m

      USE VLIDORT_INPUTS_m, Only : VLIDORT_INIT_INPUTS, VLIDORT_READ_INPUTS

      IMPLICIT NONE

!  Inputs
!  ------

      CHARACTER (LEN=*), INTENT (IN) ::   FILNAM

!  Outputs
!  -------

      TYPE(VLIDORT_Fixed_Inputs), INTENT (OUT)             :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (OUT)          :: VLIDORT_ModIn
      TYPE(VLIDORT_Fixed_LinInputs), INTENT (OUT)          :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT (OUT)       :: VLIDORT_LinModIn
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

!  Initialize standard variables
!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added (1/31/16)
!     Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.
!     Version 2.8: SS flags upgraded. 9/19/16. Argument list rearranged.
!                  VLIDORT input type structures now used for argument passing - 9/19/2017
!                  All VLIDORT fixed and modified inputs now initialized - 9/19/2017

      CALL VLIDORT_INIT_INPUTS ( &
         VLIDORT_FixIn, VLIDORT_ModIn ) !Outputs

!  Initialize linearization variables
!     Version 2.8: VLIDORT input type structures now used for argument passing - 9/19/2017
!                  All VLIDORT linearized fixed and modified inputs now initialized - 9/19/2017

      CALL VLIDORT_L_INIT_INPUTS ( &
         VLIDORT_LinFixIn, VLIDORT_LinModIn ) !Outputs

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

!  Read file error - standard variables

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
         STATUS = VLIDORT_SERIOUS
         VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
         VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
         VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
         VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS
         CLOSE(FILUNIT)
         RETURN
      ENDIF

!  Version 3.8. Introduced following routine to Read linearization inputs from File
!     -- Mirrors similar routine already in the VLIDORT code

      CALL VLIDORT_L_READ_INPUTS ( &
         VLIDORT_FixIn, VLIDORT_ModIn,       & !Inputs
         VLIDORT_LinFixIn, VLIDORT_LinModIn, & !InOut
         STATUS_SUB,                         & !Outputs
         NMESSAGES, MESSAGES, ACTIONS )        !InOut

!  Read file error - linearized variables

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

      END SUBROUTINE VLIDORT_L_INPUT_MASTER

!

      SUBROUTINE VLIDORT_L_INIT_INPUTS &
      ( VLIDORT_LinFixIn, VLIDORT_LinModIn ) !Outputs

!  VLIDORT linearized input type structures now used for argument passing - 9/19/2017
!  All VLIDORT linearized fixed and modified inputs now initialized - 9/19/2017

!  Initialises all linearized inputs for VLIDORT
!  ---------------------------------------------

!  Module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXLAYERS, ZERO
      USE VLIDORT_LinInputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT linearized input structures

      TYPE(VLIDORT_Fixed_LinInputs), INTENT (OUT)    :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT (OUT) :: VLIDORT_LinModIn

!  Initialize linearized control variables
!  =======================================

!  VLIDORT Fixed LinControl
!  ------------------------

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG    = .FALSE.
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER  = 0

      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS  = 0
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = 0
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS      = 0
      VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS       = 0

      VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES     = ' '
      VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES    = ' '

!  VLIDORT Modified LinControl
!  ---------------------------

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .FALSE.
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .FALSE.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .FALSE.

      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .FALSE.
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = .FALSE.

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .FALSE.

      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF            = .FALSE.
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF          = .FALSE.
      VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS            = .FALSE.

!  VLIDORT Fixed LinOptical
!  ------------------------

      VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT     = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT     = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT  = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_FMATRIX_UP            = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_FMATRIX_DN            = ZERO

! Finish

      RETURN
      END SUBROUTINE VLIDORT_L_INIT_INPUTS

!

      SUBROUTINE VLIDORT_LinSup_Init ( VLIDORT_LinSup )

      USE VLIDORT_LinSup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT linearized supplement input structure

      TYPE(VLIDORT_LinSup_InOut), INTENT(INOUT) :: VLIDORT_LinSup

!  Initialize VLIDORT linearized supplement inputs
!  ===============================================

      CALL VLIDORT_BRDF_LinSup_Init   ( VLIDORT_LinSup )
      CALL VLIDORT_SLEAVE_LinSup_Init ( VLIDORT_LinSup )
      CALL VLIDORT_SS_LinSup_Init     ( VLIDORT_LinSup )

!  Finish

      END SUBROUTINE VLIDORT_LinSup_Init

!

      SUBROUTINE VLIDORT_BRDF_LinSup_Init ( VLIDORT_LinSup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_LinSup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT linearized supplement input structure

      TYPE(VLIDORT_LinSup_InOut), INTENT(INOUT) :: VLIDORT_LinSup

!  Initialize VLIDORT linearized brdf supplement inputs
!  ====================================================

      VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO

      VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO

!  Finish

      END SUBROUTINE VLIDORT_BRDF_LinSup_Init

!

      SUBROUTINE VLIDORT_SLEAVE_LinSup_Init ( VLIDORT_LinSup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_LinSup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT linearized supplement input structure

      TYPE(VLIDORT_LinSup_InOut), INTENT(INOUT) :: VLIDORT_LinSup

!  Initialize VLIDORT linearized sleave supplement inputs
!  ======================================================

      VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
      VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
      VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
      VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO

!  Finish

      END SUBROUTINE VLIDORT_SLEAVE_LinSup_Init

!

SUBROUTINE VLIDORT_SS_LinSup_Init ( VLIDORT_LinSup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_LinSup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT linearized supplement input structure

      TYPE(VLIDORT_LinSup_InOut), INTENT(INOUT) :: VLIDORT_LinSup

!  Initialize VLIDORT linearized single-scatter supplement inputs
!  ==============================================================

      VLIDORT_LinSup%SS%Col%TS_COLUMNWF_SS   = ZERO
      VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB   = ZERO
      VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS = ZERO
      VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB = ZERO

      VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB = ZERO

!  Finish

END SUBROUTINE VLIDORT_SS_LinSup_Init

!

SUBROUTINE VLIDORT_L_READ_INPUTS ( &
        VLIDORT_FixIn, VLIDORT_ModIn,       & !Inputs
        VLIDORT_LinFixIn, VLIDORT_LinModIn, & !InOut
        STATUS,                             & !Outputs
        NMESSAGES, MESSAGES, ACTIONS )        !InOut

!  Version 2.8: VLIDORT input type structures now used for argument passing - 9/19/2017

!  Module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_ATMOSWFS, MAX_MESSAGES, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_INUNIT

      USE VLIDORT_AUX_m , Only : LEN_STRING, GFINDPAR, FINDPAR_ERROR

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_LinInputs_def_m

!  Inputs
!  ------

      TYPE(VLIDORT_Fixed_Inputs), INTENT(IN)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(IN) :: VLIDORT_ModIn

!  Outputs
!  -------

      TYPE(VLIDORT_Fixed_LinInputs), INTENT(INOUT)    :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT(INOUT) :: VLIDORT_LinModIn

!  Exception handling

      INTEGER, INTENT (OUT)   ::            STATUS
      INTEGER, INTENT (INOUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGES(0:MAX_MESSAGES)
      CHARACTER (LEN=*), INTENT (INOUT) ::  ACTIONS (0:MAX_MESSAGES)

!  Local variables
!  ---------------

!  General linearization flags

      LOGICAL ::            DO_COLUMN_LINEARIZATION
      LOGICAL ::            DO_PROFILE_LINEARIZATION
      LOGICAL ::            DO_ATMOS_LINEARIZATION
      LOGICAL ::            DO_SURFACE_LINEARIZATION
      !LOGICAL ::            DO_SLEAVE_WFS
      LOGICAL ::            DO_LINEARIZATION

      LOGICAL ::            DO_SIMULATION_ONLY

!  LBBF flags, New for 2p7 variables

      LOGICAL ::            DO_ATMOS_LBBF
      LOGICAL ::            DO_SURFACE_LBBF

!  Atmospheric Jacobian control

      INTEGER ::            N_TOTALPROFILE_WFS
      INTEGER ::            N_TOTALCOLUMN_WFS
      LOGICAL ::            LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER ::            LAYER_VARY_NUMBER ( MAXLAYERS )
      !INTEGER ::            N_SURFACE_WFS
      !INTEGER ::            N_SLEAVE_WFS

!  Names

      CHARACTER (LEN=31) ::  PROFILEWF_NAMES ( MAX_ATMOSWFS )
      CHARACTER (LEN=31) ::  COLUMNWF_NAMES  ( MAX_ATMOSWFS )

!  Other variables
!  ---------------

      INTEGER ::            NLAYERS
      LOGICAL ::            DO_THERMAL_TRANSONLY

      CHARACTER (LEN=9), PARAMETER :: PREFIX = 'VLIDORT -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            FILUNIT, I, NM

!  Initialize Exception handling

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

!  Proxy variables

      NLAYERS              = VLIDORT_FixIn%Cont%TS_NLAYERS
      DO_THERMAL_TRANSONLY = VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY

!  Linearization control

      PAR_STR = 'Do simulation only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SIMULATION_ONLY
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do atmospheric profile weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_PROFILE_LINEARIZATION
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do atmospheric column weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_COLUMN_LINEARIZATION
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do surface property weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_LINEARIZATION
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  New section, LBBF linearization flags. Version 2.7
!  No restrictions on usage of these options !!!!!

      PAR_STR = 'Atmospheric BB emission weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_ATMOS_LBBF
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Surface BB emission weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_LBBF
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Options currently disabled. 07 August 2014

      IF ( DO_ATMOS_LBBF .or. DO_SURFACE_LBBF ) then
         NM = NM + 1
         STATUS       = VLIDORT_SERIOUS
         MESSAGES(NM) = 'DO_ATMOS_LBBF and DO_SURFACE_LBBF not enabled 8/7/14'
         ACTIONS(NM)  = 'Turn both flags off: DO_ATMOS_LBBF, DO_SURFACE_LBBF'
         NMESSAGES    = NM
         RETURN
      ENDIF

!  LTE linearization flag. 2.6 code Replaced in Version 2.7
!      PAR_STR = 'Do LTE temperature weighting functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) DO_LTE_LINEARIZATION
!      CALL FINDPAR_ERROR &
!        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check LTE linearization. Version 2.6 code superceded.
!      IF ( DO_LTE_LINEARIZATION ) THEN
!        IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Must have THERMAL_transonly for LTE linearization'
!          ACTIONS(NM)  = 'Re-set input value'
!          STATUS  = VLIDORT_SERIOUS
!          NMESSAGES = NM
!          RETURN
!        ELSE
!          IF( DO_PROFILE_LINEARIZATION .OR. DO_COLUMN_LINEARIZATION )THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'No other kind of atmospheric linearization allowed'
!            ACTIONS(NM)  = 'Re-set input value'
!            STATUS  = VLIDORT_SERIOUS
!            NMESSAGES = NM
!            RETURN
!          ENDIF
!        ENDIF
!      ENDIF

!  Disabled Version 2.6 code
!      PAR_STR = 'Surface emission weighting functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SURFBB_LINEARIZATION
!      CALL FINDPAR_ERROR  ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Atmospheric profile weighting function input control
!  Individual layer number does not exceed whole number

      IF ( DO_PROFILE_LINEARIZATION ) THEN

        PAR_STR='Number of atmospheric profile weighting functions (total)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) N_TOTALPROFILE_WFS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  This code should be prepared, not taken from file-read
!        PAR_STR = 'Atmospheric weighting functions, layer control'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!          DO I = 1, NLAYERS
!            READ (FILUNIT,*,ERR=998) LAYER_VARY_FLAG(I), LAYER_VARY_NUMBER(I)
!          ENDDO
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      ENDIF

!  Atmospheric Bulk/column weighting function input control

       IF ( DO_COLUMN_LINEARIZATION ) THEN

         PAR_STR='Number of atmospheric column weighting functions (total)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) N_TOTALCOLUMN_WFS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
         DO I = 1, NLAYERS
           LAYER_VARY_FLAG(I)   = .TRUE.
           LAYER_VARY_NUMBER(I) = N_TOTALCOLUMN_WFS
         ENDDO
       ENDIF

!  Weighting function names

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        PAR_STR='Atmospheric profile Jacobian names (character*31)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALPROFILE_WFS
            READ (FILUNIT,'(a31)',ERR=998) PROFILEWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        PAR_STR='Atmospheric column Jacobian names (character*31)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALCOLUMN_WFS
            READ (FILUNIT,'(a31)',ERR=998) COLUMNWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Set overall linearization flags

      DO_ATMOS_LINEARIZATION = ( DO_PROFILE_LINEARIZATION .OR. DO_COLUMN_LINEARIZATION .OR. &
                                 DO_ATMOS_LBBF )

      DO_LINEARIZATION       = ( DO_ATMOS_LINEARIZATION .OR. DO_SURFACE_LINEARIZATION .OR. &
                                 DO_SURFACE_LBBF )

!  Define VLIDORT lin inputs
!  =========================

!mick mod 9/19/2017 - initializing of all lin VLIDORT type structure input variables is done
!                     in subroutine VLIDORT_L_INIT_INPUTS; only those read in from the
!                     VLIDORT config file are possibly modified here.  If conditions are
!                     applied where appropriate.

!  VLIDORT Fixed LinControl

      IF ( DO_COLUMN_LINEARIZATION ) THEN
         VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS                   = N_TOTALCOLUMN_WFS
         VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)          = LAYER_VARY_FLAG(1:NLAYERS)
         VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS)        = LAYER_VARY_NUMBER(1:NLAYERS)
         VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES(1:N_TOTALCOLUMN_WFS) = COLUMNWF_NAMES(1:N_TOTALCOLUMN_WFS)
      ENDIF
      IF ( DO_PROFILE_LINEARIZATION ) THEN
         VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS                    = N_TOTALPROFILE_WFS
         VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES(1:N_TOTALPROFILE_WFS) = PROFILEWF_NAMES(1:N_TOTALPROFILE_WFS)
      ENDIF
      !IF ( DO_SURFACE_LINEARIZATION ) &
      !   VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS = N_SURFACE_WFS
      !IF ( DO_SURFACE_LEAVING ) &
      !   VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS  = N_SLEAVE_WFS

!  VLIDORT Modified LinControl

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION

      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = DO_LINEARIZATION

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = DO_SIMULATION_ONLY

      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF            = DO_ATMOS_LBBF
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF          = DO_SURFACE_LBBF
      !IF ( DO_SURFACE_LEAVING ) &
      !   VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS         = DO_SLEAVE_WFS

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
      END SUBROUTINE VLIDORT_L_READ_INPUTS

!

      SUBROUTINE VLIDORT_L_CHECK_INPUT_DIMS_HOLD &
      ( DO_COLUMN_LINEARIZATION,  DO_PROFILE_LINEARIZATION, & ! Flags
        DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS,            & ! Flags
        N_TOTALCOLUMN_WFS, N_TOTALPROFILE_WFS,              & ! Dimension-checking
        N_SURFACE_WFS, N_SLEAVE_WFS,                        & ! Dimension-checking
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                ! Status and output

!  Check input dimensions

!  module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAX_ATMOSWFS, MAX_SURFACEWFS, &
                                 MAX_SLEAVEWFS, MAX_MESSAGES, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Flags

      LOGICAL  , intent(in) :: DO_COLUMN_LINEARIZATION
      LOGICAL  , intent(in) :: DO_PROFILE_LINEARIZATION
      LOGICAL  , intent(in) :: DO_SURFACE_LINEARIZATION
      LOGICAL  , intent(in) :: DO_SLEAVE_WFS

!  Numbers to be checked for Dimensioning

      INTEGER, intent(inout)   :: N_TOTALCOLUMN_WFS, N_TOTALPROFILE_WFS
      INTEGER, intent(inout)   :: N_SURFACE_WFS,     N_SLEAVE_WFS

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

!  Check LIDORT linearized input dimensions against maximum dimensions
!  ===================================================================

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        IF ( N_TOTALCOLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for column WFs'
         ACTIONS(NM)  = 'Action: Decrease N_TOTALCOLUMN_WFS or increase MAX_ATMOSWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        IF ( N_TOTALPROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for profile WFs'
         ACTIONS(NM)  = 'Action: Decrease N_TOTALPROFILE_WFS or increase MAX_ATMOSWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_SURFACE_LINEARIZATION ) THEN
        IF ( N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = 'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACE_WFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_SLEAVE_WFS ) THEN
        IF ( N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) =  'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  =  'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      END SUBROUTINE VLIDORT_L_CHECK_INPUT_DIMS_HOLD

!

      SUBROUTINE VLIDORT_L_CHECK_INPUT_DIMS &
        ( VLIDORT_LinFixIn, VLIDORT_LinModIn, &
          STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAX_ATMOSWFS, MAX_SURFACEWFS, &
                                 MAX_SLEAVEWFS, MAX_MESSAGES, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE VLIDORT_LinInputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_LinInputs)   , INTENT (IN) :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT (IN) :: VLIDORT_LinModIn

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

!  Check VLIDORT linearized input dimensions
!    against maximum dimensions
!  =========================================

      IF ( VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION ) THEN
        IF ( VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for column WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_TOTALCOLUMN_WFS or increase MAX_ATMOSWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION ) THEN
        IF ( VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for profile WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_TOTALPROFILE_WFS or increase MAX_ATMOSWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION ) THEN
        IF ( VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACE_WFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS ) THEN
        IF ( VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      END SUBROUTINE VLIDORT_L_CHECK_INPUT_DIMS

!

      SUBROUTINE VLIDORT_L_CHECK_INPUT ( &
        DO_SIMULATION_ONLY, N_SURFACE_WFS, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        DO_ATMOS_LINEARIZATION, DO_SURFACE_LINEARIZATION, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  1/31/21. Version 2.8.3. Add MAX_SURFACEWFS to dimensioning list

      USE VLIDORT_PARS_m, Only : MAX_MESSAGES, MAX_SURFACEWFS, VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_WARNING

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::            DO_SIMULATION_ONLY
      INTEGER, INTENT (IN) ::            N_SURFACE_WFS

      LOGICAL, INTENT (INOUT) ::         DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::         DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::         DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::         DO_SURFACE_LINEARIZATION

      INTEGER, INTENT (OUT) ::           STATUS

!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%
      INTEGER, INTENT (INOUT) ::           NMESSAGES
!   Messages and actions should be INTENT(inout)
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER (LEN=*), INTENT (INOUT) :: ACTIONS (0:MAX_MESSAGES)
!      INTEGER, INTENT (OUT) ::           NMESSAGES
!      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGES(0:MAX_MESSAGES)
!      CHARACTER (LEN=*), INTENT (OUT) :: ACTIONS (0:MAX_MESSAGES)
!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%

!  local variables

      INTEGER :: NM

!  Initialize output status

      STATUS = VLIDORT_SUCCESS

!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%
!   Messages and actions should be INTENT(inout)
!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of VLIDORT Basic Input'
!      ACTIONS(0)      = 'No Action required for this Task'
!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%

!  %%% Initialize local message count to 0, 05 March 2013
!      NM     = NMESSAGES
      NM     = 0

!  Do some checking
!  ================

!  Check formatting; should not exceed 10 entries for the number of
!  layer weighting functions. This is needed to avoid errors with
!  compilers not using variable formatting.

!      IF ( MAX_ATMOSWFS .GT. 10 ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'Max layer variations > 10, get formatting problems!'
!        ACTIONS(NM)  = 'Re-set number to avoid this'
!        STATUS = VLIDORT_SERIOUS
!      ENDIF

!  1/31/21. Version 2.8.3. Add Check on N_SURFACE_WFS

      IF ( DO_SURFACE_LINEARIZATION ) THEN
        IF ( N_SURFACE_WFS .gt. MAX_SURFACEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of surface weighting functions exceeds dimensioning'
          ACTIONS(NM)  = 'Either re-set dimensioning MAX_SURFACEWFS, or re-set input values'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Check something is being varied

      IF ( .NOT. DO_SIMULATION_ONLY ) THEN
        IF ( .NOT. DO_ATMOS_LINEARIZATION .AND. .NOT. DO_SURFACE_LINEARIZATION ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Nothing being varied or calculated'
          ACTIONS(NM)  = 'Re-set input values'
          STATUS = VLIDORT_SERIOUS
!          NMESSAGES = NM             ! @@@ Robfix 01 aug 2012, Miscount NMESSAGES
        ENDIF
      ENDIF

!  Check and fix AF flags if you want simulation only

      IF ( DO_SIMULATION_ONLY ) THEN
        IF ( DO_ATMOS_LINEARIZATION ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Input warning; simulation only ---> NO Atmos WFs'
          ACTIONS(NM)  = 'Action: VLIDORT switched off WF flag for you!'
          STATUS = VLIDORT_WARNING
          DO_PROFILE_LINEARIZATION = .FALSE.
          DO_COLUMN_LINEARIZATION  = .FALSE.
          DO_ATMOS_LINEARIZATION   = .FALSE.
        ENDIF
        IF ( DO_SURFACE_LINEARIZATION ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Input warning; simulation only ---> NO Surface WFs'
          ACTIONS(NM)  = 'Action: VLIDORT switched off WF flag for you!'
          STATUS = VLIDORT_WARNING
          DO_SURFACE_LINEARIZATION = .FALSE.
        ENDIF
      ENDIF

!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%
!  Number of messages
      NMESSAGES = NMESSAGES + NM
!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_CHECK_INPUT

!  End Module

      END MODULE vlidort_l_inputs_m

