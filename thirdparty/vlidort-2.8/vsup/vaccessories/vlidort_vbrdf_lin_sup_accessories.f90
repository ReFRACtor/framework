
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
! #            SET_VLIDORT_L_VBRDF_INPUTS                       #
! #                                                             #
! ###############################################################

!  Upgrade to VLIDORT Version 2.8
!  -----------------------------

!  Mick Mod 9/19/17. This module new.
!  Added subroutine SET_VLIDORT_L_VBRDF_INPUTS:
!    ** Subroutine to define the main VLIDORT linearized VBRDF inputs
!       by defining them using the corresponding VLIDORT linearized
!       VBRDF supplement outputs 

      MODULE vlidort_vbrdf_LinSup_accessories_m

      PRIVATE
      PUBLIC :: SET_VLIDORT_L_VBRDF_INPUTS

      CONTAINS

      SUBROUTINE SET_VLIDORT_L_VBRDF_INPUTS ( &
        VBRDF_Sup_Out, VBRDF_LinSup_Out,      & !Inputs
        VLIDORT_FixIn, VLIDORT_ModIn,         & !Inputs
        VLIDORT_LinFixIn, VLIDORT_LinModIn,   & !Inputs
        VLIDORT_Sup, VLIDORT_LinSup)            !Outputs

!  This subroutine defines the main VLIDORT VBRDF inputs using the corresponding
!  VLIDORT VBRDF supplement outputs (std & lin)

!  Use Modules

      USE VBRDF_Sup_Outputs_def_m
      USE VBRDF_LinSup_Outputs_def_m

      USE VLIDORT_PARS_m

      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_LIN_IO_DEFS_m

      USE VLIDORT_Sup_InOut_def_m
      USE VLIDORT_LinSup_InOut_def_m

      IMPLICIT NONE

!  Inputs

      TYPE(VBRDF_Sup_Outputs), INTENT(IN)          :: VBRDF_Sup_Out
      TYPE(VBRDF_LinSup_Outputs), INTENT(IN)       :: VBRDF_LinSup_Out
      TYPE(VLIDORT_Fixed_Inputs), INTENT(IN)       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(IN)    :: VLIDORT_ModIn
      TYPE(VLIDORT_Fixed_LinInputs), INTENT(IN)    :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT(IN) :: VLIDORT_LinModIn

!  Outputs

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT)    :: VLIDORT_Sup
      TYPE(VLIDORT_LinSup_InOut), INTENT(INOUT) :: VLIDORT_LinSup

!  Error output

      !LOGICAL, INTENT(OUT)                     :: FAIL
      !INTEGER, INTENT(INOUT)                   :: N_MESSAGES
      !CHARACTER (LEN=*), INTENT(INOUT)         :: MESSAGES ( MAX_MESSAGES )

!  ---------------
!  Local variables
!  ---------------

      INTEGER :: NSTOKES, N_SZANGLES, NUSERS, NAZIMS, NDISOS, NMOMS, NWFS
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

!  Set VLIDORT linearized VBRDF inputs

      IF ( VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION ) THEN
        NWFS = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS

        VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0(1:NWFS,0:NMOMS,:,1:NDISOS,1:N_SZANGLES)  = &
          VBRDF_LinSup_Out%BS_LS_BRDF_F_0 (1:NWFS,0:NMOMS,:,1:NDISOS,1:N_SZANGLES)
        VLIDORT_LinSup%BRDF%TS_LS_BRDF_F  (1:NWFS,0:NMOMS,:,1:NDISOS,1:NDISOS) = &
          VBRDF_LinSup_Out%BS_LS_BRDF_F   (1:NWFS,0:NMOMS,:,1:NDISOS,1:NDISOS)
        IF ( DO_SURFACE_EMISSION ) THEN
          VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY(1:NWFS,1:NSTOKES,1:NDISOS) = &
            VBRDF_LinSup_Out%BS_LS_EMISSIVITY (1:NWFS,1:NSTOKES,1:NDISOS)
        ENDIF

        IF ( DO_USER_STREAMS ) THEN
          VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC(1:NWFS,:,1:NUSERS,1:NAZIMS,1:N_SZANGLES) = &
            VBRDF_LinSup_Out%BS_LS_DBOUNCE_BRDFUNC (1:NWFS,:,1:NUSERS,1:NAZIMS,1:N_SZANGLES)
          VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0  (1:NWFS,0:NMOMS,:,1:NUSERS,1:N_SZANGLES)  = &
            VBRDF_LinSup_Out%BS_LS_USER_BRDF_F_0   (1:NWFS,0:NMOMS,:,1:NUSERS,1:N_SZANGLES)
          VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:NWFS,0:NMOMS,:,1:NUSERS,1:NDISOS) = &
            VBRDF_LinSup_Out%BS_LS_USER_BRDF_F     (1:NWFS,0:NMOMS,:,1:NUSERS,1:NDISOS)
          IF ( DO_SURFACE_EMISSION ) THEN
            VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY(1:NWFS,1:NSTOKES,1:NUSERS) = &
              VBRDF_LinSup_Out%BS_LS_USER_EMISSIVITY (1:NWFS,1:NSTOKES,1:NUSERS)
          ENDIF
        ENDIF
      ENDIF

      END SUBROUTINE SET_VLIDORT_L_VBRDF_INPUTS

!  Finish Module

      END MODULE vlidort_vbrdf_LinSup_accessories_m
