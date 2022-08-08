
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

      MODULE VLIDORT_LinInputs_def_m

!  This Module contains the following VLIDORT input Structures, with Intents :

!        VLIDORT_Fixed_LinControl    nested in VLIDORT_Fixed_LinInputs
!        VLIDORT_Fixed_LinOptical    nested in VLIDORT_Fixed_LinInputs
!         VLIDORT_Fixed_LinInputs    Intent(In)

!     VLIDORT_Modified_LinControl    nested in VLIDORT_Modified_LinInputs
!      VLIDORT_Modified_LinInputs    Intent(InOut)

      USE VLIDORT_PARS_m, Only : fpk, MAX_ATMOSWFS, MAX_GEOMETRIES, MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_LinControl

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, dimension(MAXLAYERS)  :: TS_LAYER_VARY_FLAG
      INTEGER, dimension(MAXLAYERS)  :: TS_LAYER_VARY_NUMBER

!  Total number of column Jacobians

      INTEGER  :: TS_N_TOTALCOLUMN_WFS

!  Total number of profile Jacobians

      INTEGER  :: TS_N_TOTALPROFILE_WFS

!  Total number of surface Jacobians

      INTEGER  :: TS_N_SURFACE_WFS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      INTEGER  :: TS_N_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Jacobian names

      CHARACTER (LEN=31), dimension(MAX_ATMOSWFS) :: TS_COLUMNWF_NAMES
      CHARACTER (LEN=31), dimension(MAX_ATMOSWFS) :: TS_PROFILEWF_NAMES

      END TYPE VLIDORT_Fixed_LinControl

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_LinOptical


!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk), dimension(MAX_ATMOSWFS, MAXLAYERS) :: TS_L_DELTAU_VERT_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS, MAXLAYERS) :: TS_L_OMEGA_TOTAL_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ) :: TS_L_GREEKMAT_TOTAL_INPUT

!  Linearized Fmatrix inputs for FO calculations
!mick fix 9/19/2017 - swapped layer & geo indices

      REAL(fpk), dimension( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, 6 ) :: TS_L_FMATRIX_UP 
      REAL(fpk), dimension( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, 6 ) :: TS_L_FMATRIX_DN 

      END TYPE VLIDORT_Fixed_LinOptical

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_LinInputs


      TYPE(VLIDORT_Fixed_LinControl)    :: Cont
      TYPE(VLIDORT_Fixed_LinOptical)    :: Optical


      END TYPE VLIDORT_Fixed_LinInputs

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_LinControl

!  Control linearization

      LOGICAL  :: TS_DO_COLUMN_LINEARIZATION
      LOGICAL  :: TS_DO_PROFILE_LINEARIZATION
      LOGICAL  :: TS_DO_ATMOS_LINEARIZATION

      LOGICAL  :: TS_DO_SURFACE_LINEARIZATION
      LOGICAL  :: TS_DO_LINEARIZATION

!  This flag moved from Fixed Lin Control, Version 2.7

      LOGICAL  :: TS_DO_SIMULATION_ONLY

!  BlackBody Jacobian Flags, Introduced March 26th 2014, Version 2.7

      LOGICAL  :: TS_DO_ATMOS_LBBF
      LOGICAL  :: TS_DO_SURFACE_LBBF

!  These two flags have been superseded in Version 2.7
!      LOGICAL  :: TS_DO_SURFBB_LINEARIZATION  !  REPLACED BY TS_DO_SURFACE_LBBF
!      LOGICAL  :: TS_DO_LTE_LINEARIZATION     !  REPLACED BY TS_DO_ATMOS_LBBF 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      LOGICAL  :: TS_DO_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      END TYPE VLIDORT_Modified_LinControl

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_LinInputs


      TYPE(VLIDORT_Modified_LinControl)    :: MCont


      END TYPE VLIDORT_Modified_LinInputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_Fixed_LinControl,    &
                VLIDORT_Fixed_LinOptical,    &
                VLIDORT_Fixed_LinInputs,     &
                VLIDORT_Modified_LinControl, &
                VLIDORT_Modified_LinInputs

   END MODULE VLIDORT_LinInputs_def_m

