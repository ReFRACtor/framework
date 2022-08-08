
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

!  1/31/21, Version 2.8.3. Restrict parameters, otherwise no changes

module VBRDF_Sup_Outputs_def_m

!  This module contains the following structures:

!  VBRDF_Sup_Outputs - Intent(In) for VLIDORT,
!                      Intent(Out) for VBRDFSup

      use vlidort_pars_m, Only : fpk, MAXMOMENTS, MAXSTOKES, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS, &
                                 MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_BRDF_KERNELS, MAX_MESSAGES

      implicit none

! #####################################################################
! #####################################################################

      type VBRDF_Sup_Outputs

!  direct bounce BRDF. 
!    Rob Fix 9/25/14. Name changed from EXACTDB --> DBOUNCE

      REAL(fpk), dimension ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS ) :: BS_DBOUNCE_BRDFUNC

!  Fourier components of BRDF, in the following order (same all threads)
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams ( Truncated direct-bounce )
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )         :: BS_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )       :: BS_BRDF_F
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )   :: BS_USER_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS ) :: BS_USER_BRDF_F

!  Emissivity

      REAL(fpk), dimension ( MAXSTOKES, MAXSTREAMS )       :: BS_EMISSIVITY
      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_STREAMS ) :: BS_USER_EMISSIVITY

!  9/27/14.  Add WSA and BSA output (one SZA only for the latter)
!  10/16/14. Added the Kernel output. Dimensioning 3/1/17

      REAL(fpk) :: BS_WSA_CALCULATED, BS_WSA_KERNELS(MAX_BRDF_KERNELS)
      REAL(fpk) :: BS_BSA_CALCULATED, BS_BSA_KERNELS(MAX_BRDF_KERNELS)

      end type VBRDF_Sup_Outputs

! #####################################################################
! #####################################################################

      TYPE VBRDF_Input_Exception_Handling


!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: BS_STATUS_INPUTREAD
      INTEGER      :: BS_NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_INPUTACTIONS


      END TYPE VBRDF_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE VBRDF_Output_Exception_Handling

!  Exception handling for Output. New code, 02 April 2014
!     Message Length should be at least 120 Characters

      INTEGER      :: BS_STATUS_OUTPUT
      INTEGER      :: BS_NOUTPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_OUTPUTMESSAGES


      END TYPE VBRDF_Output_Exception_Handling

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VBRDF_Sup_Outputs,              &
             VBRDF_Input_Exception_Handling, &
             VBRDF_Output_Exception_Handling

end module VBRDF_Sup_Outputs_def_m
