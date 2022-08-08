
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
! #            VLIDORT_UNPACK_L_THERM                           #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. No Changes.

      MODULE vlidort_l_unpack_m

      PRIVATE
      PUBLIC :: VLIDORT_UNPACK_L_THERM

      CONTAINS

      SUBROUTINE VLIDORT_UNPACK_L_THERM ( L_Therm,                                  & ! Input
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, N_TOTALATMOS_WFS, & ! Input
        L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,                       & ! Output
        L_T_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN )            ! Output

      USE VLIDORT_PARS_m, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS
 
      USE VLIDORT_LinWork_def_m

!  Input
!  -----

      TYPE(VLIDORT_LinWork_Thermal), INTENT (IN) :: L_Therm

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          N_TOTALATMOS_WFS

!  Output
!  ------

!  From THERMAL_SETUP_PLUS

      DOUBLE PRECISION, INTENT (OUT) :: L_THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DIRECT_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DIRECT_DN  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Unpack structure

              L_THERMCOEFFS ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_THERMCOEFFS ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS )
              L_DELTAU_POWER ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_DELTAU_POWER ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS )
              L_XTAU_POWER ( 1:N_PARTLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_XTAU_POWER ( 1:N_PARTLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS )
              L_TCOM1 ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_TCOM1 ( 1:NLAYERS, 1:N_THERMAL_COEFFS, 1:N_TOTALATMOS_WFS )

              L_T_DIRECT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_T_DIRECT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_TOTALATMOS_WFS )
              L_T_DIRECT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_T_DIRECT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_TOTALATMOS_WFS )

              L_T_UT_DIRECT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_T_UT_DIRECT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_TOTALATMOS_WFS )
              L_T_UT_DIRECT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_TOTALATMOS_WFS ) = &
      L_Therm%L_T_UT_DIRECT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_TOTALATMOS_WFS )

      END SUBROUTINE VLIDORT_UNPACK_L_THERM

      END MODULE vlidort_l_unpack_m

