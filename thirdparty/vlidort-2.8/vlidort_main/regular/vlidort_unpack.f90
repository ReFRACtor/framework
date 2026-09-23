
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
! #            VLIDORT_UNPACK_MISC                              #
! #            VLIDORT_UNPACK_THERM                             #
! #            VLIDORT_UNPACK_MULT                              #
! #                                                             #
! ###############################################################

      MODULE vlidort_unpack_m

!  1/31/21. Version 2.8.3.
!    -- Need additional variables SIGMA_M, SIGMA_P, ITRANS_USERM to be unpacked
!    -- Additional outputs to VLIDORT_UNPACK_MISC and VLIDORT_UNPACK_MULT

      PUBLIC :: VLIDORT_UNPACK_MISC,  &
                VLIDORT_UNPACK_THERM, &
                VLIDORT_UNPACK_MULT

      CONTAINS

      SUBROUTINE VLIDORT_UNPACK_MISC ( Misc,                       & ! Input structure to be unpacked
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS,                  & ! Input Numbers
        N_USER_VZANGLES, N_SZANGLES, NMOMENTS,                     & ! Input Numbers
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT,                    & ! unpacked Optical
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, OMEGA_GREEK,       & ! unpacked Optical
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,            & ! unpacked Trans D.O.
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, CUMTRANS,        & ! unpacked Trans User
        BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,    & ! unpacked Solar
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR, & ! unpacked Average Secant
        ITRANS_USERM, LOCAL_CSZA )                                   ! unpacked Misc

!  1/31/21. Version 2.8.3.
!    -- Argument list rearranged to follow more natural logic
!    -- ITRANS_USERM added to the unpack list

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAX_USER_STREAMS, MAXMOMENTS, MAXLAYERS, MAX_PARTLAYERS, MAXSTOKES, MAXBEAMS
      USE VLIDORT_Work_def_m

!  Input
!  -----

      TYPE(VLIDORT_Work_Miscellanous), INTENT (IN) :: Misc

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          NMOMENTS

!  Output
!  ------

!  From VLIDORT_MISCSETUPS, optical variables
!mick mod 9/19/2017 - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to output

      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  User-angle transmittances

      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (OUT) :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )

!  Solar

      INTEGER, INTENT (OUT) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL, INTENT (OUT) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Average-secant parameterization

      DOUBLE PRECISION, INTENT (OUT) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

! 1/31/21. Version 2.8.3. Green's function requires additional array ITRANS_USERM to be unpacked

      DOUBLE PRECISION, INTENT (OUT) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Unpack structure
!  ----------------

!  optical

           DELTAU_VERT ( 1:NLAYERS ) = &
      Misc%DELTAU_VERT ( 1:NLAYERS )
           PARTAU_VERT ( 1:N_PARTLAYERS ) = &
      Misc%PARTAU_VERT ( 1:N_PARTLAYERS )
           DELTAU_SLANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES ) = &
      Misc%DELTAU_SLANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES )

           LEVELS_SOLARTRANS ( 0:NLAYERS, 1:N_SZANGLES )        = &
      Misc%LEVELS_SOLARTRANS ( 0:NLAYERS, 1:N_SZANGLES )
           PARTIALS_SOLARTRANS ( 1:N_PARTLAYERS, 1:N_SZANGLES ) = &
      Misc%PARTIALS_SOLARTRANS ( 1:N_PARTLAYERS, 1:N_SZANGLES )

           OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES ) = &
      Misc%OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES )

!  Discrete ordidnate transmittances

           T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS ) = &
      Misc%T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS )
           T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS ) = &
      Misc%T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS )
           T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS ) = &
      Misc%T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS )

!  User angle sutff

           T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES ) = &
      Misc%T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES )
           T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES ) = &
      Misc%T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES )
           T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES ) = &
      Misc%T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES )
           CUMTRANS ( 1:NLAYERS, 1:N_USER_VZANGLES ) = &
      Misc%CUMTRANS ( 1:NLAYERS, 1:N_USER_VZANGLES )

!  solar stuff

           BEAM_CUTOFF ( 1:N_SZANGLES ) = &
      Misc%BEAM_CUTOFF ( 1:N_SZANGLES )
           TRANS_SOLAR_BEAM ( 1:N_SZANGLES ) = &
      Misc%TRANS_SOLAR_BEAM ( 1:N_SZANGLES )
           DO_REFLECTED_DIRECTBEAM ( 1:N_SZANGLES ) = &
      Misc%DO_REFLECTED_DIRECTBEAM ( 1:N_SZANGLES )

!  Average secant

           T_DELT_MUBAR ( 1:NLAYERS, 1:N_SZANGLES ) = &
      Misc%T_DELT_MUBAR ( 1:NLAYERS, 1:N_SZANGLES )
           T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:N_SZANGLES ) = &
      Misc%T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:N_SZANGLES )
           AVERAGE_SECANT ( 1:NLAYERS, 1:N_SZANGLES ) = &
      Misc%AVERAGE_SECANT ( 1:NLAYERS, 1:N_SZANGLES )
           INITIAL_TRANS ( 1:NLAYERS, 1:N_SZANGLES ) = &
      Misc%INITIAL_TRANS ( 1:NLAYERS, 1:N_SZANGLES )

! 1/31/21. Version 2.8.3. additional array ITRANS_USERM to be unpacked

           ITRANS_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES ) = &
      Misc%ITRANS_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES )
           LOCAL_CSZA ( 0:NLAYERS, 1:N_SZANGLES ) = &
      Misc%LOCAL_CSZA ( 0:NLAYERS, 1:N_SZANGLES )

      END SUBROUTINE VLIDORT_UNPACK_MISC

!

      SUBROUTINE VLIDORT_UNPACK_THERM ( Therm,                    & ! Input Sstructure to be unpacked
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, & ! Input Numbers
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Output
        T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )  ! Output

      USE VLIDORT_PARS_m, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_THERMAL_COEFFS
      USE VLIDORT_Work_def_m

!  Input
!  -----

      TYPE(VLIDORT_Work_Thermal), INTENT (IN) :: Therm

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES

!  Output
!  ------

!  From THERMAL_SETUP

      DOUBLE PRECISION, INTENT (OUT) :: THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS )

      DOUBLE PRECISION, INTENT (OUT) :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Unpack structure

            THERMCOEFFS ( 1:NLAYERS, 1:N_THERMAL_COEFFS ) = &
      Therm%THERMCOEFFS ( 1:NLAYERS, 1:N_THERMAL_COEFFS )
            DELTAU_POWER ( 1:NLAYERS, 1:N_THERMAL_COEFFS ) = &
      Therm%DELTAU_POWER ( 1:NLAYERS, 1:N_THERMAL_COEFFS )
            XTAU_POWER ( 1:N_PARTLAYERS, 1:N_THERMAL_COEFFS ) = &
      Therm%XTAU_POWER ( 1:N_PARTLAYERS, 1:N_THERMAL_COEFFS )
            TCOM1 ( 1:NLAYERS, 1:N_THERMAL_COEFFS ) = &
      Therm%TCOM1 ( 1:NLAYERS, 1:N_THERMAL_COEFFS )

            T_DIRECT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS ) = &
      Therm%T_DIRECT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS )
            T_DIRECT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS ) = &
      Therm%T_DIRECT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS )

            T_UT_DIRECT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS ) = &
      Therm%T_UT_DIRECT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS )
            T_UT_DIRECT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS ) = &
      Therm%T_UT_DIRECT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS )

      END SUBROUTINE VLIDORT_UNPACK_THERM

!

      SUBROUTINE VLIDORT_UNPACK_MULT ( Mult,                & ! Input structure to be unpacked
        NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, & ! Input numbers
        SIGMA_M, SIGMA_P, EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )        ! Output

!  1/31/21. Version 2.8.3. additional arrays SIGMA_M, SIGMA_P to be unpacked

      USE VLIDORT_PARS_m, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS
      USE VLIDORT_Work_def_m

!  Input
!  -----

      TYPE(VLIDORT_Work_Multiplier), INTENT (IN) :: Mult

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES

!  Output
!  ------

!  From EMULT_MASTER or EMULT_MASTER_OBSGEO

      DOUBLE PRECISION, INTENT (OUT) :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (OUT) :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Green's function requires additional arrays SIGMA_M, SIGMA_P to be unpacked

      DOUBLE PRECISION, INTENT (OUT) :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Unpack structure
!    -- 1/31/21. Version 2.8.3. arrays SIGMA_M, SIGMA_P to be unpacked

           SIGMA_M ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES ) = &
      Mult%SIGMA_M ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES )
           SIGMA_P ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES ) = &
      Mult%SIGMA_P ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_SZANGLES )

           EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES ) = &
      Mult%EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES )
           EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES ) = &
      Mult%EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:N_SZANGLES )

           UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES ) = &
      Mult%UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES )
           UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES ) = &
      Mult%UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:N_SZANGLES )

      END SUBROUTINE VLIDORT_UNPACK_MULT

      END MODULE vlidort_unpack_m
