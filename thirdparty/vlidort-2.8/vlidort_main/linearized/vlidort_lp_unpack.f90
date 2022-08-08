
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
! #            VLIDORT_UNPACK_LAP_MISC                          #
! #            VLIDORT_UNPACK_LP_MULT                           #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. No Changes.

      MODULE vlidort_lp_unpack_m

      PRIVATE
      PUBLIC :: VLIDORT_UNPACK_LAP_MISC, &
                VLIDORT_UNPACK_LP_MULT

      CONTAINS

      SUBROUTINE VLIDORT_UNPACK_LAP_MISC ( LAP_Misc,                    & ! Input
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS,                       & ! Input
        N_SZANGLES, N_USER_VZANGLES, NMOMENTS, N_TOTALPROFILE_WFS,      & ! Input
        L_DELTAU_VERT, L_OMEGA_GREEK,                                   & ! Output
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,           & ! Output
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,                 & ! Output
        LP_AVERAGE_SECANT, LP_INITIAL_TRANS,                            & ! Output
        LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,                               & ! Output
        LP_LEVELS_SOLARTRANS, LP_PARTIALS_SOLARTRANS )                    ! Output

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXBEAMS, MAXSTREAMS, &
                                 MAX_ATMOSWFS, MAX_PARTLAYERS, MAX_USER_STREAMS

      USE VLIDORT_LinWork_def_m

!  Input
!  -----

!  structure

      TYPE(VLIDORT_LinWork_Miscellanous), INTENT (IN) :: LAP_Misc

!  Control

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          N_TOTALPROFILE_WFS

!  Output From VLIDORT_LAP_MISCSETUPS
!  ----------------------------------

!  Disabled

      !DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION, INTENT (OUT) :: L_GREEKMAT_TOTAL  ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_SLANT   ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !LOGICAL, INTENT (OUT) ::          DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      !DOUBLE PRECISION, INTENT (OUT) :: L_TRUNC_FACTOR  ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized optical

      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Linearized Transmittances

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Pseudo-spherical average-secant linearization
!mick mod 9/19/2017 - added LP_LEVELS_SOLARTRANS & LP_LEVELS_SOLARTRANS to output

      DOUBLE PRECISION, INTENT (OUT) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Unpack structure

      !         L_OMEGA_TOTAL ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS ) = &
      !LAP_Misc%L_OMEGA_TOTAL ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS )
      !         L_GREEKMAT_TOTAL ( 1:N_TOTALPROFILE_WFS, 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES_SQ ) = &
      !LAP_Misc%L_GREEKMAT_TOTAL ( 1:N_TOTALPROFILE_WFS, 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES_SQ )
      !         L_DELTAU_SLANT ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES ) = &
      !LAP_Misc%L_DELTAU_SLANT ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES )
      !         DO_SCATMAT_VARIATION ( 1:NLAYERS, 1:N_TOTALPROFILE_WFS ) = &
      !LAP_Misc%DO_SCATMAT_VARIATION ( 1:NLAYERS, 1:N_TOTALPROFILE_WFS )
      !         L_TRUNC_FACTOR ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS ) = &
      !LAP_Misc%L_TRUNC_FACTOR ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS )

               L_DELTAU_VERT ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS ) = &
      LAP_Misc%L_DELTAU_VERT ( 1:N_TOTALPROFILE_WFS, 1:NLAYERS )
               L_OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%L_OMEGA_GREEK ( 0:NMOMENTS, 1:NLAYERS, 1:NSTOKES, 1:NSTOKES, 1:N_TOTALPROFILE_WFS )

               L_T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%L_T_DELT_DISORDS ( 1:NSTREAMS, 1:NLAYERS, 1:N_TOTALPROFILE_WFS )
               L_T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%L_T_DISORDS_UTDN ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALPROFILE_WFS )
               L_T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%L_T_DISORDS_UTUP ( 1:NSTREAMS, 1:N_PARTLAYERS, 1:N_TOTALPROFILE_WFS )

               L_T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%L_T_DELT_USERM ( 1:NLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS )
               L_T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%L_T_UTDN_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS )
               L_T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%L_T_UTUP_USERM ( 1:N_PARTLAYERS, 1:N_USER_VZANGLES, 1:N_TOTALPROFILE_WFS )

               LP_AVERAGE_SECANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%LP_AVERAGE_SECANT ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
               LP_INITIAL_TRANS ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%LP_INITIAL_TRANS ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
               LP_T_DELT_MUBAR ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%LP_T_DELT_MUBAR ( 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
               LP_T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%LP_T_UTDN_MUBAR ( 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
               LP_LEVELS_SOLARTRANS ( 0:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%LP_LEVELS_SOLARTRANS ( 0:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
               LP_PARTIALS_SOLARTRANS ( 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LAP_Misc%LP_PARTIALS_SOLARTRANS ( 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )

      END SUBROUTINE VLIDORT_UNPACK_LAP_MISC

!

      SUBROUTINE VLIDORT_UNPACK_LP_MULT ( LP_Mult,               & ! Input
        NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES,      & ! Input
        N_TOTALPROFILE_WFS,                                      & ! Input
        LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN ) ! Output

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAX_ATMOSWFS

      USE VLIDORT_LinWork_def_m

!  Input
!  -----

!  structure

      TYPE(VLIDORT_LinWork_Multiplier), INTENT (IN) :: LP_Mult

!  control

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      INTEGER, INTENT (IN) ::          N_TOTALPROFILE_WFS

!  Outputfrom LP_EMULT_MASTER or LP_EMULT_MASTER_OBSGEO
!  ----------------------------------------------------

      DOUBLE PRECISION, INTENT (OUT) :: LP_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: LP_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Unpack structure

              LP_EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LP_Mult%LP_EMULT_UP ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
              LP_EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LP_Mult%LP_EMULT_DN ( 1:N_USER_VZANGLES, 1:NLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )

              LP_UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LP_Mult%LP_UT_EMULT_UP ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )
              LP_UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS ) = &
      LP_Mult%LP_UT_EMULT_DN ( 1:N_USER_VZANGLES, 1:N_PARTLAYERS, 1:NLAYERS, 1:N_SZANGLES, 1:N_TOTALPROFILE_WFS )

      END SUBROUTINE VLIDORT_UNPACK_LP_MULT

      END MODULE vlidort_lp_unpack_m

