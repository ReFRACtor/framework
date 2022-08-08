
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

      MODULE VLIDORT_Work_def_m

!  This module contains the following VLIDORT work structures:
!  (n/a = not available at present). List revised, 7/5/16

!         VLIDORT_Work_RadTran    n/a
!            VLIDORT_Work_User    n/a
!           VLIDORT_Work_Solar    n/a
!         VLIDORT_Work_Thermal
!        VLIDORT_Work_Geometry    n/a
!    VLIDORT_Work_Miscellanous
!      VLIDORT_Work_Multiplier
!     VLIDORT_Work_Corrections    removed for Version 2.8
!      VLIDORT_Work_FirstOrder
!         VLIDORT_Work_Control    n/a
!         VLIDORT_Work_Optical    n/a
!            VLIDORT_Work_BRDF
!          VLIDORT_Work_SLEAVE
!          VLIDORT_Work_Output    n/a

!  1/31/21. Version 2.8.3. 
!    -- Need additional variables SIGMA_M, SIGMA_P, ITRANS_USERM to be packed/unpacked
!    -- Additional inputs  to VLIDORT_PACK_MISC   and VLIDORT_PACK_MULT
!    -- Additional outputs to VLIDORT_UNPACK_MISC and VLIDORT_UNPACK_MULT
!    -- Modification occasioned by use of  Green's function option.

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_GEOMETRIES, &
                                 MAXSTREAMS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_USER_RELAZMS,  &
                                 MAXBEAMS, MAX_THERMAL_COEFFS, MAX_DIRECTIONS, MAXSTOKES_SQ

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Thermal


      DOUBLE PRECISION :: THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )

      DOUBLE PRECISION :: TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS )

      DOUBLE PRECISION :: T_DIRECT_UP    ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DIRECT_DN    ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION :: T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS )


      END TYPE VLIDORT_Work_Thermal

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Miscellanous

!mick mod 9/19/2017 - added LEVELS_SOLARTRANS

      DOUBLE PRECISION :: DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION :: PARTAU_VERT  ( MAX_PARTLAYERS )

      DOUBLE PRECISION :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Version 2.8. these arrays are not packed (7/5/16), comment out
!      DOUBLE PRECISION :: TRUNC_FACTOR ( MAXLAYERS )
!      DOUBLE PRECISION :: FAC1 ( MAXLAYERS )

      DOUBLE PRECISION :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

      INTEGER ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      DOUBLE PRECISION :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: CUMTRANS     ( MAXLAYERS, MAX_USER_STREAMS )

      DOUBLE PRECISION :: LOCAL_CSZA   ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Need additional variable ITRANS_USERM to be packed

      DOUBLE PRECISION :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

      END TYPE VLIDORT_Work_Miscellanous

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_Multiplier


!  Version 2.8. these arrays are not packed (7/5/16), comment out
!      LOGICAL ::          EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )


      DOUBLE PRECISION :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Need additional variables SIGMA_M, SIGMA_P to be packed

      DOUBLE PRECISION :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

      END TYPE VLIDORT_Work_Multiplier

! #####################################################################
! #####################################################################

!  Version 2.8. this type structure no longer required. (7/5/16), comment out
!     --> SSCORR_NADIR ROUTINE HAS BEEN REMOVED.

!      TYPE VLIDORT_Work_Corr_Nadir
!      DOUBLE PRECISION :: ZMAT_UP ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION :: ZMAT_DN ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION :: TMS    ( MAXLAYERS )
!      DOUBLE PRECISION :: SSFDEL ( MAXLAYERS )
!      DOUBLE PRECISION :: SS_CUMSOURCE_UP &( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!      DOUBLE PRECISION :: SS_CUMSOURCE_DN ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!      END TYPE VLIDORT_Work_Corr_Nadir

! #####################################################################
! #####################################################################

!  Version 2.8. this type structure no longer required. (7/5/16), comment out
!     --> SSCORR_OUTGOING ROUTINE HAS BEEN REMOVED.

!      TYPE VLIDORT_Work_Corr_OutGo
!      DOUBLE PRECISION :: UP_MULTIPLIERS ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_MULTIPLIERS ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: UP_MULTIPLIERS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_MULTIPLIERS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: UP_LOSTRANS ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_LOSTRANS ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: UP_LOSTRANS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_LOSTRANS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: BOA_ATTN ( MAX_GEOMETRIES )
!      END TYPE VLIDORT_Work_Corr_OutGo

! #####################################################################
! #####################################################################

!  Version 2.8. this type structure no longer required. (7/5/16), comment out
!     --> DB ROUTINES HAVE BEEN REMOVED. 

!     TYPE VLIDORT_Work_Corr_DB
!      DOUBLE PRECISION :: DB_CUMSOURCE   ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!      DOUBLE PRECISION :: EXACTDB_SOURCE ( MAX_GEOMETRIES, MAXSTOKES )
!      DOUBLE PRECISION :: ATTN_DB_SAVE   ( MAX_GEOMETRIES )
!      END TYPE VLIDORT_Work_Corr_DB

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_FirstOrder


      DOUBLE PRECISION :: FO_STOKES_SS  ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DB  ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: FO_STOKES_DTA ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DTS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: FO_STOKES     ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )


      END TYPE VLIDORT_Work_FirstOrder

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_BRDF


      DOUBLE PRECISION ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      DOUBLE PRECISION ::   BRDF_F_0 ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F   ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

      DOUBLE PRECISION ::   USER_BRDF_F_0 ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_BRDF_F   ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION ::   EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )


      END TYPE VLIDORT_Work_BRDF

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Work_SLEAVE


      DOUBLE PRECISION ::   SLTERM_ISOTROPIC  ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION ::   SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      DOUBLE PRECISION ::   SLTERM_F_0      ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_SLTERM_F_0 ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )


      END TYPE VLIDORT_Work_SLEAVE 

! #####################################################################
! #####################################################################

!  Version 2.8. this type structure no longer required. (7/5/16), comment out

!      TYPE VLIDORT_Work_Corrections
!      TYPE(VLIDORT_Work_Corr_Nadir) :: Nadir
!      TYPE(VLIDORT_Work_Corr_OutGo) :: OutGo
!      TYPE(VLIDORT_Work_Corr_DB)    :: DB
!      END TYPE VLIDORT_Work_Corrections

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE

!  Cleaned up, 7/5/16 for Version 2.8

      PUBLIC :: VLIDORT_Work_Thermal,      &    
                VLIDORT_Work_Miscellanous, &    
                VLIDORT_Work_Multiplier,   &    
                VLIDORT_Work_FirstOrder,   &    
                VLIDORT_Work_BRDF,         &    
                VLIDORT_Work_SLEAVE

!  Following were defined at one point
                !VLIDORT_Work_RadTran,  &
                !VLIDORT_Work_User,     &
                !VLIDORT_Work_Solar,    &    
                !VLIDORT_Work_Geometry, &    
                !VLIDORT_Work_Corr_Nadir, &
                !VLIDORT_Work_Corr_OutGo, &
                !VLIDORT_Work_Corr_DB,    &
                !VLIDORT_Work_Control, &    
                !VLIDORT_Work_Optical, &    
                !VLIDORT_Work_Output,  & 
                !VLIDORT_Work_Corrections

      END MODULE VLIDORT_Work_def_m

