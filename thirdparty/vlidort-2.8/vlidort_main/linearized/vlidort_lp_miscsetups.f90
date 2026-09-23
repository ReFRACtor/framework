
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
! #            VLIDORT_LAP_MISCSETUPS (master)                  #
! #              VLIDORT_LP_PREPTRANS                           #
! #                                                             #
! #            LP_EMULT_MASTER (master), calling:               #
! #              LP_WHOLELAYER_EMULT_UP                         #
! #              LP_WHOLELAYER_EMULT_DN                         #
! #              LP_PARTLAYER_EMULT_UP                          #
! #              LP_PARTLAYER_EMULT_DN                          #
! #                                                             #
! #            LP_EMULT_MASTER_OBSGEO,  calling..               #
! #              LP_WHOLELAYER_EMULT_OG_UP                      #
! #              LP_WHOLELAYER_EMULT_OG_DN                      #
! #              LP_PARTLAYER_EMULT_OG_UP                       #
! #              LP_PARTLAYER_EMULT_OG_DN                       #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. No Changes here.

!  Notes for Version 2.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/10/13  - Introduce TAYLOR_LIMIT parameters in module VLIDORT_PARS, instead of "HOPITAL_TOLERANCE"
!     Rob  Fix 05/10/13  - L'Hopitals Rule replaced by Taylor series (original calculation was first term in series!)
!     Rob  Fix 02/19/14  - Final  Taylor series stuff.

      MODULE vlidort_lp_miscsetups_m

!  @@@ Rob Fix 10 May 13 - need the Taylor series routines

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_L_1

      PRIVATE
      PUBLIC  :: VLIDORT_LAP_MISCSETUPS, &
                 LP_EMULT_MASTER,        &
                 LP_EMULT_MASTER_OBSGEO

      CONTAINS

      SUBROUTINE VLIDORT_LAP_MISCSETUPS ( &
        DO_SOLAR_SOURCES, DO_USER_STREAMS, DO_DELTAM_SCALING,             & ! Input flags
        DO_PLANE_PARALLEL, DO_SOLUTION_SAVING, DO_ATMOS_LINEARIZATION,    & ! Input flags 
        NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, NBEAMS,               & ! Input numbers
        NMOMENTS, NSTOKES_SQ, N_PARTLAYERS, PARTLAYERS_LAYERIDX,          & ! Input numbers/partials
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                          & ! Input Optical props
        LV_FLAG, LV_NUMBER,                                               & ! Input Lin control 
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, DELTAU_SLANT_UNSCALED,    & ! Input derived properties
        PARTAU_SLANT_UNSCALED, LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,    & ! Input derived properties
        OMEGA_GREEK, TRUNC_FACTOR, FAC1,                                  & ! Input derived properties
        MUELLER_INDEX, QUAD_STREAMS, USER_SECANTS,                        & ! Input streams/Muller
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                   & ! Input Dis.Ord. Trans.
        BEAM_CUTOFF, T_DELT_MUBAR, T_UTDN_MUBAR, AVERAGE_SECANT,          & ! Input beam Trans.
        SOLARBEAM_ATRANS, T_DELT_USERM, T_UTDN_USERM,  T_UTUP_USERM,     & ! Input User Trans.
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, & ! Input Linearized Optical props.
        L_DELTAU_VERT, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL, L_DELTAU_SLANT,        & ! Output Linearized scaled properties
        DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, L_OMEGA_GREEK,                  & ! Output Linearized scaled properties
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, & ! Output linearized beam stuff
        LP_LEVELS_SOLARTRANS, LP_PARTIALS_SOLARTRANS,                          & ! Output linearized beam stuff
        LP_SOLARBEAM_BOATRANS, LP_SOLARBEAM_ATRANS,                            & ! Output linearized beam stuff
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,                  & ! Output linearized DisOrd Trans.
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM )                         ! Output lineraized User Trans.

!mick fix 9/19/2017 - to facilitate correction of linearized direct flux,
!                     (1) added DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED,
!                         LEVELS_SOLARTRANS, & PARTIALS_SOLARTRANS as inputs
!                     (2) added LP_LEVELS_SOLARTRANS & LP_PARTIALS_SOLARTRANS as outputs

!  4/9/19. Version 2.8.1, add BOATRANS, ATRANS and linearization thereof.

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXMOMENTS_INPUT, MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, ZERO


      USE VLIDORT_LA_MISCSETUPS_m

      IMPLICIT NONE

!  Input flags
!mick fix 9/19/2017 - added DO_SOLAR_SOURCES to input

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION

!  Input numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          NBEAMS

!  Input numbers and partials

      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTOKES_SQ
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Input linearization control

      LOGICAL, INTENT (IN) ::          LV_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LV_NUMBER ( MAXLAYERS )

!  Input bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Input optical (unscaled)

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Input scaled optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT_UNSCALED ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_SLANT_UNSCALED ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: FAC1         ( MAXLAYERS )

!  Input streams

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS ( MAX_USER_STREAMS )

!  Input discrete Ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  Input beam stuff

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Input User transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SOLARBEAM_ATRANS ( MAXBEAMS )

!  Input linearized optical (unscaled)

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Output linearized scaled opticals

      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_SLANT   ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )

      LOGICAL, INTENT (OUT) ::          DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TRUNC_FACTOR       ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_GREEK        ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Output Beam parameterization linearized

      DOUBLE PRECISION, INTENT (OUT) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: LP_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_SOLARBEAM_BOATRANS  ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_SOLARBEAM_ATRANS    ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Output linearized Discrete Ordinate tranmsittances

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Output linearized User Transmittances

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Local variables

      INTEGER :: IB, K, N, UT, Q
      DOUBLE PRECISION :: HELPTRANS

!  miscellaneous setup operations for linearized quantities

      CALL VLIDORT_LA_DELTAMSCALE ( &
        DO_DELTAM_SCALING, DO_ATMOS_LINEARIZATION,                        & ! Input flags
        NSTOKES, NLAYERS, NBEAMS, NMOMENTS, NSTOKES_SQ,                   & ! Input Numbers
        LV_FLAG, LV_NUMBER,                                               & ! Input linearization control  
        DELTAU_SLANT, TRUNC_FACTOR, FAC1,                                 & ! Input Optical (scaled)
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                          & ! Input Optical
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, & ! Input Optical linearized
        L_DELTAU_VERT, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL,                   & ! Output Optical linearized
        L_DELTAU_SLANT, DO_SCATMAT_VARIATION, L_TRUNC_FACTOR )              ! Output Optical linearized

!  Initialise single scatter albedo variational quantities

      CALL VLIDORT_LA_SSALBINIT ( &
        DO_ATMOS_LINEARIZATION, NSTOKES, NLAYERS, NMOMENTS, & ! Input flag and numbers
        LV_FLAG, LV_NUMBER, MUELLER_INDEX,  & ! Input Lin-control
        OMEGA_GREEK, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL,       & ! Inputs others
        L_OMEGA_GREEK )                                       ! Output

!  Linearization of transmittances

      CALL VLIDORT_LA_PREPTRANS ( &
        DO_SOLUTION_SAVING, DO_USER_STREAMS, NLAYERS, NSTREAMS, & ! Input flags and numbers
        N_USER_STREAMS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,      & ! Input numbers
        DELTAU_VERT, PARTAU_VERT, QUAD_STREAMS, USER_SECANTS,   & ! Input Optical and streams
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,         & ! Input Discr-Ord. Trans.
        T_DELT_USERM,   T_UTDN_USERM,   T_UTUP_USERM,           & ! Input User-strm. Trans.
        LV_FLAG, LV_NUMBER, L_DELTAU_VERT,                      & ! Input linearization
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,   & ! Output linearized Trans (DO)
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM )          ! Output linearized Trans (US)

!  Linearization of pseudo-spherical setup. Additional arguments 4/9/19 (SOLARBEAM_ATRANS)

      CALL VLIDORT_LP_PREPTRANS ( &
        DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, NLAYERS, NBEAMS,      & ! Input flags/numbers
        N_PARTLAYERS, PARTLAYERS_LAYERIDX, LV_NUMBER, DELTAU_VERT, & ! Input partial/linearization control
        PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT, AVERAGE_SECANT,  & ! Input Scaled opticals
        BEAM_CUTOFF, T_DELT_MUBAR, T_UTDN_MUBAR, SOLARBEAM_ATRANS, & ! Input Solar beam
        LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, LP_INITIAL_TRANS,        & ! Output Linearized Solar beam
        LP_AVERAGE_SECANT, LP_SOLARBEAM_ATRANS )                     ! Output Linearized Solar beam

!mick fix 9/19/2017 - added calculation of LP_LEVELS_SOLARTRANS & LP_PARTIALS_SOLARTRANS
!                     to facilitate correction of linearized direct flux (following LIDORT)

!  Whole layers 4/19/19. Add LP_SOLARBEAM_BOATRANS

      LP_LEVELS_SOLARTRANS(0:NLAYERS,1:NLAYERS,1:NBEAMS,:) = ZERO
      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
         IF ( LEVELS_SOLARTRANS(N,IB) .GT. ZERO  ) THEN
           DO K = 1, N
             IF ( LV_FLAG(K) ) THEN
               HELPTRANS = - LEVELS_SOLARTRANS(N,IB) * DELTAU_SLANT_UNSCALED(N,K,IB)
               DO Q = 1, LV_NUMBER(K)
                 LP_LEVELS_SOLARTRANS(N,K,IB,Q) = HELPTRANS * L_DELTAU_VERT_INPUT(Q,K)
                 if ( N==NLAYERS) LP_SOLARBEAM_BOATRANS(IB,K,Q) = LP_LEVELS_SOLARTRANS(N,K,IB,Q)
               ENDDO
             ENDIF
           ENDDO
         ENDIF
        ENDDO
      ENDDO

!  Partial layers

      LP_PARTIALS_SOLARTRANS(:,1:NLAYERS,1:NBEAMS,:) = zero
      IF ( N_PARTLAYERS .gt. 0 ) THEN
        DO IB = 1, NBEAMS
          do UT = 1, N_PARTLAYERS
            N = PARTLAYERS_LAYERIDX(UT)
            IF ( PARTIALS_SOLARTRANS(UT,IB) .gt. zero  ) THEN
              DO K = 1, N
                IF ( LV_FLAG(K) ) THEN
                  HELPTRANS = - PARTIALS_SOLARTRANS(UT,IB) * PARTAU_SLANT_UNSCALED(UT,K,IB)
                  DO Q = 1, LV_NUMBER(K)
                    LP_PARTIALS_SOLARTRANS(UT,K,IB,Q) = HELPTRANS * L_DELTAU_VERT_INPUT(Q,K)
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          enddo
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LAP_MISCSETUPS

!

      SUBROUTINE VLIDORT_LP_PREPTRANS ( &
        DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, NLAYERS, NBEAMS,      & ! Input flags/numbers
        N_PARTLAYERS, PARTLAYERS_LAYERIDX, LV_NUMBER, DELTAU_VERT, & ! Input partial/linearization control
        PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT, AVERAGE_SECANT,  & ! Input Scaled opticals
        BEAM_CUTOFF, T_DELT_MUBAR, T_UTDN_MUBAR, SOLARBEAM_ATRANS, & ! Input Solar beam
        LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, LP_INITIAL_TRANS,        & ! Output Linearized Solar beam
        LP_AVERAGE_SECANT, LP_SOLARBEAM_ATRANS )                     ! Output Linearized Solar beam

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Input flags
!mick fix 9/19/2017 - added DO_SOLAR_SOURCES to input

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  Input numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS

!  Input partial/linearization control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          LV_NUMBER   ( MAXLAYERS )

!  Input opticals (scaled)

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Input Beam parameterization

      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SOLARBEAM_ATRANS ( MAXBEAMS )

!  Output Beam parameterization linearized

      DOUBLE PRECISION, INTENT (OUT) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: LP_SOLARBEAM_ATRANS    ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER ::          N, Q, UT, K, IB
      DOUBLE PRECISION :: WDEL, VAR, RHO, FAC, DELT, LAMDA, LOGD

!mick fix 9/19/2017 - added this return

!  Nothing to do if no solar sources
!  ---------------------------------

      IF ( .NOT.DO_SOLAR_SOURCES ) RETURN

!mick fix 6/29/11 - initialize some outputs

      LP_T_DELT_MUBAR   = ZERO
      LP_T_UTDN_MUBAR   = ZERO
      LP_INITIAL_TRANS  = ZERO
      LP_AVERAGE_SECANT = ZERO

!  linearization of Initial transmittances
!  =======================================

!   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
            DO Q = 1, LV_NUMBER(N)
              LP_INITIAL_TRANS(N,N,IB,Q) = ZERO
            ENDDO
            IF ( N .GT. 1 ) THEN
              DO K = 1, N-1
                DO Q = 1, LV_NUMBER(K)
                  LP_INITIAL_TRANS(N,K,IB,Q) = &
                  - L_DELTAU_VERT(Q,K) * DELTAU_SLANT(N-1,K,IB)
                ENDDO
              ENDDO
            ENDIF
          ELSE
            DO K = 1, N
              DO Q = 1, LV_NUMBER(K)
                LP_INITIAL_TRANS(N,K,IB,Q) = ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LV_NUMBER(N)
                LP_AVERAGE_SECANT(N,N,IB,Q) = ZERO
              ENDDO
            ELSE
              IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
                DO Q = 1, LV_NUMBER(N)
                  LP_AVERAGE_SECANT(N,N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( DELTAU_SLANT(N,K,IB)   - DELTAU_SLANT(N-1,K,IB) ) / DELT
                  DO Q = 1, LV_NUMBER(K)
                    LP_AVERAGE_SECANT(N,K,IB,Q) = L_DELTAU_VERT(Q,K) * FAC
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LV_NUMBER(K)
                    LP_AVERAGE_SECANT(N,K,IB,Q) = ZERO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Linearization of Whole layer Transmittance factors
!  ==================================================

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

          WDEL  = T_DELT_MUBAR(N,IB)
          VAR   = - DELTAU_VERT(N) * WDEL
          LAMDA = AVERAGE_SECANT(N,IB)
          FAC   = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

          IF ( .NOT. DO_PLANE_PARALLEL ) THEN

            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LV_NUMBER(N)
                LP_T_DELT_MUBAR(N,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
              ENDDO
            ELSE
              IF  ( N.LE.BEAM_CUTOFF(IB) ) THEN
                DO Q = 1, LV_NUMBER(N)
                  RHO = LP_AVERAGE_SECANT(N,N,IB,Q)
                  LP_T_DELT_MUBAR(N,N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC + VAR * RHO
                ENDDO
                DO K = 1, N-1
                  DO Q = 1, LV_NUMBER(K)
                    RHO = LP_AVERAGE_SECANT(N,K,IB,Q)
                    LP_T_DELT_MUBAR(N,K,IB,Q) = VAR * RHO
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LV_NUMBER(K)
                    LP_T_DELT_MUBAR(N,K,IB,Q) = ZERO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

!  Plane-parallel

          ELSE IF ( DO_PLANE_PARALLEL ) THEN

            IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
              DO Q = 1, LV_NUMBER(N)
                LP_T_DELT_MUBAR(N,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
              ENDDO
              DO K = 1, N-1
                DO Q = 1, LV_NUMBER(K)
                  LP_T_DELT_MUBAR(N,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                DO Q = 1, LV_NUMBER(K)
                  LP_T_DELT_MUBAR(N,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ENDIF

!  End plane parallel vs. pseudo-spherical

          ENDIF

!  end layer and beam loops

        ENDDO
      ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  =================================================================

      DO IB = 1, NBEAMS
        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
          FAC = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

          IF ( .NOT. DO_PLANE_PARALLEL ) THEN

            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LV_NUMBER(N)
                LP_T_UTDN_MUBAR(UT,N,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
              ENDDO
            ELSE
              IF ( N .LE. BEAM_CUTOFF(IB) ) THEN
                DO Q = 1, LV_NUMBER(N)
                  RHO = LP_AVERAGE_SECANT(N,N,IB,Q)
                  LP_T_UTDN_MUBAR(UT,N,IB,Q) =  L_DELTAU_VERT(Q,N)* FAC + VAR * RHO
                ENDDO
                DO K = 1, N-1
                  DO Q = 1, LV_NUMBER(K)
                    RHO = LP_AVERAGE_SECANT(N,K,IB,Q)
                    LP_T_UTDN_MUBAR(UT,K,IB,Q) = VAR * RHO
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LV_NUMBER(K)
                    LP_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

!  Plane-parallel

          ELSE IF ( DO_PLANE_PARALLEL ) THEN

            IF ( N .LE. BEAM_CUTOFF(IB) ) THEN
              DO Q = 1, LV_NUMBER(N)
                LP_T_UTDN_MUBAR(UT,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
              ENDDO
              DO K = 1, N-1
                DO Q = 1, LV_NUMBER(K)
                  LP_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                DO Q = 1, LV_NUMBER(K)
                  LP_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ENDIF

          ENDIF

!  End optical depth and beam loops

        ENDDO
      ENDDO
     
!  4/29/19 Linearization of ATRANS.

      DO IB = 1, NBEAMS
         DO K = 1, NLAYERS
            DO Q = 1, LV_NUMBER(K)
               LOGD = - L_DELTAU_VERT(Q,K) * DELTAU_SLANT(NLAYERS,K,IB)
               LP_SOLARBEAM_ATRANS(IB,K,Q) = SOLARBEAM_ATRANS(IB) * LOGD
            ENDDO
         ENDDO
      ENDDO
   
!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LP_PREPTRANS

!

      SUBROUTINE LP_EMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, TAYLOR_ORDER,       & ! Input flags + Taylor
        NLAYERS, NBEAMS, N_USER_STREAMS, N_PARTLAYERS, LV_FLAG, LV_NUMBER, & ! Input numbers
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,       & ! Input Level control
        USER_SECANTS, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,             & ! Input optical, streams   
        BEAM_CUTOFF, T_DELT_MUBAR, EMULT_HOPRULE, SIGMA_M, SIGMA_P,        & ! Input solar beam + Multipliers
        T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, ITRANS_USERM,            & ! Input User-stream Trans.
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                      & ! Input Multipliers
        LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,              & ! Input linearized Beam stuff
        LP_T_UTDN_MUBAR, L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,   & ! Input linearized transmittanaces
        LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN )           ! Output linearized Mulitpliers

!  Linearized multipliers for the Beam source terms

!  Version 2p7 notes
!     Rob Fix 05/10/13  - Introduce Taylor series routines (first version) HOPITAL_TOLERANCE replaced by TAYLOR_SMALL.
!     Rob Fix 02/19/14  - Small numbers analysis finalized using Taylor_Order parameter, as for LIDORT

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                 MAX_PARTLAYERS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  Inputs
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

!  Numbers

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      LOGICAL, INTENT (IN) ::           LV_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LV_NUMBER ( MAXLAYERS )


!  Level control

      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

!  Streams, optical

      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Solar beam

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  user-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Multiplier Diagnostics

      LOGICAL, INTENT (IN) ::           EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Linearized Beam-parameterization

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Linearized User transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM  ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  OUTPUTS
!  =======

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (OUT) :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER :: N, UT, K, KPARS

!mick fix 7/23/2014 - initialized for packing
      LP_EMULT_UP = ZERO
      LP_EMULT_DN = ZERO
      LP_UT_EMULT_UP = ZERO
      LP_UT_EMULT_DN = ZERO

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!  Whole layer upwelling
!  ---------------------

!  Loop over all  model  layers N

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           DO K = 1, NLAYERS
             IF ( N.GE.K ) THEN
               IF ( LV_FLAG(K) ) THEN
                 KPARS = LV_NUMBER(K)
                 CALL LP_WHOLELAYER_EMULT_UP ( &
                   DO_PLANE_PARALLEL, N, NBEAMS, N_USER_STREAMS, K, KPARS, & ! Input flag/numbers
                   BEAM_CUTOFF, T_DELT_MUBAR, EMULT_UP, SIGMA_P,              & ! Input Beam stuff + Multipliers
                   T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,                & ! Input User trans.
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,      & ! Input linearized Beam stuff
                   LP_EMULT_UP )                                                ! Output Multiplier
               ENDIF
             ENDIF
           ENDDO
         ENDIF
       ENDDO

!  Partial layer upwelling
!  -----------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           DO K = 1, NLAYERS
             IF ( N.GE.K ) THEN
               IF ( LV_FLAG(K) ) THEN
                 KPARS = LV_NUMBER(K)
                 CALL LP_PARTLAYER_EMULT_UP ( &
                   DO_PLANE_PARALLEL, N, UT, NBEAMS, N_USER_STREAMS, K, KPARS,            & ! Input flag/numbers
                   BEAM_CUTOFF, T_DELT_MUBAR, UT_EMULT_UP, SIGMA_P,                       & ! Input Beam stuff, multipliers
                   T_UTUP_USERM, ITRANS_USERM, L_T_UTUP_USERM,                            & ! Input Multipliers + linearization
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, & ! Input linearized Beam stuff
                   LP_UT_EMULT_UP )                                                         ! Output Multiplier
               ENDIF
             ENDIF
           ENDDO
         ENDIF
       ENDDO

!  end upwelling

      ENDIF

!  Downwelling
!  ===========

      IF ( DO_DNWELLING ) THEN

!  Whole layer downwelling
!  -----------------------

!  Start loop over all  model  layers N

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           DO K = 1, NLAYERS
             IF ( N.GE.K ) THEN
               IF ( LV_FLAG(K) ) THEN
                 KPARS = LV_NUMBER(K)
                 CALL LP_WHOLELAYER_EMULT_DN ( &
                   DO_PLANE_PARALLEL, N, NBEAMS, N_USER_STREAMS, K, KPARS, TAYLOR_ORDER,    & ! Input flag/numbers
                   USER_SECANTS, BEAM_CUTOFF, EMULT_DN, EMULT_HOPRULE, SIGMA_M,             & ! Input Beam stuff + Multipliers
                   DELTAU_VERT, L_DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,  & ! Input User trans.
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,                    & ! Input linearized Beam stuff
                   LP_EMULT_DN )                                                              ! Output Multiplier
               ENDIF
             ENDIF
           ENDDO
         ENDIF
       ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           DO K = 1, NLAYERS
             IF ( N.GE.K ) THEN
               IF ( LV_FLAG(K) ) THEN
                 KPARS = LV_NUMBER(K)
                 CALL LP_PARTLAYER_EMULT_DN ( &
                   DO_PLANE_PARALLEL, N, UT, NBEAMS, N_USER_STREAMS, K, KPARS, TAYLOR_ORDER, & ! Input flag/numbers
                   USER_SECANTS, BEAM_CUTOFF, UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M,           & ! Input Beam stuff, multipliers
                   PARTAU_VERT, L_DELTAU_VERT, T_UTDN_USERM, ITRANS_USERM, L_T_UTDN_USERM,   & ! Input User transmittances
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_UTDN_MUBAR,                     & ! Input linearized Beam stuff
                   LP_UT_EMULT_DN )                                                            ! Output Multiplier
               ENDIF
             ENDIF
           ENDDO
         ENDIF
       ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_EMULT_MASTER

!

      SUBROUTINE LP_WHOLELAYER_EMULT_UP ( &
                   DO_PLANE_PARALLEL, N, NBEAMS, N_USER_STREAMS, K, KPARS, & ! Input flag/numbers
                   BEAM_CUTOFF, T_DELT_MUBAR, EMULT_UP, SIGMA_P,              & ! Input Beam stuff + Multipliers
                   T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,                & ! Input User trans.
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,      & ! Input linearized Beam stuff
                   LP_EMULT_UP )                                                ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          K, KPARS

!  Input beam parameterization and Multipliers

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP     ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Input user transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM   ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Input linearized Beam stuff

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output linearized multiplier

      DOUBLE PRECISION, INTENT (INOUT) :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SU, V1, V2, WDEL, UDEL
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, KPARS
             LP_EMULT_UP(UM,N,K,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              UDEL = T_DELT_USERM(N,UM)
              SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, KPARS
                V1 = -LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB)
                V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                     UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
             ENDDO
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              UDEL = T_DELT_USERM(N,UM)
              SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, KPARS
                V1 = LP_INITIAL_TRANS (N,K,IB,Q) - &
                   ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB) )
                V2 =  UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
              ENDDO
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              UDEL = T_DELT_USERM(N,UM)
              SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, KPARS
                V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                     UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_UP(UM,N,K,IB,Q) = SU * V2
              ENDDO
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, KPARS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q)
                LP_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_WHOLELAYER_EMULT_UP

!

      SUBROUTINE LP_WHOLELAYER_EMULT_DN ( &
                   DO_PLANE_PARALLEL, N, NBEAMS, N_USER_STREAMS, K, KPARS, TAYLOR_ORDER,    & ! Input flag/numbers
                   USER_SECANTS, BEAM_CUTOFF, EMULT_DN, EMULT_HOPRULE, SIGMA_M,             & ! Input Beam stuff + Multipliers
                   DELTAU_VERT, L_DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,  & ! Input User trans.
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,                    & ! Input linearized Beam stuff
                   LP_EMULT_DN )                                                              ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          K, KPARS

!  Input Taylor order and streams

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )

!  Input beam parameterization and Multipliers

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF   ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN      ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M       ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Input optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Input user transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM   ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Input linearized Beam stuff

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output linearized multiplier

      DOUBLE PRECISION, INTENT (INOUT) :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, MULT, TMEW, DELTA, L_DELTA, UDEL, L_LAM, L_MULT
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, KPARS
             LP_EMULT_DN(UM,N,K,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a). K = N. Note the use of L'Hopital's Rule flag.

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              SM   = USER_SECANTS(UM) ; MULT = EMULT_DN(UM,N,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,UM,IB) ; DELTA = DELTAU_VERT(N)
                DO Q = 1, KPARS
                  L_LAM   = LP_AVERAGE_SECANT(N,K,IB,Q)
                  L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                  CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
                  LP_EMULT_DN(UM,N,K,IB,Q) = TMEW * L_MULT
                ENDDO
              ELSE
                SD = TMEW / SIGMA_M(N,UM,IB)
                DO Q = 1, KPARS
                  V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB)
                  V2 = L_T_DELT_USERM(N,UM,Q) - LP_T_DELT_MUBAR(N,K,IB,Q)
                  LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Case (b). N > K

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              SM   = USER_SECANTS(UM) ; MULT = EMULT_DN(UM,N,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,UM,IB) ; DELTA = DELTAU_VERT(N)
                DO Q = 1, KPARS
                  L_LAM   = LP_AVERAGE_SECANT(N,K,IB,Q)
                  CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, ZERO, L_LAM, ZERO, UDEL, SM, L_mult )
                  LP_EMULT_DN(UM,N,K,IB,Q) = LP_INITIAL_TRANS(N,K,IB,Q) * MULT + TMEW * L_MULT
                ENDDO
              ELSE
                SD = TMEW / SIGMA_M(N,UM,IB)
                DO Q = 1, KPARS
                  V1 =   LP_INITIAL_TRANS(N,K,IB,Q) - &
                       ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB) )
                  V2 = - LP_T_DELT_MUBAR(N,K,IB,Q)
                  LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a) . K = N

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              SM   = USER_SECANTS(UM) ; MULT = EMULT_DN(UM,N,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,UM,IB) ; DELTA = DELTAU_VERT(N)
                DO Q = 1, KPARS
                  L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                  CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_mult )
                  LP_EMULT_DN(UM,N,K,IB,Q) = TMEW * L_MULT
                ENDDO
              ELSE
                SD = TMEW / SIGMA_M(N,UM,IB)
                DO Q = 1, KPARS
                  V2 = L_T_DELT_USERM(N,UM,Q) - LP_T_DELT_MUBAR(N,K,IB,Q)
                  LP_EMULT_DN(UM,N,K,IB,Q) = SD*V2
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Case (b). N > K

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, KPARS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q)
                LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plaen-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_WHOLELAYER_EMULT_DN

!

      SUBROUTINE LP_PARTLAYER_EMULT_UP ( &
                   DO_PLANE_PARALLEL, N, UT, NBEAMS, N_USER_STREAMS, K, KPARS,            & ! Input flag/numbers
                   BEAM_CUTOFF, T_DELT_MUBAR, UT_EMULT_UP, SIGMA_P,                       & ! Input Beam stuff, multipliers
                   T_UTUP_USERM, ITRANS_USERM, L_T_UTUP_USERM,                            & ! Input Multipliers + linearization
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, & ! Input linearized Beam stuff
                   LP_UT_EMULT_UP )                                                         ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           K, KPARS

!  Solar beam parameterization and multipliers

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF   ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_P       ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  user-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized Beam-parameterization

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Output multipliers

      DOUBLE PRECISION, INTENT (INOUT) :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SU, V1, V2, WDEL, UX_UP
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, KPARS
             LP_UT_EMULT_UP(UM,UT,K,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, KPARS
                V1 = -LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB)
                V2 =         LP_T_UTDN_MUBAR(UT,K,IB,Q) - &
                     UX_UP * LP_T_DELT_MUBAR(N, K,IB,Q) - &
                     WDEL  * L_T_UTUP_USERM(UT,UM,Q)
                LP_UT_EMULT_UP(UM,UT,K,IB,Q) = SU * V2 + &
                                               UT_EMULT_UP(UM,UT,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  ..(b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, KPARS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q) - &
                     ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB) )
                V2 =           LP_T_UTDN_MUBAR(UT,K,IB,Q) - &
                       UX_UP * LP_T_DELT_MUBAR( N,K,IB,Q)
                LP_UT_EMULT_UP(UM,UT,K,IB,Q) = SU * V2 + &
                                               UT_EMULT_UP(UM,UT,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, KPARS
                V2 =         LP_T_UTDN_MUBAR(UT,K,IB,Q) - &
                     UX_UP * LP_T_DELT_MUBAR( N,K,IB,Q) - &
                     WDEL  * L_T_UTUP_USERM(UT,UM,Q)
                LP_UT_EMULT_UP(UM,UT,K,IB,Q) = SU * V2
              ENDDO
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, KPARS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q)
                LP_UT_EMULT_UP(UM,UT,K,IB,Q) = UT_EMULT_UP(UM,UT,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_PARTLAYER_EMULT_UP

!

      SUBROUTINE LP_PARTLAYER_EMULT_DN ( &
                   DO_PLANE_PARALLEL, N, UT, NBEAMS, N_USER_STREAMS, K, KPARS, TAYLOR_ORDER, & ! Input flag/numbers
                   USER_SECANTS, BEAM_CUTOFF, UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M,           & ! Input Beam stuff, multipliers
                   PARTAU_VERT, L_DELTAU_VERT, T_UTDN_USERM, ITRANS_USERM, L_T_UTDN_USERM,   & ! Input User transmittances
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_UTDN_MUBAR,                     & ! Input linearized Beam stuff
                   LP_UT_EMULT_DN )                                                            ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           K, KPARS

!  Input Taylor order, User secants, beam cutoff

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::           BEAM_CUTOFF    ( MAXBEAMS )

!  Multiplier Diagnostics

      LOGICAL, INTENT (IN) ::           EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_M       ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  optical

      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  user-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM    ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ITRANS_USERM    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_UTDN_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )


!  Linearized Beam-parameterization

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Output Multiplier
!mick fix 9/19/2017 - changed intent from "out" to "inout"

      DOUBLE PRECISION, INTENT (INOUT) :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )


!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, MULT, TMEW, DELTA, L_DELTA, UXDN, L_LAM, L_MULT
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, KPARS
             LP_UT_EMULT_DN(UM,UT,K,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a). K = N

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              SM    = USER_SECANTS(UM) ; MULT = UT_EMULT_DN(UM,UT,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,UM,IB) ; DELTA = PARTAU_VERT(UT)
                DO Q = 1, KPARS
                  L_LAM   = LP_AVERAGE_SECANT(N,K,IB,Q)
                  L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                  CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UXDN, SM, L_mult )
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = TMEW * L_MULT
                ENDDO
              ELSE
                SD = TMEW / SIGMA_M(N,UM,IB)
                DO Q = 1, KPARS
                  V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB)
                  V2 = L_T_UTDN_USERM(UT,UM,Q) - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = SD * V2 + UT_EMULT_DN(UM,UT,IB) * V1
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Case (b). K > N

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              SM    = USER_SECANTS(UM) ; MULT = UT_EMULT_DN(UM,UT,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,UM,IB) ; DELTA = PARTAU_VERT(UT)
                DO Q = 1, KPARS
                  L_LAM   = LP_AVERAGE_SECANT(N,K,IB,Q)
                  CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, ZERO, L_LAM, ZERO, UXDN, SM, L_mult )
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = LP_INITIAL_TRANS(N,K,IB,Q) * MULT + TMEW * L_MULT
                ENDDO
              ELSE
                SD = TMEW / SIGMA_M(N,UM,IB)
                DO Q = 1, KPARS
                  V1 = LP_INITIAL_TRANS(N,K,IB,Q) - &
                     ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB) )
                  V2 = - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = SD * V2 + UT_EMULT_DN(UM,UT,IB) * V1
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              SM    = USER_SECANTS(UM) ; MULT = UT_EMULT_DN(UM,UT,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,UM,IB) ; DELTA = PARTAU_VERT(UT)
                DO Q = 1, KPARS
                  L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                  CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UXDN, SM, L_mult )
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = TMEW * L_MULT
                ENDDO
              ELSE
                SD = TMEW / SIGMA_M(N,UM,IB)
                DO Q = 1, KPARS
                  V2 = L_T_UTDN_USERM(UT,UM,Q) - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = SD * V2
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!   Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, KPARS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q)
                LP_UT_EMULT_DN(UM,UT,K,IB,Q) = UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_PARTLAYER_EMULT_DN

!

      SUBROUTINE LP_EMULT_MASTER_OBSGEO ( &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, TAYLOR_ORDER,       & ! Input flags + Taylor
        NLAYERS, NBEAMS, N_PARTLAYERS, LV_FLAG, LV_NUMBER,                 & ! Input numbers
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,       & ! Input Level control
        USER_SECANTS, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,             & ! Input optical, streams   
        BEAM_CUTOFF, T_DELT_MUBAR, EMULT_HOPRULE, SIGMA_M, SIGMA_P,        & ! Input solar beam + Multipliers
        T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, ITRANS_USERM,            & ! Input User-stream Trans.
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                      & ! Input Multipliers
        LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,              & ! Input linearized Beam stuff
        LP_T_UTDN_MUBAR, L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,   & ! Input linearized transmittanaces
        LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN )           ! Output linearized Mulitpliers

!  Linearized multipliers for the Beam source terms

!  Version 2p7 notes
!     Rob Fix 05/10/13  - Introduce Taylor series routines (first version) HOPITAL_TOLERANCE replaced by TAYLOR_SMALL.
!     Rob Fix 02/19/14  - Small numbers analysis finalized using Taylor_Order parameter, as for LIDORT

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                 MAX_PARTLAYERS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  Inputs
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

!  Numbers

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      LOGICAL, INTENT (IN) ::           LV_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LV_NUMBER ( MAXLAYERS )


!  Level control

      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

!  Streams, optical

      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Solar beam

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  user-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Multiplier Diagnostics

      LOGICAL, INTENT (IN) ::           EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Linearized Beam-parameterization

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Linearized User transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM  ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  OUTPUTS
!  =======

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (OUT) :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER :: N, UT, K, KPARS

!mick fix 7/23/2014 - initialized for packing
      LP_EMULT_UP = ZERO
      LP_EMULT_DN = ZERO
      LP_UT_EMULT_UP = ZERO
      LP_UT_EMULT_DN = ZERO

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!  Whole layer upwelling
!  ---------------------

!  Loop over all  model  layers N

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           DO K = 1, NLAYERS
             IF ( N.GE.K ) THEN
               IF ( LV_FLAG(K) ) THEN
                 KPARS = LV_NUMBER(K)
                 CALL LP_WHOLELAYER_EMULT_OG_UP ( &
                   DO_PLANE_PARALLEL, N, NBEAMS, K, KPARS,                    & ! Input flag/numbers
                   BEAM_CUTOFF, T_DELT_MUBAR, EMULT_UP, SIGMA_P,              & ! Input Beam stuff + Multipliers
                   T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,                & ! Input User trans.
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,      & ! Input linearized Beam stuff
                   LP_EMULT_UP )                                                ! Output Multiplier
               ENDIF
             ENDIF
           ENDDO
         ENDIF
       ENDDO

!  Partial layer upwelling
!  -----------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           DO K = 1, NLAYERS
             IF ( N.GE.K ) THEN
               IF ( LV_FLAG(K) ) THEN
                 KPARS = LV_NUMBER(K)
                 CALL LP_PARTLAYER_EMULT_OG_UP ( &
                   DO_PLANE_PARALLEL, N, UT, NBEAMS, K, KPARS,                           & ! Input flag/numbers
                   BEAM_CUTOFF, T_DELT_MUBAR, UT_EMULT_UP, SIGMA_P,                       & ! Input Beam stuff, multipliers
                   T_UTUP_USERM, ITRANS_USERM, L_T_UTUP_USERM,                            & ! Input Multipliers + linearization
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, & ! Input linearized Beam stuff
                   LP_UT_EMULT_UP )                                                         ! Output Multiplier
               ENDIF
             ENDIF
           ENDDO
         ENDIF
       ENDDO

!  end upwelling

      ENDIF

!  Downwelling
!  ===========

      IF ( DO_DNWELLING ) THEN

!  Whole layer downwelling
!  -----------------------

!  Start loop over all  model  layers N

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           DO K = 1, NLAYERS
             IF ( N.GE.K ) THEN
               IF ( LV_FLAG(K) ) THEN
                 KPARS = LV_NUMBER(K)
                 CALL LP_WHOLELAYER_EMULT_OG_DN ( &
                   DO_PLANE_PARALLEL, N, NBEAMS, K, KPARS, TAYLOR_ORDER,                    & ! Input flag/numbers
                   USER_SECANTS, BEAM_CUTOFF, EMULT_DN, EMULT_HOPRULE, SIGMA_M,             & ! Input Beam stuff + Multipliers
                   DELTAU_VERT, L_DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,  & ! Input User trans.
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,                    & ! Input linearized Beam stuff
                   LP_EMULT_DN )                                                              ! Output Multiplier
               ENDIF
             ENDIF
           ENDDO
         ENDIF
       ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           DO K = 1, NLAYERS
             IF ( N.GE.K ) THEN
               IF ( LV_FLAG(K) ) THEN
                 KPARS = LV_NUMBER(K)
                 CALL LP_PARTLAYER_EMULT_OG_DN ( &
                   DO_PLANE_PARALLEL, N, UT, NBEAMS, K, KPARS, TAYLOR_ORDER,                 & ! Input flag/numbers
                   USER_SECANTS, BEAM_CUTOFF, UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M,           & ! Input Beam stuff, multipliers
                   PARTAU_VERT, L_DELTAU_VERT, T_UTDN_USERM, ITRANS_USERM, L_T_UTDN_USERM,   & ! Input User transmittances
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_UTDN_MUBAR,                     & ! Input linearized Beam stuff
                   LP_UT_EMULT_DN )                                                            ! Output Multiplier
               ENDIF
             ENDIF
           ENDDO
         ENDIF
       ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_EMULT_MASTER_OBSGEO

!

      SUBROUTINE LP_WHOLELAYER_EMULT_OG_UP ( &
         DO_PLANE_PARALLEL, N, NBEAMS, K, KPARS,               & ! Input flag/numbers
         BEAM_CUTOFF, T_DELT_MUBAR, EMULT_UP, SIGMA_P,         & ! Input Beam stuff + Multipliers
         T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,           & ! Input User trans.
         LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, & ! Input linearized Beam stuff
         LP_EMULT_UP )                                           ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          K, KPARS

!  Input beam parameterization and Multipliers

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP     ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Input user transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM   ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Input linearized Beam stuff

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output linearized multiplier

      DOUBLE PRECISION, INTENT (INOUT) :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SU, V1, V2, WDEL, UDEL
      INTEGER          :: Q, IB, LUM

!  Start Beam loop
!  ===============

      LUM = 1

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO Q = 1, KPARS
           LP_EMULT_UP(LUM,N,K,IB,Q) = ZERO
         ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

          IF ( K.EQ.N ) THEN
            UDEL = T_DELT_USERM(N,IB)
            SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, KPARS
              V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,LUM,IB)
              V2 = WDEL * L_T_DELT_USERM(N,IB,Q) + &
                   UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
              LP_EMULT_UP(LUM,N,K,IB,Q) = EMULT_UP(LUM,N,IB)*V1 + SU*V2
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            UDEL = T_DELT_USERM(N,IB)
            SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, KPARS
              V1 = LP_INITIAL_TRANS (N,K,IB,Q) - &
                 ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,LUM,IB) )
              V2 =  UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
              LP_EMULT_UP(LUM,N,K,IB,Q) = EMULT_UP(LUM,N,IB)*V1 + SU*V2
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            UDEL = T_DELT_USERM(N,IB)
            SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, KPARS
              V2 = WDEL * L_T_DELT_USERM(N,IB,Q) + UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
              LP_EMULT_UP(LUM,N,K,IB,Q) =  SU * V2
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO Q = 1, KPARS
              V1 = LP_INITIAL_TRANS(N,K,IB,Q)
              LP_EMULT_UP(LUM,N,K,IB,Q) = EMULT_UP(LUM,N,IB) * V1
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_WHOLELAYER_EMULT_OG_UP

!

      SUBROUTINE LP_WHOLELAYER_EMULT_OG_DN ( &
         DO_PLANE_PARALLEL, N, NBEAMS, K, KPARS, TAYLOR_ORDER,                   & ! Input flag/numbers
         USER_SECANTS, BEAM_CUTOFF, EMULT_DN, EMULT_HOPRULE, SIGMA_M,            & ! Input Beam stuff + Multipliers
         DELTAU_VERT, L_DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM, & ! Input User trans.
         LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR,                   & ! Input linearized Beam stuff
         LP_EMULT_DN )                                                             ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL


      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          K, KPARS

!  Input Taylor order and streams

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )

!  Input beam parameterization and Multipliers

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF   ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN      ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M       ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Input optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Input user transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM   ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Input linearized Beam stuff

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output linearized multiplier

      DOUBLE PRECISION, INTENT (INOUT) :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, MULT, TMEW, DELTA, L_DELTA, UDEL, L_LAM, L_MULT
      INTEGER          :: Q, IB, LUM

!  Start Beam loop
!  ===============

      LUM = 1

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO Q = 1, KPARS
           LP_EMULT_DN(LUM,N,K,IB,Q) = ZERO
         ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a). K = N. Note the use of L'Hopital's Rule flag.

          IF ( K.EQ.N ) THEN
            SM    = USER_SECANTS(IB) ; MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UDEL = T_DELT_USERM(N,IB) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
              DO Q = 1, KPARS
                L_LAM   = LP_AVERAGE_SECANT(N,K,IB,Q)
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
                LP_EMULT_DN(LUM,N,K,IB,Q) = TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, KPARS
                V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,LUM,IB)
                V2 = L_T_DELT_USERM(N,IB,Q)-LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_DN(LUM,N,K,IB,Q) = EMULT_DN(LUM,N,IB)*V1+SD*V2
              ENDDO
            ENDIF
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            SM    = USER_SECANTS(IB) ; MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UDEL = T_DELT_USERM(N,IB) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
              DO Q = 1, KPARS
                L_LAM   = LP_AVERAGE_SECANT(N,K,IB,Q)
                CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, ZERO, L_LAM, ZERO, UDEL, SM, L_mult )
                LP_EMULT_DN(LUM,N,K,IB,Q) = TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, KPARS
                V1 =   LP_INITIAL_TRANS(N,K,IB,Q) - &
                     ( LP_AVERAGE_SECANT(N,K,IB,Q)/SIGMA_M(N,LUM,IB) )
                V2 = - LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_DN(LUM,N,K,IB,Q) = EMULT_DN(LUM,N,IB)*V1+SD*V2
              ENDDO
            ENDIF
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            SM    = USER_SECANTS(IB) ; MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UDEL = T_DELT_USERM(N,IB) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
              DO Q = 1, KPARS
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_mult )
                LP_EMULT_DN(LUM,N,K,IB,Q) = TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, KPARS
                V2 = L_T_DELT_USERM(N,IB,Q) - LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_DN(LUM,N,K,IB,Q) = SD*V2
              ENDDO
            ENDIF
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO Q = 1, KPARS
              V1 = LP_INITIAL_TRANS(N,K,IB,Q)
              LP_EMULT_DN(LUM,N,K,IB,Q) = EMULT_DN(LUM,N,IB) * V1
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plaen-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_WHOLELAYER_EMULT_OG_DN

!

      SUBROUTINE LP_PARTLAYER_EMULT_OG_UP ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, K, KPARS,                            & ! Input flag/numbers
             BEAM_CUTOFF, T_DELT_MUBAR, UT_EMULT_UP, SIGMA_P,                       & ! Input Beam stuff, multipliers
             T_UTUP_USERM, ITRANS_USERM, L_T_UTUP_USERM,                            & ! Input Multipliers + linearization
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, & ! Input linearized Beam stuff
             LP_UT_EMULT_UP )                                                         ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           K, KPARS

!  Solar beam parameterization and multipliers

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF   ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_P       ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  user-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized Beam-parameterization

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Output multipliers

      DOUBLE PRECISION, INTENT (INOUT) :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SU, V1, V2, WDEL, UX_UP
      INTEGER          :: Q, IB, LUM

!  Start Beam loop
!  ===============

      LUM = 1

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO Q = 1, KPARS
           LP_UT_EMULT_UP(LUM,UT,K,IB,Q) = ZERO
         ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case (a)

          IF ( K.EQ.N ) THEN
            UX_UP = T_UTUP_USERM(UT,IB)
            SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, KPARS
              V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,LUM,IB)
              V2 =   LP_T_UTDN_MUBAR(UT,K,IB,Q) - &
                 UX_UP * LP_T_DELT_MUBAR(N,K,IB,Q) - WDEL * L_T_UTUP_USERM(UT,IB,Q)
              LP_UT_EMULT_UP(LUM,UT,K,IB,Q) = SU * V2 + UT_EMULT_UP(LUM,UT,IB) * V1
            ENDDO
          ENDIF

!  ..(b)

          IF ( N.GT.K ) THEN
            UX_UP = T_UTUP_USERM(UT,IB)
            SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, KPARS
              V1 = LP_INITIAL_TRANS(N,K,IB,Q) - &
            ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,LUM,IB) )
              V2 = LP_T_UTDN_MUBAR(UT,K,IB,Q) - UX_UP * LP_T_DELT_MUBAR(N,K,IB,Q)
              LP_UT_EMULT_UP(LUM,UT,K,IB,Q) = SU * V2 + UT_EMULT_UP(LUM,UT,IB) * V1
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            UX_UP = T_UTUP_USERM(UT,IB)
            SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, KPARS
              V2 = LP_T_UTDN_MUBAR(UT,K,IB,Q) - &
                UX_UP * LP_T_DELT_MUBAR(N,K,IB,Q) - WDEL * L_T_UTUP_USERM(UT,IB,Q)
              LP_UT_EMULT_UP(LUM,UT,K,IB,Q) = SU * V2
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO Q = 1, KPARS
              V1 = LP_INITIAL_TRANS(N,K,IB,Q)
              LP_UT_EMULT_UP(LUM,UT,K,IB,Q) = UT_EMULT_UP(LUM,UT,IB)*V1
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_PARTLAYER_EMULT_OG_UP

!

      SUBROUTINE LP_PARTLAYER_EMULT_OG_DN ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, K, KPARS, TAYLOR_ORDER,               & ! Input flag/numbers
             USER_SECANTS, BEAM_CUTOFF, UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M,         & ! Input Beam stuff, multipliers
             PARTAU_VERT, L_DELTAU_VERT, T_UTDN_USERM, ITRANS_USERM, L_T_UTDN_USERM, & ! Input User transmittances
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_UTDN_MUBAR,                   & ! Input linearized Beam stuff
             LP_UT_EMULT_DN )                                                          ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                 MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           K, KPARS

!  Input Taylor order, user streams, beam cutoff

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::           BEAM_CUTOFF   ( MAXBEAMS )

!  Multiplier and diagnostics

      LOGICAL, INTENT (IN) ::           EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SIGMA_M       ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  optical

      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  user-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM     ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ITRANS_USERM     ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_UTDN_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS)

!  Linearized Beam-parameterization

      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Output Multiplier
!mick fix 9/19/2017 - changed intent from "out" to "inout"

      DOUBLE PRECISION, INTENT (INOUT) :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, MULT, TMEW, DELTA, L_DELTA, UXDN, L_LAM, L_MULT
      INTEGER          :: Q, IB, LUM

!  Start Beam loop
!  ===============

      LUM = 1

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO Q = 1, KPARS
           LP_UT_EMULT_DN(LUM,UT,K,IB,Q) = ZERO
         ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a), N = K

          IF ( K.EQ.N ) THEN
            SM  = USER_SECANTS(IB) ;  MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,IB) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, KPARS
                L_LAM   = LP_AVERAGE_SECANT(N,K,IB,Q)
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UXDN, SM, L_mult )
                LP_UT_EMULT_DN(LUM,UT,K,IB,Q) = TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, KPARS
                V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,LUM,IB)
                V2 = L_T_UTDN_USERM(UT,IB,Q) - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                LP_UT_EMULT_DN(LUM,UT,K,IB,Q) = SD * V2 + UT_EMULT_DN(LUM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDIF

!  Case (b), N > K. profile only

          IF ( N.GT.K ) THEN
            SM  = USER_SECANTS(IB) ;  MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,IB) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, KPARS
                L_LAM   = LP_AVERAGE_SECANT(N,K,IB,Q)
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, ZERO, L_LAM, ZERO, UXDN, SM, L_mult )
                LP_UT_EMULT_DN(LUM,UT,K,IB,Q) = LP_INITIAL_TRANS(N,K,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, KPARS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q) - &
                   ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,LUM,IB) )
                V2 = - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                LP_UT_EMULT_DN(LUM,UT,K,IB,Q) = SD * V2 + UT_EMULT_DN(LUM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a), N = K.

          IF ( K.EQ.N  ) THEN
            SM  = USER_SECANTS(IB) ;  MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,IB) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, KPARS
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UXDN, SM, L_mult )
                LP_UT_EMULT_DN(LUM,UT,K,IB,Q) = TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, KPARS
                V2 = L_T_UTDN_USERM(UT,IB,Q) - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                LP_UT_EMULT_DN(LUM,UT,K,IB,Q) = SD * V2
              ENDDO
            ENDIF
          ENDIF

!   Case (b), N > K, profile only

          IF ( N.GT.K ) THEN
            DO Q = 1, KPARS
              V1 = LP_INITIAL_TRANS(N,K,IB,Q)
              LP_UT_EMULT_DN(LUM,UT,K,IB,Q) = UT_EMULT_DN(LUM,UT,IB)*V1
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LP_PARTLAYER_EMULT_OG_DN

      END MODULE vlidort_lp_miscsetups_m

