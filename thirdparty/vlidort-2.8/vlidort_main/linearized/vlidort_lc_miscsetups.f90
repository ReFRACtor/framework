
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
! #            VLIDORT_LAC_MISCSETUPS (master)                  #
! #              VLIDORT_LC_PREPTRANS                           #
! #                                                             #
! #            LC_EMULT_MASTER (master), calling:               #
! #                LC_WHOLELAYER_EMULT_UP                       #
! #                LC_WHOLELAYER_EMULT_DN                       #
! #                LC_PARTLAYER_EMULT_UP                        #
! #                LC_PARTLAYER_EMULT_DN                        #
! #                                                             #
! #            LC_EMULT_MASTER_OBSGEO (master), calling:        #
! #                LC_WHOLELAYER_EMULT_OG_UP                    #
! #                LC_WHOLELAYER_EMULT_OG_DN                    #
! #                LC_PARTLAYER_EMULT_OG_UP                     #
! #                LC_PARTLAYER_EMULT_OG_DN                     #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. No Changes here.

!  Notes for Version 2.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/10/13  - Introduce TAYLOR_LIMIT parameters in module VLIDORT_PARS, instead of "HOPITAL_TOLERANCE"
!     Rob  Fix 05/10/13  - L'Hopitals Rule replaced by Taylor series (original calculation was first term in series!)
!     Rob  Fix 02/19/14  - Final  Taylor series stuff.

      MODULE vlidort_lc_miscsetups_m

!  @@@ Rob Fix 10 May 13 - need the Taylor series routines

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_L_1

      PRIVATE
      PUBLIC  :: VLIDORT_LAC_MISCSETUPS, &
                 LC_EMULT_MASTER,        &
                 LC_EMULT_MASTER_OBSGEO

      CONTAINS

      SUBROUTINE VLIDORT_LAC_MISCSETUPS ( &
        DO_USER_STREAMS, DO_DELTAM_SCALING, DO_PLANE_PARALLEL,            & ! Input flags
        DO_SOLUTION_SAVING,  DO_ATMOS_LINEARIZATION,                      & ! Input flags 
        NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, NBEAMS,               & ! Input numbers
        NMOMENTS, NSTOKES_SQ, N_PARTLAYERS, PARTLAYERS_LAYERIDX,          & ! Input numbers/partials
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                          & ! Input Optical props
        N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,            & ! Input Lin control 
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, DELTAU_SLANT_UNSCALED,    & ! Input derived properties
        PARTAU_SLANT_UNSCALED, LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,    & ! Input derived properties
        OMEGA_GREEK, TRUNC_FACTOR, FAC1,                                  & ! Input derived properties
        MUELLER_INDEX, QUAD_STREAMS, USER_SECANTS,                        & ! Input streams/Muller
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                   & ! Input Dis.Ord. Trans.
        BEAM_CUTOFF, T_DELT_MUBAR, T_UTDN_MUBAR, AVERAGE_SECANT,          & ! Input beam Trans.
        SOLARBEAM_ATRANS, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,       & ! Input User Trans.
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, & ! Input Linearized Optical props.
        L_DELTAU_VERT, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL, L_DELTAU_SLANT,        & ! Output Linearized scaled properties
        DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, L_OMEGA_GREEK,                   & ! Output Linearized scaled properties
        LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Output linearized beam stuff
        LC_LEVELS_SOLARTRANS,  LC_PARTIALS_SOLARTRANS,                         & ! Output linearized beam stuff
        LC_SOLARBEAM_BOATRANS, LC_SOLARBEAM_ATRANS,                            & ! Output linearized beam stuff
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,                  & ! Output linearized DisOrd Trans.
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM )                         ! Output lineraized User Trans.

!mick fix 9/19/2017 -to facilitate correction of linearized direct flux, 
!                    (1) added N_TOTALCOLUMN_WFS, DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED,
!                        LEVELS_SOLARTRANS, & PARTIALS_SOLARTRANS as inputs
!                    (2) added LC_LEVELS_SOLARTRANS & LC_PARTIALS_SOLARTRANS as outputs

!  4/9/19. Version 2.8.1, add BOATRANS, ATRANS and linearization thereof.
        
      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXMOMENTS_INPUT, MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, ZERO


      USE VLIDORT_LA_MISCSETUPS_m

      IMPLICIT NONE

!  Input flags

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

      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Input optical (unscaled)

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Input derived optical

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
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT   ( MAXLAYERS, MAXBEAMS )
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

!  Output Beam parameterization linearized. Version 2.8.1, 4/19/19 add LC_SOLARBEAM_BOATRANS,LC_SOLARBEAM_TRANS

      DOUBLE PRECISION, INTENT (OUT) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: LC_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_SOLARBEAM_BOATRANS  ( MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_SOLARBEAM_ATRANS    ( MAXBEAMS, MAX_ATMOSWFS )

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

!  Miscellaneous setup operations for linearized quantities

!  Deltam scaling of variational quantities

      CALL VLIDORT_LA_DELTAMSCALE ( &
        DO_DELTAM_SCALING, DO_ATMOS_LINEARIZATION,                        & ! Input flags
        NSTOKES, NLAYERS, NBEAMS, NMOMENTS, NSTOKES_SQ,                   & ! Input Numbers
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                               & ! Input linearization control  
        DELTAU_SLANT, TRUNC_FACTOR, FAC1,                                 & ! Input Optical (scaled)
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                          & ! Input Optical
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, & ! Input Optical linearized
        L_DELTAU_VERT, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL,                   & ! Output Optical linearized
        L_DELTAU_SLANT, DO_SCATMAT_VARIATION, L_TRUNC_FACTOR )              ! Output Optical linearized

!  Initialise single scatter albedo variational quantities

      CALL VLIDORT_LA_SSALBINIT ( &
        DO_ATMOS_LINEARIZATION, NSTOKES, NLAYERS, NMOMENTS, & ! Input flag and numbers
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, MUELLER_INDEX,  & ! Input Lin-control
        OMEGA_GREEK, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL,       & ! Inputs others
        L_OMEGA_GREEK )                                       ! Output

!  Linearization of transmittances

      CALL VLIDORT_LA_PREPTRANS ( &
        DO_SOLUTION_SAVING, DO_USER_STREAMS, NLAYERS, NSTREAMS, & ! Input flags and numbers
        N_USER_STREAMS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,      & ! Input numbers
        DELTAU_VERT, PARTAU_VERT, QUAD_STREAMS, USER_SECANTS,   & ! Input Optical and streams
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,         & ! Input Discr-Ord. Trans.
        T_DELT_USERM,   T_UTDN_USERM,   T_UTUP_USERM,           & ! Input User-strm. Trans.
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, L_DELTAU_VERT,      & ! Input linearization
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,   & ! Output linearized Trans (DO)
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM )          ! Output linearized Trans (US)

!  Linearization of pseudo-spherical setup. Additional arguments 4/9/19 (SOLARBEAM_ATRANS)

      CALL VLIDORT_LC_PREPTRANS ( &
        DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,     & ! Input Flag and numbers
        LAYER_VARY_NUMBER, DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,  & ! Input Scaled opticals
        AVERAGE_SECANT, BEAM_CUTOFF, T_DELT_MUBAR, T_UTDN_MUBAR, SOLARBEAM_ATRANS, & ! Input Solar beam
        LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, LC_INITIAL_TRANS,      & ! Output Linearized Solar beam
        LC_AVERAGE_SECANT, LC_SOLARBEAM_ATRANS )                   ! Output Linearized Solar beam

!mick fix 9/19/2017 - added calculation of LC_LEVELS_SOLARTRANS & LC_PARTIALS_SOLARTRANS
!                     to facilitate correction of linearized direct flux (following LIDORT)

!  Whole layers. 4/19/19. Add LC_SOLARBEAM_BOATRANS

      LC_SOLARBEAM_BOATRANS(1:NBEAMS,:) = ZERO
      LC_LEVELS_SOLARTRANS(0:NLAYERS,1:NBEAMS,:) = ZERO
      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
         IF ( LEVELS_SOLARTRANS(N,IB) .GT. ZERO  ) THEN
           DO Q = 1, N_TOTALCOLUMN_WFS
             HELPTRANS = ZERO
             DO K = 1, N
                HELPTRANS = HELPTRANS + DELTAU_SLANT_UNSCALED(N,K,IB)* L_DELTAU_VERT_INPUT(Q,K)
             ENDDO
             LC_LEVELS_SOLARTRANS(N,IB,Q) = - LEVELS_SOLARTRANS(N,IB) * HELPTRANS
             IF ( N == NLAYERS ) LC_SOLARBEAM_BOATRANS(IB,Q) = LC_LEVELS_SOLARTRANS(N,IB,Q)
           ENDDO
         ENDIF
        ENDDO
      ENDDO

!  Partial layers

      LC_PARTIALS_SOLARTRANS(:,1:NBEAMS,:) = ZERO
      IF ( N_PARTLAYERS .GT. 0 ) THEN
        DO IB = 1, NBEAMS
          DO UT = 1, N_PARTLAYERS
            N = PARTLAYERS_LAYERIDX(UT)
            IF ( PARTIALS_SOLARTRANS(UT,IB) .GT. ZERO  ) THEN
             DO Q = 1, N_TOTALCOLUMN_WFS
               HELPTRANS = ZERO
               DO K = 1, N
                  HELPTRANS = HELPTRANS + PARTAU_SLANT_UNSCALED(UT,K,IB)* L_DELTAU_VERT_INPUT(Q,K)
               ENDDO
               LC_PARTIALS_SOLARTRANS(UT,IB,Q) = - PARTIALS_SOLARTRANS(UT,IB) * HELPTRANS
             ENDDO
           ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LAC_MISCSETUPS

!

      SUBROUTINE VLIDORT_LC_PREPTRANS ( &
        DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,     & ! Input Flag and numbers
        LAYER_VARY_NUMBER, DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,  & ! Input Scaled opticals
        AVERAGE_SECANT, BEAM_CUTOFF, T_DELT_MUBAR, T_UTDN_MUBAR, SOLARBEAM_ATRANS, & ! Input Solar beam
        LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, LC_INITIAL_TRANS,      & ! Output Linearized Solar beam
        LC_AVERAGE_SECANT, LC_SOLARBEAM_ATRANS )                   ! Output Linearized Solar beam

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Input flag

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  Input numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS

!  Input partial/linearization control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER   ( MAXLAYERS )

!  Input opticals (scaled)

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Input Beam parameterization

      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      INTEGER         , INTENT (IN) :: BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SOLARBEAM_ATRANS ( MAXBEAMS )

!  Output Beam parameterization linearized

      DOUBLE PRECISION, INTENT (OUT) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Additional output for 2.8.1

      DOUBLE PRECISION, INTENT (OUT) :: LC_SOLARBEAM_ATRANS ( MAXBEAMS, MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER ::          N, Q, UT, K, IB
      DOUBLE PRECISION :: WDEL, VAR, RHO, FAC, DELT, LAMDA, SUM, LOGD

!  linearization of Initial transmittances
!  =======================================

!   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

      DO IB = 1, NBEAMS
        N = 1
        DO Q = 1, LAYER_VARY_NUMBER(N)
          LC_INITIAL_TRANS(N,IB,Q) = ZERO
        ENDDO
        DO N = 2, NLAYERS
          IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              SUM = ZERO  
              DO K = 1, N-1
                SUM = SUM + L_DELTAU_VERT(Q,K)*DELTAU_SLANT(N-1,K,IB)
              ENDDO
              LC_INITIAL_TRANS(N,IB,Q) = - SUM
            ENDDO
          ELSE
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_INITIAL_TRANS(N,IB,Q) = ZERO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          N = 1
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_AVERAGE_SECANT(N,IB,Q) = ZERO
          ENDDO
          DO N = 2, NLAYERS
            IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
              DELT  = DELTAU_VERT(N)
              LAMDA = AVERAGE_SECANT(N,IB)
              FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
              ENDDO
              DO K = 1, N-1
                FAC = ( DELTAU_SLANT(N,K,IB) - DELTAU_SLANT(N-1,K,IB) ) / DELT
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LC_AVERAGE_SECANT(N,IB,Q) = LC_AVERAGE_SECANT(N,IB,Q) + L_DELTAU_VERT(Q,K)*FAC
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = ZERO
              ENDDO
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
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF  ( N.LE.BEAM_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_DELT_MUBAR(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.BEAM_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  end layer and beam loops

        ENDDO
      ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  =================================================================

      DO IB = 1, NBEAMS

!  zero it

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_T_UTDN_MUBAR(UT,IB,Q) = ZERO
          ENDDO
        ENDDO

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
         FAC = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF ( N .LE. BEAM_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_UTDN_MUBAR(UT,IB,Q) = L_DELTAU_VERT(Q,N)* FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF ( N .LE. BEAM_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  End optical depth and beam loops

        ENDDO
      ENDDO

!  4/29/19 Linearization of ATRANS.

      DO IB = 1, NBEAMS
         DO Q = 1, LAYER_VARY_NUMBER(NLAYERS)
            LOGD = ZERO
            DO K = 1, NLAYERS
               LOGD = LOGD - L_DELTAU_VERT(Q,K)*DELTAU_SLANT(NLAYERS,K,IB)
            ENDDO
            LC_SOLARBEAM_ATRANS(IB,Q) = SOLARBEAM_ATRANS(IB) * LOGD
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LC_PREPTRANS

!

      SUBROUTINE LC_EMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, TAYLOR_ORDER,      & ! Input flags + Taylor
        NLAYERS, NBEAMS, N_USER_STREAMS, N_PARTLAYERS, N_COLUMNWFS,       & ! Input numbers
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input Level control
        USER_SECANTS, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,            & ! Input optical, streams   
        BEAM_CUTOFF, T_DELT_MUBAR, EMULT_HOPRULE, SIGMA_M, SIGMA_P,       & ! Input solar beam + Multipliers
        T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, ITRANS_USERM,           & ! Input User-stream Trans.
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                     & ! Input Multipliers
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,             & ! Input linearized Beam stuff
        LC_T_UTDN_MUBAR, L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,  & ! Input linearized transmittanaces
        LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN )          ! Output linearized Mulitpliers

!  Linearized multipliers for the Beam source terms

!  Version 2p7 notes
!     Rob Fix 05/10/13  - Introduce Taylor series routines (first version) HOPITAL_TOLERANCE replaced by TAYLOR_SMALL.
!     Rob Fix 02/19/14  - Small numbers analysis finalized using Taylor_Order parameter, as for LIDORT

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                 MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

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
      INTEGER, INTENT (IN) ::           N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Linearized User transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM  ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  OUTPUTS
!  =======

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (OUT) :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER :: N, UT

!mick fix 7/23/2014 - initialized for packing
      LC_EMULT_UP = ZERO
      LC_EMULT_DN = ZERO
      LC_UT_EMULT_UP = ZERO
      LC_UT_EMULT_DN = ZERO

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!  Whole layer upwelling
!  ---------------------

!  Loop over all  model  layers N
!    Profiles:  loop over all varying layers K such that K </= N
!    Columns :  K = 0

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           CALL LC_WHOLELAYER_EMULT_UP ( &
             DO_PLANE_PARALLEL, N, NBEAMS, N_USER_STREAMS, N_COLUMNWFS, & ! Input flag/numbers
             BEAM_CUTOFF, T_DELT_MUBAR, EMULT_UP, SIGMA_P,              & ! Input Beam stuff + Multipliers
             T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,                & ! Input User trans.
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,      & ! Input linearized Beam stuff
             LC_EMULT_UP )                                                ! Output Multiplier
         ENDIF
       ENDDO

!  Partial layer upwelling
!  -----------------------

!  Start loop over all partial output UT occuring in layers N
!  Start loop over all varying layers K such that K </= N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           CALL LC_PARTLAYER_EMULT_UP ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, N_USER_STREAMS, N_COLUMNWFS,         & ! Input flag/numbers
             BEAM_CUTOFF, T_DELT_MUBAR, UT_EMULT_UP, SIGMA_P,                       & ! Input Beam stuff, multipliers
             T_UTUP_USERM, ITRANS_USERM, L_T_UTUP_USERM,                            & ! Input Multipliers + linearization
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Input linearized Beam stuff
             LC_UT_EMULT_UP )                                                         ! Output Multiplier
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
           CALL LC_WHOLELAYER_EMULT_DN ( &
             DO_PLANE_PARALLEL, N, NBEAMS, N_USER_STREAMS, N_COLUMNWFS, TAYLOR_ORDER, & ! Input flag/numbers
             USER_SECANTS, BEAM_CUTOFF, EMULT_DN, EMULT_HOPRULE, SIGMA_M,             & ! Input Beam stuff + Multipliers
             DELTAU_VERT, L_DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,  & ! Input User trans.
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                    & ! Input linearized Beam stuff
             LC_EMULT_DN )                                                              ! Output Multiplier
         ENDIF
       ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           CALL LC_PARTLAYER_EMULT_DN ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, N_USER_STREAMS, N_COLUMNWFS, TAYLOR_ORDER, & ! Input flag/numbers
             USER_SECANTS, BEAM_CUTOFF, UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M,              & ! Input Beam stuff, multipliers
             PARTAU_VERT, L_DELTAU_VERT, T_UTDN_USERM, ITRANS_USERM, L_T_UTDN_USERM,      & ! Input User transmittances
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_UTDN_MUBAR,                        & ! Input linearized Beam stuff
             LC_UT_EMULT_DN )                                                               ! Output Multiplier
         ENDIF
       ENDDO

!  end downwelling

      ENDIF
  
!  Finish

      RETURN
      END SUBROUTINE LC_EMULT_MASTER

!

      SUBROUTINE LC_WHOLELAYER_EMULT_UP ( &
         DO_PLANE_PARALLEL, N, NBEAMS, N_USER_STREAMS, N_COLUMNWFS, & ! Input flag/numbers
         BEAM_CUTOFF, T_DELT_MUBAR, EMULT_UP, SIGMA_P,              & ! Input Beam stuff + Multipliers
         T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,                & ! Input User trans.
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,      & ! Input linearized Beam stuff
         LC_EMULT_UP )                                                ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output linearized multiplier

      DOUBLE PRECISION, INTENT (INOUT) :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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
           DO Q = 1, N_COLUMNWFS
             LC_EMULT_UP(UM,N,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, N_COLUMNWFS
              V1 = -LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,UM,IB)
              V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_UP(UM,N,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, N_COLUMNWFS
              V1 = LC_INITIAL_TRANS (N,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_UP(UM,N,IB,Q) = EMULT_UP(UM,N,IB)*V1 + SU * V2
            ENDDO
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_EMULT_UP

!

      SUBROUTINE LC_WHOLELAYER_EMULT_DN ( &
         DO_PLANE_PARALLEL, N, NBEAMS, N_USER_STREAMS, N_COLUMNWFS, TAYLOR_ORDER, & ! Input flag/numbers
         USER_SECANTS, BEAM_CUTOFF, EMULT_DN, EMULT_HOPRULE, SIGMA_M,             & ! Input Beam stuff + Multipliers
         DELTAU_VERT, L_DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,  & ! Input User trans.
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                    & ! Input linearized Beam stuff
         LC_EMULT_DN )                                                              ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL


      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output linearized multiplier

      DOUBLE PRECISION, INTENT (INOUT) :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, DELTA, MULT, TMEW, L_DELTA, UDEL, L_LAM, L_MULT
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. BEAM_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, N_COLUMNWFS
             LC_EMULT_DN(UM,N,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Note the use of L'Hopital's Rule flag.

         DO UM = 1, N_USER_STREAMS
           SM   = USER_SECANTS(UM) ; MULT = EMULT_DN(UM,N,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
           IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
             UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,UM,IB) ; DELTA = DELTAU_VERT(N)
             DO Q = 1, N_COLUMNWFS
               L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
               L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
               CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
               LC_EMULT_DN(UM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
             ENDDO
           ELSE
             SD = TMEW / SIGMA_M(N,UM,IB)
             DO Q = 1, N_COLUMNWFS
               V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,UM,IB)
               V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
               V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
               LC_EMULT_DN(UM,N,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
             ENDDO
           ENDIF
         ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

         DO UM = 1, N_USER_STREAMS
           SM   = USER_SECANTS(UM)   ; MULT = EMULT_DN(UM,N,IB)  ; TMEW  = ITRANS_USERM(N,UM,IB) 
           IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
             UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,UM,IB) ; DELTA = DELTAU_VERT(N)
             DO Q = 1, N_COLUMNWFS
               L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
               CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_mult )
               LC_EMULT_DN(UM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
             ENDDO
           ELSE
             SD = TMEW / SIGMA_M(N,UM,IB)
             DO Q = 1, N_COLUMNWFS
               V1 = LC_INITIAL_TRANS (N,IB,Q)
               V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
               LC_EMULT_DN(UM,N,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
             ENDDO
           ENDIF
         ENDDO

!  End clause pseudo-spherical versus plaen-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_EMULT_DN

!

      SUBROUTINE LC_PARTLAYER_EMULT_UP ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, N_USER_STREAMS, N_COLUMNWFS,         & ! Input flag/numbers
             BEAM_CUTOFF, T_DELT_MUBAR, UT_EMULT_UP, SIGMA_P,                       & ! Input Beam stuff, multipliers
             T_UTUP_USERM, ITRANS_USERM, L_T_UTUP_USERM,                            & ! Input Multipliers + linearization
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Input linearized Beam stuff
             LC_UT_EMULT_UP )                                                         ! Output Multiplier


      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Output multipliers

      DOUBLE PRECISION, INTENT (INOUT) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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
           DO Q = 1, N_COLUMNWFS
             LC_UT_EMULT_UP(UM,UT,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, N_COLUMNWFS
              V1 = -LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,UM,IB)
              V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
              V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - UX_UP * LC_T_DELT_MUBAR(N, IB,Q) - WDEL * L_T_UTUP_USERM(UT,UM,Q)
              LC_UT_EMULT_UP(UM,UT,IB,Q) = SU * V2 + UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, N_COLUMNWFS
              V1 = LC_INITIAL_TRANS (N,IB,Q)
              V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - UX_UP * LC_T_DELT_MUBAR( N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,UM,Q)
              LC_UT_EMULT_UP(UM,UT,IB,Q) = SU * V2 + UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO

!  End clause pseudo-spherical versus plaen-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_EMULT_UP

!

      SUBROUTINE LC_PARTLAYER_EMULT_DN ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, N_USER_STREAMS, N_COLUMNWFS, TAYLOR_ORDER, & ! Input flag/numbers
             USER_SECANTS, BEAM_CUTOFF, UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M,              & ! Input Beam stuff, multipliers
             PARTAU_VERT, L_DELTAU_VERT, T_UTDN_USERM, ITRANS_USERM, L_T_UTDN_USERM,      & ! Input User transmittances
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_UTDN_MUBAR,                        & ! Input linearized Beam stuff
             LC_UT_EMULT_DN )                                                               ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                 MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Output Multiplier
!mick fix 9/19/2017 - changed intent from "out" to "inout"

      DOUBLE PRECISION, INTENT (INOUT) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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
           DO Q = 1, N_COLUMNWFS
             LC_UT_EMULT_DN(UM,UT,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            SM    = USER_SECANTS(UM) ; MULT = UT_EMULT_DN(UM,UT,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,UM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, N_COLUMNWFS
                L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UXDN, SM, L_mult )
                LC_UT_EMULT_DN(UM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,UM,IB)
              DO Q = 1, N_COLUMNWFS
                V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,UM,IB)
                V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
                LC_UT_EMULT_DN(UM,UT,IB,Q) =       SD * V2 + &
                                 UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            SM    = USER_SECANTS(UM) ; MULT = UT_EMULT_DN(UM,UT,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,UM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, N_COLUMNWFS
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UXDN, SM, L_mult )
                LC_UT_EMULT_DN(UM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,UM,IB)
              DO Q = 1, N_COLUMNWFS
                V1 = LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
                LC_UT_EMULT_DN(UM,UT,IB,Q) =       SD * V2 + &
                                 UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_EMULT_DN

!

      SUBROUTINE LC_EMULT_MASTER_OBSGEO ( &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, TAYLOR_ORDER,      & ! Input flags + Taylor
        NLAYERS, NBEAMS, N_PARTLAYERS, N_COLUMNWFS,                       & ! Input numbers
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input Level control
        USER_SECANTS, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,            & ! Input optical, streams   
        BEAM_CUTOFF, T_DELT_MUBAR, EMULT_HOPRULE, SIGMA_M, SIGMA_P,       & ! Input solar beam + Multipliers
        T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, ITRANS_USERM,           & ! Input User-stream Trans.
        EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                     & ! Input Multipliers
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,             & ! Input linearized Beam stuff
        LC_T_UTDN_MUBAR, L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,  & ! Input linearized transmittanaces
        LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN )          ! Output linearized Mulitpliers

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
      INTEGER, INTENT (IN) ::           N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Linearized User transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM  ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  OUTPUTS
!  =======

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (OUT) :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER :: N, UT

!mick fix 7/23/2014 - initialized for packing
      LC_EMULT_UP = ZERO
      LC_EMULT_DN = ZERO
      LC_UT_EMULT_UP = ZERO
      LC_UT_EMULT_DN = ZERO

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!  Whole layer upwelling
!  ---------------------

!  Loop over all  model  layers N
!    Profiles:  loop over all varying layers K such that K </= N
!    Columns :  K = 0

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           CALL LC_WHOLELAYER_EMULT_OG_UP ( &
             DO_PLANE_PARALLEL, N, NBEAMS, N_COLUMNWFS,             & ! Input flag/numbers
             BEAM_CUTOFF, T_DELT_MUBAR, EMULT_UP, SIGMA_P,          & ! Input Beam stuff + Multipliers
             T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,            & ! Input User trans.
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,  & ! Input linearized Beam stuff
             LC_EMULT_UP )                                                ! Output Multiplier
         ENDIF
       ENDDO

!  Partial layer upwelling
!  -----------------------

!  Start loop over all partial output UT occuring in layers N
!  Start loop over all varying layers K such that K </= N

       DO UT = 1, N_PARTLAYERS
         N   = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           CALL LC_PARTLAYER_EMULT_OG_UP ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, N_COLUMNWFS,                         & ! Input flag/numbers
             BEAM_CUTOFF, T_DELT_MUBAR, UT_EMULT_UP, SIGMA_P,                       & ! Input Beam stuff, multipliers
             T_UTUP_USERM, ITRANS_USERM, L_T_UTUP_USERM,                            & ! Input Multipliers + linearization
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Input linearized Beam stuff
             LC_UT_EMULT_UP )                                                         ! Output Multiplier
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
!  Start loop over all varying layers K such that K </= N

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           CALL LC_WHOLELAYER_EMULT_OG_DN ( &
             DO_PLANE_PARALLEL, N, NBEAMS, N_COLUMNWFS, TAYLOR_ORDER,                & ! Input flag/numbers
             USER_SECANTS, BEAM_CUTOFF, EMULT_DN, EMULT_HOPRULE, SIGMA_M,            & ! Input Beam stuff + Multipliers
             DELTAU_VERT, L_DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM, & ! Input User trans.
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                   & ! Input linearized Beam stuff
             LC_EMULT_DN )                                                             ! Output Multiplier
         ENDIF
       ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N
!  Start loop over all varying layers K such that K </= N

       DO UT = 1, N_PARTLAYERS
         N   = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           CALL LC_PARTLAYER_EMULT_OG_DN ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, N_COLUMNWFS, TAYLOR_ORDER,            & ! Input flag/numbers
             USER_SECANTS, BEAM_CUTOFF, UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M,         & ! Input Beam stuff, multipliers
             PARTAU_VERT, L_DELTAU_VERT, T_UTDN_USERM, ITRANS_USERM, L_T_UTDN_USERM, & ! Input User transmittances
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_UTDN_MUBAR,                   & ! Input linearized Beam stuff
             LC_UT_EMULT_DN )                                                          ! Output Multiplier
         ENDIF
       ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_EMULT_MASTER_OBSGEO

!

      SUBROUTINE LC_WHOLELAYER_EMULT_OG_UP ( &
         DO_PLANE_PARALLEL, N, NBEAMS, N_COLUMNWFS,            & ! Input flag/numbers
         BEAM_CUTOFF, T_DELT_MUBAR, EMULT_UP, SIGMA_P,         & ! Input Beam stuff + Multipliers
         T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM,           & ! Input User trans.
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, & ! Input linearized Beam stuff
         LC_EMULT_UP )                                           ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output linearized multiplier

      DOUBLE PRECISION, INTENT (INOUT) :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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

         DO Q = 1, N_COLUMNWFS
           LC_EMULT_UP(LUM,N,IB,Q) = ZERO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          UDEL = T_DELT_USERM(N,IB)
          SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
          DO Q = 1, N_COLUMNWFS
            V1 = -LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,LUM,IB)
            V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
            V2 = WDEL * L_T_DELT_USERM(N,IB,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
            LC_EMULT_UP(LUM,N,IB,Q) = EMULT_UP(LUM,N,IB) * V1 + SU * V2
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          UDEL = T_DELT_USERM(N,IB)
          SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
          DO Q = 1, N_COLUMNWFS
            V1 = LC_INITIAL_TRANS (N,IB,Q)
            V2 = WDEL * L_T_DELT_USERM(N,IB,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
            LC_EMULT_UP(LUM,N,IB,Q) = EMULT_UP(LUM,N,IB)*V1 + SU * V2
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_EMULT_OG_UP

!

      SUBROUTINE LC_WHOLELAYER_EMULT_OG_DN ( &
         DO_PLANE_PARALLEL, N, NBEAMS, N_COLUMNWFS, TAYLOR_ORDER,                & ! Input flag/numbers
         USER_SECANTS, BEAM_CUTOFF, EMULT_DN, EMULT_HOPRULE, SIGMA_M,            & ! Input Beam stuff + Multipliers
         DELTAU_VERT, L_DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, L_T_DELT_USERM, & ! Input User trans.
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                   & ! Input linearized Beam stuff
         LC_EMULT_DN )                                                             ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL


      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output linearized multiplier

      DOUBLE PRECISION, INTENT (INOUT) :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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

         DO Q = 1, N_COLUMNWFS
           LC_EMULT_DN(LUM,N,IB,Q) = ZERO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          SM    = USER_SECANTS(IB) ; MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
          IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
            UDEL = T_DELT_USERM(N,IB) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
            DO Q = 1, N_COLUMNWFS
              L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
              LC_EMULT_DN(LUM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
            ENDDO
          ELSE
            SD = TMEW / SIGMA_M(N,LUM,IB)
            DO Q = 1, N_COLUMNWFS
              V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,LUM,IB)
              V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
              V2 = L_T_DELT_USERM(N,IB,Q) - LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_DN(LUM,N,IB,Q) = EMULT_DN(LUM,N,IB)*V1 + SD*V2
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

          SM    = USER_SECANTS(IB) ; MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
          IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
            UDEL = T_DELT_USERM(N,IB) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
            DO Q = 1, N_COLUMNWFS
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_mult )
              LC_EMULT_DN(LUM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
            ENDDO
          ELSE
            SD = TMEW / SIGMA_M(N,LUM,IB)
            DO Q = 1, N_COLUMNWFS
              V1 = ZERO
              V1 = LC_INITIAL_TRANS(N,IB,Q)
              V2 = L_T_DELT_USERM(N,IB,Q) - LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_DN(LUM,N,IB,Q) = EMULT_DN(LUM,N,IB)*V1 + SD*V2
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
      END SUBROUTINE LC_WHOLELAYER_EMULT_OG_DN

!

      SUBROUTINE LC_PARTLAYER_EMULT_OG_UP ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, N_COLUMNWFS,                         & ! Input flag/numbers
             BEAM_CUTOFF, T_DELT_MUBAR, UT_EMULT_UP, SIGMA_P,                       & ! Input Beam stuff, multipliers
             T_UTUP_USERM, ITRANS_USERM, L_T_UTUP_USERM,                            & ! Input Multipliers + linearization
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Input linearized Beam stuff
             LC_UT_EMULT_UP )                                                         ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Output multipliers

      DOUBLE PRECISION, INTENT (INOUT) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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

         DO Q = 1, N_COLUMNWFS
           LC_UT_EMULT_UP(LUM,UT,IB,Q) = ZERO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          UX_UP = T_UTUP_USERM(UT,IB)
          SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
          DO Q = 1, N_COLUMNWFS
            V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,LUM,IB)
            V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
            V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - UX_UP * LC_T_DELT_MUBAR(N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,IB,Q)
            LC_UT_EMULT_UP(LUM,UT,IB,Q) = SU * V2 + UT_EMULT_UP(LUM,UT,IB) * V1
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          UX_UP = T_UTUP_USERM(UT,IB)
          SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
          DO Q = 1, N_COLUMNWFS
            V1 = LC_INITIAL_TRANS(N,IB,Q)
            V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - UX_UP * LC_T_DELT_MUBAR(N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,IB,Q)
            LC_UT_EMULT_UP(LUM,UT,IB,Q) =  SU * V2 + UT_EMULT_UP(LUM,UT,IB) * V1
          ENDDO

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_EMULT_OG_UP

!

      SUBROUTINE LC_PARTLAYER_EMULT_OG_DN ( &
             DO_PLANE_PARALLEL, N, UT, NBEAMS, N_COLUMNWFS, TAYLOR_ORDER,            & ! Input flag/numbers
             USER_SECANTS, BEAM_CUTOFF, UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M,         & ! Input Beam stuff, multipliers
             PARTAU_VERT, L_DELTAU_VERT, T_UTDN_USERM, ITRANS_USERM, L_T_UTDN_USERM, & ! Input User transmittances
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_UTDN_MUBAR,                   & ! Input linearized Beam stuff
             LC_UT_EMULT_DN )                                                          ! Output Multiplier

      USE VLIDORT_PARS_M, Only : MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                 MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_COLUMNWFS

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

      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  Output Multiplier
!mick fix 9/19/2017 - changed intent from "out" to "inout"

      DOUBLE PRECISION, INTENT (INOUT) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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

         DO Q = 1, N_COLUMNWFS
           LC_UT_EMULT_DN(LUM,UT,IB,Q) = ZERO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          SM  = USER_SECANTS(IB) ;  MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
          IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
            UXDN = T_UTDN_USERM(UT,IB) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
            DO Q = 1, N_COLUMNWFS
              L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UXDN, SM, L_mult )
              LC_UT_EMULT_DN(LUM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
            ENDDO
          ELSE
            SD = TMEW / SIGMA_M(N,LUM,IB)
            DO Q = 1, N_COLUMNWFS
              V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,LUM,IB)
              V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
              V2 = L_T_UTDN_USERM(UT,IB,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
              LC_UT_EMULT_DN(LUM,UT,IB,Q) = SD * V2 + UT_EMULT_DN(LUM,UT,IB) * V1
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

          SM  = USER_SECANTS(IB) ;  MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
          IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
            UXDN = T_UTDN_USERM(UT,IB) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
            DO Q = 1, N_COLUMNWFS
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UXDN, SM, L_mult )
              LC_UT_EMULT_DN(LUM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
            ENDDO
          ELSE
            DO Q = 1, N_COLUMNWFS
              V1 = LC_INITIAL_TRANS(N,IB,Q)
              V2 = L_T_UTDN_USERM(UT,IB,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
              LC_UT_EMULT_DN(LUM,UT,IB,Q) = SD * V2 + UT_EMULT_DN(LUM,UT,IB) * V1
            ENDDO
          ENDIF

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_EMULT_OG_DN

      END MODULE vlidort_lc_miscsetups_m

