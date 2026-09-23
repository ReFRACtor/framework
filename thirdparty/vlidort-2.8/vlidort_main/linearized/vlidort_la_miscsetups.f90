
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
! #              VLIDORT_LA_DELTAMSCALE                         #
! #              VLIDORT_LA_SSALBINIT                           #
! #              VLIDORT_LA_PREPTRANS                           #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. No Changes.

      MODULE vlidort_la_miscsetups_m

      PUBLIC :: VLIDORT_LA_DELTAMSCALE, &
                VLIDORT_LA_SSALBINIT,   &
                VLIDORT_LA_PREPTRANS

      CONTAINS

      SUBROUTINE VLIDORT_LA_DELTAMSCALE ( &
        DO_DELTAM_SCALING, DO_ATMOS_LINEARIZATION,                        & ! Input flags
        NSTOKES, NLAYERS, NBEAMS, NMOMENTS, NSTOKES_SQ,                   & ! Input Numbers
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                               & ! Input linearization control  
        DELTAU_SLANT, TRUNC_FACTOR, FAC1,                                 & ! Input Optical (scaled)
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                          & ! Input Optical
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, & ! Input Optical linearized
        L_DELTAU_VERT, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL,                   & ! Output Optical linearized
        L_DELTAU_SLANT, DO_SCATMAT_VARIATION, L_TRUNC_FACTOR )              ! Output Optical linearized

!  Linearization of the deltam scaling

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAX_ATMOSWFS, MAXMOMENTS_INPUT, MAXMOMENTS, MAXLAYERS, &
                                 MAXBEAMS, MAXSTOKES_SQ, zero, one, SMALLNUM

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION

!  Input numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTOKES_SQ
      INTEGER, INTENT (IN) ::          NBEAMS

!  Input linearization control

      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input scaled optical (from delta-M scaling)

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: FAC1         ( MAXLAYERS )

!  Input Optical (original)

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Input Lnearized Optical (original)

      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  output - scaled linearization quantities

      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      LOGICAL, INTENT (OUT) ::          DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_SLANT ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )

!  local variables

      DOUBLE PRECISION :: BLD, DL, OF1, F, F1, UQ, EQ
      DOUBLE PRECISION :: T1, L_FAC1, DELS
      DOUBLE PRECISION :: FZM, ZMQ, ZLQ, UZQ_SCALE, ZQ_SCALE, UQ_SCALE
      INTEGER          :: N, Q, NM1, L, K, IB, K1, K2, O1
      LOGICAL          :: LOOP

!  Indexing

      INTEGER :: KTYPE1(4), KTYPE2(4)

      KTYPE1 = (/ 1, 6, 11, 16 /)
      KTYPE2 = (/ 2, 5, 12, 15 /)

!  slant optical thickness values
!    (Revert to the input values of deltau_slant)
!   Commented out line is wrong !!!!!!!!!!!!!!!!!!!!!!

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
!            DELS = DELTAU_SLANT(N,K,IB)/DELTAU_VERT_INPUT(K)/FAC1(K)
            DELS = DELTAU_SLANT(N,K,IB)/FAC1(K)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_DELTAU_SLANT(Q,N,K,IB) = L_DELTAU_VERT_INPUT(Q,K) * DELS
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!  Set Scattering matrix linearization flag
!   Only examine the (1,1) entry (the phase function coefficient)
!    Dimensioning bug, 21 March 2007. Care with using NM1

      K1 = 1
      IF ( DO_DELTAM_SCALING ) THEN
        NM1 = NMOMENTS+1
      ELSE
        NM1 = NMOMENTS
      ENDIF

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LOOP = .TRUE.
              L = 0
              DO WHILE (LOOP.AND.L.LT.NM1)
                L = L + 1
                LOOP = ( ABS(L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1)) .LT. 1000.0*SMALLNUM )
              ENDDO
              DO_SCATMAT_VARIATION(N,Q) = .NOT.LOOP
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  DELTAM SCALING
!  ==============

!mick fix 2/17/11 - initialize
      L_GREEKMAT_TOTAL = ZERO

!  New section added 21 December 2005 by R. Spurr
!   Linearization of the scaled delta-M inputs

      IF ( DO_DELTAM_SCALING ) THEN

        IF ( DO_ATMOS_LINEARIZATION ) THEN

          NM1 = NMOMENTS+1

!  start loop over layers

          DO N = 1, NLAYERS

!  If the layer is varying

            IF ( LAYER_VARY_FLAG(N) ) THEN

              OF1 = ( ONE - FAC1(N) ) / FAC1(N)
              F   = TRUNC_FACTOR(N)
              F1  = ONE - F

!  Start loop over varying parameters in this layer

              DO Q = 1, LAYER_VARY_NUMBER(N)

!  scale scattering matrix linearization additionally

                IF ( DO_SCATMAT_VARIATION(N,Q) ) THEN

!  set bulk property values

                  UQ  = L_OMEGA_TOTAL_INPUT(Q,N)
                  EQ  = L_DELTAU_VERT_INPUT(Q,N)
                  ZMQ = L_GREEKMAT_TOTAL_INPUT(Q,NM1,N,1)
                  FZM = F * ZMQ
                  UZQ_SCALE = ( UQ + ZMQ ) * OF1
                  ZQ_SCALE  = FZM / F1
                  L_TRUNC_FACTOR(Q,N) = FZM
                  L_OMEGA_TOTAL(Q,N)      = UQ + UZQ_SCALE - ZQ_SCALE
                  L_DELTAU_VERT(Q,N)      = EQ - UZQ_SCALE

!  do the phase function first
!Rob Fix. 1/17/14. Cannot have BLD anf F both zero.

                  K1 = 1
                  L_GREEKMAT_TOTAL(Q,0,N,K1) = ZERO
                  FZM = F * L_GREEKMAT_TOTAL_INPUT(Q,NM1,N,1)
                  DO L = 1, NMOMENTS
                    DL  = DBLE(2*L+1)
                    ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1)
                    BLD = GREEKMAT_TOTAL_INPUT(L,N,K1) / DL
                    IF ( BLD.ne.ZERO.or.F.ne.zero ) THEN
                      T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                      L_GREEKMAT_TOTAL(Q,L,N,K1) = T1 + ZQ_SCALE
                    ENDIF
                  ENDDO

!  If there is polarization.........

                  IF ( NSTOKES .GT. 1 ) THEN

!  Now do the other Type 1 elements
!Rob Fix. 1/17/14. Cannot have BLD anf F both zero.

                    DO L = 0, NMOMENTS
                      DL  = DBLE(2*L+1)
                      DO K = 2, 4
                        K1 = KTYPE1(K)
                        ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1)
                        BLD = GREEKMAT_TOTAL_INPUT(L,N,K1) / DL
                        IF ( BLD.ne.ZERO.or.F.ne.zero ) THEN
                          T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                          L_GREEKMAT_TOTAL(Q,L,N,K1) = T1 + ZQ_SCALE
                        ENDIF
                      ENDDO
                    ENDDO

!  Now do the Type 2 elements
!    Bug discovered 24 September 2007
!    Formerly, we were overwriting the Type 1 elements

                    DO L = 0, NMOMENTS
                      DL  = DBLE(2*L+1)
                      DO K = 1, 4
                        K2 = KTYPE2(K)
                        ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K2)
!  !!!! WRONG           L_GREEKMAT_TOTAL(Q,L,N,K1) = ZLQ + ZQ_SCALE
                        L_GREEKMAT_TOTAL(Q,L,N,K2) = ZLQ + ZQ_SCALE
                      ENDDO
                    ENDDO

!  Finish scaled Greek matrix linearization

                  ENDIF

!  No scattering matrix linearization

                ELSE

!  Bulk property linearization

                  UQ = L_OMEGA_TOTAL_INPUT(Q,N)
                  EQ = L_DELTAU_VERT_INPUT(Q,N)
                  L_TRUNC_FACTOR(Q,N) = ZERO
                  UQ_SCALE = UQ * OF1
                  L_OMEGA_TOTAL(Q,N) = UQ + UQ_SCALE
                  L_DELTAU_VERT(Q,N) = EQ - UQ_SCALE

!   Zero all linearized scattering matrix quantities now;

                  DO L = 0, NMOMENTS
                   DO O1 = 1, NSTOKES_SQ
                    L_GREEKMAT_TOTAL(Q,L,N,O1) = ZERO
                   ENDDO
                  ENDDO

                ENDIF

              ENDDO
            ENDIF
          ENDDO

!  Scale linearized slant path optical thickness values

          DO N = 1, NLAYERS
           DO K = 1, N

!  Buggy code discovered, 30 October 2007
!            IF ( LAYER_VARY_FLAG(K) ) THEN
!              DO Q = 1, LAYER_VARY_NUMBER(K)
!                L_FAC1 = L_TRUNC_FACTOR(Q,K) *   OMEGA_TOTAL_INPUT(K)
!     &                   + TRUNC_FACTOR(K)   * L_OMEGA_TOTAL_INPUT(Q,K)
!                DO IB = 1, NBEAMS
!                  DELS = DELTAU_SLANT(N,K,IB)/FAC1(K)
!                  L_DELTAU_SLANT(Q,N,K,IB) = - L_FAC1 * DELS +
!     &                           FAC1(K) * L_DELTAU_SLANT(Q,N,K,IB)
!                ENDDO
!              ENDDO
!            ENDIF

            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_FAC1 = L_TRUNC_FACTOR(Q,K) + TRUNC_FACTOR(K)   * L_OMEGA_TOTAL_INPUT(Q,K)
                L_FAC1 = L_FAC1 *  OMEGA_TOTAL_INPUT(K)
                DO IB = 1, NBEAMS
                  DELS = DELTAU_SLANT(N,K,IB)/FAC1(K)
                  L_DELTAU_SLANT(Q,N,K,IB) = - L_FAC1 * DELS + FAC1(K) * L_DELTAU_SLANT(Q,N,K,IB)
                ENDDO
              ENDDO
            ENDIF

          ENDDO
         ENDDO

        ENDIF

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

        IF ( DO_ATMOS_LINEARIZATION ) THEN

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_TRUNC_FACTOR(Q,N) = ZERO
                L_DELTAU_VERT(Q,N)  = L_DELTAU_VERT_INPUT(Q,N)
                L_OMEGA_TOTAL(Q,N)  = L_OMEGA_TOTAL_INPUT(Q,N)
                IF ( DO_SCATMAT_VARIATION(N,Q) ) THEN
                  DO L = 0, MAXMOMENTS
                    DO O1 = 1, NSTOKES_SQ
                      L_GREEKMAT_TOTAL(Q,L,N,O1) = L_GREEKMAT_TOTAL_INPUT(Q,L,N,O1)
                    ENDDO
                  ENDDO
                ELSE
                  DO L = 0, MAXMOMENTS
                    DO O1 = 1, NSTOKES_SQ
                      L_GREEKMAT_TOTAL(Q,L,N,O1) = ZERO
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ENDDO

        ENDIF

      ENDIF

!  Finish module

      RETURN
      END SUBROUTINE VLIDORT_LA_DELTAMSCALE

!

      SUBROUTINE VLIDORT_LA_SSALBINIT ( &
        DO_ATMOS_LINEARIZATION, NSTOKES, NLAYERS, NMOMENTS, & ! Input flag and numbers
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, MUELLER_INDEX,  & ! Input Lin-control
        OMEGA_GREEK, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL,       & ! Inputs others
        L_OMEGA_GREEK )                                       ! Output

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES_SQ, MAX_ATMOSWFS

      IMPLICIT NONE

!  Input flag and numbers

      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NMOMENTS

!  Input linearization control

      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  Other inputs

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  scaled input

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, L, Q, O1, O2, OM
      DOUBLE PRECISION :: VL

!  phase moment-weighted OMEGA and linearizations
!  Including phase function linearization

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO L = 0, NMOMENTS
            DO O1 = 1, NSTOKES
             DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              VL = L_OMEGA_TOTAL(Q,N)   + L_GREEKMAT_TOTAL(Q,L,N,OM)
              L_OMEGA_GREEK(L,N,O1,O2,Q) = OMEGA_GREEK(L,N,O1,O2) * VL
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LA_SSALBINIT

!

      SUBROUTINE VLIDORT_LA_PREPTRANS ( &
        DO_SOLUTION_SAVING, DO_USER_STREAMS, NLAYERS, NSTREAMS, & ! Input flags and numbers
        N_USER_STREAMS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,      & ! Input numbers
        DELTAU_VERT, PARTAU_VERT, QUAD_STREAMS, USER_SECANTS,   & ! Input Optical and streams
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,         & ! Input Discr-Ord. Trans.
        T_DELT_USERM,   T_UTDN_USERM,   T_UTUP_USERM,           & ! Input User-strm. Trans.
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, L_DELTAU_VERT,      & ! Input linearization
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,   & ! Output linearized Trans (DO)
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM )          ! Output linearized Trans (US)

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, &
                                 MAX_ATMOSWFS, MAX_USER_STREAMS, ZERO 

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS

!  Input numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Input optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )

!  Input streams

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS  ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Input linearization control and optical

      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Input Discrete-Ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  Input User stream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  output Linearized transmittances

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS,       MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DELT_USERM ( MAXLAYERS,       MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, Q, UT, UM, I
      DOUBLE PRECISION :: VD, VU, TRANS, UX, TRANS_D, TRANS_U
      DOUBLE PRECISION :: XT, L_TAU, L_TDEL, L_TD, L_TU, LDN, LUP

!  Linearization of discrete ordinate transmittances
!  Required for the solution saving option

!mick fix 7/23/2014 - initialized for packing
        L_T_DELT_DISORDS = ZERO
        L_T_DISORDS_UTDN = ZERO
        L_T_DISORDS_UTUP = ZERO

        L_T_DELT_USERM = ZERO
        L_T_UTDN_USERM = ZERO
        L_T_UTUP_USERM = ZERO

      IF ( DO_SOLUTION_SAVING ) THEN

!  Whole layers

        DO N = 1, NLAYERS
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TAU = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N)
            DO I = 1, NSTREAMS
              L_TDEL = - L_TAU / QUAD_STREAMS(I)
              L_T_DELT_DISORDS(I,N,Q) = T_DELT_DISORDS(I,N) * L_TDEL
            ENDDO
          ENDDO
        ENDDO

!  Partial layers

        DO UT = 1, N_PARTLAYERS
          XT  = PARTAU_VERT(UT)
          N   = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TD = - L_DELTAU_VERT(Q,N) * XT
            L_TU = - L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - XT )
            DO I = 1, NSTREAMS
              LDN = L_TD / QUAD_STREAMS(I)
              LUP = L_TU / QUAD_STREAMS(I)
              L_T_DISORDS_UTDN(I,UT,Q) = T_DISORDS_UTDN(I,UT) * LDN
              L_T_DISORDS_UTUP(I,UT,Q) = T_DISORDS_UTUP(I,UT) * LUP
            ENDDO
          ENDDO
        ENDDO

!  end solution saving option

      ENDIF

!  Linearization of Transmittance factors for User Streams
!  =======================================================

!  If no user streams, then return

      IF ( .NOT. DO_USER_STREAMS  ) RETURN

!  Whole Layer transmittance factors
!  ---------------------------------

      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS = T_DELT_USERM(N,UM) * USER_SECANTS(UM) * DELTAU_VERT(N)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_USERM(N,UM,Q) = - TRANS * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Partial Layer transmittance factors for off-grid optical depths
!  ---------------------------------------------------------------

      DO UT = 1, N_PARTLAYERS
        N = PARTLAYERS_LAYERIDX(UT)
        UX = PARTAU_VERT(UT)
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS_D = T_UTDN_USERM(UT,UM) * USER_SECANTS(UM)
            TRANS_U = T_UTUP_USERM(UT,UM) * USER_SECANTS(UM)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              VD = L_DELTAU_VERT(Q,N) * UX
              VU = L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - UX )
              L_T_UTDN_USERM(UT,UM,Q) = - TRANS_D * VD
              L_T_UTUP_USERM(UT,UM,Q) = - TRANS_U * VU
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LA_PREPTRANS

!  End module

      END MODULE vlidort_la_miscsetups_m

