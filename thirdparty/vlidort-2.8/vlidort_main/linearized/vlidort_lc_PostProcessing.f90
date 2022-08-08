
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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #                                                        #
! #            LC_WHOLELAYER_STERM_UP                      #
! #            LC_WHOLELAYER_STERM_DN                      #
! #            LC_PARTLAYER_STERM_UP                       #
! #            LC_PARTLAYER_STERM_DN                       #
! #                                                        #
! #            LC_QUAD_GFUNCMULT  (private)                #
! #            QUADCOLUMNWF_LEVEL_UP                       #
! #            QUADCOLUMNWF_LEVEL_DN                       #
! #            QUADCOLUMNWF_OFFGRID_UP                     #
! #            QUADCOLUMNWF_OFFGRID_DN                     #
! #                                                        #
! ##########################################################

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> New input arguments to the post-processing subroutines (Green's function) : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, SIGMA_P, SIGMA_M

!    ==> Other Input/Output arguments and other changes
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!           ** Sourceterm and quadrature arrays now moved to PostProcessing Module (here!)
!           ** Use post-processing mask for Observational/Lattice-Doublet distinction

      MODULE vlidort_lc_PostProcessing_m

!  1/31/21. Version 2.8.3. Taylor series routines needed here.

      USE vlidort_Taylor_m

!  1/31/21. Version 2.8.3. post-processing routines all public

      PUBLIC  :: LC_WHOLELAYER_STERM_UP , LC_WHOLELAYER_STERM_DN,  &
                 LC_PARTLAYER_STERM_UP  , LC_PARTLAYER_STERM_DN,   &
                 QUADCOLUMNWF_LEVEL_UP  , QUADCOLUMNWF_LEVEL_DN,   &
                 QUADCOLUMNWF_OFFGRID_UP, QUADCOLUMNWF_OFFGRID_DN, &
                 LC_QUAD_GFUNCMULT

      CONTAINS

!

      SUBROUTINE LC_WHOLELAYER_STERM_UP ( &
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                   & ! Input flags
             DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SOURCETERM_FLAG, FOURIER,              & ! Input flags/Fourier
             LAYER, IBEAM, NV_PARAMETERS, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK,                & ! Input numbers
             TAYLOR_ORDER, DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM,            & ! input Optical + User_streams
             BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, EMULT_UP, SIGMA_P,                     & ! input Beam stuff
             K_REAL, K_COMPLEX, LCON, MCON, UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2, UPAR_UP_1,& ! input Homogeneous/Beam-1
             ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UPAR_UP_2, PMULT_UU, PMULT_UD,         & ! Input Greens/Beam-2
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_EMULT_UP,               & ! input Linearized Beam stuff
             L_KEIGEN, NCON, PCON, L_UHOM_UPDN, L_UHOM_UPUP, L_HMULT_1, L_HMULT_2,            & ! input Linearized Homogeneous
             L_UPAR_UP_1, LC_UPAR_UP_2, L_ATERM_SAVE, L_BTERM_SAVE, L_LAYER_TSUP_UP,          & ! Input linearized PIs
             L_LAYERSOURCE )                                                                    ! OUTPUT

!  1/31/21. Version 2.8.3. Several Changes.
!    -- Introduce control flag DO_CLASSICAL_SOLUTION, and Taylor-series ordering
!    -- Introduce Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK
!    -- Introduce Green's function inputs as indicated below
!    -- Subroutine has completely new section for the Green's function treatment

!  1/31/21. Version 2.8.3. Add TAYLOR_SMALL to the parameter list

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  Inputs
!  ======

!  Control and Bookkeeping
!  -----------------------

!  flags

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES

      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG

!  1/31/21. Version 2.8.3.  Introduce control for RTE solution

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Numbers

      INTEGER, INTENT (IN) ::          IBEAM, LAYER, FOURIER
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES

!  1/31/21. Version 2.8.3.  Introduce User-stream post-processing masking (as in LIDORT0

      INTEGER, INTENT (IN) ::          N_PPSTREAMS
      INTEGER, INTENT (IN) ::          PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Several new inputs
!    -- Order of Taylor series (including terms up to EPS^n)
!    -- Local vertical optical depth and its linearization
!    -- Transmittance factors for user-defined stream angles, user secants

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER
      DOUBLE PRECISION, intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: T_DELT_USERM  ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, intent(in)  :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Main solution inputs
!  --------------------

!  1/31/21. Version 2.8.3.  Beam average-secant parameterization
!    BEAM_CUTOFF = Last layer to include Particular integral solution

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )

!  Homog Solutions bookkeeping

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  integration constamts

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Homog User solutions and multipliers

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Particular integral user solutions and beam multipliers
!    -- 1/31/21. Version 2.8.3. SIGMA_P is new addition
!  3/1/21. Bug Fix on SIGMA_P declaration

      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_P   ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. New Input: Green's function solution variables

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  1/31/21. Version 2.8.3. New Input: Green's function Multipliers

      DOUBLE PRECISION, INTENT (IN) :: PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized inputs
!  -----------------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     -- LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!     -- 1/31/21. Version 2.8.3. New additions for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) ::  LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized eigenvalues
!     -- 1/31/21. Version 2.8.3. New addition for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous solutionss and multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Particular integral inputs (Beam/Thermal)

      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_EMULT_UP( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3.  New Input: Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, M, UM, O1, Q, IB
      INTEGER ::          K, KO1, K0, K1, K2, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SPAR, H1, H2, TM, TUDA, TUUB
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!  1/31/21. Version 2.8.3. New Help variables for the Green's function solution

      DOUBLE PRECISION :: ITRANS, WDEL, WUDEL, ITRANSWDEL, MULTDN, MULTUP, EPS, DELTA, YFAC, IGAM, SM
      DOUBLE PRECISION :: L_FIRST, L_DELTA, L_KEG, L_LAM, L_GAMMA, L_WDEL, ESUM, LC_ESUM, L_TUDA, L_TUUB
      DOUBLE PRECISION :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)
      DOUBLE PRECISION :: SPARG(MAXSTOKES,MAX_ATMOSWFS)

!  Proxies

      N   = LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM
      M   = FOURIER        ! Only for debugging

!   Very important to zero both output terms (bug solved 12/29/05)
!    -- 1/31/21. Version 2.8.3. Note use of post-processing mask (et seq. in this and other subroutines)

      DO Q = 1, NV_PARAMETERS
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)  = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and..NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789
!      Version 2.8, GOTO statement removed, 19 July 2016

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Homogeneous solutions
!  =====================

!  Loop over user angles and Stokes
!    -- 1/31/21. Version 2.8.3. Use postprocessing mask

        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES

!  parameter loop

            DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LUXR  = LCON(K,N)   *   UHOM_UPDN(UM,O1,K,N)
                MUXR  = MCON(K,N)   *   UHOM_UPUP(UM,O1,K,N)
                LLUXR = LCON(K,N)   * L_UHOM_UPDN(UM,O1,K,N,Q)
                MLUXR = MCON(K,N)   * L_UHOM_UPUP(UM,O1,K,N,Q)
                NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
                PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
                H1 = ( NUXR + LLUXR ) *   HMULT_2(K,UM,N) &
                           +  LUXR    * L_HMULT_2(K,UM,N,Q)
                H2 = ( PUXR + MLUXR ) *   HMULT_1(K,UM,N) &
                           +  MUXR    * L_HMULT_1(K,UM,N,Q)
                SHOM_R = SHOM_R + H1 + H2
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

                LUXR1 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K1,N) - &
                         LCON(K2,N)   *   UHOM_UPDN(UM,O1,K2,N)
                LUXR2 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K2,N) + &
                         LCON(K2,N)   *   UHOM_UPDN(UM,O1,K1,N)
                MUXR1 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K1,N) - &
                         MCON(K2,N)   *   UHOM_UPUP(UM,O1,K2,N)
                MUXR2 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K2,N) + &
                         MCON(K2,N)   *   UHOM_UPUP(UM,O1,K1,N)

                LLUXR1 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q) - &
                         LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q)
                LLUXR2 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q) + &
                         LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q)
                MLUXR1 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q) - &
                         MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q)
                MLUXR2 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q) + &
                         MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q)

                NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) - &
                         NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
                NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) + &
                         NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
                PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) - &
                         PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
                PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) + &
                         PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

                H1 =    ( NUXR1 + LLUXR1 ) * HMULT_2(K1,UM,N) &
                      - ( NUXR2 + LLUXR2 ) * HMULT_2(K2,UM,N) &
                            +  LUXR1       * L_HMULT_2(K1,UM,N,Q) &
                            -  LUXR2       * L_HMULT_2(K2,UM,N,Q)
                H2 =    ( PUXR1 + MLUXR1 ) * HMULT_1(K1,UM,N) &
                      - ( PUXR2 + MLUXR2 ) * HMULT_1(K2,UM,N) &
                            +  MUXR1       * L_HMULT_1(K1,UM,N,Q) &
                            -  MUXR2       * L_HMULT_1(K2,UM,N,Q)

                SHOM_CR = SHOM_CR + H1 + H2

              ENDDO

!  homogeneous contribution

              L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR
            ENDDO
          ENDDO
        ENDDO

!  End scattering-only solutions

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1 if solar sources are included (taken care of earlier)
!     ----- Linearization always exists

!  1/31/21. Version 2.8.3. Note use of post-processing mask

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1 ; TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          DO Q = 1, NV_PARAMETERS
            L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + L_LAYER_TSUP_UP(UM,N,Q)*TM
          ENDDO
        ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.
!    -- 1/31/21. Version 2.8.3. This condition added

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Particular and single scatter contributions - CLASSICAL SOLUTION
!  ================================================================

!  add the MS particular solution, and if flagged, the single scatter term

!  1/31/21. Version 2.8.3. Note use of post-processing mask and classical solution contrel
!     -- No more separate observational geometry clause

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          ESUM = EMULT_UP(LUM,N,IB)
          DO Q = 1, NV_PARAMETERS
            LC_ESUM = LC_EMULT_UP(LUM,N,IB,Q)
            DO O1 = 1, NSTOKES
              SPAR = LC_UPAR_UP_2(UM,O1,N,Q) * ESUM + UPAR_UP_2(UM,O1,N) * LC_ESUM
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Particular and single scatter contributions - GREENS FUNCTION SOLUTION
!  ======================================================================

!  1/31/21. Version 2.8.3. Completely New section
!     ==> Lattice/Observational geometry, controlled by masking system
!     ==> Multiplier calculations as for LIDORT, real eigenvalues only

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Some beam particulars
!mick mod 1/5/2021 - applied neg. sign to formulations for ITRANSWDEL and L_ITRANSWDEL
!                    for consistency with WHOLELAYER_STERM_UP

         WDEL       = T_DELT_MUBAR(N,IB)
         ITRANS     = INITIAL_TRANS(N,IB)
         ITRANSWDEL = - ITRANS * WDEL
         DO Q = 1, NV_PARAMETERS
            L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
            L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
            L_ITRANSWDEL(Q) = - (ITRANS * L_WDEL + WDEL * L_ITRANS(Q))
         ENDDO

!  Start local user angle loop and initialize

         DO LUM = 1, N_PPSTREAMS
            UM    = PPSTREAM_MASK(LUM,IB)
            SM    = USER_SECANTS(UM)
            SPARG = ZERO

!  Start eigenvalue loop

            DO K = 1, K_REAL(N)

!  Downwelling multipliers

               if ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
                  EPS   = GAMMA_M(K,N)              ; DELTA = DELTAU_VERT(N)
                  WUDEL = WDEL * T_DELT_USERM(N,UM) ; YFAC  = SIGMA_P(N,LUM,IB)
                  CALL Taylor_Series_2 ( TAYLOR_ORDER, EPS, YFAC, DELTA, ONE, WUDEL, SM, MULTDN )
                  DO Q = 1, NV_PARAMETERS
                     L_KEG   = L_KEIGEN(K,N,Q)
                     L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                     L_DELTA = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N) ! Input is single normalized 
                     CALL Taylor_Series_L_2a &
                      ( TAYLOR_ORDER, EPS, YFAC, DELTA, L_DELTA, L_LAM, L_KEG, WUDEL, SM, L_MULTDN(Q) )
                       L_MULTDN(Q) = ITRANS * L_MULTDN(Q) + MULTDN * L_ITRANS(Q)
                  ENDDO
                  MULTDN = MULTDN * ITRANS
               else

!   Rob Change 1/15/21. replace MULTDN to avoid division by ATERM_SAVE
!                  MULTDN = PMULT_UD(K,UM,N) / ATERM_SAVE(K,N) ; IGAM = ONE / GAMMA_M(K,N)

                  MULTDN = PMULT_UD(K,UM,N) ; IGAM = ONE / GAMMA_M(K,N)
                  DO Q = 1, NV_PARAMETERS
                     L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(K,N,Q)
                     L_FIRST  = ITRANS * L_HMULT_2(K,UM,N,Q) + L_ITRANS(Q) * HMULT_2(K,UM,N)
                     L_MULTDN(Q) = ( L_FIRST - LC_EMULT_UP(LUM,N,IB,Q) - L_GAMMA * MULTDN ) * IGAM
                  ENDDO
               endif

!  Upwelling multipliers
!   Rob Change 1/15/21. replace MULTUP to avoid division by BTERM_SAVE
!               MULTUP = PMULT_UU(K,UM,N) / BTERM_SAVE(K,N) ; IGAM = ONE / GAMMA_P(K,N)
!mick mod 1/5/2021 - removed neg. sign from L_FIRST for consistency with WHOLELAYER_STERM_UP

               MULTUP = PMULT_UU(K,UM,N) ; IGAM = ONE / GAMMA_P(K,N)
               DO Q = 1, NV_PARAMETERS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(K,N,Q)
                  L_FIRST  = ITRANSWDEL * L_HMULT_1(K,UM,N,Q) + L_ITRANSWDEL(Q) * HMULT_1(K,UM,N)
                  L_MULTUP(Q) = ( L_FIRST + LC_EMULT_UP(LUM,N,IB,Q) - L_GAMMA * MULTUP ) * IGAM
               ENDDO

!  Complete the multipliers
!   Rob Change 1/15/21. replace this do-loop using regular derivatives (non-logarithmic)
!mick fix 1/5/2021 - changed defs of L_PMULT_DN and L_PMULT_UP here to correspond with
!                    changed defs of PMULT_UD and PMULT_UU in WHOLELAYER_STERM_UP

!               DO Q = 1, NV_PARAMETERS
!                  L_PMULT_DN(Q) = ATERM_SAVE(K,N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(K,N,Q) )
!                  L_PMULT_UP(Q) = BTERM_SAVE(K,N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(K,N,Q) )
!               ENDDO
!               DO Q = 1, NV_PARAMETERS
!                  L_PMULT_DN(Q) = ATERM_SAVE(K,N) * L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(K,N,Q)
!                  L_PMULT_UP(Q) = BTERM_SAVE(K,N) * L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(K,N,Q)
!               ENDDO
               DO Q = 1, NV_PARAMETERS
                  L_PMULT_DN(Q) = L_MULTDN(Q)
                  L_PMULT_UP(Q) = L_MULTUP(Q)
               ENDDO

!  Add the particular solution contributions
!    Rob Change 1/15/21. Now have to multiply second terms by ATERM and BTERM. Replacement follows (efficiently)
!mick fix 1/5/2021 - changed formulations for H1 and H2 due to changed defs of L_PMULT_DN and L_PMULT_UP

!               DO Q = 1, NV_PARAMETERS ; DO O1 = 1, NSTOKES
!                  H1 = UHOM_UPDN(UM,O1,K,N) * L_PMULT_DN(Q) + L_UHOM_UPDN(UM,O1,K,N,Q) * PMULT_UD(K,UM,N)
!                  H2 = UHOM_UPUP(UM,O1,K,N) * L_PMULT_UP(Q) + L_UHOM_UPUP(UM,O1,K,N,Q) * PMULT_UU(K,UM,N)
!                  SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
!               ENDDO ; ENDDO

               TUDA = PMULT_UD(K,UM,N) * ATERM_SAVE(K,N)
               TUUB = PMULT_UU(K,UM,N) * BTERM_SAVE(K,N)
!               DO Q = 1, NV_PARAMETERS ; DO O1 = 1, NSTOKES
!                  H1 = UHOM_UPDN(UM,O1,K,N) * L_PMULT_DN(Q) + L_UHOM_UPDN(UM,O1,K,N,Q) * TUDA
!                  H2 = UHOM_UPUP(UM,O1,K,N) * L_PMULT_UP(Q) + L_UHOM_UPUP(UM,O1,K,N,Q) * TUUB
!                  SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
!               ENDDO ; ENDDO
               DO Q = 1, NV_PARAMETERS
                  L_TUDA = PMULT_UD(K,UM,N) * L_ATERM_SAVE(K,N,Q) + L_PMULT_DN(Q) * ATERM_SAVE(K,N)
                  L_TUUB = PMULT_UU(K,UM,N) * L_BTERM_SAVE(K,N,Q) + L_PMULT_UP(Q) * BTERM_SAVE(K,N)
                  DO O1 = 1, NSTOKES
                    H1 = UHOM_UPDN(UM,O1,K,N) * L_TUDA + L_UHOM_UPDN(UM,O1,K,N,Q) * TUDA
                    H2 = UHOM_UPUP(UM,O1,K,N) * L_TUUB + L_UHOM_UPUP(UM,O1,K,N,Q) * TUUB
                    SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
                  ENDDO
               ENDDO

!  End eigenvalue loop

            ENDDO

!  Add Green's function result to the total

            DO Q = 1, NV_PARAMETERS
               DO O1 = 1, NSTOKES
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPARG(O1,Q)
               ENDDO
            ENDDO

!  End user stream loop

         ENDDO

!  End Green's function calculation

      ENDIF

!  If NOT operating in MS-mode only, add single scatter part
!    -- 1/31/21. Version 2.8.3. Use of post-processing mask

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          ESUM = EMULT_UP(LUM,N,IB)
          DO Q = 1, NV_PARAMETERS
            LC_ESUM = LC_EMULT_UP(LUM,N,IB,Q)
            DO O1 = 1, NSTOKES
              SPAR = L_UPAR_UP_1(UM,O1,N,Q) * ESUM + UPAR_UP_1(UM,O1,N) * LC_ESUM
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_STERM_UP

!

      SUBROUTINE LC_WHOLELAYER_STERM_DN ( &
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                   & ! Input flags
             DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SOURCETERM_FLAG, FOURIER,              & ! Input flags/Fourier
             LAYER, IBEAM, NV_PARAMETERS, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK,                & ! Input numbers
             TAYLOR_ORDER, DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM,            & ! input Optical + User_streams
             BEAM_CUTOFF, AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR, EMULT_DN, SIGMA_M,     & ! input Beam stuff
             K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, HMULT_1, HMULT_2, UPAR_DN_1,& ! input Homogeneous/Beam-1
             ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UPAR_DN_2, PMULT_DU, PMULT_DD,         & ! Input Greens/Beam-2
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_EMULT_DN,               & ! input Linearized Beam stuff
             L_KEIGEN, NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_HMULT_1, L_HMULT_2,            & ! input Linearized Homogeneous
             L_UPAR_DN_1, LC_UPAR_DN_2, L_ATERM_SAVE, L_BTERM_SAVE, L_LAYER_TSUP_DN,          & ! Input linearized PIs
             L_LAYERSOURCE )                                                                    ! OUTPUT

!  1/31/21. Version 2.8.3. Several Changes.
!    -- Introduce control flag DO_CLASSICAL_SOLUTION, and Taylor-series ordering
!    -- Introduce Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK
!    -- Introduce Green's function inputs as indicated below
!    -- Subroutine has completely new section for the Green's function treatment

!  1/31/21. Version 2.8.3. Add TAYLOR_SMALL to the parameter list

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4, TAYLOR_SMALL, LDU

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG

!  1/31/21. Version 2.8.3.  Introduce control for RTE solution

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Numbers

      INTEGER, INTENT (IN) ::          IBEAM, LAYER, FOURIER
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES

!  1/31/21. Version 2.8.3.  Introduce User-stream post-processing masking (as in LIDORT0

      INTEGER, INTENT (IN) ::          N_PPSTREAMS
      INTEGER, INTENT (IN) ::          PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Several new inputs
!    -- Order of Taylor series (including terms up to EPS^n)
!    -- Local vertical optical depth and its linearization
!    -- Transmittance factors for user-defined stream angles, user secants

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER
      DOUBLE PRECISION, intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: T_DELT_USERM  ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, intent(in)  :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Main solution inputs
!  --------------------

!  1/31/21. Version 2.8.3.  Beam average-secant parameterization
!    BEAM_CUTOFF = Last layer to include Particular integral solution

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
!  Homog solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  User solutions

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!    -- 1/31/21. Version 2.8.3. SIGMA_P is new addition
!  3/1/21. Bug Fix on SIGMA_M declaration

      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_M  ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. New Input: Green's function solution variables

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  1/31/21. Version 2.8.3. New Input: Green's function Multipliers

      DOUBLE PRECISION, INTENT (IN) :: PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized inputs
!  -----------------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     -- LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!     -- 1/31/21. Version 2.8.3. New additions for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) ::  LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized eigenvalues
!     -- 1/31/21. Version 2.8.3. New addition for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous solutionss and multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Particular integral inputs (Beam/Thermal)

      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_DN  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_EMULT_DN( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3.  New Input: Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, M, UM, O1, Q, IB
      INTEGER ::          K, KO1, K0, K1, K2, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SPAR, H1, H2, TM, TDDA, TDUB
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!  1/31/21. Version 2.8.3. New Help variables for the Green's function solution

      DOUBLE PRECISION :: ITRANS, WDEL, ITRANSWDEL, MULTDN, MULTUP, EPS, DELTA, YFAC, IGAM, SM
      DOUBLE PRECISION :: L_FIRST, L_DELTA, L_KEG, L_LAM, L_GAMMA, L_WDEL, ESUM, LC_ESUM, LAM, UDEL, &
                          L_TDDA, L_TDUB
      DOUBLE PRECISION :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)
      DOUBLE PRECISION :: SPARG(MAXSTOKES,MAX_ATMOSWFS)

!  Proxies

      N   = LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM
      M   = FOURIER   ! Only for debugging

!   Very important to zero both output terms (bug solved 12/29/05)
!    -- 1/31/21. Version 2.8.3. Note use of post-processing mask (et seq. in this and other subroutines)

      DO Q = 1, NV_PARAMETERS
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)  = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789
!      Version 2.8, GOTO statement removed, 19 July 2016

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Homogeneous solutions
!  =====================

!  Loop over user angles and Stokes

        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES

!  parameter loop

            DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LUXR  = LCON(K,N)   *   UHOM_DNDN(UM,O1,K,N)
                MUXR  = MCON(K,N)   *   UHOM_DNUP(UM,O1,K,N)
                LLUXR = LCON(K,N)   * L_UHOM_DNDN(UM,O1,K,N,Q)
                MLUXR = MCON(K,N)   * L_UHOM_DNUP(UM,O1,K,N,Q)
                NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
                PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)
                H1 = ( NUXR + LLUXR ) *   HMULT_1(K,UM,N) &
                           +  LUXR    * L_HMULT_1(K,UM,N,Q)
                H2 = ( PUXR + MLUXR ) *   HMULT_2(K,UM,N) &
                           +  MUXR    * L_HMULT_2(K,UM,N,Q)
                SHOM_R = SHOM_R + H1 + H2
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

                LUXR1 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K1,N) - &
                         LCON(K2,N)   *   UHOM_DNDN(UM,O1,K2,N)
                LUXR2 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K2,N) + &
                         LCON(K2,N)   *   UHOM_DNDN(UM,O1,K1,N)
                MUXR1 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K1,N) - &
                         MCON(K2,N)   *   UHOM_DNUP(UM,O1,K2,N)
                MUXR2 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K2,N) + &
                         MCON(K2,N)   *   UHOM_DNUP(UM,O1,K1,N)

                LLUXR1 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q) - &
                         LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q)
                LLUXR2 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q) + &
                         LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q)
                MLUXR1 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q) - &
                         MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q)
                MLUXR2 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q) + &
                         MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q)

                NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) - &
                         NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
                NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) + &
                         NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
                PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) - &
                         PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
                PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) + &
                         PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

                H1 =    ( NUXR1 + LLUXR1 ) * HMULT_1(K1,UM,N) &
                      - ( NUXR2 + LLUXR2 ) * HMULT_1(K2,UM,N) &
                            +  LUXR1       * L_HMULT_1(K1,UM,N,Q) &
                            -  LUXR2       * L_HMULT_1(K2,UM,N,Q)
                H2 =    ( PUXR1 + MLUXR1 ) * HMULT_2(K1,UM,N) &
                      - ( PUXR2 + MLUXR2 ) * HMULT_2(K2,UM,N) &
                            +  MUXR1       * L_HMULT_2(K1,UM,N,Q) &
                            -  MUXR2       * L_HMULT_2(K2,UM,N,Q)

                SHOM_CR = SHOM_CR + H1 + H2

              ENDDO

!  homogeneous contribution

              L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
          ENDDO
        ENDDO

!  End scattering clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1 if solar sources are included (taken care of earlier)
!  1/31/21. Version 2.8.3. Note use of post-processing mask

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1 ; TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          DO Q = 1, NV_PARAMETERS
            L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + L_LAYER_TSUP_DN(UM,N,Q)*TM
          ENDDO
        ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.
!    -- 1/31/21. Version 2.8.3. This condition added

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Particular and single scatter contributions - CLASSICAL SOLUTION
!  ================================================================

!  add the MS particular solution
!  1/31/21. Version 2.8.3. Note use of post-processing mask and classical solution contrel
!     -- No more separate observational geometry clause

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          ESUM = EMULT_DN(LUM,N,IB)
          DO Q = 1, NV_PARAMETERS
            LC_ESUM = LC_EMULT_DN(LUM,N,IB,Q)
            DO O1 = 1, NSTOKES
              SPAR = LC_UPAR_DN_2(UM,O1,N,Q) * ESUM + UPAR_DN_2(UM,O1,N) * LC_ESUM
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Particular and single scatter contributions - GREENS FUNCTION SOLUTION
!  ======================================================================

!  add the MS particular solution

!  1/31/21. Version 2.8.3. Completely New section
!     ==> Lattice/Observational geometry, controlled by masking system
!     ==> Multiplier calculations as for LIDORT, real eigenvalues only

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Some beam particulars
!mick fix 1/5/2021 - applied neg. sign to formulations for ITRANSWDEL and L_ITRANSWDEL
!                    for consistency with WHOLELAYER_STERM_DN

        WDEL       = T_DELT_MUBAR(N,IB)
        ITRANS     = INITIAL_TRANS(N,IB)
        ITRANSWDEL = - ITRANS * WDEL
        DO Q = 1, NV_PARAMETERS
          L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
          L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
          L_ITRANSWDEL(Q) = - (ITRANS * L_WDEL + WDEL * L_ITRANS(Q))
        ENDDO

!  start local user angle loop and initialize

        DO LUM = 1, N_PPSTREAMS
          UM    = PPSTREAM_MASK(LUM,IB)
          SM    = USER_SECANTS(UM)
          SPARG = ZERO

!  Start eigenvalue loop

         DO K = 1, K_REAL(N)

!  Downwelling multipliers. Note use of -YFAC in L_2b. Rob 01/09/14

            if ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(K,N)         ;  DELTA = DELTAU_VERT(N)
               UDEL  = T_DELT_USERM(N,UM)   ;  YFAC  = SIGMA_M(N,LUM,IB) ; LAM = AVERAGE_SECANT(N,IB)
               CALL Taylor_Series_2 ( TAYLOR_ORDER, EPS, YFAC, DELTA, UDEL, WDEL, SM, MULTDN )
               DO Q = 1, NV_PARAMETERS
                  L_KEG   = L_KEIGEN(K,N,Q)
                  L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                  L_DELTA = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N) ! Input is single normalized 
                  CALL Taylor_Series_L_2b &
                   ( TAYLOR_ORDER, EPS, -YFAC, DELTA, L_DELTA, L_LAM, L_KEG, WDEL, UDEL, SM, LAM, L_MULTDN(Q) )
                  L_MULTDN(Q) = L_MULTDN(Q) * ITRANS + MULTDN * L_ITRANS(Q)
               ENDDO
               MULTDN = MULTDN * ITRANS
            else

!   Rob Change 1/15/21. replace MULTDN to avoid division by ATERM_SAVE
!               MULTDN = PMULT_DD(K,UM,N) / ATERM_SAVE(K,N) ; IGAM = ONE / GAMMA_M(K,N)

               MULTDN = PMULT_DD(K,UM,N) ; IGAM = ONE / GAMMA_M(K,N)
               DO Q = 1, NV_PARAMETERS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(K,N,Q)
                  L_FIRST  = ITRANS * L_HMULT_1(K,UM,N,Q) + L_ITRANS(Q) * HMULT_1(K,UM,N)
                  L_MULTDN(Q) = ( L_FIRST - LC_EMULT_DN(LUM,N,IB,Q) - L_GAMMA * MULTDN ) * IGAM
               ENDDO
            endif

!  Upwelling multipliers
!   Rob Change 1/15/21. replace MULTUP to avoid division by BTERM_SAVE
!            MULTUP = PMULT_DU(K,UM,N) / BTERM_SAVE(K,N) ; IGAM = ONE / GAMMA_P(K,N)
 
            MULTUP = PMULT_DU(K,UM,N) ; IGAM = ONE / GAMMA_P(K,N)
            DO Q = 1, NV_PARAMETERS
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(K,N,Q)
               L_FIRST  = ITRANSWDEL * L_HMULT_2(K,UM,N,Q) + L_ITRANSWDEL(Q) * HMULT_2(K,UM,N)
               L_MULTUP(Q) = ( L_FIRST + LC_EMULT_DN(LUM,N,IB,Q) - L_GAMMA * MULTUP ) * IGAM
            ENDDO

!  Complete the multipliers
!   Rob Change 1/15/21. replace this do-loop using regular derivatives (non-logarithmic)
!mick fix 1/5/2021 - changed defs of L_PMULT_DN and L_PMULT_UP here to correspond with
!                    changed defs of PMULT_DD and PMULT_DU in WHOLELAYER_STERM_DN

!            DO Q = 1, NV_PARAMETERS
!               L_PMULT_DN(Q) = ATERM_SAVE(K,N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(K,N,Q) )
!               L_PMULT_UP(Q) = BTERM_SAVE(K,N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(K,N,Q) )
!            ENDDO
!            DO Q = 1, NV_PARAMETERS
!               L_PMULT_DN(Q) = ATERM_SAVE(K,N) * L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(K,N,Q)
!               L_PMULT_UP(Q) = BTERM_SAVE(K,N) * L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(K,N,Q)
!            ENDDO
            DO Q = 1, NV_PARAMETERS
               L_PMULT_DN(Q) = L_MULTDN(Q)
               L_PMULT_UP(Q) = L_MULTUP(Q)
            ENDDO

!  Add the particular solution contributions
!    Rob Change 1/15/21. Now have to multiply second terms by ATERM and BTERM. Replacement follows (efficiently)
!mick fix 1/5/2021 - changed formulations for H1 and H2 due to changed defs of L_PMULT_DN and L_PMULT_UP

!            DO Q = 1, NV_PARAMETERS ; DO O1 = 1, NSTOKES
!              H1 = UHOM_DNDN(UM,O1,K,N) * L_PMULT_DN(Q) + L_UHOM_DNDN(UM,O1,K,N,Q) * PMULT_DD(K,UM,N)
!              H2 = UHOM_DNUP(UM,O1,K,N) * L_PMULT_UP(Q) + L_UHOM_DNUP(UM,O1,K,N,Q) * PMULT_DU(K,UM,N)
!              SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
!            ENDDO ; ENDDO

            TDDA = PMULT_DD(K,UM,N) * ATERM_SAVE(K,N)
            TDUB = PMULT_DU(K,UM,N) * BTERM_SAVE(K,N)
!            DO Q = 1, NV_PARAMETERS ; DO O1 = 1, NSTOKES
!              H1 = UHOM_DNDN(UM,O1,K,N) * L_PMULT_DN(Q) + L_UHOM_DNDN(UM,O1,K,N,Q) * TDDA
!              H2 = UHOM_DNUP(UM,O1,K,N) * L_PMULT_UP(Q) + L_UHOM_DNUP(UM,O1,K,N,Q) * TDUB
!              SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
!            ENDDO ; ENDDO
            DO Q = 1, NV_PARAMETERS
              L_TDDA = PMULT_DD(K,UM,N) * L_ATERM_SAVE(K,N,Q) + L_PMULT_DN(Q) * ATERM_SAVE(K,N) 
              L_TDUB = PMULT_DU(K,UM,N) * L_BTERM_SAVE(K,N,Q) + L_PMULT_UP(Q) * BTERM_SAVE(K,N)
              DO O1 = 1, NSTOKES
                H1 = UHOM_DNDN(UM,O1,K,N) * L_TDDA + L_UHOM_DNDN(UM,O1,K,N,Q) * TDDA
                H2 = UHOM_DNUP(UM,O1,K,N) * L_TDUB + L_UHOM_DNUP(UM,O1,K,N,Q) * TDUB
                SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
              ENDDO
            ENDDO

!  End eigenvalue loop

          ENDDO

!  Add Green's function result to the total

          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPARG(O1,Q)
            ENDDO
          ENDDO

!  End user stream loop

        ENDDO

!  End Green's function calculation

      ENDIF

!  Add single scatter term if flagged
!    -- 1/31/21. Version 2.8.3. Note use of post-processing mask

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
               DO Q = 1, NV_PARAMETERS
                  SPAR = L_UPAR_DN_1(UM,O1,N,Q) *    EMULT_DN(LUM,N,IB) + &
                           UPAR_DN_1(UM,O1,N)   * LC_EMULT_DN(LUM,N,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_STERM_DN

!

      SUBROUTINE LC_PARTLAYER_STERM_UP ( &
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                & ! Input flags
           DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SOURCETERM_FLAG, FOURIER,           & ! Input flags/numbers
           IB, N, UT, NV_PARAMETERS, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,  & ! Input numbers
           PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTUP_USERM,          & ! Input Optical/Usertrans
           BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                       & ! input Beam stuff
           K_REAL, K_COMPLEX, LCON, MCON, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1,               & ! Input user-solutions
           UPAR_UP_2, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, SIGMA_P,                    & ! Input User-solutions/Greens
           ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UT_PMULT_UU, UT_PMULT_UD,           & ! Input Greens
           LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_UT_EMULT_UP,         & ! input Linearized Beam stuff
           L_KEIGEN, NCON, PCON, L_UHOM_UPDN, L_UHOM_UPUP, L_UT_HMULT_UU, L_UT_HMULT_UD, & ! input Linearized Homogeneous
           L_UPAR_UP_1, LC_UPAR_UP_2, L_ATERM_SAVE, L_BTERM_SAVE, L_LAYER_TSUP_UTUP,     & ! Input linearized PIs
           L_LAYERSOURCE )                                                                  ! OUTPUT

!  1/31/21. Version 2.8.3. Several Changes.
!    -- Introduce control flag DO_CLASSICAL_SOLUTION, and Taylor-series ordering
!    -- Introduce Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK
!    -- Introduce Green's function inputs as indicated below
!    -- Subroutine has completely new section for the Green's function treatment

!  1/31/21. Version 2.8.3. Add TAYLOR_SMALL to the parameter list

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  Inputs
!  ======

!  Control and Bookkeeping
!  -----------------------

!  flags

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES

      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG

!  1/31/21. Version 2.8.3.  Introduce control for RTE solution

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Numbers

      INTEGER, INTENT (IN) ::          N, UT, IB, FOURIER
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES

!  1/31/21. Version 2.8.3.  Introduce User-stream post-processing masking (as in LIDORT0

      INTEGER, INTENT (IN) ::          N_PPSTREAMS
      INTEGER, INTENT (IN) ::          PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Several new inputs
!    -- Order of Taylor series (including terms up to EPS^n)
!    -- Local vertical optical depth and its linearization
!    -- Transmittance factors for user-defined stream angles, user secants

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER
      DOUBLE PRECISION, INTENT(IN)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN)  :: T_UTUP_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, intent(in)  :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Main solution inputs
!  --------------------

!  1/31/21. Version 2.8.3.  Beam average-secant parameterization
!    BEAM_CUTOFF = Last layer to include Particular integral solution

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )

!  Homog Solutions bookkeeping

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  integration constamts

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Homog User solutions and multipliers

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Particular integral user solutions and beam multipliers
!    -- 1/31/21. Version 2.8.3. SIGMA_P is new addition
!  3/1/21. Bug Fix on SIGMA_P declaration

      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_1   ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_2   ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_P     ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. New Input: Green's function solution variables

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  1/31/21. Version 2.8.3. New Input: Green's function Multipliers

      DOUBLE PRECISION, INTENT (IN) :: UT_PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: UT_PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized inputs
!  -----------------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     -- LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!     -- 1/31/21. Version 2.8.3. New additions for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) ::  LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized eigenvalues
!     -- 1/31/21. Version 2.8.3. New addition for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous solutionss and multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UU  ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UD  ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Particular integral inputs (Beam/Thermal)

      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTUP   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3.  New Input: Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION, INTENT (IN) :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Output
!  ------

!  Linearized output

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER          :: M, UM, O1, Q, K, KO1, K0, K1, K2, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, TM, TUDA, TUUB
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!  1/31/21. Version 2.8.3. New Help variables for the Green's function solution

      DOUBLE PRECISION :: ITRANS, WDEL, ITRANSWDEL, EPS, DELTA, PARTA, YFAC, IGAM, UXUP, LAM
      DOUBLE PRECISION :: SM, MULTDN, MULTUP, FAC1, FAC2, MULT1, MULT2, SPAR, ESUM, LC_ESUM
      DOUBLE PRECISION :: L_FIRST, L_PARTA, L_DELTA, L_KEG, L_LAM, L_GAMMA, L_WDEL, &
                          L_MULT1, L_MULT2, L_UXUP, L_MULT, L_TUDA, L_TUUB
      DOUBLE PRECISION :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)
      DOUBLE PRECISION :: SPARG(MAXSTOKES,MAX_ATMOSWFS)

!  Proxies

      KO1 = K_REAL(N) + 1
      M = FOURIER    ! debug only

!  Very important to zero both output terms (bug solved 12/29/05)
!    -- 1/31/21. Version 2.8.3. Use postprocessing mask.

      DO Q = 1, NV_PARAMETERS
        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only
!  Version 2.8. remove GOTO Statement
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles and Stokes parameters
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES

!  parameter loop

            DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)

                LUXR  = LCON(K,N)   *   UHOM_UPDN(UM,O1,K,N)
                MUXR  = MCON(K,N)   *   UHOM_UPUP(UM,O1,K,N)
                LLUXR = LCON(K,N)   * L_UHOM_UPDN(UM,O1,K,N,Q)
                MLUXR = MCON(K,N)   * L_UHOM_UPUP(UM,O1,K,N,Q)
                NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
                PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
                H1 = ( NUXR + LLUXR ) *   UT_HMULT_UD(K,UM,UT) &
                           +  LUXR    * L_UT_HMULT_UD(K,UM,UT,Q)
                H2 = ( PUXR + MLUXR ) *   UT_HMULT_UU(K,UM,UT) &
                           +  MUXR    * L_UT_HMULT_UU(K,UM,UT,Q)
                SHOM_R = SHOM_R + H1 + H2

              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

                LUXR1 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K1,N) - &
                         LCON(K2,N)   *   UHOM_UPDN(UM,O1,K2,N)
                LUXR2 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K2,N) + &
                         LCON(K2,N)   *   UHOM_UPDN(UM,O1,K1,N)
                MUXR1 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K1,N) - &
                         MCON(K2,N)   *   UHOM_UPUP(UM,O1,K2,N)
                MUXR2 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K2,N) + &
                         MCON(K2,N)   *   UHOM_UPUP(UM,O1,K1,N)

                LLUXR1 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q) - &
                         LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q)
                LLUXR2 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q) + &
                         LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q)
                MLUXR1 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q) - &
                         MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q)
                MLUXR2 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q) + &
                         MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q)

                NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) - &
                         NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
                NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) + &
                         NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
                PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) - &
                         PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
                PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) + &
                         PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

                H1 =    ( NUXR1 + LLUXR1 ) * UT_HMULT_UD(K1,UM,UT) &
                      - ( NUXR2 + LLUXR2 ) * UT_HMULT_UD(K2,UM,UT) &
                            +  LUXR1       * L_UT_HMULT_UD(K1,UM,UT,Q) &
                            -  LUXR2       * L_UT_HMULT_UD(K2,UM,UT,Q)
                H2 =    ( PUXR1 + MLUXR1 ) * UT_HMULT_UU(K1,UM,UT) &
                      - ( PUXR2 + MLUXR2 ) * UT_HMULT_UU(K2,UM,UT) &
                            +  MUXR1       * L_UT_HMULT_UU(K1,UM,UT,Q) &
                            -  MUXR2       * L_UT_HMULT_UU(K2,UM,UT,Q)

                SHOM_CR = SHOM_CR + H1 + H2

              ENDDO

!  homogeneous contribution

              L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
          ENDDO
        ENDDO

!  End scattering clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1.0 if solar sources are included (taken care of earlier)
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1 ; TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO Q = 1, NV_PARAMETERS
            L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + L_LAYER_TSUP_UTUP(UM,UT,Q)*TM
          ENDDO
        ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.
!    -- 1/31/21. Version 2.8.3. Condition added

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Particular and single scatter contributions - CLASSICAL SOLUTION
!  ================================================================

!  add the MS particular solution

!  1/31/21. Version 2.8.3. Note use of post-processing mask and classical solution control
!     -- No more separate observational geometry clause

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          ESUM = UT_EMULT_UP(LUM,UT,IB)
          DO Q = 1, NV_PARAMETERS
            LC_ESUM = LC_UT_EMULT_UP(LUM,UT,IB,Q)
            DO O1 = 1, NSTOKES
              SPAR = LC_UPAR_UP_2(UM,O1,N,Q) * ESUM + UPAR_UP_2(UM,O1,N) * LC_ESUM
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Particular and single scatter contributions - GREENS FUNCTION SOLUTION
!  ======================================================================

!  add the MS particular solution

!  1/31/21. Version 2.8.3. Completely New section
!     ==> Lattice/Observational geometry, controlled by masking system
!     ==> Multiplier calculations as for LIDORT, real eigenvalues only

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Some beam particulars
!mick mod 1/5/2021 - applied neg. sign to formulations for ITRANSWDEL and L_ITRANSWDEL
!                    for consistency with PARTLAYER_STERM_UP

         WDEL       = T_DELT_MUBAR(N,IB)
         ITRANS     = INITIAL_TRANS(N,IB)
         ITRANSWDEL = - ITRANS * WDEL
         DO Q = 1, NV_PARAMETERS
            L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
            L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
            L_ITRANSWDEL(Q) = - (ITRANS * L_WDEL + WDEL * L_ITRANS(Q))
         ENDDO

!  start local user angle loop and initialize

         DO LUM = 1, N_PPSTREAMS
            UM    = PPSTREAM_MASK(LUM,IB)
            SM    = USER_SECANTS(UM)
            SPARG = ZERO

!  Start eigenvalue loop

            DO K = 1, K_REAL(N)

!  Downwelling multipliers
!   * Use of L_2b Taylor series applied twice with UXUP was critical. 01/09/14

               if ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
                  EPS   = GAMMA_M(K,N)      ; FAC1 = T_UTDN_MUBAR(UT,IB) ; UXUP = T_UTUP_USERM(UT,UM)
                  YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL ; LAM = LOG(ONE/WDEL)/DELTA
                  PARTA = PARTAU_VERT(UT)   ; DELTA = DELTAU_VERT(N)
                  CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTA, ZERO, FAC1, SM, MULT1 )
                  CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTA, ZERO, FAC2, SM, MULT2 )
                  MULTDN = - MULT1 + MULT2 * UXUP
                  DO Q = 1, NV_PARAMETERS
                     L_KEG   = L_KEIGEN(K,N,Q)
                     L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                     L_PARTA = L_DELTAU_VERT(Q,N) * PARTAU_VERT(UT) ! Input is single normalized 
                     L_DELTA = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N)  ! Input is single normalized 
                     L_UXUP  = - UXUP * SM * (L_DELTA - L_PARTA)
                     CALL Taylor_Series_L_2b &
                      ( TAYLOR_ORDER, EPS, -YFAC, PARTA, L_PARTA, L_LAM, L_KEG, FAC1, ZERO, SM, LAM, L_MULT1 )
                     CALL Taylor_Series_L_2b &
                      ( TAYLOR_ORDER, EPS, -YFAC, DELTA, L_DELTA, L_LAM, L_KEG, FAC2, ZERO, SM, LAM, L_MULT2 )
                     L_MULT = - L_MULT1 + L_MULT2 * UXUP + MULT2 * L_UXUP
                     L_MULTDN(Q) = L_MULT * ITRANS + MULTDN * L_ITRANS(Q)
                  ENDDO
                  MULTDN = MULTDN * ITRANS
               else

!   Rob Change 1/15/21. replace MULTDN to avoid division by ATERM_SAVE
!                  MULTDN = UT_PMULT_UD(K,UM,UT) / ATERM_SAVE(K,N) ; IGAM = ONE / GAMMA_M(K,N)

                  MULTDN = UT_PMULT_UD(K,UM,UT) ; IGAM = ONE / GAMMA_M(K,N)
                  DO Q = 1, NV_PARAMETERS
                     L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(K,N,Q)
                     L_FIRST  = ITRANS * L_UT_HMULT_UD(K,UM,UT,Q) + L_ITRANS(Q) * UT_HMULT_UD(K,UM,UT)
                     L_MULTDN(Q) = ( L_FIRST - LC_UT_EMULT_UP(LUM,UT,IB,Q) - L_GAMMA * MULTDN ) * IGAM
                  ENDDO
               endif

!  Upwelling multipliers
!   Rob Change 1/15/21. replace MULTUP to avoid division by BTERM_SAVE
!               MULTUP = UT_PMULT_UU(K,UM,UT) / BTERM_SAVE(K,N) ; IGAM = ONE / GAMMA_P(K,N)
!mick mod 1/5/2021 - removed neg. sign from L_FIRST for consistency with PARTLAYER_STERM_UP

               MULTUP = UT_PMULT_UU(K,UM,UT) ; IGAM = ONE / GAMMA_P(K,N)
               DO Q = 1, NV_PARAMETERS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(K,N,Q)
                  L_FIRST  = ITRANSWDEL * L_UT_HMULT_UU(K,UM,UT,Q) + L_ITRANSWDEL(Q) * UT_HMULT_UU(K,UM,UT)
                  L_MULTUP(Q) = ( L_FIRST + LC_UT_EMULT_UP(LUM,UT,IB,Q) - L_GAMMA * MULTUP ) * IGAM
               ENDDO

!  Complete the multipliers
!   Rob Change 1/15/21. replace this do-loop using regular derivatives (non-logarithmic)
!mick fix 1/5/2021 - changed defs of L_PMULT_DN and L_PMULT_UP here to correspond with
!                    changed defs of UT_PMULT_UD and UT_PMULT_UU in PARTLAYER_STERM_UP

!               DO Q = 1, NV_PARAMETERS
!                  L_PMULT_DN(Q) = ATERM_SAVE(K,N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(K,N,Q) )
!                  L_PMULT_UP(Q) = BTERM_SAVE(K,N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(K,N,Q) )
!               ENDDO
!               DO Q = 1, NV_PARAMETERS
!                  L_PMULT_DN(Q) = ATERM_SAVE(K,N) * L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(K,N,Q)
!                  L_PMULT_UP(Q) = BTERM_SAVE(K,N) * L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(K,N,Q)
!               ENDDO
               DO Q = 1, NV_PARAMETERS
                  L_PMULT_DN(Q) = L_MULTDN(Q)
                  L_PMULT_UP(Q) = L_MULTUP(Q)
               ENDDO

!  Add the particular solution contributions
!    Rob Change 1/15/21. Now have to multiply second terms by ATERM and BTERM. Replacement follows (efficiently)
!mick fix 1/5/2021 - changed formulations for H1 and H2 due to changed defs of L_PMULT_DN and L_PMULT_UP

!               DO Q = 1, NV_PARAMETERS ; DO O1 = 1, NSTOKES
!                  H1 = UHOM_UPDN(UM,O1,K,N) * L_PMULT_DN(Q) + L_UHOM_UPDN(UM,O1,K,N,Q) * UT_PMULT_UD(K,UM,UT)
!                  H2 = UHOM_UPUP(UM,O1,K,N) * L_PMULT_UP(Q) + L_UHOM_UPUP(UM,O1,K,N,Q) * UT_PMULT_UU(K,UM,UT)
!                  SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
!               ENDDO ; ENDDO

               TUDA = UT_PMULT_UD(K,UM,UT) * ATERM_SAVE(K,N)
               TUUB = UT_PMULT_UU(K,UM,UT) * BTERM_SAVE(K,N)
!               DO Q = 1, NV_PARAMETERS ; DO O1 = 1, NSTOKES
!                  H1 = UHOM_UPDN(UM,O1,K,N) * L_PMULT_DN(Q) + L_UHOM_UPDN(UM,O1,K,N,Q) * TUDA
!                  H2 = UHOM_UPUP(UM,O1,K,N) * L_PMULT_UP(Q) + L_UHOM_UPUP(UM,O1,K,N,Q) * TUUB
!                  SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
!               ENDDO ; ENDDO
               DO Q = 1, NV_PARAMETERS
                  L_TUDA = UT_PMULT_UD(K,UM,UT) * L_ATERM_SAVE(K,N,Q) + L_PMULT_DN(Q) * ATERM_SAVE(K,N)
                  L_TUUB = UT_PMULT_UU(K,UM,UT) * L_BTERM_SAVE(K,N,Q) + L_PMULT_UP(Q) * BTERM_SAVE(K,N)
                  DO O1 = 1, NSTOKES
                    H1 = UHOM_UPDN(UM,O1,K,N) * L_TUDA + L_UHOM_UPDN(UM,O1,K,N,Q) * TUDA
                    H2 = UHOM_UPUP(UM,O1,K,N) * L_TUUB + L_UHOM_UPUP(UM,O1,K,N,Q) * TUUB
                    SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
                  ENDDO
               ENDDO

!  End eigenvalue loop

            ENDDO

!  Add Green's function result to the total

            DO Q = 1, NV_PARAMETERS
               DO O1 = 1, NSTOKES
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPARG(O1,Q)
               ENDDO
            ENDDO

!  End user stream loop

         ENDDO

!  End Green's function solution

      ENDIF

!  Add single scatter term if flagged
!    -- 1/31/21. Version 2.8.3. Use of post-processing mask

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          ESUM = UT_EMULT_UP(LUM,UT,IB) 
          DO Q = 1, NV_PARAMETERS
            LC_ESUM = LC_UT_EMULT_UP(LUM,UT,IB,Q)
            DO O1 = 1, NSTOKES
              SPAR = L_UPAR_UP_1(UM,O1,N,Q) * ESUM + UPAR_UP_1(UM,O1,N) * LC_ESUM
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_STERM_UP

!

      SUBROUTINE LC_PARTLAYER_STERM_DN ( &
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                & ! Input flags
           DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SOURCETERM_FLAG, FOURIER,           & ! Input flags/numbers
           IB, N, UT, NV_PARAMETERS, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,  & ! Input numbers
           PARTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTDN_USERM,                       & ! Input Optical/Usertrans
           BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                       & ! input Beam stuff
           K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, UPAR_DN_1,               & ! Input user-solutions
           UPAR_DN_2, UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN, SIGMA_M,                    & ! Input User-solutions/Greens
           ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UT_PMULT_DD, UT_PMULT_DU,           & ! Input Greens
           LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_UT_EMULT_DN,         & ! input Linearized Beam stuff
           L_KEIGEN, NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_UT_HMULT_DD, L_UT_HMULT_DU, & ! input Linearized Homogeneous
           L_UPAR_DN_1, LC_UPAR_DN_2, L_ATERM_SAVE, L_BTERM_SAVE, L_LAYER_TSUP_UTDN,     & ! Input linearized PIs
           L_LAYERSOURCE )                                                                  ! OUTPUT

!  1/31/21. Version 2.8.3. Several Changes.
!    -- Introduce control flag DO_CLASSICAL_SOLUTION, and Taylor-series ordering
!    -- Introduce Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK
!    -- Introduce Green's function inputs as indicated below
!    -- Subroutine has completely new section for the Green's function treatment

!  1/31/21. Version 2.8.3. Add TAYLOR_SMALL to the parameter list

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAX_PARTLAYERS, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4, TAYLOR_SMALL, LDU

      IMPLICIT NONE

!  Inputs
!  ======

!  Control and Bookkeeping
!  -----------------------

!  flags

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES

      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG

!  1/31/21. Version 2.8.3.  Introduce control for RTE solution

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Numbers

      INTEGER, INTENT (IN) ::          N, UT, IB, FOURIER
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES

!  1/31/21. Version 2.8.3.  Introduce User-stream post-processing masking (as in LIDORT0

      INTEGER, INTENT (IN) ::          N_PPSTREAMS
      INTEGER, INTENT (IN) ::          PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  1/31/21. Version 2.8.3. Several new inputs
!    -- Order of Taylor series (including terms up to EPS^n)
!    -- Local vertical optical depth and its linearization
!    -- Transmittance factors for user-defined stream angles, user secants

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER
      DOUBLE PRECISION, INTENT(IN)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN)  :: T_UTDN_USERM  ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, intent(in)  :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Main solution inputs
!  --------------------

!  1/31/21. Version 2.8.3.  Beam average-secant parameterization
!    BEAM_CUTOFF = Last layer to include Particular integral solution

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )

!  Homog Solutions bookkeeping

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  integration constamts

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Homog User solutions and multipliers

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Particular integral user solutions and beam multipliers
!    -- 1/31/21. Version 2.8.3. SIGMA_P is new addition
!  3/1/21. Bug Fix on SIGMA_M declaration

      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1   ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2   ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_M     ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. New Input: Green's function solution variables

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  1/31/21. Version 2.8.3. New Input: Green's function Multipliers

      DOUBLE PRECISION, INTENT (IN) :: UT_PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: UT_PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized inputs
!  -----------------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     -- LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!     -- 1/31/21. Version 2.8.3. New additions for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) ::  LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized eigenvalues
!     -- 1/31/21. Version 2.8.3. New addition for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous solutionss and multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DU  ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DD  ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Particular integral inputs (Beam/Thermal)

      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTDN   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3.  New Input: Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION, INTENT (IN) :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized output

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER          :: M, UM, O1, Q, K, KO1, K0, K1, K2, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SP, H1, H2, TM, TDDA, TDUB
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!  1/31/21. Version 2.8.3. New Help variables for the Green's function solution

      DOUBLE PRECISION :: ITRANS, WDEL, ITRANSWDEL, EPS, PARTA, YFAC, IGAM, LAM
      DOUBLE PRECISION :: SM, MULTDN, MULTUP, FAC1, FAC2, SPAR, ESUM, LC_ESUM
      DOUBLE PRECISION :: L_FIRST, L_PARTA, L_KEG, L_LAM, L_GAMMA, L_WDEL, L_TDDA, L_TDUB
      DOUBLE PRECISION :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)
      DOUBLE PRECISION :: SPARG(MAXSTOKES,MAX_ATMOSWFS)

!  Proxies

      KO1 = K_REAL(N) + 1
      M = FOURIER    ! debug only

!  Very important to zero both output terms (bug solved 12/29/05)
!    -- 1/31/21. Version 2.8.3. Use postprocessing mask.

      DO Q = 1, NV_PARAMETERS
        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only
!  Version 2.8. remove GOTO Statement
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles and Stokes parameters
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES

!  parameter loop

            DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)

                LUXR  = LCON(K,N)   *   UHOM_DNDN(UM,O1,K,N)
                MUXR  = MCON(K,N)   *   UHOM_DNUP(UM,O1,K,N)
                LLUXR = LCON(K,N)   * L_UHOM_DNDN(UM,O1,K,N,Q)
                MLUXR = MCON(K,N)   * L_UHOM_DNUP(UM,O1,K,N,Q)
                NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
                PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)

                H1 = ( NUXR + LLUXR ) *   UT_HMULT_DD(K,UM,UT) &
                           +  LUXR    * L_UT_HMULT_DD(K,UM,UT,Q)
                H2 = ( PUXR + MLUXR ) *   UT_HMULT_DU(K,UM,UT) &
                           +  MUXR    * L_UT_HMULT_DU(K,UM,UT,Q)
                SHOM_R = SHOM_R + H1 + H2

              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

                LUXR1 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K1,N) - &
                         LCON(K2,N)   *   UHOM_DNDN(UM,O1,K2,N)
                LUXR2 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K2,N) + &
                         LCON(K2,N)   *   UHOM_DNDN(UM,O1,K1,N)
                MUXR1 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K1,N) - &
                         MCON(K2,N)   *   UHOM_DNUP(UM,O1,K2,N)
                MUXR2 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K2,N) + &
                         MCON(K2,N)   *   UHOM_DNUP(UM,O1,K1,N)

                LLUXR1 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q) - &
                         LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q)
                LLUXR2 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q) + &
                         LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q)
                MLUXR1 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q) - &
                         MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q)
                MLUXR2 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q) + &
                         MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q)

                NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) - &
                         NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
                NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) + &
                         NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
                PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) - &
                         PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
                PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) + &
                         PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

                H1 =    ( NUXR1 + LLUXR1 ) * UT_HMULT_DD(K1,UM,UT) &
                      - ( NUXR2 + LLUXR2 ) * UT_HMULT_DD(K2,UM,UT) &
                            +  LUXR1       * L_UT_HMULT_DD(K1,UM,UT,Q) &
                            -  LUXR2       * L_UT_HMULT_DD(K2,UM,UT,Q)
                H2 =    ( PUXR1 + MLUXR1 ) * UT_HMULT_DU(K1,UM,UT) &
                      - ( PUXR2 + MLUXR2 ) * UT_HMULT_DU(K2,UM,UT) &
                            +  MUXR1       * L_UT_HMULT_DU(K1,UM,UT,Q) &
                            -  MUXR2       * L_UT_HMULT_DU(K2,UM,UT,Q)

                SHOM_CR = SHOM_CR + H1 + H2

              ENDDO

!  homogeneous contribution

              L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
          ENDDO
        ENDDO

!  End scattering section

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1.0 if solar sources are included (taken care of earlier)
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1 ; TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO Q = 1, NV_PARAMETERS
            L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + L_LAYER_TSUP_UTDN(UM,UT,Q)*TM
          ENDDO
        ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.
!    -- 1/31/21. Version 2.8.3. This condition added

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Particular and single scatter contributions - CLASSICAL SOLUTION
!  ================================================================

!  add the MS particular solution

!  1/31/21. Version 2.8.3. Note use of post-processing mask and classical solution control
!     -- No more separate observational geometry clause

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO LUM = 1, N_PPSTREAMS
          UM   = PPSTREAM_MASK(LUM,IB)
          ESUM = UT_EMULT_DN(LUM,UT,IB)
          DO Q = 1, NV_PARAMETERS
            LC_ESUM = LC_UT_EMULT_DN(LUM,UT,IB,Q)
            DO O1 = 1, NSTOKES
              SPAR = LC_UPAR_DN_2(UM,O1,N,Q) * ESUM + UPAR_DN_2(UM,O1,N) * LC_ESUM
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Particular and single scatter contributions - GREENS FUNCTION SOLUTION
!  ======================================================================

!  add the MS particular solution

!  1/31/21. Version 2.8.3. Completely New section
!     ==> Lattice/Observational geometry, controlled by masking system
!     ==> Multiplier calculations as for LIDORT, real eigenvalues only

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Some beam particulars
!mick mod 1/5/2021 - applied neg. sign to formulations for ITRANSWDEL and L_ITRANSWDEL
!                    for consistency with PARTLAYER_STERM_DN

        WDEL       = T_DELT_MUBAR(N,IB)
        ITRANS     = INITIAL_TRANS(N,IB)
        ITRANSWDEL = - ITRANS * WDEL
        DO Q = 1, NV_PARAMETERS
          L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
          L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
          L_ITRANSWDEL(Q) = - (ITRANS * L_WDEL + WDEL * L_ITRANS(Q))
        ENDDO

!  start local user angle loop and initialize

        DO LUM = 1, N_PPSTREAMS
          UM    = PPSTREAM_MASK(LUM,IB)
          SM    = USER_SECANTS(UM)
          SPARG = ZERO

!  Start eigenvalue loop

          DO K = 1, K_REAL(N)

!  Downwelling multipliers. Note use of -YFAC in the L_2b Call. 01/09/14. Rob checked.

            if ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(K,N)      ; FAC1 = T_UTDN_USERM(UT,UM) 
               YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = T_UTDN_MUBAR(UT,IB)
               PARTA = PARTAU_VERT(UT)   ; LAM  = LOG(ONE/FAC2)/PARTA !  Need average secant here (LAM)
               CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTA, FAC1, FAC2, SM, MULTDN )
               DO Q = 1, NV_PARAMETERS
                  L_KEG   = L_KEIGEN(K,N,Q)
                  L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                  L_PARTA = L_DELTAU_VERT(Q,N) * PARTA ! Input is single normalized 
                  CALL Taylor_Series_L_2b &
                   ( TAYLOR_ORDER, EPS, -YFAC, PARTA, L_PARTA, L_LAM, L_KEG, FAC2, FAC1, SM, LAM, L_MULTDN(Q) )
                  L_MULTDN(Q) = L_MULTDN(Q) * ITRANS + MULTDN * L_ITRANS(Q)
               ENDDO
               MULTDN = MULTDN * ITRANS
            else

!   Rob Change 1/15/21. replace MULTDN to avoid division by ATERM_SAVE
!               MULTDN = UT_PMULT_DD(K,UM,UT) / ATERM_SAVE(K,N) ; IGAM = ONE / GAMMA_M(K,N)

               MULTDN = UT_PMULT_DD(K,UM,UT) ; IGAM = ONE / GAMMA_M(K,N)
               DO Q = 1, NV_PARAMETERS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(K,N,Q)
                  L_FIRST  = ITRANS * L_UT_HMULT_DD(K,UM,UT,Q) + L_ITRANS(Q) * UT_HMULT_DD(K,UM,UT)
                  L_MULTDN(Q) = ( L_FIRST - LC_UT_EMULT_DN(LUM,UT,IB,Q) - L_GAMMA * MULTDN ) * IGAM
               ENDDO
            endif

!  Upwelling multipliers
!   Rob Change 1/15/21. replace MULTUP to avoid division by BTERM_SAVE
!            MULTUP = UT_PMULT_DU(K,UM,UT) / BTERM_SAVE(K,N) ; IGAM = ONE / GAMMA_P(K,N)
!mick mod 1/5/2021 - removed neg. sign from L_FIRST for consistency with PARTLAYER_STERM_DN

            MULTUP = UT_PMULT_DU(K,UM,UT) ; IGAM = ONE / GAMMA_P(K,N)
            DO Q = 1, NV_PARAMETERS
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(K,N,Q)
               L_FIRST  = ITRANSWDEL * L_UT_HMULT_DU(K,UM,UT,Q) + L_ITRANSWDEL(Q) * UT_HMULT_DU(K,UM,UT)
               L_MULTUP(Q) = ( L_FIRST + LC_UT_EMULT_DN(LUM,UT,IB,Q) - L_GAMMA * MULTUP ) * IGAM
            ENDDO

!  Complete the multipliers
!   Rob Change 1/15/21. replace this do-loop using regular derivatives (non-logarithmic)
!mick fix 1/5/2021 - changed defs of L_PMULT_DN and L_PMULT_UP here to correspond with
!                    changed defs of UT_PMULT_DD and UT_PMULT_DU in PARTLAYER_STERM_DN

!            DO Q = 1, NV_PARAMETERS
!               L_PMULT_DN(Q) = ATERM_SAVE(K,N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(K,N,Q) )
!               L_PMULT_UP(Q) = BTERM_SAVE(K,N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(K,N,Q) )
!            ENDDO
!            DO Q = 1, NV_PARAMETERS
!               L_PMULT_DN(Q) = ATERM_SAVE(K,N) * L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(K,N,Q)
!               L_PMULT_UP(Q) = BTERM_SAVE(K,N) * L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(K,N,Q)
!            ENDDO
            DO Q = 1, NV_PARAMETERS
               L_PMULT_DN(Q) = L_MULTDN(Q)
               L_PMULT_UP(Q) = L_MULTUP(Q)
            ENDDO

!  Add the particular solution contributions
!    Rob Change 1/15/21. Now have to multiply second terms by ATERM and BTERM. Replacement follows (efficiently)
!mick fix 1/5/2021 - changed formulations for H1 and H2 due to changed defs of L_PMULT_DN and L_PMULT_UP

!            DO Q = 1, NV_PARAMETERS ; DO O1 = 1, NSTOKES
!                H1 = UHOM_DNDN(UM,O1,K,N) * L_PMULT_DN(Q) + L_UHOM_DNDN(UM,O1,K,N,Q) * UT_PMULT_DD(K,UM,UT) 
!                H2 = UHOM_DNUP(UM,O1,K,N) * L_PMULT_UP(Q) + L_UHOM_DNUP(UM,O1,K,N,Q) * UT_PMULT_DU(K,UM,UT) 
!                SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
!            ENDDO ; ENDDO

            TDDA = UT_PMULT_DD(K,UM,UT) * ATERM_SAVE(K,N)
            TDUB = UT_PMULT_DU(K,UM,UT) * BTERM_SAVE(K,N)
!            DO Q = 1, NV_PARAMETERS
!              DO O1 = 1, NSTOKES
!                H1 = UHOM_DNDN(UM,O1,K,N) * L_PMULT_DN(Q) + L_UHOM_DNDN(UM,O1,K,N,Q) * TDDA
!                H2 = UHOM_DNUP(UM,O1,K,N) * L_PMULT_UP(Q) + L_UHOM_DNUP(UM,O1,K,N,Q) * TDUB 
!                SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
!              ENDDO
!            ENDDO
            DO Q = 1, NV_PARAMETERS
              L_TDDA = UT_PMULT_DD(K,UM,UT) * L_ATERM_SAVE(K,N,Q) + L_PMULT_DN(Q) * ATERM_SAVE(K,N) 
              L_TDUB = UT_PMULT_DU(K,UM,UT) * L_BTERM_SAVE(K,N,Q) + L_PMULT_UP(Q) * BTERM_SAVE(K,N)
              DO O1 = 1, NSTOKES
                H1 = UHOM_DNDN(UM,O1,K,N) * L_TDDA + L_UHOM_DNDN(UM,O1,K,N,Q) * TDDA
                H2 = UHOM_DNUP(UM,O1,K,N) * L_TDUB + L_UHOM_DNUP(UM,O1,K,N,Q) * TDUB
                SPARG(O1,Q) = SPARG(O1,Q) + H1 + H2
              ENDDO
            ENDDO

!  End eigenvalue loop

          ENDDO

!  Add Green's function result to the total

          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPARG(O1,Q)
            ENDDO
          ENDDO

!  End user stream loop

        ENDDO

!  End Green's solution

      ENDIF

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM   = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NSTOKES
               DO Q = 1, NV_PARAMETERS
                  SP = L_UPAR_DN_1(UM,O1,N,Q) *    UT_EMULT_DN(LUM,UT,IB) + &
                         UPAR_DN_1(UM,O1,N)   * LC_UT_EMULT_DN(LUM,UT,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_STERM_DN

!

      SUBROUTINE QUADCOLUMNWF_LEVEL_UP ( &
           DO_THERMAL_TRANSONLY, NL, UTA, NV_PARAMETERS, NSTOKES,             & ! Input flag/numbers
           NSTREAMS, NLAYERS, FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,  & ! Input numbers/flux/quad/D.O.trans
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,             & ! Input Homog. solutions
           LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,                          & ! Input thermal and PI solutions
           L_T_DELT_DISORDS, L_SOLA_XPOS, L_SOLB_XNEG,  L_T_DELT_EIGEN,       & ! Linearized Homog. solutions
           L_WLOWER, L_WUPPER, NCON, PCON, L_T_WUPPER, L_BOA_THTONLY_SOURCE,  & ! Linearized thermal and PI solutions
           QATMOSWF_F )                                                         ! OUTPUT

!   1/31/21. Version 2.8.3. No Changes

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXLAYERS, &
                                 MAX_USER_LEVELS, MAX_SZANGLES, MAX_ATMOSWFS,  &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO 

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Input basic numbers

      INTEGER, INTENT (IN) ::          NL, UTA, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  Discrete Ordinate solutions

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  Linearized Input
!  ================

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS  ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Discrete Ordinate solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE  ( MAXSTREAMS, MAX_ATMOSWFS )

!  Output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION :: LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION :: TPROP, THELP, L_TPROP, L_THELP, FM

!  homogeneous and particular solution contributions SHOM and SPAR

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N  = NL + 1
      FM = FLUX_MULTIPLIER

!  For the lowest level
!  ====================

!  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY .and. NL.EQ.NLAYERS  ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
            QATMOSWF_F(Q,UTA,I,O1) = FM * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For the lowest level, scattering solution
!  -----------------------------------------

      IF ( NL .EQ. NLAYERS ) THEN

!  Stokes and streams loops

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  Parameter loop

            DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(NL)

                LXR  = LCON(K,NL)   *   SOLA_XPOS(I1,O1,K,NL)
                LLXR = LCON(K,NL)   * L_SOLA_XPOS(I1,O1,K,NL,Q)
                MLXR = MCON(K,NL)   * L_SOLB_XNEG(I1,O1,K,NL,Q)
                NXR  = NCON(K,NL,Q) *   SOLA_XPOS(I1,O1,K,NL)
                PXR  = PCON(K,NL,Q) *   SOLB_XNEG(I1,O1,K,NL)

                HOM1 = ( LLXR + NXR )  *   T_DELT_EIGEN(K,NL) &
                             +  LXR    * L_T_DELT_EIGEN(K,NL,Q)
                HOM2 = PXR + MLXR
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(NL) + 1
              DO K = 1, K_COMPLEX(NL)
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

                NXR1  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL) &
                        - NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL)
                NXR2  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL) &
                        + NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL)
                PXR1  =   PCON(K1,NL,Q) *   SOLB_XNEG(I1,O1,K1,NL) &
                        - PCON(K2,NL,Q) *   SOLB_XNEG(I1,O1,K2,NL)

                LXR1  =   LCON(K1,NL) *   SOLA_XPOS(I1,O1,K1,NL) &
                        - LCON(K2,NL) *   SOLA_XPOS(I1,O1,K2,NL)
                LXR2  =   LCON(K1,NL) *   SOLA_XPOS(I1,O1,K2,NL) &
                        + LCON(K2,NL) *   SOLA_XPOS(I1,O1,K1,NL)

                LLXR1  =   LCON(K1,NL) * L_SOLA_XPOS(I1,O1,K1,NL,Q) &
                         - LCON(K2,NL) * L_SOLA_XPOS(I1,O1,K2,NL,Q)
                LLXR2  =   LCON(K1,NL) * L_SOLA_XPOS(I1,O1,K2,NL,Q) &
                         + LCON(K2,NL) * L_SOLA_XPOS(I1,O1,K1,NL,Q)
                MLXR1  =   MCON(K1,NL) * L_SOLB_XNEG(I1,O1,K1,NL,Q) &
                         - MCON(K2,NL) * L_SOLB_XNEG(I1,O1,K2,NL,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,NL) &
                         - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,NL)
                HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,NL,Q) &
                                   - LXR2   * L_T_DELT_EIGEN(K2,NL,Q)
                HOM3CR = PXR1 + MLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR
              ENDDO

!  real part, add particular solution, complete result

              SHOM = SHOM_R + SHOM_CR
              SPAR = L_WLOWER(I1,O1,NL,Q)
              QATMOSWF_F(Q,UTA,I,O1) = FM * ( SPAR + SHOM )

!  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

!  End lowest level clause

      ENDIF

!  For other levels in the atmosphere
!  ==================================

!  For other levels, thermal transmittance only

      IF ( NL .NE. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
            THELP   =   BOA_THTONLY_SOURCE(I,O1)
            DO LAY = NLAYERS, N, -1
              L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
              L_TPROP = L_T_WUPPER(I1,LAY,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              TPROP = T_WUPPER(I1,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            QATMOSWF_F(Q,UTA,I,O1) = FM * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For other levels, scattering solution
!  -------------------------------------

      IF ( NL .NE. NLAYERS ) THEN

!  Stokes and streams loops

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  Parameter loop

            DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I1,O1,K,N)
                MXR  = MCON(K,N)   *   SOLB_XNEG(I1,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I1,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 =   NXR + LLXR
                HOM2 = ( PXR + MLXR ) *   T_DELT_EIGEN(K,N) &
                             + MXR    * L_T_DELT_EIGEN(K,N,Q)
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N) &
                        + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)

                MXR1  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K1,N) &
                        - MCON(K2,N) *   SOLB_XNEG(I1,O1,K2,N)
                MXR2  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K2,N) &
                        + MCON(K2,N) *   SOLB_XNEG(I1,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                         - LCON(K2,N) * L_SOLA_XPOS(I1,O1,K2,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K1,N,Q) &
                         - MCON(K2,N) * L_SOLB_XNEG(I1,O1,K2,N,Q)
                MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K2,N,Q) &
                         + MCON(K2,N) * L_SOLB_XNEG(I1,O1,K1,N,Q)

                HOM1CR =   ( PXR1 + MLXR1 ) *   T_DELT_EIGEN(K1,N) &
                         - ( PXR2 + MLXR2 ) *   T_DELT_EIGEN(K2,N)
                HOM2CR =             MXR1   * L_T_DELT_EIGEN(K1,N,Q) &
                                   - MXR2   * L_T_DELT_EIGEN(K2,N,Q)
                HOM3CR = NXR1 + LLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR

              ENDDO

!  real part, add particular solution, complete result

              SHOM = SHOM_R + SHOM_CR
              SPAR = L_WUPPER(I1,O1,N,Q)
              QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

!  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

!  End variability clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADCOLUMNWF_LEVEL_UP

!

      SUBROUTINE QUADCOLUMNWF_LEVEL_DN ( &
           DO_THERMAL_TRANSONLY, NL, UTA, NV_PARAMETERS, NSTOKES,             & ! Input flag/numbers
           NSTREAMS, FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,           & ! Input numbers/flux/quad/D.O.trans
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, & ! Input RT solutions
           T_WLOWER, L_T_DELT_DISORDS, L_SOLA_XPOS, L_SOLB_XNEG,              & ! Linearized Homog. solutions
           L_T_DELT_EIGEN, L_WLOWER,  NCON, PCON, L_T_WLOWER,                 & ! Linearized thermal and PI solutions
           QATMOSWF_F )                                                         ! OUTPUT

!   1/31/21. Version 2.8.3. No Changes

!  Downwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXLAYERS, &
                                 MAX_USER_LEVELS, MAX_SZANGLES, MAX_ATMOSWFS,  &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO 

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Input basic numbers

      INTEGER, INTENT (IN) ::          NL, UTA, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  Discrete Ordinate solutions

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN  ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  Linearized Input
!  ================

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS  ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Discrete Ordinate solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

!  Output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1
      DOUBLE PRECISION :: LXR, LXR1, LXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2
      DOUBLE PRECISION :: TPROP, THELP, L_TPROP, L_THELP

!  Downwelling weighting function at TOA ( or N = 0 ) is zero
!    Zero and return

      IF ( NL .EQ. 0 ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              QATMOSWF_F(Q,UTA,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For other levels in the atmosphere
!  ----------------------------------

!  Other levels, Thermal transmittance-only solution
!  Thermal transmittance solution, build from TOA downwards
!  Scattering solution, use the Discrete Ordinate solution

      IF ( NL.NE.0 .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = ZERO
            THELP   = ZERO
            DO LAY = 1, NL
              L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
              L_TPROP = L_T_WLOWER(I,LAY,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP +  L_TPROP + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              TPROP = T_WLOWER(I,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  --------------------

      IF ( NL.NE.0 ) THEN

!  Shorthand

        N = NL

!  Stokes and streams loops

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  Parameter loop

            DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_DELT_EIGEN(K,N) &
                             +  LXR   * L_T_DELT_EIGEN(K,N,Q)
                HOM2 = PXR + MLXR
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N) &
                        - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N) &
                        + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q) &
                         - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q) &
                         + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q) &
                         - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,N) &
                         - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,N)
                HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,N,Q) &
                                   - LXR2   * L_T_DELT_EIGEN(K2,N,Q)
                HOM3CR = PXR1 + MLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR

              ENDDO

!  real part, add particular solution, complete result

              SHOM = SHOM_R + SHOM_CR
              SPAR = L_WLOWER(I,O1,N,Q)
              QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

!  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

!  Finish variability clauses

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADCOLUMNWF_LEVEL_DN

!

      SUBROUTINE LC_QUAD_GFUNCMULT &
          ( IB, UT, N, NV_PARAMETERS, TAYLOR_ORDER, PARTAU_VERT, L_DELTAU_VERT,      & ! input
            BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,  & ! input
            T_UTUP_EIGEN, T_UTDN_EIGEN, ATERM_SAVE, BTERM_SAVE,                      & ! input
            K_REAL, GAMMA_M, GAMMA_P, UT_GMULT_UP, UT_GMULT_DN,                      & ! input
            LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,  & ! input
            L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,    & ! input
            LC_UT_GMULT_UP, LC_UT_GMULT_DN )                                           ! output

!  Linearization of Quadrature Green's function multiplier (off-grid only)

!  1/31/21. Version 2.8.3. Completely new subroutine
!     -- Modeled after the LIDORT Version 3.7 routine of the same nams

!  module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXSTREAMS, MAXSTOKES, MAX_PARTLAYERS, MAXEVALUES, &
                                 MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Beam index,  offgrid indices, number of parameters

      INTEGER  , intent(in)  :: IB, N, UT, NV_PARAMETERS

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Input optical properties after delta-M scaling
!  Linearized input optical properties after delta-M scaling

      DOUBLE PRECISION, intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Beam parameterization

      INTEGER         , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)
      DOUBLE PRECISION, intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, intent(in)  :: AVERAGE_SECANT ( MAXLAYERS,      MAXBEAMS )
      DOUBLE PRECISION, intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS,      MAXBEAMS )
      DOUBLE PRECISION, intent(in)  :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  transmittance factors for +/- eigenvalues at User optical depths (UTUP and UTDN)

      DOUBLE PRECISION, intent(in)  :: T_UTUP_EIGEN(MAXEVALUES,MAX_PARTLAYERS)
      DOUBLE PRECISION, intent(in)  :: T_UTDN_EIGEN(MAXEVALUES,MAX_PARTLAYERS)

!  Green's function particular integral arrays

      DOUBLE PRECISION, intent(in)  :: ATERM_SAVE(MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: BTERM_SAVE(MAXEVALUES,MAXLAYERS)

!  Number of reall eigenvalues

      INTEGER         , intent(in)  :: K_REAL(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      DOUBLE PRECISION, intent(in)  :: GAMMA_M(MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: GAMMA_P(MAXEVALUES,MAXLAYERS)

!  Green functions multipliers for off-grid optical depths

      DOUBLE PRECISION, intent(in)  :: UT_GMULT_UP(MAXEVALUES,MAX_PARTLAYERS)
      DOUBLE PRECISION, intent(in)  :: UT_GMULT_DN(MAXEVALUES,MAX_PARTLAYERS)

!  Linearized transittance factors for solar beams.

      DOUBLE PRECISION, intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      DOUBLE PRECISION, intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized transmittances, eigensolutions

      DOUBLE PRECISION, intent(in)  :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized (Positive) Eigenvalues

      DOUBLE PRECISION, intent(in)  :: L_KEIGEN(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution
!mick fix 1/5/2021 - changed 1st dim of L_BTERM_SAVE from MAXSTREAMS to MAXEVALUES

      DOUBLE PRECISION, intent(in)  :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent(in)  :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized Green functions multipliers for off-grid optical depths

      DOUBLE PRECISION, intent(inout) :: LC_UT_GMULT_UP(MAXEVALUES,MAX_PARTLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent(inout) :: LC_UT_GMULT_DN(MAXEVALUES,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER          :: Q, K
      DOUBLE PRECISION :: ZX_DN, ZX_UP, ZW, WX, WDEL, ITRANS
      DOUBLE PRECISION :: L_ZX_DN, L_ZW, L_WX(MAX_ATMOSWFS), L_WDEL(MAX_ATMOSWFS)

      DOUBLE PRECISION :: MULTDN, MULTUP, EPS, PARTA, LAM, IGAM, ITA, ITB, GUP, GDN
      DOUBLE PRECISION :: L_PARTA, L_KEG, L_LAM, L_GAMMA
      DOUBLE PRECISION :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_ITRANS(MAX_ATMOSWFS)

!  No particular solution beyond the cutoff layer.
!    [ Zero the multiplier values and exit )
!mick fix 1/5/2021 - initialize all K entries of LC_UT_GMULT_UP & LC_UT_GMULT_DN
!                    whether or not beam is cut off

      !IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
      !  DO K = 1, K_REAL(N)
      !    DO Q = 1, NV_PARAMETERS
      !      LC_UT_GMULT_UP(K,UT,Q) = ZERO
      !      LC_UT_GMULT_DN(K,UT,Q) = ZERO
      !     ENDDO
      !  ENDDO
      !  RETURN
      !ENDIF

      DO Q = 1, NV_PARAMETERS
        LC_UT_GMULT_UP(:,UT,Q) = ZERO
        LC_UT_GMULT_DN(:,UT,Q) = ZERO
      ENDDO
      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Layer constant terms

      WX     = T_UTDN_MUBAR(UT,IB)
      WDEL   = T_DELT_MUBAR(N,IB)
      ITRANS = INITIAL_TRANS(N,IB)
      DO Q = 1, NV_PARAMETERS
         L_ITRANS(Q) = LC_INITIAL_TRANS(N,IB,Q)
         L_WX(Q)     = LC_T_UTDN_MUBAR(UT,IB,Q)
         L_WDEL(Q)   = LC_T_DELT_MUBAR(N,IB,Q)
      ENDDO

!  No difference for Pseudo-spherical (average secant) or plane parallel
!   Except that LC_AVERAGE_SCANT is zero

      DO K = 1, K_REAL(N)

         ZX_DN = T_UTDN_EIGEN(K,UT)
         ZX_UP = T_UTUP_EIGEN(K,UT)
         ZW    = WDEL * ZX_UP

!  Downwelling multipliers
!   -- Have to introduce Taylor series calculation of MULTDN, as well as L_MULTDN

         IF ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
            EPS   = GAMMA_M(K,N)
            PARTA = PARTAU_VERT(UT)
            LAM   = AVERAGE_SECANT(N,IB)
            DO Q = 1, NV_PARAMETERS
               L_LAM  = LC_AVERAGE_SECANT(N,IB,Q) ; L_KEG  = L_KEIGEN(K,N,Q)
               L_PARTA = L_DELTAU_VERT(Q,N) * PARTA ! Input is single normalized
               CALL TAYLOR_SERIES_1   ( TAYLOR_ORDER, EPS, PARTA, WX, ONE, MULTDN )
               CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, PARTA, L_PARTA, L_KEG, L_LAM, WX, LAM, L_MULTDN(Q) )
            ENDDO
         ELSE
            IGAM = ONE / GAMMA_M(K,N) ; MULTDN =  ( T_UTDN_EIGEN(K,UT) - WX ) * IGAM
            DO Q = 1, NV_PARAMETERS
               L_ZX_DN = L_T_UTDN_EIGEN(K,UT,Q)
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(K,N,Q)
               L_MULTDN(Q) = ( ( L_ZX_DN - L_WX(Q) ) - MULTDN * L_GAMMA ) * IGAM
            ENDDO
         ENDIF

!  upwelling multipliers

         IGAM = ONE / GAMMA_P(K,N) ; MULTUP =  ( WX - ZW ) * IGAM
         DO Q = 1, NV_PARAMETERS
            L_ZW     = WDEL * L_T_UTUP_EIGEN(K,UT,Q) + L_WDEL(Q) * ZX_UP
            L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(K,N,Q)
            L_MULTUP(Q) = ( ( L_WX(Q) - L_ZW) - MULTUP * L_GAMMA ) * IGAM
         ENDDO

!  Complete the multipliers

         ITA = ITRANS * ATERM_SAVE(K,N); GDN = UT_GMULT_DN(K,UT)
         ITB = ITRANS * BTERM_SAVE(K,N); GUP = UT_GMULT_UP(K,UT)

!   Rob Change 1/15/21. replace this do-loop using regular derivatives (non-logarithmic)
!         DO Q = 1, NV_PARAMETERS
!            LC_UT_GMULT_DN(K,UT,Q) = ITA * L_MULTDN(Q) + GDN * ( L_ITRANS(Q) + L_ATERM_SAVE(K,N,Q) )
!            LC_UT_GMULT_UP(K,UT,Q) = ITB * L_MULTUP(Q) + GUP * ( L_ITRANS(Q) + L_BTERM_SAVE(K,N,Q) )
!         ENDDO

!mick fix 1/5/2021 - replaced GDN and GUP with MULTDN * ATERM_SAVE(K,N) and MULTUP * BTERM_SAVE(K,N), respectively (hold)
         DO Q = 1, NV_PARAMETERS
            LC_UT_GMULT_DN(K,UT,Q) = ITA * L_MULTDN(Q) + GDN * L_ITRANS(Q) + MULTDN * ITRANS * L_ATERM_SAVE(K,N,Q)  
            LC_UT_GMULT_UP(K,UT,Q) = ITB * L_MULTUP(Q) + GUP * L_ITRANS(Q) + MULTUP * ITRANS * L_BTERM_SAVE(K,N,Q)
            !LC_UT_GMULT_DN(K,UT,Q) = ITA * L_MULTDN(Q) + MULTDN * ATERM_SAVE(K,N) * L_ITRANS(Q) &
            !                       + MULTDN * ITRANS * L_ATERM_SAVE(K,N,Q)  
            !LC_UT_GMULT_UP(K,UT,Q) = ITB * L_MULTUP(Q) + MULTUP * BTERM_SAVE(K,N) * L_ITRANS(Q) &
            !                       + MULTUP * ITRANS * L_BTERM_SAVE(K,N,Q)
         ENDDO

!   End Discrete ordinate loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_QUAD_GFUNCMULT

!

      SUBROUTINE QUADCOLUMNWF_OFFGRID_UP ( &
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, DO_CLASSICAL, & ! Input flags
           N, UTA, UT, IB, NV_PARAMETERS, NSTOKES, NSTREAMS, NLAYERS, NSTKS_NSTRMS,     & ! Input numbers
           FLUX_MULTIPLIER, QUAD_STREAMS, T_UTDN_MUBAR, LC_T_UTDN_MUBAR,                & ! Input Flux/Beam/Quad
           T_DELT_DISORDS, T_DISORDS_UTUP, L_T_DELT_DISORDS, L_T_DISORDS_UTUP,          & ! Input disords
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN,         & ! Input Homogeneous 
           UT_GMULT_UP, UT_GMULT_DN, LCON, MCON, WUPPER,                                & ! Input SolarPI
           T_WUPPER, BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,                          & ! Input Thermal
           L_SOLA_XPOS, L_SOLB_XNEG, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, NCON, PCON,        & ! Input Lin. Homogeneous
           LC_UT_GMULT_UP, LC_UT_GMULT_DN, L_WUPPER, L_T_WUPPER, L_UT_T_PARTIC,         & ! Input Lin. SolarPI/Thermal
           QATMOSWF_F )                                                                   ! Output

!  1/31/21. Version 2.8.3. Rewritten to incorporate new Linearized Green's function solution
!     -- New code modeled after the LIDORT Version 3.8.1, except include DO_CLASSICAL
!     -- New arguments: multipliers UT_GMULT_UP, UT_GMULT_DN, LC_UT_GMULT_UP, LC_UT_GMULT_DN

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, &
                                 MAX_PARTLAYERS, MAX_USER_LEVELS, MAX_ATMOSWFS,          &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3. New flag for solution options

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL

!  Input basic numbers

      INTEGER, INTENT (IN) ::          N, UTA, UT, IB, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS, NLAYERS, NSTKS_NSTRMS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Beam parameterization (Initial trans, average secantm etc...)

      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS   ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP   ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Homogeneous solutions and Eigenstream transmittances

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )

!  Green's functions multipliers for off-grid optical depths
!    -- 1/31//21. Version 2.8.3. Added for the Green's function calculation

      DOUBLE PRECISION, intent (in) :: UT_GMULT_UP(MAXEVALUES,MAX_PARTLAYERS)
      DOUBLE PRECISION, intent (in) :: UT_GMULT_DN(MAXEVALUES,MAX_PARTLAYERS)

!  BVP constants and particular integrals

      DOUBLE PRECISION, INTENT (IN) :: WUPPER     ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  Linearized Input
!  ================

!  LInearized Beam parameterization

      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Homogeneous solution variables

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  LInearized PI and BVP constants

      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE  ( MAXSTREAMS, MAX_ATMOSWFS )

!  Linearized Green functions multipliers for off-grid optical depths
!    -- 1/31//21. Version 2.8.3. Added for the Green's function calculation

      DOUBLE PRECISION, intent(in) :: LC_UT_GMULT_UP(MAXEVALUES,MAX_PARTLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent(in) :: LC_UT_GMULT_DN(MAXEVALUES,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  output
!  ------

!  Quadrature-defined weighting functions

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

!  @@@@@@@ Rob fix, 2/9/11, added variables SUNP, L_SUNP

      INTEGER ::          I, I1, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM, SUNP, L_SUNP
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR, HOM4CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION :: LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION :: TPROP, THELP, L_TPROP, L_THELP, FMULT

!  short hand

      FMULT = FLUX_MULTIPLIER

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, NV_PARAMETERS
            THELP     =   BOA_THTONLY_SOURCE(I,O1)
            L_THELP   = L_BOA_THTONLY_SOURCE(I,Q)
            DO LAY = NLAYERS, N+1, -1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              L_TPROP = L_T_WUPPER(I1,LAY,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              TPROP = T_WUPPER(I1,Q) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTUP(I,UT)
            L_TPROP = L_UT_T_PARTIC(I1,UT,Q) / QUAD_STREAMS(I)
            L_THELP = L_THELP + L_TPROP + THELP * L_T_DISORDS_UTUP(I,UT,Q)
            QATMOSWF_F(Q,UTA,I,O1) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  ###########################################

!  Homogeneous
!  -----------

!  stream directions

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES

!  Parameter loop

          DO Q = 1, NV_PARAMETERS

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR  = LCON(K,N)   *   SOLA_XPOS(I1,O1,K,N)
              MXR  = MCON(K,N)   *   SOLB_XNEG(I1,O1,K,N)
              LLXR = LCON(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q)
              MLXR = MCON(K,N)   * L_SOLB_XNEG(I1,O1,K,N,Q)
              NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
              HOM1 = ( NXR + LLXR ) *   T_UTDN_EIGEN(K,UT) &
                           +  LXR   * L_T_UTDN_EIGEN(K,UT,Q)
              HOM2 = ( PXR + MLXR ) *   T_UTUP_EIGEN(K,UT) &
                           +  MXR   * L_T_UTUP_EIGEN(K,UT,Q)
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                      - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                      + NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N) &
                      - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
              PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N) &
                      + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)

              LXR1  =   LCON(K1,N) *   SOLA_XPOS(I1,O1,K1,N) &
                      - LCON(K2,N) *   SOLA_XPOS(I1,O1,K2,N)
              LXR2  =   LCON(K1,N) *   SOLA_XPOS(I1,O1,K2,N) &
                      + LCON(K2,N) *   SOLA_XPOS(I1,O1,K1,N)
              MXR1  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K1,N) &
                      - MCON(K2,N) *   SOLB_XNEG(I1,O1,K2,N)
              MXR2  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K2,N) &
                      + MCON(K2,N) *   SOLB_XNEG(I1,O1,K1,N)

              LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                       - LCON(K2,N) * L_SOLA_XPOS(I1,O1,K2,N,Q)
              LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                       + LCON(K2,N) * L_SOLA_XPOS(I1,O1,K1,N,Q)
              MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K1,N,Q) &
                       - MCON(K2,N) * L_SOLB_XNEG(I1,O1,K2,N,Q)
              MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K2,N,Q) &
                       + MCON(K2,N) * L_SOLB_XNEG(I1,O1,K1,N,Q)

              HOM1CR =   ( NXR1 + LLXR1 ) *   T_UTDN_EIGEN(K1,UT) &
                       - ( NXR2 + LLXR2 ) *   T_UTDN_EIGEN(K2,UT)
              HOM2CR =             LXR1   * L_T_UTDN_EIGEN(K1,UT,Q) &
                                 - LXR2   * L_T_UTDN_EIGEN(K2,UT,Q)
              HOM3CR =   ( PXR1 + MLXR1 ) *   T_UTUP_EIGEN(K1,UT) &
                       - ( PXR2 + MLXR2 ) *   T_UTUP_EIGEN(K2,UT)
              HOM4CR =             MXR1   * L_T_UTUP_EIGEN(K1,UT,Q) &
                                 - MXR2   * L_T_UTUP_EIGEN(K2,UT,Q)

              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR + HOM4CR

            ENDDO

!  real part

            SHOM = SHOM_R + SHOM_CR
            QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

!  Finish Q, I and O1 loops

          ENDDO
        ENDDO
      ENDDO

!  Add the linearized thermal solution  (if flagged)
!     This is the thermal solution with scattering (NOT TRANSONLY)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, NV_PARAMETERS
            SPAR = L_UT_T_PARTIC(I1,UT,Q)
            QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
          ENDDO
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN
 
!  Add Linearized Classical beam solution
!  ======================================

!   -- 1/31/21. Section now has the classical flag returned

! @@@@@@@@@@ Rob fix, 2/9/11, Replacement of commented section
!      IF ( DO_CLASSICAL_SOLUTION ) THEN
!        DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!            DO Q = 1, NV_PARAMETERS
!              SPAR = L_WUPPER(I1,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
!                       WUPPER(I1,O1,N)   * LC_T_UTDN_MUBAR(UT,IB,Q)
!              QATMOSWF_F(Q,UTA,I,O1) = &
!                QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
!            ENDDO
!          ENDDO
!        ENDDO
!      ELSE
!                P L A C E H O L D E R
!      ENDIF
! @@@@@@@@@@ Rob fix, 2/9/11, End Replaced section

      IF ( DO_CLASSICAL ) THEN
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SUNP   = WUPPER(I1,O1,N)
                L_SUNP = L_WUPPER(I1,O1,N,Q)
                IF ( O1 .eq. 1 ) THEN
                  SUNP   = SUNP - T_WUPPER(I1,N)
                  L_SUNP = L_SUNP - L_T_WUPPER(I1,N,Q)
                ENDIF
                SPAR = L_SUNP * T_UTDN_MUBAR(UT,IB) + SUNP * LC_T_UTDN_MUBAR(UT,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = L_WUPPER(I1,O1,N,Q) * T_UTDN_MUBAR(UT,IB) + WUPPER(I1,O1,N)* LC_T_UTDN_MUBAR(UT,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Linearized Green's function solution
! =====================================

!  1/31/21. Version 2.8.3. COMPLETELY NEW SECTION
!     -- Follows the LIDORT implementation closely

     IF ( .not. DO_CLASSICAL ) THEN

!  Sum up the Green's function contributions

         DO O1 = 1, NSTOKES
           DO I = 1, NSTREAMS
             DO Q = 1, NV_PARAMETERS
               SPAR = DOT_PRODUCT (    UT_GMULT_DN(1:NSTKS_NSTRMS,UT),   L_SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N,Q) ) &
                    + DOT_PRODUCT ( LC_UT_GMULT_DN(1:NSTKS_NSTRMS,UT,Q),   SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N)   ) &
                    + DOT_PRODUCT (    UT_GMULT_UP(1:NSTKS_NSTRMS,UT),   L_SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N,Q) ) &
                    + DOT_PRODUCT ( LC_UT_GMULT_UP(1:NSTKS_NSTRMS,UT,Q),   SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N)   )
               IF ( O1.gt.2) SPAR = -SPAR
               QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1)  + FMULT * SPAR
             ENDDO
           ENDDO
         ENDDO

!  End Green's function solution

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADCOLUMNWF_OFFGRID_UP

!

      SUBROUTINE QUADCOLUMNWF_OFFGRID_DN ( &
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, DO_CLASSICAL, & ! Input flags
           N, UTA, UT, IB, NV_PARAMETERS, NSTOKES, NSTREAMS, NSTKS_NSTRMS,              & ! Input numbers
           FLUX_MULTIPLIER, QUAD_STREAMS, T_UTDN_MUBAR, LC_T_UTDN_MUBAR,                    & ! Input solar beam
           T_DELT_DISORDS, T_DISORDS_UTDN, L_T_DELT_DISORDS, L_T_DISORDS_UTDN,              & ! Input disords
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN,             & ! Input Homog. solution
           UT_GMULT_UP, UT_GMULT_DN, LCON, MCON, WUPPER, T_WLOWER, T_WUPPER,                & ! Input SOlat/Thermal PIs
           L_SOLA_XPOS, L_SOLB_XNEG, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, NCON, PCON,            & ! Input Lin. Homog.
           LC_UT_GMULT_UP, LC_UT_GMULT_DN, L_WUPPER, L_T_WLOWER, L_T_WUPPER, L_UT_T_PARTIC, & ! Input Lin PIs
           QATMOSWF_F )                                                                        ! Output

!  1/31/21. Version 2.8.3. Rewritten to incorporate new Linearized Green's function solution
!     -- New code modeled after the LIDORT Version 3.8.1, except include DO_CLASSICAL
!     -- New arguments: multipliers UT_GMULT_UP, UT_GMULT_DN, LC_UT_GMULT_UP, LC_UT_GMULT_DN

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, &
                                 MAX_PARTLAYERS, MAX_USER_LEVELS, MAX_ATMOSWFS,          &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO

      implicit none

!  Inputs
!  ======

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3. New flag for solution options

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL

!  Input basic numbers

      INTEGER, INTENT (IN) ::          N, UTA, UT, IB, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS, NSTKS_NSTRMS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Beam parameterization (Initial trans, average secantm etc...)
!    -- 1/31/21. Version 2.8.3. New arguments

      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS   ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN   ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Homogeneous solutions and Eigenstream transmittances

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )

!  Green function multipliers for off-grid optical depths
!    -- 1/31//21. Version 2.8.3. Added for the Green's function calculation

      DOUBLE PRECISION, intent(in)  :: UT_GMULT_UP ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, intent(in)  :: UT_GMULT_DN ( MAXEVALUES, MAX_PARTLAYERS )

!  BVP constants and particular integrals

      DOUBLE PRECISION, INTENT (IN) :: WUPPER     ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  Linearized Input
!  ================

!  LInearized Beam parameterization

      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Homogeneous solution variables

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  LInearized PI and BVP constants

      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

!  output
!  ------

!  Quadrature-defined weighting functions

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  Linearized Green functions multipliers for off-grid optical depths
!    -- 1/31/21. Version 2.8.3. Added for the Green's function calculation

      DOUBLE PRECISION, intent(in) :: LC_UT_GMULT_UP(MAXEVALUES,MAX_PARTLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent(in) :: LC_UT_GMULT_DN(MAXEVALUES,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  @@@@@@@ Rob fix, 2/9/11, added variables SUNP, L_SUNP

      INTEGER ::          I, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM, SUNP, L_SUNP
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR, HOM4CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION :: LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION :: TPROP, THELP, L_TPROP, L_THELP, FMULT

!  Short hand

      FMULT = FLUX_MULTIPLIER

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = ZERO
            THELP   = ZERO
            DO LAY = 1, N-1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              L_TPROP = L_T_WLOWER(I,LAY,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              TPROP = T_WLOWER(I,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTDN(I,UT)
            L_TPROP = L_UT_T_PARTIC(I,UT,Q) / QUAD_STREAMS(I)
            L_THELP = L_THELP + L_TPROP + THELP * L_T_DISORDS_UTDN(I,UT,Q)
            QATMOSWF_F(Q,UTA,I,O1) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  ###########################################

!  Homogeneous
!  -----------

!  stream and Stokes directions

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

!  Parameter loop

          DO Q = 1, NV_PARAMETERS

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
              MXR  = MCON(K,N)   *   SOLB_XNEG(I,O1,K,N)
              LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
              MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
              NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              HOM1 = ( NXR + LLXR ) *   T_UTDN_EIGEN(K,UT) &
                           +  LXR   * L_T_UTDN_EIGEN(K,UT,Q)
              HOM2 = ( PXR + MLXR ) *   T_UTUP_EIGEN(K,UT) &
                           +  MXR   * L_T_UTUP_EIGEN(K,UT,Q)
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                      - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                      + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                      - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                      + PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)

              LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N) &
                      - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
              LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N) &
                      + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)
              MXR1  =   MCON(K1,N) *   SOLB_XNEG(I,O1,K1,N) &
                      - MCON(K2,N) *   SOLB_XNEG(I,O1,K2,N)
              MXR2  =   MCON(K1,N) *   SOLB_XNEG(I,O1,K2,N) &
                      + MCON(K2,N) *   SOLB_XNEG(I,O1,K1,N)

              LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q) &
                       - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
              LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q) &
                       + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
              MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q) &
                       - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)
              MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K2,N,Q) &
                       + MCON(K2,N) * L_SOLB_XNEG(I,O1,K1,N,Q)

              HOM1CR =   ( NXR1 + LLXR1 ) *   T_UTDN_EIGEN(K1,UT) &
                       - ( NXR2 + LLXR2 ) *   T_UTDN_EIGEN(K2,UT)
              HOM2CR =             LXR1   * L_T_UTDN_EIGEN(K1,UT,Q) &
                                 - LXR2   * L_T_UTDN_EIGEN(K2,UT,Q)
              HOM3CR =   ( PXR1 + MLXR1 ) *   T_UTUP_EIGEN(K1,UT) &
                       - ( PXR2 + MLXR2 ) *   T_UTUP_EIGEN(K2,UT)
              HOM4CR =             MXR1   * L_T_UTUP_EIGEN(K1,UT,Q) &
                                 - MXR2   * L_T_UTUP_EIGEN(K2,UT,Q)

              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR + HOM4CR
            ENDDO

!  real part

            SHOM = SHOM_R + SHOM_CR
            QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

!  Finish Q, I and O1 loops

          ENDDO
        ENDDO
      ENDDO

!  Add the linearized thermal solution  (if flagged)
!   THIS IS THE SOLUTION with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            SPAR = L_UT_T_PARTIC(I,UT,Q)
            QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
          ENDDO
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Add Linearized Classical beam solution
!  ======================================

!   -- 1/31/21. Section now has the classical flag returned

! @@@@@@@@@@ Rob fix, 2/9/11, Replacement of commented section
!      IF ( DO_CLASSICAL_SOLUTION ) THEN
!        DO I = 1, NSTREAMS
!          DO O1 = 1, NSTOKES
!            DO Q = 1, NV_PARAMETERS
!              SPAR = L_WUPPER(I,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
!                       WUPPER(I,O1,N)   * LC_T_UTDN_MUBAR(UT,IB,Q)
!              QATMOSWF_F(Q,UTA,I,O1) = &
!                QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
!            ENDDO
!          ENDDO
!        ENDDO
!      ELSE
!                P L A C E H O L D E R
!      ENDIF
! @@@@@@@@@@ Rob fix, 2/9/11, End Replaced section

      IF ( DO_CLASSICAL ) THEN
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SUNP   = WUPPER(I,O1,N)
                L_SUNP = L_WUPPER(I,O1,N,Q)
                IF ( O1 .eq. 1 ) THEN
                  SUNP   = SUNP - T_WUPPER(I,N)
                  L_SUNP = L_SUNP - L_T_WUPPER(I,N,Q)
                ENDIF
                SPAR = L_SUNP * T_UTDN_MUBAR(UT,IB) + SUNP * LC_T_UTDN_MUBAR(UT,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = L_WUPPER(I,O1,N,Q) * T_UTDN_MUBAR(UT,IB) + WUPPER(I,O1,N) * LC_T_UTDN_MUBAR(UT,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Linearized Green's function solution
! =====================================

!  1/31/21. Version 2.8.3. COMPLETELY NEW SECTION
!     -- Follows the LIDORT implementation closely

     IF ( .not. DO_CLASSICAL ) THEN

!  Sum up the Green's function contributions
!  -----------------------------------------

         DO O1 = 1, NSTOKES
           DO I = 1, NSTREAMS
             DO Q = 1, NV_PARAMETERS
               SPAR = DOT_PRODUCT (    UT_GMULT_UP(1:NSTKS_NSTRMS,UT),   L_SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N,Q) ) &
                    + DOT_PRODUCT ( LC_UT_GMULT_UP(1:NSTKS_NSTRMS,UT,Q),   SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N)   ) &
                    + DOT_PRODUCT (    UT_GMULT_DN(1:NSTKS_NSTRMS,UT),   L_SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N,Q) ) &
                    + DOT_PRODUCT ( LC_UT_GMULT_DN(1:NSTKS_NSTRMS,UT,Q),   SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N)   )
               IF ( O1.gt.2) SPAR = -SPAR
               QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1)  + FMULT * SPAR
             ENDDO
           ENDDO
         ENDDO

!  End Green's function solution

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADCOLUMNWF_OFFGRID_DN

!

      END MODULE vlidort_lc_PostProcessing_m
