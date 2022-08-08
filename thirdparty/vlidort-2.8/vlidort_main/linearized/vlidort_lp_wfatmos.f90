
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
! #     Top level routines--------------                   #
! #                                                        #
! #       UPUSER_PROFILEWF                                 #
! #       DNUSER_PROFILEWF                                 #
! #                                                        #
! #       ------Calling                                    #
! #                                                        #
! #            GET_LP_TOASOURCE                            #
! #            GET_LP_BOASOURCE                            #
! #                                                        #
! #       MIFLUX_PROFILEWF                                 #
! #                                                        #
! ##########################################################

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!  1/31/21. Version 2.8.3. Introduce MSST facility (linearizations)
!  1/31/21. Version 2.8.3. Quadrature and sourceterm routines moved to new module LP_POSTPROCESSING
!  1/31/21. Version 2.8.3. LPS Convergence routines moved to new module LPS_CONVERGE

      MODULE vlidort_lp_wfatmos_m

!  1/31/21. Version 2.8.3. Dependency for Quadrature and sourceterm routines

      USE vlidort_lp_PostProcessing_m

!   1/31/21. Version 2.8.3. Only the  TOA/BOA source routines are now private

      PRIVATE :: GET_LP_TOASOURCE,        GET_LP_BOASOURCE

      PUBLIC  :: UPUSER_PROFILEWF,     &
                 DNUSER_PROFILEWF,     &
                 MIFLUX_PROFILEWF

      CONTAINS

      SUBROUTINE UPUSER_PROFILEWF ( &
        DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,             & ! Input flags (RT mode)
        DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,            & ! Input flags (sources)
        DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT,    DO_INCLUDE_SURFACE,              & ! Input flags (sources)
        DO_LAMBERTIAN_SURFACE, DO_DBCORR, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,   & ! Input flags (Surface)
        DO_LAYER_SCATTERING, DO_MSSTS, FOURIER, IBEAM, ECHT_NSTOKES, NSTOKES, NSTREAMS, NLAYERS,    & ! Input flags/numbers (basic)
        N_PPSTREAMS, PPSTREAM_MASK, N_USER_LEVELS, LAYER_TO_VARY, NV_PARAMETERS,      & ! Input numbers (basic)
        LEVEL_MASK_UP, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input Level output control
        TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,          & ! Input Taylor/Optical
        FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS, USER_SECANTS, MUELLER_INDEX,     & ! Input Flux/quad/Mueller
        BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_DISORDS,       & ! Input Beam parameterization
        T_DELT_USERM, T_UTUP_USERM, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,      & ! Input UTrans/Surface.
        RF_USER_DIRECT_BEAM, SL_USERTERM, LP_TRANS_ATMOS_FINAL, CUMSOURCE_UP,         & ! Input surfradiances/Cumsource
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,            & ! Input Homogeneous Solutions
        ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_UU, PMULT_UD, SIGMA_P,        & ! Input Green's function variables
        T_WLOWER, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, HMULT_1, HMULT_2,       & ! Input Thermal/Homog
        EMULT_UP, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, UT_PMULT_UU, UT_PMULT_UD,    & ! Input multipliers  
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, L_T_DELT_DISORDS,          & ! Linearized BeamParm/DODTrans
        L_T_DELT_USERM, L_T_UTUP_USERM, L_T_WLOWER, L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP,  & ! Linearized Trans/Thermal
        L_KEIGEN, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG, L_UHOM_UPDN, L_UHOM_UPUP,    & ! Linearized Homog solutions
        L_HMULT_1, L_HMULT_2, LP_EMULT_UP, L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP, & ! Linearized multipliers
        L_UPAR_UP_1, LP_UPAR_UP_2, L_ATERM_SAVE, L_BTERM_SAVE, L_WLOWER, NCON, PCON,     & ! Linearized PI Solutions
        L_BOA_THTONLY_SOURCE, PROFILEWF_F, LP_LAYER_MSSTS_F, LP_SURF_MSSTS_F )             ! Output

!   Streamlined for Version 2.8. 7/6/16

!  1/31/21. Version 2.8.3.
!    ==> New input arguments to the post processing (Green's function) : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE, SIGMA_P
!           ** Green's function arrays ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_UU, PMULT_UD
!    ==> Other Input/Output arguments and other changes
!           ** for the multiple scatter source term linearization, use control flag DO_MSSTS
!           ** LInearizations of MSST functions (LP_LAYER_MSSTS_F, LP_SURF_MSSTS_F), now output
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!    -- (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES) number of STOKES components

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_DIRECTIONS,   &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, UPIDX

      IMPLICIT NONE

!  Inputs
!  ======

!  Input RT flags

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT

!  1/31/21. Version 2.8.3. Add DO_CLASSICAL_SOLUTION flag

      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION

!  Source flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

!  Surface flags

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_DBCORR

!  4/9/19. Replacement of DIRECTBEAM by two separate flags
      
      LOGICAL, intent(in)  ::           DO_INCLUDE_DIRECTRF
      LOGICAL, intent(in)  ::           DO_INCLUDE_DIRECTSL
!      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM

!  Layer scattering

      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  removed, Version 2.8, 7.17.16
!      LOGICAL, INTENT (IN) ::           DO_QUAD_OUTPUT

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::           DO_MSSTS

!  numbers
!    -- 1/31/21. Version 2.8.3. (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES)

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           ECHT_NSTOKES, NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           LAYER_TO_VARY, NV_PARAMETERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Level output control + Partial-layer control

      INTEGER, INTENT (IN) ::           LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3. Order of Taylor series (including terms up to EPS^n)

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER

!  1/31/21. Version 2.8.3.  Local vertical optical depths and linearization

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS ) 
      DOUBLE PRECISION, intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Flux multipliers and quadrature, secants

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, intent (in) ::  USER_SECANTS  ( MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX  ( MAXSTOKES, MAXSTOKES )

!  Main solution inputs
!  --------------------

!  1/31/21. Version 2.8.3.  Beam average-secant parameterization
!    BEAM_CUTOFF = Last layer to include Particular integral solution

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )

!  User-stream and discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM   ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Cumulative source term

      DOUBLE PRECISION, INTENT (IN) ::  CUMSOURCE_UP   ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

!  Surface reflectance variables
!    -- 1/31/21. Versdion 2.8.3. BRDF arrays are locallyl deifned for each Fourier, drop MAXMOMENTS

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR, ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F       ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F  ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  Reflected Direct beam and surface-leaving solutions

      DOUBLE PRECISION, INTENT (IN) ::  RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SL_USERTERM         ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION, INTENT (IN) ::  USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Linearized transmittance flux for water-leaving

      DOUBLE PRECISION, INTENT (IN) ::  LP_TRANS_ATMOS_FINAL (MAXBEAMS,MAXLAYERS,MAX_ATMOSWFS)
      
!  RTE solutions (Homog.)

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON    ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON    ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  User stream RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Particular integral user solutions and beam multipliers
!    -- 1/31/21. Version 2.8.3. SIGMA_P is new addition
!  3/1/21. Bug Fix on SIGMA_P declaration

      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_1   ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_2   ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP    ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIGMA_P     ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. New Input: Green's function solution variables

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  1/31/21. Version 2.8.3. New Input: Green's function Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) ::  PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) ::  UT_PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (IN) ::  UT_PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized inputs
!  -----------------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!     -- 1/31/21. Version 2.8.3. New additions for the Green's function solution

      DOUBLE PRECISION, intent(in)  ::  LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  ::  LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  ::  LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized transmittances

      DOUBLE PRECISION, INTENT (IN) ::  L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous and sdiscrete ordinate solution
!     -- 1/31/21. Version 2.8.3. L_KEIGEN added for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) ::  L_KEIGEN       ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized thermal solutions

      DOUBLE PRECISION, INTENT (IN) ::  L_T_WLOWER        ( MAXSTREAMS_2,     MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_LAYER_TSUP_UP   ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  linearized User homogeneous and beam solutions

      DOUBLE PRECISION, INTENT (IN) ::  L_UHOM_UPDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_UHOM_UPUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_UPAR_UP_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LP_UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Homog. solution multipliers

      DOUBLE PRECISION, INTENT (IN) ::  L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Beam solution multipliers

      DOUBLE PRECISION, INTENT (IN) ::  LP_EMULT_UP    ( MAX_USER_STREAMS, MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3.  New Input: Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) ::  L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  OUTPUT
!  ======

!  Linearized BOA source term for Thermal-transonly

      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_THTONLY_SOURCE ( MAXSTREAMS, MAX_ATMOSWFS )

!  Linearized output (Ostensibly output)

      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                        MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 1/31/21. Version 2.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LP_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LP_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  local variables
!  ---------------

!  Help

      LOGICAL ::          SFLAG
      INTEGER ::          M, N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, NV, NVPARS, Q, NC, UT, IB, LUM
      DOUBLE PRECISION :: L_FINAL_SOURCE

!  Local arrays

      DOUBLE PRECISION :: L_CUMUL_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_BOA_MSSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_BOA_DBSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  Proxies

      NV     = LAYER_TO_VARY
      NVPARS = NV_PARAMETERS
      M   = FOURIER
      IB  = IBEAM

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!     -- 1/31/21. Version 2.8.3. Use post-processing mask input
!    -- 1/31/21. Version 2.8.3.  (RTS 2/16/21). Zero with the true NSTOKES value (ECHT_NSTOKES)
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present); rearranged order of loops

      !IF ( DO_USER_STREAMS ) THEN
      !  DO UTA = 1, N_USER_LEVELS
      !    DO Q = 1, NV_PARAMETERS
      !      DO LUM = 1, N_PPSTREAMS
      !        UM = PPSTREAM_MASK(LUM,IB)
      !        DO O1 = 1, ECHT_NSTOKES
      !          PROFILEWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = ZERO
      !        ENDDO
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ENDIF

      IF ( DO_USER_STREAMS ) THEN
        DO O1 = 1, ECHT_NSTOKES
          DO LUM = 1, N_PPSTREAMS
            DO UTA = 1, N_USER_LEVELS
              DO Q = 1, NV_PARAMETERS
                PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,UPIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. ==> Zero MSST outputs
!rob & mick fix 1/5/2021 - added this zeroing

      DO Q = 1, NV_PARAMETERS
        DO O1 = 1, ECHT_NSTOKES
          LP_LAYER_MSSTS_F(IB,O1,1:NLAYERS,NV,Q) = ZERO
          LP_SURF_MSSTS_F (IB,O1,NV,Q)           = ZERO
        ENDDO
      ENDDO

!  Initialize post-processing recursion
!  ====================================

      IF ( DO_USER_STREAMS ) THEN

!  Get the linearized BOA source terms (diffuse and direct)
!     -- 1/31/21. Version 2.8.3. Use post-processing mask input, define BRDF arrays locally.
!mick fix 1/5/2021 - changed arg order of IBEAM & FOURIER based on GET_LP_BOASOURCE subroutine statement

        CALL GET_LP_BOASOURCE ( &
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,          & ! Input flags sources
          DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_LAMBERTIAN_SURFACE, DO_DBCORR, & ! Input flags misc
          DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,           & ! Input flags surface
          NSTOKES, NSTREAMS, NLAYERS, IBEAM, FOURIER, NV, NV_PARAMETERS,          & ! Input Numbers/Indices
          N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS, MUELLER_INDEX,  & ! input bookkeeping/quadrature
          DELTAU_SLANT, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,              & ! Input Surface
          RF_USER_DIRECT_BEAM, SL_USERTERM, LP_TRANS_ATMOS_FINAL,                 & ! Input surface radiance
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                  & ! Input Homog solutions
          LCON, MCON, T_DELT_DISORDS, T_WLOWER,                                   & ! Input Thermal and PI
          L_DELTAU_VERT, L_T_DELT_DISORDS, L_T_WLOWER,                            & ! Linearized thermal and misc.
          L_SOLA_XPOS, L_SOLB_XNEG, L_T_DELT_EIGEN, L_WLOWER, NCON, PCON,         & ! Linearized DsOr solutions
          L_BOA_THTONLY_SOURCE, L_BOA_MSSOURCE, L_BOA_DBSOURCE )                    ! OUTPUT BOA terms (linearized)

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST surface source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!

        IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) then
          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              LP_SURF_MSSTS_F(IB,O1,NV,Q) = FLUX_MULTIPLIER * L_BOA_MSSOURCE(IB,O1,Q)
            ENDDO
          ENDDO
        ENDIF

!  Set the cumulative source term equal to the BOA sum
!     -- 1/31/21. Version 2.8.3. Use post-processing mask input

        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              L_CUMUL_SOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) + L_BOA_DBSOURCE(UM,O1,Q)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            SFLAG = DO_LAYER_SCATTERING(FOURIER,N)
            NC = NLAYERS + 1 - N

!  1/31/21. Version 2.8.3. Several Changes to this Call.
!    -- Introduce control flag DO_CLASSICAL_SOLUTION, and Taylor-series TAYLOR_ORDER
!    -- Introduce Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK
!    -- Introduce Green's function inputs ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UPAR_UP_2, PMULT_UU, PMULT_UD
!    -- Introduce Other inputs LP_AVERAGE_SECANT, SIGMA_P, L_ATERM_SAVE, L_BTERM_SAVE

            CALL LP_WHOLELAYER_STERM_UP ( &
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,        & ! Input flags
             DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SFLAG, M,                   & ! Input flags/Fourier
             N, IB, NV, NVPARS, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK,               & ! Input numbers
             TAYLOR_ORDER, DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM, & ! input Optical + User_streams
             BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, EMULT_UP, SIGMA_P,          & ! input Beam stuff
             K_REAL, K_COMPLEX, LCON, MCON, UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2, UPAR_UP_1,   & ! input Homogeneous/Beam-1
             ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UPAR_UP_2, PMULT_UU, PMULT_UD,            & ! Input Greens/Beam-2
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_EMULT_UP,                  & ! input Linearized Beam stuff
             L_KEIGEN, NCON, PCON, L_UHOM_UPDN, L_UHOM_UPUP, L_HMULT_1, L_HMULT_2,               & ! input Linearized Homogeneous
             L_UPAR_UP_1, LP_UPAR_UP_2, L_ATERM_SAVE, L_BTERM_SAVE, L_LAYER_TSUP_UP,             & ! Input linearized PIs
             L_LAYER_SOURCE )                                                         ! OUTPUT

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY) THEN
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  LP_LAYER_MSSTS_F(IB,O1,N,NV,Q) = FLUX_MULTIPLIER * L_LAYER_SOURCE(IB,O1,Q)
                ENDDO
              ENDDO
            ENDIF

!  1/31/21. Version 2.8.3. Cumulative sourceterm, use Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK

            IF ( N.EQ.NV ) THEN
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                                                                      + L_T_DELT_USERM(N,UM,Q) * CUMSOURCE_UP(UM,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV ) THEN
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) + T_DELT_USERM(N,UM) * L_CUMUL_SOURCE(UM,O1,Q)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  End collection of layer sources

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT    = PARTLAYERS_OUTINDEX(UTA)
          N     = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

!  1/31/21. Version 2.8.3. Several Changes.
!    -- Introduce control flag DO_CLASSICAL_SOLUTION, and Taylor-series TAYLOR_ORDER
!    -- Introduce Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK
!    -- Introduce Green's function inputs ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UT_PMULT_UU, UT_PMULT_UD
!    -- Introduce Other inputs SIGMA_P, LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE

            CALL LP_PARTLAYER_STERM_UP ( &
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                & ! Input flags
              DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SFLAG, M, IB,                       & ! Input flags/numbers
              N, UT, NV, NV_PARAMETERS, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,  & ! Input numbers
              PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTUP_USERM,          & ! Input Optical/Usertrans
              L_T_UTUP_USERM, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,       & ! input Beam stuff
              K_REAL, K_COMPLEX, LCON, MCON, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1,               & ! Input user-solutions
              UPAR_UP_2, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, SIGMA_P,                    & ! Input User-solutions/Greens
              ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UT_PMULT_UU, UT_PMULT_UD,           & ! Input Greens
              LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_UT_EMULT_UP,         & ! input Linearized Beam stuff
              L_KEIGEN, NCON, PCON, L_UHOM_UPDN, L_UHOM_UPUP, L_UT_HMULT_UU, L_UT_HMULT_UD, & ! input Linearized Homogeneous
              L_UPAR_UP_1, LP_UPAR_UP_2, L_ATERM_SAVE, L_BTERM_SAVE, L_LAYER_TSUP_UTUP,     & ! Input linearized PIs
              L_LAYER_SOURCE )                                                                ! OUTPUT

!  1/31/21. Version 2.8.3. Cumulative sourceterm, use Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK

            IF ( N.EQ.NV ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO O1 = 1, NSTOKES
                     DO Q = 1, NV_PARAMETERS
                        L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) &
                                            +   T_UTUP_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                                            + L_T_UTUP_USERM(UT,UM,Q) *   CUMSOURCE_UP(UM,O1,NC)
                        PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                     ENDDO
                  ENDDO
               ENDDO
            ELSE IF ( N.NE.NV ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO O1 = 1, NSTOKES
                     DO Q = 1, NV_PARAMETERS
                        L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) + T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,O1,Q)
                        PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,UPIDX) =  FLUX_MULTIPLIER * L_FINAL_SOURCE
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

!  End user stream clause

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term
!    -- 1/31/21. Version 2.8.3. use Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO O1 = 1, NSTOKES
                   DO Q = 1, NV_PARAMETERS
                      L_FINAL_SOURCE =  FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,O1,Q)
!                    if (DABS(L_FINAL_SOURCE).GT.1.0d-12 ) then
                      PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,UPIDX) = L_FINAL_SOURCE
!                    endif
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

!  end level output assignation

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE UPUSER_PROFILEWF

!

      SUBROUTINE DNUSER_PROFILEWF ( &
        DO_USER_STREAMS,   DO_OBSERVATION_GEOMETRY, DO_CLASSICAL_SOLUTION,              & ! Input flags (RT mode)
        DO_SOLAR_SOURCES,  DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,               & ! Input flags (sources)
        DO_MSMODE_VLIDORT, DO_LAYER_SCATTERING,     DO_MSSTS,  FOURIER, IBEAM,          & ! Input flags (sources)
        ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, NV, NV_PARAMETERS, N_PPSTREAMS, PPSTREAM_MASK,& ! Input flags/numbers (basic)
        LEVEL_MASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,    & ! Input Level output control
        TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT, FLUX_MULT, USER_SECANTS, & ! Input Taylor/Optical/Secants
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,         & ! Input Beam parameterization
        T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN,  K_REAL, K_COMPLEX, LCON, MCON,       & ! Input User/Cumsource/BVP
        ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_DU, PMULT_DD, SIGMA_M,          & ! Input Green's function variables
        UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2, HMULT_1, HMULT_2,                   & ! Homog
        EMULT_DN, UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, UT_PMULT_DU, UT_PMULT_DD,      & ! Input multipliers  
        L_T_DELT_USERM, L_T_UTDN_USERM, LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Linearized BeamParm/Trans
        L_KEIGEN, L_UHOM_DNDN, L_UHOM_DNUP, L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN,               & ! Linearized Homog/Thermal
        L_HMULT_1, L_HMULT_2, LP_EMULT_DN, L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,      & ! Linearized multipliers
        L_UPAR_DN_1, LP_UPAR_DN_2, L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON,                    & ! Linearized PI Solutions
        PROFILEWF_F, LP_LAYER_MSSTS_F )                                                         ! Output

!   Streamlined for Version 2.8. 7/6/16

!  1/31/21. Version 2.8.3.
!    ==> New input arguments to the post processing (Green's function) : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE, SIGMA_M
!           ** Green's function arrays ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_DU, PMULT_DD
!    ==> Other Input/Output arguments and other changes
!           ** for the multiple scatter source term linearization, use control flag DO_MSSTS
!           ** LInearizations of MSST functions (LP_LAYER_MSSTS_F), now output
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!    -- (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES) number of STOKES components

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_DIRECTIONS,   &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, DNIDX

      IMPLICIT NONE

!  Input RT flags

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT

!  1/31/21. Version 2.8.3. Add DO_CLASSICAL_SOLUTION flag

      LOGICAL, INTENT (IN) ::           DO_CLASSICAL_SOLUTION

!  Source flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

!  Scattering

      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  removed, Version 2.8, 7.17.16
!      LOGICAL, INTENT (IN) ::           DO_QUAD_OUTPUT

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::           DO_MSSTS

!  numbers
!    -- 1/31/21. Version 2.8.3. (RTS 2/16/21). Use actual (ECHT_NSTOKES) and local (NSTOKES)

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           ECHT_NSTOKES, NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           NV, NV_PARAMETERS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Level output control + Partial-layer control

      INTEGER, INTENT (IN) ::           LEVEL_MASK_DN       ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3. Order of Taylor series (including terms up to EPS^n)

      INTEGER, INTENT (IN) ::           TAYLOR_ORDER

!  1/31/21. Version 2.8.3.  Local vertical optical depths and linearization

      DOUBLE PRECISION, intent(in)  ::  PARTAU_VERT   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, intent(in)  ::  DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, intent(in)  ::  L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Flux multipliers and secants

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_MULT
      DOUBLE PRECISION, intent (in) ::  USER_SECANTS  ( MAX_USER_STREAMS )

!  Main solution inputs
!  --------------------

!  1/31/21. Version 2.8.3.  Beam average-secant parameterization
!    BEAM_CUTOFF = Last layer to include Particular integral solution

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  User-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM   ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Cumulative source term

      DOUBLE PRECISION, INTENT (IN) ::  CUMSOURCE_DN   ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

!  RTE solutions (Homog.)

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON    ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON    ( MAXSTRMSTKS, MAXLAYERS )

!  User stream RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Particular integral user solutions and beam multipliers
!    -- 1/31/21. Version 2.8.3. SIGMA_M is new addition
!  3/1/21. Bug Fix on SIGMA_M declaration

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  EMULT_DN    ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) ::  SIGMA_M     ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. New Input: Green's function solution variables

      DOUBLE PRECISION, INTENT (IN) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_P    ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GAMMA_M    ( MAXEVALUES, MAXLAYERS )

!  1/31/21. Version 2.8.3. New Input: Green's function Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) ::  PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) ::  UT_PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (IN) ::  UT_PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized inputs
!  -----------------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!     -- 1/31/21. Version 2.8.3. New additions for the Green's function solution

      DOUBLE PRECISION, intent(in)  ::  LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  ::  LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  ::  LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized transmittances

      DOUBLE PRECISION, INTENT (IN) ::  L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues and BVP constants

      DOUBLE PRECISION, INTENT (IN) ::  L_KEIGEN  ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized thermal solutions

      DOUBLE PRECISION, INTENT (IN) ::  L_LAYER_TSUP_DN   ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  linearized User homogeneous and beam solutions

!     -- 1/31/21. Version 2.8.3. L_KEIGEN added for the Green's function solution

      DOUBLE PRECISION, INTENT (IN) ::  L_UHOM_DNDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_UHOM_DNUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_UPAR_DN_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LP_UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Homog. solution multipliers

      DOUBLE PRECISION, INTENT (IN) ::  L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Beam solution multipliers

      DOUBLE PRECISION, INTENT (IN) ::  LP_EMULT_DN    ( MAX_USER_STREAMS, MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3.  New Input: Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION, INTENT (IN) ::  L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) ::  L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Output
!  ======

!  Linearized output (Ostensibly output)

      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                        MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 1/31/21. Version 2.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LP_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  local variables
!  ---------------

!  help

      LOGICAL ::          SFLAG
      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, NC, UT, IB, M, LUM
      DOUBLE PRECISION :: L_FINAL_SOURCE

!  Local arrays

      DOUBLE PRECISION :: L_CUMUL_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_TOA_SOURCE  ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  Proxies

      IB  = IBEAM
      M   = FOURIER

!  Zero all Fourier component output
!     -- 1/31/21. Version 2.8.3. Use post-processing mask input
!    -- 1/31/21. Version 2.8.3.  (RTS 2/16/21). Zero with the true NSTOKES value (ECHT_NSTOKES)
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present); rearranged order of loops

      !IF ( DO_USER_STREAMS ) THEN
      !  DO UTA = 1, N_USER_LEVELS
      !    DO Q = 1, NV_PARAMETERS
      !      DO LUM = 1, N_PPSTREAMS
      !        UM = PPSTREAM_MASK(LUM,IB)
      !        DO O1 = 1, ECHT_NSTOKES
      !          PROFILEWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = ZERO
      !        ENDDO
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ENDIF

      IF ( DO_USER_STREAMS ) THEN
        DO O1 = 1, ECHT_NSTOKES
          DO LUM = 1, N_PPSTREAMS
            DO UTA = 1, N_USER_LEVELS
              DO Q = 1, NV_PARAMETERS
                PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,DNIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. ==> Zero MSST outputs
!rob & mick fix 1/5/2021 - added this zeroing

      DO Q = 1, NV_PARAMETERS
        DO O1 = 1, ECHT_NSTOKES
          LP_LAYER_MSSTS_F(IB,O1,1:NLAYERS,NV,Q) = ZERO
        ENDDO
      ENDDO

!  Initialize post-processing recursion
!  ====================================

      IF ( DO_USER_STREAMS ) THEN

!  Get the linearized TOA source terms
!     -- 1/31/21. Version 2.8.3. Use post-processing mask input

        CALL GET_LP_TOASOURCE ( &
           N_PPSTREAMS, PPSTREAM_MASK, NSTOKES, IBEAM, NV_PARAMETERS, LP_TOA_SOURCE )

!  Set the cumulative source term equal to the BOA sum
!     -- 1/31/21. Version 2.8.3. Use post-processing mask input

        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              L_CUMUL_SOURCE(UM,O1,Q) = LP_TOA_SOURCE(UM,O1,Q)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            SFLAG = DO_LAYER_SCATTERING(FOURIER,N)
            NC = N

!  1/31/21. Version 2.8.3. Several Changes to this Call.
!    -- Introduce control flag DO_CLASSICAL_SOLUTION, and Taylor-series TAYLOR_ORDER
!    -- Introduce Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK
!    -- Introduce Green's function inputs ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UPAR_DN_2, PMULT_UU, PMULT_UD
!    -- Introduce Other inputs LP_AVERAGE_SECANT, SIGMA_M, L_ATERM_SAVE, L_BTERM_SAVE

            CALL LP_WHOLELAYER_STERM_DN ( &
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,        & ! Input flags
             DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SFLAG, FOURIER,             & ! Input flags/Fourier
             N, IB, NV, NV_PARAMETERS, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK,        & ! Input numbers
             TAYLOR_ORDER, DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM, & ! input Optical + User_streams
             BEAM_CUTOFF, AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR, EMULT_DN, SIGMA_M,        & ! input Beam stuff
             K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, HMULT_1, HMULT_2, UPAR_DN_1,   & ! input Homogeneous/Beam-1
             ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UPAR_DN_2, PMULT_DU, PMULT_DD,            & ! Input Greens/Beam-2
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_EMULT_DN,                  & ! input Linearized Beam stuff
             L_KEIGEN, NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_HMULT_1, L_HMULT_2,               & ! input Linearized Homogeneous
             L_UPAR_DN_1, LP_UPAR_DN_2, L_ATERM_SAVE, L_BTERM_SAVE, L_LAYER_TSUP_DN,             & ! Input linearized PIs
             L_LAYER_SOURCE )                                                                      ! OUTPUT

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY) THEN
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  LP_LAYER_MSSTS_F(IB,O1,N,NV,Q) = FLUX_MULT * L_LAYER_SOURCE(IB,O1,Q)
                ENDDO
              ENDDO
            ENDIF

!  1/31/21. Version 2.8.3. Cumulative sourceterm, use Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK

            IF ( N.EQ.NV ) THEN
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                                                                      + L_T_DELT_USERM(N,UM,Q) * CUMSOURCE_DN(UM,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV ) THEN
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) + T_DELT_USERM(N,UM) * L_CUMUL_SOURCE(UM,O1,Q)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  End collection of layer sources

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT    = PARTLAYERS_OUTINDEX(UTA)
          N     = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

!  1/31/21. Version 2.8.3. Several Changes.
!    -- Introduce control flag DO_CLASSICAL_SOLUTION, and Taylor-series TAYLOR_ORDER
!    -- Introduce Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK
!    -- Introduce Green's function inputs ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UT_PMULT_UU, UT_PMULT_UD
!    -- Introduce Other inputs SIGMA_M, LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE

            CALL LP_PARTLAYER_STERM_DN ( &
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                & ! Input flags
             DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, SFLAG, M, IB,                       & ! Input flags/numbers
             N, UT, NV, NV_PARAMETERS, NSTOKES, N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,  & ! Input numbers
             PARTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTDN_USERM,                       & ! Input Optical/Usertrans
             BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,                       & ! input Beam stuff
             K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, UPAR_DN_1,               & ! Input user-solutions
             UPAR_DN_2, UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, SIGMA_M,                    & ! Input User-solutions/Greens
             ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, UT_PMULT_DU, UT_PMULT_DD,           & ! Input Greens
             LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_UT_EMULT_DN,         & ! input Linearized Beam stuff
             L_KEIGEN, NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_UT_HMULT_DU, L_UT_HMULT_DD, & ! input Linearized Homogeneous
             L_UPAR_DN_1, LP_UPAR_DN_2, L_ATERM_SAVE, L_BTERM_SAVE, L_LAYER_TSUP_UTDN,     & ! Input linearized PIs
             L_LAYER_SOURCE )                                                                ! OUTPUT

!  1/31/21. Version 2.8.3. Cumulative sourceterm, use Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK

            IF ( N.EQ.NV ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO O1 = 1, NSTOKES
                     DO Q = 1, NV_PARAMETERS
                        L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) &
                                            +   T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                                            + L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,O1,NC)
                        PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,DNIDX) = FLUX_MULT * L_FINAL_SOURCE
                     ENDDO
                  ENDDO
               ENDDO
            ELSE IF ( N.NE.NV ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO O1 = 1, NSTOKES
                     DO Q = 1, NV_PARAMETERS
                        L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) + T_UTDN_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,O1,Q)
                        PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,DNIDX) =  FLUX_MULT * L_FINAL_SOURCE
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

!  End user stream clause

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term
!    -- 1/31/21. Version 2.8.3. use Post-processing Mask N_PPSTREAMS, PPSTREAM_MASK

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO O1 = 1, NSTOKES
                   DO Q = 1, NV_PARAMETERS
                      L_FINAL_SOURCE =  FLUX_MULT * L_CUMUL_SOURCE(UM,O1,Q)
                      PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,DNIDX) = L_FINAL_SOURCE
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

!  End level output clause

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE DNUSER_PROFILEWF

!

      SUBROUTINE GET_LP_TOASOURCE ( &
           N_PPSTREAMS, PPSTREAM_MASK, NSTOKES, IBEAM, NV_PARAMETERS, LP_TOA_SOURCE )

!  1/31/21. Version 2.8.3. 
!    -- Use post-processing mask N_PPSTREAMS, PPSTREAM_MASK, Drop OBSERVATION geometry flag.

!  Linearized Top of the atmosphere source term

!  module, dimensions and numbers
!    -- 1/31/21. Version 2.8.3. Add MAXBEAMS

      USE VLIDORT_PARS_m, only : MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Beam, Stokes

      INTEGER  , intent(in)  :: IBEAM, NSTOKES

!  Linearization

      INTEGER  , intent(in)  :: NV_PARAMETERS

!  Output

      DOUBLE PRECISION, INTENT (OUT) :: LP_TOA_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER :: LUM, UM, Q

!  initialise TOA source function
!  ------------------------------

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IBEAM)
         DO Q = 1, NV_PARAMETERS
            LP_TOA_SOURCE(UM,1:NSTOKES,Q) = ZERO
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE GET_LP_TOASOURCE

!

      SUBROUTINE GET_LP_BOASOURCE ( &
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,          & ! Input flags sources
          DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_LAMBERTIAN_SURFACE, DO_DBCORR, & ! Input flags misc
          DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,           & ! Input flags surface
          NSTOKES, NSTREAMS, NLAYERS, IBEAM, FOURIER, NV, NV_PARAMETERS,          & ! Input numbers/indices
          N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS, MUELLER_INDEX,  & ! input bookkeeping/quadrature
          DELTAU_SLANT, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,      & ! Input Surface
          RF_USER_DIRECT_BEAM, SL_USERTERM, LP_TRANS_ATMOS_FINAL,         & ! Input surface radiance
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,          & ! Input Homog solutions
          LCON, MCON, T_DELT_DISORDS, T_WLOWER,                           & ! Input Thermal and PI
          L_DELTAU_VERT, L_T_DELT_DISORDS, L_T_WLOWER,                    & ! Linearized thermal and misc.
          L_SOLA_XPOS, L_SOLB_XNEG, L_T_DELT_EIGEN, L_WLOWER, NCON, PCON, & ! Linearized DsOr solutions
          L_BOA_THTONLY_SOURCE, L_BOA_MSSOURCE, L_BOA_DBSOURCE )            ! OUTPUT BOA terms (linearized)

!  Linearized Bottom of the atmosphere source term

!  1/31/21. Version 2.8.3. Ordering of arguments changed. 
!    -- Use of Postprocessing mask system
!    -- BRDF arrays are defined locally, no Fourier index, drop MAXMOMENTS

!  module, dimensions and numbers
!    -- 1/31/21. Version 2.8.3. Add MAXBEAMS

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE

      IMPLICIT NONE
      
!  INPUTS
!  ======

!  Flags sources

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

!  Flags RT control

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_MVOUTPUT

!  Flags surface

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_DBCORR

!  4/9/19. Replacement of DIRECTBEAM by two separate flags
      
      LOGICAL, intent(in)  ::           DO_INCLUDE_DIRECTRF
      LOGICAL, intent(in)  ::           DO_INCLUDE_DIRECTSL
!      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      
!  Numbers

      INTEGER, INTENT (IN) ::           NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

!  1/31/21. Version 2.8.3. Post-processing mask

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Bookkeeping (quadrature etc)

      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS ) 
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Surface reflectance
!    -- 1/31/21. Version 2.8.3. BRDF arrays defined locally, drop MAXMOMENTS index

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F      ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  Reflected Direct beam and surface-leaving solutions

      DOUBLE PRECISION, INTENT (IN) ::  RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SL_USERTERM         ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION, INTENT (IN) ::  USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Linearized transmittance flux for water-leaving

      DOUBLE PRECISION, INTENT (IN) ::  LP_TRANS_ATMOS_FINAL (MAXBEAMS,MAXLAYERS,MAX_ATMOSWFS)
      
!  Homogeneous solutions

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

!  Particular integrals and thermal

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Linearized inputs

      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN  ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER  ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

!  Output
!  ======

!  Quad-angle output

      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_THTONLY_SOURCE ( MAXSTREAMS, MAX_ATMOSWFS )

!  User-angle output

      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_MSSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_DBSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          M, N, NN, J, I, UM, NELEMENTS
      INTEGER ::          Q, IB, O1, O2, OM, O11
      INTEGER ::          K, KO1, K0, K1, K2, LUM

      LOGICAL ::          DO_QTHTONLY
      DOUBLE PRECISION :: DOWN   (MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION :: L_DOWN (MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)

      DOUBLE PRECISION :: REFLEC, S_REFLEC, FAC, KMULT
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION :: LXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION :: LXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION :: LXR2, NXR2, LLXR2

!  Starting section
!  ----------------

!  Fourier number, layer number

      M   = FOURIER
      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

!  Special flag
!    Version 2.8. QUAD_OUTPUT flag removed
!      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND. ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

      DO_QTHTONLY = DO_THERMAL_TRANSONLY .AND. DO_INCLUDE_MVOUTPUT

!  initialise linearized BOA source functions
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            DO Q = 1, NV_PARAMETERS
               L_BOA_MSSOURCE(UM,1:NSTOKES,Q) = ZERO
               L_BOA_DBSOURCE(UM,1:NSTOKES,Q) = ZERO
            ENDDO
         ENDDO
      ENDIF

!  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
         DO I = 1, NSTREAMS
            L_BOA_THTONLY_SOURCE(I,1:NV_PARAMETERS) = ZERO
         ENDDO
      ENDIF

!  Number of Elements
!    Only want the (1,1) component for Lambertian

      NELEMENTS = 1
      IF ( .not. DO_LAMBERTIAN_SURFACE ) NELEMENTS = NSTOKES

!  Add direct beam if flagged
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_INCLUDE_DIRECTRF.AND..NOT.DO_DBCORR .AND.DO_USER_STREAMS) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO O1 = 1, NELEMENTS
               FAC = - RF_USER_DIRECT_BEAM(UM,IB,O1) * DELTAU_SLANT(N,NV,IB)
               DO Q = 1, NV_PARAMETERS
                  L_BOA_DBSOURCE(UM,O1,Q) = L_BOA_DBSOURCE(UM,O1,Q) + L_DELTAU_VERT(Q,NV) * FAC
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Add surface-term linearization if flagged. O1 = 1 (I-component only)
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_INCLUDE_DIRECTSL .AND. DO_USER_STREAMS ) THEN
         O1 = 1
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            FAC = SL_USERTERM(UM,IB,O1)
            DO Q = 1, NV_PARAMETERS
               L_BOA_DBSOURCE(UM,O1,Q) = L_BOA_DBSOURCE(UM,O1,Q) + LP_TRANS_ATMOS_FINAL(IB,NV,Q) * FAC
            ENDDO
         ENDDO
      ENDIF

!  Exit if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Include surface only, from now on

!  1. Thermal Transmittance only
!     %%%%%%%%%%%%%%%%%%%%%%%%%%

!  Thermal transmittance solution, build from TOA downwards

      IF ( DO_THERMAL_TRANSONLY ) THEN

!  Initialise

         O1 = 1
         DO I = 1, NSTREAMS
            DOWN(I,O1) = ZERO
            L_DOWN(I,O1,1:NV_PARAMETERS) = ZERO
         ENDDO

!  Build

         DO NN = 1, NLAYERS
            IF ( NV.EQ.NN ) THEN
               DO I = 1, NSTREAMS
                  DO Q = 1, NV_PARAMETERS
                     L_DOWN(I,O1,Q) = L_DOWN(I,O1,Q) *   T_DELT_DISORDS(I,NN)     &
                                    +   DOWN(I,O1)   * L_T_DELT_DISORDS(I,NN,Q) + L_T_WLOWER(I,NN,Q)
                  ENDDO
               ENDDO
            ELSE
               DO I = 1, NSTREAMS
                  DO Q = 1, NV_PARAMETERS
                     L_DOWN(I,O1,Q) = L_DOWN(I,O1,Q) * T_DELT_DISORDS(I,NN)
                  ENDDO
               ENDDO
            ENDIF
            DO I = 1, NSTREAMS
               DOWN(I,O1) = DOWN(I,O1)*T_DELT_DISORDS(I,NN) + T_WLOWER(I,NN)
            ENDDO
         ENDDO

!  Continuation point for avoiding the next section
!  Version 2.8 remove GOTO -->        GO TO 5432

!  2. Scattering solutions
!     %%%%%%%%%%%%%%%%%%%%

      ELSE

!  Two cases:
!  (a) If  NV = N, this is also the layer that is varying --> Extras!
!  (b) If  N > NV with variations in layer NV above N

         IF ( NV .EQ. N ) THEN

!  stream and stokes loops

            DO I = 1, NSTREAMS
              DO O1 = 1, NELEMENTS

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
                    K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
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

!  real part and add particular solution

                  SHOM = SHOM_R + SHOM_CR

!  Add Particular integral linearization
!    Does not exist if thermal only and N > K, K = 0

                  SPAR = ZERO
                  IF ( DO_SOLAR_SOURCES .OR. (DO_INCLUDE_THERMEMISS .AND. N.EQ.NV) ) THEN
                    SPAR = L_WLOWER(I,O1,N,Q)
                  ENDIF

!  Final result

                  L_DOWN(I,O1,Q) = SPAR + SHOM

!  End loops over Q, O1 and I

                ENDDO
              ENDDO
            ENDDO

!  otherwise the varying layer is above the boundary layer

         ELSE IF ( NV.LT.N ) THEN

!  stream and stokes loops

            DO I = 1, NSTREAMS
              DO O1 = 1, NELEMENTS

!  Parameter loop

                DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

                  SHOM_R = ZERO
                  DO K = 1, K_REAL(N)
                    NXR = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                    PXR = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                    HOM1 = NXR * T_DELT_EIGEN(K,N)
                    HOM2 = PXR
                    SHOM_R = SHOM_R + HOM1 + HOM2
                  ENDDO

!  Complex homogeneous solutions

                  SHOM_CR = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1

                    NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &  
                            - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)  
                    NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                            + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                    PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &  
                            - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                    HOM1CR =   NXR1 * T_DELT_EIGEN(K1,N) &
                             - NXR2 * T_DELT_EIGEN(K2,N)
                    HOM2CR = PXR1
                    SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                ENDDO

!  real part homogeneous solution

                  SHOM = SHOM_R + SHOM_CR

!  Add Particular integral linearization
!    Does not exist if thermal only and N > K, K = 0

                  SPAR = ZERO
                  IF ( DO_SOLAR_SOURCES .OR. (DO_INCLUDE_THERMEMISS .AND. N.EQ.NV) ) THEN
                    SPAR = L_WLOWER(I,O1,N,Q)  
                  ENDIF

!  Final result

                  L_DOWN(I,O1,Q) = SPAR + SHOM

!  End loops over Q, O1 and I

                ENDDO
              ENDDO
            ENDDO

!  End cases

         ENDIF

!  Continuation point. Removed Version 2.8
! 5432   CONTINUE

!  End scattering vs. thermal-only

      ENDIF

!  1/31/21. Version 2.8.3. Reflectance integrand consolidated here, as is done with LIDORT
!    -- Was a bug with the THERMAL_TRANSONLY implementation

!    Set reflectance integrand  a(j).x(j).L_DOWN(-j)  Scattering solutions
!    Set reflectance integrand       a(j).L_DOWN(-j)  Thermal tranmsittance

      IF ( DO_THERMAL_TRANSONLY ) THEN
         DO O1 = 1, NELEMENTS
            DO Q = 1, NV_PARAMETERS
               DO I = 1, NSTREAMS
                  L_DOWN(I,O1,Q) = L_DOWN(I,O1,Q) * QUAD_WEIGHTS(I)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO O1 = 1, NELEMENTS
            DO Q = 1, NV_PARAMETERS
               DO I = 1, NSTREAMS
                  L_DOWN(I,O1,Q) = L_DOWN(I,O1,Q) * QUAD_STRMWTS(I)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  ###### Lambertian reflectance (same for all user-streams)
!    -- 1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
         IF ( FOURIER .EQ. 0 ) THEN
            O1 = 1 ; KMULT = SURFACE_FACTOR * ALBEDO
            DO Q = 1, NV_PARAMETERS
               REFLEC = KMULT * SUM(L_DOWN(1:NSTREAMS,O1,Q))
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  L_BOA_MSSOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) + REFLEC
               ENDDO
               IF ( DO_QTHTONLY ) THEN
                  DO I = 1, NSTREAMS
                     L_BOA_THTONLY_SOURCE(I,Q) = L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDIF

!  .. integrate with BRDF reflectance function at user angles
!  @@@@@@@@ Rob Fix, 2/9/11, DO_QTHTONLY clause was absent
!    -- 1/31/21. Version 2.8.3. Use post-processing mask
!    -- 1/31/21. Version 2.8.3. BRDF arrays defined locally, drop Fourier M index

      IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
        KMULT = SURFACE_FACTOR
        DO Q = 1, NV_PARAMETERS
          IF ( DO_USER_STREAMS ) THEN
            DO LUM = 1, N_PPSTREAMS
              UM = PPSTREAM_MASK(LUM,IB)
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + L_DOWN(J,O2,Q) * USER_BRDF_F(OM,UM,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                REFLEC = KMULT * REFLEC
                L_BOA_MSSOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) + REFLEC
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            O11 = 1
            DO I = 1, NSTREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + L_DOWN(J,O11,Q) * BRDF_F(O11,I,J)
              ENDDO
              REFLEC = KMULT * REFLEC
              L_BOA_THTONLY_SOURCE(I,Q) = L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE GET_LP_BOASOURCE

!

      SUBROUTINE MIFLUX_PROFILEWF ( &
        DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, DO_CLASSICAL,         & ! Input flags
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! Input flags
        IB, NV, NV_PARAMETERS, NSTOKES, NSTREAMS, NSTKS_NSTRMS, NLAYERS,  & ! Input numbers
        N_USER_LEVELS, N_DIRS, WHICH_DIRS, LEVEL_MASK_UP, LEVEL_MASK_DN,  & ! Level/Dir output control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,     & ! Partial layer Output control
        FLUX_MULTIPLIER, FLUX_FACTOR, FLUXVEC, LOCAL_CSZA, TAYLOR_ORDER,        & ! Input Flux/Angles/Taylor
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, PARTAU_VERT, L_DELTAU_VERT,   & ! Input quadrature/Optical
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR, & ! Beam parameterization
        T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP, K_REAL, K_COMPLEX,      & ! DODTrans/Bookkeep
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,         & ! Homog solution
        ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GMULT_UP, UT_GMULT_DN,     & ! Green's function solution
        WUPPER, LCON, MCON, T_WLOWER, T_WUPPER, BOA_THTONLY_SOURCE,             & ! SOlar/BVP/Thermal
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,                   & ! Lin. Beam parameterization
        LP_T_UTDN_MUBAR, LP_LEVELS_SOLARTRANS, LP_PARTIALS_SOLARTRANS,          & ! Lin. Beam transmittances
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, L_KEIGEN,           & ! Lin. DODTrans/Keigen
        L_SOLA_XPOS, L_SOLB_XNEG, L_T_DELT_EIGEN, L_T_UTDN_EIGEN, L_T_UTUP_EIGEN, & ! Lin. Homog solution
        L_WLOWER, L_WUPPER, NCON, PCON, L_ATERM_SAVE, L_BTERM_SAVE,               & ! Lin. PI/Greens Solutons
        L_UT_T_PARTIC, L_T_WLOWER, L_T_WUPPER, L_BOA_THTONLY_SOURCE,              & ! Lin. Thermal solutions
        MEANST_DIFFUSE_PROFWF, DNMEANST_DIRECT_PROFWF,                            & ! OUTPUT Actinic Fluxes linearized
        FLUX_DIFFUSE_PROFWF, DNFLUX_DIRECT_PROFWF )                                 ! OUTPUT Regular Fluxes linearized

!  Quadrature output at offgrid or ongrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )


! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
!mick mod 9/19/2017 - output flux variables renamed: distinguish between diffuse and direct
!mick fix 9/19/2017 - added LP_LEVELS_SOLARTRANS & LP_PARTIALS_SOLARTRANS to facilitate correction
!                     of linearized direct flux   

!  1/31/21. Version 2.8.3. Rewritten to incorporate new Linearized Green's function solution
!     -- New arguments include DO_CLASSICAL, TAYLOR_ORDER, PARTAU_VERT, L_DELTAU_VERT
!     -- New Greens function variables ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GMULT_UP, UT_GMULT_DN
!     -- Logic changed to perform adadition Partial-layer Green function multipliers initially

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAX_USER_LEVELS, MAX_SZANGLES, MAX_DIRECTIONS, MAX_ATMOSWFS,            &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, HALF, PI2, PI4, UPIDX, DNIDX

      IMPLICIT NONE

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument 

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM

!  1/31/21. Version 2.8.3. New flag for solution options

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL

!  Source flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Input basic numbers

      INTEGER, INTENT (IN) ::          IB, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS, NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  input directions

      INTEGER, INTENT (IN) ::          N_DIRS
      INTEGER, INTENT (IN) ::          WHICH_DIRS ( MAX_DIRECTIONS )

!  Level output control

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  Partial output control

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  1/31/21. Version 2.8.3. Several new inputs
!    -- Order of Taylor series (including terms up to EPS^n)
!    -- Local vertical optical depth and its linearization

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER
      DOUBLE PRECISION, INTENT(IN)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      DOUBLE PRECISION, intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Beam parameterization
!    -- 1/31/21. Version 2.8.3. Some new inputs

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN  ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )

!  Discrete Ordinate solutions

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER     ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Green's function particular integral arrays

      DOUBLE PRECISION, intent(in)  :: ATERM_SAVE(MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: BTERM_SAVE(MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: GAMMA_M(MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: GAMMA_P(MAXEVALUES,MAXLAYERS)

!  1/31/21. Version 2.8.3. Green's function multipliers for off-grid optical depths

      DOUBLE PRECISION, intent (in) :: UT_GMULT_UP(MAXEVALUES,MAX_PARTLAYERS)
      DOUBLE PRECISION, intent (in) :: UT_GMULT_DN(MAXEVALUES,MAX_PARTLAYERS)

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  Linearized Input
!  ================

!  Beam parameterization

      DOUBLE PRECISION, intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      DOUBLE PRECISION, intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: LP_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Derived Solar-beam Linearized Transmittance at all levels

      DOUBLE PRECISION, INTENT (IN) :: LP_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS  ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN  ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP  ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Discrete Ordinate solutions
!    -- 1/31/21. Version 2.8.3. L_KEIGEN added

      DOUBLE PRECISION, intent(in)  :: L_KEIGEN    ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )


!  1/31/21. Version 2.8.3. Linearizations of Saved quantities for the Green function solution
!mick fix 1/5/2021 - changed 1st dim of L_BTERM_SAVE from MAXSTREAMS to MAXEVALUES

      DOUBLE PRECISION, intent(in)  :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent(in)  :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE  ( MAXSTREAMS, MAX_ATMOSWFS )

!  Linearized output (Has to be INOUT)

      DOUBLE PRECISION, INTENT (INOUT) :: MEANST_DIFFUSE_PROFWF &
             ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: DNMEANST_DIRECT_PROFWF &
             ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_DIFFUSE_PROFWF &
             ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: DNFLUX_DIRECT_PROFWF &
             ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  local variables
!  ---------------

!  Local quadrature output (for debug)

      DOUBLE PRECISION :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  Linearized Green functions multipliers for off-grid optical depths
!    -- 1/31/21. Version 2.8.3. Added for the Green's function calculation

      DOUBLE PRECISION :: LP_UT_GMULT_UP(MAXEVALUES,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: LP_UT_GMULT_DN(MAXEVALUES,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, Q, N, O1, NLEVEL
      DOUBLE PRECISION :: SMI, SFX, L_TRANS, L_FTRANS, L_DNDIRECT_MEANST, L_DNDIRECT_FLUX

!  Mick fix 9/19/2017 - added to facilitate correction of direct flux

      DOUBLE PRECISION :: HELP, L_TRANS_SCALED, L_FTRANS_SCALED, L_DNDIRECT_MEANST_SCALED, L_DNDIRECT_FLUX_SCALED

!  First get the partial-layer multipliers (Green's function)
!  ==========================================================

!  1/31/21. Version 2.8.3. This section is completely new.

      IF ( DO_INCLUDE_MVOUTPUT .and. .not. DO_CLASSICAL ) THEN
        DO UTA = 1, N_USER_LEVELS
          NLEVEL = LEVEL_MASK_UP(UTA)
          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            CALL LP_QUAD_GFUNCMULT &
              ( IB, UT, N, NV, NV_PARAMETERS, TAYLOR_ORDER, PARTAU_VERT, L_DELTAU_VERT,  & ! Input
                BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,  & ! Input
                T_UTUP_EIGEN, T_UTDN_EIGEN, ATERM_SAVE, BTERM_SAVE,                      & ! Input
                K_REAL, GAMMA_M, GAMMA_P, UT_GMULT_UP, UT_GMULT_DN,                      & ! Input
                LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,   & ! Input
                L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_KEIGEN, L_ATERM_SAVE, L_BTERM_SAVE,    & ! Input
                LP_UT_GMULT_UP, LP_UT_GMULT_DN )                                           ! Output
          ENDIF
        ENDDO
      ENDIF

!  Get the Quadrature fields
!  =========================

!  direction loop

      DO IDIR = 1, N_DIRS
        WDIR = WHICH_DIRS(IDIR)

!  Upwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

!  1/31/21. Version 2.8.3. Rewritten to incorporate new Linearized Green's function solution
!     -- New code modeled after the LIDORT Version 3.8.1.
!     -- New arguments include DO_CLASSICAL
!     -- New Greens' multipliers UT_GMULT_UP, UT_GMULT_DN, LP_UT_GMULT_UP, LP_UT_GMULT_DN

              CALL QUADPROFILEWF_OFFGRID_UP ( &
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, DO_CLASSICAL, & ! Input flags
                N, UTA, UT, IB, NV, NV_PARAMETERS, NSTOKES, NSTREAMS, NLAYERS, NSTKS_NSTRMS, & ! Input numbers
                FLUX_MULTIPLIER, QUAD_STREAMS, T_UTDN_MUBAR, LP_T_UTDN_MUBAR,                & ! Input Flux/Beam/Quad
                T_DELT_DISORDS, T_DISORDS_UTUP, L_T_DELT_DISORDS, L_T_DISORDS_UTUP,          & ! Input disords
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN,         & ! Input Homogeneous 
                UT_GMULT_UP, UT_GMULT_DN, LCON, MCON, WUPPER,                                & ! Input SolarPI
                T_WUPPER, BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,                          & ! Input Thermal
                L_SOLA_XPOS, L_SOLB_XNEG, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, NCON, PCON,        & ! Input Lin. Homogeneous
                LP_UT_GMULT_UP, LP_UT_GMULT_DN, L_WUPPER, L_T_WUPPER, L_UT_T_PARTIC,         & ! Input Lin. SolarPI/Thermal
                QATMOSWF_F )                                                                   ! Output

            ELSE

!  1/31/21. Version 2.8.3. No Changes here

              CALL QUADPROFILEWF_LEVEL_UP ( &
                DO_THERMAL_TRANSONLY, NLEVEL, UTA, NV, NV_PARAMETERS, NSTOKES,     & ! Input flag/numbers
                NSTREAMS, NLAYERS, FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,  & ! Input numbers/flux/quad/D.O.trans
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,             & ! Input Homog. solutions
                LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,                          & ! Input thermal and PI solutions
                L_T_DELT_DISORDS, L_SOLA_XPOS, L_SOLB_XNEG,  L_T_DELT_EIGEN,       & ! Linearized Homog. solutions
                L_WLOWER, L_WUPPER, NCON, PCON, L_T_WUPPER, L_BOA_THTONLY_SOURCE,  & ! Linearized thermal and PI solutions
                QATMOSWF_F )                                                         ! OUTPUT

            ENDIF
          ENDDO
        ENDIF

!  Downwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

!  1/31/21. Version 2.8.3. Rewritten to incorporate new Linearized Green's function solution
!     -- New code modeled after the LIDORT Version 3.8.1.
!     -- New arguments include DO_CLASSICAL
!     -- New Greens' multipliers UT_GMULT_UP, UT_GMULT_DN, LP_UT_GMULT_UP, LP_UT_GMULT_DN

              CALL QUADPROFILEWF_OFFGRID_DN ( &
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, DO_CLASSICAL,     & ! Input flags
                N, UTA, UT, IB, NV, NV_PARAMETERS, NSTOKES, NSTREAMS, NSTKS_NSTRMS,              & ! Input numbers
                FLUX_MULTIPLIER, QUAD_STREAMS, T_UTDN_MUBAR, LP_T_UTDN_MUBAR,                    & ! Input solar beam
                T_DELT_DISORDS, T_DISORDS_UTDN, L_T_DELT_DISORDS, L_T_DISORDS_UTDN,              & ! Input disords
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN,             & ! Input Homog. solution
                UT_GMULT_UP, UT_GMULT_DN, LCON, MCON, WUPPER, T_WLOWER, T_WUPPER,                & ! Input SOlat/Thermal PIs
                L_SOLA_XPOS, L_SOLB_XNEG, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, NCON, PCON,            & ! Input Lin. Homog.
                LP_UT_GMULT_UP, LP_UT_GMULT_DN, L_WUPPER, L_T_WLOWER, L_T_WUPPER, L_UT_T_PARTIC, & ! Input Lin PIs
                QATMOSWF_F )                                                                        ! Output

            ELSE

!  1/31/21. Version 2.8.3. No Changes here

              CALL QUADPROFILEWF_LEVEL_DN ( &
                DO_THERMAL_TRANSONLY, NLEVEL, UTA, NV, NV_PARAMETERS, NSTOKES,     & ! Input flag/numbers
                NSTREAMS, FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,           & ! Input numbers/flux/quad/D.O.trans
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, & ! Input Homog. solutions
                T_WLOWER, L_T_DELT_DISORDS, L_SOLA_XPOS, L_SOLB_XNEG,              & ! Linearized Homog. solutions
                L_T_DELT_EIGEN, L_WLOWER, NCON, PCON, L_T_WLOWER,                  & ! Linearized thermal and PI solutions
                QATMOSWF_F )                                                         ! OUTPUT

            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux  output
!  -------------------------------

        IF ( DO_INCLUDE_MVOUTPUT ) THEN

!  Diffuse term integrated output

          DO Q = 1, NV_PARAMETERS
            DO UTA = 1, N_USER_LEVELS
              DO O1 = 1, NSTOKES
                SMI = ZERO
                SFX = ZERO
                DO I = 1, NSTREAMS
                  SMI = SMI + QUAD_WEIGHTS(I) * QATMOSWF_F(Q,UTA,I,O1)
                  SFX = SFX + QUAD_STRMWTS(I) * QATMOSWF_F(Q,UTA,I,O1)
                ENDDO
                MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR) = SMI * HALF
                FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR)   = SFX * PI2
              ENDDO
            ENDDO
          ENDDO

!  nothing further to do if no solar sources
!   Version 2.8, remove GOTO statment, 7/19/16
!          IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) GO TO 455

          IF ( DO_INCLUDE_DIRECTBEAM ) THEN

!  For the downward direction, add the direct beam contributions
!mick fix 9/19/2017 - use unscaled optical thicknesses for calculation of linearized
!                     direct flux (following LIDORT)

            IF ( WDIR .EQ. DNIDX ) THEN

!  loop over all the output optical depths

              DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

                IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
                  UT = PARTLAYERS_OUTINDEX(UTA)
                  N  = PARTLAYERS_LAYERIDX(UT)

!  For the offgrid values.......
!     .....Only contributions for layers above the PI cutoff
!    LP_INITIAL_TRANS is a logarithmic derivative

!                  IF ( N .LE. BEAM_CUTOFF(IB) ) THEN
!                    IF ( NV.LE.N ) THEN
!                      DO Q = 1, NV_PARAMETERS
!                        L_TRANS = LP_T_UTDN_MUBAR(UT,NV,IB,Q) + &
!                                  LP_INITIAL_TRANS(N,NV,IB,Q) * T_UTDN_MUBAR(UT,IB)
!                        L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
!                        DO O1 = 1, NSTOKES
!
!!  @@@ Rob fix (Forgot the Flux factor)
!!                  FTRANS = FLUXVEC(O1) * L_TRANS
!!                  FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!!                   FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * L_TRANS
!! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
!
!                          FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!                          L_DIRECT_MEANI = FTRANS / PI4
!                          L_DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IB)
!                          MINT_PROFILEWF_DIRECT(Q,NV,UTA,IB,O1)=L_DIRECT_MEANI
!                          FLUX_PROFILEWF_DIRECT(Q,NV,UTA,IB,O1)=L_DIRECT_FLUX
!                          MINT_PROFILEWF(Q,NV,UTA,IB,O1,WDIR) = &
!                                  MINT_PROFILEWF(Q,NV,UTA,IB,O1,WDIR) + L_DIRECT_MEANI
!                          FLUX_PROFILEWF(Q,NV,UTA,IB,O1,WDIR) = &
!                                  FLUX_PROFILEWF(Q,NV,UTA,IB,O1,WDIR) + L_DIRECT_FLUX
!                        ENDDO
!                      ENDDO
!                    ENDIF
!                  ENDIF

                  IF ( N .LE. BEAM_CUTOFF(IB) ) THEN
                    IF ( NV.LE.N ) THEN
                      DO Q = 1, NV_PARAMETERS

!  Linearized direct transmittances, scaled and unscaled

                        L_TRANS = LP_PARTIALS_SOLARTRANS(UT,NV,IB,Q)
                        HELP    = LP_T_UTDN_MUBAR(UT,NV,IB,Q) + &
                                  LP_INITIAL_TRANS(N,NV,IB,Q) * T_UTDN_MUBAR(UT,IB)
                        L_TRANS_SCALED = HELP * INITIAL_TRANS(N,IB)

                        DO O1 = 1, NSTOKES

!  Direct calculation with non-scaled linearized transmittances

                          L_FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
                          L_DNDIRECT_MEANST = L_FTRANS / PI4
                          L_DNDIRECT_FLUX   = L_FTRANS * LOCAL_CSZA(N,IB)
                          DNMEANST_DIRECT_PROFWF(Q,NV,UTA,IB,O1) = L_DNDIRECT_MEANST
                          DNFLUX_DIRECT_PROFWF(Q,NV,UTA,IB,O1)   = L_DNDIRECT_FLUX

!  Diffuse calculation

                          L_FTRANS_SCALED = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS_SCALED
                          L_DNDIRECT_MEANST_SCALED = L_FTRANS_SCALED / PI4
                          L_DNDIRECT_FLUX_SCALED   = L_FTRANS_SCALED * LOCAL_CSZA(N,IB)
                          MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR) = &
                                MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR) + ( L_DNDIRECT_MEANST_SCALED - L_DNDIRECT_MEANST )
                          FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR) = &
                                FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR)   + ( L_DNDIRECT_FLUX_SCALED - L_DNDIRECT_FLUX )
                        ENDDO
                      ENDDO
                    ENDIF
                  ENDIF

!  For the on-grid values
!    LP_INITIAL_TRANS is a logarithmic derivative

                ELSE
                  N = LEVEL_MASK_DN(UTA)

!                  IF ( N .LE. BEAM_CUTOFF(IB) ) THEN
!                    IF ( N.GT.0 ) THEN
!                      IF ( NV.LE.N ) THEN
!                        DO Q = 1, NV_PARAMETERS
!                          L_TRANS = LP_T_DELT_MUBAR(N,NV,IB,Q) + &
!                                    LP_INITIAL_TRANS(N,NV,IB,Q) * T_DELT_MUBAR(N,IB)
!                          L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
!                          DO O1 = 1, NSTOKES
!
!!  @@@ Rob fix (Forgot the Flux factor)
!!                  FTRANS = FLUXVEC(O1) * L_TRANS
!!                  FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!!                   FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * L_TRANS
!! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
!
!                            FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!                            L_DIRECT_MEANI = FTRANS / PI4
!                            L_DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IB)
!                            MINT_PROFILEWF_DIRECT(Q,NV,UTA,IB,O1)=L_DIRECT_MEANI
!                            FLUX_PROFILEWF_DIRECT(Q,NV,UTA,IB,O1)=L_DIRECT_FLUX
!                            MINT_PROFILEWF(Q,NV,UTA,IB,O1,WDIR) = &
!                                     MINT_PROFILEWF(Q,NV,UTA,IB,O1,WDIR) + L_DIRECT_MEANI
!                            FLUX_PROFILEWF(Q,NV,UTA,IB,O1,WDIR) = &
!                                     FLUX_PROFILEWF(Q,NV,UTA,IB,O1,WDIR) + L_DIRECT_FLUX
!                          ENDDO
!                        ENDDO
!                      ENDIF
!                    ENDIF
!                  ENDIF

                  IF ( N .LE. BEAM_CUTOFF(IB) ) THEN
                    IF ( N.GT.0 ) THEN
                      IF ( NV.LE.N ) THEN
                        DO Q = 1, NV_PARAMETERS

!  Linearized direct transmittances, scaled and unscaled

                          L_TRANS = LP_LEVELS_SOLARTRANS(N,NV,IB,Q)
                          HELP    = LP_T_DELT_MUBAR(N,NV,IB,Q) + &
                                    LP_INITIAL_TRANS(N,NV,IB,Q) * T_DELT_MUBAR(N,IB)
                          L_TRANS_SCALED = HELP * INITIAL_TRANS(N,IB)

                          DO O1 = 1, NSTOKES

!  Direct calculation with non-scaled linearized transmittances

                            L_FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
                            L_DNDIRECT_MEANST = L_FTRANS / PI4
                            L_DNDIRECT_FLUX   = L_FTRANS * LOCAL_CSZA(N,IB)
                            DNMEANST_DIRECT_PROFWF(Q,NV,UTA,IB,O1) = L_DNDIRECT_MEANST
                            DNFLUX_DIRECT_PROFWF(Q,NV,UTA,IB,O1)   = L_DNDIRECT_FLUX

!  Diffuse calculation

                            L_FTRANS_SCALED = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS_SCALED
                            L_DNDIRECT_MEANST_SCALED = L_FTRANS_SCALED / PI4
                            L_DNDIRECT_FLUX_SCALED   = L_FTRANS_SCALED * LOCAL_CSZA(N,IB)
                            MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR) = &
                                   MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR) + ( L_DNDIRECT_MEANST_SCALED - L_DNDIRECT_MEANST )
                            FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR) = &
                                   FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IB,O1,WDIR)   + ( L_DNDIRECT_FLUX_SCALED - L_DNDIRECT_FLUX )
                          ENDDO
                        ENDDO
                      ENDIF
                    ENDIF
                  ENDIF

                ENDIF

!  End UTA loop

              ENDDO

!  Finish downwelling direct contribution

            ENDIF

!  Continuation point for avoiding direct beam calculation
! 455      CONTINUE. Version 2.8, GOTO REMOVED 7/19/16. Replaced by IF clause

          ENDIF

!  Finish MV output

        ENDIF

!  end direction loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MIFLUX_PROFILEWF

!  End Module

      END MODULE vlidort_lp_wfatmos_m
