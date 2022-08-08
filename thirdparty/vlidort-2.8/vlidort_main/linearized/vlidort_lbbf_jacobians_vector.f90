
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
! #            vlidort_lbbf_jacobians_whole                     #
! #            vlidort_lbbf_jacobians_wpartials                 #
! #                                                             #
! ###############################################################

module vlidort_lbbf_jacobians_m

!  Linearization w.r.t  BB input variables.

!  LIDORT HISTORY :--
!    First   Attempt, 27 January 2014.
!    Second  Attempt, 19 March   2014. Success!  No partials             THIS ROUTINE
!    Third   Attempt, 21 March   2014. Partials Introduced.              NOT THIS ROUTINE
!    Fourth  Attempt, 25 March   2014. Thermal Transmittance only.       BOTH ROUTINES.

!  VLIDORT HISTORY :--
!    First   Attempt, 27 March    2014. CODE DISABLED for VERSION 2.7.
!    Second  Attempt, 20-21 April 2016. Success (Both routines) !!!
!       - Downwelling and Upwelling ABBFs, SBBF; Whole and Partials.
!       - Downwelling and Upwelling Flux ABBFs, SBBF; Whole and Partials.
!       -  Complex variable code not tested
!       - TRANSONLY code not tested
!       - Presented for Version 2.8, 11 July 2016

!  1/31/21. Version 2.8.3. BRDF arrays are defined locally for each Fourier

public

contains

subroutine vlidort_lbbf_jacobians_whole &
      ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
        DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
        DO_MSMODE_THERMAL, DO_POSTPROCESSING, DO_MVOUTPUT,         & ! input
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,                 & ! input
        NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
        NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & ! Input
        NTOTAL, N_SUPDIAG, N_SUBDIAG, MUELLER_INDEX,               & ! input
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                    & ! Input
        USER_STREAMS, LAYERMASK_UP, LAYERMASK_DN,                  & ! Input
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
        SURFACE_FACTOR, ALBEDO, BRDF_F, UBRDF_F,                   & ! input
        EMISSIVITY, USER_EMISSIVITY,                               & ! input
        FLUX_MULTIPLIER, DELTAUS, OMEGA_GREEK,                     & ! Input
        T_DELT_DISORDS, T_DELT_USERM,                              & ! input
        PI_XQP, PI_XQM_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
        K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,     & ! input
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
        UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
        ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
        SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
        STATUS, MESSAGE, TRACE )                                     ! Output

!  Module file of dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, MAX_USER_STREAMS,     &
                                 MAX_USER_LEVELS, MAXMOMENTS, MAXTOTAL, MAXBANDTOTAL, MAXEVALUES,        &
                                 MAXSTRMSTKS, MAXSTREAMS_2, MAXSTRMSTKS_2, MAXSTOKES_SQ, MAX_DIRECTIONS, &
                                 UPIDX, DNIDX, ZERO, HALF, ONE, TWO, PI2, PI4, VLIDORT_SERIOUS, VLIDORT_SUCCESS

!  VLIDORT module dependencies

   USE lapack_tools_m, only : DGETRF, DGETRS, DGBTRS

      implicit none

!  Subroutine input arguments
!  --------------------------

!  Master control

      LOGICAL, INTENT(IN)  :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  local control flags

      LOGICAL, INTENT(IN)  :: DO_THERMAL_TRANSONLY
      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

      LOGICAL, INTENT(IN)  :: DO_MSMODE_THERMAL
      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING
      LOGICAL, INTENT(IN)  :: DO_MVOUTPUT

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_LAMBERTIAN_SURFACE

!  Numbers

      INTEGER, INTENT(IN)  :: NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, NMOMENTS
      INTEGER, INTENT(IN)  :: NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NSTREAMS_2
      INTEGER, INTENT(IN)  :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  multiplier

      DOUBLE PRECISION, INTENT(IN)  :: FLUX_MULTIPLIER

!  output   control

      INTEGER, INTENT (IN)   :: N_USER_LEVELS
      INTEGER, INTENT (IN)   :: UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN)   :: UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  User polar directions, postprocessnd control

      DOUBLE PRECISION, INTENT(IN)  :: USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL  , INTENT(IN)  :: LAYERMASK_UP ( MAXLAYERS )
      LOGICAL  , INTENT(IN)  :: LAYERMASK_DN ( MAXLAYERS )

!  Quadrature values

      DOUBLE PRECISION, intent(in)   :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, intent(in)   :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, intent(in)   :: QUAD_STRMWTS ( MAXSTREAMS )

!  Optical properties

      DOUBLE PRECISION, INTENT(IN)   :: DELTAUS     ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN)   :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Polynomials

      DOUBLE PRECISION, INTENT(IN)   :: PI_XQP     ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT(IN)   :: PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT(IN)   :: PI_XUP     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT(IN)   :: PI_XUM     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Discrete ordinate solutions
!  ---------------------------

!  Direct solutions, stream transmittances

      DOUBLE PRECISION, intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Eigensolutions, eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER         , INTENT (IN) :: K_REAL    ( MAXLAYERS )
      INTEGER         , INTENT (IN) :: K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Matrices

      DOUBLE PRECISION, INTENT (IN) :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )

!  BVProblem stuff
!  ---------------

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER         , INTENT (IN) :: IPIVOT   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2    ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER         , INTENT (IN) :: SIPIVOT  ( MAXSTRMSTKS_2 )

!  Surface stuff
!  -------------

!  1/31/21. Version 2.8.3. BRDF arrays are defined locally for each Fourier, drop "MAXMOMENTS" dimension

      DOUBLE PRECISION, INTENT(IN)   :: SURFACE_FACTOR, ALBEDO
      DOUBLE PRECISION, intent(in)   :: BRDF_F  ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN)   :: UBRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION, intent(in)   :: EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, intent(in)   :: USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  User-angle (post-processed) solution variables
!  ----------------------------------------------

!  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION, INTENT(IN)  :: T_DELT_USERM  (MAXLAYERS,MAX_USER_STREAMS)

!  User solutions defined at user-defined stream angles

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  solution multipliers 

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Outputs
!  -------

!  Postprocessed Jacobians.
!  Outputs are all Pre-zeroed in the calling Masters

      DOUBLE PRECISION, INTENT(INOUT) :: ABBWFS_JACOBIANS &
               ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)
      DOUBLE PRECISION, INTENT(INOUT) :: SBBWFS_JACOBIANS &
               ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXSTOKES, MAX_DIRECTIONS)

!  Flux Jacobians.
!  Outputs are all Pre-zeroed in the calling Masters

      DOUBLE PRECISION, INTENT(INOUT) :: ABBWFS_FLUXES ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(INOUT) :: SBBWFS_FLUXES ( MAX_USER_LEVELS, 2, MAXSTOKES, MAX_DIRECTIONS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  LOCAL THERMAL-BBF JACOBIAN ARRAYS
!  =================================

!  Weighting function column matrices

      DOUBLE PRECISION  :: COL2_BWF  ( MAXTOTAL )
      DOUBLE PRECISION  :: SCOL2_BWF ( MAXSTRMSTKS_2 )

!  Linearized Solution constants of integration

      DOUBLE PRECISION  :: NCON(MAXSTRMSTKS,MAXLAYERS)
      DOUBLE PRECISION  :: PCON(MAXSTRMSTKS,MAXLAYERS)

!  Linearized Thermal solutions at the Upper/Lower boundary

      DOUBLE PRECISION  :: L_T_WUPPER_Gp1(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION  :: L_T_WUPPER_Gp2(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION  :: L_T_WLOWER_Gp1(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION  :: L_T_WLOWER_Gp2(MAXSTREAMS_2,MAXLAYERS)

!  Linearized Thermal layer source terms

      DOUBLE PRECISION  :: L_LAYER_TSUP_UP_Gp1(MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_UP_Gp2(MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_DN_Gp1(MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_DN_Gp2(MAX_USER_STREAMS,MAXLAYERS)

!  Linearization of Direct solutions

      DOUBLE PRECISION  :: L_T_DIRECT_UP_Gp1 (MAX_USER_STREAMS,MAXLAYERS )
      DOUBLE PRECISION  :: L_T_DIRECT_UP_Gp2 (MAX_USER_STREAMS,MAXLAYERS )
      DOUBLE PRECISION  :: L_T_DIRECT_DN_Gp1 (MAX_USER_STREAMS,MAXLAYERS )
      DOUBLE PRECISION  :: L_T_DIRECT_DN_Gp2 (MAX_USER_STREAMS,MAXLAYERS )

!  Help arrays for the discrete ordinate field

      INTEGER          :: TPIVOT ( MAXSTREAMS )
      DOUBLE PRECISION :: TMAT ( MAXSTREAMS, MAXSTREAMS ), HVEC ( MAXSTREAMS )
      DOUBLE PRECISION :: TVEC1_Gp1 ( MAXSTREAMS_2, MAXLAYERS ), TVEC1_Gp2 ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: TVEC2_Gp1 ( MAXSTREAMS_2, MAXLAYERS ), TVEC2_Gp2 ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: HVEC1_Gp1 ( MAXSTREAMS ), HVEC1_Gp2 ( MAXSTREAMS )
      DOUBLE PRECISION :: HVEC2_Gp1 ( MAXSTREAMS ), HVEC2_Gp2 ( MAXSTREAMS )
      DOUBLE PRECISION :: JVEC1_Gp1 ( MAXSTREAMS ), JVEC1_Gp2 ( MAXSTREAMS )

!  Help arrays for the post-processed field

      DOUBLE PRECISION :: T_HELP1_Gp1 ( 0:MAXMOMENTS ), T_HELP1_Gp2 ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: T_HELP2_Gp1 ( 0:MAXMOMENTS ), T_HELP2_Gp2 ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: U_TPOS1_Gp1 ( MAX_USER_STREAMS, MAXLAYERS ), U_TPOS1_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TPOS2_Gp1 ( MAX_USER_STREAMS, MAXLAYERS ), U_TPOS2_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG1_Gp1 ( MAX_USER_STREAMS, MAXLAYERS ), U_TNEG1_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG2_Gp1 ( MAX_USER_STREAMS, MAXLAYERS ), U_TNEG2_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )

!  local variables
!  ---------------

!  Reflectance integrands, BOA source terms
!    Vector change 4/27/16. Give Reflec a Stokes dimension

      DOUBLE PRECISION :: L_IDOWN(MAXSTREAMS,MAXSTOKES), DOWN(MAXSTREAMS)
      DOUBLE PRECISION :: REFLEC(MAXSTREAMS,MAXSTOKES), BBWF_QUAD(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION :: L_BOA_MSSOURCE ( MAX_USER_STREAMS, MAXSTOKES ), L_BOA_THTONLY_SOURCE (MAXSTREAMS)

!  Local layer and cumulative source terms

      DOUBLE PRECISION :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: L_CUMULSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

!  help variables

      LOGICAL   :: DO_QTHTONLY
      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC, NL, K, KO1, K0, K1, K2, LAY, O1, O2, O11, OM, JB, AA, AA1, L
      INTEGER   :: UTA, UM, N1, M, I, I1, J, J1, INFO, IR, IROW, IROW1, IROW_S, IROW1_S, IC, IC1, IC2, CMP, CM, C0

      DOUBLE PRECISION :: H1, H2, HOM1, HOM2, HOM1CR, HOM2CR
      DOUBLE PRECISION :: NXR1, NXR2, PXR1, PXR2, NUXR1, NUXR2, PUXR1, PUXR2
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SHOM, POS1, POS2, NEG1, NEG2

      DOUBLE PRECISION :: SM, TM, SF, COSMUM, SUM, FAC, FACTOR, omega1_odelt, OMEGA
      DOUBLE PRECISION :: U1, U2, D1, D2, SPAR, SPAR1, SPAR2, EMISS, L_THELP, HOMIGR
      DOUBLE PRECISION :: Udel, DelUdel, omega1, A5_xqp, A5_xqm, z1, z, zd, z1_ok, S_REFLEC
      DOUBLE PRECISION :: PN1_Gp1, PN1_Gp2, PN2_Gp1, PN2_Gp2, Sum1_Gp1, Sum1_Gp2, Sum2_Gp1, Sum2_Gp2
      CHARACTER*3 :: CI, CN, C3

      DOUBLE PRECISION :: Group1(maxlayers,2), Group2(maxlayers,2)

!  Initial section
!  ---------------

!  Exception handling

   STATUS  = VLIDORT_SUCCESS
   MESSAGE = ' '
   TRACE   = ' '

!  Proxies

   m = 0

!  Use this variable when only the nstokes = 1 component is required

   o11 = 1

!  Initial modulus = 4.pi if solar sources are included

   fac = one
   if ( do_solar_sources ) fac = pi4

!  Local flag

   DO_QTHTONLY = do_MVOUTPUT .and. DO_THERMAL_TRANSONLY

!  Control to SURFACE LBBF Section

   if ( .not. DO_ATMOS_LBBF ) go to 55

!  Group 1/2 Derivatives of TCOM1, all layers
!   Group 1, w.r.t to the upper boundary BBF
!   Group 2, w.r.t to the lower boundary BBF
!     Assumes only 2 coefficients, simplified form

   do n = 1, nlayers
      omega  = omega_greek(0,n,1,1)
      omega1 = two * ( one - omega ) ; if ( do_thermal_transonly ) omega1 = two
      omega1_odelt = omega1 / deltaus(n)
      Group1(n,1)  =   omega1
      Group1(n,2)  = - omega1_odelt
      Group2(n,1)  =   zero
      Group2(n,2)  = + omega1_odelt
   enddo

!  Linearization of Direct Term
!  ----------------------------

!  Linearization of Direct Term, Zero terms in MSMODE-only

   L_t_direct_up_Gp1 = zero
   L_t_direct_up_Gp2 = zero
   L_t_direct_dn_Gp1 = zero
   L_t_direct_dn_Gp2 = zero

!  Skip this section for the sfollowing conditions. [removed Goto 68, 9/17/16]
!   IF ( DO_POSTPROCESSING .and. DO_MSMODE_THERMAL ) go to 68
!   IF ( .not. DO_POSTPROCESSING                   ) go to 68

   IF ( DO_POSTPROCESSING .and..not.DO_MSMODE_THERMAL ) then

!  Upwelling Direct solution WHOLE-LAYER source terms

     IF ( do_upwelling ) THEN
       DO um = 1, n_user_streams
         cosmum = user_streams(um)
         do n = 1, nlayers
           if ( layermask_up(n) ) then
             Udel = t_delt_userm(n,um)
             u1 = one - Udel ; u2 = cosmum - Udel * ( cosmum + deltaus(n) )
             L_t_direct_up_Gp1(um,n) = half * ( u1 * Group1(n,1) + u2 * Group1(n,2) )
             L_t_direct_up_Gp2(um,n) = half * ( u1 * Group2(n,1) + u2 * Group2(n,2) )
           endif
         enddo
       enddo
     endif

!  Downwelling Direct solution WHOLE-LAYER source terms

     IF ( do_dnwelling ) THEN
       DO um = 1, n_user_streams
         cosmum = user_streams(um)
         do n = 1, nlayers
           if ( layermask_dn(n) ) then
             Udel = t_delt_userm(n,um)
             d1 = one - Udel ; d2 = deltaus(n) - cosmum * d1
             L_t_direct_dn_Gp1(um,n) = half * ( d1 * Group1(n,1) + d2 * Group1(n,2) )
             L_t_direct_dn_Gp2(um,n) = half * ( d1 * Group2(n,1) + d2 * Group2(n,2) )
           endif
         enddo
       enddo
     endif

!  End linearization of direct term

   endif

!  Continuation point, removed 9/17/16.
!68 continue

!  Thermal Transmittance only, quadrature solutions
!  ================================================

   if ( do_thermal_transonly ) then

     DO n = 1, nlayers
       DO aa = 1, nstreams
         aa1 = aa + nstreams
         Z = t_delt_disords(aa,n) ; zd = Z * deltaus(n) ; z1 = one - Z ; z1_ok = z1 * quad_streams(aa)
         d2 =  ( deltaus(n) - z1_ok ) ; d1 = z1
         u2 =   ( z1_ok - zd )        ; u1 = z1
!   Bugs, 4/21/16
!            L_t_wupper_Gp1(aa1,n)  = u2 * Group1(n,2) + u1 * Group1(n,1)
!            L_t_wupper_Gp2(aa1,n)  = u2 * Group2(n,2) + u1 * Group2(n,1)
!            L_t_wlower_Gp1(aa,n)   = d2 * Group1(n,2) + d1 * Group1(n,1)
!            L_t_wlower_Gp2(aa,n)   = d2 * Group2(n,2) + d1 * Group2(n,1)
         L_t_wupper_Gp1(aa1,n)  = half * ( u2 * Group1(n,2) + u1 * Group1(n,1) ) * quad_streams(aa)
         L_t_wupper_Gp2(aa1,n)  = half * ( u2 * Group2(n,2) + u1 * Group2(n,1) ) * quad_streams(aa)
         L_t_wlower_Gp1(aa,n)   = half * ( d2 * Group1(n,2) + d1 * Group1(n,1) ) * quad_streams(aa)
         L_t_wlower_Gp2(aa,n)   = half * ( d2 * Group2(n,2) + d1 * Group2(n,1) ) * quad_streams(aa)
       END DO
     END DO

!  GOTO and ENdif removed, 9/17/16
!      GO TO 74
!   endif

!  Thermal Scattering, quadrature solutions
!  =========================================

   else

!  Start Layer loop for solution derivatives
!  -----------------------------------------

!    Comment, 26 March 2014
!      We will have to go back to the originals here because not saved

     do n = 1, nlayers

!  First SOLUTION MATRIX and LU-decomposition  (Same as in Thermal CLSolution)

       DO I = 1, NSTREAMS
         TMAT(I,1:nstreams) = SAB(I,1:nstreams,O11,O11,N) * QUAD_STREAMS(I)
       ENDDO
       CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)
       IF ( INFO .GT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(C3, '(I3)' ) N
         MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
         TRACE   = 'DGETRF CALL FOR H-vector, LAYER '//C3//' LBBF_Jacobians_whole, Call 1'
         STATUS  = VLIDORT_SERIOUS ; return
       ENDIF

!  Linearized H VECTOR_1 Group-1 :  SOLUTION BY BACK-SUBSTITUTION
!  Linearized H VECTOR_1 Group-2 :  Zero

       HVEC = one
       CALL DGETRS  ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,HVEC,MAXSTREAMS,INFO)
       IF ( INFO .GT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(C3, '(I3)' ) N
         MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
         TRACE   = 'DGETRS CALL FOR H-vector, LAYER '//C3//' LBBF_Jacobians_whole, Call 2'
         STATUS  = VLIDORT_SERIOUS ; return
       ENDIF
       HVEC1_Gp1(1:nstreams) = Group1(n,1) * HVEC(1:nstreams)
       HVEC1_Gp2(1:nstreams) = zero

!  Linearized H VECTOR_2 Group-1 :  SOLUTION BY BACK-SUBSTITUTION
!  Linearized H VECTOR_2 Group-2 :  SOLUTION BY BACK-SUBSTITUTION

       HVEC2_Gp1(1:nstreams) = Group1(n,2) * HVEC(1:nstreams)
       HVEC2_Gp2(1:nstreams) = Group2(n,2) * HVEC(1:nstreams)

!  Second SOLUTION MATRIX and LU-decomposition  (Same as in Thermal CLSolution)

       DO I = 1, NSTREAMS
         TMAT(I,1:nstreams) = -DAB(I,1:nstreams,O11,O11,N)
       ENDDO
       CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)
       IF ( INFO .GT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(C3, '(I3)' ) N
         MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
         TRACE   = 'DGETRF CALL FOR J-vector, LAYER '//C3//' LBBF_Jacobians_whole, Call 3'
         STATUS  = VLIDORT_SERIOUS ; return
       ENDIF

!  Linearized J VECTOR_2 Groups 1-2 :  SOLUTION BY BACK-SUBSTITUTION

       CALL DGETRS  ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,HVEC,MAXSTREAMS,INFO)
       IF ( INFO .GT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(C3, '(I3)' ) N
         MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
         TRACE   = 'DGETRS CALL FOR H-vector, LAYER '//C3//' LBBF_Jacobians_whole, Call 4'
         STATUS  = VLIDORT_SERIOUS ; return
       ENDIF
       JVEC1_Gp1(1:nstreams) = Group1(n,2) * HVEC(1:nstreams)
       JVEC1_Gp2(1:nstreams) = Group2(n,2) * HVEC(1:nstreams)
!       if (n.eq.1) write(*,*)'Linear',L,n,HVEC1_Gp2(2), JVEC1_Gp2(2)

!  Set solution

       DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         TVEC1_Gp1(I,N)   = HALF * (HVEC1_Gp1(I) + JVEC1_Gp1(I))
         TVEC1_Gp2(I,N)   = HALF * (HVEC1_Gp2(I) + JVEC1_Gp2(I))
         TVEC1_Gp1(I1,N)  = HALF * (HVEC1_Gp1(I) - JVEC1_Gp1(I))
         TVEC1_Gp2(I1,N)  = HALF * (HVEC1_Gp2(I) - JVEC1_Gp2(I))
       ENDDO
!  This code is Wrong 4/15/16
!      DO I = 1, NSTREAMS_2
!         TVEC2_Gp1(I,N)  = HALF * HVEC2_Gp1(I)
!         TVEC2_Gp2(I,N)  = HALF * HVEC2_Gp2(I)
!      ENDDO
       DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         TVEC2_Gp1(I,N)   = HALF * HVEC2_Gp1(I)
         TVEC2_Gp2(I,N)   = HALF * HVEC2_Gp2(I)
         TVEC2_Gp1(I1,N)  = HALF * HVEC2_Gp1(I)
         TVEC2_Gp2(I1,N)  = HALF * HVEC2_Gp2(I)
       ENDDO

       DO I = 1, NSTREAMS_2
         L_T_WUPPER_Gp1(I,N) = TVEC1_Gp1(I,N)
         L_T_WUPPER_Gp2(I,N) = TVEC1_Gp2(I,N)
!         L_T_WLOWER_Gp1(I,N) = TVEC1_Gp1(I,N) + TVEC2_Gp2(I,N) * DELTAUS(N)   !   Bug, 4/20/16
         L_T_WLOWER_Gp1(I,N) = TVEC1_Gp1(I,N) + TVEC2_Gp1(I,N) * DELTAUS(N)
         L_T_WLOWER_Gp2(I,N) = TVEC1_Gp2(I,N) + TVEC2_Gp2(I,N) * DELTAUS(N)
       ENDDO

!  USER SOLUTIONS
!    4/15/16. Bad bug in the original, see commented out code

       DO L = M, NMOMENTS
         SUM1_Gp1 = ZERO ; SUM2_Gp1 = ZERO
         SUM1_Gp2 = ZERO ; SUM2_Gp2 = ZERO
         HOMIGR = half * OMEGA_GREEK(L,N,O11,O11)
         DO  J = 1, NSTREAMS
           J1 = J + NSTREAMS
           A5_xqp = QUAD_WEIGHTS(J) * PI_XQP    (L,J,O11,O11) * HOMIGR
           A5_xqm = QUAD_WEIGHTS(J) * PI_XQM_PRE(L,J,O11,O11) * HOMIGR
!           PN1_Gp1 = TVEC1_Gp1(J1,N) * A5_xqp + TVEC1_Gp1(J,N) + A5_xqm ;  SUM1_Gp1 = SUM1_Gp1 + PN1_Gp1
!           PN1_Gp2 = TVEC1_Gp2(J1,N) * A5_xqp + TVEC1_Gp2(J,N) + A5_xqm ;  SUM1_Gp2 = SUM1_Gp2 + PN1_Gp2
!           PN2_Gp1 = TVEC2_Gp1(J1,N) * A5_xqp + TVEC2_Gp1(J,N) + A5_xqm ;  SUM2_Gp1 = SUM2_Gp1 + PN2_Gp1
!           PN2_Gp2 = TVEC2_Gp2(J1,N) * A5_xqp + TVEC2_Gp2(J,N) + A5_xqm ;  SUM2_Gp2 = SUM2_Gp2 + PN2_Gp2
           PN1_Gp1 = TVEC1_Gp1(J1,N) * A5_xqp + TVEC1_Gp1(J,N) * A5_xqm ;  SUM1_Gp1 = SUM1_Gp1 + PN1_Gp1
           PN1_Gp2 = TVEC1_Gp2(J1,N) * A5_xqp + TVEC1_Gp2(J,N) * A5_xqm ;  SUM1_Gp2 = SUM1_Gp2 + PN1_Gp2
           PN2_Gp1 = TVEC2_Gp1(J1,N) * A5_xqp + TVEC2_Gp1(J,N) * A5_xqm ;  SUM2_Gp1 = SUM2_Gp1 + PN2_Gp1
           PN2_Gp2 = TVEC2_Gp2(J1,N) * A5_xqp + TVEC2_Gp2(J,N) * A5_xqm ;  SUM2_Gp2 = SUM2_Gp2 + PN2_Gp2
         ENDDO
         T_HELP1_Gp1(L) = SUM1_Gp1 ; T_HELP1_Gp2(L) = SUM1_Gp2
         T_HELP2_Gp1(L) = SUM2_Gp1 ; T_HELP2_Gp2(L) = SUM2_Gp2 
!       if (n.eq.1) write(*,*)'Linear',L,n,T_HELP1_Gp2(L)
       ENDDO

!  UPWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

       IF ( DO_UPWELLING.AND.LAYERMASK_UP(N) ) THEN
         DO UM = 1, N_USER_STREAMS
           POS1 = dot_product(T_HELP1_Gp1(M:NMOMENTS),PI_XUP(M:NMOMENTS,UM,O11,O11) )
           POS2 = dot_product(T_HELP2_Gp1(M:NMOMENTS),PI_XUP(M:NMOMENTS,UM,O11,O11) )
           U_TPOS1_Gp1(UM,N) = POS1 ; U_TPOS2_Gp1(UM,N) = POS2
           POS1 = dot_product(T_HELP1_Gp2(M:NMOMENTS),PI_XUP(M:NMOMENTS,UM,O11,O11) )
           POS2 = dot_product(T_HELP2_Gp2(M:NMOMENTS),PI_XUP(M:NMOMENTS,UM,O11,O11) )
           U_TPOS1_Gp2(UM,N) = POS1 ; U_TPOS2_Gp2(UM,N) = POS2
         ENDDO
       ENDIF

!  DNWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

       IF ( DO_DNWELLING.AND.LAYERMASK_DN(N) ) THEN
         DO UM = 1, N_USER_STREAMS
           NEG1 = dot_product(T_HELP1_Gp1(M:NMOMENTS),PI_XUM(M:NMOMENTS,UM,O11,O11) )
           NEG2 = dot_product(T_HELP2_Gp1(M:NMOMENTS),PI_XUM(M:NMOMENTS,UM,O11,O11) )
           U_TNEG1_Gp1(UM,N) = NEG1 ; U_TNEG2_Gp1(UM,N) = NEG2
           NEG1 = dot_product(T_HELP1_Gp2(M:NMOMENTS),PI_XUM(M:NMOMENTS,UM,O11,O11) )
           NEG2 = dot_product(T_HELP2_Gp2(M:NMOMENTS),PI_XUM(M:NMOMENTS,UM,O11,O11) )
           U_TNEG1_Gp2(UM,N) = NEG1 ; U_TNEG2_Gp2(UM,N) = NEG2
!       if (n.eq.1) write(*,*)'Linear',UM,n,U_TPOS1_Gp2(UM,N),U_TNEG1_Gp2(UM,N)
         ENDDO
       ENDIF

!  End layer loop

     ENDDO

!  End clause, thermal transonly vs. Scattering

   endif

!  Continuation point for thermal tranmsittance only solutions (removed, 9/17/16
!74 continue

!  Layer source term derivatives
!  =============================

!  Initialize completely, skip if no post processing

   L_LAYER_TSUP_UP_Gp1 = zero
   L_LAYER_TSUP_DN_Gp1 = zero
   L_LAYER_TSUP_UP_Gp2 = zero
   L_LAYER_TSUP_DN_Gp2 = zero

!  Postprocessing
!  ==============

   if ( do_POSTPROCESSING ) then

!   GOTO Removed 9/17/16.
!   if ( .not. do_POSTPROCESSING ) go to 69

!  Initialize with Direct term (which may be zero...)
!  --------------------------------------------------

     DO um = 1, n_user_streams
       do n = 1, nlayers
         if ( do_upwelling .and. layermask_up(n) ) then
           L_layer_tsup_up_Gp1(um,n) = fac * L_t_direct_up_Gp1(um,n)
           L_layer_tsup_up_Gp2(um,n) = fac * L_t_direct_up_Gp2(um,n)
         endif
         if ( do_dnwelling .and. layermask_dn(n) ) then
           L_layer_tsup_dn_Gp1(um,n) = fac * L_t_direct_dn_Gp1(um,n)
           L_layer_tsup_dn_Gp2(um,n) = fac * L_t_direct_dn_Gp2(um,n)
         endif
       enddo
     enddo

!  Not done if thermal Transmittance only
!   GOTO Removed 9/17/16.
!   if ( do_thermal_transonly ) go to 69

     if ( .not.do_thermal_transonly ) then

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!  --------------------------------------------------

       DO UM = 1, N_USER_STREAMS
         COSMUM = USER_STREAMS(UM)
         do n = 1, nlayers

           Udel = t_delt_userm(n,um) ;  delUdel = deltaus(n) * Udel
           u1 = one - udel ; u2 = cosmum * u1 - delUdel

!  Upwelling

           if ( do_upwelling .and. layermask_up(n) ) then
             u1 = one - udel ; u2 = cosmum * u1 - delUdel
             spar1 = u2 * U_TPOS2_Gp1(UM,N) + u1 * U_TPOS1_Gp1(UM,N)
             spar2 = u2 * U_TPOS2_Gp2(UM,N) + u1 * U_TPOS1_Gp2(UM,N)
             L_layer_tsup_up_Gp1(um,n) = L_layer_tsup_up_Gp1(um,n) + spar1 * fac
             L_layer_tsup_up_Gp2(um,n) = L_layer_tsup_up_Gp2(um,n) + spar2 * fac
           endif

!  Downwelling 

           if ( do_dnwelling .and. layermask_dn(n) ) then
!             d1 = one - udel ; d2 = - cosmum * d1 - deltaus(n)         !  BUG, Confirmed 4/20/16
             d1 = one - udel ; d2 = - cosmum * d1 + deltaus(n)
             spar1 = d2 * U_TNEG2_Gp1(UM,N) + d1 * U_TNEG1_Gp1(UM,N)
             spar2 = d2 * U_TNEG2_Gp2(UM,N) + d1 * U_TNEG1_Gp2(UM,N)
             L_layer_tsup_dn_Gp1(um,n) = L_layer_tsup_dn_Gp1(um,n) + spar1 * fac
             L_layer_tsup_dn_Gp2(um,n) = L_layer_tsup_dn_Gp2(um,n) + spar2 * fac
           endif

!  End layer and User-stream loops

         enddo
       ENDDO

!  End thermal scattering clause

     endif

!  End postprocessing clause

   endif

!  Continuation point. Removed 9/17/16.
!69 continue

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        START MAIN LOOP OVER LEVEL BBFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   DO JB = 0, NLAYERS

!  Solve BVProblem
!  ===============

!  Vector upgrade
!    -- Although the LBBF particular thermal solution is scalar, BRDF may induce polarization
!    -- Additional polarization (I/Q only) comes through the MS solutions, after solving BVP

!  Initialize

     N = JB + 1 ; N1 = JB
     COL2_BWF = zero

!  Skip BVP for tranmsittance only
!   GOTO Removed 9/17/16.
!   if ( do_thermal_transonly ) go to 75

     if ( .not.do_thermal_transonly ) then

!  Surface terms. Down = surface downwelling dependence
!    -- 1/31/21. Version 2.8.3. BRDF arrays are defined locally, drop "M" Fourier index
 
       Reflec = zero
       IF ( DO_INCLUDE_SURFACE.and.JB.ge.nlayers-1 ) THEN
         Down = zero
         if ( jb .eq. nlayers ) then
           do j = 1, nstreams
              Down(j) = L_T_WLOWER_Gp2(j,n1) * QUAD_STRMWTS(J)
           enddo
         else if (jb .eq. nlayers - 1 ) then
           do j = 1, nstreams
              Down(j) = L_T_WLOWER_Gp1(j,n) * QUAD_STRMWTS(J)
           enddo
         endif
         IF ( DO_LAMBERTIAN_SURFACE  ) THEN
           FACTOR = SURFACE_FACTOR * ALBEDO * sum(Down(1:nstreams))
           Reflec(1:nstreams,O11) = FACTOR
         ELSE
           do i = 1, nstreams
             do o2 = 1, nstokes
               om = mueller_index(o2,o11)
               FACTOR = SURFACE_FACTOR * Dot_Product(Down(1:nstreams),brdf_f(om,i,1:nstreams))
               Reflec(i,o2) = FACTOR
             enddo
           enddo
         ENDIF
       ENDIF

!  BVProblem, Special Case, N = 2

       if ( nlayers .eq. 2 ) then
         if ( JB.eq.0 ) then                 ! Correct 3/19
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1
             COL2_BWF(ir) = - L_T_WUPPER_Gp1(i,n)
           enddo
           CM = nstks_nstrms
           do i = 1, nstreams_2
             ir = nstokes*(i-1) + 1 ; ic = cm + ir
             COL2_BWF(ic)  = - L_T_WLOWER_Gp1(i,n)
           enddo
         else if ( JB.eq.1 ) then            ! Correct 3/19
           cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms ; cmp = cm + nstks_nstrms_2
           do i = 1, nstreams_2
             ir = nstokes*(i-1) + 1 ; ic = cm + ir ; ic1 = cmp + ir
             COL2_BWF(ic1)  = L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
           enddo
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1 ; ic = cmp + ir + nstks_nstrms_2 ; i1 = i + nstreams
             COL2_BWF(ir) = - L_T_WUPPER_Gp2(i,n1)
             COL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i,o11)
             if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
               do o2 = 2, nstokes
                 ic2 = ic + (o2-1)
                 COL2_BWF(ic2) = Reflec(i,o2)
               enddo
             endif
           enddo
         else if ( JB.eq.NLAYERS ) then     ! Correct 3/19
           cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms ; cmp = cm + nstks_nstrms_2
           do i = 1, nstreams_2
             ir = nstokes*(i-1) + 1 ; ic = cm + ir
             COL2_BWF(ic)  = + L_T_WUPPER_Gp1(i,n1)
           enddo
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1 ; ic = cmp + ir ; i1 = i + nstreams
             COL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i,o11)
             if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
               do o2 = 2, nstokes
                 ic2 = ic + (o2-1)
                 COL2_BWF(ic2) = Reflec(i,o2)
               enddo
             endif
           enddo
         endif
       endif

!  BVProblem for a single layer.
!   Bug 4/9/19. Special Case, N = 1 (JB = 0), N1 = 1 (JB = 1)
       
       if ( nlayers .eq. 1 ) then
         if ( JB.eq.0 ) then
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1 ; i1 = i + nstreams ; ic = ir + nstks_nstrms
             SCOL2_BWF(ir) = - L_T_WUPPER_Gp1(i,n)
             SCOL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i,o11)
             if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
               do o2 = 2, nstokes
                 ic2 = ic + o2 - 1
                 SCOL2_BWF(ic2) = Reflec(i,o2)
               enddo
             endif
           enddo
         else if ( JB.eq.1 ) then
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1 ; i1 = i + nstreams ; ic = ir + nstks_nstrms
             SCOL2_BWF(ir) = - L_T_WUPPER_Gp2(i,n1)                     ! bug fix n --> n1
             SCOL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i,o11)    ! bug fix n --> n1
             if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
               do o2 = 2, nstokes
                 ic2 = ic + o2 - 1
                 SCOL2_BWF(ic2) = Reflec(i,o2)
               enddo
             endif
           enddo
         endif
       endif

!  General Case, N > 2. Separate out the various cases

       if ( nlayers .gt. 2 ) then
         if ( JB.eq.0 ) then                 ! Correct 3/19
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1
             COL2_BWF(ir) = - L_T_WUPPER_Gp1(i,n)
           enddo
           CM = nstks_nstrms
           do i = 1, nstreams_2
             ir = nstokes*(i-1) + 1 ; ic = cm + ir
             COL2_BWF(ic)  = - L_T_WLOWER_Gp1(i,n)
           enddo
         else if ( JB.eq.1 ) then            ! Correct 3/19
           cm = nstks_nstrms
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1
             COL2_BWF(ir) = - L_T_WUPPER_Gp2(i,n1)
           enddo
           do i = 1, nstreams_2
             ir = nstokes*(i-1) + 1 ; ic = cm + ir ; ic1 = ic + nstks_nstrms_2
             COL2_BWF(ic)   = + L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
             COL2_BWF(ic1)  = - L_T_WLOWER_Gp1(i,n)
           enddo
         else if ( JB.eq.nlayers - 1 ) then  ! Correct 3/19
           cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms ; cmp = cm + nstks_nstrms_2
           do i = 1, nstreams_2
             ir = nstokes*(i-1) + 1 ; ic = cm + ir ; ic1 = cmp + ir
             COL2_BWF(ic)   = L_T_WUPPER_Gp2(i,n1)
             COL2_BWF(ic1)  = L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
           enddo
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1 ; i1 = i + nstreams ; ic = cmp + ir + nstks_nstrms_2
             COL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i,o11)
             if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
               do o2 = 2, nstokes
                 ic2 = ic + (o2-1)
                 COL2_BWF(ic2) = Reflec(i,o2)
               enddo
             endif
           enddo
         else if ( JB.eq.NLAYERS ) then     ! Correct 3/19
           cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms ; cmp = cm + nstks_nstrms_2
           do i = 1, nstreams_2
             ir = nstokes*(i-1) + 1 ; ic = cm + ir
             COL2_BWF(ic)  = + L_T_WUPPER_Gp2(i,n1)
           enddo
           do i = 1, nstreams
             ir = nstokes*(i-1) + 1 ; i1 = i + nstreams ; ic = cmp + ir
             COL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i,o11)
             if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
               do o2 = 2, nstokes
                 ic2 = ic + (o2-1)
                 COL2_BWF(ic2) = Reflec(i,o2)
               enddo
             endif
           enddo
         else                               ! Correct 3/19
           cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms
           do i = 1, nstreams_2
             ir = nstokes*(i-1) + 1 ; ic = cm + ir ; ic1 = ic + nstks_nstrms_2 ; ic2 = ic1 + nstks_nstrms_2
             COL2_BWF(ic)   = + L_T_WUPPER_Gp2(i,n1)
             COL2_BWF(ic1)  = + L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
             COL2_BWF(ic2)  = - L_T_WLOWER_Gp1(i,n)
           enddo
         endif
       endif

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

       IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

         CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_BWF, MAXTOTAL, INFO )

!  Exception handling

         IF ( INFO .LT. 0 ) THEN
           WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
           MESSAGE = 'argument i illegal value, for i = '//CI
           TRACE   = ' for Atmos BBF Level '//CN//' DGBTRS call in LBBF_Jacobians'
           STATUS  = VLIDORT_SERIOUS ; RETURN
         ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

         DO N = 1, NLAYERS
           C0 = ( N - 1 ) * NSTKS_NSTRMS_2
           KO1 = K_REAL(N) + 1
           DO K = 1, K_REAL(N)
             IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
             NCON(K,N) = COL2_BWF(C0+IROW)
             PCON(K,N) = COL2_BWF(C0+IROW1)
           ENDDO
           DO K = 1, K_COMPLEX(N)
             K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
             IROW    = K    + K_REAL(N)       ; IROW1   = IROW   + NSTKS_NSTRMS
             IROW_S  = IROW + K_COMPLEX(N)    ; IROW1_S = IROW_S + NSTKS_NSTRMS
             NCON(K1,N) = COL2_BWF(C0+IROW)   ; NCON(K2,N) = COL2_BWF(C0+IROW_S)
             PCON(K1,N) = COL2_BWF(C0+IROW1)  ; PCON(K2,N) = COL2_BWF(C0+IROW1_S)
           ENDDO
         ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

       ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

         CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
                       SCOL2_BWF, MAXSTRMSTKS_2, INFO )

!  Exception handling

         IF ( INFO .LT. 0 ) THEN
           WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
           MESSAGE = 'argument i illegal value, for i = '//CI
           TRACE   = ' for BBF Level '//CN//' DGBTRS call in 1-layer LBBF_Jacobians'
           STATUS  = VLIDORT_SERIOUS ; RETURN
         ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

         N = 1 ;  KO1 = K_REAL(N) + 1
         DO K = 1, K_REAL(N)
           IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
           NCON(K,N) = SCOL2_BWF(IROW)
           PCON(K,N) = SCOL2_BWF(IROW1)
         ENDDO
         DO K = 1, K_COMPLEX(N)
           K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
           IROW    = K    + K_REAL(N)       ; IROW1   = IROW   + NSTKS_NSTRMS
           IROW_S  = IROW + K_COMPLEX(N)    ; IROW1_S = IROW_S + NSTKS_NSTRMS
           NCON(K1,N) = SCOL2_BWF(IROW)     ; NCON(K2,N) = SCOL2_BWF(IROW_S)
           PCON(K1,N) = SCOL2_BWF(IROW1)    ; PCON(K2,N) = SCOL2_BWF(IROW1_S)
         ENDDO

!  End layer clause

        ENDIF

!  End thermal scattering clause

     endif
  
!  Continuation point. Removed 9/17/16.
!75 continue

!  Upwelling Jacobians
!  ===================

!  Skip if not applicable. GOTO 344, Removed 9/17/16
!      IF ( .not. DO_POSTPROCESSING ) GO TO 344
!      IF ( .NOT. DO_UPWELLING )      GO TO 344

     IF ( DO_POSTPROCESSING .and. DO_UPWELLING ) THEN

!  Zero BOA terms

       L_BOA_MSSOURCE       = zero
       L_BOA_THTONLY_SOURCE = zero

!  Skip BOA terms if no surface. Rmoved GOTO, 9/17/16
!       IF ( .not. DO_INCLUDE_SURFACE ) go to 76

       IF ( DO_INCLUDE_SURFACE ) THEN

!  Reflected  downwelling solution 
!  -------------------------------

!   Distinguish between thermal transmittance only, and scattered solution

         IF ( DO_THERMAL_TRANSONLY ) THEN
           L_IDOWN = zero ; O1 = 1
           do n = 1, nlayers
             L_IDOWN(1:nstreams,o1) = L_IDOWN(1:nstreams,o1) * T_DELT_DISORDS(1:nstreams,N)
             if ( JB.eq.n )     L_IDOWN(1:nstreams,o1) = L_IDOWN(1:nstreams,o1) + L_T_WLOWER_GP2(1:nstreams,N)
             if ( JB.eq.n - 1 ) L_IDOWN(1:nstreams,o1) = L_IDOWN(1:nstreams,o1) + L_T_WLOWER_GP1(1:nstreams,N)
           enddo
           L_IDOWN(1:nstreams,o1) = L_IDOWN(1:nstreams,o1) * quad_weights(1:nstreams)
         ELSE
           N = NLAYERS
           do i = 1, nstreams
             do o1 = 1, nstokes
!  ...Zero contributions 
               SPAR = zero ; SHOM_R = zero ; SHOM_CR = zero
!  ..thermal particular integral; only scalar component
               if ( o1.eq.1.and.JB.eq.nlayers )     SPAR = L_T_WLOWER_GP2(i,N)
               if ( o1.eq.1.and.JB.eq.nlayers - 1 ) SPAR = L_T_WLOWER_GP1(i,N)
!  .. Homogeneous Real part
               DO K = 1, K_REAL(N)
                 HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N)* T_DELT_EIGEN(K,N)
                 HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)
                 SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO
!  .. Homogeneous Complex part
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                 K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                 NXR1    = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                 NXR2    = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                 HOM2CR  = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                 HOM1CR =  NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
                 SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO
!  .. Gather solution
               SHOM = SHOM_R + SHOM_CR
               L_IDOWN(I,O1) = ( SPAR + SHOM ) * QUAD_STRMWTS(I)
             ENDDO
           ENDDO
         ENDIF

!  BOA MS source terms
!  -------------------

!  Polarized calculation for the BRDF surface.
!    -- 1/31/21. Version 2.8.3. BRDF arrays are defined locally, drop "M" Fourier index

         IF ( DO_LAMBERTIAN_SURFACE ) THEN
           FACTOR = SURFACE_FACTOR * ALBEDO * SUM(L_IDOWN(1:nstreams,o11))
           L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) = FACTOR
           IF ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:NSTREAMS) =  FACTOR
         ELSE
           DO UM = 1, N_USER_STREAMS
             do o1 = 1, nstokes
               FACTOR = ZERO
               DO J = 1, NSTREAMS
                 S_REFLEC = ZERO
                 DO O2 = 1, NSTOKES
                   OM = MUELLER_INDEX(O1,O2)
                   S_REFLEC = S_REFLEC + L_IDOWN(J,O2) * UBRDF_F(OM,UM,J)
                 ENDDO
                 FACTOR = FACTOR + S_REFLEC
               ENDDO
               L_BOA_MSSOURCE(UM,o1) = SURFACE_FACTOR * FACTOR
             enddo
           ENDDO
           IF ( DO_QTHTONLY ) THEN
             DO I = 1, NSTREAMS
               FACTOR = DOT_PRODUCT(L_IDOWN(1:nstreams,O11),BRDF_F(O11,I,1:NSTREAMS))
               L_BOA_THTONLY_SOURCE(I) =  SURFACE_FACTOR * FACTOR
             ENDDO
           ENDIF
         ENDIF

!  End surface inclusion

        ENDIF

!  continuation point. removed, 9/17/16
!76    continue

!  Recursion Loop for linearized Post-processed Jacobians (upwelling)
!  ------------------------------------------------------------------

!  Set the cumulative source term equal to BOA values

       L_CUMULSOURCE = zero
       DO UM = 1, N_USER_STREAMS
         L_CUMULSOURCE(UM,1:nstokes) = L_BOA_MSSOURCE(UM,1:nstokes) 
       ENDDO

!  initialise cumulative source term loop

       NC  = 0
       NUT = 0
       NSTART = NLAYERS
       NUT_PREV = NSTART + 1

!  loop over all output optical depths

       DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

         NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
         NUT = NLEVEL + 1

!  Cumulative lop to include source terms, to layer NUT

         DO N = NSTART, NUT, -1
           NC = NLAYERS + 1 - N

!  Homogeneous Solution Contributions (only present with scattered light)
   
           L_LAYERSOURCE = zero
           if ( .not. do_thermal_transonly ) then
             DO UM = 1, N_USER_STREAMS
               DO o1 = 1, nstokes
!  ....Real homogeneous solutions
                 SHOM_R = ZERO
                 DO K = 1, K_REAL(N)
                   H1 = NCON(K,N) * UHOM_UPDN(UM,O1,K,N) * HMULT_2(K,UM,N)
                   H2 = PCON(K,N) * UHOM_UPUP(UM,O1,K,N) * HMULT_1(K,UM,N)
                   SHOM_R = SHOM_R + H1 + H2
                 ENDDO
!  ....Complex homogeneous solutions
                 SHOM_CR = ZERO
                 DO K = 1, K_COMPLEX(N)
                   K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                   NUXR1 = NCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                   NUXR2 = NCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                   PUXR1 = PCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                   PUXR2 = PCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                   H1 = NUXR1 * HMULT_2(K1,UM,N) - NUXR2 * HMULT_2(K2,UM,N)
                   H2 = PUXR1 * HMULT_1(K1,UM,N) - PUXR2 * HMULT_1(K2,UM,N)
                   SHOM_CR = SHOM_CR + H1 + H2
                 ENDDO
                 L_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR
               ENDDO
             ENDDO
           endif

!  Add thermal emission source terms (Includes direct and diffuse). Unpolarized
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     -----Only present for layers N adjacent to the level JB that is varying

           o1 = 1 ; TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
           if ( N.eq.JB + 1 ) then
             DO UM = 1, N_USER_STREAMS
               L_LAYERSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + L_LAYER_TSUP_UP_Gp1(UM,N)*TM
             ENDDO
           else if ( N.eq.JB ) then
             DO UM = 1, N_USER_STREAMS
               L_LAYERSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + L_LAYER_TSUP_UP_Gp2(UM,N)*TM
             ENDDO
           endif

!  Add to Linearized cumulative source sterm

           DO UM = 1, N_USER_STREAMS
             DO o1 = 1, nstokes
               L_CUMULSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,o1)
             ENDDO
           ENDDO

!  End layer recursion loop

         ENDDO

!  User-defined stream output, just set to the cumulative source term

         DO UM = 1, N_USER_STREAMS
           DO o1 = 1, nstokes
             ABBWFS_JACOBIANS(UTA,UM,JB,O1,UPIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,O1)
           ENDDO
         ENDDO

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT

!  end loop over optical depth

       ENDDO

!  End Upwelling Postprocessing clause

     endif

!  continuation point. removed 9/17/16
!344   continue

!  Downwelling Jacobians
!  =====================

!  Skip if not applicable. GOTO 34,5 Removed 9/17/16
!      IF ( .not. DO_POSTPROCESSING ) GO TO 345
!      IF ( .NOT. DO_DNWELLING )      GO TO 345

     IF ( DO_POSTPROCESSING .and. DO_DNWELLING ) THEN

!  Initialize post-processing recursion
!  Set the cumulative source term equal to TOA values

       L_CUMULSOURCE = zero

!  Recursion Loop for linearized Post-processed Jacobians (Downwelling)
!  --------------------------------------------------------------------

!  initialise cumulative source term loop

       NC  = 0
       NUT = 0
       NSTART = 1
       NUT_PREV = NSTART - 1

!  loop over all output optical depths

       DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

         NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
 
         NUT = NLEVEL
         DO N = NSTART, NUT
           NC = N

!  Homogeneous Solution Contributions (only present with scattered light)
   
           L_LAYERSOURCE = zero
           if ( .not. do_thermal_transonly ) then
             DO UM = 1, N_USER_STREAMS
               DO o1 = 1, nstokes
!  ....Real homogeneous solutions
                 SHOM_R = ZERO
                 DO K = 1, K_REAL(N)
                   H1 = NCON(K,N) * UHOM_DNDN(UM,O1,K,N) * HMULT_1(K,UM,N)
                   H2 = PCON(K,N) * UHOM_DNUP(UM,O1,K,N) * HMULT_2(K,UM,N)
                   SHOM_R = SHOM_R + H1 + H2
                 ENDDO
!  ....Complex homogeneous solutions
                 SHOM_CR = ZERO
                 DO K = 1, K_COMPLEX(N)
                   K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                   NUXR1 = NCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
                   NUXR2 = NCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
                   PUXR1 = PCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
                   PUXR2 = PCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_DnUP(UM,O1,K1,N)
                   H1 = NUXR1 * HMULT_1(K1,UM,N) - NUXR2 * HMULT_1(K2,UM,N)
                   H2 = PUXR1 * HMULT_2(K1,UM,N) - PUXR2 * HMULT_2(K2,UM,N)
                   SHOM_CR = SHOM_CR + H1 + H2
                 ENDDO
                 L_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR
               ENDDO
             ENDDO
           endif

!  Add thermal emission source terms (Includes direct and diffuse). Unpolarized
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     -----Only present for layers N adjacent to the level JB that is varying

           o1 = 1 ; TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
           if ( N.eq.JB + 1 ) then
             DO UM = 1, N_USER_STREAMS
               L_LAYERSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + L_LAYER_TSUP_DN_Gp1(UM,N)*TM
             ENDDO
           else if ( N.eq.JB ) then
             DO UM = 1, N_USER_STREAMS
               L_LAYERSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + L_LAYER_TSUP_DN_Gp2(UM,N)*TM
             ENDDO
           endif

!  Add to Linearized cumulative source sterm

           DO UM = 1, N_USER_STREAMS
             DO o1 = 1, nstokes
               L_CUMULSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,o1)
             ENDDO
           ENDDO

!  End layer loop

         ENDDO

!  User-defined stream output, just set to the cumulative source term

         DO UM = 1, N_USER_STREAMS
           DO o1 = 1, nstokes
             ABBWFS_JACOBIANS(UTA,UM,JB,O1,DNIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,O1)
           ENDDO
         ENDDO

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT

!  end loop over optical depth

       ENDDO

!  End downwelling postprocessing

     endif

!  continuation point. Removed 9/17/16
!345   continue

!  Flux Jacobians
!  ==============

!  Upwelling FLux output
!  ---------------------

      if ( DO_MVOUTPUT .and. DO_UPWELLING ) THEN
         DO UTA = 1, N_USER_LEVELS
            NL = UTAU_LEVEL_MASK_UP(UTA) ; N = NL + 1 ; BBWF_QUAD = zero ; o11 = 1

!  quadrature field at bottom level

            IF ( NL .EQ. NLAYERS  ) THEN
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     BBWF_QUAD(I,o11) = FLUX_MULTIPLIER * L_BOA_THTONLY_SOURCE(I)
                  enddo
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS ; spar = zero
                     DO O1 = 1, nstokes
                        SHOM_R = ZERO ; SHOM_CR = zero
                        DO K = 1, K_REAL(NL)
                           HOM1 = NCON(K,NL) * SOLA_XPOS(I1,O1,K,NL) * T_DELT_EIGEN(K,NL)
                           HOM2 = PCON(K,NL) * SOLB_XNEG(I1,O1,K,NL)
                           SHOM_R = SHOM_R + HOM1 + HOM2
                        ENDDO
                        DO K = 1, K_COMPLEX(NL)
                           K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                           NXR1   = NCON(K1,NL) * SOLA_XPOS(I1,O1,K1,NL) - NCON(K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
                           NXR2   = NCON(K1,NL) * SOLA_XPOS(I1,O1,K2,NL) + NCON(K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
                           HOM2CR = PCON(K1,NL) * SOLB_XNEG(I1,O1,K1,NL) - PCON(K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
                           HOM1CR = NXR1 * T_DELT_EIGEN(K1,NL) - NXR2 * T_DELT_EIGEN(K2,NL)
                           SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                        ENDDO
                        BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                     ENDDO
                     if ( JB.eq.NL-1 ) SPAR = L_T_WLOWER_Gp1(I1,NL)
                     if ( JB.eq.NL )   SPAR = L_T_WLOWER_Gp2(I1,NL)
                     BBWF_QUAD(I,o11) = BBWF_QUAD(I,o11) + FLUX_MULTIPLIER * SPAR 
                  ENDDO
               endif

!  Quadrature field other levels

            ELSE
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     L_THELP = L_BOA_THTONLY_SOURCE(I)
                     DO LAY = NLAYERS, N, -1
                        SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                        if ( JB.eq.LAY-1 ) SPAR = L_T_WUPPER_Gp1(I1,LAY)
                        if ( JB.eq.LAY )   SPAR = L_T_WUPPER_Gp2(I1,LAY)
                        L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                     ENDDO
                     BBWF_QUAD(I,o11) = FLUX_MULTIPLIER * L_THELP
                  enddo
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS ; spar = zero
                     DO O1 = 1, nstokes
                        SHOM_R = ZERO ; SHOM_CR = zero
                        DO K = 1, K_REAL(N)
                           HOM1 = NCON(K,N) * SOLA_XPOS(I1,O1,K,N)
                           HOM2 = PCON(K,N) * SOLB_XNEG(I1,O1,K,N) * T_DELT_EIGEN(K,N)
                           SHOM_R = SHOM_R + HOM1 + HOM2
                        ENDDO
                        DO K = 1, K_COMPLEX(N)
                           K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                           HOM1CR = NCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                           PXR1   = PCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - PCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                           PXR2   = PCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
                           HOM2CR = PXR1 * T_DELT_EIGEN(K1,N) - PXR2 * T_DELT_EIGEN(K2,N)
                           SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                        ENDDO
                        BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                     ENDDO
                     if ( JB.eq.N-1 ) SPAR = L_T_WUPPER_Gp1(I1,N)
                     if ( JB.eq.N )   SPAR = L_T_WUPPER_Gp2(I1,N)
                     BBWF_QUAD(I,o11) = BBWF_QUAD(I,o11) + FLUX_MULTIPLIER * SPAR 
                  ENDDO
               ENDIF
            ENDIF

!  Set fluxes

            do o1 = 1, nstokes
               SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_WEIGHTS(1:NSTREAMS))
               SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS))
               ABBWFS_FLUXES(UTA,1,JB,O1,UPIDX) = HALF * SM
               ABBWFS_FLUXES(UTA,2,JB,O1,UPIDX) = PI2  * SF
            enddo

!  End level output loop

         ENDDO
      ENDIF

!  Downwelling FLux output
!  -----------------------

      if ( DO_MVOUTPUT .and. DO_DNWELLING ) THEN
         DO UTA = 1, N_USER_LEVELS
            NL = UTAU_LEVEL_MASK_DN(UTA) ; N = NL ; BBWF_QUAD = ZERO ; o11 = 1

            IF ( NL .NE. 0  ) THEN
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     L_THELP = ZERO
                     DO LAY = 1, NL
                        SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                        if ( JB.eq.LAY-1 ) SPAR = L_T_WLOWER_Gp1(I,LAY)
                        if ( JB.eq.LAY )   SPAR = L_T_WLOWER_Gp2(I,LAY)
                        L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                     ENDDO
                     BBWF_QUAD(I,o11) = FLUX_MULTIPLIER * L_THELP
                  enddo
               else
                  DO I = 1, NSTREAMS
                     spar = zero
                     DO O1 = 1, nstokes
                        SHOM_R = ZERO ; SHOM_CR = zero
                        DO K = 1, K_REAL(N)
                           HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N) * T_DELT_EIGEN(K,N)
                           HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)
                           SHOM_R = SHOM_R + HOM1 + HOM2
                        ENDDO
                        DO K = 1, K_COMPLEX(N)
                           K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                           NXR1   = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                           NXR2   = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                           HOM2CR = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                           HOM1CR = NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
                           SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                        ENDDO
                        BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                     ENDDO
                     if ( JB.eq.N-1 ) SPAR = L_T_WLOWER_Gp1(I,N)
                     if ( JB.eq.N )   SPAR = L_T_WLOWER_Gp2(I,N)
                     BBWF_QUAD(I,o11) = BBWF_QUAD(I,o11) + FLUX_MULTIPLIER * SPAR 
                  ENDDO
               ENDIF
            ENDIF

!  Set fluxes

            do o1 = 1, nstokes
               SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_WEIGHTS(1:NSTREAMS))
               SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS))
               ABBWFS_FLUXES(UTA,1,JB,O1,DNIDX) = HALF * SM
               ABBWFS_FLUXES(UTA,2,JB,O1,DNIDX) = PI2  * SF
            enddo

!  End level output loop

         ENDDO
      ENDIF

!  End loop over BBWFS

   enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      SURFACE LBBF JACOBIANS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Continuation point. This is the only one left in, 9/17/16.

55 continue

!  Finish if not required

   if ( .not. DO_SURFACE_LBBF ) RETURN

!  Solve BVProblem
!  ===============

!  Skip BVP for thermal Transmittance only
!   GOTO Removed 9/17/16.
!   if ( do_thermal_transonly ) go to 79

   if ( .not.do_thermal_transonly ) then

!  Initialize and Set the Column vector

     if ( nlayers .gt.1 ) then
       COL2_BWF = zero ; C0 = NLAYERS * NSTKS_NSTRMS_2 - NSTKS_NSTRMS
       DO I = 1, NSTREAMS
         DO O1 = 1, nstokes
           IR = NSTOKES*(I-1) + O1 ; CM = C0 + IR
           COL2_BWF(CM) = EMISSIVITY(O1,I)
         ENDDO
       ENDDO
     else
       SCOL2_BWF = zero ; C0 = NSTKS_NSTRMS
       DO I = 1, NSTREAMS
         DO O1 = 1, nstokes
           IR = NSTOKES*(I-1) + O1 ; CM = C0 + IR
           SCOL2_BWF(CM) = EMISSIVITY(O1,I)
         ENDDO
       ENDDO
     endif

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

     IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

       CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_BWF, MAXTOTAL, INFO )

!  Exception handling

       IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = ' for Surface BBF, DGBTRS call in LBBF_Jacobians'
         STATUS  = VLIDORT_SERIOUS ; RETURN
       ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

       DO N = 1, NLAYERS
         C0 = ( N - 1 ) * NSTKS_NSTRMS_2
         KO1 = K_REAL(N) + 1
         DO K = 1, K_REAL(N)
           IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
           NCON(K,N) = COL2_BWF(C0+IROW)
           PCON(K,N) = COL2_BWF(C0+IROW1)
         ENDDO
         DO K = 1, K_COMPLEX(N)
           K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
           IROW    = K    + K_REAL(N)       ; IROW1   = IROW   + NSTKS_NSTRMS
           IROW_S  = IROW + K_COMPLEX(N)    ; IROW1_S = IROW_S + NSTKS_NSTRMS
           NCON(K1,N) = COL2_BWF(C0+IROW)   ; NCON(K2,N) = COL2_BWF(C0+IROW_S)
           PCON(K1,N) = COL2_BWF(C0+IROW1)  ; PCON(K2,N) = COL2_BWF(C0+IROW1_S)
         ENDDO
       ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

     ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

       CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
                       SCOL2_BWF, MAXSTRMSTKS_2, INFO )

!  Exception handling

       IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = ' for Surface BBF, DGBTRS call in 1-layer LBBF_Jacobians'
         STATUS  = VLIDORT_SERIOUS ; RETURN
       ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

       N = 1 ;  KO1 = K_REAL(N) + 1
       DO K = 1, K_REAL(N)
         IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
         NCON(K,N) = SCOL2_BWF(IROW)
         PCON(K,N) = SCOL2_BWF(IROW1)
       ENDDO
       DO K = 1, K_COMPLEX(N)
         K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
         IROW    = K    + K_REAL(N)       ; IROW1   = IROW   + NSTKS_NSTRMS
         IROW_S  = IROW + K_COMPLEX(N)    ; IROW1_S = IROW_S + NSTKS_NSTRMS
         NCON(K1,N) = SCOL2_BWF(IROW)     ; NCON(K2,N) = SCOL2_BWF(IROW_S)
         PCON(K1,N) = SCOL2_BWF(IROW1)    ; PCON(K2,N) = SCOL2_BWF(IROW1_S)
       ENDDO

!  End choice over number of layers

     ENDIF

!  Clause for ending BVP calculation

   endif

!  Continuation point for skipping BVP (removed, 9/17/16)
!79 continue

!  Upwelling Jacobians
!  ===================

!  Skip if not applicable. GOTO 544, Removed 9/17/16
!      IF ( .not. DO_POSTPROCESSING ) GO TO 544
!      IF ( .NOT. DO_UPWELLING )      GO TO 544

   IF ( DO_POSTPROCESSING .and. DO_UPWELLING ) THEN

!  BOA source terms

     L_BOA_MSSOURCE = zero
     IF ( DO_INCLUDE_SURFACE ) THEN
       N = NLAYERS
       L_IDOWN              = zero ! L_Down is zero for thermal transmittance only
       L_BOA_THTONLY_SOURCE = zero ! Only non-zero if thermal transmittance only (and DO_QTHTONLY)

       if ( .not. do_thermal_transonly ) then
         do i = 1, nstreams
           do o1 = 1, nstokes
             SHOM_R = zero ; SHOM_CR = zero
!  .. Homogeneous Real part
             DO K = 1, K_REAL(N)
               HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N)* T_DELT_EIGEN(K,N)
               HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)
               SHOM_R = SHOM_R + HOM1 + HOM2
             ENDDO
!  .. Homogeneous Complex part
             KO1 = K_REAL(N) + 1
             DO K = 1, K_COMPLEX(N)
               K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
               NXR1    = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
               NXR2    = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
               HOM2CR  = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
               HOM1CR =  NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
               SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
             ENDDO
!  .. Gather solution
             SHOM = SHOM_R + SHOM_CR
             L_IDOWN(I,O1) = SHOM *QUAD_STRMWTS(I)
           ENDDO
         ENDDO
       endif

!  Surface Source Term (Polarized calculation for BRDF surfaces)
!   R. Spurr, Bug Fix, 4/9/19. User_Emissivity should not be used when MS-only is true.
!    -- 1/31/21. Version 2.8.3. BRDF arrays are defined locally, drop "M" Fourier index

       IF ( DO_LAMBERTIAN_SURFACE ) THEN
         FACTOR = SURFACE_FACTOR * ALBEDO * SUM(L_IDOWN(1:nstreams,o11)) ; EMISS = one - ALBEDO
         L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) = FACTOR
!  Bug   L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) = L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) + EMISS
         IF (.not.DO_MSMODE_THERMAL) L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) = L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) + EMISS
         if ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:nstreams) = EMISS
       ELSE
         DO UM = 1, N_USER_STREAMS
           do o1 = 1, nstokes
             FACTOR = ZERO
             DO J = 1, NSTREAMS
               S_REFLEC = ZERO
               DO O2 = 1, NSTOKES
                 OM = MUELLER_INDEX(O1,O2)
                 S_REFLEC = S_REFLEC + L_IDOWN(J,O2) * UBRDF_F(OM,UM,J)
               ENDDO
               FACTOR = FACTOR + S_REFLEC
             ENDDO
!  Bug       L_BOA_MSSOURCE(UM,o1) = SURFACE_FACTOR * FACTOR + USER_EMISSIVITY(O1,UM)
             IF (.not.DO_MSMODE_THERMAL) L_BOA_MSSOURCE(UM,o1) = SURFACE_FACTOR * FACTOR + USER_EMISSIVITY(O1,UM)
           enddo
         ENDDO
         if ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:nstreams) = EMISSIVITY(o11,1:nstreams)
       ENDIF

!  End surface inclusion

     ENDIF

!  Upwelling post-processing recursion
!  -----------------------------------

!  Initialize

     L_CUMULSOURCE = zero
     DO UM = 1, N_USER_STREAMS
       L_CUMULSOURCE(UM,1:nstokes) = L_BOA_MSSOURCE(UM,1:nstokes) 
     ENDDO

!  Recursion loop

     NC  = 0;  NUT = 0
     NSTART = NLAYERS ; NUT_PREV = NSTART + 1
     DO UTA = N_USER_LEVELS, 1, -1
       NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
       NUT = NLEVEL + 1
       DO N = NSTART, NUT, -1
         NC = NLAYERS + 1 - N
         DO UM = 1, N_USER_STREAMS
           DO O1 = 1, nstokes
             SHOM_R = ZERO ; SHOM_CR = zero
             IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!  ....Real homogeneous solutions
               DO K = 1, K_REAL(N)
                 H1 = NCON(K,N) * UHOM_UPDN(UM,O1,K,N) * HMULT_2(K,UM,N)
                 H2 = PCON(K,N) * UHOM_UPUP(UM,O1,K,N) * HMULT_1(K,UM,N)
                 SHOM_R = SHOM_R + H1 + H2
               ENDDO
!  ....Complex homogeneous solutions
               DO K = 1, K_COMPLEX(N)
                 K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                 NUXR1 = NCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                 NUXR2 = NCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                 PUXR1 = PCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                 PUXR2 = PCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                 H1 = NUXR1 * HMULT_2(K1,UM,N) - NUXR2 * HMULT_2(K2,UM,N)
                 H2 = PUXR1 * HMULT_1(K1,UM,N) - PUXR2 * HMULT_1(K2,UM,N)
                 SHOM_CR = SHOM_CR + H1 + H2
               ENDDO
             ENDIF
             L_CUMULSOURCE(UM,O1) = SHOM_R + SHOM_CR + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,O1)
           ENDDO
         ENDDO
       ENDDO
       DO UM = 1, N_USER_STREAMS
         SBBWFS_JACOBIANS(UTA,UM,1:nstokes,UPIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,1:nstokes)
       ENDDO
       IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
       NUT_PREV = NUT
     ENDDO

!  end postprocessing upwelling clause

   endif

!  continuation point, removed, 9/17/16.
!544   continue

!  Downwelling Jacobians
!  =====================

!  Skip if not applicable. GOTO 545, Removed 9/17/16
!      IF ( .not. DO_POSTPROCESSING ) GO TO 545
!      IF ( .NOT. DO_UPWELLING )      GO TO 545

   IF ( DO_POSTPROCESSING .and. DO_DNWELLING ) THEN

!  Downwelling post-processing recursion

     L_CUMULSOURCE = zero
     NC  = 0 ; NUT = 0
     NSTART = 1 ;  NUT_PREV = NSTART - 1
     DO UTA = 1, N_USER_LEVELS
       NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
       NUT = NLEVEL
       DO N = NSTART, NUT
         NC = N
         DO UM = 1, N_USER_STREAMS
           DO O1 = 1, nstokes
             SHOM_R = ZERO ; SHOM_CR = zero
             IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!  ....Real homogeneous solutions
               DO K = 1, K_REAL(N)
                 H1 = NCON(K,N) * UHOM_DNDN(UM,O1,K,N) * HMULT_1(K,UM,N)
                 H2 = PCON(K,N) * UHOM_DNUP(UM,O1,K,N) * HMULT_2(K,UM,N)
                 SHOM_R = SHOM_R + H1 + H2
               ENDDO
!  ....Complex homogeneous solutions
               DO K = 1, K_COMPLEX(N)
                 K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                 NUXR1 = NCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
                 NUXR2 = NCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
                 PUXR1 = PCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
                 PUXR2 = PCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
                 H1 = NUXR1 * HMULT_1(K1,UM,N) - NUXR2 * HMULT_1(K2,UM,N)
                 H2 = PUXR1 * HMULT_2(K1,UM,N) - PUXR2 * HMULT_2(K2,UM,N)
                 SHOM_CR = SHOM_CR + H1 + H2
               ENDDO
             ENDIF
             L_CUMULSOURCE(UM,O1) = SHOM_R + SHOM_CR + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,O1)
           ENDDO
         ENDDO
       ENDDO
       DO UM = 1, N_USER_STREAMS
         SBBWFS_JACOBIANS(UTA,UM,1:nstokes,DNIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,1:nstokes)
       ENDDO
       IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
       NUT_PREV = NUT
     ENDDO

!  End postprocessing downwelling

   endif

!  continuation point, removed 9/17/16.
!545   continue

!  Flux Jacobians
!  ==============

!  Upwelling FLux output

   if ( DO_MVOUTPUT .and. DO_UPWELLING ) THEN
      DO UTA = 1, N_USER_LEVELS
         NL = UTAU_LEVEL_MASK_UP(UTA) ; N = NL + 1 ; BBWF_QUAD = zero
         IF ( NL .EQ. NLAYERS  ) THEN
            if ( do_thermal_transonly ) then
               BBWF_QUAD(1:nstreams,o11) = FLUX_MULTIPLIER * L_BOA_THTONLY_SOURCE(1:nstreams)
            else
               DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  DO O1 = 1, nstokes
                     SHOM_R = ZERO ; SHOM_CR = zero
                     DO K = 1, K_REAL(NL)
                        HOM1 = NCON(K,NL) * SOLA_XPOS(I1,O1,K,NL) * T_DELT_EIGEN(K,NL)
                        HOM2 = PCON(K,NL) * SOLB_XNEG(I1,O1,K,NL)
                        SHOM_R = SHOM_R + HOM1 + HOM2
                     ENDDO
                     DO K = 1, K_COMPLEX(NL)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        NXR1   = NCON(K1,NL) * SOLA_XPOS(I1,O1,K1,NL) - NCON(K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
                        NXR2   = NCON(K1,NL) * SOLA_XPOS(I1,O1,K2,NL) + NCON(K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
                        HOM2CR = PCON(K1,NL) * SOLB_XNEG(I1,O1,K1,NL) - PCON(K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
                        HOM1CR = NXR1 * T_DELT_EIGEN(K1,NL) - NXR2 * T_DELT_EIGEN(K2,NL)
                        SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                     ENDDO
                     BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                  ENDDO
               ENDDO
            endif
         ELSE
            if ( do_thermal_transonly ) then
               DO I = 1, NSTREAMS
                  SHOM = L_BOA_THTONLY_SOURCE(I)
                  DO LAY = NLAYERS, N, -1
                    SHOM = SHOM * T_DELT_DISORDS(I,LAY)
                  ENDDO
                  BBWF_QUAD(I,O11) = FLUX_MULTIPLIER * SHOM
               ENDDO
            else
               DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  DO O1 = 1, nstokes
                     SHOM_R = ZERO ; SHOM_CR = zero
                     DO K = 1, K_REAL(N)
                        HOM1 = NCON(K,N) * SOLA_XPOS(I1,O1,K,N)
                        HOM2 = PCON(K,N) * SOLB_XNEG(I1,O1,K,N) * T_DELT_EIGEN(K,N)
                        SHOM_R = SHOM_R + HOM1 + HOM2
                     ENDDO
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        HOM1CR = NCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                        PXR1   = PCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - PCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                        PXR2   = PCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
                        HOM2CR = PXR1 * T_DELT_EIGEN(K1,N) - PXR2 * T_DELT_EIGEN(K2,N)
                        SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                     ENDDO
                     BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                  ENDDO
               ENDDO
            endif
         ENDIF

!  Assign fluxes

         do o1 = 1, nstokes
            SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_WEIGHTS(1:NSTREAMS))
            SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS))
            SBBWFS_FLUXES(UTA,1,O1,UPIDX) = HALF * SM
            SBBWFS_FLUXES(UTA,2,O1,UPIDX) = PI2  * SF
         enddo
      ENDDO
   ENDIF

!  Downwelling FLux output. Nothing for the transmittance-only case.
!  -----------------------------------------------------------------

   if ( DO_MVOUTPUT .and. DO_DNWELLING ) THEN
      DO UTA = 1, N_USER_LEVELS
         NL = UTAU_LEVEL_MASK_DN(UTA) ; N = NL ;  BBWF_QUAD = ZERO
         IF ( NL .NE. 0 .and. .not. do_thermal_transonly  ) THEN
            DO I = 1, NSTREAMS
               DO O1 = 1, nstokes
                  SHOM_R = ZERO ; SHOM_CR = zero
                  DO K = 1, K_REAL(N)
                     HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N) * T_DELT_EIGEN(K,N)
                     HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)
                     SHOM_R = SHOM_R + HOM1 + HOM2
                  ENDDO
                  DO K = 1, K_COMPLEX(N)
                     K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                     NXR1   = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                     NXR2   = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                     HOM2CR = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                     HOM1CR = NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
                     SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                  ENDDO
                  BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
               ENDDO
            ENDDO
         ENDIF

!  Assign fluxes

         do o1 = 1, nstokes
            SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_WEIGHTS(1:NSTREAMS))
            SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS))
            SBBWFS_FLUXES(UTA,1,O1,DNIDX) = HALF * SM
            SBBWFS_FLUXES(UTA,2,O1,DNIDX) = PI2  * SF
         enddo

      ENDDO
   ENDIF

!  FINISH

   return
end subroutine vlidort_lbbf_jacobians_whole

!  

subroutine vlidort_lbbf_jacobians_wpartials &
      ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
        DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
        DO_MSMODE_THERMAL, DO_POSTPROCESSING, DO_MVOUTPUT,         & ! input
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,                 & ! input
        NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
        NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & ! Input
        NTOTAL, N_SUPDIAG, N_SUBDIAG,  MUELLER_INDEX,              & ! input
        N_PARTLAYERS, PARTLAYERS_LAYERIDX,                         & ! Input
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                    & ! Input
        USER_STREAMS, LAYERMASK_UP, LAYERMASK_DN,                  & ! Input
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
        SURFACE_FACTOR, ALBEDO, BRDF_F, UBRDF_F,                   & ! input
        EMISSIVITY, USER_EMISSIVITY,  FLUX_MULTIPLIER,             & ! input
        DELTAUS, PARTAUS, OMEGA_GREEK,                             & ! Input
        T_DELT_DISORDS, T_UTDN_DISORDS, T_UTUP_DISORDS,            & ! input
        PI_XQP, PI_XQM_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,                   & ! input
        T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                  & ! input
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                  & ! Input
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
        UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,        & ! Input
        UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
        ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
        SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
        STATUS, MESSAGE, TRACE )                                     ! Output

!  Module file of dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, MAX_USER_STREAMS,     &
                                 MAX_USER_LEVELS, MAXMOMENTS, MAXTOTAL, MAXBANDTOTAL, MAXEVALUES,        &
                                 MAXSTRMSTKS, MAXSTREAMS_2, MAXSTRMSTKS_2, MAXSTOKES_SQ, MAX_DIRECTIONS, &
                                 UPIDX, DNIDX, ZERO, HALF, ONE, TWO, PI2, PI4, VLIDORT_SERIOUS, VLIDORT_SUCCESS

!  VLIDORT module dependencies

   USE lapack_tools_m, only : DGETRF, DGETRS, DGBTRS

      implicit none

!  Subroutine input arguments
!  --------------------------

!  Master control

      LOGICAL, INTENT(IN)  :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  local control flags

      LOGICAL, INTENT(IN)  :: DO_THERMAL_TRANSONLY
      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

      LOGICAL, INTENT(IN)  :: DO_MSMODE_THERMAL
      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING
      LOGICAL, INTENT(IN)  :: DO_MVOUTPUT

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_LAMBERTIAN_SURFACE

!  Numbers

      INTEGER, INTENT(IN)  :: NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, NMOMENTS
      INTEGER, INTENT(IN)  :: NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NSTREAMS_2
      INTEGER, INTENT (IN) :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  multiplier

      DOUBLE PRECISION, INTENT(IN)  :: FLUX_MULTIPLIER

!  Partials control

      INTEGER  , intent(in)  :: N_PARTLAYERS
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)

!  output   control

      INTEGER, INTENT (IN)   :: N_USER_LEVELS
      INTEGER, INTENT (IN)   :: UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN)   :: UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  User polar directions, postprocessnd control

      DOUBLE PRECISION, INTENT(IN)  :: USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL         , INTENT(IN)  :: LAYERMASK_UP ( MAXLAYERS )
      LOGICAL         , INTENT(IN)  :: LAYERMASK_DN ( MAXLAYERS )

!  Quadrature values

      DOUBLE PRECISION, intent(in)   :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, intent(in)   :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, intent(in)   :: QUAD_STRMWTS ( MAXSTREAMS )

!  Optical properties

      DOUBLE PRECISION, INTENT(IN)   :: DELTAUS     ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN)   :: PARTAUS     ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT(IN)   :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Polynomials

      DOUBLE PRECISION, INTENT(IN)   :: PI_XQP     ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT(IN)   :: PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT(IN)   :: PI_XUP     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT(IN)   :: PI_XUM     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Discrete ordinate solutions
!  ---------------------------

!  Direct solutions, stream transmittances

      DOUBLE PRECISION, intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: T_UTDN_DISORDS(MAXSTREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, intent(in)  :: T_UTUP_DISORDS(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigensolutions, eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      INTEGER         , INTENT (IN) :: K_REAL    ( MAXLAYERS )
      INTEGER         , INTENT (IN) :: K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Matrices

      DOUBLE PRECISION, INTENT (IN) :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )

!  BVProblem stuff
!  ---------------

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER         , INTENT (IN) :: IPIVOT   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2    ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER         , INTENT (IN) :: SIPIVOT  ( MAXSTRMSTKS_2 )

!  Surface stuff
!  -------------

!  1/31/21. Version 2.8.3. BRDF arrays are defined locally, drop "MAXMOMENTS" dimension

      DOUBLE PRECISION, INTENT(IN)   :: SURFACE_FACTOR, ALBEDO
      DOUBLE PRECISION, intent(in)   :: BRDF_F  ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN)   :: UBRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION, intent(in)   :: EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, intent(in)   :: USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  User-angle (post-processed) solution variables
!  ----------------------------------------------

!  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION, INTENT(IN)  :: T_DELT_USERM  (MAXLAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION, intent(in)  :: T_UTUP_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS)
      DOUBLE PRECISION, intent(in)  :: T_UTDN_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS)

!  User solutions defined at user-defined stream angles

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  solution multipliers 

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: UT_HMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, intent(in)  :: UT_HMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, intent(in)  :: UT_HMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION, intent(in)  :: UT_HMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Outputs
!  -------

!  Postprocessed Jacobians.
!  Outputs are all Pre-zeroed in the calling Masters

      DOUBLE PRECISION, INTENT(INOUT) :: ABBWFS_JACOBIANS &
               ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)
      DOUBLE PRECISION, INTENT(INOUT) :: SBBWFS_JACOBIANS &
               ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXSTOKES, MAX_DIRECTIONS)

!  Flux Jacobians.
!  Outputs are all Pre-zeroed in the calling Masters

      DOUBLE PRECISION, INTENT(INOUT) :: ABBWFS_FLUXES ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(INOUT) :: SBBWFS_FLUXES ( MAX_USER_LEVELS, 2, MAXSTOKES, MAX_DIRECTIONS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  LOCAL THERMAL-BBF JACOBIAN ARRAYS
!  =================================

!  Weighting function column matrices

      DOUBLE PRECISION  :: COL2_BWF  ( MAXTOTAL )
      DOUBLE PRECISION  :: SCOL2_BWF ( MAXSTRMSTKS_2 )

!  Linearized Solution constants of integration

      DOUBLE PRECISION  :: NCON(MAXSTRMSTKS,MAXLAYERS)
      DOUBLE PRECISION  :: PCON(MAXSTRMSTKS,MAXLAYERS)

!  Linearized Thermal solutions at the Upper/Lower boundary

      DOUBLE PRECISION  :: L_T_WUPPER_Gp1(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION  :: L_T_WUPPER_Gp2(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION  :: L_T_WLOWER_Gp1(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION  :: L_T_WLOWER_Gp2(MAXSTREAMS_2,MAXLAYERS)

!  Linearized Thermal layer source terms

      DOUBLE PRECISION  :: L_LAYER_TSUP_UP_Gp1(MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_UP_Gp2(MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_DN_Gp1(MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_DN_Gp2(MAX_USER_STREAMS,MAXLAYERS)

!  Linearization of Direct solutions

      DOUBLE PRECISION  :: L_T_DIRECT_UP_Gp1 (MAX_USER_STREAMS,MAXLAYERS )
      DOUBLE PRECISION  :: L_T_DIRECT_UP_Gp2 (MAX_USER_STREAMS,MAXLAYERS )
      DOUBLE PRECISION  :: L_T_DIRECT_DN_Gp1 (MAX_USER_STREAMS,MAXLAYERS )
      DOUBLE PRECISION  :: L_T_DIRECT_DN_Gp2 (MAX_USER_STREAMS,MAXLAYERS )

!  Help arrays for the discrete ordinate field

      INTEGER          :: TPIVOT ( MAXSTREAMS )
      DOUBLE PRECISION :: TMAT ( MAXSTREAMS, MAXSTREAMS ), HVEC ( MAXSTREAMS )
      DOUBLE PRECISION :: TVEC1_Gp1 ( MAXSTREAMS_2, MAXLAYERS ), TVEC1_Gp2 ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: TVEC2_Gp1 ( MAXSTREAMS_2, MAXLAYERS ), TVEC2_Gp2 ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: HVEC1_Gp1 ( MAXSTREAMS ), HVEC1_Gp2 ( MAXSTREAMS )
      DOUBLE PRECISION :: HVEC2_Gp1 ( MAXSTREAMS ), HVEC2_Gp2 ( MAXSTREAMS )
      DOUBLE PRECISION :: JVEC1_Gp1 ( MAXSTREAMS ), JVEC1_Gp2 ( MAXSTREAMS )

!  Help arrays for the post-processed field

      DOUBLE PRECISION :: T_HELP1_Gp1 ( 0:MAXMOMENTS ), T_HELP1_Gp2 ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: T_HELP2_Gp1 ( 0:MAXMOMENTS ), T_HELP2_Gp2 ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: U_TPOS1_Gp1 ( MAX_USER_STREAMS, MAXLAYERS ), U_TPOS1_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TPOS2_Gp1 ( MAX_USER_STREAMS, MAXLAYERS ), U_TPOS2_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG1_Gp1 ( MAX_USER_STREAMS, MAXLAYERS ), U_TNEG1_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG2_Gp1 ( MAX_USER_STREAMS, MAXLAYERS ), U_TNEG2_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )

!  Linearized Partial-layer quantities
!  -----------------------------------

!  Linearization of Direct solutions

      DOUBLE PRECISION  :: L_T_UT_DIRECT_UP_Gp1 ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION  :: L_T_UT_DIRECT_UP_Gp2 ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION  :: L_T_UT_DIRECT_DN_Gp1 ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION  :: L_T_UT_DIRECT_DN_Gp2 ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearization of Partial-layer solutions

      DOUBLE PRECISION  :: L_ut_t_partic_Gp1(MAXSTREAMS_2,MAX_PARTLAYERS)
      DOUBLE PRECISION  :: L_ut_t_partic_Gp2(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Partial layer sources

      DOUBLE PRECISION  :: L_LAYER_TSUP_UTUP_Gp1(MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_UTUP_Gp2(MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_UTDN_Gp1(MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION  :: L_LAYER_TSUP_UTDN_Gp2(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Reflectance integrands, BOA source terms
!    Vector change 4/27/16. Give Reflec a Stokes dimension

      DOUBLE PRECISION :: L_IDOWN(MAXSTREAMS,MAXSTOKES), DOWN(MAXSTREAMS)
      DOUBLE PRECISION :: REFLEC(MAXSTREAMS,MAXSTOKES), BBWF_QUAD(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION :: L_BOA_MSSOURCE ( MAX_USER_STREAMS, MAXSTOKES ), L_BOA_THTONLY_SOURCE (MAXSTREAMS)

!  Local layer and cumulative source terms

      DOUBLE PRECISION :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: L_CUMULSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

!  help variables

      LOGICAL   :: DO_QTHTONLY
      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC, NL, K, KO1, K0, K1, K2, LAY, O1, O2, O11, OM, JB, AA, AA1, L
      INTEGER   :: UTA, UT, UM, N1, M, I, I1, J, J1, INFO, IR, IROW, IROW1, IROW_S, IROW1_S, IC, IC1, IC2, CMP, CM, C0

      DOUBLE PRECISION :: H1, H2, HOM1, HOM2, HOM1CR, HOM2CR
      DOUBLE PRECISION :: NXR1, NXR2, PXR1, PXR2, NUXR1, NUXR2, PUXR1, PUXR2
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SHOM, POS1, POS2, NEG1, NEG2
      DOUBLE PRECISION :: Uxup, Uxdn, zup, zdn, zd1, zu1, zd1_ok, quad, cmdel
      DOUBLE PRECISION :: NXR_CR, NXR_CI, PXR_CR, PXR_CI                           ! Added 4/21/16

      DOUBLE PRECISION :: SM, TM, SF, COSMUM, SUM, FAC, FACTOR, omega1_odelt, omega
      DOUBLE PRECISION :: U1, U2, D1, D2, SPAR, SPAR1, SPAR2, EMISS, L_THELP, HOMIGR
      DOUBLE PRECISION :: Udel, DelUdel, omega1, A5_xqp, A5_xqm, z1, z, zd, z1_ok, S_REFLEC, FINAL_SOURCE
      DOUBLE PRECISION :: PN1_Gp1, PN1_Gp2, PN2_Gp1, PN2_Gp2, Sum1_Gp1, Sum1_Gp2, Sum2_Gp1, Sum2_Gp2
      CHARACTER*3 :: CI, CN, C3

      DOUBLE PRECISION :: Group1(maxlayers,2), Group2(maxlayers,2)

!  Initial section
!  ---------------

!  Exception handling

   STATUS  = VLIDORT_SUCCESS
   MESSAGE = ' '
   TRACE   = ' '

!  Proxies

   m  = 0

!  Use this variable when only the nstokes = 1 component is required

   o11 = 1

!  Initial modulus = 4.pi if solar sources are included

   fac = one
   if ( do_solar_sources ) fac = pi4

!  Local flag

   DO_QTHTONLY = do_MVOUTPUT .and. DO_THERMAL_TRANSONLY

!  Control to SURFACE LBBF Section

   if ( .not. DO_ATMOS_LBBF ) go to 55

!  Group 1/2 Derivatives of TCOM1, all layers
!   Group 1, w.r.t to the upper boundary BBF
!   Group 2, w.r.t to the lower boundary BBF
!     Assumes only 2 coefficients, simplified form

   do n = 1, nlayers
      omega = omega_greek(0,n,1,1)
      omega1 = two * ( one - omega ) ; if ( do_thermal_transonly ) omega1 = two
      omega1_odelt = omega1 / deltaus(n)
      Group1(n,1)  =   omega1
      Group1(n,2)  = - omega1_odelt
      Group2(n,1)  =   zero
      Group2(n,2)  = + omega1_odelt
   enddo

!  Linearization of Direct Term
!  ----------------------------

!  Zero the terms  first, then skip if in MSMODE-only of Luxes-only
!  Linearization of Direct Term, Zero terms in MSMODE-only

   L_t_direct_up_Gp1 = zero ; L_t_direct_up_Gp2 = zero
   L_t_direct_dn_Gp1 = zero ; L_t_direct_dn_Gp2 = zero
   L_t_ut_direct_up_Gp1 = zero ; L_t_ut_direct_up_Gp2 = zero
   L_t_ut_direct_dn_Gp1 = zero ; L_t_ut_direct_dn_Gp2 = zero

   IF ( DO_POSTPROCESSING .and. DO_MSMODE_THERMAL ) go to 68
   IF ( .not. DO_POSTPROCESSING                   ) go to 68

!  Upwelling Direct solution source terms
!  --------------------------------------

   IF ( do_upwelling ) THEN
      DO um = 1, n_user_streams
         cosmum = user_streams(um)

!  Whole layer....

         do n = 1, nlayers
            if ( layermask_up(n) ) then
               Udel = t_delt_userm(n,um)
               u1 = one - Udel ; u2 = cosmum - Udel * ( cosmum + deltaus(n) )
               L_t_direct_up_Gp1(um,n) = half * ( u1 * Group1(n,1) + u2 * Group1(n,2) )
               L_t_direct_up_Gp2(um,n) = half * ( u1 * Group2(n,1) + u2 * Group2(n,2) )
            endif
         enddo

!  Partial layer...

         IF ( n_partlayers.ne.0 ) THEN
            DO ut = 1, n_PARTLAYERS
               Uxup = t_utup_userm(ut,um)
               n  = partlayers_layeridx(ut) ; cmdel = cosmum + deltaus(n)
               u1 = one - Uxup ; u2 = partaus(ut) + cosmum - Uxup * cmdel
               L_t_ut_direct_up_Gp1(um,ut) = half * ( u1 * Group1(n,1) + u2 * Group1(n,2) )
               L_t_ut_direct_up_Gp2(um,ut) = half * ( u1 * Group2(n,1) + u2 * Group2(n,2) )
            enddo
         endif

      enddo
   endif

!  Downwelling Direct solution source terms
!  ----------------------------------------

   IF ( do_dnwelling ) THEN
      DO um = 1, n_user_streams
         cosmum = user_streams(um)

!  Whole layer

         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               Udel = t_delt_userm(n,um)
               d1 = one - Udel ; d2 = deltaus(n) - cosmum * d1
               L_t_direct_dn_Gp1(um,n) = half * ( d1 * Group1(n,1) + d2 * Group1(n,2) )
               L_t_direct_dn_Gp2(um,n) = half * ( d1 * Group2(n,1) + d2 * Group2(n,2) )
            endif
         enddo

!  Partial layer...

         IF ( n_partlayers.ne.0 ) THEN
            DO ut = 1, n_PARTLAYERS
               Uxdn = t_utdn_userm(ut,um)
               n  = partlayers_layeridx(ut)
               d1 = one - Uxdn ; d2 = partaus(ut) - cosmum * d1
               L_t_ut_direct_dn_Gp1(um,ut) = half * ( d1 * Group1(n,1) + d2 * Group1(n,2) )
               L_t_ut_direct_dn_Gp2(um,ut) = half * ( d1 * Group2(n,1) + d2 * Group2(n,2) )
            enddo
         endif

!  Finish

      enddo
   endif

!  Continuation point when Linearization of direct term not required

68 continue

!  Thermal Transmittance only, quadrature solutions
!  ================================================

   if ( do_thermal_transonly ) then

!   Whole layer

      DO n = 1, nlayers
         DO aa = 1, nstreams
            aa1 = aa + nstreams
            Z = t_delt_disords(aa,n) ; zd = Z * deltaus(n) ; z1 = one - Z ; z1_ok = z1 * quad_streams(aa)
            d2 =  ( deltaus(n) - z1_ok ) ; d1 = z1
            u2 =   ( z1_ok - zd )        ; u1 = z1
!   Bugs, 4/21/16
!            L_t_wupper_Gp1(aa1,n)  = u2 * Group1(n,2) + u1 * Group1(n,1)
!            L_t_wupper_Gp2(aa1,n)  = u2 * Group2(n,2) + u1 * Group2(n,1)
!            L_t_wlower_Gp1(aa,n)   = d2 * Group1(n,2) + d1 * Group1(n,1)
!            L_t_wlower_Gp2(aa,n)   = d2 * Group2(n,2) + d1 * Group2(n,1)
            L_t_wupper_Gp1(aa1,n)  = half * ( u2 * Group1(n,2) + u1 * Group1(n,1) ) * quad_streams(aa)
            L_t_wupper_Gp2(aa1,n)  = half * ( u2 * Group2(n,2) + u1 * Group2(n,1) ) * quad_streams(aa)
            L_t_wlower_Gp1(aa,n)   = half * ( d2 * Group1(n,2) + d1 * Group1(n,1) ) * quad_streams(aa)
            L_t_wlower_Gp2(aa,n)   = half * ( d2 * Group2(n,2) + d1 * Group2(n,1) ) * quad_streams(aa)
         END DO
      END DO

!  Partial layer

      if ( DO_MVOUTPUT .and. n_PARTLAYERS .gt. 0 ) then
         DO ut = 1, n_PARTLAYERS
            n  = partlayers_layeridx(ut)
            do aa = 1, nstreams
               aa1 = aa + nstreams
               Zup  = t_utup_disords(aa,ut) ; Zdn  = t_utdn_disords(aa,ut)
               Zu1 = one - Zup ; zd1 = one - Zdn ; zd1_ok =  zd1 * quad_streams(aa)
               d1 = zd1 ; d2 =  ( partaus(ut) - zd1_ok )
               u1 = zu1 ; u2 = ( zd1_ok + partaus(ut) - deltaus(n) * zup)
!   Bugs, 4/21/16
!               L_ut_t_partic_Gp1(aa1,ut)  = u2 * Group1(n,2) + u1 * Group1(n,1)
!               L_ut_t_partic_Gp2(aa1,ut)  = u2 * Group2(n,2) + u1 * Group2(n,1)
!               L_ut_t_partic_Gp1(aa,ut)   = d2 * Group1(n,2) + d1 * Group1(n,1)
!               L_ut_t_partic_Gp2(aa,ut)   = d2 * Group2(n,2) + d1 * Group2(n,1)
               L_ut_t_partic_Gp1(aa1,ut)  = half * ( u2 * Group1(n,2) + u1 * Group1(n,1) ) * quad_streams(aa)
               L_ut_t_partic_Gp2(aa1,ut)  = half * ( u2 * Group2(n,2) + u1 * Group2(n,1) ) * quad_streams(aa)
               L_ut_t_partic_Gp1(aa,ut)   = half * ( d2 * Group1(n,2) + d1 * Group1(n,1) ) * quad_streams(aa)
               L_ut_t_partic_Gp2(aa,ut)   = half * ( d2 * Group2(n,2) + d1 * Group2(n,1) ) * quad_streams(aa)
            enddo
         enddo
      endif

      GO TO 74
   endif

!  Start Layer loop for solution derivatives
!  =========================================

!  Comment, 26 March 2014
!    We will have to go back to the originals here because not saved

   do n = 1, nlayers

!  First SOLUTION MATRIX and LU-decomposition  (Same as in Thermal CLSolution)

      DO I = 1, NSTREAMS
         TMAT(I,1:nstreams) = SAB(I,1:nstreams,O11,O11,N) * QUAD_STREAMS(I)
      ENDDO
      CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)
      IF ( INFO .GT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(C3, '(I3)' ) N
         MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
         TRACE   = 'DGETRF CALL FOR H-vector, LAYER '//C3//' LBBF_Jacobians_whole, Call 1'
         STATUS  = VLIDORT_SERIOUS ; return
      ENDIF

!  Linearized H VECTOR_1 Group-1 :  SOLUTION BY BACK-SUBSTITUTION
!  Linearized H VECTOR_1 Group-2 :  Zero

      HVEC = one
      CALL DGETRS  ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,HVEC,MAXSTREAMS,INFO)
      IF ( INFO .GT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(C3, '(I3)' ) N
         MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
         TRACE   = 'DGETRS CALL FOR H-vector, LAYER '//C3//' LBBF_Jacobians_whole, Call 2'
         STATUS  = VLIDORT_SERIOUS ; return
      ENDIF
      HVEC1_Gp1(1:nstreams) = Group1(n,1) * HVEC(1:nstreams)
      HVEC1_Gp2(1:nstreams) = zero

!  Linearized H VECTOR_2 Group-1 :  SOLUTION BY BACK-SUBSTITUTION
!  Linearized H VECTOR_2 Group-2 :  SOLUTION BY BACK-SUBSTITUTION

      HVEC2_Gp1(1:nstreams) = Group1(n,2) * HVEC(1:nstreams)
      HVEC2_Gp2(1:nstreams) = Group2(n,2) * HVEC(1:nstreams)

!  Second SOLUTION MATRIX and LU-decomposition  (Same as in Thermal CLSolution)

      DO I = 1, NSTREAMS
         TMAT(I,1:nstreams) = -DAB(I,1:nstreams,O11,O11,N)
      ENDDO
      CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)
      IF ( INFO .GT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(C3, '(I3)' ) N
         MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
         TRACE   = 'DGETRF CALL FOR J-vector, LAYER '//C3//' LBBF_Jacobians_whole, Call 3'
         STATUS  = VLIDORT_SERIOUS ; return
      ENDIF

!  Linearized J VECTOR_2 Groups 1-2 :  SOLUTION BY BACK-SUBSTITUTION

      CALL DGETRS  ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,HVEC,MAXSTREAMS,INFO)
      IF ( INFO .GT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(C3, '(I3)' ) N
         MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
         TRACE   = 'DGETRS CALL FOR H-vector, LAYER '//C3//' LBBF_Jacobians_whole, Call 4'
         STATUS  = VLIDORT_SERIOUS ; return
      ENDIF
      JVEC1_Gp1(1:nstreams) = Group1(n,2) * HVEC(1:nstreams)
      JVEC1_Gp2(1:nstreams) = Group2(n,2) * HVEC(1:nstreams)

!  Set solution

      DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         TVEC1_Gp1(I,N)   = HALF * (HVEC1_Gp1(I) + JVEC1_Gp1(I))
         TVEC1_Gp2(I,N)   = HALF * (HVEC1_Gp2(I) + JVEC1_Gp2(I))
         TVEC1_Gp1(I1,N)  = HALF * (HVEC1_Gp1(I) - JVEC1_Gp1(I))
         TVEC1_Gp2(I1,N)  = HALF * (HVEC1_Gp2(I) - JVEC1_Gp2(I))
      ENDDO
!  This code is Wrong 4/15/16
!      DO I = 1, NSTREAMS_2
!         TVEC2_Gp1(I,N)  = HALF * HVEC2_Gp1(I)
!         TVEC2_Gp2(I,N)  = HALF * HVEC2_Gp2(I)
!      ENDDO
      DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         TVEC2_Gp1(I,N)   = HALF * HVEC2_Gp1(I)
         TVEC2_Gp2(I,N)   = HALF * HVEC2_Gp2(I)
         TVEC2_Gp1(I1,N)  = HALF * HVEC2_Gp1(I)
         TVEC2_Gp2(I1,N)  = HALF * HVEC2_Gp2(I)
      ENDDO

      DO I = 1, NSTREAMS_2
         L_T_WUPPER_Gp1(I,N) = TVEC1_Gp1(I,N)
         L_T_WUPPER_Gp2(I,N) = TVEC1_Gp2(I,N)
!         L_T_WLOWER_Gp1(I,N) = TVEC1_Gp1(I,N) + TVEC2_Gp2(I,N) * DELTAUS(N)   !   Bug 4/20/16
         L_T_WLOWER_Gp1(I,N) = TVEC1_Gp1(I,N) + TVEC2_Gp1(I,N) * DELTAUS(N)
         L_T_WLOWER_Gp2(I,N) = TVEC1_Gp2(I,N) + TVEC2_Gp2(I,N) * DELTAUS(N)
      ENDDO

!  USER SOLUTIONS
!    4/15/16. Bad bug in the original, see commented out code

      DO L = M, NMOMENTS
         SUM1_Gp1 = ZERO ; SUM2_Gp1 = ZERO
         SUM1_Gp2 = ZERO ; SUM2_Gp2 = ZERO
         HOMIGR = half * OMEGA_GREEK(L,N,O11,O11)
         DO  J = 1, NSTREAMS
            J1 = J + NSTREAMS
            A5_xqp = QUAD_WEIGHTS(J) * PI_XQP    (L,J,O11,O11) * HOMIGR
            A5_xqm = QUAD_WEIGHTS(J) * PI_XQM_PRE(L,J,O11,O11) * HOMIGR
!            PN1_Gp1 = TVEC1_Gp1(J1,N) * A5_xqp + TVEC1_Gp1(J,N) + A5_xqm ;  SUM1_Gp1 = SUM1_Gp1 + PN1_Gp1
!            PN1_Gp2 = TVEC1_Gp2(J1,N) * A5_xqp + TVEC1_Gp2(J,N) + A5_xqm ;  SUM1_Gp2 = SUM1_Gp2 + PN1_Gp2
!            PN2_Gp1 = TVEC2_Gp1(J1,N) * A5_xqp + TVEC2_Gp1(J,N) + A5_xqm ;  SUM2_Gp1 = SUM2_Gp1 + PN2_Gp1
!            PN2_Gp2 = TVEC2_Gp2(J1,N) * A5_xqp + TVEC2_Gp2(J,N) + A5_xqm ;  SUM2_Gp2 = SUM2_Gp2 + PN2_Gp2
            PN1_Gp1 = TVEC1_Gp1(J1,N) * A5_xqp + TVEC1_Gp1(J,N) * A5_xqm ;  SUM1_Gp1 = SUM1_Gp1 + PN1_Gp1
            PN1_Gp2 = TVEC1_Gp2(J1,N) * A5_xqp + TVEC1_Gp2(J,N) * A5_xqm ;  SUM1_Gp2 = SUM1_Gp2 + PN1_Gp2
            PN2_Gp1 = TVEC2_Gp1(J1,N) * A5_xqp + TVEC2_Gp1(J,N) * A5_xqm ;  SUM2_Gp1 = SUM2_Gp1 + PN2_Gp1
            PN2_Gp2 = TVEC2_Gp2(J1,N) * A5_xqp + TVEC2_Gp2(J,N) * A5_xqm ;  SUM2_Gp2 = SUM2_Gp2 + PN2_Gp2
         ENDDO
         T_HELP1_Gp1(L) = SUM1_Gp1 ; T_HELP1_Gp2(L) = SUM1_Gp2
         T_HELP2_Gp1(L) = SUM2_Gp1 ; T_HELP2_Gp2(L) = SUM2_Gp2 
      ENDDO

!  UPWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

      IF ( DO_UPWELLING.AND.LAYERMASK_UP(N) ) THEN
         DO UM = 1, N_USER_STREAMS
            POS1 = dot_product(T_HELP1_Gp1(M:NMOMENTS),PI_XUP(M:NMOMENTS,UM,O11,O11) )
            POS2 = dot_product(T_HELP2_Gp1(M:NMOMENTS),PI_XUP(M:NMOMENTS,UM,O11,O11) )
            U_TPOS1_Gp1(UM,N) = POS1 ; U_TPOS2_Gp1(UM,N) = POS2
            POS1 = dot_product(T_HELP1_Gp2(M:NMOMENTS),PI_XUP(M:NMOMENTS,UM,O11,O11) )
            POS2 = dot_product(T_HELP2_Gp2(M:NMOMENTS),PI_XUP(M:NMOMENTS,UM,O11,O11) )
            U_TPOS1_Gp2(UM,N) = POS1 ; U_TPOS2_Gp2(UM,N) = POS2
         ENDDO
      ENDIF

!  DNWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

      IF ( DO_DNWELLING.AND.LAYERMASK_DN(N) ) THEN
         DO UM = 1, N_USER_STREAMS
            NEG1 = dot_product(T_HELP1_Gp1(M:NMOMENTS),PI_XUM(M:NMOMENTS,UM,O11,O11) )
            NEG2 = dot_product(T_HELP2_Gp1(M:NMOMENTS),PI_XUM(M:NMOMENTS,UM,O11,O11) )
            U_TNEG1_Gp1(UM,N) = NEG1 ; U_TNEG2_Gp1(UM,N) = NEG2
            NEG1 = dot_product(T_HELP1_Gp2(M:NMOMENTS),PI_XUM(M:NMOMENTS,UM,O11,O11) )
            NEG2 = dot_product(T_HELP2_Gp2(M:NMOMENTS),PI_XUM(M:NMOMENTS,UM,O11,O11) )
            U_TNEG1_Gp2(UM,N) = NEG1 ; U_TNEG2_Gp2(UM,N) = NEG2
         ENDDO
      ENDIF

!  End layer loop

   ENDDO

!  Linearization of Partial Green's function
!  -----------------------------------------

   if ( DO_MVOUTPUT .and. n_PARTLAYERS .gt. 0 ) then

!  start loop over offgrid optical depths

      DO ut = 1, n_PARTLAYERS
         n  = partlayers_layeridx(ut)

!  Simple rule here.

         DO i = 1, nstreams_2
            L_ut_t_partic_Gp1(i,ut) = TVEC1_Gp1(I,N) + TVEC2_Gp1(I,N) * PARTAUS(UT)
            L_ut_t_partic_Gp2(i,ut) = TVEC1_Gp2(I,N) + TVEC2_Gp2(I,N) * PARTAUS(UT)
         enddo

!  Following placeholder code commented out 4/21/16
!  upwelling solutions
!         IF ( do_upwelling ) THEN
!            DO i = 1, nstreams
!               i1 = i + nstreams
!               spar1 = zero ;  spar2 = zero
!  PLACEHOLDER  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    HERE
!               L_ut_t_partic_Gp1(i1,ut) = spar1
!               L_ut_t_partic_Gp2(i1,ut) = spar2
!            END DO
!         END IF
!  Downwelling solutions
!         IF ( do_dnwelling ) THEN
!            DO i = 1, nstreams
!               i1 = i + nstreams
!               spar1 = zero ;  spar2 = zero
!  PLACEHOLDER  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    HERE
!               L_ut_t_partic_Gp1(i,ut) = spar1
!               L_ut_t_partic_Gp2(i,ut) = spar2
!            END DO
!         END IF

! finish off-grid solutions

      END DO
   END IF

!  Layer source term derivatives
!  =============================

!  Continuation point for thermal tranmsittance only solutions

74 continue

!  Initialize completely, skip if no post processing

   L_LAYER_TSUP_UP_Gp1   = zero ; L_LAYER_TSUP_DN_Gp1   = zero
   L_LAYER_TSUP_UP_Gp2   = zero ; L_LAYER_TSUP_DN_Gp2   = zero
   L_LAYER_TSUP_UTUP_Gp1 = zero ; L_LAYER_TSUP_UTDN_Gp1 = zero
   L_LAYER_TSUP_UTUP_Gp2 = zero ; L_LAYER_TSUP_UTDN_Gp2 = zero

   if ( .not. do_POSTPROCESSING ) go to 69

!  Initialize with Direct term (which may be zero...)
!  --------------------------------------------------

   DO um = 1, n_user_streams
      do n = 1, nlayers
         if ( do_upwelling.and.layermask_up(n) ) then
            L_layer_tsup_up_Gp1(um,n) = fac * L_t_direct_up_Gp1(um,n)
            L_layer_tsup_up_Gp2(um,n) = fac * L_t_direct_up_Gp2(um,n)
         endif
         if ( do_dnwelling .and. layermask_dn(n) ) then
            L_layer_tsup_dn_Gp1(um,n) = fac * L_t_direct_dn_Gp1(um,n)
            L_layer_tsup_dn_Gp2(um,n) = fac * L_t_direct_dn_Gp2(um,n)
         endif
      enddo
      if ( do_upwelling .and. n_partlayers .ne. 0 ) then
         do ut = 1, n_PARTLAYERS
            L_layer_tsup_utup_Gp1(um,ut) = fac * L_t_ut_direct_up_Gp1(um,ut)
            L_layer_tsup_utup_Gp2(um,ut) = fac * L_t_ut_direct_up_Gp2(um,ut)
         enddo
      endif
      if ( do_dnwelling .and. n_partlayers .ne. 0 ) then
         do ut = 1, n_PARTLAYERS
            L_layer_tsup_utdn_Gp1(um,ut) = fac * L_t_ut_direct_dn_Gp1(um,ut)
            L_layer_tsup_utdn_Gp2(um,ut) = fac * L_t_ut_direct_dn_Gp2(um,ut)
         enddo
      endif
   enddo

!  done if thermal Transmittance only

   if ( do_thermal_transonly ) go to 69

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!  --------------------------------------------------

   DO UM = 1, N_USER_STREAMS

      COSMUM = USER_STREAMS(UM)

!  Whole-layer loop ------------>

      do n = 1, nlayers

         Udel = t_delt_userm(n,um) ;  delUdel = deltaus(n) * Udel
         u1 = one - udel ; u2 = cosmum * u1 - delUdel

!  Upwelling

         if ( do_upwelling .and. layermask_up(n) ) then
            u1 = one - udel ; u2 = cosmum * u1 - delUdel
            spar1 = u2 * U_TPOS2_Gp1(UM,N) + u1 * U_TPOS1_Gp1(UM,N)
            spar2 = u2 * U_TPOS2_Gp2(UM,N) + u1 * U_TPOS1_Gp2(UM,N)
            L_layer_tsup_up_Gp1(um,n) = L_layer_tsup_up_Gp1(um,n) + spar1 * fac
            L_layer_tsup_up_Gp2(um,n) = L_layer_tsup_up_Gp2(um,n) + spar2 * fac
         endif

!  Downwelling 

         if ( do_dnwelling .and. layermask_dn(n) ) then
!            d1 = one - udel ; d2 = - cosmum * d1 - deltaus(n)        ! Bug 4/20/16
            d1 = one - udel ; d2 = - cosmum * d1 + deltaus(n)
            spar1 = d2 * U_TNEG2_Gp1(UM,N) + d1 * U_TNEG1_Gp1(UM,N)
            spar2 = d2 * U_TNEG2_Gp2(UM,N) + d1 * U_TNEG1_Gp2(UM,N)
            L_layer_tsup_dn_Gp1(um,n) = L_layer_tsup_dn_Gp1(um,n) + spar1 * fac
            L_layer_tsup_dn_Gp2(um,n) = L_layer_tsup_dn_Gp2(um,n) + spar2 * fac
         endif

!  End whole layer loop

      enddo

!  Partial layer loop. Checked 4/21/16

      if ( n_partlayers .ne. 0 ) then
         do ut = 1, n_PARTLAYERS
            n = partlayers_layeridx(ut)

!  Upwelling

            if ( do_upwelling ) then
               Uxup = t_utup_userm(ut,um) ; u1 = one - Uxup ; delUdel = partaus(ut) - deltaus(n) * Uxup
               u2 = cosmum * u1 + delUdel
!               u2 = cosmum * u1 - delUdel     BUG 4/21/16
               spar1 = u2 * U_TPOS2_Gp1(UM,N) + u1 * U_TPOS1_Gp1(UM,N)
               spar2 = u2 * U_TPOS2_Gp2(UM,N) + u1 * U_TPOS1_Gp2(UM,N)
               L_layer_tsup_utup_Gp1(um,ut) = L_layer_tsup_utup_Gp1(um,ut) + spar1 * fac
               L_layer_tsup_utup_Gp2(um,ut) = L_layer_tsup_utup_Gp2(um,ut) + spar2 * fac
!if (ut.eq.1) write(*,*)'Linear',ut,n,um,L_layer_tsup_utup_Gp1(um,ut),L_layer_tsup_utup_Gp2(um,ut)
            endif

!  Dnwelling

            if ( do_dnwelling ) then
               Uxdn = t_utdn_userm(ut,um) ; U1 = one - Uxdn
               u2 = - cosmum*U1 + partaus(ut) 
               spar1 = u2 * U_TNEG2_Gp1(UM,N) + u1 * U_TNEG1_Gp1(UM,N)
               spar2 = u2 * U_TNEG2_Gp2(UM,N) + u1 * U_TNEG1_Gp2(UM,N)
!               spar1 = u2 * U_TPOS2_Gp1(UM,N) + u1 * U_TPOS1_Gp1(UM,N) ! Bug 4/21/16
!               spar2 = u2 * U_TPOS2_Gp2(UM,N) + u1 * U_TPOS1_Gp2(UM,N) ! Bug 4/21/16
               L_layer_tsup_utdn_Gp1(um,ut) = L_layer_tsup_utdn_Gp1(um,ut) + spar1 * fac
               L_layer_tsup_utdn_Gp2(um,ut) = L_layer_tsup_utdn_Gp2(um,ut) + spar2 * fac
            endif

!  End partial layer loop

         enddo
      endif

!  End user-stream loop

   ENDDO

!  Continuation point

69 continue

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        START MAIN LOOP OVER LEVEL BBFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   DO JB = 0, NLAYERS

!  Solve BVProblem
!  ===============

!  Vector upgrade
!    -- Although the LBBF particular thermal solution is scalar, BRDF may induce polarization
!    -- Additional polarization (I/Q only) comes through the MS solutions, after solving BVP

!  Initialize

      N = JB + 1 ; N1 = JB
      COL2_BWF = zero

!  Skip BVP for tranmsittance only

      if ( do_thermal_transonly ) go to 75

!  Surface terms. Down = surface downwelling dependence
!    -- 1/31/21. Version 2.8.3. BRDF arrays are defined locally, drop "M" Fourier index

      Reflec = zero
      IF ( DO_INCLUDE_SURFACE.and.JB.ge.nlayers-1 ) THEN
        Down = zero
        if ( jb .eq. nlayers ) then
           do j = 1, nstreams
              Down(j) = L_T_WLOWER_Gp2(j,n1) * QUAD_STRMWTS(J)
           enddo
        else if (jb .eq. nlayers - 1 ) then
           do j = 1, nstreams
              Down(j) = L_T_WLOWER_Gp1(j,n) * QUAD_STRMWTS(J)
           enddo
        endif
         IF ( DO_LAMBERTIAN_SURFACE  ) THEN
            FACTOR = SURFACE_FACTOR * ALBEDO * sum(Down(1:nstreams))
            Reflec(1:nstreams,O11) = FACTOR
         ELSE
            do i = 1, nstreams
               do o2 = 1, nstokes
                  om = mueller_index(o2,o11)
                  FACTOR = SURFACE_FACTOR * Dot_Product(Down(1:nstreams),brdf_f(om,i,1:nstreams))
                  Reflec(i,o2) = FACTOR
               enddo
            enddo
         ENDIF
      ENDIF

!  BVProblem, Special Case, N = 2

      if ( nlayers .eq. 2 ) then
         if ( JB.eq.0 ) then                 ! Correct 3/19
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1
               COL2_BWF(ir) = - L_T_WUPPER_Gp1(i,n)
            enddo
            CM = nstks_nstrms
            do i = 1, nstreams_2
               ir = nstokes*(i-1) + 1 ; ic = cm + ir
               COL2_BWF(ic)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.1 ) then            ! Correct 3/19
            cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms ; cmp = cm + nstks_nstrms_2
            do i = 1, nstreams_2
               ir = nstokes*(i-1) + 1 ; ic = cm + ir ; ic1 = cmp + ir
               COL2_BWF(ic1)  = L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1 ; ic = cmp + ir + nstks_nstrms_2 ; i1 = i + nstreams
               COL2_BWF(ir) = - L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i,o11)
               if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
                 do o2 = 2, nstokes
                    ic2 = ic + (o2-1)
                    COL2_BWF(ic2) = Reflec(i,o2)
                 enddo
               endif
            enddo
         else if ( JB.eq.NLAYERS ) then     ! Correct 3/19
            cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms ; cmp = cm + nstks_nstrms_2
            do i = 1, nstreams_2
               ir = nstokes*(i-1) + 1 ; ic = cm + ir
               COL2_BWF(ic)  = + L_T_WUPPER_Gp1(i,n1)
            enddo
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1 ; ic = cmp + ir ; i1 = i + nstreams
               COL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i,o11)
               if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
                 do o2 = 2, nstokes
                    ic2 = ic + (o2-1)
                    COL2_BWF(ic2) = Reflec(i,o2)
                 enddo
               endif
            enddo
         endif
      endif

!  BVProblem, Single layer case
!   Bug 4/9/19. Special Case, N = 1 (JB = 0), N1 = 1 (JB = 1)

      if ( nlayers .eq. 1 ) then
         if ( JB.eq.0 ) then
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1 ; i1 = i + nstreams ; ic = ir + nstks_nstrms
               SCOL2_BWF(ir) = - L_T_WUPPER_Gp1(i,n)
               SCOL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i,o11)
               if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
                 do o2 = 2, nstokes
                    ic2 = ic + o2 - 1
                    SCOL2_BWF(ic2) = Reflec(i,o2)
                 enddo
               endif
            enddo
         else if ( JB.eq.1 ) then
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1 ; i1 = i + nstreams ; ic = ir + nstks_nstrms
               SCOL2_BWF(ir) = - L_T_WUPPER_Gp2(i,n1)                    ! Bug fix n --> n1 
               SCOL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i,o11)   ! Bug fix n --> n1
               if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
                 do o2 = 2, nstokes
                    ic2 = ic + o2 - 1
                    SCOL2_BWF(ic2) = Reflec(i,o2)
                 enddo
               endif
            enddo
         endif
      endif

!  General Case, N > 2. Separate out the various cases

      if ( nlayers .gt. 2 ) then
         if ( JB.eq.0 ) then                 ! Correct 3/19
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1
               COL2_BWF(ir) = - L_T_WUPPER_Gp1(i,n)
            enddo
            CM = nstks_nstrms
            do i = 1, nstreams_2
               ir = nstokes*(i-1) + 1 ; ic = cm + ir
               COL2_BWF(ic)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.1 ) then            ! Correct 3/19
            cm = nstks_nstrms
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1
               COL2_BWF(ir) = - L_T_WUPPER_Gp2(i,n1)
            enddo
            do i = 1, nstreams_2
               ir = nstokes*(i-1) + 1 ; ic = cm + ir ; ic1 = ic + nstks_nstrms_2
               COL2_BWF(ic)   = + L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
               COL2_BWF(ic1)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.nlayers - 1 ) then  ! Correct 3/19
            cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms ; cmp = cm + nstks_nstrms_2
            do i = 1, nstreams_2
               ir = nstokes*(i-1) + 1 ; ic = cm + ir ; ic1 = cmp + ir
               COL2_BWF(ic)   = L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic1)  = L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1 ; i1 = i + nstreams ; ic = cmp + ir + nstks_nstrms_2
               COL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i,o11)
               if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
                 do o2 = 2, nstokes
                    ic2 = ic + (o2-1)
                    COL2_BWF(ic2) = Reflec(i,o2)
                 enddo
               endif
            enddo
         else if ( JB.eq.NLAYERS ) then     ! Correct 3/19
            cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms ; cmp = cm + nstks_nstrms_2
            do i = 1, nstreams_2
               ir = nstokes*(i-1) + 1 ; ic = cm + ir
               COL2_BWF(ic)  = + L_T_WUPPER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               ir = nstokes*(i-1) + 1 ; i1 = i + nstreams ; ic = cmp + ir
               COL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i,o11)
               if ( nstokes.gt.1.and..not.do_lambertian_surface ) then
                 do o2 = 2, nstokes
                    ic2 = ic + (o2-1)
                    COL2_BWF(ic2) = Reflec(i,o2)
                 enddo
               endif
            enddo
         else                               ! Correct 3/19
            cm = JB * nstks_nstrms_2 - 3 * nstks_nstrms
            do i = 1, nstreams_2
               ir = nstokes*(i-1) + 1 ; ic = cm + ir ; ic1 = ic + nstks_nstrms_2 ; ic2 = ic1 + nstks_nstrms_2
               COL2_BWF(ic)   = + L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic1)  = + L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
               COL2_BWF(ic2)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         endif
      endif

!  debug  BVP linearization. 19 March 2014
!      if ( jb.eq.m) then
!         do n = 1, ntotal
!            write(24,*)n,COL2_BWF(n),COL2_BWF(n)
!         enddo
!      endif

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

         CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_BWF, MAXTOTAL, INFO )

!  Exception handling

         IF ( INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = ' for Atmos BBF Level '//CN//' DGBTRS call in LBBF_Jacobians'
            STATUS  = VLIDORT_SERIOUS ; RETURN
         ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

         DO N = 1, NLAYERS
            C0 = ( N - 1 ) * NSTKS_NSTRMS_2
            KO1 = K_REAL(N) + 1
            DO K = 1, K_REAL(N)
               IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
               NCON(K,N) = COL2_BWF(C0+IROW)
               PCON(K,N) = COL2_BWF(C0+IROW1)
            ENDDO
            DO K = 1, K_COMPLEX(N)
               K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
               IROW    = K    + K_REAL(N)       ; IROW1   = IROW   + NSTKS_NSTRMS
               IROW_S  = IROW + K_COMPLEX(N)    ; IROW1_S = IROW_S + NSTKS_NSTRMS
               NCON(K1,N) = COL2_BWF(C0+IROW)   ; NCON(K2,N) = COL2_BWF(C0+IROW_S)
               PCON(K1,N) = COL2_BWF(C0+IROW1)  ; PCON(K2,N) = COL2_BWF(C0+IROW1_S)
            ENDDO
         ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

         CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
                       SCOL2_BWF, MAXSTRMSTKS_2, INFO )

!  Exception handling

         IF ( INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = ' for BBF Level '//CN//' DGBTRS call in 1-layer LBBF_Jacobians'
            STATUS  = VLIDORT_SERIOUS ; RETURN
         ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

         N = 1 ;  KO1 = K_REAL(N) + 1
         DO K = 1, K_REAL(N)
            IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
            NCON(K,N) = SCOL2_BWF(IROW)
            PCON(K,N) = SCOL2_BWF(IROW1)
         ENDDO
         DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
            IROW    = K    + K_REAL(N)       ; IROW1   = IROW   + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)    ; IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON(K1,N) = SCOL2_BWF(IROW)     ; NCON(K2,N) = SCOL2_BWF(IROW_S)
            PCON(K1,N) = SCOL2_BWF(IROW1)    ; PCON(K2,N) = SCOL2_BWF(IROW1_S)
         ENDDO

!  End BV Problem

      ENDIF

!  Upwelling Jacobians
!  ===================

!  Continuation point for thermal transmittance only

75    continue

!  Skip if not applicable

      IF ( .not. DO_POSTPROCESSING ) GO TO 344
      IF ( .NOT. DO_UPWELLING )      GO TO 344

!  Skip BOA terms if no surface

      L_BOA_MSSOURCE       = zero
      L_BOA_THTONLY_SOURCE = zero
      IF ( .not. DO_INCLUDE_SURFACE ) go to 76

!  Reflected  downwelling solution 
!  -------------------------------

!   Distinguish between thermal transmittance only, and scattered solution

      IF ( DO_THERMAL_TRANSONLY ) THEN
         L_IDOWN = zero ; O1 = 1
         do n = 1, nlayers
            L_IDOWN(1:nstreams,o1) = L_IDOWN(1:nstreams,o1) * T_DELT_DISORDS(1:nstreams,N)
            if ( JB.eq.n )     L_IDOWN(1:nstreams,o1) = L_IDOWN(1:nstreams,o1) + L_T_WLOWER_GP2(1:nstreams,N)
            if ( JB.eq.n - 1 ) L_IDOWN(1:nstreams,o1) = L_IDOWN(1:nstreams,o1) + L_T_WLOWER_GP1(1:nstreams,N)
         enddo
         L_IDOWN(1:nstreams,o1) = L_IDOWN(1:nstreams,o1) * quad_weights(1:nstreams)
      ELSE
         N = NLAYERS
         do i = 1, nstreams
            do o1 = 1, nstokes
!  ...Zero contributions 
               SPAR = zero ; SHOM_R = zero ; SHOM_CR = zero
!  ..thermal particular integral; scalar component only
               if ( o1.eq.1.and.JB.eq.nlayers )     SPAR = L_T_WLOWER_GP2(i,N)
               if ( o1.eq.1.and.JB.eq.nlayers - 1 ) SPAR = L_T_WLOWER_GP1(i,N)
!  .. Homogeneous Real part
               DO K = 1, K_REAL(N)
                HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N)* T_DELT_EIGEN(K,N)
                HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO
!  .. Homogeneous Complex part
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NXR1    = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                  NXR2    = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                  HOM2CR  = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                  HOM1CR =  NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
                  SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO
!  .. Gather solution
               SHOM = SHOM_R + SHOM_CR
               L_IDOWN(I,O1) = ( SPAR + SHOM ) * QUAD_STRMWTS(I)
            ENDDO
         ENDDO
      ENDIF

!  BOA MS source terms
!  -------------------

!  Polarized calculation for the BRDF surface.
!    -- 1/31/21. Version 2.8.3. BRDF arrays are defined locally, drop "M" Fourier index

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
         FACTOR = SURFACE_FACTOR * ALBEDO * SUM(L_IDOWN(1:nstreams,o11))
         L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) = FACTOR
         IF ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:NSTREAMS) =  FACTOR
      ELSE
         DO UM = 1, N_USER_STREAMS
            do o1 = 1, nstokes
               FACTOR = ZERO
               DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                     OM = MUELLER_INDEX(O1,O2)
                     S_REFLEC = S_REFLEC + L_IDOWN(J,O2) * UBRDF_F(OM,UM,J)
                  ENDDO
                  FACTOR = FACTOR + S_REFLEC
               ENDDO
               L_BOA_MSSOURCE(UM,o1) = SURFACE_FACTOR * FACTOR
            enddo
         ENDDO
         IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              FACTOR = DOT_PRODUCT(L_IDOWN(1:nstreams,O11),BRDF_F(O11,I,1:NSTREAMS))
              L_BOA_THTONLY_SOURCE(I) =  SURFACE_FACTOR * FACTOR
            ENDDO
         ENDIF
      ENDIF

!  continuation point

76    continue

!  Recursion Loop for linearized Post-processed Jacobians (upwelling)
!  ------------------------------------------------------------------

!  Set the cumulative source term equal to BOA values

      L_CUMULSOURCE = zero
      DO UM = 1, N_USER_STREAMS
         L_CUMULSOURCE(UM,1:nstokes) = L_BOA_MSSOURCE(UM,1:nstokes) 
      ENDDO

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

         NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
         NUT = NLEVEL + 1

!  Cumulative lop to include source terms, to layer NUT

         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

!  Homogeneous Solution Contributions (only present with scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  DO o1 = 1, nstokes
!  ....Real homogeneous solutions
                     SHOM_R = ZERO
                     DO K = 1, K_REAL(N)
                        H1 = NCON(K,N) * UHOM_UPDN(UM,O1,K,N) * HMULT_2(K,UM,N)
                        H2 = PCON(K,N) * UHOM_UPUP(UM,O1,K,N) * HMULT_1(K,UM,N)
                        SHOM_R = SHOM_R + H1 + H2
                     ENDDO
!  ....Complex homogeneous solutions
                     SHOM_CR = ZERO
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        NUXR1 = NCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                        NUXR2 = NCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                        PUXR1 = PCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                        PUXR2 = PCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                        H1 = NUXR1 * HMULT_2(K1,UM,N) - NUXR2 * HMULT_2(K2,UM,N)
                        H2 = PUXR1 * HMULT_1(K1,UM,N) - PUXR2 * HMULT_1(K2,UM,N)
                        SHOM_CR = SHOM_CR + H1 + H2
                     ENDDO
                     L_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR
                  ENDDO
               ENDDO
            endif

!  Add thermal emission source terms (Includes direct and diffuse). Unpolarized
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     -----Only present for layers N adjacent to the level JB that is varying

            o1 = 1 ; TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + L_LAYER_TSUP_UP_Gp1(UM,N)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + L_LAYER_TSUP_UP_Gp2(UM,N)*TM
               ENDDO
            endif

!  Add to Linearized cumulative source sterm

            DO UM = 1, N_USER_STREAMS
               DO o1 = 1, nstokes
                  L_CUMULSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,o1)
               ENDDO
            ENDDO

!  End layer recursion loop

         ENDDO

!  User-defined stream output
!  --------------------------

!  Offgrid output

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

   !  Homogeneous Solution Contributions (only present with scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  DO o1 = 1, nstokes
!  ....Real homogeneous solutions
                     SHOM_R = ZERO
                     DO K = 1, K_REAL(N)
                        H1 = NCON(K,N) * UHOM_UPDN(UM,O1,K,N) * UT_HMULT_UD(K,UM,UT)
                        H2 = PCON(K,N) * UHOM_UPUP(UM,O1,K,N) * UT_HMULT_UU(K,UM,UT)
                        SHOM_R = SHOM_R + H1 + H2
                     ENDDO
!  ....Complex homogeneous solutions
                     SHOM_CR = ZERO
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        NUXR1 = NCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                        NUXR2 = NCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                        PUXR1 = PCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                        PUXR2 = PCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                        H1 = NUXR1 * UT_HMULT_UD(K1,UM,UT) - NUXR2 * UT_HMULT_UD(K2,UM,UT)
                        H2 = PUXR1 * UT_HMULT_UU(K1,UM,UT) - PUXR2 * UT_HMULT_UU(K2,UM,UT)
                        SHOM_CR = SHOM_CR + H1 + H2
                     ENDDO
                     L_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR
                  ENDDO
               ENDDO
            endif

!  Add thermal emission term (direct and diffuse)

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM,o11) = L_LAYERSOURCE(UM,o11) + L_LAYER_TSUP_UTUP_Gp1(UM,UT)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM,o11) = L_LAYERSOURCE(UM,o11) + L_LAYER_TSUP_UTUP_Gp2(UM,UT)*TM
               ENDDO
            endif

!  Final answer

            DO UM = 1, N_USER_STREAMS
               DO o1 = 1, nstokes
                  FINAL_SOURCE = L_LAYERSOURCE(UM,O1) + T_UTUP_USERM(UT,UM) * L_CUMULSOURCE(UM,O1)
                  ABBWFS_JACOBIANS(UTA,UM,JB,O1,UPIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
               ENDDO
            ENDDO

!  Ongrid output, just set to the cumulative source term

         ELSE
            DO UM = 1, N_USER_STREAMS
               DO o1 = 1, nstokes
                  ABBWFS_JACOBIANS(UTA,UM,JB,O1,UPIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,O1)
               ENDDO
            ENDDO
         ENDIF

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  continuation point

344   continue

!  Downwelling Jacobians
!  =====================

!  Skip if not applicable

      IF ( .not. DO_POSTPROCESSING ) GO TO 345
      IF ( .NOT. DO_DNWELLING )      GO TO 345

!  Initialize post-processing recursion
!  Set the cumulative source term equal to TOA values

      L_CUMULSOURCE = zero

!  Recursion Loop for linearized Post-processed Jacobians (downwelling)
!  --------------------------------------------------------------------

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

         NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
 
         NUT = NLEVEL
         DO N = NSTART, NUT
            NC = N

!  Homogeneous Solution Contributions (only present with scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  DO o1 = 1, nstokes
!  ....Real homogeneous solutions
                     SHOM_R = ZERO
                     DO K = 1, K_REAL(N)
                        H1 = NCON(K,N) * UHOM_DNDN(UM,O1,K,N) * HMULT_1(K,UM,N)
                        H2 = PCON(K,N) * UHOM_DNUP(UM,O1,K,N) * HMULT_2(K,UM,N)
                        SHOM_R = SHOM_R + H1 + H2
                     ENDDO
!  ....Complex homogeneous solutions
                     SHOM_CR = ZERO
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        NUXR1 = NCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
                        NUXR2 = NCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
                        PUXR1 = PCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
                        PUXR2 = PCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_DnUP(UM,O1,K1,N)
                        H1 = NUXR1 * HMULT_1(K1,UM,N) - NUXR2 * HMULT_1(K2,UM,N)
                        H2 = PUXR1 * HMULT_2(K1,UM,N) - PUXR2 * HMULT_2(K2,UM,N)
                        SHOM_CR = SHOM_CR + H1 + H2
                     ENDDO
                     L_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR
                  ENDDO
               ENDDO
            endif
 
!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     -----Only if adjacent layers to the level that is varying

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one / PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM,o11) = L_LAYERSOURCE(UM,o11) + L_LAYER_TSUP_DN_Gp1(UM,N)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM,o11) = L_LAYERSOURCE(UM,o11) + L_LAYER_TSUP_DN_Gp2(UM,N)*TM
               ENDDO
            endif

!  Add to Linearized cumulative source sterm

            DO UM = 1, N_USER_STREAMS
               do o1 = 1, nstokes
                  L_CUMULSOURCE(UM,o1) = L_LAYERSOURCE(UM,o1) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,o1)
               ENDDO
            ENDDO

!  End layer loop

         ENDDO

!  User-defined stream output
!  --------------------------

!  Offgrid output

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

   !  Homogeneous Solution Contributions (only present with scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  DO o1 = 1, nstokes
!  ....Real homogeneous solutions
                     SHOM_R = ZERO
                     DO K = 1, K_REAL(N)
!                        H1 = NCON(K,N) * UHOM_UPDN(UM,O1,K,N) * UT_HMULT_DD(K,UM,UT)  ! BUG 4/21/16
!                        H2 = PCON(K,N) * UHOM_UPUP(UM,O1,K,N) * UT_HMULT_DU(K,UM,UT)  ! BUG 4/21/16
                        H1 = NCON(K,N) * UHOM_DNDN(UM,O1,K,N) * UT_HMULT_DD(K,UM,UT)
                        H2 = PCON(K,N) * UHOM_DNUP(UM,O1,K,N) * UT_HMULT_DU(K,UM,UT)
                        SHOM_R = SHOM_R + H1 + H2
                     ENDDO
!  ....Complex homogeneous solutions
                     SHOM_CR = ZERO
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        NUXR1 = NCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
                        NUXR2 = NCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
                        PUXR1 = PCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
                        PUXR2 = PCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
                        H1 = NUXR1 * UT_HMULT_DD(K1,UM,UT) - NUXR2 * UT_HMULT_DD(K2,UM,UT)
                        H2 = PUXR1 * UT_HMULT_DU(K1,UM,UT) - PUXR2 * UT_HMULT_DU(K2,UM,UT)
                        SHOM_CR = SHOM_CR + H1 + H2
                     ENDDO
                     L_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR
                  ENDDO
               ENDDO
            endif

!  Add thermal emission term (direct and diffuse)

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM,o11) = L_LAYERSOURCE(UM,o11) + L_LAYER_TSUP_UTDN_Gp1(UM,UT)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM,o11) = L_LAYERSOURCE(UM,o11) + L_LAYER_TSUP_UTDN_Gp2(UM,UT)*TM
               ENDDO
            endif

!  Final answer

            DO UM = 1, N_USER_STREAMS
               DO o1 = 1, nstokes
                  FINAL_SOURCE = L_LAYERSOURCE(UM,O1) + T_UTDN_USERM(UT,UM) * L_CUMULSOURCE(UM,O1)
                  ABBWFS_JACOBIANS(UTA,UM,JB,O1,DNIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
               ENDDO
            ENDDO

!  Ongrid output, just set to the cumulative source term

         ELSE
            DO UM = 1, N_USER_STREAMS
               DO o1 = 1, nstokes
                  ABBWFS_JACOBIANS(UTA,UM,JB,O1,DNIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,O1)
               ENDDO
            ENDDO
         ENDIF

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  continuation point

345   continue

!  Flux Jacobians
!  ==============

!  Upwelling FLux output
!  ---------------------

      if ( DO_MVOUTPUT .and. DO_UPWELLING ) THEN
         DO UTA = 1, N_USER_LEVELS
            BBWF_QUAD = ZERO ; O11 = 1

!  Partial layer output for linearized Quadrature field
!  ----------------------------------------------------

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
               UT = PARTLAYERS_OUTINDEX(UTA)
               N  = PARTLAYERS_LAYERIDX(UT)
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     QUAD = QUAD_STREAMS(I)
                     L_THELP = L_BOA_THTONLY_SOURCE(I)
                     DO LAY = NLAYERS, N+1, -1
                        SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                        if ( JB.eq.LAY-1 ) SPAR = L_T_WUPPER_Gp1(I1,LAY)
                        if ( JB.eq.LAY )   SPAR = L_T_WUPPER_Gp2(I1,LAY)
                        L_THELP = L_THELP + ( SPAR / QUAD )
                     ENDDO
                     SPAR = zero ; L_THELP = L_THELP * T_UTUP_DISORDS(I,UT)
                     if ( JB.eq.N-1 ) SPAR = L_UT_T_PARTIC_Gp1(I1,UT)
                     if ( JB.eq.N )   SPAR = L_UT_T_PARTIC_Gp2(I1,UT)
!                     BBWF_QUAD(I,O11) = FLUX_MULTIPLIER * L_THELP       ! Bug 4/21/16
                     BBWF_QUAD(I,O11) = FLUX_MULTIPLIER * ( L_THELP + (SPAR/QUAD) )
                  ENDDO
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS ; spar = zero
                     DO O1 = 1, nstokes
                        SHOM_R = ZERO ; SHOM_CR = zero
                        DO K = 1, K_REAL(N)
!   Buggy code 4/21/16
!                           HOM1 = NCON(K,N) * SOLA_XPOS(I1,O1,K,N) * T_UTUP_EIGEN(K,UT)
!                           HOM2 = PCON(K,N) * SOLB_XNEG(I1,O1,K,N)
                           HOM1 = NCON(K,N) * SOLA_XPOS(I1,O1,K,N) * T_UTDN_EIGEN(K,UT)
                           HOM2 = PCON(K,N) * SOLB_XNEG(I1,O1,K,N) * T_UTUP_EIGEN(K,UT)
                           SHOM_R = SHOM_R + HOM1 + HOM2
                        ENDDO

                        DO K = 1, K_COMPLEX(N)
                           K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1

!  Buggy code replaced, 4/21/16, but new stuff is untested
!                           NXR1   = NCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
!                           NXR2   = NCON(K1,N) * SOLA_XPOS(I1,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I1,O1,K1,N)
!                           HOM2CR = PCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
!                           HOM1CR = NXR1 * T_UTUP_EIGEN(K1,UT) - NXR2 * T_UTUP_EIGEN(K2,UT)

!  replacement, not tested yet....
                           NXR_CR = NCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                           NXR_CI = NCON(K1,N) * SOLA_XPOS(I1,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I1,O1,K1,N)
                           PXR_CR = PCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
                           PXR_CI = PCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
                           HOM1CR =  NXR_CR * T_UTDN_EIGEN(K1,UT) - NXR_CI * T_UTDN_EIGEN(K2,UT)
                           HOM2CR =  PXR_CR * T_UTUP_EIGEN(K1,UT) - PXR_CI * T_UTUP_EIGEN(K2,UT)

                           SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                        ENDDO
                        BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                     ENDDO
                     if ( JB.eq.N-1 ) SPAR = L_UT_T_PARTIC_Gp1(I1,UT)
                     if ( JB.eq.N )   SPAR = L_UT_T_PARTIC_Gp2(I1,UT)
                     BBWF_QUAD(I,o11) = BBWF_QUAD(I,o11) + FLUX_MULTIPLIER * SPAR 
                  ENDDO
               endif
            ENDIF

!  Level-boundary output of linearized quadrature field
!  ----------------------------------------------------

            IF ( .not.PARTLAYERS_OUTFLAG(UTA) ) THEN
               NL = UTAU_LEVEL_MASK_UP(UTA) ; N = NL + 1

!  quadrature field at bottom level

               IF ( NL .EQ. NLAYERS  ) THEN
                  if ( do_thermal_transonly ) then
                     DO I = 1, NSTREAMS
                        BBWF_QUAD(I,o11) = FLUX_MULTIPLIER * L_BOA_THTONLY_SOURCE(I)
                     enddo
                  else
                     DO I = 1, NSTREAMS
                        I1 = I + NSTREAMS ; spar = zero
                        DO O1 = 1, nstokes
                           SHOM_R = ZERO ; SHOM_CR = zero
                           DO K = 1, K_REAL(NL)
                              HOM1 = NCON(K,NL) * SOLA_XPOS(I1,O1,K,NL) * T_DELT_EIGEN(K,NL)
                              HOM2 = PCON(K,NL) * SOLB_XNEG(I1,O1,K,NL)
                              SHOM_R = SHOM_R + HOM1 + HOM2
                           ENDDO
                           DO K = 1, K_COMPLEX(NL)
                              K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                              NXR1   = NCON(K1,NL) * SOLA_XPOS(I1,O1,K1,NL) - NCON(K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
                              NXR2   = NCON(K1,NL) * SOLA_XPOS(I1,O1,K2,NL) + NCON(K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
                              HOM2CR = PCON(K1,NL) * SOLB_XNEG(I1,O1,K1,NL) - PCON(K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
                              HOM1CR = NXR1 * T_DELT_EIGEN(K1,NL) - NXR2 * T_DELT_EIGEN(K2,NL)
                              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                           ENDDO
                           BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                        ENDDO
                        if ( JB.eq.NL-1 ) SPAR = L_T_WLOWER_Gp1(I1,NL)
                        if ( JB.eq.NL )   SPAR = L_T_WLOWER_Gp2(I1,NL)
                        BBWF_QUAD(I,o11) = BBWF_QUAD(I,o11) + FLUX_MULTIPLIER * SPAR 
                     ENDDO
                  endif

!  Quadrature field other levels

               ELSE
                  if ( do_thermal_transonly ) then
                     DO I = 1, NSTREAMS
                        I1 = I + NSTREAMS
                        L_THELP = L_BOA_THTONLY_SOURCE(I)
                        DO LAY = NLAYERS, N, -1
                           SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                           if ( JB.eq.LAY-1 ) SPAR = L_T_WUPPER_Gp1(I1,LAY)
                           if ( JB.eq.LAY )   SPAR = L_T_WUPPER_Gp2(I1,LAY)
                           L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                        ENDDO
                        BBWF_QUAD(I,o11) = FLUX_MULTIPLIER * L_THELP
                     enddo
                  else
                     DO I = 1, NSTREAMS
                        I1 = I + NSTREAMS ; spar = zero
                        DO O1 = 1, nstokes
                           SHOM_R = ZERO ; SHOM_CR = zero
                           DO K = 1, K_REAL(N)
                              HOM1 = NCON(K,N) * SOLA_XPOS(I1,O1,K,N)
                              HOM2 = PCON(K,N) * SOLB_XNEG(I1,O1,K,N) * T_DELT_EIGEN(K,N)
                              SHOM_R = SHOM_R + HOM1 + HOM2
                           ENDDO
                           DO K = 1, K_COMPLEX(N)
                              K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                              HOM1CR = NCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                              PXR1   = PCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - PCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                              PXR2   = PCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
                              HOM2CR = PXR1 * T_DELT_EIGEN(K1,N) - PXR2 * T_DELT_EIGEN(K2,N)
                              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                           ENDDO
                           BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                        ENDDO
                        if ( JB.eq.N-1 ) SPAR = L_T_WUPPER_Gp1(I1,N)
                        if ( JB.eq.N )   SPAR = L_T_WUPPER_Gp2(I1,N)
                        BBWF_QUAD(I,o11) = BBWF_QUAD(I,o11) + FLUX_MULTIPLIER * SPAR 
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF

!  Set fluxes

            do o1 = 1, nstokes
               SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_WEIGHTS(1:NSTREAMS))
               SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS))
               ABBWFS_FLUXES(UTA,1,JB,O1,UPIDX) = HALF * SM
               ABBWFS_FLUXES(UTA,2,JB,O1,UPIDX) = PI2  * SF
            enddo

!  End output level loop

         ENDDO
      ENDIF

!  Downwelling FLux output
!  -----------------------

      if ( DO_MVOUTPUT .and. DO_DNWELLING ) THEN
         DO UTA = 1, N_USER_LEVELS
            BBWF_QUAD = ZERO

!  Partial layer output for linearized Quadrature field
!  ----------------------------------------------------

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
               UT = PARTLAYERS_OUTINDEX(UTA)
               N  = PARTLAYERS_LAYERIDX(UT)
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     L_THELP = ZERO
                     DO LAY = 1, N - 1
                        SPAR = ZERO ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                        if ( JB.eq.LAY-1 ) SPAR = L_T_WLOWER_Gp1(I,LAY)
                        if ( JB.eq.LAY )   SPAR = L_T_WLOWER_Gp2(I,LAY)
                        L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                     ENDDO
                     SPAR = ZERO ; L_THELP = L_THELP * T_UTDN_DISORDS(I,UT)
                     if ( JB.eq.N-1 ) SPAR = L_UT_T_PARTIC_Gp1(I,UT)
                     if ( JB.eq.N )   SPAR = L_UT_T_PARTIC_Gp2(I,UT)
                     L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                     BBWF_QUAD(I,O11) = FLUX_MULTIPLIER * L_THELP
                  enddo
               else
                  DO I = 1, NSTREAMS
                     spar = zero
                     DO O1 = 1, nstokes
                        SHOM_R = ZERO ; SHOM_CR = zero
                        DO K = 1, K_REAL(N)
                           HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N) * T_UTDN_EIGEN(K,UT)
!                           HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)       ! Bug 4/21/16
                           HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N) * T_UTUP_EIGEN(K,UT)
                           SHOM_R = SHOM_R + HOM1 + HOM2
                        ENDDO
                        DO K = 1, K_COMPLEX(N)
                           K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1

!  Buggy Code 4/21/16, replaced but not tested

!                           NXR1   = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
!                           NXR2   = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
!                           HOM2CR = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
!                           HOM1CR = NXR1 * T_UTDN_EIGEN(K1,UT) - NXR2 * T_UTDN_EIGEN(K2,UT)

!  replacement, not tested yet....
                           NXR_CR = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                           NXR_CI = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                           PXR_CR = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                           PXR_CI = PCON(K1,N) * SOLB_XNEG(I,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I,O1,K1,N)
                           HOM1CR =  NXR_CR * T_UTDN_EIGEN(K1,UT) - NXR_CI * T_UTDN_EIGEN(K2,UT)
                           HOM2CR =  PXR_CR * T_UTUP_EIGEN(K1,UT) - PXR_CI * T_UTUP_EIGEN(K2,UT)

                           SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                        ENDDO
                        BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                     ENDDO
                     if ( JB.eq.N-1 ) SPAR = L_UT_T_PARTIC_Gp1(I,UT)
                     if ( JB.eq.N )   SPAR = L_UT_T_PARTIC_Gp2(I,UT)
                     BBWF_QUAD(I,o11) = BBWF_QUAD(I,o11) + FLUX_MULTIPLIER * SPAR 
                  ENDDO
               endif
            ENDIF

!  Level-boundary output of linearized quadrature field
!  ----------------------------------------------------

            IF ( .not.PARTLAYERS_OUTFLAG(UTA) ) THEN
               NL = UTAU_LEVEL_MASK_DN(UTA) ; N = NL

!  Quadrature field at other levels than top

               IF ( NL .NE. 0  ) THEN
                  if ( do_thermal_transonly ) then
                     DO I = 1, NSTREAMS
                        L_THELP = ZERO
                        DO LAY = 1, NL
                           SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                           if ( JB.eq.LAY-1 ) SPAR = L_T_WLOWER_Gp1(I,LAY)
                           if ( JB.eq.LAY )   SPAR = L_T_WLOWER_Gp2(I,LAY)
                           L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                        ENDDO
                        BBWF_QUAD(I,o11) = FLUX_MULTIPLIER * L_THELP
                     enddo
                  else
                     DO I = 1, NSTREAMS
                        spar = zero
                        DO O1 = 1, nstokes
                           SHOM_R = ZERO ; SHOM_CR = zero
                           DO K = 1, K_REAL(N)
                              HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N) * T_DELT_EIGEN(K,N)
                              HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)
                              SHOM_R = SHOM_R + HOM1 + HOM2
                           ENDDO
                           DO K = 1, K_COMPLEX(N)
                              K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                              NXR1   = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                              NXR2   = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                              HOM2CR = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                              HOM1CR = NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
                              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                           ENDDO
                           BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                        ENDDO
                        if ( JB.eq.N-1 ) SPAR = L_T_WLOWER_Gp1(I,N)
                        if ( JB.eq.N )   SPAR = L_T_WLOWER_Gp2(I,N)
                        BBWF_QUAD(I,o11) = BBWF_QUAD(I,o11) + FLUX_MULTIPLIER * SPAR 
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF

!  Set Fluxes

            do o1 = 1, nstokes
               SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_WEIGHTS(1:NSTREAMS))
               SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS))
               ABBWFS_FLUXES(UTA,1,JB,O1,DNIDX) = HALF * SM
               ABBWFS_FLUXES(UTA,2,JB,O1,DNIDX) = PI2  * SF
            enddo

!  End output level loop

         ENDDO
      ENDIF

!  End loop over BBWFS

   enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      SURFACE LBBF JACOBIANS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Continuation point

55 continue

!  Finish if not required

   if ( .not. DO_SURFACE_LBBF ) RETURN

!  Solve BVProblem
!  ===============

!  Skip BVP for thermal Transmittance only

   if ( do_thermal_transonly ) go to 79

!  Initialize and Set the Column vector

   if ( nlayers .gt.1 ) then
      COL2_BWF = zero ; C0 = NLAYERS * NSTKS_NSTRMS_2 - NSTKS_NSTRMS
      DO I = 1, NSTREAMS
         DO O1 = 1, nstokes
            IR = NSTOKES*(I-1) + O1 ; CM = C0 + IR
            COL2_BWF(CM) = EMISSIVITY(O1,I)
         ENDDO
      ENDDO
   else
      SCOL2_BWF = zero ; C0 = NSTKS_NSTRMS
      DO I = 1, NSTREAMS
         DO O1 = 1, nstokes
            IR = NSTOKES*(I-1) + O1 ; CM = C0 + IR
            SCOL2_BWF(CM) = EMISSIVITY(O1,I)
         ENDDO
      ENDDO
   endif

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

   IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

      CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_BWF, MAXTOTAL, INFO )

!  Exception handling

      IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = ' for Surface BBF, DGBTRS call in LBBF_Jacobians'
         STATUS  = VLIDORT_SERIOUS ; RETURN
      ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

      DO N = 1, NLAYERS
         C0 = ( N - 1 ) * NSTKS_NSTRMS_2
         KO1 = K_REAL(N) + 1
         DO K = 1, K_REAL(N)
            IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
            NCON(K,N) = COL2_BWF(C0+IROW)
            PCON(K,N) = COL2_BWF(C0+IROW1)
         ENDDO
         DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
            IROW    = K    + K_REAL(N)       ; IROW1   = IROW   + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)    ; IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON(K1,N) = COL2_BWF(C0+IROW)   ; NCON(K2,N) = COL2_BWF(C0+IROW_S)
            PCON(K1,N) = COL2_BWF(C0+IROW1)  ; PCON(K2,N) = COL2_BWF(C0+IROW1_S)
         ENDDO
      ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

   ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

      CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
                       SCOL2_BWF, MAXSTRMSTKS_2, INFO )

!  Exception handling

      IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = ' for Surface BBF, DGBTRS call in 1-layer LBBF_Jacobians'
         STATUS  = VLIDORT_SERIOUS ; RETURN
      ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

      N = 1 ;  KO1 = K_REAL(N) + 1
      DO K = 1, K_REAL(N)
         IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
         NCON(K,N) = SCOL2_BWF(IROW)
         PCON(K,N) = SCOL2_BWF(IROW1)
      ENDDO
      DO K = 1, K_COMPLEX(N)
         K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
         IROW    = K    + K_REAL(N)       ; IROW1   = IROW   + NSTKS_NSTRMS
         IROW_S  = IROW + K_COMPLEX(N)    ; IROW1_S = IROW_S + NSTKS_NSTRMS
         NCON(K1,N) = SCOL2_BWF(IROW)     ; NCON(K2,N) = SCOL2_BWF(IROW_S)
         PCON(K1,N) = SCOL2_BWF(IROW1)    ; PCON(K2,N) = SCOL2_BWF(IROW1_S)
      ENDDO

!  End choice over number of layers

   ENDIF

!  Continuation point for skipping BVP

79 continue

!  Upwelling Jacobians
!  ===================

!  Skip if not applicable

   IF ( .not. DO_POSTPROCESSING ) GO TO 544
   IF ( .NOT. DO_UPWELLING )      GO TO 544

!  BOA source terms

   L_BOA_MSSOURCE = zero
   IF ( DO_INCLUDE_SURFACE ) THEN

      N = NLAYERS
      L_IDOWN              = zero ! L_Down is zero for thermal transmittance only
      L_BOA_THTONLY_SOURCE = zero ! Only non-zero if thermal transmittance only (and DO_QTHTONLY)

      if ( .not. do_thermal_transonly ) then
         do i = 1, nstreams
            do o1 = 1, nstokes
               SHOM_R = zero ; SHOM_CR = zero
!  .. Homogeneous Real part
               DO K = 1, K_REAL(N)
                HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N)* T_DELT_EIGEN(K,N)
                HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO
!  .. Homogeneous Complex part
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NXR1    = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                  NXR2    = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                  HOM2CR  = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                  HOM1CR =  NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
                  SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO
!  .. Gather solution
               SHOM = SHOM_R + SHOM_CR
               L_IDOWN(I,O1) = SHOM *QUAD_STRMWTS(I)
            ENDDO
         ENDDO
      endif

!  Surface Source Term. Polarized calculation for the BRDF surface.
!   R. Spurr, Bug Fix, 4/9/19. User_Emissivity should not be used when MS-only is true.
!    -- 1/31/21. Version 2.8.3. BRDF arrays are defined locally, drop "M" Fourier index

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
         FACTOR = SURFACE_FACTOR * ALBEDO * SUM(L_IDOWN(1:nstreams,o11)) ; EMISS = one - ALBEDO
         L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) = FACTOR
!  Bug   L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) = L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) + EMISS
         IF (.not.DO_MSMODE_THERMAL) L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) = L_BOA_MSSOURCE(1:N_USER_STREAMS,o11) + EMISS
         if ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:nstreams) = EMISS
      ELSE
         DO UM = 1, N_USER_STREAMS
            do o1 = 1, nstokes
               FACTOR = ZERO
               DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                     OM = MUELLER_INDEX(O1,O2)
                     S_REFLEC = S_REFLEC + L_IDOWN(J,O2) * UBRDF_F(OM,UM,J)
                  ENDDO
                  FACTOR = FACTOR + S_REFLEC
               ENDDO
!  Bug         L_BOA_MSSOURCE(UM,o1) = SURFACE_FACTOR * FACTOR + USER_EMISSIVITY(O1,UM)
               IF (.not.DO_MSMODE_THERMAL) L_BOA_MSSOURCE(UM,o1) = SURFACE_FACTOR * FACTOR + USER_EMISSIVITY(O1,UM)
               L_BOA_MSSOURCE(UM,o1) = SURFACE_FACTOR * FACTOR + USER_EMISSIVITY(O1,UM)
            enddo
         ENDDO
         if ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:nstreams) = EMISSIVITY(o11,1:nstreams)
      ENDIF

   ENDIF

!  Upwelling post-processing recursion
!  -----------------------------------

!  Initialie

   DO UM = 1, N_USER_STREAMS
      L_CUMULSOURCE(UM,1:nstokes) = L_BOA_MSSOURCE(UM,1:nstokes) 
   ENDDO


   NC  = 0;  NUT = 0
   NSTART = NLAYERS ; NUT_PREV = NSTART + 1

!  Start user level output llop

   DO UTA = N_USER_LEVELS, 1, -1

      NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
      NUT = NLEVEL + 1

!  Layer cumulative terms

      DO N = NSTART, NUT, -1
         NC = NLAYERS + 1 - N
         DO UM = 1, N_USER_STREAMS
            DO O1 = 1, nstokes
               SHOM_R = ZERO ; SHOM_CR = zero
               IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!  ....Real homogeneous solutions
                  DO K = 1, K_REAL(N)
                     H1 = NCON(K,N) * UHOM_UPDN(UM,O1,K,N) * HMULT_2(K,UM,N)
                     H2 = PCON(K,N) * UHOM_UPUP(UM,O1,K,N) * HMULT_1(K,UM,N)
                     SHOM_R = SHOM_R + H1 + H2
                  ENDDO
!  ....Complex homogeneous solutions
                  DO K = 1, K_COMPLEX(N)
                     K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                     NUXR1 = NCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                     NUXR2 = NCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                     PUXR1 = PCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                     PUXR2 = PCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                     H1 = NUXR1 * HMULT_2(K1,UM,N) - NUXR2 * HMULT_2(K2,UM,N)
                     H2 = PUXR1 * HMULT_1(K1,UM,N) - PUXR2 * HMULT_1(K2,UM,N)
                     SHOM_CR = SHOM_CR + H1 + H2
                  ENDDO
               ENDIF
               L_CUMULSOURCE(UM,O1) = SHOM_R + SHOM_CR + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,O1)
            ENDDO
         ENDDO
      ENDDO

!  Offgrid (partial) output : Need to evaulate extra term
!  Layer-boundary    output : just set to the cumulative source term

      IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
         DO UM = 1, N_USER_STREAMS
            DO O1 = 1, nstokes
               SHOM_R = ZERO ; SHOM_CR = zero
               IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!  ....Real homogeneous solutions
                  DO K = 1, K_REAL(N)
                     H1 = NCON(K,N) * UHOM_UPDN(UM,O1,K,N) * UT_HMULT_UD(K,UM,UT) 
                     H2 = PCON(K,N) * UHOM_UPUP(UM,O1,K,N) * UT_HMULT_UU(K,UM,UT) 
                     SHOM_R = SHOM_R + H1 + H2
                  ENDDO
!  ....Complex homogeneous solutions
                  DO K = 1, K_COMPLEX(N)
                     K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                     NUXR1 = NCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                     NUXR2 = NCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                     PUXR1 = PCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                     PUXR2 = PCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                     H1 = NUXR1 * UT_HMULT_UD(K1,UM,UT) - NUXR2 * UT_HMULT_UD(K2,UM,UT)
                     H2 = PUXR1 * UT_HMULT_UU(K1,UM,UT) - PUXR2 * UT_HMULT_UU(K2,UM,UT)
                     SHOM_CR = SHOM_CR + H1 + H2
                  ENDDO
               ENDIF
               L_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR 
               FINAL_SOURCE = L_LAYERSOURCE(UM,O1) + T_UTUP_USERM(UT,UM) * L_CUMULSOURCE(UM,O1)
               SBBWFS_JACOBIANS(UTA,UM,O1,UPIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO
         ENDDO
      ELSE
         DO UM = 1, N_USER_STREAMS
            SBBWFS_JACOBIANS(UTA,UM,1:nstokes,UPIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,1:nstokes)
         ENDDO
      ENDIF

!  End output level loop

      IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
      NUT_PREV = NUT
   ENDDO

!  continuation point

544   continue

!  Downwelling Jacobians
!  =====================

!  Skip if not applicable

   IF ( .not. DO_POSTPROCESSING ) GO TO 545
   IF ( .NOT. DO_DNWELLING )      GO TO 545

!  Downwelling post-processing recursion
!  -------------------------------------

!  Initialize

   L_CUMULSOURCE = zero
   NC  = 0 ; NUT = 0
   NSTART = 1 ;  NUT_PREV = NSTART - 1

!  Start output level loop

   DO UTA = 1, N_USER_LEVELS
      NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
      NUT = NLEVEL

!  Layer cumulative terms

      DO N = NSTART, NUT
         NC = N
         DO UM = 1, N_USER_STREAMS
            DO O1 = 1, nstokes
               SHOM_R = ZERO ; SHOM_CR = zero
               IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!  ....Real homogeneous solutions
                  DO K = 1, K_REAL(N)
                     H1 = NCON(K,N) * UHOM_DNDN(UM,O1,K,N) * HMULT_1(K,UM,N)
                     H2 = PCON(K,N) * UHOM_DNUP(UM,O1,K,N) * HMULT_2(K,UM,N)
                     SHOM_R = SHOM_R + H1 + H2
                  ENDDO
!  ....Complex homogeneous solutions
                  DO K = 1, K_COMPLEX(N)
                     K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                     NUXR1 = NCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
                     NUXR2 = NCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
                     PUXR1 = PCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
                     PUXR2 = PCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
                     H1 = NUXR1 * HMULT_1(K1,UM,N) - NUXR2 * HMULT_1(K2,UM,N)
                     H2 = PUXR1 * HMULT_2(K1,UM,N) - PUXR2 * HMULT_2(K2,UM,N)
                     SHOM_CR = SHOM_CR + H1 + H2
                  ENDDO
               ENDIF
               L_CUMULSOURCE(UM,O1) = SHOM_R + SHOM_CR + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,O1)
            ENDDO
         ENDDO
      ENDDO

!  Offgrid (partial) output : Need to evaulate extra term
!  Layer-boundary    output : just set to the cumulative source term

      IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
         DO UM = 1, N_USER_STREAMS
            DO O1 = 1, NSTOKES
               SHOM_R = ZERO ; SHOM_CR = zero
               IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!  ....Real homogeneous solutions
                  DO K = 1, K_REAL(N)
                     H1 = NCON(K,N) * UHOM_DNDN(UM,O1,K,N) * UT_HMULT_DD(K,UM,UT) 
                     H2 = PCON(K,N) * UHOM_DNUP(UM,O1,K,N) * UT_HMULT_DU(K,UM,UT) 
                     SHOM_R = SHOM_R + H1 + H2
                  ENDDO
!  ....Complex homogeneous solutions
                  DO K = 1, K_COMPLEX(N)
                     K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                     NUXR1 = NCON(K1,N) * UHOM_DNDN(UM,O1,K1,N) - NCON(K2,N) * UHOM_DNDN(UM,O1,K2,N)
                     NUXR2 = NCON(K1,N) * UHOM_DNDN(UM,O1,K2,N) + NCON(K2,N) * UHOM_DNDN(UM,O1,K1,N)
                     PUXR1 = PCON(K1,N) * UHOM_DNUP(UM,O1,K1,N) - PCON(K2,N) * UHOM_DNUP(UM,O1,K2,N)
                     PUXR2 = PCON(K1,N) * UHOM_DNUP(UM,O1,K2,N) + PCON(K2,N) * UHOM_DNUP(UM,O1,K1,N)
                     H1 = NUXR1 * UT_HMULT_DD(K1,UM,UT) - NUXR2 * UT_HMULT_DD(K2,UM,UT)
                     H2 = PUXR1 * UT_HMULT_DU(K1,UM,UT) - PUXR2 * UT_HMULT_DU(K2,UM,UT)
                     SHOM_CR = SHOM_CR + H1 + H2
                  ENDDO
               ENDIF
               L_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR 
               FINAL_SOURCE = L_LAYERSOURCE(UM,O1) + T_UTDN_USERM(UT,UM) * L_CUMULSOURCE(UM,O1)
               SBBWFS_JACOBIANS(UTA,UM,O1,DNIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO
         ENDDO
      ELSE
         DO UM = 1, N_USER_STREAMS
            SBBWFS_JACOBIANS(UTA,UM,1:nstokes,DNIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,1:nstokes)
         ENDDO
      ENDIF

!  End output level loop

      IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
      NUT_PREV = NUT
   ENDDO

!  continuation point

545   continue

!  Flux Jacobians
!  ==============

!  Upwelling FLux output
!  ---------------------

   if ( DO_MVOUTPUT .and. DO_UPWELLING ) THEN
      DO UTA = 1, N_USER_LEVELS
         BBWF_QUAD = ZERO

!  Partial layer output for linearized Quadrature field

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            if ( do_thermal_transonly ) then
               DO I = 1, NSTREAMS
                  SHOM = L_BOA_THTONLY_SOURCE(I)
                  DO LAY = NLAYERS,  N + 1, -1
                     SHOM = SHOM * T_DELT_DISORDS(I,LAY)
                  ENDDO
                  SHOM = SHOM * T_UTUP_DISORDS(I,UT)
                  BBWF_QUAD(I,O11) = FLUX_MULTIPLIER * SHOM
               ENDDO
            else
               DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  DO O1 = 1, nstokes
                     SHOM_R = ZERO ; SHOM_CR = zero
                     DO K = 1, K_REAL(N)
!                        HOM1 = NCON(K,N) * SOLA_XPOS(I1,O1,K,N)   !  BUG 4/21/16 ????
                        HOM1 = NCON(K,N) * SOLA_XPOS(I1,O1,K,N) * T_UTDN_EIGEN(K,UT)
                        HOM2 = PCON(K,N) * SOLB_XNEG(I1,O1,K,N) * T_UTUP_EIGEN(K,UT)
                        SHOM_R = SHOM_R + HOM1 + HOM2
                     ENDDO
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1

!  Buggy code 4/21/16
!                        HOM1CR = NCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
!                        PXR1   = PCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - PCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
!                        PXR2   = PCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
!                        HOM2CR = PXR1 * T_UTUP_EIGEN(K1,UT) - PXR2 * T_UTUP_EIGEN(K2,UT)

!  replacement, not tested yet....
                        NXR_CR = NCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                        NXR_CI = NCON(K1,N) * SOLA_XPOS(I1,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I1,O1,K1,N)
                        PXR_CR = PCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
                        PXR_CI = PCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
                        HOM1CR =  NXR_CR * T_UTDN_EIGEN(K1,UT) - NXR_CI * T_UTDN_EIGEN(K2,UT)
                        HOM2CR =  PXR_CR * T_UTUP_EIGEN(K1,UT) - PXR_CI * T_UTUP_EIGEN(K2,UT)

                        SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                     ENDDO
                     BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                  ENDDO
               ENDDO
            endif
         ENDIF

!  Level-boundary output of linearized quadrature field

         IF ( .not.PARTLAYERS_OUTFLAG(UTA) ) THEN
            NL = UTAU_LEVEL_MASK_UP(UTA) ; N = NL + 1
            IF ( NL .EQ. NLAYERS  ) THEN
               if ( do_thermal_transonly ) then
                  BBWF_QUAD(1:nstreams,o11) = FLUX_MULTIPLIER * L_BOA_THTONLY_SOURCE(1:nstreams)
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     DO O1 = 1, nstokes
                        SHOM_R = ZERO ; SHOM_CR = zero
                        DO K = 1, K_REAL(NL)
                           HOM1 = NCON(K,NL) * SOLA_XPOS(I1,O1,K,NL) * T_DELT_EIGEN(K,NL)
                           HOM2 = PCON(K,NL) * SOLB_XNEG(I1,O1,K,NL)
                           SHOM_R = SHOM_R + HOM1 + HOM2
                        ENDDO
                        DO K = 1, K_COMPLEX(NL)
                           K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                           NXR1   = NCON(K1,NL) * SOLA_XPOS(I1,O1,K1,NL) - NCON(K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
                           NXR2   = NCON(K1,NL) * SOLA_XPOS(I1,O1,K2,NL) + NCON(K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
                           HOM2CR = PCON(K1,NL) * SOLB_XNEG(I1,O1,K1,NL) - PCON(K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
                           HOM1CR = NXR1 * T_DELT_EIGEN(K1,NL) - NXR2 * T_DELT_EIGEN(K2,NL)
                           SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                        ENDDO
                        BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                     ENDDO
                  ENDDO
               endif
            ELSE
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     SHOM = L_BOA_THTONLY_SOURCE(I)
                     DO LAY = NLAYERS, N, -1
                       SHOM = SHOM * T_DELT_DISORDS(I,LAY)
                     ENDDO
                     BBWF_QUAD(I,O11) = FLUX_MULTIPLIER * SHOM
                  ENDDO
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     DO O1 = 1, nstokes
                        SHOM_R = ZERO ; SHOM_CR = zero
                        DO K = 1, K_REAL(N)
                           HOM1 = NCON(K,N) * SOLA_XPOS(I1,O1,K,N)
                           HOM2 = PCON(K,N) * SOLB_XNEG(I1,O1,K,N) * T_DELT_EIGEN(K,N)
                           SHOM_R = SHOM_R + HOM1 + HOM2
                        ENDDO
                        DO K = 1, K_COMPLEX(N)
                           K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                           HOM1CR = NCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                           PXR1   = PCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - PCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                           PXR2   = PCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
                           HOM2CR = PXR1 * T_DELT_EIGEN(K1,N) - PXR2 * T_DELT_EIGEN(K2,N)
                           SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                        ENDDO
                        BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                     ENDDO
                  ENDDO
               endif

            ENDIF
         ENDIF

!  Assign fluxes

         do o1 = 1, nstokes
            SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_WEIGHTS(1:NSTREAMS))
            SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS))
            SBBWFS_FLUXES(UTA,1,O1,UPIDX) = HALF * SM
            SBBWFS_FLUXES(UTA,2,O1,UPIDX) = PI2  * SF
         enddo

!   End loop

      ENDDO
   ENDIF

!  Downwelling FLux output
!  -----------------------

   if ( DO_MVOUTPUT .and. DO_DNWELLING ) THEN
      DO UTA = 1, N_USER_LEVELS
         BBWF_QUAD = ZERO

!  Partial layer output for linearized Quadrature field

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            IF ( .not. do_thermal_transonly  ) THEN
               DO I = 1, NSTREAMS
                  DO O1 = 1, nstokes
                     SHOM_R = ZERO ; SHOM_CR = zero
                     DO K = 1, K_REAL(N)
                        HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N) * T_UTDN_EIGEN(K,UT)
!                        HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)         !  Bug 4/21/16
                        HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N) * T_UTUP_EIGEN(K,UT)
                        SHOM_R = SHOM_R + HOM1 + HOM2
                     ENDDO
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1

!  Buggy code, 4/21/16
!                        NXR1   = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
!                        NXR2   = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
!                        HOM2CR = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
!                        HOM1CR = NXR1 * T_UTDN_EIGEN(K1,UT) - NXR2 * T_UTDN_EIGEN(K2,UT)

!  replacement, not tested yet....
                        NXR_CR = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                        NXR_CI = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                        PXR_CR = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                        PXR_CI = PCON(K1,N) * SOLB_XNEG(I,O1,K2,N) + PCON(K2,N) * SOLB_XNEG(I,O1,K1,N)
                        HOM1CR =  NXR_CR * T_UTDN_EIGEN(K1,UT) - NXR_CI * T_UTDN_EIGEN(K2,UT)
                        HOM2CR =  PXR_CR * T_UTUP_EIGEN(K1,UT) - PXR_CI * T_UTUP_EIGEN(K2,UT)

                        SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                     ENDDO
                     BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                  ENDDO
               ENDDO
            ENDIF
         ENDIF

!  Level-boundary output of linearized quadrature field

         IF ( .not.PARTLAYERS_OUTFLAG(UTA) ) THEN
            NL = UTAU_LEVEL_MASK_DN(UTA) ; N = NL
            IF ( NL .NE. 0 .and. .not. do_thermal_transonly  ) THEN
               DO I = 1, NSTREAMS
                  DO O1 = 1, nstokes
                     SHOM_R = ZERO ; SHOM_CR = zero
                     DO K = 1, K_REAL(N)
                        HOM1 = NCON(K,N) * SOLA_XPOS(I,O1,K,N) * T_DELT_EIGEN(K,N)
                        HOM2 = PCON(K,N) * SOLB_XNEG(I,O1,K,N)
                        SHOM_R = SHOM_R + HOM1 + HOM2
                     ENDDO
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        NXR1   = NCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                        NXR2   = NCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + NCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                        HOM2CR = PCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                        HOM1CR = NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
                        SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
                     ENDDO
                     BBWF_QUAD(I,o1) = FLUX_MULTIPLIER * ( SHOM_R + SHOM_CR )
                  ENDDO
               ENDDO
            ENDIF
         ENDIF

!  Integrate field to get Fluxes

         do o1 = 1, nstokes
            SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_WEIGHTS(1:NSTREAMS))
            SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS))
            SBBWFS_FLUXES(UTA,1,O1,DNIDX) = HALF * SM
            SBBWFS_FLUXES(UTA,2,O1,DNIDX) = PI2  * SF
         enddo

!  End output loop

      ENDDO
   ENDIF

!  FINISH

   return
end subroutine vlidort_lbbf_jacobians_wpartials

!  End module

end module vlidort_lbbf_jacobians_m
