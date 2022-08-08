
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
! #            VLIDORT_LC_MediaProps                            #
! #                                                             #
! ###############################################################

MODULE VLIDORT_LC_MediaProps_m

!  1/31/21. Version 2.8.3. No Changes.

!  Mark II NOTES (4/28/19)
!  =======================

!    ** Here, results are for only User angles, but TOA-UP and BOA-DN fluxes are recorded
!    ** Module is a simplified and streamlined version of the Mark I code.
!    ** By reciprocity, we have ALBMED_FLUXES(2) = TRNMED_FLUXES(1)
!    ** TRNMED_FLUXES(2) = Diffuse surface backscatter, the Sbterm in planetary problem
  
!  ORIGINAL (MARK I) NOTES (8/12/17)
!  =================================

!  Module for Computing Medium Albedos and Transmissivities AND COLUMN linearization.
!  for Isotropic illumination sources at TOA/BOA of unit magnitude
!    -- first introduced by R. Spurr 8/11/17. for LIDORT 3.7, 9/28/18 for VLIDORT 2.8
!   USER_ANGLE OUTPUT is reverse-logic from the DISORT output
!   QUAD_ANGLE OUTPUT is the same logic as from DISORT.

!  Parameter types

   USE VLIDORT_PARS_m, Only : fpk

!  Back-subsitution routine

   USE vlidort_bvproblem_m, Only : BVP_BACKSUB_1

public

contains

subroutine VLIDORT_LC_MediaProps &
       ( DO_USER_STREAMS, DO_ALBTRN_MEDIA, DO_COLUMN_WFS, NSTOKES, NLAYERS,      & ! Input Control
         NSTREAMS, N_USER_STREAMS, NSTREAMS_2, N_COLUMN_WFS, NSTKS_NSTRMS,       & ! Input Control
         NSTKS_NSTRMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, QUAD_STRMWTS, DELTAU_VERT,& ! Input quad/deltaus
         K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,              & ! Input Homog solutions
         USER_STREAMS, T_DELT_USERM, UHOM_UPDN, UHOM_UPUP,                   & ! Input User solutions
         HMULT_1, HMULT_2, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                 & ! Input Multipliers, BVP
         L_DELTAU_VERT, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,            & ! Input Lin solutions
         L_T_DELT_USERM, L_UHOM_UPDN, L_UHOM_UPUP, L_HMULT_1, L_HMULT_2,     & ! Input Lin solutions
         ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,             & ! output Main
         LC_ALBMED_USER, LC_ALBMED_FLUXES, LC_TRNMED_USER, LC_TRNMED_FLUXES, & ! output Linearized
         STATUS, MESSAGE, TRACE_1, TRACE_2 )                                   ! Output Exception handling

!  Mark II NOTES (4/28/19)
!  =======================

!    ** Here, results are for only User angles, but TOA-UP and BOA-DN fluxes are recorded
!    ** Module is a simplified and streamlined version of the Mark I code.
!    ** By reciprocity, we have ALBMED_FLUXES(2) = TRNMED_FLUXES(1)
!    ** TRNMED_FLUXES(2) = Diffuse surface backscatter, the Sbterm in planetary problem
  
!  ORIGINAL (MARK I) NOTES (10/23/18)
!  =================================
  
!  Module created by R. Spurr 10/23/18. (For LIDORT, created 8/15/17).
!    ** Reproduces the Results from DISORT for the ibcnd = 1 condition. 
!    ** Here, results are for both User and Quad angles, in 1 call (DISORT requires 2 calls)
!    ** Linearizations checked by offline unit testing.

!  Revision by R. Spurr, 30 October 2018
!    ** METHOD 1: Full discrete-ordinate RT treatment for both illumination problems with non-zero albedo present
!    ** METHOD 2: Introduction of spherical albedo/transmittance (Planetary problem, nonzero albedo). Modeled After DISORT code
!    ** Simplification of the QUAD/FLUX contributions. Old code removed.
!    ** Major reduction in input arguments. Introduced the METHOD index

!  module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXSTOKES, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                                 MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAXBANDTOTAL, MAXTOTAL, &
                                 MAX_ATMOSWFS, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, ONE, TWO

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Linearization control

      LOGICAL  , intent(in)  :: DO_COLUMN_WFS

!  Control for the two illumination problems (1 = TOA, 2 = BOA)

      LOGICAL  , intent(in)  :: DO_ALBTRN_MEDIA(2)

!  Control integers

      INTEGER  , intent(in)  :: NSTOKES, NSTREAMS, N_USER_STREAMS, NLAYERS, N_COLUMN_WFS

!  Bookkeeping control integers

      INTEGER  , intent(in)  :: NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG

!  Quadratures and User streams

      DOUBLE PRECISION, intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, intent(in)  :: USER_STREAMS ( MAX_USER_STREAMS )

!   Input optical properties

      DOUBLE PRECISION, intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!   Homogeneous RTE solutions

      INTEGER         , INTENT (IN) :: K_REAL    ( MAXLAYERS )
      INTEGER         , INTENT (IN) :: K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!   Transmittance factors for user-defined stream angles

      DOUBLE PRECISION, intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Homogeneous RTE solutions defined at user-defined stream angles
!   Dont need the downwelling solutions

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multiplier variables

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Matrix, Band-matrix

      DOUBLE PRECISION, intent(in)  :: SMAT2   (MAXSTRMSTKS_2, MAXSTRMSTKS_2)
      DOUBLE PRECISION, intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT (MAXSTRMSTKS_2)

!  Linearized quantities
!  ---------------------

!   Transmittance factors for user-defined polar directions

      DOUBLE PRECISION, intent(in) :: L_T_DELT_USERM   ( MAXLAYERS,  MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized optical property inputs

      DOUBLE PRECISION, intent(in) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized Eigenvector solutions

      DOUBLE PRECISION, intent(in) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in) :: L_SOLA_XPOS    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in) :: L_SOLB_XNEG    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized user solutions

      DOUBLE PRECISION, intent(in) :: L_UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in) :: L_UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized multipliers

      DOUBLE PRECISION, intent(in) :: L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in) :: L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Outputs
!  -------

!  Special Media-property output. -- Introduced 8/11/17 R. Spurr. Pre-initialized
!     ** Output for User-angles, also fluxes

      DOUBLE PRECISION, intent (inout) :: ALBMED_USER   ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, intent (inout) :: TRNMED_USER   ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, intent (inout) :: ALBMED_FLUXES ( MAXSTOKES, 2 )
      DOUBLE PRECISION, intent (inout) :: TRNMED_FLUXES ( MAXSTOKES, 2 )

!  Linearized special media-property output

      DOUBLE PRECISION, intent (inout) :: LC_ALBMED_USER   ( MAXSTOKES, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent (inout) :: LC_TRNMED_USER   ( MAXSTOKES, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent (inout) :: LC_ALBMED_FLUXES ( MAXSTOKES, 2, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent (inout) :: LC_TRNMED_FLUXES ( MAXSTOKES, 2, MAX_ATMOSWFS )

!  Exception handling

      INTEGER      , intent(out)   :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE_1, TRACE_2

!  Local Variables
!  ===============

!  Column vector for solving BCs

      DOUBLE PRECISION :: COL2  (MAXTOTAL,1)
      DOUBLE PRECISION :: SCOL2 (MAXSTRMSTKS_2,1)

!  Solution constants of integration, and related quantities

      DOUBLE PRECISION :: LCON(MAXSTRMSTKS,MAXLAYERS)
      DOUBLE PRECISION :: MCON(MAXSTRMSTKS,MAXLAYERS)

!  Linearizations

      DOUBLE PRECISION :: COL2_WF  (MAXTOTAL,MAX_ATMOSWFS)
      DOUBLE PRECISION :: SCOL2_WF (MAXSTRMSTKS_2,MAX_ATMOSWFS)

      DOUBLE PRECISION :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Cumulative sources

      DOUBLE PRECISION :: CUMSOURCE_USER (MAXSTOKES, MAX_USER_STREAMS,0:MAXLAYERS), TOTAL_USERTRN(MAX_USER_STREAMS)

!  other variables

      INTEGER          :: UM, I, I1, O1, Q, IR, IROW, N, N1, NC, K, KO1, K0, K1, K2
      INTEGER          :: OFFSET, C0, CM, STATUS_SUB

      DOUBLE PRECISION :: SHOM, H1, H2, H3, H4, H_R, H_CR, LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI
      DOUBLE PRECISION :: TOTAL_OD, DELTA, CNEG, CPOS, T1, T1R, T1I, T2, T2R, T2I, ISOFLUX
      
      DOUBLE PRECISION :: L_TOTAL_OD(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_SOURCE(MAX_ATMOSWFS), L_SHOM_R(MAX_ATMOSWFS), L_SHOM_CR(MAX_ATMOSWFS)
      DOUBLE PRECISION :: NUXR, PUXR, LLUXR, MLUXR, L_HOM_R, L_HOM_CR, L_HOMU, L_HOMD, L_HOM
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2
      
      DOUBLE PRECISION :: ALBMED_QUAD ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: TRNMED_QUAD ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: LC_ALBMED_QUAD ( MAXSTOKES, MAXSTREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_TRNMED_QUAD ( MAXSTOKES, MAXSTREAMS, MAX_ATMOSWFS )

!  Start of Code
!  =============

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      
!  Unit flux at TOA or BOA

      ISOFLUX = ONE

!  BugFix, 8/12/19. initialize these terms, otherwise lose reproducibility!

       LC_ALBMED_USER = ZERO
       LC_TRNMED_USER = ZERO

!  Total OD for direct-source transmittance

      TOTAL_OD = SUM(DELTAU_VERT(1:NLAYERS)) ; L_TOTAL_OD = zero
      DO N = 1, NLAYERS
         DELTA = DELTAU_VERT(N)
         DO Q = 1, N_COLUMN_WFS
            L_TOTAL_OD(Q) =  L_TOTAL_OD(Q) + L_DELTAU_VERT(Q,N) * DELTA 
         ENDDO
      ENDDO

!  Total Transmittances

      IF ( DO_USER_STREAMS ) THEN
         DO UM = 1, N_USER_STREAMS
            TOTAL_USERTRN(UM) = exp ( - TOTAL_OD / USER_STREAMS(UM) )
         ENDDO
      ENDIF
      
!  First Problem (Isotropic Source at TOP). ALBEDO of MEDIA
!  ========================================================

      IF ( DO_ALBTRN_MEDIA(1) ) THEN
         
!  -- set up Column for solution vector (the "B" as in AX=B)

         IF ( NLAYERS .eq. 1 ) then
            SCOL2 = zero
            DO I = 1, NSTREAMS
               IROW = NSTOKES*(I-1)+1 ; SCOL2(IROW,1) = ISOFLUX
            ENDDO
         ELSE
            COL2 = zero
            DO I = 1, NSTREAMS
               IROW = NSTOKES*(I-1)+1 ; COL2(IROW,1)  = ISOFLUX
            ENDDO
         ENDIF

!  --Solve using backsubstitution.

         CALL BVP_BACKSUB_1 &
           ( NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,           &
             NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX, &
             BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2, SCOL2,   &
             LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )         ! Output

!  -- exception handling

         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            TRACE_2 = 'Call from BVP_BACKSUB_1, TOA illumination media problem'
            STATUS = VLIDORT_SERIOUS ; return
         ENDIF

!  -- Upwelling output at User streams.

         IF ( DO_USER_STREAMS ) THEN
            DO N = NLAYERS, 1, -1
               DO UM = 1, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                     H_R = ZERO
                     DO K = 1, K_REAL(N)
                        LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
                        MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
                        H_R = H_R + LUXR * HMULT_2(K,UM,N) + MUXR * HMULT_1(K,UM,N)
                     ENDDO
                     H_CR = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                        LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                        MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                        MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                        H_CR = H_CR + LUX_CR * HMULT_2(K1,UM,N) - LUX_CI * HMULT_2(K2,UM,N) + &
                                      MUX_CR * HMULT_1(K1,UM,N) - MUX_CI * HMULT_1(K2,UM,N)
                     ENDDO
                     SHOM = H_R + H_CR
                     ALBMED_USER(O1,UM) = SHOM + T_DELT_USERM(N,UM)*ALBMED_USER(O1,UM)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

!  -- Upwelling output at User streams.

         IF ( DO_USER_STREAMS ) THEN
            NC = 0 ; CUMSOURCE_USER(:,:,NC) = ZERO
            DO N = NLAYERS, 1, -1
               NC = NLAYERS + 1 - N
               DO UM = 1, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                     H_R = ZERO
                     DO K = 1, K_REAL(N)
                        LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
                        MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
                        H_R = H_R + LUXR * HMULT_2(K,UM,N) + MUXR * HMULT_1(K,UM,N)
                     ENDDO
                     H_CR = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                        LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                        LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                        MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                        MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                        H_CR = H_CR + LUX_CR * HMULT_2(K1,UM,N) - LUX_CI * HMULT_2(K2,UM,N) + &
                                      MUX_CR * HMULT_1(K1,UM,N) - MUX_CI * HMULT_1(K2,UM,N)
                     ENDDO
                     SHOM = H_R + H_CR
                     CUMSOURCE_USER(O1,UM,NC) = SHOM + T_DELT_USERM(N,UM) * CUMSOURCE_USER(O1,UM,NC-1)
                  ENDDO
               ENDDO
            ENDDO
            DO O1 = 1, NSTOKES
               ALBMED_USER(O1,1:N_USER_STREAMS) = CUMSOURCE_USER(O1,1:N_USER_STREAMS,NLAYERS)
            ENDDO
         ENDIF

       
!  Fluxes for this problem.
!    -- TOA flux is spherical transmittance
   
         N = 1
         DO I = 1, NSTREAMS
            I1 = I+ NSTREAMS
            DO O1 = 1, NSTOKES
               H_R = ZERO
               DO K = 1, K_REAL(N)
                  LUXR = LCON(K,N) * SOLA_XPOS(I1,O1,K,N)
                  MUXR = MCON(K,N) * SOLB_XNEG(I1,O1,K,N)
                  H_R = H_R + LUXR + MUXR * T_DELT_EIGEN(K,N)
               ENDDO
               H_CR = ZERO ; KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                  LUX_CR = LCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - LCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                  LUX_CI = LCON(K1,N) * SOLA_XPOS(I1,O1,K2,N) + LCON(K2,N) * SOLA_XPOS(I1,O1,K1,N)
                  MUX_CR = MCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - MCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
                  MUX_CI = MCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + MCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
                  H_CR = H_CR + LUX_CR + MUX_CR * T_DELT_EIGEN(K1,N) - MUX_CI * T_DELT_EIGEN(K2,N)
               ENDDO
               ALBMED_QUAD(O1,I) = H_R + H_CR
            ENDDO 
         ENDDO
         DO O1 = 1, NSTOKES
            ALBMED_FLUXES(O1,1) = TWO * DOT_PRODUCT ( ALBMED_QUAD(O1,1:NSTREAMS), QUAD_STRMWTS(1:NSTREAMS) )
         ENDDO
         
         N = NLAYERS
         DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
               H_R = ZERO
               DO K = 1, K_REAL(N)
                  LUXR = LCON(K,N) * SOLA_XPOS(I,O1,K,N)
                  MUXR = MCON(K,N) * SOLB_XNEG(I,O1,K,N)
                  H_R = H_R + LUXR * T_DELT_EIGEN(K,N) + MUXR
               ENDDO
               H_CR = ZERO ; KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                  LUX_CR = LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                  LUX_CI = LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                  MUX_CR = MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                  MUX_CI = MCON(K1,N) * SOLB_XNEG(I,O1,K2,N) + MCON(K2,N) * SOLB_XNEG(I,O1,K1,N)
                  H_CR = H_CR + MUX_CR + LUX_CR * T_DELT_EIGEN(K1,N) - LUX_CI * T_DELT_EIGEN(K2,N)
               ENDDO
               ALBMED_QUAD(O1,I) = H_R + H_CR
            ENDDO            
         ENDDO
         DO O1 = 1, NSTOKES
            ALBMED_FLUXES(O1,2) = TWO * DOT_PRODUCT ( ALBMED_QUAD(O1,1:NSTREAMS), QUAD_STRMWTS(1:NSTREAMS) )
         ENDDO
         
!  LINEARIZATION
!  -------------

!  Linearized Boundary value problem ==> 4 different conditions
         
         IF ( DO_COLUMN_WFS ) THEN
            COL2_WF = zero

!   I. First condition at TOA, top layer
!   ------------------------------------

            N = 1
            DO I = 1, NSTREAMS
               IR = NSTOKES*(I-1)
               DO O1 = 1, NSTOKES
                  IROW = IR + O1
                  DO Q = 1, N_COLUMN_WFS
!...Real
                     L_HOM_R  = ZERO
                     DO K = 1, K_REAL(N)
                        CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
                        CNEG = T_DELT_EIGEN(K,N) * L_SOLB_XNEG(I,O1,K,N,Q) + L_T_DELT_EIGEN(K,N,Q) * SOLB_XNEG(I,O1,K,N)
                        T1 = LCON(K,N) * CPOS ; T2 = MCON(K,N) * CNEG
                        L_HOM_R = L_HOM_R + T1 + T2
                     ENDDO
!...Complex
                     L_HOM_CR  = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                        K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                        T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                             L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
                        T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                             - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                             - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                        T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                             + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                             + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                        T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                        L_HOM_CR = L_HOM_CR + T1 + T2
                     ENDDO
!  ....Final contribution and end first condition loops
                     COL2_WF(IROW,Q) = - L_HOM_R - L_HOM_CR
                  ENDDO
               ENDDO
            ENDDO

!  II. Second Cndition Intermediate Levels
!  ---------------------------------------

            DO N = 1, NLAYERS - 1
               N1 = N + 1
               C0  = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
               DO I = 1, NSTREAMS_2
                 IR = NSTOKES*(I-1)
                 DO O1 = 1, NSTOKES
                   IROW = IR + O1
                   CM = C0 + IROW
                   DO Q = 1, N_COLUMN_WFS
!...Real
                     L_HOM_R  = ZERO
                     DO K = 1, K_REAL(N)
                       CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
                       CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                            L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                       T1 = LCON(K,N) * CPOS ; T2 = MCON(K,N) * CNEG
                       L_HOM_R = L_HOM_R + T1 + T2
                     ENDDO
!...Complex
                     L_HOM_CR  = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                       K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                       T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                            L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
                       T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                            - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                            + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                            - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                       T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                            + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                            + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                            + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                       T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                       L_HOM_CR = L_HOM_CR + T1 + T2
                     ENDDO
                     L_HOMU = L_HOM_R + L_HOM_CR     ! first half
!...Real
                     L_HOM_R  = ZERO
                     DO K = 1, K_REAL(N)
                       CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
                       CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                            L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                       T1 = LCON(K,N) * CPOS ; T2 = MCON(K,N) * CNEG
                       L_HOM_R = L_HOM_R + T1 + T2
                     ENDDO
!...Complex
                     L_HOM_CR  = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                       K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                       T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                            L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
                       T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                            - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                            - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                       T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                            + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                            + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                       T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                       L_HOM_CR = L_HOM_CR + T1 + T2
                     ENDDO
!  ...Final contribution and end third condition loops
                     L_HOMD = L_HOM_R + L_HOM_CR     ! Second half
                     L_HOM    = L_HOMU - L_HOMD
                     COL2_WF(CM,Q) = L_HOM
                   ENDDO
                 ENDDO
               ENDDO
            ENDDO

!  III. Third condition, Boundary layer
!  ------------------------------------

            N = NLAYERS
            C0 = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_COLUMN_WFS
!...Real
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    CPOS = T_DELT_EIGEN(K,N) * L_SOLA_XPOS(I1,O1,K,N,Q) + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                    CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                    T1 = LCON(K,N) * CPOS ; T2 = MCON(K,N) * CNEG
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO
!...Complex
                  L_HOM_CR  = ZERO ; KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                    T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) - T_DELT_EIGEN(K2,N)     * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)   - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                    T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) + T_DELT_EIGEN(K2,N)     * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)   + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = L_SOLB_XNEG(I1,O1,K1,N,Q)
                    T2I = L_SOLB_XNEG(I1,O1,K2,N,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO
                  COL2_WF(CM,Q) = - L_HOM_R - L_HOM_CR
!  ...Final contribution and end Fourth condition loops
                ENDDO
              ENDDO
            ENDDO

!  -- Copy for the one-layer case
!mick mod 8/20/2019 - do vector assignment

            IF ( NLAYERS .EQ. 1 ) THEN
              !DO I = 1, NSTKS_NSTRMS
              !  SCOL2_WF(I,1:N_COLUMN_WFS) = COL2_WF(I,1:N_COLUMN_WFS)
              !ENDDO
              SCOL2_WF(1:NSTKS_NSTRMS,1:N_COLUMN_WFS) = COL2_WF(1:NSTKS_NSTRMS,1:N_COLUMN_WFS)
            ENDIF

!  --Solve using backsubstitution.
!mick fix 1/23/2018 - trimmed 1st dim of COL2_WF --> COL2 and SCOL2_WF --> SCOL2 and added IF block
!Rob fix  10/25/18  - bug copying COL2_WF to COL2
!mick fix 8/20/2019 - trimmed 1st dim of SCOL2_WF --> SCOL2 only and added IF block

            DO Q = 1, N_COLUMN_WFS
               !COL2(1:NTOTAL,1)  = COL2_WF(1:NTOTAL,Q) ; SCOL2(:,1) = SCOL2_WF(:,Q)
               IF ( NLAYERS .EQ. 1 ) THEN
                 SCOL2(1:NSTKS_NSTRMS,1) = SCOL2_WF(1:NSTKS_NSTRMS,Q)
               ELSE
                 COL2(:,1) = COL2_WF(:,Q)
               ENDIF
               CALL BVP_BACKSUB_1 &
                 ( NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,           & ! Inputs
                   NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX, & ! Inputs
                   BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2, SCOL2,   & ! Inputs
                   NCON(:,:,Q), PCON(:,:,Q), STATUS_SUB, MESSAGE, TRACE_1 )   ! Output
               IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  TRACE_2 = 'Call from BVP_BACKSUB_1, TOA illumination media problem, Linearized ALBMED'
                  STATUS = VLIDORT_SERIOUS ; return
               ENDIF
            ENDDO

!  -- Upwelling output at User streams.

            IF ( DO_USER_STREAMS ) THEN
               DO N = NLAYERS, 1, - 1
                 NC = NLAYERS + 1 - N
                 DO UM = 1, N_USER_STREAMS
                   DO O1 = 1, NSTOKES
!...Real
                     L_SHOM_R = ZERO
                     DO K = 1, K_REAL(N)
                       DO Q = 1, N_COLUMN_WFS
                         NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
                         PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
                         H3 = NUXR * HMULT_2(K,UM,N) 
                         H4 = PUXR * HMULT_1(K,UM,N)
                         LUXR  = LCON(K,N)   *   UHOM_UPDN(UM,O1,K,N)
                         MUXR  = MCON(K,N)   *   UHOM_UPUP(UM,O1,K,N)
                         LLUXR = LCON(K,N)   * L_UHOM_UPDN(UM,O1,K,N,Q)
                         MLUXR = MCON(K,N)   * L_UHOM_UPUP(UM,O1,K,N,Q)
                         H1 = LLUXR * HMULT_2(K,UM,N) + LUXR * L_HMULT_2(K,UM,N,Q)
                         H2 = MLUXR * HMULT_1(K,UM,N) + MUXR * L_HMULT_1(K,UM,N,Q)
                         L_SHOM_R(Q) = L_SHOM_R(Q) + H1 + H2 + H3 + H4
                       ENDDO
                     ENDDO
!...Complex
                     L_SHOM_CR = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                       K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                       DO Q = 1, N_COLUMN_WFS
                         NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) - NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
                         NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) + NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
                         PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) - PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
                         PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) + PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)
                         H3 =  NUXR1 * HMULT_2(K1,UM,N) - NUXR2 * HMULT_2(K2,UM,N) 
                         H4 =  PUXR1 * HMULT_1(K1,UM,N) - PUXR2 * HMULT_1(K2,UM,N) 
                         LUXR1 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K1,N) - LCON(K2,N)   *   UHOM_UPDN(UM,O1,K2,N)
                         LUXR2 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K2,N) + LCON(K2,N)   *   UHOM_UPDN(UM,O1,K1,N)
                         MUXR1 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K1,N) - MCON(K2,N)   *   UHOM_UPUP(UM,O1,K2,N)
                         MUXR2 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K2,N) + MCON(K2,N)   *   UHOM_UPUP(UM,O1,K1,N)
                         LLUXR1 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q) - LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q)
                         LLUXR2 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q) + LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q)
                         MLUXR1 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q) - MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q)
                         MLUXR2 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q) + MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q)
                         H1 =   LLUXR1 * HMULT_2(K1,UM,N) + LUXR1 * L_HMULT_2(K1,UM,N,Q) &
                             -  LLUXR2 * HMULT_2(K2,UM,N) - LUXR2 * L_HMULT_2(K2,UM,N,Q)
                         H2 =   MLUXR1 * HMULT_1(K1,UM,N) + MUXR1 * L_HMULT_1(K1,UM,N,Q) &
                              - MLUXR2 * HMULT_1(K2,UM,N) - MUXR2 * L_HMULT_1(K2,UM,N,Q)
                         L_SHOM_CR(Q) = L_SHOM_CR(Q) + H1 + H2 + H3 + H4
                       ENDDO
                     ENDDO
!  ...Final contribution
                     DO Q = 1, N_COLUMN_WFS
                       L_SOURCE(Q) = L_SHOM_R(Q) + L_SHOM_CR(Q)
                       LC_ALBMED_USER(O1,UM,Q) =  L_SOURCE(Q) + T_DELT_USERM(N,UM) * LC_ALBMED_USER(O1,UM,Q)
                       LC_ALBMED_USER(O1,UM,Q) = LC_ALBMED_USER(O1,UM,Q) + L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_USER(O1,UM,NC-1)
                     ENDDO
!  ...Finish loops
                   ENDDO
                 ENDDO
               ENDDO
            ENDIF

!  -- Linearized Flux problem #1
!      ** Use upwelling discrete-ordinate solution at TOA (N = 1)

            N = 1
            DO I = 1, NSTREAMS
               I1 = I + NSTREAMS
               DO O1 = 1, NSTOKES
!...Real solution
                 L_SHOM_R = ZERO
                 DO K = 1, K_REAL(N)
                   DO Q = 1, N_COLUMN_WFS
                     NUXR  = NCON(K,N,Q) * SOLA_XPOS(I1,O1,K,N)
                     PUXR  = PCON(K,N,Q) * SOLB_XNEG(I1,O1,K,N) * T_DELT_EIGEN(K,N)
                     L_SHOM_R(Q) = L_SHOM_R(Q) + NUXR + PUXR
                     LUXR  = LCON(K,N)   *   SOLA_XPOS(I1,O1,K,N) ; LLUXR = LCON(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q)
                     MUXR  = MCON(K,N)   *   SOLB_XNEG(I1,O1,K,N) ; MLUXR = MCON(K,N)   * L_SOLB_XNEG(I1,O1,K,N,Q)
                     H2 = MLUXR * T_DELT_EIGEN(K,N) + MUXR * L_T_DELT_EIGEN(K,N,Q)
                     L_SHOM_R(Q) = L_SHOM_R(Q) + LLUXR + H2
                   ENDDO
                 ENDDO
!...Complex solution
                 L_SHOM_CR = ZERO ; KO1 = K_REAL(N) + 1
                 DO K = 1, K_COMPLEX(N)
                   K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                   DO Q = 1, N_COLUMN_WFS
                     NUXR1 =  NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                     PUXR1 =  PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N) - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                     PUXR2 =  PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
                     H2 =  PUXR1 * T_DELT_EIGEN(K1,N) - PUXR2 * T_DELT_EIGEN(K2,N)
                     L_SHOM_CR(Q) = L_SHOM_CR(Q) + NUXR1 + H2
                     LUXR1 =  LCON(K1,N)   *   SOLA_XPOS(I1,O1,K1,N) - LCON(K2,N)   *   SOLA_XPOS(I1,O1,K2,N)
                     LUXR2 =  LCON(K1,N)   *   SOLA_XPOS(I1,O1,K2,N) + LCON(K2,N)   *   SOLA_XPOS(I1,O1,K1,N)
                     MUXR1 =  MCON(K1,N)   *   SOLB_XNEG(I1,O1,K1,N) - MCON(K2,N)   *   SOLB_XNEG(I1,O1,K2,N)
                     MUXR2 =  MCON(K1,N)   *   SOLB_XNEG(I1,O1,K2,N) + MCON(K2,N)   *   SOLB_XNEG(I1,O1,K1,N)
                     LLUXR1 = LCON(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) - LCON(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
                     MLUXR1 = MCON(K1,N)   * L_SOLB_XNEG(I1,O1,K1,N,Q) - MCON(K2,N)   * L_SOLB_XNEG(I1,O1,K2,N,Q)
                     MLUXR2 = MCON(K1,N)   * L_SOLB_XNEG(I1,O1,K2,N,Q) + MCON(K2,N)   * L_SOLB_XNEG(I1,O1,K1,N,Q)
                     H2 =   MLUXR1 * T_DELT_EIGEN(K1,N) + MUXR1 * L_T_DELT_EIGEN(K1,N,Q) &
                          - MLUXR2 * T_DELT_EIGEN(K2,N) - MUXR2 * L_T_DELT_EIGEN(K2,N,Q)
                     L_SHOM_CR(Q) = L_SHOM_CR(Q) + LLUXR1 + H2
                   ENDDO
                 ENDDO
!  ...Final contribution
                 DO Q = 1, N_COLUMN_WFS
                   L_SOURCE(Q) = L_SHOM_R(Q) + L_SHOM_CR(Q)
                   LC_ALBMED_QUAD(O1,I,Q) =  L_SOURCE(Q)
                 ENDDO
               ENDDO
            ENDDO
!  ...FLux calculation
            DO Q = 1, N_COLUMN_WFS
               DO O1 = 1, NSTOKES
                  LC_ALBMED_FLUXES(O1,1,Q) = TWO * DOT_PRODUCT ( LC_ALBMED_QUAD(O1,1:NSTREAMS,Q), QUAD_STRMWTS(1:NSTREAMS) )
               ENDDO
            ENDDO
           
!  -- Linearized Flux problem #2
!      ** Use Downwelling discrete-ordinate solution at BOA (N = NLAYERS)

            N = NLAYERS
            DO I = 1, NSTREAMS
               DO O1 = 1, NSTOKES
!...Real solution
                 L_SHOM_R = ZERO
                 DO K = 1, K_REAL(N)
                   DO Q = 1, N_COLUMN_WFS
                     NUXR  = NCON(K,N,Q) * SOLA_XPOS(I,O1,K,N) * T_DELT_EIGEN(K,N)
                     PUXR  = PCON(K,N,Q) * SOLB_XNEG(I,O1,K,N)
                     L_SHOM_R(Q) = L_SHOM_R(Q) + NUXR + PUXR
                     LUXR  = LCON(K,N) * SOLA_XPOS(I,O1,K,N) ; LLUXR = LCON(K,N) * L_SOLA_XPOS(I,O1,K,N,Q)
                     MUXR  = MCON(K,N) * SOLB_XNEG(I,O1,K,N) ; MLUXR = MCON(K,N) * L_SOLB_XNEG(I,O1,K,N,Q)
                     H2 = LLUXR * T_DELT_EIGEN(K,N) + LUXR * L_T_DELT_EIGEN(K,N,Q)
                     L_SHOM_R(Q) = L_SHOM_R(Q) + MLUXR + H2
                   ENDDO
                 ENDDO
!...Complex solution
                 L_SHOM_CR = ZERO ; KO1 = K_REAL(N) + 1
                 DO K = 1, K_COMPLEX(N)
                   K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                   DO Q = 1, N_COLUMN_WFS
                     NUXR1 =  NCON(K1,N,Q) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N,Q) * SOLA_XPOS(I,O1,K2,N)
                     PUXR1 =  PCON(K1,N,Q) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N,Q) * SOLB_XNEG(I,O1,K2,N)
                     PUXR2 =  PCON(K1,N,Q) * SOLB_XNEG(I,O1,K2,N) + PCON(K2,N,Q) * SOLB_XNEG(I,O1,K1,N)
                     H2 =  NUXR1 * T_DELT_EIGEN(K1,N) - NUXR2 * T_DELT_EIGEN(K2,N)
                     L_SHOM_CR(Q) = L_SHOM_CR(Q) + PUXR1 + H2
                     LUXR1 =  LCON(K1,N)   *   SOLA_XPOS(I,O1,K1,N) - LCON(K2,N)   *   SOLA_XPOS(I,O1,K2,N)
                     LUXR2 =  LCON(K1,N)   *   SOLA_XPOS(I,O1,K2,N) + LCON(K2,N)   *   SOLA_XPOS(I,O1,K1,N)
                     MUXR1 =  MCON(K1,N)   *   SOLB_XNEG(I,O1,K1,N) - MCON(K2,N)   *   SOLB_XNEG(I,O1,K2,N)
                     MUXR2 =  MCON(K1,N)   *   SOLB_XNEG(I,O1,K2,N) + MCON(K2,N)   *   SOLB_XNEG(I,O1,K1,N)
                     LLUXR1 = LCON(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) - LCON(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
                     MLUXR1 = MCON(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) - MCON(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
                     MLUXR2 = MCON(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) + MCON(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
                     H2 =   LLUXR1 * T_DELT_EIGEN(K1,N) + LUXR1 * L_T_DELT_EIGEN(K1,N,Q) &
                          - LLUXR2 * T_DELT_EIGEN(K2,N) - LUXR2 * L_T_DELT_EIGEN(K2,N,Q)
                     L_SHOM_CR(Q) = L_SHOM_CR(Q) + MLUXR1 + H2
                   ENDDO
                 ENDDO
!  ...Final contribution
                 DO Q = 1, N_COLUMN_WFS
                   L_SOURCE(Q) = L_SHOM_R(Q) + L_SHOM_CR(Q)
                   LC_ALBMED_QUAD(O1,I,Q) =  L_SOURCE(Q)
                 ENDDO
               ENDDO
            ENDDO
!  ...FLux calculation
            DO Q = 1, N_COLUMN_WFS
               DO O1 = 1, NSTOKES
                  LC_ALBMED_FLUXES(O1,2,Q) = TWO * DOT_PRODUCT ( LC_ALBMED_QUAD(O1,1:NSTREAMS,Q), QUAD_STRMWTS(1:NSTREAMS) )
               ENDDO
            ENDDO
           
!  End linearization

         ENDIF
     
!  End first problem

      ENDIF


!  ===================================================================
!  Second Problem (Isotropic Source at BOTTOM). TRANSMITTANCE of MEDIA
!  ===================================================================

      IF ( DO_ALBTRN_MEDIA(2) ) THEN

!  -- set up Column for solution vector (the "B" as in AX=B)

         IF ( NLAYERS .eq. 1 ) then
            SCOL2 = zero ; OFFSET = NSTKS_NSTRMS
            DO I = 1, NSTREAMS
               IROW = OFFSET + NSTOKES*(I-1) + 1
               SCOL2(IROW, 1) = ISOFLUX
            ENDDO
         ELSE
            COL2 = zero  ; OFFSET = NTOTAL - NSTKS_NSTRMS 
            DO I = 1, NSTREAMS
              IROW = OFFSET + NSTOKES*(I-1) + 1
              COL2(IROW, 1) = ISOFLUX
            ENDDO
         ENDIF

!  --Solve using backsubstitution.

         CALL BVP_BACKSUB_1 &
            ( NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,           &
              NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX, &
              BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2, SCOL2,   &
              LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )         ! Output

!  -- exception handling

         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
           TRACE_2 = 'Call from BVP_BACKSUB_1, BOA illumination media problem #2, regular TRNMED'
           STATUS = VLIDORT_SERIOUS ; return
         ENDIF

!  -- Upwelling output at User streams.
!      ( Illumination is at the BOTTOM of the medium)
!      ( Don't forget to add the direct transmittance of the illuminated bottom surface )

         IF ( DO_USER_STREAMS ) THEN
           NC = 0 ; CUMSOURCE_USER(:,:,NC) = ZERO
           DO N = NLAYERS, 1, -1
             NC = NLAYERS + 1 - N
             DO UM = 1, N_USER_STREAMS
               DO O1 = 1, NSTOKES
                 H_R = ZERO
                 DO K = 1, K_REAL(N)
                   LUXR = LCON(K,N) * UHOM_UPDN(UM,O1,K,N)
                   MUXR = MCON(K,N) * UHOM_UPUP(UM,O1,K,N)
                   H_R = H_R + LUXR * HMULT_2(K,UM,N) + MUXR * HMULT_1(K,UM,N)
                 ENDDO
                 H_CR = ZERO ; KO1 = K_REAL(N) + 1
                 DO K = 1, K_COMPLEX(N)
                   K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                   LUX_CR = LCON(K1,N) * UHOM_UPDN(UM,O1,K1,N) - LCON(K2,N) * UHOM_UPDN(UM,O1,K2,N)
                   LUX_CI = LCON(K1,N) * UHOM_UPDN(UM,O1,K2,N) + LCON(K2,N) * UHOM_UPDN(UM,O1,K1,N)
                   MUX_CR = MCON(K1,N) * UHOM_UPUP(UM,O1,K1,N) - MCON(K2,N) * UHOM_UPUP(UM,O1,K2,N)
                   MUX_CI = MCON(K1,N) * UHOM_UPUP(UM,O1,K2,N) + MCON(K2,N) * UHOM_UPUP(UM,O1,K1,N)
                   H_CR = H_CR + LUX_CR * HMULT_2(K1,UM,N) - LUX_CI * HMULT_2(K2,UM,N) + &
                                 MUX_CR * HMULT_1(K1,UM,N) - MUX_CI * HMULT_1(K2,UM,N)
                 ENDDO
                 SHOM = H_R + H_CR
                 CUMSOURCE_USER(O1,UM,NC) = SHOM + T_DELT_USERM(N,UM) * CUMSOURCE_USER(O1,UM,NC-1)
               ENDDO
             ENDDO
           ENDDO
           DO UM = 1, N_USER_STREAMS
             DO O1 = 1, NSTOKES
               TRNMED_USER(O1,UM) = CUMSOURCE_USER(O1,UM,NLAYERS)
             ENDDO
             TRNMED_USER(1,UM) = TRNMED_USER(1,UM) + ISOFLUX * TOTAL_USERTRN(UM)
           ENDDO
         ENDIF

!  Fluxes for this problem.
!    -- BOA flux is spherical albedo

         N = 1
         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
               H_R = ZERO
               DO K = 1, K_REAL(N)
                  LUXR = LCON(K,N) * SOLA_XPOS(I1,O1,K,N)
                  MUXR = MCON(K,N) * SOLB_XNEG(I1,O1,K,N)
                  H_R = H_R + LUXR + MUXR * T_DELT_EIGEN(K,N)
               ENDDO
               H_CR = ZERO ; KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                  LUX_CR = LCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) - LCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
                  LUX_CI = LCON(K1,N) * SOLA_XPOS(I1,O1,K2,N) + LCON(K2,N) * SOLA_XPOS(I1,O1,K1,N)
                  MUX_CR = MCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) - MCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
                  MUX_CI = MCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) + MCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
                  H_CR = H_CR + LUX_CR + MUX_CR * T_DELT_EIGEN(K1,N) - MUX_CI * T_DELT_EIGEN(K2,N)
               ENDDO
               TRNMED_QUAD(O1,I) = H_R + H_CR
            ENDDO 
         ENDDO
         DO O1 = 1, NSTOKES
            TRNMED_FLUXES(O1,1) = TWO * DOT_PRODUCT ( TRNMED_QUAD(O1,1:NSTREAMS), QUAD_STRMWTS(1:NSTREAMS) )
         ENDDO

         N = NLAYERS
         DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
               H_R = ZERO
               DO K = 1, K_REAL(N)
                  LUXR = LCON(K,N) * SOLA_XPOS(I,O1,K,N)
                  MUXR = MCON(K,N) * SOLB_XNEG(I,O1,K,N)
                  H_R = H_R + LUXR * T_DELT_EIGEN(K,N) + MUXR
               ENDDO
               H_CR = ZERO ; KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                  LUX_CR = LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                  LUX_CI = LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                  MUX_CR = MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                  MUX_CI = MCON(K1,N) * SOLB_XNEG(I,O1,K2,N) + MCON(K2,N) * SOLB_XNEG(I,O1,K1,N)
                  H_CR = H_CR + MUX_CR + LUX_CR * T_DELT_EIGEN(K1,N) - LUX_CI * T_DELT_EIGEN(K2,N)
               ENDDO
               TRNMED_QUAD(O1,I) = H_R + H_CR
            ENDDO 
         ENDDO
         DO O1 = 1, NSTOKES
            TRNMED_FLUXES(O1,2) = TWO * DOT_PRODUCT ( TRNMED_QUAD(O1,1:NSTREAMS), QUAD_STRMWTS(1:NSTREAMS) )
         ENDDO
         
!  LINEARIZATION
!  -------------

!  Linearized Boundary value problem ==> 4 different conditions
         
         IF ( DO_COLUMN_WFS ) THEN
            COL2_WF = zero

!   I. First condition at TOA, top layer
!   ------------------------------------

            N = 1
            DO I = 1, NSTREAMS
               IR = NSTOKES*(I-1)
               DO O1 = 1, NSTOKES
                  IROW = IR + O1
                  DO Q = 1, N_COLUMN_WFS

!...Real
                     L_HOM_R  = ZERO
                     DO K = 1, K_REAL(N)
                       CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
                       CNEG = T_DELT_EIGEN(K,N) * L_SOLB_XNEG(I,O1,K,N,Q) + L_T_DELT_EIGEN(K,N,Q) * SOLB_XNEG(I,O1,K,N)
                       T1 = LCON(K,N) * CPOS ; T2 = MCON(K,N) * CNEG
                       L_HOM_R = L_HOM_R + T1 + T2
                     ENDDO
!...Complex
                     L_HOM_CR  = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                       K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                       T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                            L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
                       T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                            - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                            + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                            - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                       T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                            + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                            + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                            + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                       T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                       L_HOM_CR = L_HOM_CR + T1 + T2
                     ENDDO
!  ...Final contribution, and end first condition loops
                     COL2_WF(IROW,Q) = - L_HOM_R - L_HOM_CR     
                  ENDDO
               ENDDO
            ENDDO

!  II. Second Condition Intermediate Levels
!  ----------------------------------------

!  Bug 8/12/19
!    First  half LHOMU uses layer N1 solutions.  Bug was that we were using N instead of N1
!    Second half LHOMD uses layer N  solutions


            DO N = 1, NLAYERS - 1
               N1 = N + 1
               C0  = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
               DO I = 1, NSTREAMS_2
                 IR = NSTOKES*(I-1)
                 DO O1 = 1, NSTOKES
                   IROW = IR + O1
                   CM = C0 + IROW
                   DO Q = 1, N_COLUMN_WFS
!...Real (Bug fix 8/12/19. Now using N1 instead of N)
                     L_HOM_R  = ZERO
                     DO K = 1, K_REAL(N1)
                       CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
                       CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + &
                            L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
                       T1 = LCON(K,N1) * CPOS ; T2 = MCON(K,N1) * CNEG
                       L_HOM_R = L_HOM_R + T1 + T2
                     ENDDO
!...Complex (Bug fix 8/12/19. Now using N1 instead of N)
                     L_HOM_CR  = ZERO ; KO1 = K_REAL(N1) + 1
                     DO K = 1, K_COMPLEX(N1)
                       K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                       T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) - &
                            L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
                       T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                            - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                            + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1) &
                            - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
                       T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                            + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                            + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1) &
                            + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
                       T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
                       L_HOM_CR = L_HOM_CR + T1 + T2
                     ENDDO
                     L_HOMU = L_HOM_R + L_HOM_CR     ! first half
!...Real
                     L_HOM_R  = ZERO
                     DO K = 1, K_REAL(N)
                       CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
                       CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                            L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                       T1 = LCON(K,N) * CPOS ; T2 = MCON(K,N) * CNEG
                       L_HOM_R = L_HOM_R + T1 + T2
                     ENDDO
!...Complex
                     L_HOM_CR  = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                       K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                       T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                            L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
                       T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                            - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                            - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                       T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                            + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                            + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                       T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                       L_HOM_CR = L_HOM_CR + T1 + T2
                     ENDDO
!  ...Final contribution and end third condition loops
                     L_HOMD = L_HOM_R + L_HOM_CR     ! Second half
                     L_HOM    = L_HOMU - L_HOMD
                     COL2_WF(CM,Q) = L_HOM
                   ENDDO
                 ENDDO
               ENDDO
            ENDDO

!  III. Third condition, Boundary layer
!  ------------------------------------

            N = NLAYERS
            C0 = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_COLUMN_WFS
!...Real
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    CPOS = T_DELT_EIGEN(K,N) * L_SOLA_XPOS(I1,O1,K,N,Q) + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                    CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                    T1 = LCON(K,N) * CPOS ; T2 = MCON(K,N) * CNEG
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO
!...Complex
                  L_HOM_CR  = ZERO ; KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                    T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) - T_DELT_EIGEN(K2,N)     * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)   - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                    T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) + T_DELT_EIGEN(K2,N)     * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)   + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = L_SOLB_XNEG(I1,O1,K1,N,Q)
                    T2I = L_SOLB_XNEG(I1,O1,K2,N,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO
                  COL2_WF(CM,Q) = - L_HOM_R - L_HOM_CR
!  ...Final contribution and end Fourth condition loops
                ENDDO
              ENDDO
            ENDDO

!  -- Copy for the one-layer case
!mick mod 8/20/2019 - do vector assignment

            IF ( NLAYERS .EQ. 1 ) THEN
              !DO I = 1, NSTKS_NSTRMS
              !  SCOL2_WF(I,1:N_COLUMN_WFS) = COL2_WF(I,1:N_COLUMN_WFS)
              !ENDDO
              SCOL2_WF(1:NSTKS_NSTRMS,1:N_COLUMN_WFS) = COL2_WF(1:NSTKS_NSTRMS,1:N_COLUMN_WFS)
            ENDIF

!  --Solve using backsubstitution.
!mick fix 1/23/2018 - trimmed 1st dim of COL2_WF --> COL2 and SCOL2_WF --> SCOL2 and added IF block
!Rob fix  10/25/18  - bug copying COL2_WF to COL2
!mick fix 8/20/2019 - trimmed 1st dim of SCOL2_WF --> SCOL2 only and added IF block

            DO Q = 1, N_COLUMN_WFS
               !COL2(1:NTOTAL,1)  = COL2_WF(1:NTOTAL,Q) ; SCOL2(:,1) = SCOL2_WF(:,Q)
               IF ( NLAYERS .EQ. 1 ) THEN
                 SCOL2(1:NSTKS_NSTRMS,1) = SCOL2_WF(1:NSTKS_NSTRMS,Q)
               ELSE
                 COL2(:,1) = COL2_WF(:,Q)
               ENDIF
               CALL BVP_BACKSUB_1 &
                 ( NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,           & ! Inputs
                   NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX, & ! Inputs
                   BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2, SCOL2,   & ! Inputs
                   NCON(:,:,Q), PCON(:,:,Q), STATUS_SUB, MESSAGE, TRACE_1 )   ! Output
               IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  TRACE_2 = 'Call from BVP_BACKSUB_1, TOA illumination media problem, Linearized TRNMED'
                  STATUS = VLIDORT_SERIOUS ; return
               ENDIF
            ENDDO
            
!  -- Upwelling output at User streams.

            IF ( DO_USER_STREAMS ) THEN
               DO N = NLAYERS, 1, - 1
                 NC = NLAYERS + 1 - N
                 DO UM = 1, N_USER_STREAMS
                   DO O1 = 1, NSTOKES
!...Real
                     L_SHOM_R = ZERO
                     DO K = 1, K_REAL(N)
                       DO Q = 1, N_COLUMN_WFS
                         NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
                         PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
                         H1 = NUXR * HMULT_2(K,UM,N) 
                         H2 = PUXR * HMULT_1(K,UM,N)
                         L_SHOM_R(Q) = L_SHOM_R(Q) + H1 + H2
                         LUXR  = LCON(K,N)   *   UHOM_UPDN(UM,O1,K,N)
                         MUXR  = MCON(K,N)   *   UHOM_UPUP(UM,O1,K,N)
                         LLUXR = LCON(K,N)   * L_UHOM_UPDN(UM,O1,K,N,Q)
                         MLUXR = MCON(K,N)   * L_UHOM_UPUP(UM,O1,K,N,Q)
                         H1 = LLUXR * HMULT_2(K,UM,N) + LUXR * L_HMULT_2(K,UM,N,Q)
                         H2 = MLUXR * HMULT_1(K,UM,N) + MUXR * L_HMULT_1(K,UM,N,Q)
                         L_SHOM_R(Q) = L_SHOM_R(Q) + H1 + H2
                       ENDDO
                     ENDDO
!...Complex
                     L_SHOM_CR = ZERO ; KO1 = K_REAL(N) + 1
                     DO K = 1, K_COMPLEX(N)
                       K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                       DO Q = 1, N_COLUMN_WFS
                         NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) - NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
                         NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) + NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
                         PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) - PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
                         PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) + PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)
                         H1 =  NUXR1 * HMULT_2(K1,UM,N) - NUXR2 * HMULT_2(K2,UM,N) 
                         H2 =  PUXR1 * HMULT_1(K1,UM,N) - PUXR2 * HMULT_1(K2,UM,N) 
                         L_SHOM_CR(Q) = L_SHOM_CR(Q) + H1 + H2
                         LUXR1 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K1,N) - LCON(K2,N)   *   UHOM_UPDN(UM,O1,K2,N)
                         LUXR2 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K2,N) + LCON(K2,N)   *   UHOM_UPDN(UM,O1,K1,N)
                         MUXR1 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K1,N) - MCON(K2,N)   *   UHOM_UPUP(UM,O1,K2,N)
                         MUXR2 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K2,N) + MCON(K2,N)   *   UHOM_UPUP(UM,O1,K1,N)
                         LLUXR1 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q) - LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q)
                         LLUXR2 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q) + LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q)
                         MLUXR1 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q) - MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q)
                         MLUXR2 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q) + MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q)
                         H1 =   LLUXR1 * HMULT_2(K1,UM,N) + LUXR1 * L_HMULT_2(K1,UM,N,Q) &
                             -  LLUXR2 * HMULT_2(K2,UM,N) - LUXR2 * L_HMULT_2(K2,UM,N,Q)
                         H2 =   MLUXR1 * HMULT_1(K1,UM,N) + MUXR1 * L_HMULT_1(K1,UM,N,Q) &
                              - MLUXR2 * HMULT_1(K2,UM,N) - MUXR2 * L_HMULT_1(K2,UM,N,Q)
                         L_SHOM_CR(Q) = L_SHOM_CR(Q) + H1 + H2
                       ENDDO
                     ENDDO
!  ...Final contribution and end loops
                     DO Q = 1, N_COLUMN_WFS
                       L_SOURCE(Q) = L_SHOM_R(Q) + L_SHOM_CR(Q)
                       LC_TRNMED_USER(O1,UM,Q) =  L_SOURCE(Q) + T_DELT_USERM(N,UM) * LC_TRNMED_USER(O1,UM,Q)
                       LC_TRNMED_USER(O1,UM,Q) = LC_TRNMED_USER(O1,UM,Q) + L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_USER(O1,UM,NC-1)
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO
!  ...Add direct tranmsittance of surface flux.
               O1 = 1
               DO UM = 1, N_USER_STREAMS
                 DO Q = 1, N_COLUMN_WFS
                   LC_TRNMED_USER(O1,UM,Q) = LC_TRNMED_USER(O1,UM,Q) - TOTAL_USERTRN(UM) * L_TOTAL_OD(Q) / USER_STREAMS(UM)
                 ENDDO
               ENDDO
            ENDIF

!  -- Linearized Flux problem #1
!   -- Use upwelling discrete-ordinate solution at TOA (N = 1)

            N = 1
            DO I = 1, NSTREAMS
               I1 = I + NSTREAMS
               DO O1 = 1, NSTOKES
!...Real solution
                 L_SHOM_R = ZERO
                 DO K = 1, K_REAL(N)
                   DO Q = 1, N_COLUMN_WFS
                     NUXR  = NCON(K,N,Q) * SOLA_XPOS(I1,O1,K,N)
                     PUXR  = PCON(K,N,Q) * SOLB_XNEG(I1,O1,K,N) * T_DELT_EIGEN(K,N)
                     L_SHOM_R(Q) = L_SHOM_R(Q) + NUXR + PUXR
                     LUXR  = LCON(K,N)   *   SOLA_XPOS(I1,O1,K,N) ; LLUXR = LCON(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q)
                     MUXR  = MCON(K,N)   *   SOLB_XNEG(I1,O1,K,N) ; MLUXR = MCON(K,N)   * L_SOLB_XNEG(I1,O1,K,N,Q)
                     H2 = MLUXR * T_DELT_EIGEN(K,N) + MUXR * L_T_DELT_EIGEN(K,N,Q)
                     L_SHOM_R(Q) = L_SHOM_R(Q) + LLUXR + H2
                   ENDDO
                 ENDDO
!...Complex solution
                 L_SHOM_CR = ZERO ; KO1 = K_REAL(N) + 1
                 DO K = 1, K_COMPLEX(N)
                   K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                   DO Q = 1, N_COLUMN_WFS
                     NUXR1 =  NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                     PUXR1 =  PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N) - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                     PUXR2 =  PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N) + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
                     H2 =  PUXR1 * T_DELT_EIGEN(K1,N) - PUXR2 * T_DELT_EIGEN(K2,N)
                     L_SHOM_CR(Q) = L_SHOM_CR(Q) + NUXR1 + H2
                     LUXR1 =  LCON(K1,N)   *   SOLA_XPOS(I1,O1,K1,N) - LCON(K2,N)   *   SOLA_XPOS(I1,O1,K2,N)
                     LUXR2 =  LCON(K1,N)   *   SOLA_XPOS(I1,O1,K2,N) + LCON(K2,N)   *   SOLA_XPOS(I1,O1,K1,N)
                     MUXR1 =  MCON(K1,N)   *   SOLB_XNEG(I1,O1,K1,N) - MCON(K2,N)   *   SOLB_XNEG(I1,O1,K2,N)
                     MUXR2 =  MCON(K1,N)   *   SOLB_XNEG(I1,O1,K2,N) + MCON(K2,N)   *   SOLB_XNEG(I1,O1,K1,N)
                     LLUXR1 = LCON(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) - LCON(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
                     MLUXR1 = MCON(K1,N)   * L_SOLB_XNEG(I1,O1,K1,N,Q) - MCON(K2,N)   * L_SOLB_XNEG(I1,O1,K2,N,Q)
                     MLUXR2 = MCON(K1,N)   * L_SOLB_XNEG(I1,O1,K2,N,Q) + MCON(K2,N)   * L_SOLB_XNEG(I1,O1,K1,N,Q)
                     H2 =   MLUXR1 * T_DELT_EIGEN(K1,N) + MUXR1 * L_T_DELT_EIGEN(K1,N,Q) &
                          - MLUXR2 * T_DELT_EIGEN(K2,N) - MUXR2 * L_T_DELT_EIGEN(K2,N,Q)
                     L_SHOM_CR(Q) = L_SHOM_CR(Q) + LLUXR1 + H2
                   ENDDO
                 ENDDO
!  ...Final contribution and end loops
                 DO Q = 1, N_COLUMN_WFS
                   L_SOURCE(Q) = L_SHOM_R(Q) + L_SHOM_CR(Q)
                   LC_TRNMED_QUAD(O1,I,Q) =  L_SOURCE(Q)
                 ENDDO
               ENDDO
            ENDDO
!  FLux calculation
            DO Q = 1, N_COLUMN_WFS
               DO O1 = 1, NSTOKES
                  LC_TRNMED_FLUXES(O1,1,Q) = TWO * DOT_PRODUCT ( LC_TRNMED_QUAD(O1,1:NSTREAMS,Q), QUAD_STRMWTS(1:NSTREAMS) )
               ENDDO
            ENDDO

!  -- Linearized Flux problem #2
!      ** Use downwelling discrete-ordinate solution at BOA (N = NLAYERS)

            N = NLAYERS
            DO I = 1, NSTREAMS
               DO O1 = 1, NSTOKES
!...Real solution
                 L_SHOM_R = ZERO
                 DO K = 1, K_REAL(N)
                   DO Q = 1, N_COLUMN_WFS
                     NUXR  = NCON(K,N,Q) * SOLA_XPOS(I,O1,K,N) * T_DELT_EIGEN(K,N)
                     PUXR  = PCON(K,N,Q) * SOLB_XNEG(I,O1,K,N)
                     L_SHOM_R(Q) = L_SHOM_R(Q) + NUXR + PUXR
                     LUXR  = LCON(K,N) * SOLA_XPOS(I,O1,K,N) ; LLUXR = LCON(K,N) * L_SOLA_XPOS(I,O1,K,N,Q)
                     MUXR  = MCON(K,N) * SOLB_XNEG(I,O1,K,N) ; MLUXR = MCON(K,N) * L_SOLB_XNEG(I,O1,K,N,Q)
                     H2 = LLUXR * T_DELT_EIGEN(K,N) + LUXR * L_T_DELT_EIGEN(K,N,Q)
                     L_SHOM_R(Q) = L_SHOM_R(Q) + MLUXR + H2
                   ENDDO
                 ENDDO
!...Complex solution
                 L_SHOM_CR = ZERO ; KO1 = K_REAL(N) + 1
                 DO K = 1, K_COMPLEX(N)
                   K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                   DO Q = 1, N_COLUMN_WFS
                     NUXR1 =  NCON(K1,N,Q) * SOLA_XPOS(I,O1,K1,N) - NCON(K2,N,Q) * SOLA_XPOS(I,O1,K2,N)
                     PUXR1 =  PCON(K1,N,Q) * SOLB_XNEG(I,O1,K1,N) - PCON(K2,N,Q) * SOLB_XNEG(I,O1,K2,N)
                     PUXR2 =  PCON(K1,N,Q) * SOLB_XNEG(I,O1,K2,N) + PCON(K2,N,Q) * SOLB_XNEG(I,O1,K1,N)
                     H2 =  NUXR1 * T_DELT_EIGEN(K1,N) - NUXR2 * T_DELT_EIGEN(K2,N)
                     L_SHOM_CR(Q) = L_SHOM_CR(Q) + PUXR1 + H2
                     LUXR1 =  LCON(K1,N)   *   SOLA_XPOS(I,O1,K1,N) - LCON(K2,N)   *   SOLA_XPOS(I,O1,K2,N)
                     LUXR2 =  LCON(K1,N)   *   SOLA_XPOS(I,O1,K2,N) + LCON(K2,N)   *   SOLA_XPOS(I,O1,K1,N)
                     MUXR1 =  MCON(K1,N)   *   SOLB_XNEG(I,O1,K1,N) - MCON(K2,N)   *   SOLB_XNEG(I,O1,K2,N)
                     MUXR2 =  MCON(K1,N)   *   SOLB_XNEG(I,O1,K2,N) + MCON(K2,N)   *   SOLB_XNEG(I,O1,K1,N)
                     LLUXR1 = LCON(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) - LCON(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
                     MLUXR1 = MCON(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) - MCON(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
                     MLUXR2 = MCON(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) + MCON(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
                     H2 =   LLUXR1 * T_DELT_EIGEN(K1,N) + LUXR1 * L_T_DELT_EIGEN(K1,N,Q) &
                          - LLUXR2 * T_DELT_EIGEN(K2,N) - LUXR2 * L_T_DELT_EIGEN(K2,N,Q)
                     L_SHOM_CR(Q) = L_SHOM_CR(Q) + MLUXR1 + H2
                   ENDDO
                 ENDDO
!  ...Final contribution
                 DO Q = 1, N_COLUMN_WFS
                   L_SOURCE(Q) = L_SHOM_R(Q) + L_SHOM_CR(Q)
                   LC_TRNMED_QUAD(O1,I,Q) =  L_SOURCE(Q)
                 ENDDO
!  ...Finish O1 (Nstokes) and I (nstreams) loops
               ENDDO
            ENDDO
!  ...FLux calculation
            DO Q = 1, N_COLUMN_WFS
               DO O1 = 1, NSTOKES
                  LC_TRNMED_FLUXES(O1,2,Q) = TWO * DOT_PRODUCT ( LC_TRNMED_QUAD(O1,1:NSTREAMS,Q), QUAD_STRMWTS(1:NSTREAMS) )
               ENDDO
            ENDDO

!  End linearization

         ENDIF
        
!  End second problem
         
      ENDIF

!  done

      return
end subroutine VLIDORT_LC_MediaProps

!  End module

END MODULE VLIDORT_LC_MediaProps_m


