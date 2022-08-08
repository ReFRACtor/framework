
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
! #            VLIDORT_MediaProps                               #
! #                                                             #
! ###############################################################

MODULE VLIDORT_MediaProps_m

!  1/31/21. Version 2.8.3. NO CHANGES

!  Mark II NOTES (4/28/19)
!  =======================

!    ** Here, results are for only User angles, but TOA-UP and BOA-DN fluxes are recorded
!    ** Module is a simplified and streamlined version of the Mark I code.
!    ** By reciprocity, we have ALBMED_FLUXES(2) = TRNMED_FLUXES(1)
!    ** TRNMED_FLUXES(2) = Diffuse surface backscatter, the Sbterm in planetary problem
  
!  ORIGINAL (MARK I) NOTES (8/12/17)
!  =================================

!  Module for Computing Medium Albedos and Transmissivities
!  for Isotropic illumination sources at TOA/BOA of unit magnitude
!    -- introduced by R. Spurr 8/11/17 for LIDORT 3.7, 9/28/18 for VLIDORT 2.8
!   USER_ANGLE OUTPUT is reverse-logic from the DISORT output
!   QUAD_ANGLE OUTPUT is the same logic as from DISORT.

!  Parameter types

   USE VLIDORT_PARS_m, only : fpk

!  Back-subsitution routine

   USE vlidort_bvproblem_m, Only : BVP_BACKSUB_1

public

contains

subroutine VLIDORT_MediaProps &
     ( DO_USER_STREAMS, DO_ALBTRN_MEDIA, NSTOKES,                       & ! Input control
       NLAYERS, NSTREAMS, N_USER_STREAMS, NSTKS_NSTRMS, NSTKS_NSTRMS_2, & ! Input control numbers
       NTOTAL, N_SUBDIAG, N_SUPDIAG, QUAD_STRMWTS, DELTAU_VERT,         & ! Input quad/deltaus
       K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,           & ! Input Homog solutions
       USER_STREAMS, T_DELT_USERM, UHOM_UPDN, UHOM_UPUP,                & ! Input User solutions
       HMULT_1, HMULT_2, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,              & ! Input Multipliers, BVP
       ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,          & ! Output Main
       STATUS, MESSAGE, TRACE_1, TRACE_2 )                                ! Output Exception handling
  
!  Mark II NOTES (4/28/19)
!  =======================

!    ** Here, results are for only User angles, but TOA-UP and BOA-DN fluxes are recorded
!    ** Module is a simplified and streamlined version of the Mark I code.
!    ** By reciprocity, we have ALBMED_FLUXES(2) = TRNMED_FLUXES(1)
!    ** TRNMED_FLUXES(2) = Diffuse surface backscatter, the Sbterm in planetary problem
  
!  ORIGINAL (MARK I) NOTES (8/12/17)
!  =================================
  
!  Module first created by R. Spurr 8/12/17 (LIDORT), VLIDORT 9/30/18, 10/22/18.
!    ** Reproduces the Results from DISORT for the ibcnd = 1 condition. NSTOKES = 1
!    ** Here, results are for both User and Quad angles, in 1 call (DISORT requires 2 calls)


!  module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXSTOKES, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                                 MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAXBANDTOTAL, MAXTOTAL, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, ONE, TWO

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flag

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Control for the two illumination problems (1 = TOA, 2 = BOA)

      LOGICAL  , intent(in)  :: DO_ALBTRN_MEDIA(2)
      
!  Control integers

      INTEGER  , intent(in)  :: NSTOKES, NSTREAMS, N_USER_STREAMS, NLAYERS

!  Bookkeeping control integers

      INTEGER  , intent(in)  :: NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG

!  Quadratures and User streams

      DOUBLE PRECISION, intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, intent(in)  :: USER_STREAMS ( MAX_USER_STREAMS )

!   Input optical properties

      DOUBLE PRECISION, intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!  Homogeneous RTE solutions

      INTEGER         , INTENT (IN) :: K_REAL    ( MAXLAYERS )
      INTEGER         , INTENT (IN) :: K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, intent(in)  :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!   Transmittance factors for user-defined stream angles

      DOUBLE PRECISION, intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Homogeneous RTE solutions defined at user-defined stream angles
!   - Only the upwelling are required

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

!  Outputs
!  -------

!  Special Media-property output. -- Introduced 8/11/17 R. Spurr. Pre-initialized
!     ** Output for User-angles, also fluxes

      DOUBLE PRECISION, intent (inout) :: ALBMED_USER ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, intent (inout) :: TRNMED_USER ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, intent (inout) :: ALBMED_FLUXES ( MAXSTOKES, 2 )
      DOUBLE PRECISION, intent (inout) :: TRNMED_FLUXES ( MAXSTOKES, 2 )

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

!  other variables

      INTEGER          :: UM, I, O1, I1, IROW, N, K, KO1, K0, K1, K2
      INTEGER          :: OFFSET, STATUS_SUB

      DOUBLE PRECISION :: SHOM, H_R, H_CR, LUXR, MUXR, LUX_CR, LUX_CI, MUX_CR, MUX_CI
      DOUBLE PRECISION :: ISOFLUX, TOTAL_OD
      DOUBLE PRECISION :: ALBMED_QUAD ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: TRNMED_QUAD ( MAXSTOKES, MAXSTREAMS )

!  Start of Code
!  ============

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  Total OD for direct-source transmittance

      TOTAL_OD = SUM(DELTAU_VERT(1:NLAYERS)) 

!  Unit flux at TOA or BOA

      ISOFLUX = ONE
      
!  First Problem (Isotropic Source at TOA). ALBEDO of MEDIA
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

!  End first problem

      ENDIF

!  Second Problem (Isotropic Source at BOTTOM). TRANSMITTANCE of MEDIA
!  ===================================================================

      IF ( DO_ALBTRN_MEDIA(2) ) THEN
         
!  --set up Column for solution vector (the "B" as in AX=B)

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
             LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )       ! Output

!  -- exception handling

         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            TRACE_2 = 'Call from BVP_BACKSUB_1, BOA illumination media problem #2'
            STATUS = VLIDORT_SERIOUS ; return
         ENDIF

!  User field TOP OF THE ATMOSPHERE UPWELLING
!     ** Add direct transmittance illuminated bottom surface = ISOFLUX * EXP(-TOTAL_OD/USER_STREAMS(UM))

         IF ( DO_USER_STREAMS ) THEN

!  TOA upwelling output by recursion (source function integration)

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
                     TRNMED_USER(O1,UM) = SHOM + T_DELT_USERM(N,UM)*TRNMED_USER(O1,UM)
                  ENDDO
               ENDDO
            ENDDO
            DO UM = 1, N_USER_STREAMS
               TRNMED_USER(1,UM) = TRNMED_USER(1,UM) + ISOFLUX * EXP(-TOTAL_OD/USER_STREAMS(UM))
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
      ENDIF

!  done

      return
end subroutine VLIDORT_MediaProps

!  End module

END MODULE VLIDORT_MediaProps_m


