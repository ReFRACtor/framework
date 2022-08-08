
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
! #       Subroutines in this Module                            #
! #                                                             #
! # Homogeneous RTE Solutions                                   #
! #                                                             #
! #              VLIDORT_QHOM_SOLUTION                          #
! #              VLIDORT_QHOM_NORMS                             #
! #              VLIDORT_UHOM_SOLUTION                          #
! #                                                             #
! #              HMULT_MASTER (master)                          #
! #                WHOLELAYER_HMULT                             #
! #                PARTLAYER_HMULT_UP                           #
! #                PARTLAYER_HMULT_DN                           #
! #                                                             #
! # Particular Integral Solutions (Classical and Greens)        #
! #                                                             #
! #              VLIDORT_CBEAM_SOLUTION                         #
! #              VLIDORT_CUSER_SOLUTION                         #
! #              VLIDORT_GBEAM_SOLUTION                         #
! #              VLIDORT_GUSER_SOLUTION                         #
! #                                                             #
! ###############################################################

!  Notes for Version 2.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix 05/10/13  - HSINGO defined using TAYLOR_SMALL
!     Rob Fix 02/19/14  - ZETA_M = RHO_M for the degenerate case....
!                         otherwise ZETA_M = ONE / RHO_M
!                         COMPLEX HSINGO removed.

!  1/31/21. Version 2.8.3. Major Changes for Green's function option
!    ==> New Modules VLIDORT_QHOM_NORMS (needed for the Green's function solution)
!    ==> New Modules VLIDORT_GBEAM_SOLUTION, VLIDORT_GUSER_SOLUTION
!    ==> Classical Modules renamed QBEAM --> CBEAM, UBEAM --> CUSER
!    ==> TOLERANCE is now an input for the ASYMTX call
!    ==> HMULT_MASTER moved here from former module vlidort_multipliers.f90
!    ==> Private routines WHOLELAYER_HMULT, PARTLAYER_HMULT_UP, PARTLAYER_HMULT_DN

      MODULE vlidort_solutions_m

      PUBLIC :: VLIDORT_QHOM_SOLUTION,  &
                VLIDORT_QHOM_NORMS,     &
                VLIDORT_UHOM_SOLUTION,  &
                HMULT_MASTER,           &
                VLIDORT_CBEAM_SOLUTION, &
                VLIDORT_CUSER_SOLUTION, &
                VLIDORT_GBEAM_SOLUTION, &
                VLIDORT_GUSER_SOLUTION

      PRIVATE :: WHOLELAYER_HMULT, PARTLAYER_HMULT_UP, PARTLAYER_HMULT_DN

      CONTAINS

      SUBROUTINE VLIDORT_QHOM_SOLUTION ( &
        DO_SOLUTION_SAVING, GIVEN_LAYER, FOURIER, TOLERANCE,                  & ! Flag/Fourier/Layer/Tolerance
        NSTOKES, NSTREAMS, N_USER_LEVELS, NMOMENTS, NSTKS_NSTRMS,             & ! Input Control integers
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,         & ! Input partial control
        QUAD_STREAMS, QUAD_HALFWTS, DO_REAL_EIGENSOLVER, DO_LAYSCAT, & ! Quadrature and bookkeeping
        DELTAU_VERT, PARTAU_VERT, OMEGA_GREEK, PI_XQP, PI_XQM, PI_XQM_PRE,    & ! Input optical and PI matrices
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                       & ! Input transmittances (discrete Ords.)
        SAB, DAB, EIGENMAT_SAVE, REAL_KSQ, IMAG_KSQ,                          & ! Output Eigenproblem
        LEFT_EVEC, RITE_EVEC, EIGENDEGEN, EIGENMASK_R, EIGENMASK_C,           & ! Output Eigenproblem
        FWD_SUMVEC, FWD_DIFVEC, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,     & ! Output Homogeneous solutions
        K_REAL, K_COMPLEX, KEIGEN, KEIGEN_CSQ, SOLA_XPOS, SOLB_XNEG,          & ! Output Homogeneous solutions
        STATUS, MESSAGE, TRACE )                                                ! Exception handling

!  Numerical solution of Eigenproblem.

!  1/31/21. Version 2.8.3.
!   --- TOLERANCE is now an input to the ASYMTX call
!   --- KEIGEN_CSQ now has an extra "LAYERS" dimension

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_LEVELS,   &
                                 MAXSTREAMS_2, MAXSTRMSTKS, MAXEVALUES, MAXMOMENTS, VLIDORT_SERIOUS,  &
                                 VLIDORT_SUCCESS, ZERO, ONE, HALF, MINUS_ONE, MAX_TAU_QPATH

      USE VLIDORT_AUX_m , only : ASYMTX 
      USE LAPACK_TOOLS_m, only : DGEEV

      IMPLICIT NONE

!  Input flag, Fourier and given layer

      LOGICAL, INTENT (IN) ::           DO_SOLUTION_SAVING
      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER

!  1/31/21. Version 2.8.3. TOLERANCE is now an input to the ASYMTX call

      DOUBLE PRECISION, INTENT(IN) ::   TOLERANCE

!  input control numbers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS

!  Partial output control

      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  quadrature arrays

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )

!  Bookkeeping - real eigenvalues and layer scattering

      LOGICAL, INTENT (IN) ::           DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_LAYSCAT ( 0:MAXMOMENTS, MAXLAYERS )

!  Optical 

      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Pi matrices

      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  output
!  ------

!  Eigenproblem

      DOUBLE PRECISION, INTENT (INOUT) ::  SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  EIGENMAT_SAVE ( MAXEVALUES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) ::  REAL_KSQ    ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  IMAG_KSQ    ( MAXSTRMSTKS )

      DOUBLE PRECISION, INTENT (INOUT) ::  LEFT_EVEC   ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  RITE_EVEC   ( MAXSTRMSTKS, MAXSTRMSTKS )

      LOGICAL, INTENT (INOUT) ::           EIGENDEGEN  ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::           EIGENMASK_R ( MAXEVALUES )
      INTEGER, INTENT (INOUT) ::           EIGENMASK_C ( MAXEVALUES )

!  sum and difference vectors

      DOUBLE PRECISION, INTENT (INOUT) ::  FWD_SUMVEC ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) ::  FWD_DIFVEC ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  Stream transmittances

      DOUBLE PRECISION, INTENT (INOUT) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Eigensolution bookkeeping

      INTEGER, INTENT (INOUT) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (INOUT) ::           K_COMPLEX ( MAXLAYERS )

!  Eigenvalues
!    -- 1/31/21. Version 2.8.3. KEIGEN_CSQ now has an extra "LAYERS" dimension

      DOUBLE PRECISION, INTENT (INOUT) ::  KEIGEN     ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  KEIGEN_CSQ ( MAXEVALUES, MAXLAYERS )

!  Homogeneous RTE solutions

      DOUBLE PRECISION, INTENT (INOUT) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  local matrix for eigenvalue computation

      DOUBLE PRECISION :: EIGENMAT ( MAXSTRMSTKS, MAXSTRMSTKS )

!  output from Eigenpackage module ASYMTX
!    Following output now stored in commons:
!      DOUBLE PRECISION REAL_KSQ(MAXSTRMSTKS)
!      DOUBLE PRECISION RITE_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)

!  Additional diagnostic output from ASYMTX, and Work space.
!    Dimension of scratch array WK increased to 4 times the basic dimens
!    R. Spurr, RTS, December 2005: More work space needed to avoid seg-f
!    Internal dimensioning in ASYMTX is the same.

!     DOUBLE PRECISION ::  WK(16*MAXSTRMSTKS)
      DOUBLE PRECISION ::  WK(4*MAXSTRMSTKS)
      INTEGER ::           IER
      LOGICAL ::           ASYMTX_FAILURE
      CHARACTER (LEN=3) :: CN, CI

!  output for Eigenpackage DGEEV ( LAPACK)
!    Following now stored in commons...............
!      DOUBLE PRECISION REAL_KSQ(MAXSTRMSTKS),IMAG_KSQ(MAXSTRMSTKS)
!      DOUBLE PRECISION LEFT_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)
!      DOUBLE PRECISION RITE_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)

!  Additional diagnostic output and Work space from DGEEV.
!    Dimension LWORK increased from 4 to 64 times the basic dimension
!    J. Kujanpaa, FMI, August 2005: More work space for DGEEV needed.

      INTEGER, PARAMETER :: LWORK = 64*MAXSTRMSTKS
      DOUBLE PRECISION ::   WORK ( LWORK )
      INTEGER ::            DGEEV_INFO

!  optimization of Sum and Difference matrices SAB and DAB
!  ...requires intermediate storage of these two arrays.
!   Timing tests by J. Kujanpaa, FMI, August 2005.

      DOUBLE PRECISION :: &
          H2PARR ( MAXSTOKES, 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS ), &
          H2MARR ( MAXSTOKES, 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS )

!  Miscellaneous local variables

      INTEGER ::          I, J, I1, L, N, NA, UTA, UT, NMINST
      INTEGER ::          M, AA, AA1, K, K_R, K_C, K_AA, KO1, K0, K1, K2
      INTEGER ::          O1, O2, O3, IR, IROW, JC, JCOL

      DOUBLE PRECISION :: DP, DM, FAC, SC_R, XINV, HELP
      DOUBLE PRECISION :: NORM_R, TAU_DN, TAU_UP, FWD_H1
      DOUBLE PRECISION :: RKSQ, IKSQ, KMUT, KVAL
      DOUBLE PRECISION :: CARG, SARG, ARG, MOD_KVAL
      DOUBLE PRECISION :: HELP_CR, HELP_CI
      DOUBLE PRECISION :: FWD_H1_CR, FWD_H1_CI
      DOUBLE PRECISION :: SC_CR, SC_CI, MODKRT

!  Enhancement. 2/3/16, 5/9/16 using symbolic dimensions 
      DOUBLE PRECISION :: FPD(MAXSTREAMS,MAXSTOKES) , FND(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION :: FPD1(MAXSTREAMS,MAXSTOKES), FPD2(MAXSTREAMS,MAXSTOKES), &
                          FND1(MAXSTREAMS,MAXSTOKES), FND2(MAXSTREAMS,MAXSTOKES)

!  Enhancement using dynamic dimensions 7/15/15, generated errors
!      DOUBLE PRECISION :: FPD(NSTREAMS,NSTOKES) , FND(NSTREAMS,NSTOKES) 
!      DOUBLE PRECISION :: FPD1(NSTREAMS,NSTOKES), FPD2(NSTREAMS,NSTOKES), &
!                          FND1(NSTREAMS,NSTOKES), FND2(NSTREAMS,NSTOKES)

!  Variables replaced by enhancements
!      DOUBLE PRECISION :: XSQ_R, YSQ_R, H1, H2, FWD_R, FWD_I, &
!                          H1P, H1M, H2P, H2M, FWD_H2, FWD_H2_CR, FWD_H2_CI

!  original complex variables were:
!      COMPLEX*16      FWD_H1_C, FWD_H2_C, ADJ_H1_C, ADJ_H2_C
!      COMPLEX*16      HELP_C, SC_C, H1_C, H2_C

!  Code start
!  ----------

!  initialise exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Initialization

      KO1 = 0

!  initialise status

      STATUS = VLIDORT_SUCCESS

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  Flipper number

      NMINST = MIN(2,NSTOKES)

!  If there is no scattering in this layer, we can set
!  the solution vectors directly and move on.
!    Just a transmitttance case.

      IF ( DO_SOLUTION_SAVING ) THEN
       IF ( .NOT. DO_LAYSCAT(M,N) ) THEN
        K_COMPLEX(N) = 0
        K_REAL(N)    = NSTKS_NSTRMS

!  Enhancement # 1. S. Quesada 7/15/15, installed 2/3/16
        SOLA_XPOS(1:NSTREAMS,                   1:NSTOKES, 1:K_REAL(N), N) = ZERO
        SOLA_XPOS(1+NSTREAMS:NSTREAMS+NSTREAMS, 1:NSTOKES, 1:K_REAL(N), N) = ZERO
        SOLB_XNEG(1+NSTREAMS:NSTREAMS+NSTREAMS, 1:NSTOKES, 1:K_REAL(N), N) = ZERO
        SOLB_XNEG(1:NSTREAMS,                   1:NSTOKES, 1:K_REAL(N), N) = ZERO

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XINV = ONE / QUAD_STREAMS(I)
          DO O1 = 1, NMINST
            K = NSTOKES*(I-1) + O1
            KEIGEN(K,N) = XINV
            SOLA_XPOS(I,O1,K,N)  = ONE
            SOLB_XNEG(I1,O1,K,N) = ONE
          ENDDO
          DO O1 = 3, NSTOKES
            K = NSTOKES*(I-1) + O1
            KEIGEN(K,N) = XINV
            SOLA_XPOS(I,O1,K,N)  = MINUS_ONE
            SOLB_XNEG(I1,O1,K,N) = MINUS_ONE
          ENDDO
        ENDDO

        GO TO 3456
       ENDIF
      ENDIF

!  Scattering solutions
!  ====================

!  Construct Eigenmatrix
!  ---------------------

!  zero the eigensolution matrix
!   Not needed. Redundancy test by J. Kujanpaa, FMI, August 2005.
!      DO I = 1, NSTKS_NSTRMS
!        DO J = 1, NSTKS_NSTRMS
!          EIGENMAT(I,J) = ZERO
!        ENDDO
!      ENDDO

!  Develop Sum and Difference matrices, part I
!   Introduction of holding arrays by J. Kujanpaa, FMI, August 2005

!  Enhancement # 2. S. Quesada 7/15/15, installed 2/3/16
       DO J = 1, NSTREAMS
       DO O2 = 1, NSTOKES
         DO L = M, NMOMENTS
          DO O3 = 1, NSTOKES
           H2PARR(O3,L,O2,J) = sum( OMEGA_GREEK(L,N,O3,1:NSTOKES)*PI_XQP    (L,J,1:NSTOKES,O2) )
           H2MARR(O3,L,O2,J) = sum( OMEGA_GREEK(L,N,O3,1:NSTOKES)*PI_XQM_PRE(L,J,1:NSTOKES,O2) )
          ENDDO
         ENDDO
        ENDDO
      ENDDO

!  Develop Sum and Difference matrices, part II of calculation
!   Use of holding arrays by J. Kujanpaa, FMI, August 2005

!  Enhancement # 3. S. Quesada 7/15/15, installed 2/3/16
      DO I = 1, NSTREAMS
       XINV = - ONE / QUAD_STREAMS(I)
       DO J = 1, NSTREAMS
        FAC = XINV * QUAD_HALFWTS(J)
        DO O1 = 1, NSTOKES
         DO O2 = 1, NSTOKES
          DP = ZERO
          DM = ZERO
          DO L = M, NMOMENTS
           DP = DP + sum( PI_XQP(L,I,O1,1:NSTOKES)*H2PARR(1:NSTOKES,L,O2,J) )
           DM = DM + sum( PI_XQP(L,I,O1,1:NSTOKES)*H2MARR(1:NSTOKES,L,O2,J) )
          ENDDO
          SAB(I,J,O1,O2,N) = FAC * ( DP + DM )
          DAB(I,J,O1,O2,N) = FAC * ( DP - DM )
         ENDDO
        ENDDO
       ENDDO
       DO O1 = 1, NSTOKES
        SAB(I,I,O1,O1,N) = SAB(I,I,O1,O1,N) - XINV
        DAB(I,I,O1,O1,N) = DAB(I,I,O1,O1,N) - XINV
       ENDDO
      ENDDO

!  Note. In some circumstances, H2PARR = H2MARR
!        Then DP = Dm in the above code, and DAB is Diagonal.
!        Then there are pair-degenerate eigenvalues, equal to the
!        discrete ordinate directions.

!  Compute Eigenmatrix

!  Enhancement # 4. S. Quesada 7/15/15, installed 2/3/16
      DO I = 1, NSTREAMS
       IR = NSTOKES*(I-1)
       DO J = 1, NSTREAMS
        JC = NSTOKES*(J-1)
        DO O1 = 1, NSTOKES
         IROW = IR + O1
         DO O2 = 1, NSTOKES
          JCOL = JC + O2
          EIGENMAT(IROW,JCOL) = sum( DAB(I,1:NSTREAMS,O1,1:NSTOKES,N) * SAB(1:NSTREAMS,J,1:NSTOKES,O2,N) )
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!  save Eigenmatrix (original is destroyed, want to use later)
!    Looping order changed by J. Kujanpaa, FMI, August 2005

!  Enhancement # 5. S. Quesada 7/15/15, installed 2/3/16
      EIGENMAT_SAVE(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS,N) = EIGENMAT(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS)

!  Debug for testing new eigensolver
!  7 August 2006. Result --> NEW LAPACK ROUTINES NOT FASTER

!      IF ( N.EQ.1) then
!        OPEN(88,file='esolver.inp',status='unknown')
!      ENDIF
!      write(88,'(2i5,l2)')n,NSTKS_NSTRMS,DO_REAL_EIGENSOLVER(M,N)
!      DO I = 1, NSTKS_NSTRMS
!        write(88,'(1p32e20.10)')(eigenmat(i,j),j=1,NSTKS_NSTRMS)
!      enddo
!      if ( n.eq.nlayers)then
!        close(88)
!        pause'end debug'
!      ENDIF

!  Eigensolver package
!  -------------------

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Here is the DISORT ASYMTX package, as used in LIDORT scalar codes
!    This is sufficient for Real Symmetric Eigenmatrices.

!  1/31/21. Version 2.8.3. TOLERANCE is now an input for this call.

       CALL  ASYMTX &
            ( TOLERANCE, EIGENMAT, NSTKS_NSTRMS, MAXSTRMSTKS, MAXSTRMSTKS, &
              RITE_EVEC, REAL_KSQ, IER, WK, MESSAGE, ASYMTX_FAILURE )

!  Exception handling 1

       IF ( ASYMTX_FAILURE  ) THEN
         WRITE(CN,'(I3)')N
         TRACE   ='ASYMTX error in VLIDORT_QHOM_SOLUTION, Layer='//CN
         STATUS  = VLIDORT_SERIOUS
         RETURN
       ENDIF

!  Exception handling 2

       IF ( IER.GT.0 ) THEN
         WRITE(CI,'(I3)')IER
         WRITE(CN,'(I3)')N
         MESSAGE = 'eigenvalue '//CI//' has not converged'
         TRACE   = 'ASYMTX error in VLIDORT_QHOM_SOLUTION, Layer ='//CN
         STATUS  = VLIDORT_SERIOUS
         RETURN
       ENDIF

!  renormalize - this is vital for ASMTX output

!  Enhancement # 6. S. Quesada 7/15/15, installed 2/3/16
       DO AA = 1, NSTKS_NSTRMS
        NORM_R = SQRT( sum( RITE_EVEC(1:NSTKS_NSTRMS,AA)*RITE_EVEC(1:NSTKS_NSTRMS,AA) ) )
        RITE_EVEC(1:NSTKS_NSTRMS,AA) = RITE_EVEC(1:NSTKS_NSTRMS,AA) / NORM_R
       ENDDO

      ELSE

!  Complex roots : Must use DGEEV (LAPACK module)
!  DGEEV is 1.7-2.0 times slower than ASYMTX because look for complex ro
!  Also first call to DGEEV is 20 times slower.

       CALL DGEEV &
            ( 'V', 'V', NSTKS_NSTRMS, EIGENMAT, &
              MAXSTRMSTKS, REAL_KSQ, IMAG_KSQ, &
              LEFT_EVEC, MAXSTRMSTKS, RITE_EVEC, MAXSTRMSTKS, &
              WORK, LWORK, DGEEV_INFO )

!  Exception handling

       IF ( DGEEV_INFO .NE. 0 ) THEN
        WRITE(CN,'(I3)')N
        TRACE = 'DGEEV call, layer '//CN//' in VLIDORT_QHOM_SOLUTION'
        IF ( DGEEV_INFO .LT. 0 ) THEN
         WRITE(CI,'(I3)')IABS(DGEEV_INFO)
         MESSAGE='Argument # '//CI//' had an illegal value'
        ELSE
         MESSAGE='QR algorithm in DGEEV failed to compute all e-values'
        ENDIF
        STATUS = VLIDORT_SERIOUS
        RETURN
       ENDIF

      ENDIF

!  Sort real and complex eigenvalues; get masks
!  --------------------------------------------

!  Initialise counts

       K_R  = 0
       K_C  = 0
       K_AA = 0

!  all eigenvalues are real if ASYMTX is used

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
        K_R = K_R + 1
        EIGENMASK_R(K_R) = AA
       ENDDO
      ENDIF

!  Degeneracy mask

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
         KVAL = SQRT(REAL_KSQ(AA))
         EIGENDEGEN(AA,N) = .FALSE.
         DO K = 1, NSTKS_NSTRMS
           IF (K.NE.AA.AND..NOT.EIGENDEGEN(AA,N)) THEN
             KMUT= SQRT(REAL_KSQ(K))
             EIGENDEGEN(AA,N) = (ABS(KMUT-KVAL).LT.1.0D-04)
           ENDIF
         ENDDO
       ENDDO
      ENDIF

!  For the DGEEV eigensolver--
!    real eigenvalues if IMAG_KSQ = 0, otherwise complex in pairs

      IF ( .NOT.DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
        IF ( IMAG_KSQ(AA) .EQ. ZERO ) THEN
          K_R = K_R + 1
          EIGENMASK_R(K_R) = AA
        ELSE
          K_AA = K_AA + 1
          IF ( K_AA .EQ. 1 ) THEN
            K_C = K_C + 1
            EIGENMASK_C(K_C) = AA
          ELSE IF ( K_AA .EQ. 2 ) THEN
            K_AA = 0
          ENDIF
        ENDIF
       ENDDO
      ENDIF

!  save number of eigenvalues

      K_REAL(N)    = K_R
      K_COMPLEX(N) = K_C

!  Real Solutions
!  --------------

!  start loop over real eigenvalues

      DO K = 1, K_REAL(N)

!  get eigenvector entry

        AA = EIGENMASK_R(K)

!  eigenvalue and separation constants (save the eigenvalues)

        KEIGEN(K,N) = SQRT(REAL_KSQ(AA))
        SC_R    = ONE / KEIGEN(K,N)

!  set Right eigenvector norms
!    This piece of code should not be necessary as
!    (a) DGEEV norms are all = 1 upon output
!    (b) ASMTYX output has been renormalized to unity

! Enhancement # 7. S. Quesada 7/15/15, installed 2/3/16
        NORM_R = SQRT( sum( RITE_EVEC(1:NSTKS_NSTRMS,AA)*RITE_EVEC(1:NSTKS_NSTRMS,AA) ) )

!  Forward Sum-vector = normalized Right eigenvector
! Enhancement # 8. S. Quesada 7/15/15, installed 2/3/16

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          FWD_SUMVEC(I,1:NSTOKES,K) = RITE_EVEC(IR+1:IR+NSTOKES, AA)/NORM_R
        ENDDO

!  Debug (important)

!        if ( do_debug_write ) then
!         if (m.eq.1.and.n.gt.4) then
!          write(27,'(2I3,1pe20.12)') N,K,KEIGEN(K,N)
!          DO I = 1, NSTREAMS
!           DO O1 = 1, NSTOKES
!            write(27,'(4i3,1pe20.12)')N,k,i,o1,FWD_SUMVEC(I,O1,K)
!           ENDDO
!          ENDDO
!         ENDIF
!        endif

!  Find Forward difference vectors (Siewert's notation)
! Enhancement # 9. S. Quesada 7/15/15, installed 2/3/16

        DO I = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          FWD_H1 = sum( SAB (I,1:NSTREAMS,O1,1:NSTOKES,N)*FWD_SUMVEC(1:NSTREAMS,1:NSTOKES,K) )
          FWD_DIFVEC(I,O1,K) = FWD_H1 * SC_R
         ENDDO
        ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions
! Enhancement # 10. S. Quesada 7/15/15, installed 2/3/16
!mick fix 9/19/2017 - added DO loop and if block for SOLA_XPOS & SOLB_XNEG

        XINV = HALF
        FPD(1:NSTREAMS, 1:NSTOKES) = XINV * ( FWD_SUMVEC(1:NSTREAMS, 1:NSTOKES, K) + FWD_DIFVEC(1:NSTREAMS, 1:NSTOKES, K) )
        FND(1:NSTREAMS, 1:NSTOKES) = XINV * ( FWD_SUMVEC(1:NSTREAMS, 1:NSTOKES, K) - FWD_DIFVEC(1:NSTREAMS, 1:NSTOKES, K) )
        
        SOLA_XPOS(1:NSTREAMS, 1:NSTOKES, K, N)  = FPD(1:NSTREAMS, 1:NSTOKES)
        SOLB_XNEG(1:NSTREAMS, 1:NSTOKES, K, N)  = FND(1:NSTREAMS, 1:NSTOKES)

        !SOLA_XPOS(1+NSTREAMS:NSTREAMS+NSTREAMS, 1:2,       K, N)  = FND(1:NSTREAMS, 1:2)
        !SOLB_XNEG(1+NSTREAMS:NSTREAMS+NSTREAMS, 1:2,       K, N)  = FPD(1:NSTREAMS, 1:2)

        !SOLA_XPOS(1+NSTREAMS:NSTREAMS+NSTREAMS, 3:NSTOKES, K, N)  = - FND(1:NSTREAMS, 3:NSTOKES)
        !SOLB_XNEG(1+NSTREAMS:NSTREAMS+NSTREAMS, 3:NSTOKES, K, N)  = - FPD(1:NSTREAMS, 3:NSTOKES)

        DO O1 = 1, NSTOKES
          IF ( O1 .LE. 2 ) THEN
            SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS, O1, K, N)  =   FND(1:NSTREAMS, O1)
            SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS, O1, K, N)  =   FPD(1:NSTREAMS, O1)
          ELSE
            SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS, 3:NSTOKES, K, N)  = - FND(1:NSTREAMS, 3:NSTOKES)
            SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS, 3:NSTOKES, K, N)  = - FPD(1:NSTREAMS, 3:NSTOKES)
            EXIT
          ENDIF
        ENDDO

!  end loop over real eigenvalues

      ENDDO

!  Complex solutions
!  -----------------

!  Offsets

      KO1 = K_REAL(N) + 1

!  start loop over complex eigenvalues

      DO K = 1, K_COMPLEX(N)

!  Bookkeeping indices

        K0 = 2 * ( K - 1 )
        K1 = KO1 + K0
        K2 = K1  + 1

!  get eigenvector entry

        AA  = EIGENMASK_C(K)
        AA1 = AA + 1

!  set eigenvalue and separation constants for one of pair of Conjugates
!   ---> (the complex conjugate follows from this information)

        RKSQ = REAL_KSQ(AA) * REAL_KSQ(AA)
        IKSQ = IMAG_KSQ(AA) * IMAG_KSQ(AA)
        MOD_KVAL = SQRT(RKSQ + IKSQ)
        ARG = HALF * ATAN ( IMAG_KSQ(AA) / REAL_KSQ(AA) )
        CARG = COS(ARG)
        SARG = SIN(ARG)

!  original complex --------------------------------------------------
!        KEIGEN_C(K,N) = DSQRT(MOD_KVAL) * CMPLX ( CARG, SARG )
!        SC_C  = ONE_C / KEIGEN_C(K,N)
!  -------------------------------------------------------------------

!  replacement complex
!    -- 1/31/21. Version 2.8.3. KEIGEN_CSQ now has an extra "LAYERS" dimension

        KEIGEN_CSQ(K,N) = MOD_KVAL
        MODKRT          = SQRT(MOD_KVAL)
        KEIGEN(K1,N)    = MODKRT * CARG
        KEIGEN(K2,N)    = MODKRT * SARG
        SC_CR =   KEIGEN(K1,N) / MOD_KVAL
        SC_CI = - KEIGEN(K2,N) / MOD_KVAL

!  set Right eigenvector norms
!   ---> Follows convention of way that eigenvectors are fomatted
!        upon output from the LAPACK module DGEEV.

! Enhancement # 11. S. Quesada 7/15/15, installed 2/3/16
        NORM_R = SQRT( sum( RITE_EVEC(1:NSTKS_NSTRMS,AA)  * RITE_EVEC(1:NSTKS_NSTRMS,AA) + &
                            RITE_EVEC(1:NSTKS_NSTRMS,AA1) * RITE_EVEC(1:NSTKS_NSTRMS,AA1) ) )

!  Forward Sum-vector = normalized Right eigenvector

! Enhancement # 12. S. Quesada 7/15/15, installed 2/3/16
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          FWD_SUMVEC(I,1:NSTOKES,K1) = RITE_EVEC(IR+1:IR+NSTOKES,AA) /NORM_R
          FWD_SUMVEC(I,1:NSTOKES,K2) = RITE_EVEC(IR+1:IR+NSTOKES,AA1)/NORM_R
        ENDDO

!  Debug (important)

!        if ( do_debug_write ) then
!         if (m.eq.0.and.n.eq.6)then
!          write(90,'(2i3,1p2e20.10)')
!     &            m,k,KEIGEN(K1,N),KEIGEN(K2,N)
!          DO I = 1, NSTREAMS
!           DO O1 = 1, NSTOKES
!             write(90,'(4i3,1p4e20.10)')
!     &      m,k,i,o1,FWD_SUMVEC(I,O1,K1),FWD_SUMVEC(I,O1,K2)
!           ENDDO
!          ENDDO
!         ENDIF
!        endif

!  Find Forward difference vectors (Siewert's notation)

! Enhancement # 13. S. Quesada 7/15/15, installed 2/3/16
        DO O1 = 1, NSTOKES
         DO I = 1, NSTREAMS
          FWD_H1_CR = sum( SAB(I,1:NSTREAMS,O1,1:NSTOKES,N) * FWD_SUMVEC(1:NSTREAMS,1:NSTOKES,K1) )
          FWD_H1_CI = sum( SAB(I,1:NSTREAMS,O1,1:NSTOKES,N) * FWD_SUMVEC(1:NSTREAMS,1:NSTOKES,K2) )
          FWD_DIFVEC(I,O1,K1) = FWD_H1_CR*SC_CR - FWD_H1_CI*SC_CI
          FWD_DIFVEC(I,O1,K2) = FWD_H1_CR*SC_CI + FWD_H1_CI*SC_CR
         ENDDO
        ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

! Enhancement # 14. S. Quesada 7/15/15, installed 2/3/16
        XINV = HALF
        FPD1(1:NSTREAMS, 1:NSTOKES) = XINV * ( FWD_SUMVEC(1:NSTREAMS, 1:NSTOKES, K1) + FWD_DIFVEC(1:NSTREAMS, 1:NSTOKES, K1) )
        FPD2(1:NSTREAMS, 1:NSTOKES) = XINV * ( FWD_SUMVEC(1:NSTREAMS, 1:NSTOKES, K2) + FWD_DIFVEC(1:NSTREAMS, 1:NSTOKES, K2) )
        FND1(1:NSTREAMS, 1:NSTOKES) = XINV * ( FWD_SUMVEC(1:NSTREAMS, 1:NSTOKES, K1) - FWD_DIFVEC(1:NSTREAMS, 1:NSTOKES, K1) )
        FND2(1:NSTREAMS, 1:NSTOKES) = XINV * ( FWD_SUMVEC(1:NSTREAMS, 1:NSTOKES, K2) - FWD_DIFVEC(1:NSTREAMS, 1:NSTOKES, K2) )

        SOLA_XPOS(1:NSTREAMS,                   1:NSTOKES, K1, N)  = FPD1(1:NSTREAMS, 1:NSTOKES)
        SOLB_XNEG(1:NSTREAMS,                   1:NSTOKES, K1, N)  = FND1(1:NSTREAMS, 1:NSTOKES)
        SOLA_XPOS(1:NSTREAMS,                   1:NSTOKES, K2, N)  = FPD2(1:NSTREAMS, 1:NSTOKES)
        SOLB_XNEG(1:NSTREAMS,                   1:NSTOKES, K2, N)  = FND2(1:NSTREAMS, 1:NSTOKES)

        SOLA_XPOS(1+NSTREAMS:NSTREAMS+NSTREAMS, 1:2,       K1, N)  = FND1(1:NSTREAMS, 1:2)
        SOLB_XNEG(1+NSTREAMS:NSTREAMS+NSTREAMS, 1:2,       K1, N)  = FPD1(1:NSTREAMS, 1:2)
        SOLA_XPOS(1+NSTREAMS:NSTREAMS+NSTREAMS, 1:2,       K2, N)  = FND2(1:NSTREAMS, 1:2)
        SOLB_XNEG(1+NSTREAMS:NSTREAMS+NSTREAMS, 1:2,       K2, N)  = FPD2(1:NSTREAMS, 1:2)

        SOLA_XPOS(1+NSTREAMS:NSTREAMS+NSTREAMS, 3:NSTOKES, K1, N)  = - FND1(1:NSTREAMS, 3:NSTOKES)
        SOLB_XNEG(1+NSTREAMS:NSTREAMS+NSTREAMS, 3:NSTOKES, K1, N)  = - FPD1(1:NSTREAMS, 3:NSTOKES)
        SOLA_XPOS(1+NSTREAMS:NSTREAMS+NSTREAMS, 3:NSTOKES, K2, N)  = - FND2(1:NSTREAMS, 3:NSTOKES)
        SOLB_XNEG(1+NSTREAMS:NSTREAMS+NSTREAMS, 3:NSTOKES, K2, N)  = - FPD2(1:NSTREAMS, 3:NSTOKES)

!  end eigenvalue loop

      ENDDO

!  control point for avoiding eigensolver
!  ======================================

 3456 CONTINUE

!  Eigenstream transmittance factors for whole layer
!  -------------------------------------------------

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then linearized transmittances are
!  linearized discrete ordinate transmittances.

      IF ( DO_SOLUTION_SAVING .AND. &
                   .NOT.DO_LAYSCAT(M,N) ) THEN

! Enhancement # 15. S. Quesada 7/15/15, installed 2/3/16
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          T_DELT_EIGEN(IR+1:IR+NSTOKES,N) = T_DELT_DISORDS(I,N)
        ENDDO

!  Otherwise compute them as normal

      ELSE

        DO K = 1, K_REAL(N)
          HELP = KEIGEN(K,N) * DELTAU_VERT(N)
          IF ( HELP .GT. MAX_TAU_QPATH ) THEN
            T_DELT_EIGEN(K,N) = ZERO
          ELSE
            T_DELT_EIGEN(K,N) = EXP(-HELP)
          ENDIF
        ENDDO

        DO K = 1, K_COMPLEX(N)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          HELP = DELTAU_VERT(N) * KEIGEN(K1,N)
          IF ( HELP .GT. MAX_TAU_QPATH ) THEN
            T_DELT_EIGEN(K1,N) = ZERO
            T_DELT_EIGEN(K2,N) = ZERO
          ELSE
            HELP_CR = EXP ( - HELP)
            HELP_CI = DELTAU_VERT(N) * KEIGEN(K2,N)
            T_DELT_EIGEN(K1,N) =   HELP_CR * COS ( HELP_CI )
            T_DELT_EIGEN(K2,N) = - HELP_CR * SIN ( HELP_CI )
          ENDIF
        ENDDO

      ENDIF

!      if ( m.eq.0.and.n.eq.20 ) then
!       DO K = 1, min(6,K_REAL(N))
!         WRITE(*,'(2i4,1p6e24.12)')m,k,T_DELT_EIGEN(K,N)
!       ENDDO
!      ENDIF

!  Eigenstream transmittance factors for partial layers
!  ----------------------------------------------------

! Loop over user optical depths

      DO UTA = 1, N_USER_LEVELS
       IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
        UT = PARTLAYERS_OUTINDEX(UTA)
        NA = PARTLAYERS_LAYERIDX(UT)

!  If the layer of occurence = given layer then do the calculation

         IF ( NA .EQ. N ) THEN

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

          IF ( DO_SOLUTION_SAVING .AND. &
               .NOT.DO_LAYSCAT(M,N) ) THEN

! Enhancement # 16. S. Quesada 7/15/15, installed 2/3/16
           DO I = 1, NSTREAMS
             IR = NSTOKES*(I-1)
             T_UTDN_EIGEN(IR+1:IR+NSTOKES,UT) = T_DISORDS_UTDN(I,UT)
             T_UTUP_EIGEN(IR+1:IR+NSTOKES,UT) = T_DISORDS_UTUP(I,UT)
           ENDDO

!  Otherwise compute them as Eigenstream transmittance factors

          ELSE

!  partial layer optical depths

           TAU_DN = PARTAU_VERT(UT)
           TAU_UP = DELTAU_VERT(N) - TAU_DN

!  Real eigenvalues

           DO K = 1, K_REAL(N)
            HELP = KEIGEN(K,N) * TAU_DN
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTDN_EIGEN(K,UT) = ZERO
            ELSE
              T_UTDN_EIGEN(K,UT) = EXP(-HELP)
            ENDIF
            HELP = KEIGEN(K,N) * TAU_UP
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTUP_EIGEN(K,UT) = ZERO
            ELSE
              T_UTUP_EIGEN(K,UT) = EXP(-HELP)
            ENDIF
           ENDDO

!  Complex eigenvalues

           DO K = 1, K_COMPLEX(N)
            K0 = 2 * ( K - 1 )
            K1 = KO1 + K0
            K2 = K1  + 1
            HELP = TAU_DN * KEIGEN(K1,N)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTDN_EIGEN(K1,UT) = ZERO
              T_UTDN_EIGEN(K2,UT) = ZERO
            ELSE
              HELP_CR = EXP ( - HELP )
              HELP_CI = TAU_DN * KEIGEN(K2,N)
              T_UTDN_EIGEN(K1,UT) =   HELP_CR * COS ( HELP_CI )
              T_UTDN_EIGEN(K2,UT) = - HELP_CR * SIN ( HELP_CI )
            ENDIF
            HELP = TAU_UP * KEIGEN(K1,N)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTUP_EIGEN(K1,UT) = ZERO
              T_UTUP_EIGEN(K2,UT) = ZERO
            ELSE
              HELP_CR = DEXP ( - HELP )
              HELP_CI = TAU_UP * KEIGEN(K2,N)
              T_UTUP_EIGEN(K1,UT) =   HELP_CR * COS ( HELP_CI )
              T_UTUP_EIGEN(K2,UT) = - HELP_CR * SIN ( HELP_CI )
            ENDIF
           ENDDO

!  End clause to use Discrete ordinates or not

          ENDIF

!  end loop over off-grid optical depths

         ENDIF
        ENDIF
      ENDDO

!  Debug

!        if ( m.lt.3.and.n.gt.0) then
!        write(97,'(2i5)')M,N
!        DO K = 1, NSTREAMS
!          WRITE(97,'(I3,1p6e17.9)')K,KEIGEN(K,N),
!     *          (SOLA_XPOS(I,1,K,N),I=1,NSTREAMS,2)
!        ENDDO
!        endif
!c        if (n.eq.26)pause
!      if (n.eq.6) pause'end hom'

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_QHOM_SOLUTION

!

      SUBROUTINE VLIDORT_QHOM_NORMS &
         ( FOURIER, NSTOKES, NSTREAMS, NLAYERS, QUAD_STRMWTS, & ! Input numbers
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! Input solutions
           NORM_SAVED )                                         ! Output

!  Eigenproblem, solution norms for Green's function

!  1/31/21. Version 2.8.3. New for this version (programmed after the LIDORT routine of the same name)

!  module,  dimensions and numbers

      USE VLIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXEVALUES, MAXSTOKES, MAXLAYERS, ZERO

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  FOURIER component

      INTEGER  , intent(in)  :: FOURIER

!  Number of stokes elements, discrete-ordinate streams streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NSTOKES
      INTEGER  , intent(in)  :: NLAYERS

!  Quadrature

      DOUBLE PRECISION, intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Eigensolution obookkeeping

      INTEGER, INTENT (IN) ::    K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::    K_COMPLEX ( MAXLAYERS )

!  Homogeneous RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  outputs
!  -------

!  Saved quantities for the Green function solution (Normalizations)
!    --INDICES REVERSED FROM LIDORT values

      DOUBLE PRECISION, intent(inout) :: NORM_SAVED ( MAXEVALUES, MAXLAYERS)

!  Local variables
!  ---------------

!  Miscellaneous local variables

      INTEGER          :: N, K, O1, J, KO1, K0, K1, k2
      DOUBLE PRECISION :: T1, T2, NORM
      DOUBLE PRECISION :: T1_CR, T1_CI, T2_CR, T2_CI, NORM_CR, NORM_CI
    
!  For all layers, save the norms

      DO N = 1, NLAYERS
         KO1 = K_REAL(N) + 1

!  Loop over real eigensolutions

         DO K = 1, K_REAL(N)
            NORM = ZERO
            DO J = 1, NSTREAMS
               T1 = zero ; T2 = zero
               DO O1 = 1, NSTOKES
                  T1 = T1 + SOLA_XPOS(J,O1,K,N) * SOLA_XPOS(J,O1,K,N)
                  T2 = T2 + SOLB_XNEG(J,O1,K,N) * SOLB_XNEG(J,O1,K,N)         
               ENDDO
               NORM = NORM + QUAD_STRMWTS(J) * ( T1 - T2 )
            ENDDO
            NORM_SAVED(K,N) = NORM
         ENDDO

!  start loop over complex eigenvalues
!     SUGGESTED CODE with the Adjoint solutions

         DO K = 1, K_COMPLEX(N)
            K0 = 2 * ( K - 1 ) ; K1 = KO1 + K0 ; K2 = K1  + 1
            NORM_CR = ZERO ; NORM_CI = zero
            DO J = 1, NSTREAMS
               T1_CR = zero ; T2_CR = zero
               T1_CI = zero ; T2_CI = zero
               DO O1 = 1, NSTOKES
!                  T1_CR = T1_CR + ADJA_XPOS(J,O1,K1,N) * ADJA_XPOS(J,O1,K1,N) &
!                                - ADJA_XPOS(J,O1,K2,N) * ADJA_XPOS(J,O1,K2,N)
!                  T2_CR = T2_CR + ADJB_XNEG(J,O1,K1,N) * ADJB_XNEG(J,O1,K1,N) &
!                                - ADJB_XNEG(J,O1,K2,N) * ADJB_XNEG(J,O1,K2,N)
!                  T1_CI = T1_CI + ADJA_XPOS(J,O1,K1,N) * ADJA_XPOS(J,O1,K2,N) &
!                                + ADJA_XPOS(J,O1,K2,N) * ADJA_XPOS(J,O1,K1,N)
!                  T2_CI = T2_CI + ADJB_XNEG(J,O1,K1,N) * ADJB_XNEG(J,O1,K1,N) &
!                                + ADJB_XNEG(J,O1,K2,N) * ADJB_XNEG(J,O1,K2,N)
               ENDDO
               NORM_CR = NORM_CR + QUAD_STRMWTS(J) * ( T1_CR - T2_CR )
               NORM_CI = NORM_CI + QUAD_STRMWTS(J) * ( T1_CI - T2_CI )
            ENDDO
            NORM_SAVED(K1,N) = NORM_CR
            NORM_SAVED(K2,N) = NORM_CI
         ENDDO

!  End Layer Loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE VLIDORT_QHOM_NORMS

!

SUBROUTINE VLIDORT_UHOM_SOLUTION ( &
        DO_UPWELLING, DO_DNWELLING, GIVEN_LAYER, FOURIER,            & ! Input flags/indices
        NSTOKES, NSTREAMS, N_USER_STREAMS, NMOMENTS, DO_LAYSCAT,     & ! Input numbers/flag
        QUAD_HALFWTS, USER_SECANTS, OMEGA_GREEK,                     & ! Input quadratures + optical
        PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE, PI_XUM_POST, PI_XUP_PRE, & ! Input PI Matrices
        K_REAL, K_COMPLEX, KEIGEN, KEIGEN_CSQ, SOLA_XPOS, SOLB_XNEG, & ! Input Homog. RTE solutions
        UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                  & ! Output user Homog
        HSINGO, HELPSTOKES, ZETA_M, ZETA_P )                           ! Output user Homog

!  Small-number Taylor series expansions, added 30 October 2007
!  @@@ Rob Fix 10 May 13 - Redefine HSINGO, supercedes above comment

!  1/31/21. Version 2.8.3. 
!    -- KEIGEN_CSQ now has an extra "LAYERS" dimension
!    -- UM loops always start with 1 now, not LOCAL_UM_START (argument dropped)

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, MAXMOMENTS, &
                                 MAXSTREAMS_2, MAXEVALUES, ZERO, ONE, TWO, TAYLOR_SMALL                   

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING

!  removed, 7/5/16
!      LOGICAL, INTENT (IN) ::           DO_DEBUG_WRITE
!      LOGICAL, INTENT (IN) ::           DO_FDTEST

!  Input indices/numbers

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS

!  bookkeeping. 1/31/21. Version 2.8.3. Drop LOCAL_UM_START

      LOGICAL, INTENT (IN) ::           DO_LAYSCAT ( 0:MAXMOMENTS, MAXLAYERS )

!  Quadrature and Optical

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS  ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK   ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS  ( MAX_USER_STREAMS ) removed 7/5/16

!  PI matrices

      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM_POST ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP_PRE  ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Homogeneous solutions
!   -- 1/31/21. Version 2.8.3. KEIGEN_CSQ now has an extra "LAYERS" dimension

      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  KEIGEN     ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  KEIGEN_CSQ ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  OUTPUT
!  ======

!  User-angle homogeneous solutions

      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Auxiliary quantities for the User solutions

      LOGICAL, INTENT (INOUT) ::          HSINGO ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: HELPSTOKES ( 0:MAXMOMENTS, MAXEVALUES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Local variables
!  ---------------

      INTEGER ::          UM, J, L, N, M, K, O1, KO1, K0, K1, K2

      DOUBLE PRECISION :: GAUX_R  ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: GAUX_CR ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: GAUX_CI ( 0:MAXMOMENTS, MAXSTOKES )

!  for real solutions

      DOUBLE PRECISION :: SM, SP, SPS, A5
      DOUBLE PRECISION :: RHO_P, RHO_M, SECMUI

!  for complex solutions

      DOUBLE PRECISION :: RHO_P_CR, RHO_M_CR, MODULUS_P, MODULUS_M
      DOUBLE PRECISION :: A_USER, A_HELP, B_USER, B_HELP
!      DOUBLE PRECISION :: SN, SNS, HPC, HMC, SNCR, SNCI,  SNSCR, SNSCI, SGCR, SGCI, SUMG
      DOUBLE PRECISION :: SMCR, SPCR, SPSCR, SMCI, SPCI, SPSCI

!  Notes for Version 2.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix 05/10/13  - HSINGO defined using TAYLOR_SMALL
!     Rob Fix 02/19/14  - ZETA_M = RHO_M for the degenerate case....
!                         otherwise ZETA_M = ONE / RHO_M
!                         COMPLEX HSINGO removed.

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  Initialization

      K   = 0

!  Zeta constants
!  ==============

!  Real Zeta constants (always required)

      DO UM = 1, N_USER_STREAMS
        SECMUI = USER_SECANTS(UM)
        DO K = 1, K_REAL(N)
          RHO_P = SECMUI + KEIGEN(K,N)
          RHO_M = SECMUI - KEIGEN(K,N)
          ZETA_P(K,UM,N) = ONE / RHO_P
          HSINGO(K,UM,N) = ( ABS(RHO_M) .LT. TAYLOR_SMALL )
          IF ( HSINGO(K,UM,N) ) then
            ZETA_M(K,UM,N) = RHO_M
          ELSE
            ZETA_M(K,UM,N) = ONE / RHO_M
          ENDIF
        ENDDO
      ENDDO

!  Complex Zeta constants
!   -- 1/31/21. Version 2.8.3. KEIGEN_CSQ now has an extra "LAYERS" dimension

      IF ( K_COMPLEX(N) .GT. 0 ) THEN
        KO1 = K_REAL(N) + 1
        DO UM = 1, N_USER_STREAMS
          A_USER = USER_SECANTS(UM) * USER_SECANTS(UM)
          B_USER = TWO * USER_SECANTS(UM)
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            RHO_P_CR = USER_SECANTS(UM) + KEIGEN(K1,N)
            RHO_M_CR = USER_SECANTS(UM) - KEIGEN(K1,N)
            A_HELP = A_USER + KEIGEN_CSQ(K,N)            ! Add layer index, Version 2.8.3.
            B_HELP = B_USER * KEIGEN(K1,N)
            MODULUS_P = A_HELP + B_HELP
            MODULUS_M = A_HELP - B_HELP
            ZETA_P(K1,UM,N) = RHO_P_CR / MODULUS_P
            ZETA_P(K2,UM,N) = -KEIGEN(K2,N) / MODULUS_P
            ZETA_M(K1,UM,N) = RHO_M_CR / MODULUS_M
            ZETA_M(K2,UM,N) = KEIGEN(K2,N) / MODULUS_M
          ENDDO
        ENDDO
      ENDIF

!  If there is no scattering avoid these solutions
!  Warning - Solutions are zeroed in the scalar code before returning
!  @@@ Rob Fix. Zero them here too.

!mick fix 2/1/2011 - modified IF to initialize solutions
!      IF ( .NOT. DO_LAYSCAT(M,N) ) RETURN

      IF ( .NOT. DO_LAYSCAT(M,N) ) THEN
       IF ( DO_DNWELLING ) THEN
          UHOM_DNDN(:,:,:,N) = ZERO
          UHOM_DNUP(:,:,:,N) = ZERO
       ENDIF
       IF ( DO_UPWELLING ) THEN
          UHOM_UPDN(:,:,:,N) = ZERO
          UHOM_UPUP(:,:,:,N) = ZERO
       ENDIF
       HELPSTOKES(M,:,:)  = ZERO
       RETURN
      ENDIF

!  User defined solutions (Real)
!  -----------------------------

      DO K = 1, K_REAL(N)

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors.
!    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec

! Enhancement # 17. S. Quesada 7/15/15, installed 2/3/16

       DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
         SPS = ZERO
         DO J = 1, NSTREAMS
          A5 = QUAD_HALFWTS(J)
          SP = sum( PI_XQP    (L,J,O1,1:NSTOKES)*A5*SOLA_XPOS(J,1:NSTOKES,K,N) )
          SM = sum( PI_XQM_PRE(L,J,O1,1:NSTOKES)*A5*SOLB_XNEG(J,1:NSTOKES,K,N) )
          SPS = SPS + SP + SM
         ENDDO
         HELPSTOKES(L,K,O1) = SPS
        ENDDO
        DO O1 = 1, NSTOKES
         GAUX_R(L,O1) = sum( OMEGA_GREEK(L,N,O1,1:NSTOKES) * HELPSTOKES(L,K,1:NSTOKES) )
        ENDDO
       ENDDO

!  Now sum over all harmonic contributions (downwelling)

       IF ( DO_DNWELLING ) THEN
! Enhancement # 18. S. Quesada 7/15/15, installed 2/3/16
        DO O1 = 1, NSTOKES
         DO UM = 1, N_USER_STREAMS
          UHOM_DNDN(UM,O1,K,N) = sum( PI_XUP     (M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_R(M:NMOMENTS,1:NSTOKES) )
          UHOM_DNUP(UM,O1,K,N) = sum( PI_XUM_POST(M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_R(M:NMOMENTS,1:NSTOKES) )
         ENDDO
        ENDDO
       ENDIF

!  Now sum over all harmonic contributions (Upwelling)

       IF ( DO_UPWELLING ) THEN
! Enhancement # 19. S. Quesada 7/15/15, installed 2/3/16
        DO O1 = 1, NSTOKES
         DO UM = 1, N_USER_STREAMS
          UHOM_UPDN(UM,O1,K,N) = sum( PI_XUM    (M:NMOMENTS,UM,O1,1:NSTOKES)  * GAUX_R(M:NMOMENTS,1:NSTOKES) )
          UHOM_UPUP(UM,O1,K,N) = sum( PI_XUP_PRE(M:NMOMENTS,UM,O1,1:NSTOKES)  * GAUX_R(M:NMOMENTS,1:NSTOKES) )
         ENDDO
        ENDDO
       ENDIF

!  Debug (important)
!        IF ( DO_DEBUG_WRITE .AND. DO_FDTEST) THEN
!          DO UM = 1, N_USER_STREAMS ;  DO O1 = 1, NSTOKES
!            WRITE(94,'(4I3,1P4E20.10)')M,K,UM,O1, UHOM_UPUP(UM,O1,K,N), UHOM_UPDN(UM,O1,K,N), &
!               UHOM_DNUP(UM,O1,K,N), UHOM_DNDN(UM,O1,K,N)
!          ENDDO ; ENDDO
!        ENDIF

!  end loop over real eigensolutions

      ENDDO

!  User defined solutions (Complex)
!  -------------------------------

!  offset

      KO1 = K_REAL(N) + 1

      DO K = 1, K_COMPLEX(N)

       K0 = 2 * K - 2
       K1 = KO1 + K0
       K2 = K1  + 1

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors.
!    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec

! Enhancement # 20. S. Quesada 7/15/15, installed 2/3/16

       DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
         SPSCR = ZERO
         SPSCI = ZERO
         DO J = 1, NSTREAMS
          SPCR = sum( PI_XQP(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J) * SOLA_XPOS(J,1:NSTOKES,K1,N) )
          SPCI = sum( PI_XQP(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J) * SOLA_XPOS(J,1:NSTOKES,K2,N) )
          SMCR = sum( PI_XQM_PRE(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J) * SOLB_XNEG(J,1:NSTOKES,K1,N) )
          SMCI = sum( PI_XQM_PRE(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J) * SOLB_XNEG(J,1:NSTOKES,K2,N) )
          SPSCR = SPSCR + SPCR + SMCR
          SPSCI = SPSCI + SPCI + SMCI
         ENDDO
         HELPSTOKES(L,K1,O1) = SPSCR
         HELPSTOKES(L,K2,O1) = SPSCI
        ENDDO
        DO O1 = 1, NSTOKES
         GAUX_CR(L,O1) = sum( OMEGA_GREEK(L,N,O1,1:NSTOKES) * HELPSTOKES(L,K1,1:NSTOKES) )
         GAUX_CI(L,O1) = sum( OMEGA_GREEK(L,N,O1,1:NSTOKES) * HELPSTOKES(L,K2,1:NSTOKES) )
        ENDDO
       ENDDO
   
!  Now sum over all harmonic contributions (downwelling)

       IF ( DO_DNWELLING ) THEN
! Enhancement # 21. S. Quesada 7/15/15, installed 2/3/16
        DO O1 = 1, NSTOKES
         DO UM = 1, N_USER_STREAMS
          UHOM_DNDN(UM,O1,K1,N) = sum( PI_XUP     (M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_CR(M:NMOMENTS,1:NSTOKES) )
          UHOM_DNUP(UM,O1,K1,N) = sum( PI_XUM_POST(M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_CR(M:NMOMENTS,1:NSTOKES) )
          UHOM_DNDN(UM,O1,K2,N) = sum( PI_XUP     (M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_CI(M:NMOMENTS,1:NSTOKES) )
          UHOM_DNUP(UM,O1,K2,N) = sum( PI_XUM_POST(M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_CI(M:NMOMENTS,1:NSTOKES) )
         ENDDO
        ENDDO
       ENDIF

!  Now sum over all harmonic contributions (upwelling)

       IF ( DO_UPWELLING ) THEN
! Enhancement # 22. S. Quesada 7/15/15, installed 2/3/16
         DO O1 = 1, NSTOKES
         DO UM = 1, N_USER_STREAMS
          UHOM_UPDN(UM,O1,K1,N) = sum( PI_XUM    (M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_CR(M:NMOMENTS,1:NSTOKES) )
          UHOM_UPUP(UM,O1,K1,N) = sum( PI_XUP_PRE(M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_CR(M:NMOMENTS,1:NSTOKES) )
          UHOM_UPDN(UM,O1,K2,N) = sum( PI_XUM    (M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_CI(M:NMOMENTS,1:NSTOKES) )
          UHOM_UPUP(UM,O1,K2,N) = sum( PI_XUP_PRE(M:NMOMENTS,UM,O1,1:NSTOKES) * GAUX_CI(M:NMOMENTS,1:NSTOKES) )
         ENDDO
        ENDDO
      ENDIF

!  Debug (important)
!        if ( do_debug_write.and.do_fdtest) then
!          DO UM = 1, N_USER_STREAMS ; DO O1 = 1, NSTOKES
!            write(96,'(4i3,1p4e20.10)')m,k,um,o1,UHOM_UPUP(UM,O1,K2,N),UHOM_UPDN(UM,O1,K2,N),&
!                  UHOM_DNUP(UM,O1,K2,N), UHOM_DNDN(UM,O1,K2,N)
!          enddo ; enddo
!        endif

!  end loop over complex eigensolutions

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_UHOM_SOLUTION

!

      SUBROUTINE HMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING,                                   & ! flags
        NLAYERS, N_USER_STREAMS, N_USER_LEVELS, TAYLOR_ORDER,         & ! Basic numbers
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                       & ! whole-layer control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! partial-layer control
        USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                       & ! secants, optical
        T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                     & ! User-stream transmittances
        T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                     & ! eigenstreamTransmittances
        K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,                    & ! RTE Eigen solution
        HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD ) ! OUTPUT Multipliers

!  1/31/21. Version 2.8.3. 
!    -- UM loops always start with 1 now, not LOCAL_UM_START (argument dropped)

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS,  MAX_USER_STREAMS, MAX_USER_LEVELS, MAXEVALUES

      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING

!  Numbers

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           TAYLOR_ORDER      ! 2p7 Taylor-order argument added

!  whole-layer control

      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

!  Partial-layer control

      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  secants and optical

      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT  ( MAX_PARTLAYERS )

!  User stream transm.

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenstream trans.

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Eigensolution variables

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           HSINGO ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  OUTPUT
!  ======

!  whole-layer multipliers

      DOUBLE PRECISION, INTENT (OUT) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Partial-layer multipliers

      DOUBLE PRECISION, INTENT (OUT) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Local variables
!  ---------------

      INTEGER          :: N, UT, UTA

!  whole layer multipliers. Version 2p7 Taylor-order argument added
!    -- 1/31/21. Version 2.8.3.  UM loops always start with 1 now, not LOCAL_UM_START (argument dropped)

      CALL WHOLELAYER_HMULT ( &
        NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,      & ! Numbers, 2p7 Taylor-order argument added
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,     & ! Layer control
        USER_SECANTS, DELTAU_VERT,                  & ! streams, optical
        T_DELT_USERM, T_DELT_EIGEN,                 & ! Transmittances
        K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,  & ! Eigensolution stuff
        HMULT_1, HMULT_2 )

!  partial layer multipliers
!   @@@ Rob Fix 10 May 13 - DELTAU_VERT no longer required from Downwelling
!    -- 1/31/21. Version 2.8.3.  UM loops always start with 1 now, not LOCAL_UM_START (argument dropped)

      DO UTA = 1, N_USER_LEVELS

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

          IF ( DO_UPWELLING ) THEN
            CALL PARTLAYER_HMULT_UP ( &
              N, UT, N_USER_STREAMS, TAYLOR_ORDER,                    & ! Numbers, 2p7 Taylor-order argument added
              USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                 & ! streams, optical
              T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTUP_USERM, & ! Transmittances
              K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,              & ! Eigensolution stuff
              UT_HMULT_UU, UT_HMULT_UD )
          END IF

          IF ( DO_DNWELLING ) THEN
            CALL PARTLAYER_HMULT_DN ( &
              N, UT, N_USER_STREAMS, TAYLOR_ORDER, USER_SECANTS, PARTAU_VERT, & ! streams, optical
              T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTDN_USERM,         & ! Transmittances
              K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,                      & ! Eigensolution stuff
              UT_HMULT_DD, UT_HMULT_DU )
          END IF

        ENDIF

      ENDDO

!  Finish

!  @@@ Rob Fix 10/13. This is a useful place to pause when debugging
!      pause'Finish HMULT master'

      RETURN
      END SUBROUTINE HMULT_MASTER

!

      SUBROUTINE WHOLELAYER_HMULT ( &
        NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,      & ! Numbers, 2p7 Taylor-order argument added
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,     & ! Layer control
        USER_SECANTS, DELTAU_VERT,                  & ! streams, optical
        T_DELT_USERM, T_DELT_EIGEN,                 & ! Transmittances
        K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,  & ! Eigensolution stuff
        HMULT_1, HMULT_2 )

!  1/31/21. Version 2.8.3.
!    -- UM loops always start with 1 now, not LOCAL_UM_START (argument dropped)

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_USER_STREAMS, MAX_USER_LEVELS, MAXEVALUES, ONE

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_1

      IMPLICIT NONE

!  Numbers

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           TAYLOR_ORDER      ! 2p7 Taylor-order argument added

!  whole-layer control

      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

!  secants and optical

      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT  ( MAXLAYERS )


!  User stream + eigenstream transm.

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )


!  Eigensolution variables

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           HSINGO ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  OUTPUT
!  ======

!  whole-layer multipliers

      DOUBLE PRECISION, INTENT (OUT) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Local variables
!  ---------------

      INTEGER          :: UM, K, N, KO1, K0, K1, K2
      DOUBLE PRECISION :: UDEL, SM, EPS
      DOUBLE PRECISION :: ZDEL_R, ZUDEL_R, THETA_1_R, THETA_2_R
      DOUBLE PRECISION :: ZDEL_CR,    ZDEL_CI
      DOUBLE PRECISION :: THETA_1_CR, THETA_1_CI
      DOUBLE PRECISION :: THETA_2_CR, THETA_2_CI

!  Start loops over layers and user-streams
!   Only done if layers are flagged

!  Version 2p7 notes
!     Rob Fix 05/10/13  - Introduce Taylor series routines (first version). Included Complex series
!     Rob Fix 02/19/14  - Small numbers analysis finalized using Taylor_Order parameter, as for LIDORT
!                          Complex Taylor routine dropped, not needed.

      DO N = 1, NLAYERS

        IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
          DO UM = 1, N_USER_STREAMS

!  setup

            UDEL = T_DELT_USERM(N,UM)
            SM = USER_SECANTS(UM)

!  Get the real multipliers
!  Version 2p7, 5/10/13, 2/19/14 - Call Taylor series routine

            DO K = 1, K_REAL(N)
              ZDEL_R    = T_DELT_EIGEN(K,N)
              ZUDEL_R   = ZDEL_R * UDEL
              THETA_2_R = ONE    - ZUDEL_R
              THETA_1_R = ZDEL_R - UDEL
              IF ( HSINGO(K,UM,N) )THEN
                EPS = ZETA_M(K,UM,N)
                CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAU_VERT(N), UDEL, SM, HMULT_1(K,UM,N) )
!                HMULT_1(K,UM,N) = SM * UDEL * DELTAU_VERT(N)   (Old code, was zero order)
              ELSE
                HMULT_1(K,UM,N) = SM * THETA_1_R * ZETA_M(K,UM,N)
              ENDIF
              HMULT_2(K,UM,N) = SM * THETA_2_R * ZETA_P(K,UM,N)
            ENDDO

!  Get the complex multipliers
!    Version 2p7. No Taylor-series degeneracy in the complex case.

            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              ZDEL_CR    = T_DELT_EIGEN(K1,N)
              ZDEL_CI    = T_DELT_EIGEN(K2,N)
              THETA_2_CR = ONE - ZDEL_CR * UDEL
              THETA_2_CI =     - ZDEL_CI * UDEL
              THETA_1_CR = ZDEL_CR - UDEL
              THETA_1_CI = ZDEL_CI
              HMULT_1(K1,UM,N) = SM * &
                   ( THETA_1_CR * ZETA_M(K1,UM,N) - &
                     THETA_1_CI * ZETA_M(K2,UM,N) )
              HMULT_1(K2,UM,N) = SM * &
                   ( THETA_1_CR * ZETA_M(K2,UM,N) + &
                     THETA_1_CI * ZETA_M(K1,UM,N) )
              HMULT_2(K1,UM,N) = SM * &
                   ( THETA_2_CR * ZETA_P(K1,UM,N) - &
                     THETA_2_CI * ZETA_P(K2,UM,N) )
              HMULT_2(K2,UM,N) = SM * &
                   ( THETA_2_CR * ZETA_P(K2,UM,N) + &
                     THETA_2_CI * ZETA_P(K1,UM,N) )
           ENDDO

!  Start loops over layers and user-streams

          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE WHOLELAYER_HMULT

!

      SUBROUTINE PARTLAYER_HMULT_UP ( &
          N, UT, N_USER_STREAMS, TAYLOR_ORDER,                    & ! Numbers, 2p7 Taylor-order argument added
          USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                 & ! streams, optical
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTUP_USERM, & ! Transmittances
          K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,              & ! Eigensolution stuff
          UT_HMULT_UU, UT_HMULT_UD )

!  1/31/21. Version 2.8.3.
!    -- UM loops always start with 1 now, not LOCAL_UM_START (argument dropped)

!  Partial layer INTEGRATED homogeneous solution MULTIPLIERS (Upwelling)

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS,  MAX_USER_STREAMS, MAX_USER_LEVELS, MAXEVALUES

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_1

      IMPLICIT NONE

!  INPUTS
!  ======

!  Numbers

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           TAYLOR_ORDER      ! 2p7 Taylor-order argument added

!  secants and optical

      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT  ( MAX_PARTLAYERS )

!  User stream transm.

      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenstream Trans.

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Eigensolution variables

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           HSINGO ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  OUTPUT
!  ======

!  Partial-layer multipliers
!mick fix 9/19/2017 - changed intent from "out" back to "inout"

      DOUBLE PRECISION, INTENT (INOUT) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Local variables
!  ---------------

      INTEGER          :: UM, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: UX_UP, SM, DX, EPS

!  real multipliers

      DOUBLE PRECISION :: ZDEL_R, ZX_UP_R, ZX_DN_R
      DOUBLE PRECISION :: THETA_DN_R, THETA_UP_R

!  complex multipliers

      DOUBLE PRECISION :: THETA_DN_CR, THETA_UP_CR
      DOUBLE PRECISION :: THETA_DN_CI, THETA_UP_CI
      DOUBLE PRECISION :: ZDEL_CR, ZX_UP_CR, ZX_DN_CR
      DOUBLE PRECISION :: ZDEL_CI, ZX_UP_CI, ZX_DN_CI

!  Partial layer multipliers

      DO UM = 1, N_USER_STREAMS

!  set-up

        UX_UP = T_UTUP_USERM(UT,UM)
        SM    = USER_SECANTS(UM)

!  real multipliers

        DO K = 1, K_REAL(N)
          ZDEL_R     = T_DELT_EIGEN(K,N)
          ZX_UP_R    = T_UTUP_EIGEN(K,UT)
          ZX_DN_R    = T_UTDN_EIGEN(K,UT)
          THETA_DN_R = ZX_DN_R - ZDEL_R * UX_UP
          THETA_UP_R = ZX_UP_R - UX_UP
          UT_HMULT_UD(K,UM,UT) = SM * THETA_DN_R * ZETA_P(K,UM,N)
          IF ( HSINGO(K,UM,N) )THEN
            DX  = DELTAU_VERT(N) - PARTAU_VERT(UT)
            EPS = ZETA_M(K,UM,N)
            CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DX, UX_UP, SM, UT_HMULT_UU(K,UM,UT) )
!            UT_HMULT_UU(K,UM,UT) = SM * UX_UP * DX  (Old code, was zero order)
          ELSE
            UT_HMULT_UU(K,UM,UT) = SM * THETA_UP_R * ZETA_M(K,UM,N)
          ENDIF
        ENDDO

!  Complex multipliers

        KO1 = K_REAL(N) + 1

        DO K = 1, K_COMPLEX(N)

          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1

          ZDEL_CR     = T_DELT_EIGEN(K1,N)
          ZDEL_CI     = T_DELT_EIGEN(K2,N)
          ZX_UP_CR    = T_UTUP_EIGEN(K1,UT)
          ZX_UP_CI    = T_UTUP_EIGEN(K2,UT)
          ZX_DN_CR    = T_UTDN_EIGEN(K1,UT)
          ZX_DN_CI    = T_UTDN_EIGEN(K2,UT)

          THETA_DN_CR = ZX_DN_CR - ZDEL_CR * UX_UP
          THETA_DN_CI = ZX_DN_CI - ZDEL_CI * UX_UP
          THETA_UP_CR = ZX_UP_CR - UX_UP
          THETA_UP_CI = ZX_UP_CI

          UT_HMULT_UD(K1,UM,UT) = SM * &
                   ( THETA_DN_CR * ZETA_P(K1,UM,N) - &
                     THETA_DN_CI * ZETA_P(K2,UM,N) )
          UT_HMULT_UD(K2,UM,UT) = SM * &
                   ( THETA_DN_CR * ZETA_P(K2,UM,N) + &
                     THETA_DN_CI * ZETA_P(K1,UM,N) )
          UT_HMULT_UU(K1,UM,UT) = SM * &
                   ( THETA_UP_CR * ZETA_M(K1,UM,N) - &
                     THETA_UP_CI * ZETA_M(K2,UM,N) )
          UT_HMULT_UU(K2,UM,UT) = SM * &
                   ( THETA_UP_CR * ZETA_M(K2,UM,N) + &
                     THETA_UP_CI * ZETA_M(K1,UM,N) )

        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE PARTLAYER_HMULT_UP

!

      SUBROUTINE PARTLAYER_HMULT_DN ( &
          N, UT, N_USER_STREAMS, TAYLOR_ORDER,                    & ! Numbers, 2p7 Taylor-order argument added
          USER_SECANTS, PARTAU_VERT,                              & ! streams, optical
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_UTDN_USERM, & ! Transmittances
          K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,              & ! Eigensolution stuff
          UT_HMULT_DD, UT_HMULT_DU )

!  Partial layer INTEGRATED homogeneous solution MULTIPLIERS (Downwelling)

!  1/31/21. Version 2.8.3.
!    -- UM loops always start with 1 now, not LOCAL_UM_START (argument dropped)

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_1

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS,  MAX_USER_STREAMS, MAX_USER_LEVELS, MAXEVALUES

      IMPLICIT NONE

!  INPUTS
!  ======

!  Numbers

      INTEGER, INTENT (IN) ::           N, UT
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           TAYLOR_ORDER      ! 2p7 Taylor-order argument added

!  secants and optical

      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT  ( MAX_PARTLAYERS )

!  User stream transm.

      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenstream Trans.

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Eigensolution variables

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           HSINGO ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  OUTPUT
!  ======

!  Partial-layer multipliers
!mick fix 9/19/2017 - changed intent from "out" back to "inout"

      DOUBLE PRECISION, INTENT (INOUT) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Local variables
!  ---------------

      INTEGER          :: UM, K, KO1, K0, K1, K2
      DOUBLE PRECISION :: UX_DN, SM, DX, EPS

!  real multipliers

      DOUBLE PRECISION :: ZDEL_R, ZX_UP_R, ZX_DN_R
      DOUBLE PRECISION :: THETA_DN_R, THETA_UP_R

!  complex multipliers

      DOUBLE PRECISION :: THETA_DN_CR, THETA_UP_CR
      DOUBLE PRECISION :: THETA_DN_CI, THETA_UP_CI
      DOUBLE PRECISION :: ZDEL_CR, ZX_UP_CR, ZX_DN_CR
      DOUBLE PRECISION :: ZDEL_CI, ZX_UP_CI, ZX_DN_CI

!  Partial layer multipliers

      DO UM = 1, N_USER_STREAMS

!  set-up

        UX_DN = T_UTDN_USERM(UT,UM)
        SM    = USER_SECANTS(UM)

!  real multipliers

        DO K = 1, K_REAL(N)
          ZDEL_R     = T_DELT_EIGEN(K,N)
          ZX_UP_R    = T_UTUP_EIGEN(K,UT)
          ZX_DN_R    = T_UTDN_EIGEN(K,UT)
          THETA_DN_R = ZX_DN_R - UX_DN
          THETA_UP_R = ZX_UP_R - ZDEL_R * UX_DN
          IF ( HSINGO(K,UM,N) ) THEN
            DX  =  PARTAU_VERT(UT)
            EPS = ZETA_M(K,UM,N)
            CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DX, UX_DN, SM, UT_HMULT_DD(K,UM,UT) )
          ELSE
            UT_HMULT_DD(K,UM,UT) = SM * THETA_DN_R * ZETA_M(K,UM,N)
          ENDIF
          UT_HMULT_DU(K,UM,UT) = SM * THETA_UP_R * ZETA_P(K,UM,N)
        ENDDO

!  Complex multipliers

        KO1 = K_REAL(N) + 1

        DO K = 1, K_COMPLEX(N)

          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1

          ZDEL_CR     = T_DELT_EIGEN(K1,N)
          ZDEL_CI     = T_DELT_EIGEN(K2,N)
          ZX_UP_CR    = T_UTUP_EIGEN(K1,UT)
          ZX_UP_CI    = T_UTUP_EIGEN(K2,UT)
          ZX_DN_CR    = T_UTDN_EIGEN(K1,UT)
          ZX_DN_CI    = T_UTDN_EIGEN(K2,UT)

          THETA_DN_CR = ZX_DN_CR - UX_DN
          THETA_DN_CI = ZX_DN_CI
          THETA_UP_CR = ZX_UP_CR - ZDEL_CR * UX_DN
          THETA_UP_CI = ZX_UP_CI - ZDEL_CI * UX_DN

          UT_HMULT_DD(K1,UM,UT) = SM * &
                   ( THETA_DN_CR * ZETA_M(K1,UM,N) - &
                     THETA_DN_CI * ZETA_M(K2,UM,N) )
          UT_HMULT_DD(K2,UM,UT) = SM * &
                   ( THETA_DN_CR * ZETA_M(K2,UM,N) + &
                     THETA_DN_CI * ZETA_M(K1,UM,N) )
          UT_HMULT_DU(K1,UM,UT) = SM * &
                   ( THETA_UP_CR * ZETA_P(K1,UM,N) - &
                     THETA_UP_CI * ZETA_P(K2,UM,N) )
          UT_HMULT_DU(K2,UM,UT) = SM * &
                   ( THETA_UP_CR * ZETA_P(K2,UM,N) + &
                     THETA_UP_CI * ZETA_P(K1,UM,N) )

        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE PARTLAYER_HMULT_DN

!

      SUBROUTINE VLIDORT_CBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM,                               & ! Input indices
        NSTOKES, NSTREAMS, NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS,     & ! Input numbers
        FLUX_FACTOR, DFLUX, QUAD_STREAMS, DO_LAYSCAT, OMEGA_GREEK, & ! Input Flux/Qauds/Optical
        BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT,  & ! Input beam attenuation
        PI_XQP, PI_XQM, PI_X0P, PI_XQM_POST,                       & ! Input PI matrices
        SAB, DAB, EIGENMAT_SAVE,                                   & ! Input matrices from RTE
        QSUMVEC_SAVE, QDIFVEC_SAVE, QVEC_SAVE, QDIF_SAVE,          & ! Output beam auxiliary vectors
        QMAT_SAVE, QPIVOT, BVEC, WUPPER, WLOWER,                   & ! Output matrix and Beam solutions
        STATUS, MESSAGE, TRACE )                                     ! Exception handling

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!  Linear Matrix algebra.

!  1/31/21. Version 2.8.3.
!    -- Subroutine renamed (QBEAM --> CBEAM). Argument list cleaned up.

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXBEAMS, MAXSTREAMS, MAXMOMENTS, MAXSTREAMS_2,  &
                                 MAXEVALUES, MAXSTRMSTKS, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, HALF, PI4

      USE LAPACK_TOOLS_m, only : DGETRF, DGETRS

      IMPLICIT NONE

!  INPUTS
!  ======

!  flags, removed 7/5/16
!      LOGICAL, INTENT (IN) ::           DO_DEBUG_WRITE
!      LOGICAL, INTENT (IN) ::           DO_FDTEST

!  indices

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS

!  Flux and quadrature
!      DOUBLE PRECISION, INTENT (IN) ::  DMAT  ( MAXSTOKES, MAXSTOKES )   ! removed from argument list 7/5/16

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )

!  Input bookkeeping + optical

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::           DO_LAYSCAT  ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Beam in average-secant parameterization

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  PI Matrices

      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Eigenproblem set-ups from RTE

      DOUBLE PRECISION, INTENT (IN) ::  SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  EIGENMAT_SAVE ( MAXEVALUES, MAXEVALUES, MAXLAYERS )

!  OUTPUT
!  ======

!  Auxiliary solution vectors

      DOUBLE PRECISION, INTENT (INOUT) ::  QSUMVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  QDIFVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  QVEC_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  QDIF_SAVE ( MAXSTRMSTKS )

!  Matrix and pivot

      DOUBLE PRECISION, INTENT (INOUT) ::  QMAT_SAVE ( MAXSTRMSTKS, MAXSTRMSTKS )
      INTEGER, INTENT (INOUT) ::           QPIVOT    ( MAXSTRMSTKS )

!  Particular integral solutions

      DOUBLE PRECISION, INTENT (INOUT) ::  BVEC   ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  help variables

      INTEGER ::          I, J, I1, L, N, M, O1, O2, IR, JR, IROW
      DOUBLE PRECISION :: SUM, GSUM, INV_X0SQ, SECBAR, F1
      DOUBLE PRECISION :: WSUM, WDIF, TRANS1, TRANS2, WVEC_D, WVEC_U
      DOUBLE PRECISION :: S_TPA, S_TMA, S_HELP, HELP
      DOUBLE PRECISION :: HELP_QFUNC ( MAXLAYERS, 0:MAXMOMENTS, MAXSTOKES )

!  Local error handling

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI, C3
      CHARACTER (LEN=2) :: CB

!  Initial section
!  ---------------

!  Initialize status

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  initialise indices

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  OR No scattering in this layer then:
!    ---> Zero the boundary layer values and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) .OR. .NOT. DO_LAYSCAT(M,N) ) THEN

! Enhancement # 23. S. Quesada 7/15/15, installed 2/3/16
        WUPPER(1:NSTREAMS_2,1:NSTOKES,N) = ZERO
        WLOWER(1:NSTREAMS_2,1:NSTOKES,N) = ZERO

!mick fix
        QSUMVEC_SAVE = ZERO
        QDIFVEC_SAVE = ZERO
        QVEC_SAVE = ZERO
        QDIF_SAVE = ZERO
        QMAT_SAVE = ZERO
        QPIVOT = 0
        BVEC(1:NSTREAMS_2,1:NSTOKES,N) = ZERO

        RETURN
      ENDIF

!  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR

!  Initialise matrix and column vector for solution
!  Changed at FMI 28 september 2004

! Enhancement # 24. S. Quesada 7/15/15, installed 2/3/16
      QMAT_SAVE(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS) = ZERO
      QVEC_SAVE(1:NSTKS_NSTRMS) = ZERO

!  Set up sum and difference vectors for Beam source terms
!  ( sum vector may be required again in linearization )
!  Auxiliary matrix for Q functions

! Enhancement # 25. S. Quesada 7/15/15, installed 2/3/16
      DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
          GSUM = ZERO
          DO O2 = 1, NSTOKES
            GSUM = GSUM + OMEGA_GREEK(L,N,O1,O2) * sum( PI_X0P(L,IBEAM,N,O2,1:NSTOKES) * DFLUX(1:NSTOKES) )
          ENDDO
          HELP_QFUNC(N,L,O1) = GSUM
        ENDDO
      ENDDO

! Enhancement # 26. S. Quesada 7/15/15, installed 2/3/16
      DO O1 = 1, NSTOKES
        DO I = 1, NSTREAMS
          S_TPA = sum( PI_XQP(M:NMOMENTS,I,O1,1:NSTOKES)      * HELP_QFUNC(N,M:NMOMENTS,1:NSTOKES) )
          S_TMA = sum( PI_XQM_POST(M:NMOMENTS,I,O1,1:NSTOKES) * HELP_QFUNC(N,M:NMOMENTS,1:NSTOKES) )
          QSUMVEC_SAVE(I,O1) = F1 * ( S_TPA + S_TMA ) / QUAD_STREAMS(I)
          QDIFVEC_SAVE(I,O1) = F1 * ( S_TPA - S_TMA ) / QUAD_STREAMS(I)
        ENDDO
      ENDDO

!  solution matrix for the reduced problem
!  ( matrix should be saved in the LU decomposition form)
!  Changed at FMI 28 september

! Enhancement # 27. S. Quesada 7/15/15, installed 2/3/16
      QMAT_SAVE(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS) = EIGENMAT_SAVE(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS,N)

      DO I = 1, NSTKS_NSTRMS
        QMAT_SAVE(I,I) = QMAT_SAVE(I,I) - INV_X0SQ
      ENDDO

!  RHS vector for the reduced problem
!  ( this vector will be the answer after the linear algebra solution,
!    and may be needed again if there is linearization )

! Enhancement # 28. S. Quesada 7/15/15, installed 2/3/16
      DO O1 = 1, NSTOKES
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          IROW = IR + O1
          QVEC_SAVE(IROW) = sum( DAB(I,1:NSTREAMS,O1,1:NSTOKES,N)*QSUMVEC_SAVE(1:NSTREAMS,1:NSTOKES) ) + QDIFVEC_SAVE(I,O1) * SECBAR
        ENDDO
      ENDDO

!  L-U decomposition of the solution matrix

      CALL DGETRF &
            ( NSTKS_NSTRMS, NSTKS_NSTRMS, &
              QMAT_SAVE, MAXSTRMSTKS, QPIVOT, INFO )

      IF ( INFO .GT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO
        WRITE(CB, '(I2)' ) IBEAM
        WRITE(C3, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRF call in VLIDORT_CBEAM_SOLUTION, Beam/layer'// &
            ' numbers = '//CB//', '//C3
        STATUS  = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  Solution of reduced problem by back-substitution

      CALL DGETRS &
          ( 'N', NSTKS_NSTRMS, 1, &
             QMAT_SAVE, MAXSTRMSTKS, QPIVOT, &
             QVEC_SAVE, MAXSTRMSTKS, INFO )

      IF ( INFO .LT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO
        WRITE(CB, '(I2)' ) IBEAM
        WRITE(C3, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRS call in VLIDORT_CBEAM_SOLUTION, Beam/layer'// &
            ' numbers = '//CB//', '//C3
        STATUS  = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  Assigning beam particular integral solution vector

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        IR = ( I - 1 ) * NSTOKES
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          HELP = ZERO
          DO J = 1, NSTREAMS
            JR = ( J - 1 ) * NSTOKES
! Enhancement # 29. S. Quesada 7/15/15, installed 2/3/16
            S_HELP = sum( SAB(I,J,O1,1:NSTOKES,N)*QVEC_SAVE(JR+1:JR+NSTOKES) )
            HELP = HELP + S_HELP
          ENDDO
          QDIF_SAVE(IROW) = ( HELP - QSUMVEC_SAVE(I,O1) ) / SECBAR
          WSUM = QVEC_SAVE(IROW)
          WDIF = QDIF_SAVE(IROW)
          WVEC_D = HALF * ( WSUM + WDIF )
          WVEC_U = HALF * ( WSUM - WDIF )
          BVEC(I,O1,N)  = WVEC_D
          BVEC(I1,O1,N) = WVEC_U
          IF ( O1 .GT. 2 ) THEN
            BVEC(I1,O1,N) = - WVEC_U
          ENDIF
!          if(i.eq.4.and.n.eq.24)write(38,'(4i3,1p2e20.10)')
!     &           m,ibeam,i,O1,BVEC(I,O1,N)
        ENDDO
      ENDDO

!  Old code: solution assignation

!      DO I = 1, NSTREAMS
!        I1 = I + NSTREAMS
!        DO O1 = 1, NSTOKES
!          BVEC(I,O1,N) = WVEC(I,O1,N)
!          HELP = ZERO
!          DO O2 = 1, NSTOKES
!            HELP = HELP + DMAT(O1,O2) * WVEC(I1,O2,N)
!          ENDDO
!          BVEC(I1,O1,N) = HELP
!        ENDDO
!      ENDDO

!  Values at the layer boundaries
!  (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1

! Enhancement # 30. S. Quesada 7/15/15, installed 2/3/16
      WUPPER(1:NSTREAMS_2, 1:NSTOKES, N) = BVEC(1:NSTREAMS_2, 1:NSTOKES, N) * TRANS1
      WLOWER(1:NSTREAMS_2, 1:NSTOKES, N) = BVEC(1:NSTREAMS_2, 1:NSTOKES, N) * TRANS2

!  debug

!      IF ( DO_DEBUG_WRITE ) THEN
!       J = 97
!       IF ( DO_FDTEST ) J = 98
!       IF ( N.EQ.14.and.m.eq.1) THEN
!        DO I = 1, NSTREAMS*2
!         write(*,'(4i4,1p4e17.9)')IBEAM,M,N,I,BVEC(I,1,N),
!     &      INITIAL_TRANS(N,IBEAM), T_DELT_MUBAR(N,IBEAM)
!        ENDDO
!       ENDIF
!       IF ( DO_FDTEST. and. M.EQ.1) PAUSE
!      ENDIF

!  debug check solution Forward problem - Works
!    removed 2 July. See Older version, labeled _debug

!      IF ( N.le.6.and.m.eq.0) THEN
!       DO I = 1,40
!        write(91,'(3i3,1p4e17.7)')M,N,I,(BVEC(I,o1,N),o1=1,4)
!       ENDDO
!      ENDIF

!      if (n.eq.6) pause 'end beam'

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CBEAM_SOLUTION

!

      SUBROUTINE VLIDORT_CUSER_SOLUTION ( &
        DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,                & ! Input flags
        GIVEN_LAYER, FOURIER, IBEAM, NSTOKES, NSTREAMS, N_USER_STREAMS,     & ! Input numbers
        NMOMENTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                   & ! Input bookkeeping
        DO_LAYSCAT, FLUX_FACTOR, DFLUX, BEAM_CUTOFF, QUAD_HALFWTS, & ! Input flux/quadrature/Bookkeeping
        OMEGA_GREEK, PI_XQP, PI_XUP, PI_XUM, PI_X0P, PI_XQM_PRE, BVEC,      & ! Input optical/PI/Beam-PI
        HELPSTOKES_BEAM, UPAR_DN_1, UPAR_DN_2, UPAR_UP_1, UPAR_UP_2 )         ! Output user solutions

!  1/31/21. Version 2.8.3.
!    -- Subroutine renamed (UBEAM --> CUSER). Argument list cleaned up.
!    -- UM loops always start with 1 now, not LOCAL_UM_START (argument dropped)

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, &
                                 MAXBEAMS, MAXMOMENTS, MAXSTREAMS_2, ZERO, PI4

      IMPLICIT NONE

!  INPUTS
!  ======

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY

!  Input numbers

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS

!  Input bookkeeping 

      INTEGER, INTENT (IN) ::           NMOMENTS
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )

      LOGICAL, INTENT (IN) ::           DO_LAYSCAT ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      INTEGER, INTENT (IN) ::           BEAM_CUTOFF  ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )

!  Input Optical

      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Input PI matrices

      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  input Particular integral vector

      DOUBLE PRECISION, INTENT (IN) ::  BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  OUTPUT - User solutions
!  =======================

      DOUBLE PRECISION, INTENT (INOUT) :: HELPSTOKES_BEAM ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Local variables
!  ---------------

      LOGICAL ::          DO_LAYER
      INTEGER ::          UM, L, N, M, O1, O2, J, J1
      DOUBLE PRECISION :: F1, SP, SM, SPS
      DOUBLE PRECISION :: HELP_Q1 ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: HELP_Q2 ( 0:MAXMOMENTS, MAXSTOKES )

      DOUBLE PRECISION, DIMENSION(MAXSTOKES)  ::  SGW =  (/ 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)

!  Variables removed by enhancement
!      DOUBLE PRECISION :: SUM1, SUM2, T1, T2, A5, SPT

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer.
!  OR No scattering in this layer
!    --> Zero the boundary layer values and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) .OR. .NOT. DO_LAYSCAT(M,N) ) THEN
        IF ( DO_UPWELLING ) THEN
! Enhancement # 31. S. Quesada 7/15/15, installed 2/3/16
          UPAR_UP_1(1:N_USER_STREAMS, 1:NSTOKES, N) = ZERO
          UPAR_UP_2(1:N_USER_STREAMS, 1:NSTOKES, N) = ZERO
        ENDIF
        IF ( DO_DNWELLING ) THEN
! Enhancement # 32. S. Quesada 7/15/15, installed 2/3/16
          UPAR_DN_1(1:N_USER_STREAMS, 1:NSTOKES, N) = ZERO
          UPAR_DN_2(1:N_USER_STREAMS, 1:NSTOKES, N) = ZERO
        ENDIF
        RETURN
      ENDIF

!  existence flag

      DO_LAYER = ( STERM_LAYERMASK_UP(N) .OR. &
                   STERM_LAYERMASK_DN(N) )

!  first function

      IF ( DO_LAYER ) THEN
! Enhancement # 33. S. Quesada 7/15/15, installed 2/3/16
       DO O1 = 1, NSTOKES
        DO L = M, NMOMENTS
          SPS = ZERO
          DO O2 = 1, NSTOKES
            SPS = SPS + OMEGA_GREEK(L,N,O1,O2) * sum( PI_X0P(L,IBEAM,N,O2,1:NSTOKES) * DFLUX(1:NSTOKES) )
          ENDDO
          HELP_Q1(L,O1) = F1 * SPS
        ENDDO
      ENDDO
     ENDIF

!  For each moment, do inner sum over computational angles
!    Equivalent to Gvec = GVECTOR. Eq 34. GAUx = Greekmat.Gvec

      IF ( DO_LAYER ) THEN
! Enhancement # 34. S. Quesada 7/15/15, installed 2/3/16
        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           J1 = J + NSTREAMS
           SP = sum( PI_XQP    (L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J)*BVEC(J,1:NSTOKES,N) )
           SM = sum( SGW(1:NSTOKES)*PI_XQM_PRE(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J)*BVEC(J1,1:NSTOKES,N) )
           SPS = SPS + SP + SM
          ENDDO
          HELPSTOKES_BEAM(L,O1) = SPS
         ENDDO
         DO O1 = 1, NSTOKES
          HELP_Q2(L,O1) = sum( OMEGA_GREEK(L,N,O1,1:NSTOKES) * HELPSTOKES_BEAM(L,1:NSTOKES) )
         ENDDO
        ENDDO
      ENDIF

!  Now sum over all harmonic contributions (Upwelling)
!    Direct and integarated contributions

      IF ( DO_UPWELLING ) THEN
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
! Enhancement # 35. S. Quesada 7/15/15, installed 2/3/16
         DO O1 = 1, NSTOKES
          DO UM = 1, N_USER_STREAMS
           UPAR_UP_1(UM,O1,N) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,UM,O1,1:NSTOKES) )
           UPAR_UP_2(UM,O1,N) = sum( HELP_Q2(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,UM,O1,1:NSTOKES) )
          ENDDO
         ENDDO
        ELSE
! Enhancement # 36. S. Quesada 7/15/15, installed 2/3/16
         DO O1 = 1, NSTOKES
          UPAR_UP_1(IBEAM,O1,N) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,IBEAM,O1,1:NSTOKES) )
          UPAR_UP_2(IBEAM,O1,N) = sum( HELP_Q2(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,IBEAM,O1,1:NSTOKES) )
         ENDDO
        ENDIF
      ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct and integarated contributions

      IF ( DO_DNWELLING ) THEN
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
! Enhancement # 37. S. Quesada 7/15/15, installed 2/3/16
         DO O1 = 1, NSTOKES
          DO UM = 1, N_USER_STREAMS
           UPAR_DN_1(UM,O1,N) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES) )
           UPAR_DN_2(UM,O1,N) = sum( HELP_Q2(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES) )
          ENDDO
         ENDDO
        ELSE
! Enhancement # 38. S. Quesada 7/15/15, installed 2/3/16
         DO O1 = 1, NSTOKES
          UPAR_DN_1(IBEAM,O1,N) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,IBEAM,O1,1:NSTOKES) )
          UPAR_DN_2(IBEAM,O1,N) = sum( HELP_Q2(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,IBEAM,O1,1:NSTOKES) )
         ENDDO
        ENDIF
      ENDIF

!  debug stokes = 1
!      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
!        DO UM = 1, N_USER_STREAMS
!          write(*,*)n,um,UPAR_UP_1(UM,1,N),UPAR_UP_2(UM,1,N)
!        ENDDO
!        DO UM = 1, N_USER_STREAMS
!          write(*,*)n,um,UPAR_DN_1(UM,1,N),UPAR_DN_2(UM,1,N)
!        ENDDO
!      ELSE
!        write(*,*)n,ibeam,UPAR_UP_1(IBEAM,1,N),UPAR_UP_2(IBEAM,1,N)
!        write(*,*)n,ibeam,UPAR_DN_1(IBEAM,1,N),UPAR_DN_2(IBEAM,1,N)
!      ENDIF
!      PAUSE

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CUSER_SOLUTION

!

      SUBROUTINE VLIDORT_GBEAM_SOLUTION &
        ( FOURIER, IBEAM, GIVEN_LAYER, TAYLOR_ORDER,                      & ! Input Control
          NSTOKES, NSTREAMS, NSTREAMS_2, NSTKS_NSTRMS, NMOMENTS,          & ! Input computational
          FLUX_FACTOR, DFLUX, QUAD_WTS, DO_LAYSCAT,                       & ! Input bookkeeping
          BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,       & ! Input average-secant
          DELTAUS, OMEGA_GREEK, PI_XQP, PI_XQM, PI_X0P, PI_XQM_POST,      & ! Input Optical/PI
          KEIGEN, KEIGEN_CSQ, K_REAL, K_COMPLEX,                          & ! Input Eigenvalues                      
          SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, NORM_SAVED,                 & ! Input Eigensolutions
          DPI, DMI, ATERM_SAVE, BTERM_SAVE, CFUNC, DFUNC,                 & ! Output
          GAMMA_M, GAMMA_P, GFUNC_UP, GFUNC_DN, WUPPER, WLOWER )            ! Output

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXBEAMS, MAXSTREAMS, MAXMOMENTS, MAXSTREAMS_2,  &
                                 MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4, TWO, TAYLOR_SMALL, SDU

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_1

!  1/31/21.  Version 2.8.3. COMPLETELY NEW SUBROUTINE
!    -- Green's function solar-beam particular integral, one layer only.
!    -- Added output DPI, DMI (no longer local).
!    -- includes layer dimension on KEIGEN_CSQ

      IMPLICIT NONE

!  subroutine input arguments
!  ==========================

!  Given layer index, Beam index and Fourier number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Numbers

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2, NSTKS_NSTRMS
      INTEGER  , intent(in)  :: NMOMENTS, NSTOKES

!  Flux and quadrature

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WTS ( MAXSTREAMS )

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYSCAT (0:MAXMOMENTS, MAXLAYERS)

!  Beam in average-secant parameterization
!    BEAM_CUTOFF = Last layer to include Particular integral solution

      INTEGER, INTENT (IN) ::           BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Input Optical depths required for Taylor-series limiting cases

      DOUBLE PRECISION, intent(in)  :: DELTAUS (MAXLAYERS)

!  PI Matrices, OMEGA_GREEK

      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Eigensolution bookkeeping

      INTEGER, INTENT (in) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (in) ::           K_COMPLEX ( MAXLAYERS )

!  Eigenvalues

      DOUBLE PRECISION, intent(in)  :: KEIGEN     (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: KEIGEN_CSQ (MAXEVALUES,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      DOUBLE PRECISION, intent(in)  :: T_DELT_EIGEN (MAXEVALUES,MAXLAYERS)

!  Eigensolutions and normalization factors

      DOUBLE PRECISION, intent(in)  :: SOLA_XPOS (MAXSTREAMS_2,MAXSTOKES,MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: SOLB_XNEG (MAXSTREAMS_2,MAXSTOKES,MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, intent(in)  :: NORM_SAVED (MAXEVALUES,MAXLAYERS)

!  subroutine output arguments
!  ===========================

!  Green function solution
!  ***********************

!mick fix 6/29/11 - change most outputs from "out" to "inout"

!  Saved quantities for the Green function solution

      DOUBLE PRECISION, intent(inout) :: ATERM_SAVE (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(inout) :: BTERM_SAVE (MAXEVALUES,MAXLAYERS)

!  Layer C and D functions

      DOUBLE PRECISION, intent(inout) :: CFUNC (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(inout) :: DFUNC (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(inout) :: GFUNC_UP (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(inout) :: GFUNC_DN (MAXEVALUES,MAXLAYERS)

!  Holding arrays for multiplier coefficients

      DOUBLE PRECISION, intent(inout) :: GAMMA_M (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(inout) :: GAMMA_P (MAXEVALUES,MAXLAYERS)

!  General beam solutions at the boundaries

      DOUBLE PRECISION, intent(inout) :: WUPPER (MAXSTREAMS_2,MAXSTOKES,MAXLAYERS)
      DOUBLE PRECISION, intent(inout) :: WLOWER (MAXSTREAMS_2,MAXSTOKES,MAXLAYERS)

!  Must be output now, for use in the linearizations

      DOUBLE PRECISION, intent(inout) :: DMI (MAXSTREAMS,MAXSTOKES), DPI (MAXSTREAMS,MAXSTOKES)

!  local variables
!  ===============

!  Arrays (status unclear, may be output)
!      DOUBLE PRECISION  :: AGM (MAXEVALUES,MAXLAYERS), BGP (MAXEVALUES,MAXLAYERS)

      INTEGER             :: L, I, I1, M, N, NS1, NS2, KO1, K, K1, K2, K0, O1, O2
      DOUBLE PRECISION    :: EPS, CONST, SECBAR, A_USER, B_USER, WDEL, ZDEL_R, ZWDEL_R, F1
      DOUBLE PRECISION    :: RHO_P_CR, RHO_M_CR, A_HELP, B_HELP, MODULUS_P, MODULUS_M, GSUM, CSUM, DSUM
      DOUBLE PRECISION    :: ZDEL_CR, ZDEL_CI, THETA_2_CR, THETA_2_CI, THETA_1_CR, THETA_1_CI

      DOUBLE PRECISION    :: HELP1 (NSTREAMS), HELP2 (NSTREAMS)
      DOUBLE PRECISION    :: HELP_QFUNC(0:MAXMOMENTS,MAXSTOKES)

!  initialise indices

      N = GIVEN_LAYER
      M = FOURIER
!mick eff 3/22/2017
      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  No particular solution beyond the cutoff layer, Or no scattering in this layer...
!  ... Zero the boundary layer values and exit. Mick eff 3/22/2017

!  initialise indices

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  OR No scattering in this layer then:
!    ---> Zero the boundary layer values and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) .OR. .NOT.DO_LAYSCAT(M,N) ) THEN
        WUPPER(1:NSTREAMS_2,1:NSTOKES,N) = ZERO
        WLOWER(1:NSTREAMS_2,1:NSTOKES,N) = ZERO
        RETURN
      ENDIF

!  constants for the layer

      SECBAR = AVERAGE_SECANT(N,IBEAM)
      CONST  = INITIAL_TRANS (N,IBEAM)
      WDEL   = T_DELT_MUBAR  (N,IBEAM)
      
!  3. Optical depth integrations for the discrete ordinate solution
!     =============================================================

!  Real-values 

      DO K = 1, K_REAL(N)
         ZDEL_R  = T_DELT_EIGEN(K,N)
         ZWDEL_R = ZDEL_R * WDEL
         GAMMA_P(K,N) = SECBAR + KEIGEN(K,N)
         GAMMA_M(K,N) = SECBAR - KEIGEN(K,N) ; EPS = GAMMA_M(K,N)
         IF ( ABS(EPS) .LT. TAYLOR_SMALL ) THEN
            CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAUS(N), WDEL, CONST, CFUNC(K,N)   )
         ELSE
            CFUNC(K,N) = CONST * ( ZDEL_R - WDEL ) / GAMMA_M(K,N)
         ENDIF
         DFUNC(K,N)  = CONST * ( ONE - ZWDEL_R ) / GAMMA_P(K,N)
      ENDDO

!  Complex values 

      IF ( K_COMPLEX(N) .GT. 0 ) THEN
         KO1 = K_REAL(N) + 1
         A_USER = SECBAR * SECBAR
         B_USER = TWO * SECBAR
         DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
            RHO_P_CR  = SECBAR + KEIGEN(K1,N)    ; RHO_M_CR  = SECBAR - KEIGEN(K1,N)
            A_HELP    = A_USER + KEIGEN_CSQ(K,N) ; B_HELP    = B_USER * KEIGEN(K1,N)
            MODULUS_P = A_HELP + B_HELP          ; MODULUS_M = A_HELP - B_HELP
            GAMMA_P(K1,N) = RHO_P_CR / MODULUS_P ; GAMMA_P(K2,N) = -KEIGEN(K2,N) / MODULUS_P
            GAMMA_M(K1,N) = RHO_M_CR / MODULUS_M ; GAMMA_M(K2,N) =  KEIGEN(K2,N) / MODULUS_M
            ZDEL_CR    = T_DELT_EIGEN(K1,N)   ; ZDEL_CI    = T_DELT_EIGEN(K2,N)
            THETA_2_CR = ONE - ZDEL_CR * WDEL ; THETA_2_CI =   - ZDEL_CI * WDEL
            THETA_1_CR = ZDEL_CR - WDEL       ; THETA_1_CI = ZDEL_CI
            CFUNC(K1,N) = CONST * ( THETA_1_CR * GAMMA_M(K1,N) - THETA_1_CI * GAMMA_M(K2,N) )
            CFUNC(K2,N) = CONST * ( THETA_1_CR * GAMMA_M(K2,N) + THETA_1_CI * GAMMA_M(K1,N) )
            DFUNC(K1,N) = CONST * ( THETA_2_CR * GAMMA_P(K1,N) - THETA_2_CI * GAMMA_P(K2,N) )
            DFUNC(K2,N) = CONST * ( THETA_2_CR * GAMMA_P(K2,N) + THETA_2_CI * GAMMA_P(K1,N) )
         ENDDO
      ENDIF

!  4. Form quantities independent of optical depth
!     ============================================

!  Set up help arrays (independent of eigenvector). Mick eff 3/22/2017
!mick fix 1/5/2021 - moved defining of DPI & DMI from inside first set of loops below to after it
!                    (so all needed elements of the help array HELP_QFUNC are defined first)

      DO O1 = 1, NSTOKES
        DO L = M, NMOMENTS
          GSUM = ZERO
          DO O2 = 1, NSTOKES
            GSUM = GSUM + OMEGA_GREEK(L,N,O1,O2) * sum( PI_X0P(L,IBEAM,N,O2,1:NSTOKES) * DFLUX(1:NSTOKES) )
          ENDDO
          HELP_QFUNC(L,O1) = GSUM
        ENDDO
      ENDDO
      DO O1 = 1, NSTOKES
        DO I = 1, NSTREAMS
          DPI(I,O1) = F1 * sum( PI_XQP(M:NMOMENTS,I,O1,1:NSTOKES)      * HELP_QFUNC(M:NMOMENTS,1:NSTOKES) )
          DMI(I,O1) = F1 * sum( PI_XQM_POST(M:NMOMENTS,I,O1,1:NSTOKES) * HELP_QFUNC(M:NMOMENTS,1:NSTOKES) )
        ENDDO
      ENDDO

!  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE. Mick eff 3/22/2017
!   --Real Valued solution

      DO K = 1, K_REAL(N)
        DO I = 1, NSTREAMS
           HELP1(I) = dot_product(DPI(I,1:nstokes),SOLA_XPOS(I,1:nstokes,K,N)) &
                    + dot_product(DMI(I,1:nstokes),SOLB_XNEG(I,1:nstokes,K,N))
           HELP2(I) = dot_product(DMI(I,1:nstokes),SOLA_XPOS(I,1:nstokes,K,N)) &
                    + dot_product(DPI(I,1:nstokes),SOLB_XNEG(I,1:nstokes,K,N))
        ENDDO
        ATERM_SAVE(K,N) = DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS),HELP1 ) / NORM_SAVED(K,N)
        BTERM_SAVE(K,N) = DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS),HELP2 ) / NORM_SAVED(K,N)
      ENDDO

!  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE. Mick eff 3/22/2017
!   -- Complex Valued solution is a PLACEHOLDER In Version 2.8.3

!  5. Green function multipliers
!     ==========================

      DO I = 1, NSTKS_NSTRMS
         CSUM = CFUNC(I,N) ; GFUNC_DN(I,N) = CSUM * ATERM_SAVE(I,N)
         DSUM = DFUNC(I,N) ; GFUNC_UP(I,N) = DSUM * BTERM_SAVE(I,N)
      ENDDO

!  6. Set particular integral from Green function expansion
!     =====================================================

!  Particular integrals at lower and upper boundaries.
!   Only works for NSTOKES = 1 or 3.

      DO O1 = 1, NSTOKES
         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            WUPPER(I, O1,N) = DOT_PRODUCT ( GFUNC_UP(1:NSTKS_NSTRMS,N),SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N) )
            WUPPER(I1,O1,N) = DOT_PRODUCT ( GFUNC_UP(1:NSTKS_NSTRMS,N),SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N) )
            WLOWER(I1,O1,N) = DOT_PRODUCT ( GFUNC_DN(1:NSTKS_NSTRMS,N),SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N) )
            WLOWER(I, O1,N) = DOT_PRODUCT ( GFUNC_DN(1:NSTKS_NSTRMS,N),SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N) )
            IF ( O1.gt.2) then
              WLOWER(I1,O1,N) = - WLOWER(I1,O1,N)
              WUPPER(I1,O1,N) = - WUPPER(I1,O1,N)
            ENDIF
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_GBEAM_SOLUTION

!

      SUBROUTINE VLIDORT_GUSER_SOLUTION ( &
            DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,      & ! Input flags
            LAYER, FOURIER, IBEAM, NSTOKES, N_USER_STREAMS, NMOMENTS, & ! Input numbers
            DO_LAYSCAT, BEAM_CUTOFF, FLUX_FACTOR,                     & ! Input Bookkeeping
            DFLUX, OMEGA_GREEK, PI_XUP, PI_XUM, PI_X0P,               & ! Input optical/PI/flux
            UPAR_DN_1, UPAR_UP_1 )                                      ! Output user solutions

!  1/31/21.  Version 2.8.3. COMPLETELY NEW SUBROUTINE
!    -- Green's function solar-beam particular integral, one layer only, User solutions

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, &
                                 MAXBEAMS, MAXMOMENTS, MAXSTREAMS_2, ZERO, PI4

      IMPLICIT NONE

!  INPUTS
!  ======

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY

!  Input numbers

      INTEGER, INTENT (IN) ::           LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS

!  Input bookkeeping 

      LOGICAL, INTENT (IN) ::           DO_LAYSCAT ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      INTEGER, INTENT (IN) ::           BEAM_CUTOFF  ( MAXBEAMS )

!  Input Optical

      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Input PI matrices

      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  OUTPUT - User solutions
!  =======================

      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Local variables
!  ---------------

      INTEGER ::          UM, L, N, M, O1, O2
      DOUBLE PRECISION :: F1, SPS
      DOUBLE PRECISION :: HELP_Q1 ( 0:MAXMOMENTS, MAXSTOKES )

!  Scattering solutions
!  ====================

!  Layer and Fourier

      N = LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer.
!  OR No scattering in this layer
!    --> Zero the boundary layer values and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) .OR. .NOT. DO_LAYSCAT(M,N) ) THEN
        IF ( DO_UPWELLING ) UPAR_UP_1(1:N_USER_STREAMS, 1:NSTOKES, N) = ZERO
        IF ( DO_DNWELLING ) UPAR_DN_1(1:N_USER_STREAMS, 1:NSTOKES, N) = ZERO
        RETURN
      ENDIF

!  For each moment do inner sum over computational angles

      DO O1 = 1, NSTOKES
        DO L = M, NMOMENTS
          SPS = ZERO
          DO O2 = 1, NSTOKES
            SPS = SPS + OMEGA_GREEK(L,N,O1,O2) * sum( PI_X0P(L,IBEAM,N,O2,1:NSTOKES) * DFLUX(1:NSTOKES) )
          ENDDO
          HELP_Q1(L,O1) = F1 * SPS
        ENDDO
      ENDDO

!  Now sum over all harmonic contributions (Upwelling)
!    Direct contribution Only. Observational Geometry, tie to the IBEAM index.

      IF ( DO_UPWELLING ) THEN
         IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO O1 = 1, NSTOKES
               DO UM = 1, N_USER_STREAMS
                  UPAR_UP_1(UM,O1,N) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,UM,O1,1:NSTOKES) )
               ENDDO
            ENDDO
         ELSE
            DO O1 = 1, NSTOKES
               UPAR_UP_1(IBEAM,O1,N) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,IBEAM,O1,1:NSTOKES) )
            ENDDO
         ENDIF
      ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct contribution Only. Observational Geometry, tie to the IBEAM index.

      IF ( DO_DNWELLING ) THEN
         IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO O1 = 1, NSTOKES
               DO UM = 1, N_USER_STREAMS
                  UPAR_DN_1(UM,O1,N) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES) )
               ENDDO
            ENDDO
         ELSE
! Enhancement # 38. S. Quesada 7/15/15, installed 2/3/16
            DO O1 = 1, NSTOKES
               UPAR_DN_1(IBEAM,O1,N) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,IBEAM,O1,1:NSTOKES) )
            ENDDO
         ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_GUSER_SOLUTION

!  End module

      END MODULE vlidort_solutions_m

