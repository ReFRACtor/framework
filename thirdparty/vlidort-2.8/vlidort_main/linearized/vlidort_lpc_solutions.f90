
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
! #              VLIDORT_L_QHOM_SOLUTION                        #
! #              VLIDORT_L_QHOM_NORMS                           #
! #              VLIDORT_L_UHOM_SOLUTION                        #
! #                                                             #
! #              L_HMULT_MASTER                                 #
! #                                                             #
! #              VLIDORT_L_GBEAM_SOLUTION                       #
! #              VLIDORT_L_GUSER_SOLUTION                       #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. Major Changes for Green's function option
!    ==> New Module VLIDORT_L_QHOM_NORMS (needed for the Green's function solution)
!    ==> New Modules VLIDORT_L_GBEAM_SOLUTION, VLIDORT_L_GUSER_SOLUTION

!  Notes for Version 2.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix 05/10/13  - HSINGO defined using TAYLOR_SMALL
!     Rob Fix 02/19/14  - ZETA_M = RHO_M for the degenerate case....
!                         otherwise ZETA_M = ONE / RHO_M
!                         COMPLEX HSINGO removed.
!     Rob Fix 02/09/16  - Do-loop optimatization by S, Quesada implemeted

      MODULE vlidort_lpc_solutions_m

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_L_1

      PRIVATE
      PUBLIC :: VLIDORT_L_QHOM_SOLUTION , &
                VLIDORT_L_QHOM_NORMS    , &
                VLIDORT_L_UHOM_SOLUTION , &
                L_HMULT_MASTER          , &
                VLIDORT_L_GBEAM_SOLUTION, &
                VLIDORT_L_GUSER_SOLUTION

      CONTAINS

      SUBROUTINE VLIDORT_L_QHOM_SOLUTION ( &
        DO_SOLUTION_SAVING, GIVEN_LAYER, FOURIER,                         & ! Input FlagL_QG and Indices
        NSTOKES, NSTREAMS, N_PARTLAYERS, DO_VARY, N_PARAMETERS, NMOMENTS, & ! Input Numbers+Lin.Control
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, PARTLAYERS_LAYERIDX,    & ! Input Numbers
        DO_REAL_EIGENSOLVER, DO_LAYER_SCAT, QUAD_STREAMS, QUAD_HALFWTS,   & ! Input bookkeeping and quadrature
        DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT, L_OMEGA_GREEK,           & ! Input Optical
        PI_XQP, PI_XQM_PRE, SAB, DAB, EIGENMAT_SAVE, REAL_KSQ, IMAG_KSQ,  & ! Input Homog solution
        LEFT_EVEC, RITE_EVEC, EIGENDEGEN, EIGENMASK_R, EIGENMASK_C,       & ! Input Homog solution
        K_REAL, K_COMPLEX, KEIGEN, KEIGEN_CSQ, FWD_SUMVEC, FWD_DIFVEC,    & ! Input Homog solution
        T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                         & ! Input Homog solution
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,             & ! Input linearized Homog Trans
        L_SAB, L_DAB, L_EIGENMAT, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG,     & ! Output lineaerized Homog. Sol.
        L_T_DELT_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,                   & ! Output linearized Homog Trans
        STATUS, MESSAGE, TRACE )                                            ! Exception handling

!  Linearization of the homogeneous solutions.

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, MAXMOMENTS, &
                                 MAXSTREAMS_2, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTRMSTKS_21, MAXSTRMSTKS_P1,   &
                                 MAXEVALUES, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, HALF, ONE, TWO

      USE LAPACK_TOOLS_m, Only : DGETRS, DGETRF


      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING

!  indices

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER

!  Linearization control

      LOGICAL, INTENT (IN) ::          DO_VARY
      INTEGER, INTENT (IN) ::          N_PARAMETERS

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS

      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )

!  Bookkeeping

      LOGICAL, INTENT (IN) ::          DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCAT       ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Pi matrices and eignevalue problem arrays

      DOUBLE PRECISION, INTENT (IN)    :: PI_XQP     ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN)    :: PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: EIGENMAT_SAVE ( MAXEVALUES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: REAL_KSQ  ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: IMAG_KSQ  ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: LEFT_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: RITE_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )

      LOGICAL, INTENT (INOUT) ::          EIGENDEGEN  ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          EIGENMASK_R ( MAXEVALUES )
      INTEGER, INTENT (INOUT) ::          EIGENMASK_C ( MAXEVALUES )

!  Homogeneous solution results
!    -- 1/31/21. Version 2.8.3. KEIGEN_CSQ now has an extra "LAYERS" dimension

      INTEGER, INTENT (INOUT) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: KEIGEN     ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: KEIGEN_CSQ ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: FWD_SUMVEC ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) :: FWD_DIFVEC ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  homog. solution transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Linearized discrete ordinate tranmsittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  OUTPUT
!  ======

!  Linearized  eigenmatrices

      DOUBLE PRECISION, INTENT (INOUT) ::  L_SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_EIGENMAT ( MAXEVALUES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous solution

      DOUBLE PRECISION, INTENT (INOUT) ::  L_KEIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous transmittances

      DOUBLE PRECISION, INTENT (INOUT) ::  L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Exception (Strictly out)

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ===============

!  local matrices for eigenvalue computation
!  -----------------------------------------

!  For the real calculation based on DGEEV

      DOUBLE PRECISION :: HMAT1(MAXSTRMSTKS,MAXSTRMSTKS)
      DOUBLE PRECISION :: HVEC1(MAXSTRMSTKS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: HVEC_AUX(MAXSTRMSTKS,MAX_ATMOSWFS)

!  For the complex calculations based on DGEEV

      DOUBLE PRECISION :: HMAT2(MAXSTRMSTKS_2,MAXSTRMSTKS_2)
      DOUBLE PRECISION :: HVEC2(MAXSTRMSTKS_2,MAX_ATMOSWFS)

      DOUBLE PRECISION :: HMATRED(MAXSTRMSTKS_21,MAXSTRMSTKS_21)
      DOUBLE PRECISION :: HVECRED(MAXSTRMSTKS_21,MAX_ATMOSWFS)

      DOUBLE PRECISION :: HVEC_CR(MAXSTRMSTKS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: HVEC_CI(MAXSTRMSTKS,MAX_ATMOSWFS)

      DOUBLE PRECISION :: ADJN_CR(MAXSTRMSTKS)
      DOUBLE PRECISION :: ADJN_CI(MAXSTRMSTKS)

!  For the real calculation based on ASYMTX

      DOUBLE PRECISION ::  HMAT1A(MAXSTRMSTKS_P1,MAXSTRMSTKS_P1)
      DOUBLE PRECISION ::  HVEC1A(MAXSTRMSTKS_P1,MAX_ATMOSWFS)

!  Local linearized vectors

      DOUBLE PRECISION ::  L_FWD_DIFVEC_R(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION ::  L_FWD_SUMVEC_R(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)

      DOUBLE PRECISION ::  L_FWD_DIFVEC_CR(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION ::  L_FWD_SUMVEC_CR(MAXSTREAMS,MAXSTOKES)

      DOUBLE PRECISION ::  L_FWD_DIFVEC_CI(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION ::  L_FWD_SUMVEC_CI(MAXSTREAMS,MAXSTOKES)

!  Pivoting arrays for LAPACK linear-algebra solvers

      INTEGER ::          LAPACK_IPIV1  (MAXSTRMSTKS)
      INTEGER ::          LAPACK_IPIV21 (MAXSTRMSTKS_21)
      INTEGER ::          LAPACK_IPIV1A (MAXSTRMSTKS_P1)

!  Miscellaneous local variables

      INTEGER ::          I, J, I1, JC, IR, JCOL, IROW, IROW1, ILAST
      INTEGER ::          L, N, M, Q, O1, O2, O3, NA, UT
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: DP, DM, KSQD, KVAL, KVAL2, XINV, FAC
      DOUBLE PRECISION :: TBAS, TBAS_CR, TBAS_CI, LCR, LCI, H1P, H1M
      DOUBLE PRECISION :: TAU_UP, TAU_DN
      DOUBLE PRECISION :: TBASUP_CR, TBASUP_CI, TBASDN_CR, TBASDN_CI
      INTEGER ::          LAPACK_INFO, NSTKS_NSTRMS_21, NSTMSTKS_P1

!  additional variables

      INTEGER           :: AA, AA1, TC, TR, I_NONTRIVIAL
      DOUBLE PRECISION  :: SR, SCR, SCI, HR, HCR, HCI, MODULUS
      DOUBLE PRECISION  :: KVAL_CR, KVAL_CI, KVAL2_CR, KVAL2_CI
      DOUBLE PRECISION  :: MOD_KVAL, SEPCON_CR, SEPCON_CI
      DOUBLE PRECISION  :: DENOM_CR, DENOM_CI
      CHARACTER (LEN=3) :: CN, CI

!  These variables declared for the do-loop optimization. 2/9/16, 5/9/16 Using Symbolic dimensions.

!      DOUBLE PRECISION  :: LFPD1, LFPD2, LFND1, LFND2, LFPD, LFND, removed 2/9/16
!      DOUBLE PRECISION  :: LFPD(NSTREAMS,NSTOKES), LFND(NSTREAMS,NSTOKES)
!      DOUBLE PRECISION  :: LFPD1(NSTREAMS,NSTOKES), LFPD2(NSTREAMS,NSTOKES),
!      DOUBLE PRECISION  :: LFND1(NSTREAMS,NSTOKES), LFND2(NSTREAMS,NSTOKES)

      DOUBLE PRECISION :: VAR(MAX_ATMOSWFS)
      DOUBLE PRECISION  :: LFPD (MAXSTREAMS,MAXSTOKES), LFND (MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION  :: LFPD1(MAXSTREAMS,MAXSTOKES), LFPD2(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION  :: LFND1(MAXSTREAMS,MAXSTOKES), LFND2(MAXSTREAMS,MAXSTOKES)

!  Logical flags, Introduced in Version 2.8 to avoid GOTO statements

      LOGICAL :: LOOP
      LOGICAL :: LINEARIZE_EMATRIX
      LOGICAL :: DO_LINEARIZE_ASYMTX
      LOGICAL :: DO_LINEARIZE_DGEEV

!  Start of code
!  -------------

!  initialise status

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Layer and Fourier number

      N = GIVEN_LAYER
      M = FOURIER

!  initialize flags. New for Version 2.8

      LINEARIZE_EMATRIX   = .true.
      DO_LINEARIZE_ASYMTX = .false.
      DO_LINEARIZE_DGEEV  = .false.

!  nothing to do if not varying

      IF ( .NOT. DO_VARY ) RETURN

!  For solution saving option with no scattering,
!  linearized solution vectors are all zero
!  Exit the routine by going to continuation point.
!    Version 2.8 remove GOTO 3456 statement, replace by IF control.
!  Enhancement # 1, 2/9/16

      IF ( DO_SOLUTION_SAVING ) THEN
        IF ( .NOT. DO_LAYER_SCAT(M,N) ) THEN
          L_SOLA_XPOS(1:NSTREAMS_2, 1:NSTOKES, 1:K_REAL(N), N, 1:N_PARAMETERS)  = ZERO
          L_SOLB_XNEG(1:NSTREAMS_2, 1:NSTOKES, 1:K_REAL(N), N, 1:N_PARAMETERS)  = ZERO
          L_KEIGEN(1:K_REAL(N),N,1:N_PARAMETERS)  = ZERO
          LINEARIZE_EMATRIX = .false.   !!! GOTO 3456
        ENDIF
      ENDIF

!  Linearize the Eigenmatrix
!  =========================

!  start clause. IF statement introduced for Version 2.8

      IF ( LINEARIZE_EMATRIX ) then
 
!  set up linearizations of SAB and DAB. All real variables.
!   Thoroughly debugged, November 2005
!  Enhancement # 2, 2/9/16

        DO I = 1, NSTREAMS
          XINV = - ONE / QUAD_STREAMS(I)
          DO J = 1, NSTREAMS
            FAC =  QUAD_HALFWTS(J) * XINV
            DO Q = 1, N_PARAMETERS
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  DP = ZERO ; DM = ZERO
                  DO L = M, NMOMENTS
                    H1P = ZERO ; H1M = ZERO
                    DO O3 = 1, NSTOKES
                      H1P = H1P + PI_XQP(L,I,O1,O3) * &
                                sum( L_OMEGA_GREEK(L,N,O3,1:NSTOKES,Q) * PI_XQP    (L,J,1:NSTOKES,O2) )
                      H1M = H1M + PI_XQP(L,I,O1,O3) * &
                                sum( L_OMEGA_GREEK(L,N,O3,1:NSTOKES,Q) * PI_XQM_PRE(L,J,1:NSTOKES,O2) )
                    ENDDO
                    DP = DP + H1P ; DM = DM + H1M
                  ENDDO
                  L_SAB(I,J,O1,O2,N,Q) = FAC * ( DP + DM )
                  L_DAB(I,J,O1,O2,N,Q) = FAC * ( DP - DM )
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  set up linearized eigenmatrices, saved for all layers. Real variables
!   Thoroughly debugged, November 2005
!  Enhancement # 3, 2/9/16

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO J = 1, NSTREAMS
              JC = NSTOKES*(J-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                DO O2 = 1, NSTOKES
                  JCOL = JC + O2
                  L_EIGENMAT(IROW,JCOL,N,Q) = &
                             sum( L_DAB(I,1:NSTREAMS,O1,1:NSTOKES,N,Q) *   SAB(1:NSTREAMS,J,1:NSTOKES,O2,N) &
                                  + DAB(I,1:NSTREAMS,O1,1:NSTOKES,N)   * L_SAB(1:NSTREAMS,J,1:NSTOKES,O2,N,Q)  )
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Linearization of Real Eigensolutions
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  start eigenvalue loop

        DO K = 1, K_REAL(N)

!  some basic quantities

          AA    = EIGENMASK_R(K)
          KVAL  = KEIGEN(K,N)
          KVAL2 = 2.0D0 * KVAL
          KSQD  = KVAL * KVAL

!   REAL SOLUTIONS ONLY (ASYMTX linearization)
!   ==========================================

          IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Count the number of entries in the eigenvector ( = TC)
!    If the eigenvector is trivial (one entry, TC = 1).....!Rob fix 4/12/12
!      ......Zero the results and skip linearization for this eigenvalue
!      ......Version 2.8, Using flag instead of GOTO statement
!  Enhancement # 4, 2/9/16
!  Enhancement # 5, 2/9/16

            DO_LINEARIZE_ASYMTX = .true.
            TC = count( RITE_EVEC(1:NSTKS_NSTRMS,AA).NE.ZERO )
            IF ( TC.EQ.1 .and. NSTKS_NSTRMS.gt.1 ) THEN
              HVEC1A(1:NSTKS_NSTRMS+1, 1:N_PARAMETERS) = ZERO
              DO_LINEARIZE_ASYMTX = .false.   !! GOTO 4555
            ENDIF

!  Nothing to do also if there is degeneracy
!    NEW CODE by R. Spurr 22 August 2008.
!  Enhancement # 6, 2/9/16

            IF ( DO_LINEARIZE_ASYMTX .and. EIGENDEGEN(K,N) ) THEN
              HVEC1A(1:NSTKS_NSTRMS+1, 1:N_PARAMETERS) = ZERO
              DO_LINEARIZE_ASYMTX = .false.   !! GOTO 4555
            ENDIF

!  Version 2.8, Only proceed to linearize if flag set.

            IF ( DO_LINEARIZE_ASYMTX ) THEN

!  initialise solution matrix HMAT1A (important to do this!)
!    HMAT1A is destroyed upon passage through the LAPACK linear modules
!  Enhancement # 7, 2/9/16

              HMAT1A(1:MAXSTRMSTKS_P1, 1:MAXSTRMSTKS_P1) = ZERO

!  Determine solution matrix HMAT (A in the matrix equation A.x = B)
!  Enhancement # 8, 2/9/16

              NSTMSTKS_P1 = NSTKS_NSTRMS + 1
              HMAT1A(1:NSTKS_NSTRMS,2:NSTKS_NSTRMS+1) = - EIGENMAT_SAVE(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS,N)
              DO I = 1, NSTKS_NSTRMS
                HMAT1A(I,I+1) = HMAT1A(I,I+1) + KSQD
              ENDDO
              HMAT1A(1:NSTKS_NSTRMS,1)             = TWO * KVAL * RITE_EVEC(1:NSTKS_NSTRMS,AA)
              HMAT1A(NSTMSTKS_P1,2:NSTKS_NSTRMS+1) = RITE_EVEC(1:NSTKS_NSTRMS,AA)
              HMAT1A(NSTMSTKS_P1,1)                = ZERO

!        HMAT1A(NSTMSTKS_P1,1) = 1.879d-15   ! should technically be zero

!  solution column vectors (B in the matrix equation A.x = B).
!    (One for each parameter to be varied)
!  Enhancement # 9, 2/9/16

              DO Q = 1, N_PARAMETERS
                DO I = 1, NSTKS_NSTRMS
                  HVEC1A(I,Q) = sum ( L_EIGENMAT(I,1:NSTKS_NSTRMS,N,Q) * RITE_EVEC(1:NSTKS_NSTRMS,AA) )
                ENDDO
                HVEC1A(NSTMSTKS_P1,Q) = ZERO
              ENDDO

!#######################################################################

!  Find any degenerate columns J, and replace the (J,J) entry with
!   a small number. This is a quick fix that ensures that the LAPACK
!   routine does not fall over because of degeneracy.
!    --This problem is only present in the vector case.
!  Fiddle code
!        DO J = 1, NSTKS_NSTRMS
!          NC = 0
!          DO I = 1,  NSTKS_NSTRMS
!           IF ( HMAT1A(I,J+1).EQ.ZERO) NC = NC + 1
!          ENDDO
!          IF ( NC .EQ. NSTKS_NSTRMS ) THEN
!            do I = 1, NSTKS_NSTRMS
!               HMAT1A(J,I+1) = 1.237d-15*DBLE(i)
!            enddo
!            if ( HMAT1A(J+1,J+1).EQ.ZERO)
!     &    HMAT1A(J+1,J+1) = 1.056d-15
!          ENDIF
!        ENDDO
!  ultimate fiddle
!        DO J = 1, NSTMSTKS_P1
!           Do I = 1, NSTMSTKS_P1
!             IF (HMAT1A(I,J).EQ.ZERO )HMAT1A(I,J) = 1.0d-15*(I+J)
!           ENDDO
!        ENDDO
!  solve System for the linearized sum-vector + eigenvalue, together
!     Solve matrix system  using LAPACK modules DGETRF and DEGTRS
!        IF (m.eq.2.and.n.eq.18) THEN
!           WRITE(45, *) m, k, N, NSTMSTKS_P1, LAPACK_INFO
!           WRITE(45, '(17D24.12)')
!     &          ((HMAT1A(j, i), i = 1, NSTMSTKS_P1), j = 1, NSTMSTKS_P1
!c           if(lapack_info.ne.0)STOP '2,18,9'
!        ENDIF

!#######################################################################

!  Resume normal code

              CALL DGETRF  ( NSTMSTKS_P1, NSTMSTKS_P1, HMAT1A, &
                       MAXSTRMSTKS_P1,  LAPACK_IPIV1A, LAPACK_INFO )

!  Exception handling 1, on DGETRF

              IF ( LAPACK_INFO .GT. 0 ) THEN
                WRITE(CI, '(I3)' ) LAPACK_INFO
                WRITE(CN, '(I3)' ) N
                MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
                TRACE   = 'DGETRF call # 1 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
                STATUS  = VLIDORT_SERIOUS
                RETURN
              ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
                WRITE(CI, '(I3)' ) LAPACK_INFO
                WRITE(CN, '(I3)' ) N
                MESSAGE = 'argument i illegal value, for i = '//CI
                TRACE   = 'DGETRF call # 1 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
                STATUS  = VLIDORT_SERIOUS
                RETURN
              ENDIF

              CALL DGETRS ('N', NSTMSTKS_P1, N_PARAMETERS, HMAT1A, &
                     MAXSTRMSTKS_P1, LAPACK_IPIV1A,  HVEC1A, &
                     MAXSTRMSTKS_P1, LAPACK_INFO )

!  Exception handling 1, on DGETRS

              IF ( LAPACK_INFO .LT. 0 ) THEN
                WRITE(CI, '(I3)' ) LAPACK_INFO
                WRITE(CN, '(I3)' ) N
                MESSAGE = 'argument i illegal value, for i = '//CI
                TRACE   = 'DGETRS call # 1 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
                STATUS  = VLIDORT_SERIOUS
                RETURN
              ENDIF

!   End linearization clause (replaces GOTO 4555 statement)

            ENDIF

!  Continuation point if linearized skipped
! 4555   CONTINUE

!  Start loop over varying parameters, to assign answers
!  Enhancement # 10, 2/9/16

            L_KEIGEN(K,N,1:N_PARAMETERS) = HVEC1A(1,1:N_PARAMETERS)
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              L_FWD_SUMVEC_R(I,1:NSTOKES,1:N_PARAMETERS) = HVEC1A(IR+1+1:IR+NSTOKES+1, 1:N_PARAMETERS)
            ENDDO

!  Debug (important).  Tested 28 December 2005.
!        DO Q = 1, N_PARAMETERS
!         IF ( DO_DEBUG_WRITE ) THEN
!          IF (Q.EQ.1.and.m.eq.2.and.n.eq.1) THEN
!           write(37,'(2i3,L3,1p2e20.12)') M,K,EIGENDEGEN(K,N),KEIGEN(K,N),L_KEIGEN(K,N,Q)
!           DO I = 1, NSTREAMS
!            DO O1 = 1, NSTOKES
!                 write(37,'(4i3,1p2e20.12)') m,k,i,o1,RITE_EVEC(I,k ),L_FWD_SUMVEC_R(I,O1,Q)
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDIF
!        ENDDO

!   REAL SOLUTIONS (DGEEV linearization)
!   ====================================

          ELSE

!  initialise

            I_NONTRIVIAL = 0

!  Count the number of entries in the eigenvector ( = TC)
!    If the eigenvector is trivial (one entry, TC = 1).....!Rob fix 4/12/12
!      ......Zero the results and skip linearization for this eigenvalue
!      ......Version 2.8, Using flag instead of GOTO statement
!  Enhancement # 11, 2/9/16

            DO_LINEARIZE_DGEEV = .true.
            TC = count( RITE_EVEC(1:NSTKS_NSTRMS,AA).NE.ZERO )
            IF ( TC.EQ.1  ) THEN
              L_KEIGEN(K, N,       1:N_PARAMETERS) = ZERO
              HVEC1(1:NSTKS_NSTRMS,1:N_PARAMETERS) = ZERO
              DO_LINEARIZE_DGEEV = .false.   !! GOTO 4556
            ENDIF

!  Version 2.8, Only proceed to linearize if flag set.

            IF ( DO_LINEARIZE_DGEEV ) THEN

!  Establish the first non-trivial row
!      * This is necessary to avoid underflow problems
!      * Version 2.8, remove GOTO 3442 construction...............
!          DO JCOL = 1, NSTKS_NSTRMS
!            if ( DABS(RITE_EVEC(JCOL,AA)).GT.1.0D-10 ) THEN
!              TR = JCOL ; GO TO 3442
!            ENDIF
!          ENDDO
! 3442     CONTINUE

              TR = 0 ; JCOL = 0 ; LOOP = .true.
              DO WHILE ( LOOP .and. JCOL .lt. NSTKS_NSTRMS )
                JCOL = JCOl + 1
                LOOP = ( ABS(RITE_EVEC(JCOL,AA)).LT.1.0D-10 )
              ENDDO
              TR = JCOL

!  Use the adjoint method to get L(k) = <X^,L(G)X> / 2k <X^,X>
!    Enhancement # 13, 2/9/16
!    Enhancement # 14, 2/9/16

              SR = sum( LEFT_EVEC(1:NSTKS_NSTRMS,AA) * RITE_EVEC(1:NSTKS_NSTRMS,AA) )
              MODULUS = ONE / KVAL2 / SR
              DO Q = 1, N_PARAMETERS
                DO IROW = 1, NSTKS_NSTRMS
                  HVEC_AUX(IROW,Q) = sum( L_EIGENMAT(IROW,1:NSTKS_NSTRMS,N,Q) * RITE_EVEC(1:NSTKS_NSTRMS,AA) )
                ENDDO
                L_KEIGEN(K,N,Q) = sum( LEFT_EVEC(1:NSTKS_NSTRMS,AA) * HVEC_AUX(1:NSTKS_NSTRMS,Q) ) * MODULUS
              ENDDO

!  Now solve for the linearized eigenvector L(X) from A.L(X) = B
!     Determine matrix HMAT (=A) and vectors HVEC (B)

!  Initialise first, as HMAT is destroyed after LAPACK routines
!    Enhancement # 15, 2/9/16

              HMAT1(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS) = ZERO

!  set up the matrix HMAT
!    Use the linearization of the Eigen-equation for all but one
!    of the rows. For the FIRST row, use the linearization of the
!    normalization condition on the eigenvectors.
 !  Enhancement # 16, 2/9/16
 !  Enhancement # 17, 2/9/16

              DO IROW = 1, NSTKS_NSTRMS
                IF ( IROW.NE.TR ) THEN
                  HMAT1(IROW,1:NSTKS_NSTRMS) = - EIGENMAT_SAVE(IROW,1:NSTKS_NSTRMS,N)
                  HMAT1(IROW,IROW) = HMAT1(IROW,IROW) + KSQD
                ELSE
                  HMAT1(IROW,1:NSTKS_NSTRMS) = RITE_EVEC(1:NSTKS_NSTRMS,AA)
                ENDIF
              ENDDO

!  set up the HVEC:
!    Use the linearization of the Eigen-equation for all but one
!    of the rows. For the FIRST row, use the linearization of the
!    normalization condition on the eigenvectors.

              DO Q = 1, N_PARAMETERS
                I_NONTRIVIAL = 0
                DO IROW = 1, NSTKS_NSTRMS
                  IF ( IROW.NE.TR ) THEN
                    HVEC1(IROW,Q) = HVEC_AUX(IROW,Q) - KVAL2 * L_KEIGEN(K,N,Q) * RITE_EVEC(IROW,AA)
                    IF ( DABS(HVEC1(IROW,Q)).LT.1.0D-12)THEN
                      HVEC1(IROW,Q) = ZERO
                    ELSE
                      I_NONTRIVIAL = I_NONTRIVIAL + 1
                    ENDIF
                  ELSE
                    HVEC1(IROW,Q) = ZERO
                  ENDIF
                ENDDO
              ENDDO

!  Solve matrix system  using LAPACK modules DGETRF and DEGTRS
!  Only do if vector has non-trivial entries

              IF ( I_NONTRIVIAL.GT.0) THEN

!  DGETRF call

                CALL DGETRF  ( NSTKS_NSTRMS, NSTKS_NSTRMS, HMAT1, &
                               MAXSTRMSTKS,  LAPACK_IPIV1, LAPACK_INFO )

!  Exception handling 2, DGTRF

                IF ( LAPACK_INFO .GT. 0 ) THEN
                  WRITE(CI, '(I3)' ) LAPACK_INFO
                  WRITE(CN, '(I3)' ) N
                  MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
                  TRACE   = 'DGETRF call # 2 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
                  STATUS  = VLIDORT_SERIOUS
                  RETURN
                ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
                  WRITE(CI, '(I3)' ) LAPACK_INFO
                  WRITE(CN, '(I3)' ) N
                  MESSAGE = 'argument i illegal value, for i = '//CI
                  TRACE   = 'DGETRF call # 2 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
                  STATUS  = VLIDORT_SERIOUS
                  RETURN
                ENDIF

!  DGETRS Call 

                CALL DGETRS ('N', NSTKS_NSTRMS, N_PARAMETERS, HMAT1, &
                           MAXSTRMSTKS, LAPACK_IPIV1,       HVEC1, &
                           MAXSTRMSTKS, LAPACK_INFO )

!  Exception handling 2, DGETRS

                IF ( LAPACK_INFO .LT. 0 ) THEN
                  WRITE(CI, '(I3)' ) LAPACK_INFO
                  WRITE(CN, '(I3)' ) N
                  MESSAGE = 'argument i illegal value, for i = '//CI
                  TRACE   = 'DGETRS call # 2 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
                  STATUS  = VLIDORT_SERIOUS
                  RETURN
                ENDIF

              ENDIF

!   End linearization clause (replaces GOTO 4556 statement)

            ENDIF

!  Continuation point if linearized skipped
! 4556   CONTINUE

!  Start loop over varying parameters, and assign
!    linearized sum vector = linearized eigensolution
!  Enhancement # 18, 2/9/16

            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              L_FWD_SUMVEC_R(I,1:NSTOKES,1:N_PARAMETERS) = HVEC1(IR+1:IR+NSTOKES,1:N_PARAMETERS)
            ENDDO

!  Debug (important).  Tested 28 December 2005.
!        DO Q = 1, N_PARAMETERS
!         IF ( DO_DEBUG_WRITE ) THEN
!          DO Q = 1, N_PARAMETERS
!          IF (Q.EQ.1.and.m.eq.0.and.n.eq.1) THEN
!           write(37,'(2i3,1p2e20.12)') M,K,KEIGEN(K,N),L_KEIGEN(K,N,Q)
!           DO I = 1, NSTREAMS
!            DO O1 = 1, NSTOKES
!                 write(37,'(4i3,1p2e20.12)') m,k,i,o1,FWD_SUMVEC(I,O1,K),L_FWD_SUMVEC_R(I,O1,Q)
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDIF
!        ENDDO

!  End Clause ASYMTX vs. DGEEV

          ENDIF

!  Resume General Code for Real-Solution linearizations
!  ===================================================

!  Start loop over varying parameters

          DO Q = 1, N_PARAMETERS

!  linearized difference vector
!  Enhancement # 19, 2/9/16

            DO O1 = 1, NSTOKES
              DO I = 1, NSTREAMS
                SR = sum ( - L_SAB(I,1:NSTREAMS,O1,1:NSTOKES,N,Q) *   FWD_SUMVEC  (1:NSTREAMS,1:NSTOKES,K) &
                           -   SAB(I,1:NSTREAMS,O1,1:NSTOKES,N)   * L_FWD_SUMVEC_R(1:NSTREAMS,1:NSTOKES,Q) )
                HR = ( - SR - L_KEIGEN(K,N,Q) * FWD_DIFVEC(I,O1,K) )
                L_FWD_DIFVEC_R(I,O1) = HR / KVAL
              ENDDO
            ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions
!  Enhancement # 20, 2/9/16
!mick fix 9/19/2017 - added DO loop and if block for L_SOLA_XPOS & L_SOLB_XNEG

            LFPD(1:NSTREAMS,1:NSTOKES) = HALF * &
                ( L_FWD_SUMVEC_R(1:NSTREAMS,1:NSTOKES,Q)+L_FWD_DIFVEC_R(1:NSTREAMS,1:NSTOKES) )
            LFND(1:NSTREAMS,1:NSTOKES) = HALF * &
                ( L_FWD_SUMVEC_R(1:NSTREAMS,1:NSTOKES,Q)-L_FWD_DIFVEC_R(1:NSTREAMS,1:NSTOKES) )

            L_SOLA_XPOS(1:NSTREAMS,1:NSTOKES,K,N,Q) = LFPD(1:NSTREAMS,1:NSTOKES)
            L_SOLB_XNEG(1:NSTREAMS,1:NSTOKES,K,N,Q) = LFND(1:NSTREAMS,1:NSTOKES)

            !L_SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS,1:2,K,N,Q) = LFND(1:NSTREAMS,1:2)
            !L_SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS,1:2,K,N,Q) = LFPD(1:NSTREAMS,1:2)
         
            !L_SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS,3:NSTOKES,K,N,Q) = - LFND(1:NSTREAMS,3:NSTOKES)
            !L_SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS,3:NSTOKES,K,N,Q) = - LFPD(1:NSTREAMS,3:NSTOKES)

            DO O1 = 1, NSTOKES
              IF ( O1 .LE. 2 ) THEN
                L_SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS,O1,K,N,Q) = LFND(1:NSTREAMS,O1)
                L_SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS,O1,K,N,Q) = LFPD(1:NSTREAMS,O1)
              ELSE
                L_SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS,3:NSTOKES,K,N,Q) = - LFND(1:NSTREAMS,3:NSTOKES)
                L_SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS,3:NSTOKES,K,N,Q) = - LFPD(1:NSTREAMS,3:NSTOKES)
                EXIT
              ENDIF
            ENDDO

!  End parameter loop

          ENDDO

!  End real eigenvalue loop

        ENDDO

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Linearization of COMPLEX Eigensolutions
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Only if flagged

        IF ( .NOT. DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Offset

          KO1 = K_REAL(N) + 1

!  start loop over eigensolutions

          DO K = 1, K_COMPLEX(N)

!  Bookkeeping indices

            K0 = 2 * ( K - 1 )
            K1 = KO1 + K0
            K2 = K1  + 1

!  Basic quantities
!    -- 1/31/21. Version 2.8.3. KEIGEN_CSQ now has an extra "LAYERS" dimension

            AA  = EIGENMASK_C(K)
            AA1 = AA + 1
            KVAL_CR  = KEIGEN(K1,N)
            KVAL_CI  = KEIGEN(K2,N)
            KVAL2_CR = TWO * KVAL_CR
            KVAL2_CI = TWO * KVAL_CI
            MOD_KVAL = KEIGEN_CSQ(K,N)
            SEPCON_CR =   KEIGEN(K1,N) / MOD_KVAL
            SEPCON_CI = - KEIGEN(K2,N) / MOD_KVAL

!  Establish the zero imaginary part of one eigenvector component
!      * Version 2.8, remove GOTO 3443 construction...............
!          TR = 0
!          DO JCOL = 1, NSTKS_NSTRMS
!            if ( RITE_EVEC(JCOL,AA1).EQ.ZERO ) THEN
!              TR = JCOL + NSTKS_NSTRMS ; GO TO 3443
!            ENDIF
!          ENDDO
! 3443     CONTINUE

            TR = 0 ; JCOL = 0 ; LOOP = .true.
            DO WHILE ( LOOP .and. JCOL .lt. NSTKS_NSTRMS )
              JCOL = JCOl + 1
              LOOP = ( RITE_EVEC(JCOL,AA1).NE.ZERO )
            ENDDO
            TR = JCOL

!  Use the adjoint theory to get L(k) from 2kL(k) = <X^,L(G)X> / <X^,X>
!      MUST USE THE COMPLEX CONJUGATE OF THE LEFT EIGENVECTOR.
!  --------------- Adjoint  norm
!  Enhancement # 21, 2/9/16

            SCR = sum( + LEFT_EVEC(1:NSTKS_NSTRMS,AA)  * RITE_EVEC(1:NSTKS_NSTRMS,AA) &
                       + LEFT_EVEC(1:NSTKS_NSTRMS,AA1) * RITE_EVEC(1:NSTKS_NSTRMS,AA1) )
            SCI = sum( - LEFT_EVEC(1:NSTKS_NSTRMS,AA1) * RITE_EVEC(1:NSTKS_NSTRMS,AA) &
                       + LEFT_EVEC(1:NSTKS_NSTRMS,AA)  * RITE_EVEC(1:NSTKS_NSTRMS,AA1) )
            ADJN_CR(1:NSTKS_NSTRMS) = + KVAL2_CR * RITE_EVEC(1:NSTKS_NSTRMS,AA) &
                                      - KVAL2_CI * RITE_EVEC(1:NSTKS_NSTRMS,AA1)
            ADJN_CI(1:NSTKS_NSTRMS) = + KVAL2_CI * RITE_EVEC(1:NSTKS_NSTRMS,AA) &
                                      + KVAL2_CR * RITE_EVEC(1:NSTKS_NSTRMS,AA1)

            DENOM_CR = KVAL2_CR * SCR - KVAL2_CI * SCI
            DENOM_CI = KVAL2_CR * SCI + KVAL2_CI * SCR
            MODULUS  = ONE / ( DENOM_CR*DENOM_CR + DENOM_CI * DENOM_CI )

!  ---------------- Adjoint numerator, and final computation
!  Enhancement # 22, 2/9/16

            DO Q = 1, N_PARAMETERS
              DO IROW = 1, NSTKS_NSTRMS
                HVEC_CR(IROW,Q) = sum( L_EIGENMAT(IROW,1:NSTKS_NSTRMS,N,Q) * RITE_EVEC(1:NSTKS_NSTRMS,AA) )
                HVEC_CI(IROW,Q) = sum( L_EIGENMAT(IROW,1:NSTKS_NSTRMS,N,Q) * RITE_EVEC(1:NSTKS_NSTRMS,AA1) )
              ENDDO
              HCR = sum( + LEFT_EVEC(1:NSTKS_NSTRMS,AA)  * HVEC_CR(1:NSTKS_NSTRMS,Q) &
                         + LEFT_EVEC(1:NSTKS_NSTRMS,AA1) * HVEC_CI(1:NSTKS_NSTRMS,Q) )
              HCI = sum( - LEFT_EVEC(1:NSTKS_NSTRMS,AA1) * HVEC_CR(1:NSTKS_NSTRMS,Q) &
                         + LEFT_EVEC(1:NSTKS_NSTRMS,AA)  * HVEC_CI(1:NSTKS_NSTRMS,Q) )
              L_KEIGEN(K1,N,Q) = ( HCR*DENOM_CR + HCI*DENOM_CI )*MODULUS
              L_KEIGEN(K2,N,Q) = ( HCI*DENOM_CR - HCR*DENOM_CI )*MODULUS
            ENDDO

!  Now solve for the linearized eigenvector L(X) from A.L(X) = B
!     Determine matrix HMAT (=A) and vectors HVEC (B)

!  ---------- Initialise first, as HMAT is destroyed after LAPACK routin
!  Enhancement # 23, 2/9/16

            HMAT2(1:NSTKS_NSTRMS_2,1:NSTKS_NSTRMS_2) = ZERO

!  ---------- set up matrix for Full solution
!  Enhancement # 24, 2/9/16

            HMAT2(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS)   = - EIGENMAT_SAVE(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS,N)
            HMAT2(NSTKS_NSTRMS+1:NSTKS_NSTRMS+NSTKS_NSTRMS,NSTKS_NSTRMS+1:NSTKS_NSTRMS+NSTKS_NSTRMS) = & 
                                                     - EIGENMAT_SAVE(1:NSTKS_NSTRMS,1:NSTKS_NSTRMS,N)
            DO I = 1, NSTKS_NSTRMS
              I1 = I + NSTKS_NSTRMS
              HMAT2(I,I)   = HMAT2(I,I)   + REAL_KSQ(AA)
              HMAT2(I1,I1) = HMAT2(I1,I1) + REAL_KSQ(AA)
              HMAT2(I,I1)  = - IMAG_KSQ(AA)
              HMAT2(I1,I)  = + IMAG_KSQ(AA)
            ENDDO

!  ---------- set up the vectors for full solution
!  Enhancement # 25, 2/9/16

            DO I = 1, NSTKS_NSTRMS
              I1 = I + NSTKS_NSTRMS
              HVEC2(I ,1:N_PARAMETERS) = HVEC_CR(I,1:N_PARAMETERS) - &
                               ( + L_KEIGEN(K1,N,1:N_PARAMETERS) * ADJN_CR(I) &
                                 - L_KEIGEN(K2,N,1:N_PARAMETERS) * ADJN_CI(I) )
              HVEC2(I1,1:N_PARAMETERS) = HVEC_CI(I,1:N_PARAMETERS) - & 
                               ( + L_KEIGEN(K1,N,1:N_PARAMETERS) * ADJN_CI(I) &
                                 + L_KEIGEN(K2,N,1:N_PARAMETERS) * ADJN_CR(I))
            ENDDO

!  ---------- Strike out TR row/column to get reduced matrix and vectors
!  Enhancement # 26, 2/9/16

!mick mod
        !NSTKS_NSTRMS_2  = 2 * NSTKS_NSTRMS
            NSTKS_NSTRMS_21 = 2 * NSTKS_NSTRMS - 1
            HMATRED(1:TR-1, 1:TR-1)              = HMAT2(1:TR-1, 1:TR-1)
            HMATRED(1:TR-1, TR:NSTKS_NSTRMS_2-1) = HMAT2(1:TR-1, TR+1:NSTKS_NSTRMS_2)
            HVECRED(1:TR-1, 1:N_PARAMETERS)      = HVEC2(1:TR-1, 1:N_PARAMETERS)

            HMATRED(TR:NSTKS_NSTRMS_2-1, 1:TR-1)              = HMAT2(TR+1:NSTKS_NSTRMS_2, 1:TR-1)
            HMATRED(TR:NSTKS_NSTRMS_2-1, TR:NSTKS_NSTRMS_2-1) = HMAT2(TR+1:NSTKS_NSTRMS_2, TR+1:NSTKS_NSTRMS_2)
            HVECRED(TR:NSTKS_NSTRMS_2-1, 1:N_PARAMETERS)      = HVEC2(TR+1:NSTKS_NSTRMS_2, 1:N_PARAMETERS)

!  ----------  Replace the last row by linearized normalization constrai
!                ( IF (TR.LT.NSTKS_NSTRMS_2) THEN  possible exception...
!                Should not be necessary to flag this case )
!!  Enhancement # 27, 2/9/16

            ILAST =  NSTKS_NSTRMS_21
            HMATRED(ILAST, 1:NSTKS_NSTRMS)      = RITE_EVEC(1:NSTKS_NSTRMS,AA)
            HMATRED(ILAST, NSTKS_NSTRMS+1:TR-1) = RITE_EVEC(1:TR-1-NSTKS_NSTRMS, AA1)
            HMATRED(ILAST, TR:NSTKS_NSTRMS_2-1) = RITE_EVEC(TR+1-NSTKS_NSTRMS:NSTKS_NSTRMS_2-NSTKS_NSTRMS, AA1)
            HVECRED(ILAST, 1:N_PARAMETERS)      = ZERO

!  ---------- Solve matrix system using LAPACK modules DGETRF and DEGTRS

            CALL DGETRF  ( NSTKS_NSTRMS_21, NSTKS_NSTRMS_21, HMATRED, &
                           MAXSTRMSTKS_21,  LAPACK_IPIV21,   LAPACK_INFO )

!  Exception handling 3

            IF ( LAPACK_INFO .GT. 0 ) THEN
              WRITE(CI, '(I3)' ) LAPACK_INFO
              WRITE(CN, '(I3)' ) N
              MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
              TRACE   = 'DGETRF call # 3 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
              STATUS  = VLIDORT_SERIOUS
              RETURN
            ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
              WRITE(CI, '(I3)' ) LAPACK_INFO
              WRITE(CN, '(I3)' ) N
              MESSAGE = 'argument i illegal value, for i = '//CI
              TRACE   = 'DGETRF call # 3 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
              STATUS  = VLIDORT_SERIOUS
              RETURN
            ENDIF

            CALL DGETRS ('N', NSTKS_NSTRMS_21, N_PARAMETERS, HMATRED, &
                         MAXSTRMSTKS_21, LAPACK_IPIV21,      HVECRED, &
                         MAXSTRMSTKS_21, LAPACK_INFO )

!  Exception handling 3

            IF ( LAPACK_INFO .LT. 0 ) THEN
              WRITE(CI, '(I3)' ) LAPACK_INFO
              WRITE(CN, '(I3)' ) N
              MESSAGE = 'argument i illegal value, for i = '//CI
              TRACE   = 'DGETRS call # 3 in VLIDORT_L_QHOM_SOLUTION,'// ' layer # '//CN
              STATUS  = VLIDORT_SERIOUS
              RETURN
            ENDIF

!  Start loop over varying parameters

            DO Q = 1, N_PARAMETERS

!  linearized sum vector = linearized eigensolution

              DO I = 1, NSTREAMS
                IR = NSTOKES*(I-1)
                DO O1 = 1, NSTOKES
                  IROW  = IR + O1
                  IROW1 = IROW + NSTKS_NSTRMS
                  L_FWD_SUMVEC_CR(I,O1) = HVECRED(IROW,Q)
                  IF ( IROW1.LT.TR)THEN
                    L_FWD_SUMVEC_CI(I,O1) = HVECRED(IROW1,Q)
                  ELSE IF (IROW1.EQ.TR) THEN
                    L_FWD_SUMVEC_CI(I,O1) = ZERO
                  ELSE
                    L_FWD_SUMVEC_CI(I,O1) = HVECRED(IROW1-1,Q)
                  ENDIF
                ENDDO
              ENDDO

!  Debug code (important)
!         IF (DO_DEBUG_WRITE) THEN
!          IF ( Q.EQ.1.AND.M.EQ.0.AND.N.EQ.1) THEN
!           WRITE(91,'(2I3,1P4E20.10)')
!     &            M,K,KEIGEN(K1,N),KEIGEN(K2,N),
!     &            L_KEIGEN(K1,N,Q),L_KEIGEN(K2,N,Q)
!           DO I = 1, NSTREAMS
!            DO O1 = 1, NSTOKES
!              WRITE(91,'(4I3,1P4E20.10)')
!     &        M,K,I,O1,FWD_SUMVEC(I,O1,K1),FWD_SUMVEC(I,O1,K2),
!     &         L_FWD_SUMVEC_CR(I,O1),L_FWD_SUMVEC_CI(I,O1)
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDIF

!  Determine the linearized difference vector
!  Enhancement # 28, 2/9/16

              DO O1 = 1, NSTOKES
                DO I = 1, NSTREAMS
                  SCR = sum( - L_SAB(I,1:NSTREAMS,O1,1:NSTOKES,N,Q) *   FWD_SUMVEC(1:NSTREAMS,1:NSTOKES,K1) &
                             -   SAB(I,1:NSTREAMS,O1,1:NSTOKES,N)   * L_FWD_SUMVEC_CR(1:NSTREAMS,1:NSTOKES) )
                  SCI = sum( - L_SAB(I,1:NSTREAMS,O1,1:NSTOKES,N,Q) *   FWD_SUMVEC(1:NSTREAMS,1:NSTOKES,K2) &
                            -   SAB(I,1:NSTREAMS,O1,1:NSTOKES,N)   * L_FWD_SUMVEC_CI(1:NSTREAMS,1:NSTOKES) )
                  HCR = - SCR - L_KEIGEN(K1,N,Q) * FWD_DIFVEC(I,O1,K1) &
                              + L_KEIGEN(K2,N,Q) * FWD_DIFVEC(I,O1,K2)
                  HCI = - SCI - L_KEIGEN(K2,N,Q) * FWD_DIFVEC(I,O1,K1) &
                              - L_KEIGEN(K1,N,Q) * FWD_DIFVEC(I,O1,K2)
                  L_FWD_DIFVEC_CR(I,O1) = HCR*SEPCON_CR - HCI*SEPCON_CI
                  L_FWD_DIFVEC_CI(I,O1) = HCI*SEPCON_CR + HCR*SEPCON_CI
                ENDDO
              ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions
!  Enhancement # 29, 2/9/16

              LFND1(1:NSTREAMS,1:NSTOKES) = HALF * &
                    (L_FWD_SUMVEC_CR(1:NSTREAMS,1:NSTOKES)-L_FWD_DIFVEC_CR(1:NSTREAMS,1:NSTOKES))
              LFPD1(1:NSTREAMS,1:NSTOKES) = HALF * &
                    (L_FWD_SUMVEC_CR(1:NSTREAMS,1:NSTOKES)+L_FWD_DIFVEC_CR(1:NSTREAMS,1:NSTOKES))
              LFND2(1:NSTREAMS,1:NSTOKES) = HALF * &
                    (L_FWD_SUMVEC_CI(1:NSTREAMS,1:NSTOKES)-L_FWD_DIFVEC_CI(1:NSTREAMS,1:NSTOKES))
              LFPD2(1:NSTREAMS,1:NSTOKES) = HALF * &
                    (L_FWD_SUMVEC_CI(1:NSTREAMS,1:NSTOKES)+L_FWD_DIFVEC_CI(1:NSTREAMS,1:NSTOKES))

              L_SOLA_XPOS(1:NSTREAMS,1:NSTOKES,K1,N,Q)  = LFPD1(1:NSTREAMS,1:NSTOKES)
              L_SOLB_XNEG(1:NSTREAMS,1:NSTOKES,K1,N,Q)  = LFND1(1:NSTREAMS,1:NSTOKES)
              L_SOLA_XPOS(1:NSTREAMS,1:NSTOKES,K2,N,Q)  = LFPD2(1:NSTREAMS,1:NSTOKES)
              L_SOLB_XNEG(1:NSTREAMS,1:NSTOKES,K2,N,Q)  = LFND2(1:NSTREAMS,1:NSTOKES)

              L_SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS,1:2,K1,N,Q)  = LFND1(1:NSTREAMS,1:2)
              L_SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS,1:2,K1,N,Q)  = LFPD1(1:NSTREAMS,1:2)
              L_SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS,1:2,K2,N,Q)  = LFND2(1:NSTREAMS,1:2)
              L_SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS,1:2,K2,N,Q)  = LFPD2(1:NSTREAMS,1:2)

              L_SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS,3:NSTOKES,K1,N,Q)  = - LFND1(1:NSTREAMS,3:NSTOKES)
              L_SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS,3:NSTOKES,K1,N,Q)  = - LFPD1(1:NSTREAMS,3:NSTOKES)
              L_SOLA_XPOS(NSTREAMS+1:NSTREAMS+NSTREAMS,3:NSTOKES,K2,N,Q)  = - LFND2(1:NSTREAMS,3:NSTOKES)
              L_SOLB_XNEG(NSTREAMS+1:NSTREAMS+NSTREAMS,3:NSTOKES,K2,N,Q)  = - LFPD2(1:NSTREAMS,3:NSTOKES)

!  End parameter loop

            ENDDO

!  End complex eigenvalue loop

          ENDDO

!  End clause for not using the real eigensolver

        ENDIF

!  End Eigensolver linearization clause

      ENDIF

!  Version 2.8, GOTO Statement 3456 removed.
!  control point for avoiding eigensolver
! 3456 CONTINUE

!  Linearized Eigenstream transmittance factors for whole layer
!  ------------------------------------------------------------

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then linearized transmittances are
!  linearized discrete ordinate transmittances.

      IF ( DO_SOLUTION_SAVING .AND. .NOT.DO_LAYER_SCAT(M,N) ) THEN

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              L_T_DELT_EIGEN(IROW,N,Q) = L_T_DELT_DISORDS(I,N,Q)
            ENDDO
          ENDDO
        ENDDO

!  Otherwise compute them as normal

      ELSE

!  Real
!  Enhancement # 30, 2/9/16

        DO K = 1, K_REAL(N)
          TBAS = - T_DELT_EIGEN(K,N) * DELTAU_VERT(N)
          L_T_DELT_EIGEN(K,N,1:N_PARAMETERS) = TBAS * &
                    ( L_KEIGEN(K,N,1:N_PARAMETERS) + KEIGEN(K,N)*L_DELTAU_VERT(1:N_PARAMETERS,N) )
        ENDDO

!  Complex

        DO K = 1, K_COMPLEX(N)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          TBAS_CR = - T_DELT_EIGEN(K1,N) *  DELTAU_VERT(N)
          TBAS_CI = - T_DELT_EIGEN(K2,N) *  DELTAU_VERT(N)
          DO Q = 1, N_PARAMETERS
            LCR = L_KEIGEN(K1,N,Q) + KEIGEN(K1,N) * L_DELTAU_VERT(Q,N)
            LCI = L_KEIGEN(K2,N,Q) + KEIGEN(K2,N) * L_DELTAU_VERT(Q,N)
            L_T_DELT_EIGEN(K1,N,Q) = TBAS_CR * LCR - TBAS_CI * LCI
            L_T_DELT_EIGEN(K2,N,Q) = TBAS_CI * LCR + TBAS_CR * LCI
          ENDDO
        ENDDO

      ENDIF

!  debug
!      if ( m.eq.0.and.n.eq.20 ) then
!       DO K = 1, min(6,K_REAL(N))
!         WRITE(*,'(2i4,1p6e24.12)')m,k,L_T_DELT_EIGEN(K,N,1)
!       ENDDO
!      ENDIF

!  Eigenstream transmittance factors for partial layers
!  ----------------------------------------------------

!  Done irrespective of solution saving
!    Loop over user optical depths

      DO UT = 1, N_PARTLAYERS
        NA  = PARTLAYERS_LAYERIDX(UT)

!  If the layer of occurrence is the given layer

        IF ( NA .EQ. N ) THEN

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

          IF ( DO_SOLUTION_SAVING .AND. .NOT.DO_LAYER_SCAT(M,N) ) THEN

!  Enhancement # 31, 2/9/16
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                L_T_UTDN_EIGEN(IROW,UT,1:N_PARAMETERS)=L_T_DISORDS_UTDN(I,UT,1:N_PARAMETERS)
                L_T_UTUP_EIGEN(IROW,UT,1:N_PARAMETERS)=L_T_DISORDS_UTUP(I,UT,1:N_PARAMETERS)
              ENDDO
            ENDDO

!  Otherwise, Compute the Eigenstream transmittance factors

          ELSE

!  Partial layer optical depths

            TAU_DN = PARTAU_VERT(UT)
            TAU_UP = DELTAU_VERT(N) - TAU_DN

!  Real eigenvalues
!  Enhancement # 32, 2/9/16

            DO K = 1, K_REAL(N)
              VAR(1:N_PARAMETERS) = L_KEIGEN(K,N,1:N_PARAMETERS) + KEIGEN(K,N) * L_DELTAU_VERT(1:N_PARAMETERS,N)
              L_T_UTDN_EIGEN(K,UT,1:N_PARAMETERS) = ( - TAU_DN * T_UTDN_EIGEN(K,UT) ) * VAR(1:N_PARAMETERS)
              L_T_UTUP_EIGEN(K,UT,1:N_PARAMETERS) = ( - TAU_UP * T_UTUP_EIGEN(K,UT) ) * VAR(1:N_PARAMETERS)
            ENDDO

!  Complex eigenvalues
!    Bug fixed 16 December 2005

            DO K = 1, K_COMPLEX(N)
              K0 = 2 * ( K - 1 )
              K1 = KO1 + K0
              K2 = K1  + 1
              TBASUP_CR = - TAU_UP * T_UTUP_EIGEN(K1,UT)
              TBASUP_CI = - TAU_UP * T_UTUP_EIGEN(K2,UT)
              TBASDN_CR = - TAU_DN * T_UTDN_EIGEN(K1,UT)
              TBASDN_CI = - TAU_DN * T_UTDN_EIGEN(K2,UT)
              DO Q = 1, N_PARAMETERS
                LCR = L_KEIGEN(K1,N,Q) + KEIGEN(K1,N) * L_DELTAU_VERT(Q,N)
                LCI = L_KEIGEN(K2,N,Q) + KEIGEN(K2,N) * L_DELTAU_VERT(Q,N)
                L_T_UTDN_EIGEN(K1,UT,Q) = TBASDN_CR*LCR - TBASDN_CI*LCI
                L_T_UTDN_EIGEN(K2,UT,Q) = TBASDN_CI*LCR + TBASDN_CR*LCI
                L_T_UTUP_EIGEN(K1,UT,Q) = TBASUP_CR*LCR - TBASUP_CI*LCI
                L_T_UTUP_EIGEN(K2,UT,Q) = TBASUP_CI*LCR + TBASUP_CR*LCI
              ENDDO
            ENDDO

!  End clause for using Discrete ordinate values or not

          ENDIF

!  End loop over partial layer depths

        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_QHOM_SOLUTION

!

      SUBROUTINE VLIDORT_L_QHOM_NORMS &
          ( FOURIER, NSTOKES, NSTREAMS, NLAYERS, QUAD_STRMWTS,      & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, K_REAL, K_COMPLEX,  & ! Input
            SOLA_XPOS, SOLB_XNEG, L_SOLA_XPOS, L_SOLB_XNEG,         & ! Input
            L_NORM_SAVED )                                            ! Output

!  1/31/21. Version 2.8.3. Introduced this new subroutine, based on LIDORT code.
!   ==> Eigenproblem, linearized solution norms for Green's function
!   ==> only for Real-value solutions (NSTOKES < 4). Complex solutions are a PLACEHOLDER
!   ==> No efficiency applied to the coding

!  module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXEVALUES, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, ZERO, TWO

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  FOURIER component (debug only)

      INTEGER  , intent(in)  :: FOURIER

!  Number of stokes elements, discrete-ordinate streams streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NSTOKES
      INTEGER  , intent(in)  :: NLAYERS

!  Quadrature

      DOUBLE PRECISION, intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Linearization control

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Eigensolution obookkeeping

      INTEGER, INTENT (IN) ::    K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::    K_COMPLEX ( MAXLAYERS )

!  Eigenvector solutions

      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Linearized Eigenvector solutions

      DOUBLE PRECISION, INTENT (IN) ::  L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  outputs
!  -------

!  Linearized norms for the Green function solution

      DOUBLE PRECISION, INTENT(out) :: L_NORM_SAVED ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER          :: N, K, O1, J, KO1, K0, K1, K2, Q, M
      DOUBLE PRECISION :: T1, T2, NORM
      DOUBLE PRECISION :: T1_CR, T1_CI, T2_CR, T2_CI, NORM_CR, NORM_CI

!  debug

      M = FOURIER

!  For all layers, save the norms

      DO N = 1, NLAYERS
        KO1 = K_REAL(N) + 1
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)

!  Loop over real eigensolutions

            DO K = 1, K_REAL(N)
              NORM = ZERO
              DO J = 1, NSTREAMS
                T1 = zero ; T2 = zero
                DO O1 = 1, NSTOKES
                  T1 = T1 + L_SOLA_XPOS(J,O1,K,N,Q) * SOLA_XPOS(J,O1,K,N)
                  T2 = T2 + L_SOLB_XNEG(J,O1,K,N,Q) * SOLB_XNEG(J,O1,K,N)         
                ENDDO
                NORM = NORM + QUAD_STRMWTS(J) * ( T1 - T2 )
              ENDDO
              L_NORM_SAVED(K,N,Q) = TWO * NORM
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
!                  T1_CR = T1_CR + L_ADJA_XPOS(J,O1,K1,N,Q) * ADJA_XPOS(J,O1,K1,N) &
!                                - L_ADJA_XPOS(J,O1,K2,N,Q) * ADJA_XPOS(J,O1,K2,N)
!                  T2_CR = T2_CR + L_ADJB_XNEG(J,O1,K1,N,Q) * ADJB_XNEG(J,O1,K1,N) &
!                                - L_ADJB_XNEG(J,O1,K2,N,Q) * ADJB_XNEG(J,O1,K2,N)
!                  T1_CI = T1_CI + L_ADJA_XPOS(J,O1,K1,N,Q) * ADJA_XPOS(J,O1,K2,N) &
!                                + L_ADJA_XPOS(J,O1,K2,N,Q) * ADJA_XPOS(J,O1,K1,N)
!                  T2_CI = T2_CI + L_ADJB_XNEG(J,O1,K1,N,Q) * ADJB_XNEG(J,O1,K1,N) &
!                                + L_ADJB_XNEG(J,O1,K2,N,Q) * ADJB_XNEG(J,O1,K2,N)
                ENDDO
                NORM_CR = NORM_CR + QUAD_STRMWTS(J) * ( T1_CR - T2_CR )
                NORM_CI = NORM_CI + QUAD_STRMWTS(J) * ( T1_CI - T2_CI )
              ENDDO
              L_NORM_SAVED(K1,N,Q) = TWO * NORM_CR
              L_NORM_SAVED(K2,N,Q) = TWO * NORM_CI
            ENDDO

!  End variations and Layer Loop

          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_QHOM_NORMS

!

      SUBROUTINE VLIDORT_L_UHOM_SOLUTION ( &
        DO_UPWELLING, DO_DNWELLING, DO_DEBUG_WRITE,              & ! Input flags
        GIVEN_LAYER, FOURIER, NSTOKES, NSTREAMS, N_USER_STREAMS, & ! Input numbers
        DO_VARY, N_PARAMETERS, NMOMENTS,                         & ! Input Lin-control, numbers
        QUAD_HALFWTS, USER_SECANTS, DO_LAYER_SCAT,               & ! Input streams/scat
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                  & ! Input bookkeeping
        OMEGA_GREEK, L_OMEGA_GREEK, PI_XQP, PI_XQM_PRE,          & ! Input optical + PI-quad
        PI_XUP, PI_XUM, PI_XUM_POST, PI_XUP_PRE,                 & ! Input PI-User
        K_REAL, K_COMPLEX, KEIGEN, HELPSTOKES,                   & ! Input eigensolution
        UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,              & ! Input User solutions
        ZETA_M, ZETA_P, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG,      & ! Input Zeta + linearized Homog
        L_UHOM_DNDN, L_UHOM_DNUP, L_UHOM_UPDN, L_UHOM_UPUP,      & ! OUTPUT Linearized user solutions
        L_ZETA_M, L_ZETA_P )                                       ! OUTPUT linearized Zeta

!    -- 1/31/21. Version 2.8.3.  LOCAL_UM_START (argument dropped)

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS, MAX_USER_STREAMS, &
                                 MAXMOMENTS, MAXEVALUES, MAXSTREAMS_2, MAXSTRMSTKS, ZERO, ONE, TWO

      IMPLICIT NONE

!  INPUTS
!  ======

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_DEBUG_WRITE

!  Input linearization control

      LOGICAL, INTENT (IN) ::          DO_VARY
      INTEGER, INTENT (IN) ::          N_PARAMETERS

!  Input numbers

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          NMOMENTS

!  Input streams, scattering

      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS  ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Input bookkeeping

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCAT ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )

!  Input optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK   ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Input Pi matrices

      DOUBLE PRECISION, INTENT (IN) :: PI_XQP     ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) :: PI_XUP      ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM      ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM_POST ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUP_PRE  ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Eigenproblem solution

      INTEGER, INTENT (INOUT) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: KEIGEN     ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: HELPSTOKES ( 0:MAXMOMENTS, MAXEVALUES, MAXSTOKES )

!  User solutions

      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Linearized Eigensolution
 
      DOUBLE PRECISION, INTENT (INOUT) :: L_KEIGEN    ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )


!  OUTPUT
!  ======

!  Linearized user solutions

      DOUBLE PRECISION, INTENT (INOUT) :: L_UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Zeta

      DOUBLE PRECISION, INTENT (INOUT) :: L_ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: L_HELPSTOKES_R  ( MAXSTOKES )
      DOUBLE PRECISION :: L_HELPSTOKES_CR ( MAXSTOKES )
      DOUBLE PRECISION :: L_HELPSTOKES_CI ( MAXSTOKES )

      DOUBLE PRECISION :: L_GAUX_R  ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: L_GAUX_CR ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: L_GAUX_CI ( 0:MAXMOMENTS, MAXSTOKES )

      LOGICAL ::          DOLAYER_PP
      INTEGER ::          UM, J, L, N, M, Q, O1
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPS, SP, SM, SPSCR, SPSCI, SPCR, SPCI, SMCR, SMCI
      DOUBLE PRECISION :: KISQ, RHO_P_CR, RHO_M_CR, MP, MM

!  Enhancement variables, 2/9/16, 5/9/16 using symbolic dimensions
!      DOUBLE PRECISION :: LMP(N_PARAMETERS), LMM(N_PARAMETERS)
      DOUBLE PRECISION :: LMP(MAX_ATMOSWFS), LMM(MAX_ATMOSWFS)

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  general post-processing flag

      DOLAYER_PP = ( STERM_LAYERMASK_UP(N) .OR. STERM_LAYERMASK_DN(N) )

!  Linearized Zeta constants
!  =========================

      IF ( DOLAYER_PP ) THEN

!  start User angle loop

        DO UM = 1, N_USER_STREAMS
          SM = USER_SECANTS(UM)

!  Linearized Real Zeta constants (always required)
!  Enhancement # 33, 2/9/16

          DO K = 1, K_REAL(N)
             L_ZETA_P(K,UM,N,1:N_PARAMETERS) = - L_KEIGEN(K,N,1:N_PARAMETERS) * ZETA_P(K,UM,N) * ZETA_P(K,UM,N)
             L_ZETA_M(K,UM,N,1:N_PARAMETERS) = + L_KEIGEN(K,N,1:N_PARAMETERS) * ZETA_M(K,UM,N) * ZETA_M(K,UM,N)
          ENDDO

!  linearization of Complex Zeta constants
!  Enhancement # 34, 2/9/16

          IF ( K_COMPLEX(N) .GT. 0 ) THEN
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              KISQ = KEIGEN(K2,N)*KEIGEN(K2,N)
              RHO_P_CR = SM + KEIGEN(K1,N)
              RHO_M_CR = SM - KEIGEN(K1,N)
              MP = ONE / ( RHO_P_CR*RHO_P_CR + KISQ )
              MM = ONE / ( RHO_M_CR*RHO_M_CR + KISQ )
              LMP(1:N_PARAMETERS)    = TWO * &
                   ( L_KEIGEN(K2,N,1:N_PARAMETERS)*KEIGEN(K2,N) + RHO_P_CR * L_KEIGEN(K1,N,1:N_PARAMETERS) )
              LMM(1:N_PARAMETERS)    = TWO * &
                   ( L_KEIGEN(K2,N,1:N_PARAMETERS)*KEIGEN(K2,N) - RHO_M_CR * L_KEIGEN(K1,N,1:N_PARAMETERS) )
              L_ZETA_P(K1,UM,N,1:N_PARAMETERS) = - MP * &
                   ( ZETA_P(K1,UM,N)*LMP(1:N_PARAMETERS) - L_KEIGEN(K1,N,1:N_PARAMETERS) )
              L_ZETA_P(K2,UM,N,1:N_PARAMETERS) = - MP * &
                  ( ZETA_P(K2,UM,N)*LMP(1:N_PARAMETERS) + L_KEIGEN(K2,N,1:N_PARAMETERS) )
              L_ZETA_M(K1,UM,N,1:N_PARAMETERS) = - MM * &
                  ( ZETA_M(K1,UM,N)*LMM(1:N_PARAMETERS) + L_KEIGEN(K1,N,1:N_PARAMETERS) )
              L_ZETA_M(K2,UM,N,1:N_PARAMETERS) = - MM * &
                  ( ZETA_M(K2,UM,N)*LMM(1:N_PARAMETERS) - L_KEIGEN(K2,N,1:N_PARAMETERS) )
            ENDDO
          ENDIF

!  End zeta linearization

        ENDDO
      ENDIF

!  Eigenvector interpolation to user-defined angles
!  ================================================

!  Only do this if both flags are set,
!    Or if there is a solution for scattering in this layer

      IF ( .NOT.DOLAYER_PP .OR. .NOT. DO_VARY .OR. .NOT. DO_LAYER_SCAT(M,N) ) RETURN

!  For each parameter

      DO Q = 1, N_PARAMETERS

!  User defined solutions (Real)
!  -----------------------------

        DO K = 1, K_REAL(N)

!  For each moment, do inner sum over computational angles
!  for the positive and negative linearized eigenvectors.
!  Linearized the product with OMEGA_GREEK ---> L_GAUX
!  Enhancement # 35, 2/9/16

          DO L = M, NMOMENTS
            DO O1 = 1, NSTOKES
              SPS = ZERO
              DO J = 1, NSTREAMS
                SP = sum( QUAD_HALFWTS(J)*PI_XQP    (L,J,O1,1:NSTOKES)*L_SOLA_XPOS(J,1:NSTOKES,K,N,Q) )
                SM = sum( QUAD_HALFWTS(J)*PI_XQM_PRE(L,J,O1,1:NSTOKES)*L_SOLB_XNEG(J,1:NSTOKES,K,N,Q) )
                SPS = SPS + SP + SM
              ENDDO
              L_HELPSTOKES_R(O1) = SPS
            ENDDO
            DO O1 = 1, NSTOKES
              L_GAUX_R(L,O1) = sum( OMEGA_GREEK(L,N,O1,1:NSTOKES)   * L_HELPSTOKES_R(1:NSTOKES) &
                                + L_OMEGA_GREEK(L,N,O1,1:NSTOKES,Q) *   HELPSTOKES(L,K,1:NSTOKES) )
            ENDDO
          ENDDO

!  Now sum over all harmonic contributions (downwelling)
!  Enhancement # 36, 2/9/16

          IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N)) THEN
            DO O1 = 1, NSTOKES
              DO UM = 1, N_USER_STREAMS
                L_UHOM_DNDN(UM,O1,K,N,Q) = sum( PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES)      * L_GAUX_R(M:NMOMENTS,1:NSTOKES) )
                L_UHOM_DNUP(UM,O1,K,N,Q) = sum( PI_XUM_POST(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_R(M:NMOMENTS,1:NSTOKES) )
              ENDDO
            ENDDO
          ENDIF

!  Now sum over all harmonic contributions (Upwelling)
!  Enhancement # 37, 2/9/16

          IF ( DO_UPWELLING .AND. STERM_LAYERMASK_UP(N)) THEN
            DO O1 = 1, NSTOKES
              DO UM = 1, N_USER_STREAMS
                L_UHOM_UPDN(UM,O1,K,N,Q) = sum( PI_XUM    (M:NMOMENTS,UM,O1,1:NSTOKES)  * L_GAUX_R(M:NMOMENTS,1:NSTOKES) )
                L_UHOM_UPUP(UM,O1,K,N,Q) = sum( PI_XUP_PRE(M:NMOMENTS,UM,O1,1:NSTOKES)  * L_GAUX_R(M:NMOMENTS,1:NSTOKES) )
              ENDDO
            ENDDO
          ENDIF

!  Debug (important)

          IF ( DO_DEBUG_WRITE.AND.Q.EQ.2) THEN
            DO UM = 1, N_USER_STREAMS
               DO O1 = 1, NSTOKES
                 WRITE(95,'(4I3,1P8E20.10)')M,K,UM,O1, &
                   L_UHOM_UPUP(UM,O1,K,N,Q), L_UHOM_UPDN(UM,O1,K,N,Q), &
                   L_UHOM_DNUP(UM,O1,K,N,Q), L_UHOM_DNDN(UM,O1,K,N,Q), &
                   UHOM_UPUP(UM,O1,K,N), UHOM_UPDN(UM,O1,K,N), &
                   UHOM_DNUP(UM,O1,K,N), UHOM_DNDN(UM,O1,K,N)
              ENDDO
            ENDDO
          ENDIF

!  end loop over real eigensolutions

        ENDDO

!  User defined solutions (Complex)
!  -------------------------------

!  offset

        KO1 = K_REAL(N) + 1

!  start loop over complex eigenvalues

        DO K = 1, K_COMPLEX(N)

          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors.
!    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec
!  Enhancement # 38, 2/9/16

          DO L = M, NMOMENTS
            DO O1 = 1, NSTOKES
              SPSCR = ZERO
              SPSCI = ZERO
              DO J = 1, NSTREAMS
                SPCR = sum( PI_XQP(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J) * L_SOLA_XPOS(J,1:NSTOKES,K1,N,Q) )
                SPCI = sum( PI_XQP(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J) * L_SOLA_XPOS(J,1:NSTOKES,K2,N,Q) )
                SMCR = sum( PI_XQM_PRE(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J) * L_SOLB_XNEG(J,1:NSTOKES,K1,N,Q) )
                SMCI = sum( PI_XQM_PRE(L,J,O1,1:NSTOKES)*QUAD_HALFWTS(J) * L_SOLB_XNEG(J,1:NSTOKES,K2,N,Q) )
                SPSCR = SPSCR + SPCR + SMCR
                SPSCI = SPSCI + SPCI + SMCI
              ENDDO
              L_HELPSTOKES_CR(O1) = SPSCR
              L_HELPSTOKES_CI(O1) = SPSCI
            ENDDO
            DO O1 = 1, NSTOKES
              L_GAUX_CR(L,O1) = sum( OMEGA_GREEK(L,N,O1,1:NSTOKES)     * L_HELPSTOKES_CR(1:NSTOKES) &
                                   + L_OMEGA_GREEK(L,N,O1,1:NSTOKES,Q) *   HELPSTOKES(L,K1,1:NSTOKES) )
              L_GAUX_CI(L,O1) = sum( OMEGA_GREEK(L,N,O1,1:NSTOKES)     * L_HELPSTOKES_CI(1:NSTOKES) &
                                   + L_OMEGA_GREEK(L,N,O1,1:NSTOKES,Q) *   HELPSTOKES(L,K2,1:NSTOKES) )
            ENDDO
          ENDDO

!  Now sum over all harmonic contributions (downwelling)
!  Enhancement # 39, 2/9/16

          IF ( DO_DNWELLING ) THEN
            DO O1 = 1, NSTOKES
              DO UM = 1, N_USER_STREAMS
                L_UHOM_DNDN(UM,O1,K1,N,Q) = sum( PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_CR(M:NMOMENTS,1:NSTOKES) )
                L_UHOM_DNUP(UM,O1,K1,N,Q) = sum( PI_XUM_POST(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_CR(M:NMOMENTS,1:NSTOKES) )
                L_UHOM_DNDN(UM,O1,K2,N,Q) = sum( PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_CI(M:NMOMENTS,1:NSTOKES) )
                L_UHOM_DNUP(UM,O1,K2,N,Q) = sum( PI_XUM_POST(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_CI(M:NMOMENTS,1:NSTOKES) )
              ENDDO
            ENDDO
          ENDIF

!  Now sum over all harmonic contributions (upwelling)
!  Enhancement # 40, 2/9/16

          IF ( DO_UPWELLING ) THEN
            DO O1 = 1, NSTOKES
              DO UM = 1, N_USER_STREAMS
                L_UHOM_UPDN(UM,O1,K1,N,Q) = sum( PI_XUM(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_CR(M:NMOMENTS,1:NSTOKES) )
                L_UHOM_UPUP(UM,O1,K1,N,Q) = sum( PI_XUP_PRE(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_CR(M:NMOMENTS,1:NSTOKES) )
                L_UHOM_UPDN(UM,O1,K2,N,Q) = sum( PI_XUM(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_CI(M:NMOMENTS,1:NSTOKES) )
                L_UHOM_UPUP(UM,O1,K2,N,Q) = sum( PI_XUP_PRE(M:NMOMENTS,UM,O1,1:NSTOKES) * L_GAUX_CI(M:NMOMENTS,1:NSTOKES) )
              ENDDO
            ENDDO
          ENDIF

!  Debug (important)

          IF ( DO_DEBUG_WRITE.AND.Q.EQ.2) THEN
            DO UM = 1, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                 WRITE(97,'(4I3,1P8E20.10)')M,K,UM,O1, &
                   L_UHOM_UPUP(UM,O1,K2,N,Q), L_UHOM_UPDN(UM,O1,K2,N,Q), &
                   L_UHOM_DNUP(UM,O1,K2,N,Q), L_UHOM_DNDN(UM,O1,K2,N,Q), &
                   UHOM_UPUP(UM,O1,K2,N), UHOM_UPDN(UM,O1,K2,N), &
                   UHOM_DNUP(UM,O1,K2,N), UHOM_DNDN(UM,O1,K2,N)
              ENDDO
            ENDDO
          ENDIF

!  end loop over complex eigensolutions

        ENDDO

!  End loop over parameters

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_UHOM_SOLUTION

!

      SUBROUTINE L_HMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, NLAYERS, N_PARTLAYERS,    & !  Flags/Taylor/Numbers
        N_USER_STREAMS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                 & ! Numbers + Lin-control
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,  PARTLAYERS_LAYERIDX,       & ! Output control
        USER_SECANTS, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,              & ! sstreams/Optical
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                           & ! Transmittances User
        T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                           & ! Transmittances Homog.
        K_REAL, K_COMPLEX, HSINGO, HMULT_1, HMULT_2, ZETA_M, ZETA_P,        & ! Multipliers and Zetas
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,                     & ! Linearized Transmittances User
        L_T_DELT_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,                     & ! LinearizedTransmittances Homog.
        L_KEIGEN, L_ZETA_M, L_ZETA_P,                                       & ! Linearized Zeta/eigenvalues
        L_HMULT_1, L_HMULT_2, L_UT_HMULT_UU,                                & ! OUTPUT - Linearized Multipliers 
        L_UT_HMULT_UD, L_UT_HMULT_DU, L_UT_HMULT_DD )                         ! OUTPUT - Linearized Multipliers 

!  1/31/21. Version 2.8.3.  LOCAL_UM_START (argument dropped). UM loops now start with 1

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAX_USER_STREAMS, MAX_ATMOSWFS, MAXMOMENTS, MAXEVALUES,     &
                                 MAXSTREAMS_2, MAXSTRMSTKS, ZERO, ONE

      IMPLICIT NONE

!  INPUTS
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING

!  Numbers

      INTEGER, INTENT (IN) ::          TAYLOR_ORDER  ! 2p7 Taylor-order argument added
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Linearization control

      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  streams

      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Output Level/Partial control

      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT   ( MAX_PARTLAYERS )  ! @@@ Rob Fix 5/10/13 added
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigensolution variables

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Multipliers

      LOGICAL, INTENT (IN) ::          HSINGO  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: ZETA_M  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: ZETA_P  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Linearized transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized solution variables

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )

!  OUTPUT
!  ======

!  Linearized whole-layer multipliers

      DOUBLE PRECISION, INTENT (OUT) :: L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized partial-layer multipliers

      DOUBLE PRECISION, INTENT (OUT) :: L_UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

!  Control

      LOGICAL ::          DO_RTSOL_VARY ( MAXLAYERS )
      INTEGER ::          NPARAMS_VARY  ( MAXLAYERS )

!  integers

      INTEGER ::          N, UT, UM, K, Q, KO1, K0, K1, K2

!  real variables

      DOUBLE PRECISION :: UDEL, ZDEL_R, L_UDEL, L_ZDEL_R
      DOUBLE PRECISION :: L_T2, L_T1, HOM1, HOM2, SM, EPS, DELTA, L_DELTA, L_MULT
      DOUBLE PRECISION :: UX_UP, L_UX_UP(MAX_ATMOSWFS)
      DOUBLE PRECISION :: UX_DN, L_UX_DN(MAX_ATMOSWFS)
      DOUBLE PRECISION :: ZX_UP_R, ZX_DN_R
      DOUBLE PRECISION :: THETA_DN_R, THETA_UP_R
      DOUBLE PRECISION :: L_ZX_UP_R, L_ZX_DN_R
      DOUBLE PRECISION :: L_THETA_DN_R, L_THETA_UP_R

!  variables for the complex multiplier linearization

      DOUBLE PRECISION :: ZDEL_CR, ZDEL_CI
      DOUBLE PRECISION :: THETA_1_CR, L_THETA_1_CR
      DOUBLE PRECISION :: THETA_2_CR, L_THETA_2_CR
      DOUBLE PRECISION :: THETA_1_CI, L_THETA_1_CI
      DOUBLE PRECISION :: THETA_2_CI, L_THETA_2_CI
      DOUBLE PRECISION :: THETA_DN_CR, THETA_UP_CR
      DOUBLE PRECISION :: THETA_DN_CI, THETA_UP_CI
      DOUBLE PRECISION :: ZX_UP_CR, ZX_DN_CR
      DOUBLE PRECISION :: ZX_UP_CI, ZX_DN_CI
      DOUBLE PRECISION :: L_THETA_DN_CR, L_THETA_UP_CR
      DOUBLE PRECISION :: L_THETA_DN_CI, L_THETA_UP_CI
      DOUBLE PRECISION :: L_ZDEL_CR, L_ZX_UP_CR, L_ZX_DN_CR
      DOUBLE PRECISION :: L_ZDEL_CI, L_ZX_UP_CI, L_ZX_DN_CI

!  Local control
!  Enhancement # 41, 2/9/16

      DO_RTSOL_VARY(1:NLAYERS) = LAYER_VARY_FLAG  (1:NLAYERS)
      NPARAMS_VARY (1:NLAYERS) = LAYER_VARY_NUMBER(1:NLAYERS)

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!   Only done if layers are flagged

      DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
          IF ( DO_RTSOL_VARY(N) ) THEN
            DO UM = 1, N_USER_STREAMS

!  setup

              UDEL = T_DELT_USERM(N,UM)
              SM   = USER_SECANTS(UM)

!  linearized Real multipliers
!    Chain Rule: similar to scalar code.

              DO K = 1, K_REAL(N)
                ZDEL_R    = T_DELT_EIGEN(K,N)
                DO Q = 1, NPARAMS_VARY(N)
                  L_ZDEL_R = L_T_DELT_EIGEN(K,N,Q)
                  L_UDEL   = L_T_DELT_USERM(N,UM,Q)
                  IF ( HSINGO(K,UM,N) ) THEN
                    EPS      = ZETA_M(K,UM,N) ; DELTA   = DELTAU_VERT(N)
                    L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                    CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, &
                                           L_KEIGEN(K,N,Q), ZERO, UDEL, SM, L_MULT )
                    L_HMULT_1(K,UM,N,Q) = SM * L_MULT
!  10 May 13       CALL TAYLOR_SERIES_L_1 ( EPS, DELTAU_VERT(N), SM, ONE, UDEL, &
!  10 May 13                                L_KEIGEN(K,N,Q), L_DEL, L_HMULT_1(K,UM,N,Q) )
!   Old            HOM1 = L_DELTAU_VERT(Q,N)*(ONE-SM) - &
!   Old                  HALF * L_KEIGEN(K,N,Q) * DELTAU_VERT(N)
!   Old            L_HMULT_1(K,UM,N,Q) = HMULT_1(K,UM,N) * HOM1
                  ELSE
                    L_T1 = L_ZDEL_R - L_UDEL
                    HOM1 =   L_KEIGEN(K,N,Q)*HMULT_1(K,UM,N) + SM*L_T1
                    L_HMULT_1(K,UM,N,Q) = ZETA_M(K,UM,N) * HOM1
                  ENDIF
                  L_T2 = - ZDEL_R * L_UDEL - L_ZDEL_R * UDEL
                  HOM2 = - L_KEIGEN(K,N,Q)*HMULT_2(K,UM,N) + SM*L_T2
                  L_HMULT_2(K,UM,N,Q) = ZETA_P(K,UM,N) * HOM2
                ENDDO
              ENDDO

!  Linearized Complex multipliers
!    Chain rule using real and complex parts.

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
                DO Q = 1, NPARAMS_VARY(N)
                  L_ZDEL_CR = L_T_DELT_EIGEN(K1,N,Q)
                  L_ZDEL_CI = L_T_DELT_EIGEN(K2,N,Q)
                  L_UDEL    = L_T_DELT_USERM(N,UM,Q)
                  L_THETA_2_CR = - L_ZDEL_CR * UDEL - ZDEL_CR * L_UDEL
                  L_THETA_2_CI = - L_ZDEL_CI * UDEL - ZDEL_CI * L_UDEL
                  L_THETA_1_CR = L_ZDEL_CR - L_UDEL
                  L_THETA_1_CI = L_ZDEL_CI

                  L_HMULT_1(K1,UM,N,Q) = SM * &
                   (   THETA_1_CR * L_ZETA_M(K1,UM,N,Q) + &
                     L_THETA_1_CR *   ZETA_M(K1,UM,N)   - &
                       THETA_1_CI * L_ZETA_M(K2,UM,N,Q) - &
                     L_THETA_1_CI *   ZETA_M(K2,UM,N)   )
                  L_HMULT_1(K2,UM,N,Q) = SM * &
                   (   THETA_1_CR * L_ZETA_M(K2,UM,N,Q) + &
                     L_THETA_1_CR *   ZETA_M(K2,UM,N)   + &
                       THETA_1_CI * L_ZETA_M(K1,UM,N,Q) + &
                     L_THETA_1_CI *   ZETA_M(K1,UM,N)   )

                  L_HMULT_2(K1,UM,N,Q) = SM * &
                   (   THETA_2_CR * L_ZETA_P(K1,UM,N,Q) + &
                     L_THETA_2_CR *   ZETA_P(K1,UM,N)   - &
                       THETA_2_CI * L_ZETA_P(K2,UM,N,Q) - &
                     L_THETA_2_CI *   ZETA_P(K2,UM,N)   )
                  L_HMULT_2(K2,UM,N,Q) = SM * &
                   (   THETA_2_CR * L_ZETA_P(K2,UM,N,Q) + &
                     L_THETA_2_CR *   ZETA_P(K2,UM,N)   + &
                       THETA_2_CI * L_ZETA_P(K1,UM,N,Q) + &
                     L_THETA_2_CI *   ZETA_P(K1,UM,N)   )

                ENDDO
              ENDDO

!  End loops over user angles and layers

            ENDDO
          ENDIF
        ENDIF
      ENDDO

!  partial layer multipliers
!  -------------------------

      DO UT = 1, N_PARTLAYERS
        N  = PARTLAYERS_LAYERIDX(UT)

!  UPWELLING
!  ---------

        IF ( DO_UPWELLING.AND.DO_RTSOL_VARY(N) ) THEN

!  start code with loop over user angles

          DO UM = 1, N_USER_STREAMS

!  set-up
!  Enhancement # 42, 2/9/16

            UX_UP = T_UTUP_USERM(UT,UM)
            SM    = USER_SECANTS(UM)
            L_UX_UP(1:NPARAMS_VARY(N)) = L_T_UTUP_USERM(UT,UM,1:NPARAMS_VARY(N))

!  Linearization of Real multipliers
!    Use Chain rule. Similar to scalar case.

            DO K = 1, K_REAL(N)
              ZDEL_R     = T_DELT_EIGEN(K,N)
              ZX_UP_R    = T_UTUP_EIGEN(K,UT)
              ZX_DN_R    = T_UTDN_EIGEN(K,UT)
              THETA_DN_R = ZX_DN_R - ZDEL_R * UX_UP
              THETA_UP_R = ZX_UP_R - UX_UP
              DO Q = 1, NPARAMS_VARY(N)
                L_ZDEL_R     = L_T_DELT_EIGEN(K,N,Q)
                L_ZX_UP_R    = L_T_UTUP_EIGEN(K,UT,Q)
                L_ZX_DN_R    = L_T_UTDN_EIGEN(K,UT,Q)
                L_THETA_DN_R = L_ZX_DN_R - L_ZDEL_R * UX_UP - ZDEL_R * L_UX_UP(Q)
                L_THETA_UP_R = L_ZX_UP_R - L_UX_UP(Q)
                L_UT_HMULT_UD(K,UM,UT,Q) = SM * &
                    ( L_THETA_DN_R *   ZETA_P(K,UM,N) + THETA_DN_R * L_ZETA_P(K,UM,N,Q) )
                IF ( HSINGO(K,UM,N) )THEN
                  EPS   = ZETA_M(K,UM,N) ; DELTA   = DELTAU_VERT(N) - PARTAU_VERT(UT)
                  L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA
                  CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, &
                                           L_KEIGEN(K,N,Q), ZERO, UX_UP, SM, L_MULT )
                  L_UT_HMULT_UU(K,UM,UT,Q)= SM * L_MULT
                ELSE
                  L_UT_HMULT_UU(K,UM,UT,Q) = SM * &
                     ( L_THETA_UP_R *   ZETA_M(K,UM,N) + THETA_UP_R * L_ZETA_M(K,UM,N,Q) )
                ENDIF
              ENDDO
            ENDDO

!  Linearization of Complex multipliers, Use Chain rule
!  Enhancement # 43, 2/9/16

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

              L_UT_HMULT_UD(K1,UM,UT,1:NPARAMS_VARY(N)) = SM * &
                  ( L_T_UTDN_EIGEN(K1,UT,1:NPARAMS_VARY(N)) &
                  - L_T_DELT_EIGEN(K1,N,1:NPARAMS_VARY(N)) * UX_UP - ZDEL_CR * L_UX_UP(1:NPARAMS_VARY(N)) *   ZETA_P(K1,UM,N)   + &
                      THETA_DN_CR * L_ZETA_P(K1,UM,N,1:NPARAMS_VARY(N)) - &
                    L_T_UTDN_EIGEN(K2,UT,1:NPARAMS_VARY(N)) &
                  - L_T_DELT_EIGEN(K2,N,1:NPARAMS_VARY(N)) * UX_UP - ZDEL_CI * L_UX_UP(1:NPARAMS_VARY(N)) *   ZETA_P(K2,UM,N)   - &
                      THETA_DN_CI * L_ZETA_P(K2,UM,N,1:NPARAMS_VARY(N)) )
              L_UT_HMULT_UD(K2,UM,UT,1:NPARAMS_VARY(N)) = SM * &
                  ( L_T_UTDN_EIGEN(K1,UT,1:NPARAMS_VARY(N)) &
                  - L_T_DELT_EIGEN(K1,N,1:NPARAMS_VARY(N)) * UX_UP - ZDEL_CR * L_UX_UP(1:NPARAMS_VARY(N)) *   ZETA_P(K2,UM,N)   + &
                      THETA_DN_CR * L_ZETA_P(K2,UM,N,1:NPARAMS_VARY(N)) + &
                    L_T_UTDN_EIGEN(K2,UT,1:NPARAMS_VARY(N)) &
                  - L_T_DELT_EIGEN(K2,N,1:NPARAMS_VARY(N)) * UX_UP - ZDEL_CI * L_UX_UP(1:NPARAMS_VARY(N)) *   ZETA_P(K1,UM,N)   + &
                      THETA_DN_CI * L_ZETA_P(K1,UM,N,1:NPARAMS_VARY(N)) )
              L_UT_HMULT_UU(K1,UM,UT,1:NPARAMS_VARY(N)) = SM * &
                  ( L_T_UTUP_EIGEN(K1,UT,1:NPARAMS_VARY(N)) - L_UX_UP(1:NPARAMS_VARY(N)) *   ZETA_M(K1,UM,N)   + &
                      THETA_UP_CR * L_ZETA_M(K1,UM,N,1:NPARAMS_VARY(N)) - &
                    L_T_UTUP_EIGEN(K2,UT,1:NPARAMS_VARY(N)) *   ZETA_M(K2,UM,N)   - &
                      THETA_UP_CI * L_ZETA_M(K2,UM,N,1:NPARAMS_VARY(N)) )
              L_UT_HMULT_UU(K2,UM,UT,1:NPARAMS_VARY(N)) = SM * &
                  ( L_T_UTUP_EIGEN(K1,UT,1:NPARAMS_VARY(N)) - L_UX_UP(1:NPARAMS_VARY(N)) *   ZETA_M(K2,UM,N)   + &
                      THETA_UP_CR * L_ZETA_M(K2,UM,N,1:NPARAMS_VARY(N)) + &
                    L_T_UTUP_EIGEN(K2,UT,1:NPARAMS_VARY(N)) *   ZETA_M(K1,UM,N)   + &
                      THETA_UP_CI * L_ZETA_M(K1,UM,N,1:NPARAMS_VARY(N)) )
            ENDDO

!  Finish loop over user angles

          ENDDO

!  End upwelling clause

        ENDIF

!  DOWNWELLING
!  -----------

        IF ( DO_DNWELLING.AND.DO_RTSOL_VARY(N) ) THEN

!  start code with loop over user angles

          DO UM = 1, N_USER_STREAMS

!  set-up
!  Enhancement # 44, 2/9/16


            UX_DN = T_UTDN_USERM(UT,UM)
            SM    = USER_SECANTS(UM)
            L_UX_DN(1:NPARAMS_VARY(N)) = L_T_UTDN_USERM(UT,UM,1:NPARAMS_VARY(N))

!  Linearization of Real multipliers
!    Use Chain rule. Similar to scalar case.

            DO K = 1, K_REAL(N)
              ZDEL_R     = T_DELT_EIGEN(K,N)
              ZX_UP_R    = T_UTUP_EIGEN(K,UT)
              ZX_DN_R    = T_UTDN_EIGEN(K,UT)
              THETA_DN_R = ZX_DN_R - UX_DN
              THETA_UP_R = ZX_UP_R - ZDEL_R * UX_DN
              DO Q = 1, NPARAMS_VARY(N)
                L_ZDEL_R     = L_T_DELT_EIGEN(K,N,Q)
                L_ZX_UP_R    = L_T_UTUP_EIGEN(K,UT,Q)
                L_ZX_DN_R    = L_T_UTDN_EIGEN(K,UT,Q)
                L_THETA_UP_R = L_ZX_UP_R - L_ZDEL_R * UX_DN - ZDEL_R * L_UX_DN(Q)
                L_THETA_DN_R = L_ZX_DN_R - L_UX_DN(Q)
                IF ( HSINGO(K,UM,N) )THEN
                  EPS      = ZETA_M(K,UM,N) ; DELTA   = PARTAU_VERT(UT)
                  L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA  ! Input is Single normalized (partau)
                  CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, &
                                           L_KEIGEN(K,N,Q), ZERO, UX_DN, SM, L_MULT )
                  L_UT_HMULT_DD(K,UM,UT,Q) = SM * L_MULT
                ELSE
                  L_UT_HMULT_DD(K,UM,UT,Q) = SM * &
                       ( L_THETA_DN_R *   ZETA_M(K,UM,N) + THETA_DN_R * L_ZETA_M(K,UM,N,Q) )
                ENDIF
                L_UT_HMULT_DU(K,UM,UT,Q) = SM * &
                     ( L_THETA_UP_R *   ZETA_P(K,UM,N) + THETA_UP_R * L_ZETA_P(K,UM,N,Q) )
              ENDDO
            ENDDO

!  Linearization of Complex multipliers
!    Use Chain rule

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

              THETA_UP_CR = ZX_UP_CR - ZDEL_CR * UX_DN
              THETA_UP_CI = ZX_UP_CI - ZDEL_CI * UX_DN
              THETA_DN_CR = ZX_DN_CR - UX_DN
              THETA_DN_CI = ZX_DN_CI

              DO Q = 1, NPARAMS_VARY(N)
                L_ZDEL_CR     = L_T_DELT_EIGEN(K1,N,Q)
                L_ZDEL_CI     = L_T_DELT_EIGEN(K2,N,Q)
                L_ZX_UP_CR    = L_T_UTUP_EIGEN(K1,UT,Q)
                L_ZX_UP_CI    = L_T_UTUP_EIGEN(K2,UT,Q)
                L_ZX_DN_CR    = L_T_UTDN_EIGEN(K1,UT,Q)
                L_ZX_DN_CI    = L_T_UTDN_EIGEN(K2,UT,Q)

                L_THETA_UP_CR = L_ZX_UP_CR - L_ZDEL_CR * UX_DN - ZDEL_CR * L_UX_DN(Q)
                L_THETA_UP_CI = L_ZX_UP_CI - L_ZDEL_CI * UX_DN - ZDEL_CI * L_UX_DN(Q)
                L_THETA_DN_CR = L_ZX_DN_CR - L_UX_DN(Q)
                L_THETA_DN_CI = L_ZX_DN_CI

                L_UT_HMULT_DD(K1,UM,UT,Q) = SM * &
                   ( L_THETA_DN_CR *   ZETA_M(K1,UM,N)   + &
                       THETA_DN_CR * L_ZETA_M(K1,UM,N,Q) - &
                     L_THETA_DN_CI *   ZETA_M(K2,UM,N)   - &
                       THETA_DN_CI * L_ZETA_M(K2,UM,N,Q) )
                L_UT_HMULT_DD(K2,UM,UT,Q) = SM * &
                   ( L_THETA_DN_CR *   ZETA_M(K2,UM,N)   + &
                       THETA_DN_CR * L_ZETA_M(K2,UM,N,Q) + &
                     L_THETA_DN_CI *   ZETA_M(K1,UM,N)   + &
                       THETA_DN_CI * L_ZETA_M(K1,UM,N,Q) )

                L_UT_HMULT_DU(K1,UM,UT,Q) = SM * &
                   ( L_THETA_UP_CR *   ZETA_P(K1,UM,N)   + &
                       THETA_UP_CR * L_ZETA_P(K1,UM,N,Q) - &
                     L_THETA_UP_CI *   ZETA_P(K2,UM,N)   - &
                       THETA_UP_CI * L_ZETA_P(K2,UM,N,Q) )
                L_UT_HMULT_DU(K2,UM,UT,Q) = SM * &
                   ( L_THETA_UP_CR *   ZETA_P(K2,UM,N)   + &
                       THETA_UP_CR * L_ZETA_P(K2,UM,N,Q) + &
                     L_THETA_UP_CI *   ZETA_P(K1,UM,N)   + &
                       THETA_UP_CI * L_ZETA_P(K1,UM,N,Q) )

              ENDDO
            ENDDO

!  Finish loop over user angles

          ENDDO

!  end of downwelling clause

        ENDIF

!  End off-grid loop

      ENDDO

!  Finish

! @@@ Rob Fix 5/10/13 - This is a good place to pause for debug
!      pause'End L_HMULT_MASTER'

      RETURN
      END SUBROUTINE L_HMULT_MASTER

!

      SUBROUTINE VLIDORT_L_GBEAM_SOLUTION &
         ( FOURIER, IBEAM, GIVEN_LAYER, NSTOKES, NSTREAMS,        & ! input Numbers
           NSTKS_NSTRMS, NMOMENTS, DOVARY, N_PARAMETERS,          & ! input Numbers
           DO_LAYSCAT, FLUX_FACTOR, DFLUX, QUAD_WTS, BEAM_CUTOFF, & ! Input Bookkeeping
           PI_XQP, PI_XQM, PI_X0P, PI_XQM_POST,                   & ! Input GSF,
           L_OMEGA_GREEK, K_REAL, SOLA_XPOS, SOLB_XNEG,           & ! Input Optical/Eigensolutions
           NORM_SAVED, ATERM_SAVE, BTERM_SAVE, DMI, DPI,          & ! Input Green's function stuff
           L_NORM_SAVED, L_SOLA_XPOS, L_SOLB_XNEG,                & ! Input linearized solutions
           L_ATERM_SAVE, L_BTERM_SAVE )                             ! Output

!  1/31/21. New routine: Linearization part of Green's function solution (non-multipliers)
!    -- Real-valued Eigensolutions only (Nstokes = 1 or 3)

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXBEAMS, MAXSTREAMS, MAXMOMENTS, MAXSTREAMS_2,  &
                                 MAXEVALUES, MAXSTRMSTKS, MAX_ATMOSWFS, ZERO, ONE, PI4


      IMPLICIT NONE

!  subroutine input arguments
!  ==========================

!  Given layer index, Beam index and Fourier number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Numbers

      INTEGER  , intent(in)  :: NSTREAMS, NSTKS_NSTRMS
      INTEGER  , intent(in)  :: NMOMENTS, NSTOKES

!  Linearization control

      LOGICAL  , intent(in)  :: DOVARY
      INTEGER  , intent(in)  :: N_PARAMETERS

!  Flux and quadrature

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WTS ( MAXSTREAMS )

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYSCAT (0:MAXMOMENTS, MAXLAYERS)

!  Last layer to include Particular integral solution

      INTEGER  , INTENT (IN) :: BEAM_CUTOFF ( MAXBEAMS )

!  PI Matrices, Linearized OMEGA_GREEK

      DOUBLE PRECISION, INTENT (IN) ::  L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Eigensolution bookkeeping.
!   -- real-valued eigensolutions only (K_COMPLEX not required here (at the moment))

      INTEGER, INTENT (in) ::          K_REAL ( MAXLAYERS )

!  Eigensolutions and normalization factors

      DOUBLE PRECISION, intent(in) ::  SOLA_XPOS (MAXSTREAMS_2,MAXSTOKES,MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, intent(in) ::  SOLB_XNEG (MAXSTREAMS_2,MAXSTOKES,MAXEVALUES,MAXLAYERS)

!  Linearized Eigenvector solutions and normalizations

      DOUBLE PRECISION, intent(in) ::  L_SOLA_XPOS (MAXSTREAMS_2,MAXSTOKES,MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent(in) ::  L_SOLB_XNEG (MAXSTREAMS_2,MAXSTOKES,MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Input Green's function normalizations, and Aterm/Bterm vectors

      DOUBLE PRECISION, intent(in) ::  NORM_SAVED ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, intent(in) ::  ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, intent(in) ::  BTERM_SAVE ( MAXEVALUES, MAXLAYERS )

!  Input Greens function DMI and DPI variables

      DOUBLE PRECISION, intent(in) :: DMI (MAXSTREAMS,MAXSTOKES), DPI (MAXSTREAMS,MAXSTOKES)

!  Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION, intent(in) :: L_NORM_SAVED (MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  subroutine output arguments
!  ===========================

!  Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION, intent(inout) :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent(inout) :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  linearizations of component variables

      DOUBLE PRECISION :: L_DMI(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_DPI(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)

!  help variables

      LOGICAL          :: DO_FIRST
      INTEGER          :: K, O1, O2, L, I, M, N, Q, IB
      DOUBLE PRECISION :: GSUM, F1, JK
      DOUBLE PRECISION :: HELP1 (NSTREAMS), HELP2 (NSTREAMS)
      DOUBLE PRECISION :: HELP_QFUNC(0:MAXMOMENTS,MAXSTOKES)

!  initialise indices

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

!  Check existence
!  ===============

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. BEAM_CUTOFF(IB) ) .AND. DO_LAYSCAT(M,N)

!  If there is nothing varying or if there is no solution, zero the output and exit

      IF ( .NOT.DOVARY .OR. .NOT.DO_FIRST ) THEN
        DO Q = 1, N_PARAMETERS
          L_ATERM_SAVE(1:NSTKS_NSTRMS,N,Q) = ZERO
          L_BTERM_SAVE(1:NSTKS_NSTRMS,N,Q) = ZERO
        ENDDO
        RETURN
      ENDIF

!  Form quantities independent of optical depth
!  ============================================

!  Set up linearizations of help arrays (independent of eigenvector)
!mick fix 1/13/2021 - moved defining of L_DPI & L_DMI from inside first set of loops below to after it
!                     (so all needed elements of the help array HELP_QFUNC are defined first)

      DO Q = 1, N_PARAMETERS
        DO O1 = 1, NSTOKES
          DO L = M, NMOMENTS
            GSUM = ZERO
            DO O2 = 1, NSTOKES
              GSUM = GSUM + L_OMEGA_GREEK(L,N,O1,O2,Q) * sum( PI_X0P(L,IBEAM,N,O2,1:NSTOKES) * DFLUX(1:NSTOKES) )
            ENDDO
            HELP_QFUNC(L,O1) = GSUM
          ENDDO
        ENDDO
        DO O1 = 1, NSTOKES
          DO I = 1, NSTREAMS
            L_DPI(I,O1,Q) = F1 * sum( PI_XQP(M:NMOMENTS,I,O1,1:NSTOKES)      * HELP_QFUNC(M:NMOMENTS,1:NSTOKES) )
            L_DMI(I,O1,Q) = F1 * sum( PI_XQM_POST(M:NMOMENTS,I,O1,1:NSTOKES) * HELP_QFUNC(M:NMOMENTS,1:NSTOKES) )
          ENDDO
        ENDDO
      ENDDO

!  linearize quantities independent of TAU (L_ATERM_SAVE, L_BTERM_SAVE)
!  ==> Real-valued solutions only

      DO Q = 1, N_PARAMETERS
        DO K = 1, K_REAL(N)

!  A-term derivatives
!  ------------------

          DO I = 1, NSTREAMS
            HELP1(I) = dot_product(L_DPI(I,1:nstokes,Q),  SOLA_XPOS(I,1:nstokes,K,N)) &
                     + dot_product(DPI  (I,1:nstokes)  ,L_SOLA_XPOS(I,1:nstokes,K,N,Q))
            HELP2(I) = dot_product(L_DMI(I,1:nstokes,Q),  SOLB_XNEG(I,1:nstokes,K,N)) &
                     + dot_product(DMI  (I,1:nstokes)  ,L_SOLB_XNEG(I,1:nstokes,K,N,Q))
          ENDDO

!  Rob Chages 1/15/21. Make L_ATERM a normal derivative

          JK = DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS), (HELP1+HELP2) )
          L_ATERM_SAVE(K,N,Q) = ( JK - ATERM_SAVE(K,N) * L_NORM_SAVED(K,N,Q) ) / NORM_SAVED(K,N)

!  here is the old logarithmic derivative
!          L_ATERM_SAVE(K,N,Q) = ( ( DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS), (HELP1+HELP2) ) / ATERM_SAVE(K,N) ) &
!                                   - L_NORM_SAVED(K,N,Q) ) / NORM_SAVED(K,N)

!  B-term derivatives
!  ------------------

          DO I = 1, NSTREAMS
            HELP1(I) = dot_product(L_DMI(I,1:nstokes,Q),  SOLA_XPOS(I,1:nstokes,K,N)) &
                     + dot_product(DMI  (I,1:nstokes)  ,L_SOLA_XPOS(I,1:nstokes,K,N,Q))
            HELP2(I) = dot_product(L_DPI(I,1:nstokes,Q),  SOLB_XNEG(I,1:nstokes,K,N)) &
                     + dot_product(DPI  (I,1:nstokes)  ,L_SOLB_XNEG(I,1:nstokes,K,N,Q))
          ENDDO

!  Rob Chages 1/15/21. Make L_BTERM a normal derivative

          JK = DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS), (HELP1+HELP2) )
          L_BTERM_SAVE(K,N,Q) = ( JK - BTERM_SAVE(K,N) * L_NORM_SAVED(K,N,Q) ) / NORM_SAVED(K,N)

!  here is the old logarithmic derivative
!          L_BTERM_SAVE(K,N,Q) = ( ( DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS), (HELP1+HELP2) ) / BTERM_SAVE(K,N) ) &
!                                   - L_NORM_SAVED(K,N,Q) ) / NORM_SAVED(K,N)

!  end eigenvalue and parameter loops

        ENDDO
      ENDDO

!  Complex Valued solution - PLACEHOLDER
!  -------------------------------------

!  CODE NOT READY

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_GBEAM_SOLUTION

!

      SUBROUTINE VLIDORT_L_GUSER_SOLUTION &
          ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,  & ! Input Flags
            LAYER, FOURIER, IBEAM, NSTOKES, N_USER_STREAMS,       & ! Input numbers
            NMOMENTS, DOVARY, N_PARAMETERS,                       & ! Input Numbers
            DO_LAYER_SCATTERING, BEAM_CUTOFF, FLUX_FACTOR, DFLUX, & ! Input flux/Bookkeeping
            L_OMEGA_GREEK, PI_XUP, PI_XUM, PI_X0P,                & ! Input optical/PI
            L_UPAR_DN_1, L_UPAR_UP_1 )                              ! Output

!  1/31/21. Version 2.8.3. New routine: Linearization Green's function user solution
!        -- No reference to LOCAL_UM_START.

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_STREAMS, &
                                 MAXBEAMS, MAX_ATMOSWFS, MAXMOMENTS, MAXSTREAMS_2, ZERO, PI4

      implicit none

!  subroutine input arguments
!  --------------------------

!  Direction flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Obervational Geometry flag

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Input numbers

      INTEGER, INTENT (IN) ::           LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS

!  Variation control

      LOGICAL, intent(in)  ::           DOVARY
      INTEGER, intent(in)  ::           N_PARAMETERS

!  Input bookkeeping 

      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      INTEGER, INTENT (IN) ::           BEAM_CUTOFF  ( MAXBEAMS )

!  Input Optical

      DOUBLE PRECISION, INTENT (IN) ::  L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Input PI matrices

      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Subroutine output arguments
!  ---------------------------

!  Linearized Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION, INTENT (INOUT) :: L_UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      LOGICAL   :: DO_FIRST
      INTEGER   :: UM, L, N, M, Q, IB, O1, O2
      DOUBLE PRECISION :: F1, SPS
      DOUBLE PRECISION :: HELP_Q1 ( 0:MAXMOMENTS, MAXSTOKES )
!  Layer and Fourier

      N = LAYER
      M = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

!  Check existence
!  ---------------

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. BEAM_CUTOFF(IB) ) .AND. DO_LAYER_SCATTERING(M,N)

!  If no solution or no variation, zero output and exit

      IF ( .NOT.DOVARY .OR. .NOT.DO_FIRST ) THEN
        DO Q = 1, N_PARAMETERS
          IF ( DO_UPWELLING ) L_UPAR_UP_1(1:N_USER_STREAMS, 1:NSTOKES,N,Q) = ZERO
          IF ( DO_DNWELLING ) L_UPAR_DN_1(1:N_USER_STREAMS, 1:NSTOKES,N,Q) = ZERO
        ENDDO
        RETURN
      ENDIF

!  start parameter loop

      DO Q = 1, N_PARAMETERS

!  For each moment do inner sum over computational angles

        DO O1 = 1, NSTOKES
          DO L = M, NMOMENTS
            SPS = ZERO
            DO O2 = 1, NSTOKES
              SPS = SPS + L_OMEGA_GREEK(L,N,O1,O2,Q) * sum( PI_X0P(L,IBEAM,N,O2,1:NSTOKES) * DFLUX(1:NSTOKES) )
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
                L_UPAR_UP_1(UM,O1,N,Q) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,UM,O1,1:NSTOKES) )
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              L_UPAR_UP_1(IBEAM,O1,N,Q) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,IBEAM,O1,1:NSTOKES) )
            ENDDO
          ENDIF
        ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct contribution Only. Observational Geometry, tie to the IBEAM index.

        IF ( DO_DNWELLING ) THEN
          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO O1 = 1, NSTOKES
              DO UM = 1, N_USER_STREAMS
                L_UPAR_DN_1(UM,O1,N,Q) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES) )
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              L_UPAR_DN_1(IBEAM,O1,N,Q) = sum( HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,IBEAM,O1,1:NSTOKES) )
            ENDDO
          ENDIF
        ENDIF

!  end parameter loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_GUSER_SOLUTION

!  End module

      END MODULE vlidort_lpc_solutions_m

