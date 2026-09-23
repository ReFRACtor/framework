
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
! #            L_BVP_BACKSUB                                    #
! #            L_BVP_SURFACE_SETUP                              #
! #                                                             #
! ###############################################################

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO construction with IF block
!     * Use  VAR_INDEX and N_WEIGHTFUNCS instead

!  1/31/21. Version 2.8.3.
!   -- Bookkeeping changes (BRDF_F array defined locally, each Fourier)

      MODULE vlidort_lpc_bvproblem_m

      PRIVATE
      PUBLIC :: L_BVP_BACKSUB, L_BVP_SURFACE_SETUP

      CONTAINS

      SUBROUTINE L_BVP_BACKSUB ( &
        VAR_INDEX, N_WEIGHTFUNCS, NLAYERS, NTOTAL,   & ! Numbers
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS,          & ! Numbers
        NSTKS_NSTRMS_2, K_REAL, K_COMPLEX,           & ! Numbers
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Input BVP matrices
        COL2_WF, SCOL2_WF,                           & ! Input BVP vectors
        NCON, PCON, STATUS, MESSAGE, TRACE )           ! Output and excpetion handling

!  1/31/21. Version 2.8.3. No Changes here

!  Solves the linearized boundary value problem.

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_ATMOSWFS, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS
      USE LAPACK_TOOLS_m, Only : DGBTRS, DGETRS

      IMPLICIT NONE

!  Linearization control

      INTEGER, INTENT (IN) ::          VAR_INDEX
      INTEGER, INTENT (IN) ::          N_WEIGHTFUNCS

!  numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  Bookkeeping

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  BVP matrices

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2    ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT  ( MAXSTRMSTKS_2 )

!  Integration constants. Not needed  Version 2.8
!      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
!      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Column Vectors

      DOUBLE PRECISION, INTENT (INOUT) :: COL2_WF  ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  Linearizaed integration constants

      DOUBLE PRECISION, INTENT (OUT) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local variables
!  ---------------

      INTEGER ::           C0, LAY, K, Q, KO1, K0, K1, K2
      INTEGER ::           IROW, IROW1, IROW_S, IROW1_S
      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI, CN

!  Intialize status

      STATUS = VLIDORT_SUCCESS

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

!  Version 2.8. 6/28/16. revised exception handling

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_WEIGHTFUNCS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              COL2_WF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          if ( VAR_INDEX .eq. 0 ) then
            TRACE   = 'Column_Wfs ; DGBTRS call in L_BVP_BACKSUB'
          else
            WRITE(CN, '(I3)' ) VAR_INDEX
            TRACE   = 'Profile_Wfs with varying layer = '//CN// '; DGBTRS call in L_BVP_BACKSUB'
          endif
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON and PCON, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTKS_NSTRMS_2
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            DO Q = 1, N_WEIGHTFUNCS
              NCON(K,LAY,Q) = COL2_WF(C0+IROW,Q)
              PCON(K,LAY,Q) = COL2_WF(C0+IROW1,Q)
            ENDDO
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(LAY)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(LAY)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            DO Q = 1, N_WEIGHTFUNCS
              NCON(K1,LAY,Q) = COL2_WF(C0+IROW,   Q)
              NCON(K2,LAY,Q) = COL2_WF(C0+IROW_S, Q)
              PCON(K1,LAY,Q) = COL2_WF(C0+IROW1,  Q)
              PCON(K2,LAY,Q) = COL2_WF(C0+IROW1_S,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS &
           ( 'N', NTOTAL, N_WEIGHTFUNCS, SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              SCOL2_WF, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE = 'Atmos_Wfs for 1-layer: DGETRS call in L_BVP_BACKSUB'
          STATUS = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set linearized integration constants NCON and PCON, 1 layer

        LAY = 1
        KO1 = K_REAL(LAY) + 1
        DO K = 1, K_REAL(LAY)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          DO Q = 1, N_WEIGHTFUNCS
            NCON(K,LAY,Q) = SCOL2_WF(IROW, Q)
            PCON(K,LAY,Q) = SCOL2_WF(IROW1,Q)
          ENDDO
        ENDDO
        DO K = 1, K_COMPLEX(LAY)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(LAY)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(LAY) + K_COMPLEX(LAY)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          DO Q = 1, N_WEIGHTFUNCS
            NCON(K1,LAY,Q) = SCOL2_WF(IROW,    Q)
            NCON(K2,LAY,Q) = SCOL2_WF(IROW_S,  Q)
            PCON(K1,LAY,Q) = SCOL2_WF(IROW1,   Q)
            PCON(K2,LAY,Q) = SCOL2_WF(IROW1_S, Q)
          ENDDO
        ENDDO

      ENDIF

!  debug------------------------------------------

!       if ( fourier .EQ.3 .and.ibeam.eq.1 ) then
!        do n = 1, nlayers
!          write(76,'(a,2i3,1p2e18.10)')
!     &     'hey',n,variation_index,ncon(3,n,1),pcon(3,n,1)
!        enddo
!       endif

!        if ( do_debug_write ) then
!        if ( fourier.eq.1 ) then
!         Q = 1
!         DO LAY = 1, NLAYERS
!          KO1 = K_REAL(LAY) + 1
!          DO K = 1, K_REAL(LAY)
!           write(51,'(4i3,1p4e20.10)')FOURIER,IBEAM,K,LAY,
!     &                LCON(K,LAY),  MCON(K,LAY),
!     &                NCON(K,LAY,Q),PCON(K,LAY,Q)
!          ENDDO
!          DO K = 1, K_COMPLEX(LAY)
!           K0 = 2*K - 2
!           K1 = KO1 + K0
!           K2 = K1  + 1
!           write(51,'(4i3,1p8e20.10)')FOURIER,IBEAM,K,LAY,
!     &                LCON(K1,LAY), MCON(K1,LAY),
!     &                LCON(K2,LAY), MCON(K2,LAY),
!     &                NCON(K1,LAY,Q),PCON(K1,LAY,Q),
!     &                NCON(K2,LAY,Q),PCON(K2,LAY,Q)
!          ENDDO
!         ENDDO
!        endif
!        ENDIF

!  finish

      RETURN
      END SUBROUTINE L_BVP_BACKSUB

!

      SUBROUTINE L_BVP_SURFACE_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BOUNDARY, & ! Flags
        IBEAM, FOURIER, N_WEIGHTFUNCS, NSTOKES, NSTREAMS, NLAYERS,    & ! Numbers
        SURFACE_FACTOR, QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,      & ! Surface input
        NSTKS_NSTRMS, MUELLER_INDEX, K_REAL, K_COMPLEX,               & ! bookkeeping
        SOLA_XPOS, T_DELT_EIGEN, L_T_DELT_EIGEN,                      & ! RT Solutions
        L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,                           & ! Linearized RT Solutions
        R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )                               ! Output reflected solutions

!  1/31/21. Version 2.8.3. Bookkeeping changes (BRDF/SLEAVE arrays defined locally, each Fourier)

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTREAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTOKES_SQ, ZERO

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          MODIFIED_BOUNDARY

!  Numbers

      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          N_WEIGHTFUNCS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Surface inputs
!    -- 1/31/21. Version 2.8.3. BRDF array defined locally, drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  Bookkeeping

      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  Solutions

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

!  Linearized solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output - linearized surface-reflected solutions

      DOUBLE PRECISION, INTENT (OUT) :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  Local variables
!  ---------------

!  help arrays

      DOUBLE PRECISION :: PV_W ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  Other variables

      INTEGER ::          I, J, O1, IB, Q, N, M
      INTEGER ::          K, KO1, K0, K1, K2, MO2(MAXSTOKES)

      DOUBLE PRECISION :: H1R, H1I, REFL_P, REFL_B, REFL_M, KMULT
      DOUBLE PRECISION :: H1, H2, HP, HM

      DOUBLE PRECISION :: H1_CR,   H2_CR,   H1_CI,   H2_CI
      DOUBLE PRECISION :: H1_S_CR, H2_S_CR, H1_S_CI, H2_S_CI

!  Initial section
!  ---------------

!  Beam index and Fourier component

      IB = IBEAM
      M = FOURIER

!  Always zero the result to start

      DO I = 1, NSTREAMS
       DO O1 = 1, NSTOKES
        DO Q = 1, N_WEIGHTFUNCS
          R2_L_BEAM(I,O1,Q)    = ZERO
          DO K = 1, NSTKS_NSTRMS
            R2_L_HOMP(I,O1,K,Q) = ZERO
            R2_L_HOMM(I,O1,K,Q) = ZERO
          ENDDO
        ENDDO
       ENDDO
      ENDDO

!  Return if no surface contributions

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if Fourier component > 0 (Lambertian)

      IF ( DO_LAMBERTIAN_SURFACE .and. M.gt.0 ) RETURN

!  Set up Auxiliary arrays
!  -----------------------

!  Last layer

      N = NLAYERS

!  Need Beam arrays regardless
!   Always something from layers above the PBL and also from PBL

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          DO Q = 1, N_WEIGHTFUNCS
            PV_W(J,O1,Q) = L_WLOWER(J,O1,N,Q) * QUAD_STRMWTS(J)
          ENDDO
        ENDDO
      ENDDO

!   Modified boundary condition
!     This applies when PBL is varying layer.
!       Require homogeneous solutions linearizations

      IF ( MODIFIED_BOUNDARY ) THEN

!  start loops

        DO J = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, N_WEIGHTFUNCS

!  real homogeneous solution contributions

           DO K = 1, K_REAL(N)
            H1 = L_SOLA_XPOS(J,O1,K,N,Q) *   T_DELT_EIGEN(K,N) + &
                   SOLA_XPOS(J,O1,K,N)   * L_T_DELT_EIGEN(K,N,Q)
            H2 = L_SOLB_XNEG(J,O1,K,N,Q)
            HV_P(J,O1,K,Q) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K,Q) = QUAD_STRMWTS(J)*H2
           ENDDO

!  Complex homogeneous solution contributions

           KO1 = K_REAL(N) + 1
           DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            H1R = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K1,N)   &
                - L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K2,N)   &
                +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K1,N,Q) &
                -   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K2,N,Q)
            H1I = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K2,N)   &
                + L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K1,N)   &
                +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K2,N,Q) &
                +   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K1,N,Q)
            HV_P(J,O1,K1,Q) = QUAD_STRMWTS(J)* H1R
            HV_P(J,O1,K2,Q) = QUAD_STRMWTS(J)* H1I
            HV_M(J,O1,K1,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K1,N,Q)
            HV_M(J,O1,K2,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K2,N,Q)
           ENDDO

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  End modified BCL4 condition

      ENDIF

!  Lambertian condition
!  ====================

      if ( DO_LAMBERTIAN_SURFACE ) THEN

!  reflection

        KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO

!  only 1 component

        O1 = 1

!  Integrated Downward reflection (Calculation)
!  --------------------------------------------

        DO Q = 1, N_WEIGHTFUNCS

!  Particular solution (only for the first Stokes component)

          REFL_B = KMULT * SUM(PV_W(1:NSTREAMS,O1,Q))
          R2_L_BEAM(1:NSTREAMS,O1,Q) = REFL_B

!  Homogeneous solutions for the modified condition

          IF ( MODIFIED_BOUNDARY ) THEN

!  Homogeneous real solutions

            DO K = 1, K_REAL(NLAYERS)
              REFL_P = KMULT * SUM(HV_P(1:NSTREAMS,O1,K,Q))
              REFL_M = KMULT * SUM(HV_M(1:NSTREAMS,O1,K,Q))
              R2_L_HOMP(1:NSTREAMS,O1,K,Q) = REFL_P
              R2_L_HOMM(1:NSTREAMS,O1,K,Q) = REFL_M
            ENDDO

!  Homogeneous complex solutions

            KO1 = K_REAL(NLAYERS) + 1
            DO K = 1, K_COMPLEX(NLAYERS)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              H1_CR = KMULT * SUM(HV_P(1:NSTREAMS,O1,K1,Q))
              H1_CI = KMULT * SUM(HV_P(1:NSTREAMS,O1,K2,Q))
              H2_CR = KMULT * SUM(HV_M(1:NSTREAMS,O1,K1,Q))
              H2_CI = KMULT * SUM(HV_M(1:NSTREAMS,O1,K2,Q))
              R2_L_HOMP(1:NSTREAMS,O1,K1,Q) = H1_CR
              R2_L_HOMP(1:NSTREAMS,O1,K2,Q) = H1_CI
              R2_L_HOMM(1:NSTREAMS,O1,K1,Q) = H2_CR
              R2_L_HOMM(1:NSTREAMS,O1,K2,Q) = H2_CI
            ENDDO

!  End modified boundary condition clause

          ENDIF

!  end parameter loop

        ENDDO

!  BRDF surface condition
!  ======================

      ELSE

!  start loops

        DO Q = 1, N_WEIGHTFUNCS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES

!  Particular solution
!    -- 1/31/21. Version 2.8.3. BRDF array defined locally, drop FOURIER index "M"

              MO2(1:NSTOKES) = MUELLER_INDEX(O1,1:NSTOKES)
              REFL_B = ZERO
              DO J = 1, NSTREAMS
                H1 = DOT_PRODUCT(PV_W(J,1:NSTOKES,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                REFL_B = REFL_B + H1
              ENDDO
              REFL_B = REFL_B * SURFACE_FACTOR
              R2_L_BEAM(I,O1,Q) = R2_L_BEAM(I,O1,Q) + REFL_B

!  Homogeneous solutions for the modified condition

              IF ( MODIFIED_BOUNDARY ) THEN

!  Homogeneous real solutions
!    -- 1/31/21. Version 2.8.3. BRDF array defined locally, drop FOURIER index "M"

                DO K = 1, K_REAL(NLAYERS)
                  REFL_P = ZERO
                  REFL_M = ZERO
                  DO J = 1, NSTREAMS
                    HP = DOT_PRODUCT(HV_P(J,1:NSTOKES,K,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    HM = DOT_PRODUCT(HV_M(J,1:NSTOKES,K,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    REFL_P = REFL_P + HP
                    REFL_M = REFL_M + HM
                  ENDDO
                  R2_L_HOMP(I,O1,K,Q) =  REFL_P * SURFACE_FACTOR
                  R2_L_HOMM(I,O1,K,Q) =  REFL_M * SURFACE_FACTOR
                ENDDO

!  homogeneous complex solutions
!    -- 1/31/21. Version 2.8.3. BRDF array defined locally, drop FOURIER index "M"

                KO1 = K_REAL(NLAYERS) + 1
                DO K = 1, K_COMPLEX(NLAYERS)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  H1_CR = ZERO
                  H2_CR = ZERO
                  H1_CI = ZERO
                  H2_CI = ZERO
                  DO J = 1, NSTREAMS
                    H1_S_CR = DOT_PRODUCT(HV_P(J,1:NSTOKES,K1,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    H2_S_CR = DOT_PRODUCT(HV_M(J,1:NSTOKES,K1,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    H1_S_CI = DOT_PRODUCT(HV_P(J,1:NSTOKES,K2,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    H2_S_CI = DOT_PRODUCT(HV_M(J,1:NSTOKES,K2,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    H1_CR = H1_CR + H1_S_CR
                    H2_CR = H2_CR + H2_S_CR
                    H1_CI = H1_CI + H1_S_CI
                    H2_CI = H2_CI + H2_S_CI
                  ENDDO
                  R2_L_HOMP(I,O1,K1,Q) = H1_CR * SURFACE_FACTOR
                  R2_L_HOMM(I,O1,K1,Q) = H2_CR * SURFACE_FACTOR
                  R2_L_HOMP(I,O1,K2,Q) = H1_CI * SURFACE_FACTOR
                  R2_L_HOMM(I,O1,K2,Q) = H2_CI * SURFACE_FACTOR
                ENDDO

!  End modified condition

              ENDIF

!  end parameter and stream and Stokes loops

            ENDDO
          ENDDO
        ENDDO

!  End clause BRDF vs. Lambertian

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_BVP_SURFACE_SETUP

      END MODULE vlidort_lpc_bvproblem_m
