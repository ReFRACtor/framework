
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
! #              VLIDORT_LC_CBEAM_SOLUTION                      #
! #              VLIDORT_LC_CUSER_SOLUTION                      #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. No Changes except naming.
!    ==> Classical Modules renamed QBEAM --> CBEAM, UBEAM --> CUSER

      MODULE vlidort_lc_solutions_m

      PRIVATE
      PUBLIC :: VLIDORT_LC_CBEAM_SOLUTION, VLIDORT_LC_CUSER_SOLUTION

      CONTAINS

      SUBROUTINE VLIDORT_LC_CBEAM_SOLUTION ( &
        DO_PLANE_PARALLEL, DOVARY, GIVEN_LAYER, FOURIER, IBEAM,  & ! Input flags/numbers
        NSTOKES, NSTREAMS, NMOMENTS, N_PARAMETERS,               & ! Input numbers
        NSTREAMS_2, NSTKS_NSTRMS, DO_LAYER_SCATTERING,           & ! Input bookkeeping
        QUAD_STREAMS, FLUX_FACTOR, BEAM_CUTOFF, AVERAGE_SECANT,  & ! Input quads/Beam
        DFLUX, PI_XQP, PI_X0P, PI_XQM_POST, SAB, DAB,            & ! Input Eigenproblem
        QSUMVEC_SAVE, QDIFVEC_SAVE, QVEC_SAVE, QDIF_SAVE,        & ! Input solution
        QMAT_SAVE, QPIVOT, L_SAB, L_DAB, L_EIGENMAT,             & ! Input solution, lin Eigen.
        L_OMEGA_GREEK, LC_AVERAGE_SECANT,                        & ! Input optical + beam
        LC_BVEC, STATUS, MESSAGE, TRACE )                          ! Output + status

!  linearized values of the classical particular solution.

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS, MAXMOMENTS, MAXSTREAMS_2, &
                                 MAXSTRMSTKS, MAXEVALUES, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, HALF, TWO, PI4

      USE LAPACK_TOOLS_m, Only : DGETRS

      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DOVARY

!  Indices

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          N_PARAMETERS

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS

!  quadrature and bookkeeping

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Flux and beam parameterization

      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      INTEGER, INTENT (IN) ::          BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Pi matrix stuff

      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Eigenmatrix components

      DOUBLE PRECISION, INTENT (IN) :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )

!  Saved solution variables

      DOUBLE PRECISION, INTENT (IN) :: QSUMVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: QDIFVEC_SAVE ( MAXSTREAMS, MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) :: QVEC_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (IN) :: QDIF_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (IN) :: QMAT_SAVE ( MAXSTRMSTKS, MAXSTRMSTKS )
      INTEGER, INTENT (IN) ::          QPIVOT    ( MAXSTRMSTKS )

!  Linearized  values (strictly in)

      DOUBLE PRECISION, INTENT (IN) :: L_SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_EIGENMAT ( MAXEVALUES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized output (Ostensibly out)

      DOUBLE PRECISION, INTENT (INOUT) ::  LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Exception

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  Saved vector

      DOUBLE PRECISION ::  QSAVE ( MAXSTRMSTKS )

!  linearization arrays

      DOUBLE PRECISION ::  L_QVEC     ( MAXSTRMSTKS, MAX_ATMOSWFS )
      DOUBLE PRECISION ::  L_QSUMVEC  ( MAXSTREAMS,  MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION ::  L_QDIFVEC  ( MAXSTREAMS,  MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION ::  HELP_TEMPO ( 0:MAXMOMENTS, MAXSTOKES )

!  help variables

      DOUBLE PRECISION ::  HELP, WSUM, WDIF, SECBAR
      DOUBLE PRECISION ::  TM1, TM2, TM3, L_WVEC_D, L_WVEC_U
      DOUBLE PRECISION ::  SUM, GSUM, S_TPA, S_TMA, F1

      INTEGER          ::  I, J, I1, L, N, M, INFO, IB, Q
      INTEGER          ::  O1, O2, IR, JC, IROW, JCOL
      LOGICAL          ::  DO_FIRST
      CHARACTER (LEN=3) :: CN, CI

!  Start of code
!  -------------

!  initialise Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Set layer, fourier, beam index

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

!  safety

      SECBAR = ZERO

!  set up driving vector (using saved results)
!  This must be done regardless of whether layer N is varying or not.
!  Enhancement # 1, 6/27/16

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        QSAVE(IR+1:IR+NSTOKES) = QDIFVEC_SAVE(I,1:NSTOKES) + TWO*AVERAGE_SECANT(N,IB)*QVEC_SAVE(IR+1:IR+NSTOKES)
      ENDDO

!  Linearization for layer N is in two parts:
!    1A. Linearization due to variations in Layer N itself (Profiles)
!    1B. Linearization due to columns
!    2. Linearization due to variations in layers K < N (Profiles)

!  Part 1.
!  =======

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. BEAM_CUTOFF(IB) ) .AND. DO_LAYER_SCATTERING(M,N)

!  No particular solution beyond the cutoff layer.
!  Or no scattering in this layer
!    [ Zero the boundary layer values and Return (There is no Part 2 here) ]
!  Enhancement # 2, 6/27/16

      IF (.NOT. DOVARY .OR. .NOT. DO_FIRST  ) THEN
        LC_BVEC(1:NSTREAMS_2, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
        RETURN
      ENDIF

!  solar zenith cosine for this layer

      SECBAR = AVERAGE_SECANT(N,IB)

!  For each varying parameter

      DO Q = 1, N_PARAMETERS

!  Auxiliary matrix for Q functions
!  Enhancement # 3, 6/27/16

        DO O1 = 1, NSTOKES
          DO L = M, NMOMENTS
            GSUM = ZERO
            DO O2 = 1, NSTOKES
              GSUM = GSUM + L_OMEGA_GREEK(L,N,O1,O2,Q) * sum( PI_X0P(L,IBEAM,N,O2,1:NSTOKES) * DFLUX(1:NSTOKES) )
            ENDDO
            HELP_TEMPO(L,O1) = GSUM
          ENDDO
        ENDDO

!  Set up linearized sum and difference vectors for Beam source terms
!  Enhancement # 4, 6/27/16

        DO O1 = 1, NSTOKES
          DO I = 1, NSTREAMS
            S_TPA = sum( PI_XQP(M:NMOMENTS,I,O1,1:NSTOKES)      * HELP_TEMPO(M:NMOMENTS,1:NSTOKES) )
            S_TMA = sum( PI_XQM_POST(M:NMOMENTS,I,O1,1:NSTOKES) * HELP_TEMPO(M:NMOMENTS,1:NSTOKES) )
            L_QSUMVEC(I,O1,Q) = F1 * ( S_TPA + S_TMA ) / QUAD_STREAMS(I)
            L_QDIFVEC(I,O1,Q) = F1 * ( S_TPA - S_TMA ) / QUAD_STREAMS(I)
          ENDDO
        ENDDO

!   setup linearized RHS vector
!  ( use results from the original solution )
!    WARNING. BE CAREFUL OF SIGNS on TM1, TM2, TM3
!     slightly different from the scalar model case.
!  Enhancement # 5, 6/27/16

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            HELP = ZERO
            DO J = 1, NSTREAMS
              JC = NSTOKES*(J-1)
              HELP = HELP + sum( & 
                   - L_EIGENMAT(IROW,JC+1:JC+NSTOKES,N,Q) * QVEC_SAVE(JC+1:JC+NSTOKES) + &
                   DAB(I,J,O1,1:NSTOKES,N)   * L_QSUMVEC(J,1:NSTOKES,Q) + &
                   L_DAB(I,J,O1,1:NSTOKES,N,Q) *   QSUMVEC_SAVE(J,1:NSTOKES) )
            ENDDO
            L_QVEC(IROW,Q)  = HELP + L_QDIFVEC(I,O1,Q) * SECBAR

          ENDDO
        ENDDO

!  end parameter loop

      ENDDO

!  additional terms for the quasi-spherical case
!  ( layers greater than one )

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( N.GT.1 ) THEN
          DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                HELP = QSAVE(IROW) * LC_AVERAGE_SECANT(N,IB,Q)
                L_QVEC(IROW,Q) = L_QVEC(IROW,Q) + HELP

              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Solve problem by back substitution for all of the RHS vectors
!  ( uses L_U decomposition of QMAT_SAVE from original solution )

      CALL DGETRS &
              ('N',NSTKS_NSTRMS,N_PARAMETERS,QMAT_SAVE, &
                MAXSTRMSTKS,QPIVOT,L_QVEC,MAXSTRMSTKS,INFO)

!  Exception handling 1

      IF ( INFO .LT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO
        WRITE(CN, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRS call # 1 in VLIDORT_LC_QBEAM_SOLUTION,'// &
                  ' Beam P.I. linearization (N,N), layer # '//CN
        STATUS  = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  assign solutions for the quasi-spherical case, N > 1
!    Note change of sign with HELP, this is because SAB as
!    defined in VLIDORT is the negative of SAB in the scalar model.

      IF ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 ) THEN

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            TM3 = - QDIF_SAVE(IROW) / SECBAR
            DO Q = 1, N_PARAMETERS
              HELP = ZERO
              DO J = 1, NSTREAMS
                JC = NSTOKES*(J-1)
                DO O2 = 1, NSTOKES
                  JCOL = JC + O2
                  TM1 =   SAB(I,J,O1,O2,N)   * L_QVEC(JCOL,Q)
                  TM2 = L_SAB(I,J,O1,O2,N,Q) *   QVEC_SAVE(JCOL)
                  HELP = HELP + TM1 + TM2
                ENDDO
              ENDDO
              WSUM = L_QVEC(IROW,Q)
              TM2 = ( HELP - L_QSUMVEC(I,O1,Q) ) / SECBAR
              WDIF = LC_AVERAGE_SECANT(N,IB,Q) * TM3 + TM2
              L_WVEC_D = HALF * ( WSUM + WDIF )
              L_WVEC_U = HALF * ( WSUM - WDIF )
              LC_BVEC(I,O1,N,Q)  = L_WVEC_D
              LC_BVEC(I1,O1,N,Q) = L_WVEC_U
              IF ( O1 .GT. 2 ) LC_BVEC(I1,O1,N,Q) = - L_WVEC_U
            ENDDO
          ENDDO
        ENDDO

!  assign solutions for plane/parallel & quasi-spherical case N = 1

      ELSE

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_PARAMETERS
              HELP = ZERO
              DO J = 1, NSTREAMS
                JC = NSTOKES*(J-1)
                DO O2 = 1, NSTOKES
                  JCOL = JC + O2
                  TM1 =   SAB(I,J,O1,O2,N)   * L_QVEC(JCOL,Q)
                  TM2 = L_SAB(I,J,O1,O2,N,Q) *   QVEC_SAVE(JCOL)
                  HELP = HELP - TM1 - TM2
                ENDDO
              ENDDO
              WSUM = L_QVEC(IROW,Q)
              WDIF = ( - HELP - L_QSUMVEC(I,O1,Q) ) / SECBAR
              L_WVEC_D = HALF * ( WSUM + WDIF )
              L_WVEC_U = HALF * ( WSUM - WDIF )
              LC_BVEC(I,O1,N,Q)  = L_WVEC_D
              LC_BVEC(I1,O1,N,Q) = L_WVEC_U
              IF ( O1 .GT. 2 ) LC_BVEC(I1,O1,N,Q) = - L_WVEC_U
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Old code: Toggle to LC_BVEC
!      DO Q = 1, N_PARAMETERS
!       DO I = 1, NSTREAMS
!        I1 = I + NSTREAMS
!        DO O1 = 1, NSTOKES
!         LC_BVEC(I,O1,N,Q) = L_WVEC(I,O1,N,LVARY,Q)
!         HELP = ZERO
!         DO O2 = 1, NSTOKES
!          HELP = HELP + DMAT(O1,O2) * L_WVEC(I1,O2,N,LVARY,Q)
!         ENDDO
!         LC_BVEC(I1,O1,N,Q) = HELP
!        ENDDO
!       ENDDO
!      ENDDO

!  Part 2.
!  =======

!  No part 2 (Profile Jacobians only). 
!  Version 2.8, goto statement removed.

!  Continuation point
! 2222 CONTINUE

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LC_CBEAM_SOLUTION

!

      SUBROUTINE VLIDORT_LC_CUSER_SOLUTION ( &
        DO_OBSERVATION_GEOMETRY, DO_UPWELLING, DO_DNWELLING,         & ! input flags
        DO_VARY, GIVEN_LAYER, FOURIER, IBEAM,                        & ! Input flags/indices
        NSTOKES, NSTREAMS, NMOMENTS, N_USER_STREAMS, N_PARAMETERS,   & ! Input numbers
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DO_LAYER_SCATTERING, & ! Bookkeeping
        DFLUX, FLUX_FACTOR, QUAD_HALFWTS,                            & ! Bookkeeping and quads
        PI_XQP, PI_XUP, PI_XUM, PI_X0P, PI_XQM_PRE, BEAM_CUTOFF,     & ! Pi Matrices and cutoff
        OMEGA_GREEK, HELPSTOKES_BEAM, LC_BVEC, L_OMEGA_GREEK,        & ! solutions and optical
        L_UPAR_DN_1, L_UPAR_UP_1, LC_UPAR_DN_2, LC_UPAR_UP_2 )         ! output

!  1/31/21. Version 2.8.3. Drop LOCAL_UM_START (replaced by 1)

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAXMOMENTS, MAXSTREAMS_2, MAXSTRMSTKS, ZERO, PI4

      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_VARY

!  numbers

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_PARAMETERS

!  Bookkeeping
!    -- 1/31/21. Version 2.8.3. LOCAL_UM_START dropped

      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )

!  Flux and quadratures

      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )

!  Pi Matrices

      DOUBLE PRECISION, INTENT (IN) :: PI_XQP     ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUP     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_X0P     ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Beam stuff

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: HELPSTOKES_BEAM ( 0:MAXMOMENTS, MAXSTOKES )

!  Optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK   ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Linearized Beam vectors (Ostensibly in)

      DOUBLE PRECISION, INTENT (INOUT) :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  OUTPUT
!  ======

!  Linearized User solutions (Ostensibly OUTPUT)

      DOUBLE PRECISION, INTENT (INOUT) :: L_UPAR_DN_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UPAR_UP_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LC_UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LC_UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

!  help

      LOGICAL ::          DO_LAYER, DO_FIRST
      INTEGER ::          IB, UM, L, N, M, O1, O2, J, J1, Q
      DOUBLE PRECISION :: H, F1, SP, SM, SPS, SGW ( MAXSTOKES )

!  local arrays

      DOUBLE PRECISION :: L_HELP_Q1 ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: L_HELP_Q2 ( 0:MAXMOMENTS, MAXSTOKES )

      DOUBLE PRECISION :: L_HELPSTOKES_BEAM ( MAXSTOKES )

!  Data

      SGW = (/ 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)

!  Proxies

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM

!  initialize

      DO_LAYER = .FALSE.
      F1 = FLUX_FACTOR / PI4

!  Part 1. Linearizations for quantities varying in Layer N
!  ========================================================

!  Check existence
!  ---------------

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. BEAM_CUTOFF(IB) ) .AND. DO_LAYER_SCATTERING(M,N)

!  If no solution or no variation, zero output and RETURN
!     Enhancement # 6, 6/27/16
!     Enhancement # 7, 6/27/16
!     Versionrsion 2.8, remove GOTO statement and return instead (NO PART 2)
!    -- 1/31/21. Version 2.8.3. LOCAL_UM_START replaced by 1

      IF ( .NOT. DO_VARY .OR. .NOT. DO_FIRST ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          IF ( DO_UPWELLING ) THEN
            L_UPAR_UP_1 (1:N_USER_STREAMS, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
            LC_UPAR_UP_2(1:N_USER_STREAMS, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            L_UPAR_DN_1 (1:N_USER_STREAMS, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
            LC_UPAR_DN_2(1:N_USER_STREAMS, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
          ENDIF
        ELSE
          IF ( DO_UPWELLING ) THEN
            L_UPAR_UP_1 (IB, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
            LC_UPAR_UP_2(IB, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            L_UPAR_DN_1 (IB, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
            LC_UPAR_DN_2(IB, 1:NSTOKES, N, 1:N_PARAMETERS) = ZERO
          ENDIF
        ENDIF
        RETURN             !! GO TO 2222
      ENDIF

!  existence flag

      DO_LAYER = ( STERM_LAYERMASK_UP(N) .OR. STERM_LAYERMASK_DN(N) )

!  start parameter loop

      DO Q = 1, N_PARAMETERS

!  first function (linearization of primary scattering of beam)
!  Enhancement # 8, 6/27/16

        IF ( DO_LAYER ) THEN
          DO O1 = 1, NSTOKES
            DO L = M, NMOMENTS
              SPS = ZERO
              DO O2 = 1, NSTOKES
                 SPS = SPS + L_OMEGA_GREEK(L,N,O1,O2,Q) * sum( PI_X0P(L,IBEAM,N,O2,1:NSTOKES) * DFLUX(1:NSTOKES) )
              ENDDO
              L_HELP_Q1(L,O1) = SPS
            ENDDO
          ENDDO
        ENDIF

!  For each moment, do inner sum over computational angles
!   (Linearization of Diffuse scattering of beam)
!  Enhancement # 9, 6/27/16

        IF ( DO_LAYER ) THEN
          DO L = M, NMOMENTS
            DO O1 = 1, NSTOKES
              SPS = ZERO
              DO J = 1, NSTREAMS
                J1 = J + NSTREAMS
                H = QUAD_HALFWTS(J)
                SP = sum ( H  *                  PI_XQP    (L,J,O1,1:NSTOKES)*LC_BVEC(J ,1:NSTOKES,N,Q) )
                SM = sum ( H  * SGW(1:NSTOKES) * PI_XQM_PRE(L,J,O1,1:NSTOKES)*LC_BVEC(J1,1:NSTOKES,N,Q) )
                SPS = SPS + SP + SM
              ENDDO
              L_HELPSTOKES_BEAM(O1) = SPS
            ENDDO
            DO O1 = 1, NSTOKES
              L_HELP_Q2(L,O1) = sum( OMEGA_GREEK(L,N,O1,1:NSTOKES) * L_HELPSTOKES_BEAM(1:NSTOKES) + &
                                   L_OMEGA_GREEK(L,N,O1,1:NSTOKES,Q) * HELPSTOKES_BEAM(L,1:NSTOKES) )
            ENDDO
          ENDDO
        ENDIF

!  Now sum over all harmonic contributions (Upwelling)
!    Direct and integrated contributions
!  Enhancement # 10, 6/27/16
!  Enhancement # 11, 6/27/16
!    -- 1/31/21. Version 2.8.3. LOCAL_UM_START replaced by 1

        IF ( DO_UPWELLING .AND. STERM_LAYERMASK_UP(N) ) THEN
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO O1 = 1, NSTOKES
              DO UM = 1, N_USER_STREAMS
                L_UPAR_UP_1 (UM,O1,N,Q) = F1 * sum( L_HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,UM,O1,1:NSTOKES) )
                LC_UPAR_UP_2(UM,O1,N,Q) = sum( L_HELP_Q2(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,UM,O1,1:NSTOKES) )
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              L_UPAR_UP_1 (IB,O1,N,Q) = F1 * sum( L_HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,IB,O1,1:NSTOKES) )
              LC_UPAR_UP_2(IB,O1,N,Q) = sum( L_HELP_Q2(M:NMOMENTS,1:NSTOKES)*PI_XUM(M:NMOMENTS,IB,O1,1:NSTOKES) )
            ENDDO
          ENDIF
        ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct and integrated contributions
!     Enhancement # 12, 6/27/16
!     Enhancement # 13, 6/27/16
!    -- 1/31/21. Version 2.8.3. LOCAL_UM_START replaced by 1

        IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N) ) THEN
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO O1 = 1, NSTOKES
              DO UM = 1, N_USER_STREAMS
                L_UPAR_DN_1 (UM,O1,N,Q) = F1 * sum( L_HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES) )
                LC_UPAR_DN_2(UM,O1,N,Q) = sum( L_HELP_Q2(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,UM,O1,1:NSTOKES) )
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              L_UPAR_DN_1 (IB,O1,N,Q) = F1 * sum( L_HELP_Q1(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,IB,O1,1:NSTOKES) )
              LC_UPAR_DN_2(IB,O1,N,Q) = sum( L_HELP_Q2(M:NMOMENTS,1:NSTOKES)*PI_XUP(M:NMOMENTS,IB,O1,1:NSTOKES) )
            ENDDO
          ENDIF
        ENDIF

!  end parameter loop

      ENDDO

!  Part 2. Version 2.8, removed.
!  Continuation point
! 2222 CONTINUE

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LC_CUSER_SOLUTION

      END MODULE vlidort_lc_solutions_m

