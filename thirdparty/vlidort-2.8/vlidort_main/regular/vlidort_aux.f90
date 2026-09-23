
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
! #  This is VLIDORT_AUX.f90. This is a collection of utility   #
! #  modules compiled from a variety of sources. This is the    #
! #  only part of the VLIDORT package not written by the Author.#
! #  The subroutines in VLIDORT_AUX are listed below with their #
! #  source of origin (order of appearance).                    #
! #                                                             #
! #      lcstring:              E. Mikusch, DLR 1994            #
! #      len_string:            E. Mikusch, DLR 1994            #
! #      vlidort_read_error:    M. Christi, 2017                #
! #      vlidort_write_status:  R.J.D. Spurr, RTS, 2005, 2010   #
! #      asymtx:                adapted by RJDS from the        #
! #                             DISORT module ASYMTX, 1988      #
! #      getquad2               M. Christi, 2017                #
! #      issort                 M. Christi, 2017                #
! #      rssort                 M. Christi, 2017                #
! #      rqsort_idx             M. Christi, 2017                #
! #      cfplgarr:              T.P. Kurosu, 1997 (adapted      #
! #                               Numerical Recipes).           #
! #      xfindpar:              J. Lavagnino, SAO, 1991         #
! #      gfindpar:              J. Lavagnino, SAO, 1991         #
! #      FINDPAR_ERROR          R. Spurr, 2001                  #
! #                                                             #
! ###############################################################

!  Version 2.8. LEN_STRING and LCSTRING have been removed from everywhere except XFINDPAR. 7/5/16
!   --> use Trim and Adjustl functions in Fortran 90 wherever possible.

!  1/31/21. Version 2.8.3. Changed the call to ASYMTX to include tolerance variable

      MODULE vlidort_aux_m

!  Parameter types

      USE VLIDORT_pars_m

      CONTAINS

      SUBROUTINE LCSTRING ( STRING )

!  Convert a string to lower-case.  (ASCII character set assumed)

      IMPLICIT NONE

!  Modified argument

      CHARACTER * (*)     STRING

!  Local variables

      INTEGER         CCC, III, OFFSET

      OFFSET = ICHAR ( 'a' ) - ICHAR ( 'A' )
      DO 10 III = 1, LEN ( STRING )
        CCC = ICHAR ( STRING ( III : III ) )
        IF ( CCC .GE. ICHAR ( 'A' ) .AND. CCC .LE. ICHAR ( 'Z' ) ) THEN
           STRING ( III : III ) = CHAR ( CCC + OFFSET )
        END IF
10    CONTINUE

      END SUBROUTINE LCSTRING

!

      INTEGER FUNCTION LEN_STRING(STR)

!       ----------------------------------------------------------------
!       Length of string.
!
!       This function determines the length of the string contained in
!       the character variable. Therefore it searches the last non-blank
!       character and returnes the position of this value as the length.
!       If a character variable contains blanks only, 0 is returned.
!        2 Feb 1994, Eberhard Mikusch
!       ----------------------------------------------------------------

      IMPLICIT NONE

!     Input argument

      CHARACTER * (*) STR

!     Get last non-blank character

      LEN_STRING = LEN(STR)
      DO WHILE ((LEN_STRING .GT.0) .AND. (STR(LEN_STRING:LEN_STRING) .EQ. ' '))
        LEN_STRING = LEN_STRING - 1
      END DO
      RETURN

      END FUNCTION LEN_STRING

!

      SUBROUTINE VLIDORT_READ_ERROR ( ERRORFILE, VLIDORT_InputStatus )

!  Module, dimensions and numbers

      USE VLIDORT_pars_m, only : VLIDORT_ERRUNIT
      USE VLIDORT_Outputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  VLIDORT file inputs status

      TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = VLIDORT_ERRUNIT

!  Write VLIDORT configuration file read errors to VLIDORT error file 

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*)' FATAL:   Wrong input from VLIDORT input file-read'
      WRITE(W,*)'  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)')'    ** Number of messages = ',VLIDORT_InputStatus%TS_NINPUTMESSAGES
      DO N = 1, VLIDORT_InputStatus%TS_NINPUTMESSAGES
        WRITE(W,'(A,I3,A,A)')'Message # ',N,' : ',&
          ADJUSTL(TRIM(VLIDORT_InputStatus%TS_INPUTMESSAGES(N)))
        WRITE(W,'(A,I3,A,A)')'Action  # ',N,' : ',&
          ADJUSTL(TRIM(VLIDORT_InputStatus%TS_INPUTACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Read-input fail: Look at file ' // TRIM(ERRORFILE)
      STOP 'Program stopped in call to VLIDORT_READ_ERROR'

      END SUBROUTINE VLIDORT_READ_ERROR

!

      SUBROUTINE VLIDORT_WRITE_STATUS &
         ( ERRORFILE, ERRORUNIT, OPENFILEFLAG, VLIDORT_Status )

!  Upgraded, 13 October 2010

!mick mod 9/19/2017 - added variable VLIDORT_ABORT
!                   - enhanced internal error handling for fatal cases

      USE VLIDORT_PARS_m, Only : MAX_MESSAGES,    VLIDORT_SUCCESS, &
                                 VLIDORT_SERIOUS, VLIDORT_WARNING
      USE VLIDORT_Outputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file and unit (will only opened if overall status not successfu

      CHARACTER (LEN=*), INTENT(IN) :: ERRORFILE
      INTEGER, INTENT(IN)           :: ERRORUNIT

      LOGICAL, INTENT(INOUT)        :: OPENFILEFLAG

!  VLIDORT status

      TYPE(VLIDORT_Exception_Handling), INTENT(IN) :: VLIDORT_Status

!  Local variables

!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NCHECKMESSAGES
      CHARACTER (LEN=120) :: CKMESSAGES (0:MAX_MESSAGES)
      CHARACTER (LEN=120) :: ACTIONS (0:MAX_MESSAGES)

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER ::             STATUS_CALCULATION
      CHARACTER (LEN=120) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  Help

      INTEGER ::           N, W
      LOGICAL ::           STATUS_OVERALL, VLIDORT_ABORT

!  Copy inputs to local variables

      STATUS_INPUTCHECK            = VLIDORT_Status%TS_STATUS_INPUTCHECK
      NCHECKMESSAGES               = VLIDORT_Status%TS_NCHECKMESSAGES
      CKMESSAGES(0:NCHECKMESSAGES) = VLIDORT_Status%TS_CHECKMESSAGES(0:NCHECKMESSAGES)
      ACTIONS(0:NCHECKMESSAGES)    = VLIDORT_Status%TS_ACTIONS(0:NCHECKMESSAGES)

!mick fix 2/4/2013  - added IF condition
!mick mod 9/19/2017 - added ELSE section here; removed previous IF block
!                     that followed this one

      IF (STATUS_INPUTCHECK .NE. VLIDORT_SERIOUS) THEN
        STATUS_CALCULATION = VLIDORT_Status%TS_STATUS_CALCULATION

        MESSAGE            = VLIDORT_Status%TS_MESSAGE
        TRACE_1            = VLIDORT_Status%TS_TRACE_1
        TRACE_2            = VLIDORT_Status%TS_TRACE_2
        TRACE_3            = VLIDORT_Status%TS_TRACE_3

        STATUS_OVERALL = ( STATUS_INPUTCHECK  .NE. VLIDORT_SUCCESS ) .or. &
                         ( STATUS_CALCULATION .NE. VLIDORT_SUCCESS )
      ELSE
        STATUS_OVERALL = .true.
      ENDIF

!  Open if flagged

      IF ( .not.OPENFILEFLAG .and. STATUS_OVERALL ) THEN
        OPENFILEFLAG = .true.
        OPEN(ERRORUNIT, FILE = ERRORFILE, STATUS = 'unknown' )
      ENDIF

!  Define some local variables

      W = ERRORUNIT
      VLIDORT_ABORT = .false.

!  Check status and write messages

      IF ( STATUS_INPUTCHECK .EQ. VLIDORT_WARNING ) THEN
        WRITE(W,*)' WARNING: Suspicious input from VLIDORT input check'
        WRITE(W,*)'  ------ VLIDORT has executed with internal defaults'
        WRITE(W,*)'  ------ Here are the messages and actions '
        write(W,'(A,I3)')'    ** Number of messages = ',NCHECKMESSAGES
        DO N = 1, NCHECKMESSAGES
          write(W,'(A,I3,A,A)')'Message # ',N,' : ',Adjustl(Trim(CKMESSAGES(N)))
          write(W,'(A,I3,A,A)')'Action  # ',N,' : ',Adjustl(Trim(ACTIONS(N)))
        ENDDO
      ELSE IF ( STATUS_INPUTCHECK .EQ. VLIDORT_SERIOUS ) THEN
        WRITE(W,*)' FATAL:   Wrong input from VLIDORT input check'
        WRITE(W,*)'  ------ VLIDORT will not execute'
        WRITE(W,*)'  ------ Here are the messages and actions '
        write(W,'(A,I3)')'    ** Number of messages = ',NCHECKMESSAGES
        DO N = 1, NCHECKMESSAGES
          write(W,'(A,I3,A,A)')'Message # ',N,' : ',Adjustl(Trim(CKMESSAGES(N)))
          write(W,'(A,I3,A,A)')'Action  # ',N,' : ',Adjustl(Trim(ACTIONS(N)))
        ENDDO

        !Set abort flag and alert the user at the screen
        VLIDORT_ABORT = .true.
        WRITE(*,*)
        WRITE(*,*)'VLIDORT input abort'
      ENDIF

!  Model calculation status

      IF ( STATUS_INPUTCHECK .NE. VLIDORT_SERIOUS ) THEN
        IF ( STATUS_CALCULATION .NE. VLIDORT_SUCCESS ) THEN
          WRITE(W,*)' FATAL: VLIDORT execution failure '
          WRITE(W,*)'  -----Here is the message and traces : '
          WRITE(W,'(A,A)')' Message : ', Adjustl(Trim(MESSAGE))
          WRITE(W,'(A,A)')' Trace 1 : ', Adjustl(Trim(TRACE_1))
          WRITE(W,'(A,A)')' Trace 2 : ', Adjustl(Trim(TRACE_2))
          WRITE(W,'(A,A)')' Trace 3 : ', Adjustl(Trim(TRACE_3))

          !Set abort flag and alert the user at the screen
          VLIDORT_ABORT = .true.
          WRITE(*,*)
          WRITE(*,*)'VLIDORT calculation abort'
        ENDIF
      ENDIF

!  If fatal error, abort

      IF ( VLIDORT_ABORT ) THEN
        CLOSE(ERRORUNIT)
        WRITE(*,*)
        WRITE(*,'(/1x,a)')'Main program failed to execute, ' // &
          'error messages in ' // trim(ERRORFILE)
        STOP 'Program stopped in call to VLIDORT_WRITE_STATUS'
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_WRITE_STATUS

!

      SUBROUTINE  ASYMTX( TOL, AAD, M, IA, IEVEC, &
                            EVECD, EVALD, IER, WKD, &
                            MESSAGE, BAD_STATUS )

!  1/31/21. Version 2.8.3. Add the tolerance variable

!    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

!       Solves eigenfunction problem for real asymmetric matrix
!       for which it is known a priori that the eigenvalues are real.

!       This is an adaptation of a subroutine EIGRF in the IMSL
!       library to use real instead of complex arithmetic, accounting
!       for the known fact that the eigenvalues and eigenvectors in
!       the discrete ordinate solution are real.  Other changes include
!       putting all the called subroutines in-line, deleting the
!       performance index calculation, updating many DO-loops
!       to Fortran77, and in calculating the machine precision
!       TOL instead of specifying it in a data statement.

!       EIGRF is based primarily on EISPACK routines.  The matrix is
!       first balanced using the parlett-reinsch algorithm.  Then
!       the Martin-Wilkinson algorithm is applied.

!       References:
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!             Sources and Development of Mathematical Software,
!             Prentice-Hall, Englewood Cliffs, NJ
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!             Clarendon Press, Oxford

!   I N P U T    V A R I A B L E S:

!        AAD  :  input asymmetric matrix, destroyed after solved
!        M    :  order of  A
!       IA    :  first dimension of  A
!    IEVEC    :  first dimension of  EVECD

!   O U T P U T    V A R I A B L E S:

!       EVECD :  (unnormalized) eigenvectors of  A
!                   ( column J corresponds to EVALD(J) )

!       EVALD :  (unordered) eigenvalues of  A ( dimension at least M )

!       IER   :  if .NE. 0, signals that EVALD(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.

!  include file for dimensioning

      USE VLIDORT_PARS_m, Only : MAXSTRMSTKS, ZERO, ONE

      implicit none

      LOGICAL            BAD_STATUS
      CHARACTER*(*)      MESSAGE

!   S C R A T C H   V A R I A B L E S:

!       WKD    :  WORK AREA ( DIMENSION AT LEAST 2*M )
!+---------------------------------------------------------------------+

!  input/output arguments
!    - 1/31/21. Version 2.8.3. Add the tolerance variable

      INTEGER              M, IA, IEVEC, IER
      DOUBLE PRECISION TOL, &
               AAD(MAXSTRMSTKS,MAXSTRMSTKS), WKD(4*MAXSTRMSTKS), &
               EVALD(MAXSTRMSTKS), EVECD(MAXSTRMSTKS,MAXSTRMSTKS)

!  local variables (explicit declaration

      LOGICAL           NOCONV, NOTLAS
      INTEGER              I, J, L, K, KKK, LLL
      INTEGER              N, N1, N2, IN, LB, KA, II
      DOUBLE PRECISION  C1, C2, C3, C4, C5, C6
      DATA             C1 / 0.4375D0 /
      DATA               C2 / 0.5D0 /
      DATA               C3 / 0.75D0 /
      DATA              C4 / 0.95D0 /
      DATA              C5 / 16.0D0 /
      DATA              C6 / 256.0D0 /
      DOUBLE PRECISION  DISCRI, SGN, RNORM, W, F, G, H, P, Q, R
      DOUBLE PRECISION  REPL, COL, ROW, SCALE, T, X, Z, S, Y, UU, VV

      IER = 0
      BAD_STATUS = .FALSE.
      MESSAGE = ' '

!  1/31/21. Version 2.8.3. Removed this section
!---4/18/08. Use 1.0d-12 (1.0d-24 has convergence issues in Rayleigh lay
!      TOL = 0.0000001
!      TOL = 1.0D-12
!      TOL = 1.0D-24
!      TOL = 1.0D-20

!       Here change to bypass D1MACH:
!        TOL = D1MACH(4)

      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M ) THEN
        MESSAGE = 'ASYMTX--bad input variable(s)'
        BAD_STATUS = .TRUE.
        RETURN
      ENDIF
!                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES
      IF ( M.EQ.1 )  THEN
         EVALD(1) = AAD(1,1)
         EVECD(1,1) = 1.0D0
         RETURN
      ELSE IF ( M.EQ.2 )  THEN
         DISCRI = ( AAD(1,1) - AAD(2,2) )**2 + 4.0D0*AAD(1,2)*AAD(2,1)
         IF ( DISCRI.LT.ZERO ) THEN
           MESSAGE = 'ASYMTX--COMPLEX EVALS IN 2X2 CASE'
           BAD_STATUS = .TRUE.
           RETURN
         ENDIF
         SGN = ONE
         IF ( AAD(1,1).LT.AAD(2,2) )  SGN = - ONE
         EVALD(1) = 0.5D0*( AAD(1,1) + AAD(2,2) + SGN*DSQRT(DISCRI) )
         EVALD(2) = 0.5D0*( AAD(1,1) + AAD(2,2) - SGN*DSQRT(DISCRI) )
         EVECD(1,1) = ONE
         EVECD(2,2) = ONE
         IF ( AAD(1,1).EQ.AAD(2,2) .AND. &
               (AAD(2,1).EQ.ZERO.OR.AAD(1,2).EQ.ZERO) ) THEN
            RNORM = DABS(AAD(1,1))+DABS(AAD(1,2))+ &
                      DABS(AAD(2,1))+DABS(AAD(2,2))
            W = TOL * RNORM
            EVECD(2,1) = AAD(2,1) / W
            EVECD(1,2) = - AAD(1,2) / W
         ELSE
            EVECD(2,1) = AAD(2,1) / ( EVALD(1) - AAD(2,2) )
            EVECD(1,2) = AAD(1,2) / ( EVALD(2) - AAD(1,1) )
         ENDIF
         RETURN
      END IF
!                                        ** INITIALIZE OUTPUT VARIABLES
      DO 20 I = 1, M
         EVALD(I) = ZERO
         DO 10 J = 1, M
            EVECD(I,J) = ZERO
10       CONTINUE
         EVECD(I,I) = ONE
20    CONTINUE
!                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
!                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
!                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                  ** AND PUSH THEM DOWN
      RNORM = ZERO
      L  = 1
      K  = M

30    KKK = K
         DO 70  J = KKK, 1, -1
            ROW = ZERO
            DO 40 I = 1, K
               IF ( I.NE.J ) ROW = ROW + DABS( AAD(J,I) )
40          CONTINUE
            IF ( ROW.EQ.ZERO ) THEN
               WKD(K) = J
               IF ( J.NE.K ) THEN
                  DO 50 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,K)
                     AAD(I,K) = REPL
50                CONTINUE
                  DO 60 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(K,I)
                     AAD(K,I) = REPL
60                CONTINUE
               END IF
               K = K - 1
               GO TO 30
            END IF
70       CONTINUE
!                                     ** SEARCH FOR COLUMNS ISOLATING AN
!                                       ** EIGENVALUE AND PUSH THEM LEFT
80    LLL = L
         DO 120 J = LLL, K
            COL = ZERO
            DO 90 I = L, K
               IF ( I.NE.J ) COL = COL + DABS( AAD(I,J) )
90          CONTINUE
            IF ( COL.EQ.ZERO ) THEN
               WKD(L) = J
               IF ( J.NE.L ) THEN
                  DO 100 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,L)
                     AAD(I,L) = REPL
100               CONTINUE
                  DO 110 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(L,I)
                     AAD(L,I) = REPL
110               CONTINUE
               END IF
               L = L + 1
               GO TO 80
            END IF
120      CONTINUE
!                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K
      DO 130 I = L, K
         WKD(I) = ONE
130   CONTINUE

140   NOCONV = .FALSE.
         DO 200 I = L, K
            COL = ZERO
            ROW = ZERO
            DO 150 J = L, K
               IF ( J.NE.I ) THEN
                  COL = COL + DABS( AAD(J,I) )
                  ROW = ROW + DABS( AAD(I,J) )
               END IF
150         CONTINUE
            F = ONE
            G = ROW / C5
            H = COL + ROW
160         IF ( COL.LT.G ) THEN
               F   = F * C5
               COL = COL * C6
               GO TO 160
            END IF
            G = ROW * C5
170         IF ( COL.GE.G ) THEN
               F   = F / C5
               COL = COL / C6
               GO TO 170
            END IF
!                                                         ** NOW BALANCE
            IF ( (COL+ROW)/F .LT. C4*H ) THEN
               WKD(I)  = WKD(I) * F
               NOCONV = .TRUE.
               DO 180 J = L, M
                  AAD(I,J) = AAD(I,J) / F
180            CONTINUE
               DO 190 J = 1, K
                  AAD(J,I) = AAD(J,I) * F
190            CONTINUE
            END IF
200      CONTINUE

      IF ( NOCONV ) GO TO 140
!                                  ** IS -A- ALREADY IN HESSENBERG FORM?
      IF ( K-1 .LT. L+1 ) GO TO 350
!                                   ** TRANSFER -A- TO A HESSENBERG FORM
      DO 290 N = L+1, K-1
         H        = ZERO
         WKD(N+M) = ZERO
         SCALE    = ZERO
!                                                        ** SCALE COLUMN
         DO 210 I = N, K
            SCALE = SCALE + DABS(AAD(I,N-1))
210      CONTINUE
         IF ( SCALE.NE.ZERO ) THEN
            DO 220 I = K, N, -1
               WKD(I+M) = AAD(I,N-1) / SCALE
               H = H + WKD(I+M)**2
220         CONTINUE
            G = - SIGN( DSQRT(H), WKD(N+M) )
            H = H - WKD(N+M) * G
            WKD(N+M) = WKD(N+M) - G
!                                                 ** FORM (I-(U*UT)/H)*A
            DO 250 J = N, M
               F = ZERO
               DO 230  I = K, N, -1
                  F = F + WKD(I+M) * AAD(I,J)
230            CONTINUE
               DO 240 I = N, K
                  AAD(I,J) = AAD(I,J) - WKD(I+M) * F / H
240            CONTINUE
250         CONTINUE
!                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 280 I = 1, K
               F = ZERO
               DO 260  J = K, N, -1
                  F = F + WKD(J+M) * AAD(I,J)
260            CONTINUE
               DO 270 J = N, K
                  AAD(I,J) = AAD(I,J) - WKD(J+M) * F / H
270            CONTINUE
280         CONTINUE
            WKD(N+M)  = SCALE * WKD(N+M)
            AAD(N,N-1) = SCALE * G
         END IF
290   CONTINUE

      DO 340  N = K-2, L, -1
         N1 = N + 1
         N2 = N + 2
         F  = AAD(N1,N)
         IF ( F.NE.ZERO ) THEN
            F  = F * WKD(N1+M)
            DO 300 I = N2, K
               WKD(I+M) = AAD(I,N)
300         CONTINUE
            IF ( N1.LE.K ) THEN
               DO 330 J = 1, M
                  G = ZERO
                  DO 310 I = N1, K
                     G = G + WKD(I+M) * EVECD(I,J)
310               CONTINUE
                  G = G / F
                  DO 320 I = N1, K
                     EVECD(I,J) = EVECD(I,J) + G * WKD(I+M)
320               CONTINUE
330            CONTINUE
            END IF
         END IF
340   CONTINUE

350   CONTINUE
      N = 1
      DO 370 I = 1, M
         DO 360 J = N, M
            RNORM = RNORM + DABS(AAD(I,J))
360      CONTINUE
         N = I
         IF ( I.LT.L .OR. I.GT.K ) EVALD(I) = AAD(I,I)
370   CONTINUE
      N = K
      T = ZERO
!                                         ** SEARCH FOR NEXT EIGENVALUES
380   IF ( N.LT.L ) GO TO 530
      IN = 0
      N1 = N - 1
      N2 = N - 2
!                          ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
390   CONTINUE
      DO 400 I = L, N
         LB = N+L - I
         IF ( LB.EQ.L ) GO TO 410
         S = DABS( AAD(LB-1,LB-1) ) + DABS( AAD(LB,LB) )
         IF ( S.EQ.ZERO ) S = RNORM
         IF ( DABS(AAD(LB,LB-1)) .LE. TOL*S ) GO TO 410
400   CONTINUE

410   X = AAD(N,N)
      IF ( LB.EQ.N ) THEN
!                                        ** ONE EIGENVALUE FOUND
         AAD(N,N)  = X + T
         EVALD(N) = AAD(N,N)
         N = N1
         GO TO 380
      END IF

      Y = AAD(N1,N1)
      W = AAD(N,N1) * AAD(N1,N)
      IF ( LB.EQ.N1 ) THEN
!                                        ** TWO EIGENVALUES FOUND
         P = (Y-X) * C2
         Q = P**2 + W
         Z = DSQRT( DABS(Q) )
         AAD(N,N) = X + T
         X = AAD(N,N)
         AAD(N1,N1) = Y + T
!                                        ** REAL PAIR
         Z = P + SIGN(Z,P)
         EVALD(N1) = X + Z
         EVALD(N)  = EVALD(N1)
         IF ( Z.NE.ZERO ) EVALD(N) = X - W / Z
         X = AAD(N,N1)
!                                  ** EMPLOY SCALE FACTOR IN CASE
!                                  ** X AND Z ARE VERY SMALL
         R = SQRT( X*X + Z*Z )
         P = X / R
         Q = Z / R
!                                             ** ROW MODIFICATION
         DO 420 J = N1, M
            Z = AAD(N1,J)
            AAD(N1,J) = Q * Z + P * AAD(N,J)
            AAD(N,J)  = Q * AAD(N,J) - P * Z
420      CONTINUE
!                                             ** COLUMN MODIFICATION
         DO 430 I = 1, N
            Z = AAD(I,N1)
            AAD(I,N1) = Q * Z + P * AAD(I,N)
            AAD(I,N)  = Q * AAD(I,N) - P * Z
430      CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
         DO 440 I = L, K
            Z = EVECD(I,N1)
            EVECD(I,N1) = Q * Z + P * EVECD(I,N)
            EVECD(I,N)  = Q * EVECD(I,N) - P * Z
440      CONTINUE

         N = N2
         GO TO 380
      END IF

      IF ( IN.EQ.30 ) THEN
!                    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR
!                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE
         IER = N
         GO TO 670
      END IF
!                                                          ** FORM SHIFT
      IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
         T = T + X
         DO 450 I = L, N
            AAD(I,I) = AAD(I,I) - X
450      CONTINUE
         S = DABS(AAD(N,N1)) + DABS(AAD(N1,N2))
         X = C3 * S
         Y = X
         W = - C1 * S**2
      END IF

      IN = IN + 1
!                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS

      DO 460 J = LB, N2
         I = N2+LB - J
         Z = AAD(I,I)
         R = X - Z
         S = Y - Z
         P = ( R * S - W ) / AAD(I+1,I) + AAD(I,I+1)
         Q = AAD(I+1,I+1) - Z - R - S
         R = AAD(I+2,I+1)
         S = DABS(P) + DABS(Q) + DABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF ( I.EQ.LB ) GO TO 470
         UU = DABS( AAD(I,I-1) ) * ( DABS(Q) + DABS(R) )
         VV = DABS(P)*(DABS(AAD(I-1,I-1))+DABS(Z)+DABS(AAD(I+1,I+1)))
         IF ( UU .LE. TOL*VV ) GO TO 470
460   CONTINUE

470   CONTINUE
      AAD(I+2,I) = ZERO
      DO 480 J = I+3, N
         AAD(J,J-2) = ZERO
         AAD(J,J-3) = ZERO
480   CONTINUE

!             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N

      DO 520 KA = I, N1
         NOTLAS = KA.NE.N1
         IF ( KA.EQ.I ) THEN
            S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
            IF ( LB.NE.I ) AAD(KA,KA-1) = - AAD(KA,KA-1)
         ELSE
            P = AAD(KA,KA-1)
            Q = AAD(KA+1,KA-1)
            R = ZERO
            IF ( NOTLAS ) R = AAD(KA+2,KA-1)
            X = DABS(P) + DABS(Q) + DABS(R)
            IF ( X.EQ.ZERO ) GO TO 520
            P = P / X
            Q = Q / X
            R = R / X
            S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
            AAD(KA,KA-1) = - S * X
         END IF
         P = P + S
         X = P / S
         Y = Q / S
         Z = R / S
         Q = Q / P
         R = R / P
!                                                    ** ROW MODIFICATION
         DO 490 J = KA, M
            P = AAD(KA,J) + Q * AAD(KA+1,J)
            IF ( NOTLAS ) THEN
               P = P + R * AAD(KA+2,J)
               AAD(KA+2,J) = AAD(KA+2,J) - P * Z
            END IF
            AAD(KA+1,J) = AAD(KA+1,J) - P * Y
            AAD(KA,J)   = AAD(KA,J)   - P * X
490      CONTINUE
!                                                 ** COLUMN MODIFICATION
         DO 500 II = 1, MIN0(N,KA+3)
            P = X * AAD(II,KA) + Y * AAD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * AAD(II,KA+2)
               AAD(II,KA+2) = AAD(II,KA+2) - P * R
            END IF
            AAD(II,KA+1) = AAD(II,KA+1) - P * Q
            AAD(II,KA)   = AAD(II,KA) - P
500      CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
         DO 510 II = L, K
            P = X * EVECD(II,KA) + Y * EVECD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * EVECD(II,KA+2)
               EVECD(II,KA+2) = EVECD(II,KA+2) - P * R
            END IF
            EVECD(II,KA+1) = EVECD(II,KA+1) - P * Q
            EVECD(II,KA)   = EVECD(II,KA) - P
510      CONTINUE

520   CONTINUE
      GO TO 390
!                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR
530   CONTINUE
      IF ( RNORM.NE.ZERO ) THEN
         DO 560  N = M, 1, -1
            N2 = N
            AAD(N,N) = ONE
            DO 550  I = N-1, 1, -1
               W = AAD(I,I) - EVALD(N)
               IF ( W.EQ.ZERO ) W = TOL * RNORM
               R = AAD(I,N)
               DO 540 J = N2, N-1
                  R = R + AAD(I,J) * AAD(J,N)
540            CONTINUE
               AAD(I,N) = - R / W
               N2 = I
550         CONTINUE
560      CONTINUE
!                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS

         DO 580 I = 1, M
            IF ( I.LT.L .OR. I.GT.K ) THEN
               DO 570 J = I, M
                  EVECD(I,J) = AAD(I,J)
570            CONTINUE
            END IF
580      CONTINUE
!                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF ( K.NE.0 ) THEN
            DO 600  J = M, L, -1
               DO 600 I = L, K
                  Z = ZERO
                  DO 590 N = L, MIN0(J,K)
                     Z = Z + EVECD(I,N) * AAD(N,J)
590               CONTINUE
                  EVECD(I,J) = Z
600         CONTINUE
         END IF

      END IF

      DO 620 I = L, K
         DO 620 J = 1, M
            EVECD(I,J) = EVECD(I,J) * WKD(I)
620   CONTINUE
!                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO 640  I = L-1, 1, -1
         J = INT(WKD(I))
         IF ( I.NE.J ) THEN
            DO 630 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
630         CONTINUE
         END IF
640   CONTINUE

      DO 660 I = K+1, M
         J = INT(WKD(I))
         IF ( I.NE.J ) THEN
            DO 650 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
650         CONTINUE
         END IF
660   CONTINUE
!
  670 CONTINUE

      RETURN
      END SUBROUTINE  ASYMTX

!

      SUBROUTINE GETQUAD2(A,B,N,ROOTS,WGTS)

!  Computes N roots and weights for Gauss-Legendre quadrature on the interval (a,b)

      IMPLICIT NONE

!  Limits of interval

      REAL(FPK), INTENT(IN)  :: A, B

!  Dimension

      INTEGER, INTENT(IN) :: N

!  Quadrature roots and weights

      REAL(FPK), INTENT(OUT) :: ROOTS(N), WGTS(N)

!  Local variables

      INTEGER   :: I, M, N2, NM1
      REAL(FPK) :: IR, MR, NR
      REAL(FPK) :: MIDPT, SFAC
      REAL(FPK) :: DLP_DX, LP, LPM1, LPM2, X, XOLD, XX

!  Threshold for Newton's Method

      REAL(FPK), PARAMETER :: QEPS = 1.0D-13

!  Since roots are symmetric about zero on the interval (-1,1), split the interval
!  in half and only work on the lower half of the interval (-1,0).

      N2 = INT((N + 1)/2)
      NR = REAL(N,FPK)

!  Define the shift [midpoint of (a,b)] and scale factor to later move roots from
!  the interval (-1,1) to the interval (a,b)

      MIDPT = HALF*(B + A)
      SFAC  = HALF*(B - A)

      DO M = 1, N2

!  Find current root of the related Nth order Legendre Polynomial on (-1,0) by Newton's
!  Method using two Legendre Polynomial recurrence relations (e.g. see Abramowitz &
!  Stegan (1972))

         !Define starting point [ after Tricomi (1950) ]
         MR = REAL(M,FPK)
         XX = PIE*(MR - QUARTER)/(NR + HALF)
         X  = (ONE - (NR - ONE)/(8.0_FPK*NR**3) &
             - ONE/(384.0_FPK*NR**4)*(39.0_FPK - 28.0_FPK/SIN(XX)**2))*COS(XX)

         !Use Newton's Method
         DO 
            LPM1 = ZERO ; LP = ONE
            DO I = 1, N
               IR = REAL(I,FPK) ; LPM2 = LPM1 ; LPM1 = LP
               LP = ((TWO*IR - ONE)*X*LPM1 - (IR - ONE)*LPM2)/IR
            ENDDO
            DLP_DX = NR*(X*LP - LPM1)/(X**2 - ONE)
            XOLD = X ; X = XOLD - LP/DLP_DX
            IF (ABS(X-XOLD) <= QEPS) EXIT
         ENDDO

!  Shift and scale the current root (and its symmetric counterpart) from the interval (-1,1)
!  to the interval (a,b).  Define their related weights (e.g. see Abramowitz & Stegan (1972)).
!  Note:
!  If (1) N is even or (2) N is odd and M /= N2, then ROOTS(M) and ROOTS(NM1) are unique.
!  If N is odd and M = N2, then M = NM1 and ROOTS(M) = ROOTS(NM1) are one and the same root.

         !On interval lower half: (a,midpt)
         ROOTS(M)   = MIDPT - SFAC*X
         WGTS(M)    = (TWO*SFAC)/((ONE - X**2)*DLP_DX**2)

         !On interval upper half: (midpt,b)
         NM1 = N - M + 1
         ROOTS(NM1) = MIDPT + SFAC*X
         WGTS (NM1) = WGTS(M)

      ENDDO

      END SUBROUTINE GETQUAD2

!

      SUBROUTINE ISSORT(N,ILIST)

!  Sorts a 1-D array of integers in ascending order using a simple selection
!  ("interchange") sort if necessary.  Can handle cases where elements in the
!  array appear more than once.

      IMPLICIT NONE

!  Dimension

      INTEGER, INTENT(IN) :: N

!  List to sort

      INTEGER, INTENT(INOUT) :: ILIST(N)

!  Local variables

      INTEGER :: I, ITEMP, J, SMALL_INDEX
      LOGICAL :: IN_ORDER

!  Do prechecks

      !(1) Case where N = 1
      IF (N == 1) RETURN

      !(2) Case where list is already in order
      IN_ORDER = .TRUE.
      DO I = 1, N-1
        IF (ILIST(I) > ILIST(I+1)) THEN
          IN_ORDER = .FALSE.
          EXIT
        ENDIF
      ENDDO
      IF (IN_ORDER) RETURN

!  Sort

      DO I = 1, N-1
        SMALL_INDEX = I
        DO J = I+1, N
          IF (ILIST(J) < ILIST(SMALL_INDEX)) SMALL_INDEX = J
        ENDDO
        IF (SMALL_INDEX /= I) THEN
          ITEMP              = ILIST(I)
          ILIST(I)           = ILIST(SMALL_INDEX)
          ILIST(SMALL_INDEX) = ITEMP
        ENDIF
      ENDDO

      END SUBROUTINE ISSORT

!

      SUBROUTINE RSSORT(N,RLIST)

!  Sorts a 1-D array of reals in ascending order using a simple selection
!  ("interchange") sort if necessary.  Can handle cases where elements in the
!  array appear more than once.

      IMPLICIT NONE

!  Dimension

      INTEGER, INTENT(IN) :: N

!  List to sort

      REAL(FPK), INTENT(INOUT) :: RLIST(N)

!  Local variables

      INTEGER   :: I, J, SMALL_INDEX
      REAL(FPK) :: RTEMP
      LOGICAL   :: IN_ORDER

!  Do prechecks

      !(1) Case where N = 1
      IF (N == 1) RETURN

      !(2) Case where list is already in order
      IN_ORDER = .TRUE.
      DO I = 1, N-1
        IF (RLIST(I) > RLIST(I+1)) THEN
          IN_ORDER = .FALSE.
          EXIT
        ENDIF
      ENDDO
      IF (IN_ORDER) RETURN

!  Sort

      DO I = 1, N-1
        SMALL_INDEX = I
        DO J = I+1, N
          IF (RLIST(J) < RLIST(SMALL_INDEX)) SMALL_INDEX = J
        ENDDO
        IF (SMALL_INDEX /= I) THEN
          RTEMP              = RLIST(I)
          RLIST(I)           = RLIST(SMALL_INDEX)
          RLIST(SMALL_INDEX) = RTEMP
        ENDIF
      ENDDO

      END SUBROUTINE RSSORT

!

      RECURSIVE SUBROUTINE RQSORT_IDX(M,RLIST,N,IDX)

!  Takes a 1-D input array of reals RLIST and returns a 1-D array of integers IDX
!  whose entries are the indices of RLIST in an order needed to put the elements
!  of RLIST in ascending order.  Uses a "quicksort" sort if necessary.  Can handle
!  cases where elements in RLIST appear more than once.

!  Notes: * Before initial call, output IDX must be initialized with the integers
!           [1,2,...,N].
!         * On initial call, M = N.

      IMPLICIT NONE

!  Dimension of list

      INTEGER, INTENT(IN) :: M

!  List whose indices may require sorting

      REAL(FPK), INTENT(IN) :: RLIST(M)

!  Dimension of sorted index sublist

      INTEGER, INTENT(IN) :: N

!  Sorted index sublist

      INTEGER, INTENT(INOUT) :: IDX(N)

!  Local variables

      INTEGER   :: ITEMP, LE, NSUB1, NSUB2, PIVOT, RE

!  Do precheck

      !Case where N = 1
      IF (N == 1) RETURN

!  Sort

!  (1) Partition the current list

      LE = 1 ; RE = N
      PIVOT = INT((1 + N)/2)
      DO WHILE (LE < RE) !If >= , fall out of loop

        !Moving right, look for a key >= the pivot
        DO WHILE (RLIST(IDX(LE)) <  RLIST(IDX(PIVOT)))
          LE = LE + 1
        ENDDO

        !Moving left, look for a key < the pivot
        DO WHILE ((RLIST(IDX(RE)) >= RLIST(IDX(PIVOT))) .AND. (LE < RE))
          RE = RE - 1
        ENDDO

        IF (LE < RE) THEN
          !Case where LE < PIVOT < RE or LE < RE < PIVOT:
          !     swap IDX(LE) & IDX(RE)
          ITEMP   = IDX(LE)
          IDX(LE) = IDX(RE)
          IDX(RE) = ITEMP

          !Case where the LE = PIVOT < RE:
          !     (1) swap IDX(LE) & IDX(RE) (done above)
          !     (2) set PIVOT = RE         (done here)
          IF (PIVOT == LE) PIVOT = RE
        ELSE
          !Case where LE = RE < PIVOT:
          !     (1) swap IDX(LE) & IDX(PIVOT) (done here)
          !     (2) set LE = LE + 1           (done below)
          IF (LE < PIVOT) THEN
            ITEMP      = IDX(LE)
            IDX(LE)    = IDX(PIVOT)
            IDX(PIVOT) = ITEMP
          ENDIF

          !Case where LE = PIVOT = RE or LE = RE = PIVOT:
          !     set LE = LE + 1
          LE = LE + 1
        ENDIF
      ENDDO

      !Note: pivot ends up on the right side of the left sublist

!  (2) Sort the two sublists

      DO 
        NSUB1 = RE
        NSUB2 = N - LE + 1
        IF (NSUB2 == 0) THEN
          !Case where only the left sublist exists
          LE = LE - 1
          RE = RE - 1
        ELSE
          EXIT
        ENDIF
      ENDDO

      !For safety
      IF (NSUB1 == 0 .OR. NSUB2 == 0) THEN
        WRITE(*,*)
        WRITE(*,*) 'NSUB1 = ',NSUB1
        WRITE(*,*) 'NSUB2 = ',NSUB2
        WRITE(*,*) 'ERROR IN RQSORT_IDX! EITHER NSUB1 OR NSUB2 = 0'
        WRITE(*,*)
        STOP 'Safety-first stoppage in RQSORT_IDX'
      ENDIF

      CALL RQSORT_IDX(M,RLIST,NSUB1,IDX(1:RE))
      CALL RQSORT_IDX(M,RLIST,NSUB2,IDX(LE:N))

      END SUBROUTINE RQSORT_IDX

!

      SUBROUTINE CFPLGARR (MAXMOM, L, M, X, CFPLG)

      implicit none

!-------- This Routine Calculates  -----------

! ###########################################################
! #                                                         #
! #          [(l - m)! / (l + m)!]^{1/2} P(l,m)(x)          #
! #                                                         #
! ###########################################################

!     Based on the PLG routine from the numerical receipes,
!     this revised routine calculates all legendre polynomials
!     from Plg(m,m) to Plg(l,m).
!     The maximum order for the Plg is OMAX


!     ---------------
!     Input Variables
!     ---------------

      Integer            maxmom
      Integer            l, m       ! order of legendre polynomial
      Double Precision   x          ! abscissa


!     ---------------
!     Local Variables
!     ---------------

      Integer            i, ll

      Double Precision   fact, fll, fmm, fmmp1, somx2

!     ---------------
!     Output Variables
!     ---------------

      Double Precision   cfplg (0:maxmom) ! array of Plgs

!     --- Compute F(M,M)  ---

      fmm = 1.0d0

      If (m .gt. 0) Then
         somx2 = Dsqrt( (1.0d0-x) * (1.0d0+x) )
         fact = 1.0d0
         Do i = 1, m
            fmm = -fmm*fact*somx2 / Dsqrt( DBLE( (2*i)*(2*i-1) ) )
            fact = fact + 2.0d0
         End Do
      End If

      cfplg(m) = fmm

!     --- Compute F(M+1,M)  ---

      If (l .gt. m) Then
         fmmp1 = x * Dsqrt( DBLE( 2*m+1 ) ) * fmm
         cfplg(m+1) = fmmp1

!     --- Compute F(L,M) ---

         If ( l .gt. (m+1) ) Then

            Do ll = m+2, l
               fll = (fmmp1 * x * DBLE(2*ll-1) - &
                    fmm * DBLE(ll+m-1) &
                    * Dsqrt( (DBLE( ll-m-1 ) ) / DBLE( ll+m-1) ) ) &
                    / Dsqrt( DBLE( (ll+m)*(ll-m) ) )
               fmm   = fmmp1
               fmmp1 = fll
               cfplg(ll) = fmmp1
            End Do

         End If

      End If

      Return
      END SUBROUTINE CFPLGARR

!  FINDPAR input subroutines.

!----- The basic idea

!  FINDPAR is a poor man's RDPAR.  It implements some of the
! central features of RDPAR --- the ability to find a named set
! of data values in an input file, the relaxation of the need for
! data in the input file to appear in the order in which it's
! used, the provision for comments in the input file --- in
! (almost) standard FORTRAN 77.

!  Among the features of RDPAR that are not implemented in FINDPAR:
! all redirection of input (from one parameter name to another, from
! file to terminal, from file to command line); freer-format data
! specification on input; tracing options.

!----- Using FINDPAR

!  To use FINDPAR, you open your input file using ordinary FORTRAN
! procedures.  You will also use ordinary FORTRAN to read data
! values.  The only difference is that you can use the FINDPAR call
! at any point to position the file to a point after a line containing
! a parameter name.  These parameter names help to document what the
! input values are for; they also don't have to appear in the order in
! which the program uses them.  An input file might look like this:

! ! Parameters for July 16 run.

! NLAYERS
! 5

! LAYER PARAMETERS
! 5,1013,296
! .001,.001,.001,.001,.001,.001,.001,.001

! MOLECULE FLAGS
! t,t,t,t,t,t,t,t

!  By convention, ! or a space introduces a comment line.  A
! parameter name appears on a line of its own, followed by the data
! lines to be read using FORTRAN statements.  Here is a sample of
! FORTRAN code to read this parameter file:

!      OPEN (UNIT = IN_UNIT, FILE = INPUT, STATUS = 'OLD')
!      CALL FINDPAR ( IN_UNIT, 'NLAYERS' )
!      READ (IN_UNIT, *) NLAYERS
!      CALL FINDPAR ( IN_UNIT, 'LAYER PARAMETERS' )
!      READ (IN_UNIT, *) Z, PRESS, TEMP
!      READ (IN_UNIT, *) BRO, CLO, HCHO, NO2, O2, O3, OCLO, SO2
!      CALL FINDPAR ( IN_UNIT, 'MOLECULE FLAGS' )
!      READ (IN_UNIT, *) IF_BRO, IF_CLO, IF_HCHO, IF_NO2, IF_O2, IF_O3,
!     $ IF_OCLO, IF_SO2
!      CLOSE (UNIT = IN_UNIT)

!----- Details of parameter file format

! The parameter file consists of a sequence of comments and
! parameter-name/data chunks, in any order (except as noted below).

! Comments are essentially anything that doesn't happen to get matched
! as a parameter name or read as data.  However, by convention comments
! are introduced by an exclamation point, to make absolutely sure that
! they don't get mistaken for parameter names.  It is also conventional
! to ``comment out'' a parameter name that you don't want to match by
! simply indenting it one space.

! Parameter names can theoretically be almost any sequence of
! characters.  However, FINDPAR ignores trailing spaces both in the
! parameter name you give it and in lines it reads from the parameter
! file.  And, by convention, parameter names don't begin with spaces.
! Parameter names may contain embedded spaces.  Parameter name matching
! is case-insensitive but otherwise exact: if you use odd things
! like control characters in a parameter name they have to match
! exactly, too, so it's usually best to avoid control characters and
! tabs in parameter names.

! The data lines that follow a line with a parameter name are ordinary
! FORTRAN input lines.  FINDPAR has no influence on how they're handled,
! and all the ordinary rules about FORTRAN READ statements apply.

! Comments may not appear on the same lines as parameter names or data,
! and they can't appear within parameter name/data chunks.  When FINDPAR
! locates a parameter name it positions the file at the line following
! that name; your program will then use ordinary FORTRAN READ statements
! to read the data, and will probably be confused by comments that appea
! at that point.

! There's a maximum length for input lines and parameter names, set
! by the parameter MAXLINE in the code.

! By convention, a parameter name that ends in ? precedes an
! input line that contains a single logical variable (T or F).

!----- Repeated parameters

! Sometimes you want to read several parameters with the same name:
! for example, parameters for each layer of the atmosphere, where there
! might be a large number of layers.  In this case the order of
! parameters in the file does matter, since you usually just want to
! put the numbers in the file in some natural order rather than giving
! a number in each input line specifying where it appears in the input
! order.

! The normal FINDPAR approach is to organize these parameters as
! follows: begin with some parameter that appears only once in the file
! (the number of layers, for example); then follow it with the instances
! of the repeated parameter and associated data.  The only-once
! parameter can even be a dummy, with no data following; its importance
! is that reading it gets the parameter file positioned at the right
! point.

! The other problem here is knowing when a list of repeated parameters
! is over.  FINDPAR doesn't have any of the tricks for doing this that
! RDPAR has; you must either know exactly how many instances of the
! repeated parameter there are going to be (by reading some once-only
! parameter that tells you), or else you need some special numbers to
! flag the end of the list (zeros, negative numbers, etc.).  You can't
! just keep looking for the same parameter name indefinitely, because
! FINDPAR will just rewind the file for you and loop through its
! contents forever.

!----- Errors and XFINDPAR

! If FINDPAR can't find a parameter, it prints an error message
! and returns.  Your READ statements following the FINDPAR call are
! very likely to run into errors of their own at this point, since
! you'll be positioned at some random point in the file (usually at
! the beginning, in this version, but that isn't true in all cases).

! If you want to do something more intelligent about such errors,
! you can use XFINDPAR, which is a function returning a logical value:
! .true. if the parameter name was found, .false. if not; XFINDPAR
! doesn't display FINDPAR's error message when a parameter name cannot
! be found.

! Both FINDPAR and XFINDPAR can run into ordinary FORTRAN input errors
! as they read the file: no special action is taken on these---the
! system's default action, whatever that is, occurs.

!       9/10/91         John Lavagnino

!  Find a parameter name in a parameter file.

!        SUBROUTINE FINDPAR ( NUNIT, PARNAME )

!  Input arguments.

!        INTEGER         NUNIT
!        CHARACTER * (*) PARNAME

!  Local variables.

!        LOGICAL         XFINDPAR
!        EXTERNAL        XFINDPAR

!  No local variables need to be SAVEd.

!***********************************************************************

!        IF ( .NOT. XFINDPAR ( NUNIT, PARNAME ) ) THEN
!          WRITE (*, *) 'FINDPAR error'
!          WRITE (*, *) '   Unit number ', NUNIT
!          WRITE (*, *) '   Parameter name ', PARNAME
!        END IF

!        END
!
!  Find a parameter name in a paramter file: return .true. if found,
! .false. if not.

      LOGICAL FUNCTION XFINDPAR ( NUNIT, PARNAME )

        implicit none

!  Parameters.

        INTEGER         MAXLINE
        PARAMETER       ( MAXLINE = 132 )

!  Input arguments.

        INTEGER         NUNIT
        CHARACTER * (*) PARNAME

!  Local variables.

        INTEGER                 III
        INTEGER                 PARLEN
        INTEGER                 START
        INTEGER                 END
        LOGICAL                 ENDSEEN
        CHARACTER * (MAXLINE)   LINE
        CHARACTER * (MAXLINE)   NAME

!  No local variables need to be SAVEd.

!***********************************************************************

!  Determine the length of the parameter name.

        DO 10 III = LEN ( PARNAME ), 1, -1
          IF ( PARNAME ( III : III ) .NE. ' ' ) THEN
            PARLEN = III
            GO TO 20
          END IF
10      CONTINUE

!  If we get to here, then name contains nothing but blanks.
!  We just return, claiming success, in such a case.

        GO TO 500

!  If we get here, then there's a non-null parameter name; but it
!  might still be too long, in which case we always return
!  signaling failure, since we couldn't ever find such a name in the
!  file.

20      CONTINUE
        IF ( PARLEN  .GT.  MAXLINE ) GO TO 400

!  Convert the name to lower-case.

        NAME = PARNAME ( 1 : PARLEN )
        CALL LCSTRING ( NAME ( 1 : PARLEN ) )

!  Top of main loop.

        ENDSEEN = .FALSE.

100     CONTINUE

          LINE = ' '
          READ ( UNIT = NUNIT, FMT = '(A)', END = 200 ) LINE
          CALL LCSTRING ( LINE ( 1 : PARLEN ) )
          IF ( LINE ( 1 : PARLEN ) .NE. NAME ( 1 : PARLEN ) ) &
                GO TO 100

          START = PARLEN + 1
          END = MAXLINE
          IF ( START  .GT.  END ) GO TO 500
          IF ( LINE ( START : END ) .EQ. ' ' ) GO TO 500
          GO TO 100

!  End-of-file branch.

200       CONTINUE
          REWIND ( UNIT = NUNIT )
          IF ( ENDSEEN ) THEN
            GO TO 400
          ELSE
            ENDSEEN = .TRUE.
            GO TO 100
          END IF

!  End of loop: failure, no parameter name found.

400     CONTINUE
        XFINDPAR = .FALSE.
        RETURN

!  End of loop: successful location of parameter name.

500     CONTINUE
        XFINDPAR = .TRUE.
        RETURN

      END FUNCTION XFINDPAR

!

      LOGICAL FUNCTION GFINDPAR (NUNIT, PREFIX, &
                                          ERROR, PARNAME)

!  Find a parameter name in a parameter file, with optional added
! prefix on parameter name.  Error messages are displayed but
! instead of stopping the program a logical flag is returned.
!  Returns .TRUE. if some form of the parameter name was found;
! .FALSE. if not.
! ERROR and GFINDPAR are initialised for a successful search.
! ERROR is just the negation of GFINDPAR. Error messages are removed.

        implicit none

!  Input arguments.  If PREFIX is something other than a bunch of
! blanks, it will be added onto PARNAME, with a blank to separate
! it, for the FINDPAR call.  If that FINDPAR call fails, then
! PARNAME without prefix is tried instead, and only if that call fails
! is there an error message.  This means that unprefixed versions of
! the parameter name act as defaults.
!  If PREFIX is a bunch of blanks, this works much as FINDPAR does.
!  (This is assuming that PREFIX has no trailing blanks.  Probably
! it should be checking for them and trimming them if necessary.)

        INTEGER         NUNIT
        CHARACTER * (*) PREFIX
        CHARACTER * (*) PARNAME

!  Modified argument.  If all goes well, this is not changed.  If
! the parameter name can't be found, it's set to .TRUE.

        LOGICAL         ERROR

!  Local variables.

        CHARACTER * 132 LINE

!        LOGICAL         XFINDPAR
!        EXTERNAL        XFINDPAR

!***********************************************************************

!  Initialise return values.

        GFINDPAR = .TRUE.
        ERROR    = .FALSE.

!  Compose search line for first search.

        IF ( PREFIX  .EQ.  ' ' ) THEN
          LINE = PARNAME
        ELSE
          LINE = PREFIX // ' ' // PARNAME
        END IF

!  Try the parameter name with prefix.

        IF (XFINDPAR (NUNIT, LINE)) THEN
          RETURN
        END IF

!  If that doesn't work, try it without the prefix, if the prefix
! was non-null.

        IF (PREFIX  .NE.  ' ') THEN
          IF (XFINDPAR (NUNIT, PARNAME)) THEN
            RETURN
          END IF
        END IF

!  If we got here, we just couldn't find the parameter name.

        GFINDPAR = .FALSE.
        ERROR    = .TRUE.

      END FUNCTION GFINDPAR

!

      SUBROUTINE FINDPAR_ERROR &
           (  ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  include file

      USE VLIDORT_PARS_m, only : MAX_MESSAGES, VLIDORT_SERIOUS

      implicit none

!  subroutine input arguments

      LOGICAL       ERROR
      CHARACTER*(*) PAR_STR

!  subroutine Output or In/Out arguments

      INTEGER       STATUS
      INTEGER       NM
      CHARACTER*(*) MESSAGES(0:MAX_MESSAGES), ACTIONS(0:MAX_MESSAGES)

!  local variables

!      INTEGER       LEN_STRING
!      EXTERNAL      LEN_STRING

      IF ( ERROR ) THEN
        NM = NM + 1
        STATUS = VLIDORT_SERIOUS
        MESSAGES(NM) = 'Cannot find string: '//Adjustl(Trim(PAR_STR))
        ACTIONS(NM)  = 'Check Spelling of String in input file'
      ENDIF

!  finish

      RETURN
      END SUBROUTINE FINDPAR_ERROR

!  End module

      END MODULE vlidort_aux_m

