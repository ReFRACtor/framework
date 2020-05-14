Module pca_auxiliaries_m

!  Collection of Numerical routines

use iso_c_binding

private
public PCA_ASYMTX, PCA_Ranker, PCA_LINTP2, PCA_DGETRF, PCA_DGETRS

!  precision parameters

      INTEGER, PARAMETER :: sp     = selected_real_kind(6)
      INTEGER, PARAMETER :: dp     = selected_real_kind(15)

contains

! ###############################################################
! #                                                             #
! #  pca_auxiliaries.f90. Collection of Numerical Routines.     #
! #                                                             #
! #     PCA_Ranker   (same as Indexx1 from Numerical Recipes)   #  
! #     PCA_ASYMTX   (adapted from LIDORT routine)              #
! #     PCA_LINTP2   (Linear interpolation)                     #
! #                                                             #
! #   utility  modules compiled from LAPACK:                    #
! #                                                             #
! #      PCA_dgetrf:                  LAPACK                    #
! #      PCA_dgetrs:                  LAPACK                    #
! #      PCA_dlaswp:                  LAPACK                    #
! #      PCA_ilaenv:                  LAPACK                    #
! #      PCA_lsame:                   LAPACK                    #
! #      PCA_xerbla:                  LAPACK                    #
! #      PCA_dgetf2:                  LAPACK                    #
! #      PCA_dcopy:                   LAPACK                    #
! #      PCA_dgemm:                   LAPACK                    #
! #      PCA_dgemv:                   LAPACK                    #
! #      PCA_dger:                    LAPACK                    #
! #      PCA_dscal:                   LAPACK                    #
! #      PCA_dswap:                   LAPACK                    #
! #      PCA_dtrsm:                   LAPACK                    #
! #      PCA_idamax:                  LAPACK                    #
! #                                                             #
! ###############################################################

SUBROUTINE PCA_ASYMTX( AAD, M, IA, IEVEC, IEVEC2, TOL, EVECD, EVALD, IER, WKD, &
                       MESSAGE_LEN, MESSAGE_OUT, BAD_STATUS ) bind(C)

   implicit none

   real(kind=dp), parameter :: zero = 0.0d0
   real(kind=dp), parameter :: one  = 1.0d0

!    =======  D O U B L E    P R E ! I S I O N    V E R S I O N  ======

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
!    IEVE!    :  first dimension of  EVECD

!   O U T P U T    V A R I A B L E S:

!       EVECD :  (unnormalized) eigenvectors of  A 
!                   ( column J corresponds to EVALD(J) )

!       EVALD :  (unordered) eigenvalues of  A ( dimension at least M )

!       IER   :  if .NE. 0, signals that EVALD(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.

      LOGICAL(c_bool)       , intent(out) :: BAD_STATUS
      integer(c_int), intent(in) :: message_len
      character(kind=c_char), intent(out) :: MESSAGE_OUT(MESSAGE_LEN+1)

!   S ! R A T ! H   V A R I A B L E S:

!       WKD    :  WORK AREA ( DIMENSION AT LEAST 2*M )

!  input/output arguments

      INTEGER(c_int), intent(in)         :: M, IA, IEVEC, IEVEC2
      INTEGER(c_int), intent(out)        :: IER

      REAL(kind=c_double), intent(in)    :: TOL

      REAL(kind=c_double), intent(inout) :: AAD(IEVEC,IEVEC)
      REAL(kind=c_double), intent(inout) :: WKD(IEVEC2)

      REAL(kind=c_double), intent(out)   :: EVALD(IEVEC)
      REAL(kind=c_double), intent(out)   :: EVECD(IEVEC,IEVEC)

!  local variables (explicit declaration

      LOGICAL          :: NOCONV, NOTLAS
      INTEGER          :: I, J, L, K, KKK, LLL
      INTEGER          :: N, N1, N2, IN, LB, KA, II
      REAL(kind=dp)        :: C1, C2, C3, C4, C5, C6
      DATA               C1 / 0.4375D0 /
      DATA               C2 / 0.5D0 /
      DATA               C3 / 0.75D0 /
      DATA               C4 / 0.95D0 /
      DATA               C5 / 16.0D0 /
      DATA               C6 / 256.0D0 /
!      REAL(kind=dp)        :: TOL, DISCRI, SGN, RNORM, W, F, G, H, P, Q, R
      REAL(kind=dp)        :: DISCRI, SGN, RNORM, W, F, G, H, P, Q, R
      REAL(kind=dp)        :: REPL, COL, ROW, SCALE, T, X, Z, S, Y, UU, VV

!  these variables by R. Spurr, Removes GOTO statements
!     Introduced 08 Feburary 2010, RT SOLUTIONS Inc.

      LOGICAL          :: B70, B120, B200, B380, B390, B400, B460
      LOGICAL          :: TBLOCK, HESSENBERG_FORM, B490, B520

      character(kind=c_char, len=message_len) :: message_lcl
      integer :: len_idx

!  output status

      IER         = 0
      BAD_STATUS  = .FALSE.
      MESSAGE_LCL = ' '

!       Here change to bypass D1MACH:

!      TOL = 0.0000001
!       TOL = 1.0D-6
!       TOL = 1.0D-12
!       TOL = 1.0D-48
!        TOL = D1MACH(4)

      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M ) THEN
        MESSAGE_LCL = 'ASYMTX--bad input variable(s)'
        do len_idx = 1, message_len
          message_out(len_idx) = message_lcl(len_idx:len_idx)
        end do

        BAD_STATUS = .TRUE.
        RETURN
      ENDIF
!                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES

      IF ( M.EQ.1 )  THEN
         EVALD(1) = AAD(1,1)
         EVECD(1,1) = ONE
         RETURN
      ELSE IF ( M.EQ.2 )  THEN
         DISCRI = ( AAD(1,1) - AAD(2,2) )**2 + 4.0D0*AAD(1,2)*AAD(2,1)
         IF ( DISCRI.LT.ZERO ) THEN
           MESSAGE_LCL = 'ASYMTX--COMPLEX EVALS IN 2X2 CASE'
           do len_idx = 1, message_len
             message_out(len_idx) = message_lcl(len_idx:len_idx)
           end do

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

      TBLOCK = .TRUE.    
      DO WHILE (TBLOCK)

         KKK = K 
         B70 = .TRUE.
         J = KKK + 1

         DO WHILE (B70.and.TBLOCK)

            J = J - 1
            TBLOCK = (J.GT.1)

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
               B70 = .false.
  
            ENDIF

         ENDDO
      ENDDO

!                                     ** SEARCH FOR COLUMNS ISOLATING AN
!                                       ** EIGENVALUE AND PUSH THEM LEFT

      TBLOCK = .TRUE.    
      DO WHILE (TBLOCK)

         LLL = L
         B120 = .TRUE.
         J = LLL - 1

         DO WHILE (B120.and.TBLOCK)

            J = J + 1
            TBLOCK = (J.LT.K)

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
               B120 = .false.
             ENDIF

           ENDDO

       ENDDO

!      do i = 1, 4
!       write(*,'(4f10.6)')(AAD(I,J),J=1,4)
!      enddo
!       pause'after 120'

!                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K

      DO 130 I = L, K
         WKD(I) = ONE
130   CONTINUE

      B200 = .true.

      DO WHILE ( B200 )

         noconv = .false.

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

            DO WHILE ( COL .LT. G )
               F   = F * C5
               COL = COL * C6
            ENDDO

            G = ROW * C5
            DO WHILE ( COL .GE. G )
               F   = F / C5
               COL = COL / C6
            ENDDO

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

         B200 = ( noconv.eqv..true.)

      ENDDO

!      do i = 1, 4
!       write(*,'(4f10.6)')(AAD(I,J),J=1,4)
!      enddo
!       pause'before 350'

!  ** IS -A- ALREADY IN HESSENBERG FORM?

      HESSENBERG_FORM = ( K-1 .LT. L+1 )

      IF ( .not. HESSENBERG_FORM ) THEN

!   ** TRANSFER -A- TO A HESSENBERG FORM

         DO 290 N = L+1, K-1
            H        = ZERO
            WKD(N+M) = ZERO
            SCALE    = ZERO
!                                                        ** SCALE COLUMN
            DO 210 I = N, K
               SCALE = SCALE + DABS(AAD(I,N-1))
210         CONTINUE
            IF ( SCALE.NE.ZERO ) THEN
               DO 220 I = K, N, -1
                  WKD(I+M) = AAD(I,N-1) / SCALE
                  H = H + WKD(I+M)**2
220            CONTINUE
               G = - SIGN( DSQRT(H), WKD(N+M) )
               H = H - WKD(N+M) * G
               WKD(N+M) = WKD(N+M) - G
!                                                 ** FORM (I-(U*UT)/H)*A
               DO 250 J = N, M
                  F = ZERO
                  DO 230  I = K, N, -1
                     F = F + WKD(I+M) * AAD(I,J)
230               CONTINUE
                  DO 240 I = N, K
                     AAD(I,J) = AAD(I,J) - WKD(I+M) * F / H
240               CONTINUE
250            CONTINUE
!                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
               DO 280 I = 1, K
                  F = ZERO
                  DO 260  J = K, N, -1
                     F = F + WKD(J+M) * AAD(I,J)
260               CONTINUE
                  DO 270 J = N, K
                     AAD(I,J) = AAD(I,J) - WKD(J+M) * F / H
270               CONTINUE
280            CONTINUE
               WKD(N+M)  = SCALE * WKD(N+M)
               AAD(N,N-1) = SCALE * G
            END IF
290      CONTINUE

         DO 340  N = K-2, L, -1
            N1 = N + 1
            N2 = N + 2
            F  = AAD(N1,N)
            F  = F * WKD(N1+M)
            IF (F.NE.ZERO) THEN
!               F  = F * WKD(N1+M)
               DO 300 I = N2, K
                  WKD(I+M) = AAD(I,N)
300            CONTINUE
               IF ( N1.LE.K ) THEN
                  DO 330 J = 1, M
                     G = ZERO
                     DO 310 I = N1, K
                        G = G + WKD(I+M) * EVECD(I,J)
310                  CONTINUE
                     G = G / F
                     DO 320 I = N1, K
                        EVECD(I,J) = EVECD(I,J) + G * WKD(I+M)
320                  CONTINUE
330               CONTINUE
               END IF
            END IF
340      CONTINUE

!  End clause for conversion to Hessenberg Form

      ENDIF

!      do i = 1, 4
!       write(*,'(4f10.6)')(AAD(I,J),J=1,4)
!      enddo
!      do i = 1, 4
!       write(*,'(4f10.6)')(EVECD(I,J),J=1,4)
!      enddo
!       pause'before 350'

!  Set first Eigenvalues

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

! #################################################################
! #################################################################
! ################################################################

!    ** SEARCH FOR NEXT EIGENVALUES

      B380 = ( N.GE.L )
      DO WHILE (B380 )

         IN = 0
         N1 = N - 1
         N2 = N - 2

!  ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
 
         B390 = .true.
         DO WHILE ( B390 )

! 400  Block

            I = L - 1
            B400 = .true.
            do while ( B400 .and. I.lt.N )
               I = I + 1
               LB = N+L - I
               IF ( LB.EQ.L ) B400 = .false.
               if ( B400 ) then
                 S = DABS( AAD(LB-1,LB-1) ) + DABS( AAD(LB,LB) )
                 IF ( S.EQ.ZERO ) S = RNORM
                 IF ( DABS(AAD(LB,LB-1)) .LE. TOL*S ) B400 = .false.
               endif
            enddo

!   ** ONE EIGENVALUE FOUND

            X = AAD(N,N)      
            IF ( LB.EQ.N ) THEN
               AAD(N,N)  = X + T
               EVALD(N) = AAD(N,N)
               N = N1
               B390 = .false.
               B380 = ( N.GE.L )
            END IF

!  Other wise ** TWO EIGENVALUES FOUND

            IF ( B390 ) THEN

               Y = AAD(N1,N1)
               W = AAD(N,N1) * AAD(N1,N)

              IF ( LB.EQ.N1 ) THEN
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
420               CONTINUE
!                                             ** COLUMN MODIFICATION
                  DO 430 I = 1, N
                     Z = AAD(I,N1)
                     AAD(I,N1) = Q * Z + P * AAD(I,N)
                     AAD(I,N)  = Q * AAD(I,N) - P * Z
430               CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
                  DO 440 I = L, K
                     Z = EVECD(I,N1)
                     EVECD(I,N1) = Q * Z + P * EVECD(I,N)
                     EVECD(I,N)  = Q * EVECD(I,N) - P * Z
440               CONTINUE

                  N = N2
                  B390 = .false.
                  B380 = ( N.GE.L )
               END IF

            ENDIF

!  Carry on

            if ( B390 ) THEN

!    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR & return
!    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE

               IF ( IN.EQ.30 ) THEN
                  IER = N
                  RETURN
               END IF

!     ** FORM SHIFT

               IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
                   T = T + X
                   DO 450 I = L, N
                      AAD(I,I) = AAD(I,I) - X
450                CONTINUE
                   S = DABS(AAD(N,N1)) + DABS(AAD(N1,N2))
                   X = C3 * S
                   Y = X
                   W = - C1 * S**2
                END IF

                IN = IN + 1

!  ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS

                B460 = .true.
                J = LB - 1 
                DO WHILE (B460.and.J.LT.N2)
                   J = J + 1
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
                   IF ( I.EQ.LB ) B460 = .false.
                   if ( B460 ) then
                      UU = DABS( AAD(I,I-1) ) * ( DABS(Q) + DABS(R) )
                      VV = DABS(P) * &
                       (DABS(AAD(I-1,I-1))+DABS(Z)+DABS(AAD(I+1,I+1)))
                     IF ( UU .LE. TOL*VV ) B460 = .false.
                   endif
                ENDDO

!  B470 was here

               AAD(I+2,I) = ZERO
               DO 480 J = I+3, N
                  AAD(J,J-2) = ZERO
                  AAD(J,J-3) = ZERO
480            CONTINUE

!   ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N

               KA = I - 1
               B520 = .true.

               DO WHILE (B520.and.KA.lt.N1)

                  KA = KA + 1 
                  NOTLAS = KA.NE.N1
                  B490 = .true.

                  IF ( KA.EQ.I ) THEN
                     S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
                     IF ( LB.NE.I ) AAD(KA,KA-1) = - AAD(KA,KA-1)
                  ELSE
                     P = AAD(KA,KA-1)
                     Q = AAD(KA+1,KA-1)
                     R = ZERO
                     IF ( NOTLAS ) R = AAD(KA+2,KA-1)
                     X = DABS(P) + DABS(Q) + DABS(R)
                     IF ( X.EQ.ZERO ) THEN
                        B490 = .false.
                     ELSE
                        P = P / X
                        Q = Q / X
                        R = R / X
                        S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
                        AAD(KA,KA-1) = - S * X
                     END IF
                  ENDIF

!  Only do reminder if set

                  IF ( B490 ) then

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
490                  CONTINUE
!                                                 ** COLUMN MODIFICATION
                     DO 500 II = 1, MIN0(N,KA+3)
                        P = X * AAD(II,KA) + Y * AAD(II,KA+1)
                        IF ( NOTLAS ) THEN
                           P = P + Z * AAD(II,KA+2)
                           AAD(II,KA+2) = AAD(II,KA+2) - P * R
                        END IF
                        AAD(II,KA+1) = AAD(II,KA+1) - P * Q
                        AAD(II,KA)   = AAD(II,KA) - P
500                  CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
                     DO 510 II = L, K
                        P = X * EVECD(II,KA) + Y * EVECD(II,KA+1)
                        IF ( NOTLAS ) THEN
                           P = P + Z * EVECD(II,KA+2)
                           EVECD(II,KA+2) = EVECD(II,KA+2) - P * R
                        END IF
                        EVECD(II,KA+1) = EVECD(II,KA+1) - P * Q
                        EVECD(II,KA)   = EVECD(II,KA) - P
510                  CONTINUE

!  B490 clause

                  ENDIF

!  End B520 do while

               ENDDO

!  Clause for B390

            ENDIF

!  Finish loop for B390

         ENDDO

!  Finish Loop for B380

      ENDDO

!      do i = 1, 4
!       write(*,'(5F14.7)')(AAD(I,J),J=1,4), evald(i)
!      enddo
!      pause'after evals'

!  ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR

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
            DO 610  J = M, L, -1
               DO 600 I = L, K
                  Z = ZERO
                  DO 590 N = L, MIN0(J,K)
                     Z = Z + EVECD(I,N) * AAD(N,J)
590               CONTINUE
                  EVECD(I,J) = Z
600            CONTINUE
610         CONTINUE
         END IF

      END IF


!  Set egienvectors

      DO 625 I = L, K
         DO 620 J = 1, M
            EVECD(I,J) = EVECD(I,J) * WKD(I)
620      CONTINUE
625   CONTINUE



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
!  670 CONTINUE    ! Removed

!  Finish

      RETURN
END SUBROUTINE PCA_ASYMTX

!

subroutine PCA_Ranker(n,arrin,indx) bind(C)

! sorting the array indices in ascending order of array values
!   Numerical recipes routine, renamed

   implicit none

!  inputs

   integer(c_int)      :: n
   real(kind=c_double) :: arrin(n)

!  outputs

   integer(c_int)      :: indx(n)

!  local variables

   integer      :: i,j,l,ir,indxt
   real(kind=dp) :: q

      do 11 j = 1,n
        indx(j) = j
   11 continue
      if (n .eq. 1) return
        l = n/2+1
        ir = n
   10 continue
      if (l .gt. 1) then
        l = l-1
        indxt = indx(l)
        q = arrin(indxt)
      else
        indxt = indx(ir)
        q = arrin(indxt)
        indx(ir) = indx(1)
        ir = ir-1
        if (ir .eq. 1) then
          indx(1) = indxt
          return
        endif
      endif
      i = l
      j = l+l
   20 if (j .le. ir) then
        if (j .lt. ir) then
          if (arrin(indx(j)) .lt. arrin(indx(j+1))) j = j+1
        endif
        if (q .lt. arrin(indx(j))) then
          indx(i) = indx(j)
          i = j
          j = j+j
        else
          j = ir+1
        endif
        go to 20
      endif
      indx(i)=indxt
      go to 10

    return
end subroutine PCA_Ranker

!


  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER            ::  n
      REAL(kind=dp)       :: yp1,ypn,x(n),y(n),y2(n)
      INTEGER, PARAMETER :: NMAX=6094
      INTEGER            :: i,k
      REAL(kind=dp)       :: p,qn,sig,un,u(NMAX)

!  Check dimension

      if (n.gt.NMAX) then
        write(*,*)'Error in spline routine: too small NMAX =',NMAX
        stop
      endif

      if (yp1.gt..99e30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.0d0*( (y(i+1)-y(i)) / (x(i+1)-x(i)) - &
                     (y(i)-y(i-1)) / (x(i)-x(i-1))   &
                   ) / (x(i+1)-x(i-1)) - sig*u(i-1)  &
            )/p
      enddo

      if (ypn.gt..99d30) then
        qn=0.0d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo

      return
  END SUBROUTINE SPLINE


  SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER            ::  n
      REAL(kind=dp)       :: x, y, xa(n),ya(n),y2a(n)
      INTEGER            ::  k,khi,klo
      REAL(kind=dp)       :: a,b,h

      klo=1
      khi=n

 1    if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
!       if (h.eq.0.0d0) pause
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+ &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

      return
  END SUBROUTINE SPLINT


  SUBROUTINE PCA_LINTP2 (N1,X1,Y1,N2,X2,Y2)

      INTEGER        ::     N1,N2
      REAL(kind=dp)   :: X1(N1),Y1(N1),X2(N2),Y2(N2)
      INTEGER        :: INV, N2P, J, I2, I1, J2, I
      REAL(kind=dp)  :: X, YP, XP, S

      IF(N1.LE.0.OR.N2.LE.0.OR.N1.GT.32767.OR.N2.GT.32767)GOTO 9
!
! --- Checking the condition for inverting arrays ---
!
      INV=0
      IF (N1.EQ.1.OR.N2.EQ.1) GOTO 500
      IF ((X1(N1)-X1(1))*(X2(N2)-X2(1))) 300,300,500
!
! --- Inversion of new grid ---
!

      N2P=0
  300 CONTINUE
      INV=1
      N2P=N2/2
      DO 301 J=1,N2P
      XP=X2(J)
      X2(J)=X2(N2-J+1)
      X2(N2-J+1)=XP
  301 CONTINUE
!
! --- Main block ---
!
  500 IF (N1.EQ.1) GOTO 7
      S=DSIGN(1.0D0,X1(N1)-X1(1))
      I2=1
      I1=2
    1 IF((X2(I2)-X1(1))*S) 2,2,3
    2 Y2(I2)=Y1(1)
      I2=I2+1
      IF(I2.GT.N2) GOTO 999
      GOTO 1
    3 IF((X2(I2)-X1(I1))*S) 4,4,5
    4 X=(X2(I2)-X1(I1-1))/(X1(I1)-X1(I1-1))
      Y2(I2)=Y1(I1)*X+Y1(I1-1)*(1.0D0-X)
      I2=I2+1
      IF(I2.GT.N2) GOTO 999
      GOTO 3
    5 I1=I1+1
      IF(I1.LE.N1) GOTO 3
      DO 6 J2=I2,N2
    6 Y2(J2)=Y1(N1)
      GOTO 999
    7 DO 8 I=1,N2
    8 Y2(I)=Y1(1)
  999 CONTINUE
!
! --- Checking the condition for back inversion ---
!
      IF (INV.NE.1) GOTO 1000
!
! --- Back inversion ---
!
      DO 302 J=1,N2P
      XP=X2(J)
      YP=Y2(J)
      X2(J)=X2(N2-J+1)
      Y2(J)=Y2(N2-J+1)
      X2(N2-J+1)=XP
      Y2(N2-J+1)=YP
  302 CONTINUE
!
! --- Exit block ---
!
 1000 RETURN
    9 PRINT 100,N1,N2
  100 FORMAT(10X,'  Error in subroutine LINTP2 !',2I15)
      STOP '  Error in subroutine LINTP2 !'

  END SUBROUTINE PCA_LINTP2
    

SUBROUTINE PCA_DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER        ::    INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER        ::    IPIV( * )
      REAL(kind=dp)      ::   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL(kind=dp)      :: array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(kind=dp)      ::   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            :: I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      !EXTERNAL           PCA_DGEMM, PCA_DGETF2, PCA_DLASWP, PCA_DTRSM, PCA_XERBLA
!     ..
!     .. External Functions ..
!      INTEGER            :: PCA_ILAENV
      !EXTERNAL           PCA_ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PCA_XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Determine the block size for this environment.
!
      NB = PCA_ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL PCA_DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL PCA_DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL PCA_DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL PCA_DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                            IPIV, 1 )
!
!              Compute block row of U.
!
               CALL PCA_DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                           N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
                           LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL PCA_DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                              N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,    &
                              A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),  &
                              LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DGETRF
!
END SUBROUTINE PCA_DGETRF

SUBROUTINE PCA_DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER*1        :: TRANS
      INTEGER            :: INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(kind=dp)      ::   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) REAL(kind=dp)      :: array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) REAL(kind=dp)      :: array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(kind=dp)      ::   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL           ::  NOTRAN
!     ..
!     .. External Functions ..
!      LOGICAL            :: PCA_LSAME
      !EXTERNAL           PCA_LSAME
!     ..
!     .. External Subroutines ..
      !EXTERNAL           PCA_DLASWP, PCA_DTRSM, PCA_XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = PCA_LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.PCA_LSAME( TRANS, 'T' ) .AND. .NOT. &
          PCA_LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PCA_XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )  RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL PCA_DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL PCA_DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                     ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL PCA_DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                     NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL PCA_DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                     ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL PCA_DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
                     A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL PCA_DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of DGETRS
!
END SUBROUTINE PCA_DGETRS

!


SUBROUTINE PCA_DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER         :: INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(kind=dp)      ::   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) REAL(kind=dp)      :: array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER           :: I, IP, IX
!     ..
!     .. External Subroutines ..
      !EXTERNAL           PCA_DSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.EQ.0 ) RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I ) &
               CALL PCA_DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I ) &
               CALL PCA_DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I ) &
               CALL PCA_DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of DLASWP
!
END SUBROUTINE PCA_DLASWP


INTEGER FUNCTION PCA_ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    :: NAME, OPTS
      INTEGER            :: ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPE! for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPE!   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL

!  COmpiler dummies (R. Spurr)

      integer            nn3
      character*(1)      aa1
!     ..
!     .. Executable Statements ..
!

      nn3 = n3
      aa1 = opts(1:1)

!     Formerly a listed Goto: 
!       GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
      IF ( ISPEC.EQ.1 .or. ISPEC.eq.2 .or. ISPEC.eq.3 ) THEN

!  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      PCA_ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) )    &
                  SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )  &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )  RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
!     Formerly a listed Goto:  GO TO ( 110, 200, 300 ) ISPEC
!
!  110 CONTINUE
      IF ( ISPEC .EQ. 1 ) THEN
!
!     ISPE! = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      PCA_ILAENV = NB
      RETURN
!
!  200 CONTINUE
      ELSE IF ( ISPEC .EQ. 2 ) THEN
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.  &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      PCA_ILAENV = NBMIN
      RETURN
!
!  300 CONTINUE
      ELSE IF ( ISPEC .EQ. 3 ) THEN
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.  &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      PCA_ILAENV = NX
      RETURN
!
!     END clause ISPEC = 1, 2 or 3
!
      ENDIF
!
!  400 CONTINUE
      ELSE IF ( ISPEC .EQ. 4 ) THEN
!
!     ISPE! = 4:  number of shifts (used by xHSEQR)
!
      PCA_ILAENV = 6
      RETURN
!
!  500 CONTINUE
      ELSE IF ( ISPEC .EQ. 5 ) THEN
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      PCA_ILAENV = 2
      RETURN
!
!  600 CONTINUE 
      ELSE IF ( ISPEC .EQ. 6 ) THEN
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      PCA_ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
!  700 CONTINUE
      ELSE IF ( ISPEC .EQ. 7 ) THEN
!
!     ISPEC = 7:  number of processors (not used)
!
      PCA_ILAENV = 1
      RETURN
!
!  800 CONTINUE
      ELSE IF ( ISPEC .EQ. 8 ) THEN
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      PCA_ILAENV = 50
      RETURN

!
      ELSE
!
!     Invalid value of ISPEC
!
      PCA_ILAENV = -1
      RETURN
!
!     End of ILAENV
!
      ENDIF

!
!     End of ILAENV
!
END FUNCTION PCA_ILAENV

!

LOGICAL FUNCTION PCA_LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER         ::   INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      PCA_LSAME = CA.EQ.CB
      IF( PCA_LSAME )  RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDI! machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDI! is assumed - ZCODE is the EBCDI! code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.  &
             INTA.GE.145 .AND. INTA.LE.153 .OR.  &
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.  &
             INTB.GE.145 .AND. INTB.LE.153 .OR.  &
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      PCA_LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
END FUNCTION PCA_LSAME

SUBROUTINE PCA_XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        :: SRNAME
      INTEGER            :: INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
END SUBROUTINE PCA_XERBLA

!

SUBROUTINE PCA_DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(kind=dp)      ::   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL(kind=dp)      :: array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(kind=dp)      ::   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER         ::   J, JP
!     ..
!     .. External Functions ..
!      INTEGER         ::   PCA_IDAMAX
      !EXTERNAL             PCA_IDAMAX
!     ..
!     .. External Subroutines ..
      !EXTERNAL           PCA_DGER, PCA_DSCAL, PCA_DSWAP, PCA_XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PCA_XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + PCA_IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )  &
               CALL PCA_DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )  &
               CALL PCA_DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL PCA_DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,  &
                       A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of DGETF2
!
END SUBROUTINE PCA_DGETF2

subroutine PCA_dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      REAL(kind=dp)     ::  dx(*),dy(*)
      integer       :: i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return

      if(incx.eq.1.and.incy.eq.1)then

        m = mod(n,7)
        if( m .ne. 0 ) then
          do 30 i = 1,m
            dy(i) = dx(i)
   30     continue
          if( n .lt. 7 ) return
        endif

        mp1 = m + 1
        do 50 i = mp1,n,7
          dy(i) = dx(i)
          dy(i + 1) = dx(i + 1)
          dy(i + 2) = dx(i + 2)
          dy(i + 3) = dx(i + 3)
          dy(i + 4) = dx(i + 4)
          dy(i + 5) = dx(i + 5)
          dy(i + 6) = dx(i + 6)
   50   continue

      else
!
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do 10 i = 1,n
          dy(iy) = dx(ix)
          ix = ix + incx
          iy = iy + incy
   10   continue

      endif

!  Finish

      return
end subroutine PCA_dcopy

SUBROUTINE PCA_DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
                         BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1    ::   TRANSA, TRANSB
      INTEGER        ::   M, N, K, LDA, LDB, LDC
      REAL(kind=dp)      ::   ALPHA, BETA
!     .. Array Arguments ..
      REAL(kind=dp)      ::   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     ! := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and ! are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  ! an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(kind=dp)      :: array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL(kind=dp)      :: array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then ! need not be set on input.
!           Unchanged on exit.
!
!  !      - REAL(kind=dp)      :: array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  ! must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case ! need not be set on entry.
!           On exit, the array  !  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*! ).
!
!  LD!    - INTEGER.
!           On entry, LD! specifies the first dimension of ! as declared
!           in  the  calling  (sub)  program.   LD!  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!      LOGICAL         :: PCA_LSAME
      !EXTERNAL           PCA_LSAME
!     .. External Subroutines ..
      !EXTERNAL           PCA_XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL        ::   NOTA, NOTB
      INTEGER        ::   I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL(kind=dp)      ::   TEMP
!     .. Parameters ..
      REAL(kind=dp)      ::   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA  = PCA_LSAME( TRANSA, 'N' )
      NOTB  = PCA_LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.   &
               ( .NOT.PCA_LSAME( TRANSA, 'C' ) ).AND.   &
               ( .NOT.PCA_LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.   &
               ( .NOT.PCA_LSAME( TRANSB, 'C' ) ).AND.   &
               ( .NOT.PCA_LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL PCA_XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.   &
          ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )   &
         RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  ! := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  ! := alpha!A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  ! := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  ! := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
END SUBROUTINE PCA_DGEMM

SUBROUTINE PCA_DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,   &
                         BETA, Y, INCY )
!     .. Scalar Arguments ..
      REAL(kind=dp)      ::   ALPHA, BETA
      INTEGER        ::   INCX, INCY, LDA, M, N
      CHARACTER*1    ::   TRANS
!     .. Array Arguments ..
      REAL(kind=dp)      ::   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(kind=dp)      :: array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - REAL(kind=dp)      :: array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL(kind=dp)      :: array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(kind=dp)      ::   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
      REAL(kind=dp)      ::   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY

!     .. External Functions ..
!      LOGICAL        ::  PCA_LSAME
      !EXTERNAL           PCA_LSAME
!     .. External Subroutines ..
      !EXTERNAL           PCA_XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.PCA_LSAME( TRANS, 'N' ).AND.   &
               .NOT.PCA_LSAME( TRANS, 'T' ).AND.   &
               .NOT.PCA_LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL PCA_XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.   &
          ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )   &
         RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( PCA_LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO ) RETURN
      IF( PCA_LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
END SUBROUTINE PCA_DGEMV

!

SUBROUTINE PCA_DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      REAL(kind=dp)      ::   ALPHA
      INTEGER        ::   INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      REAL(kind=dp)      ::   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL(kind=dp)      :: array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL(kind=dp)      :: array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL(kind=dp)      :: array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(kind=dp)   , PARAMETER   :: ZERO = 0.0D+0
!     .. Local Scalars ..
      REAL(kind=dp)      ::  TEMP
      INTEGER        ::  I, INFO, IX, J, JY, KX

!     .. External Subroutines ..
!      !EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL PCA_XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )  RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of DGER  .
!
END SUBROUTINE PCA_DGER

!

subroutine PCA_dscal(n,da,dx,incx)

!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return

!        code for increment not equal to 1

      if(incx.ne.1) then
        nincx = n*incx
        do 10 i = 1,nincx,incx
          dx(i) = da*dx(i)
   10   continue
        return
      endif

!  clean-up loop

      m = mod(n,5)
      if( m .ne. 0 ) then
        do 30 i = 1,m
          dx(i) = da*dx(i)
   30   continue
        if( n .lt. 5 ) return
      endif

      mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return

end subroutine PCA_dscal

subroutine PCA_dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return

!       code for both increments equal to 1

      if(incx.eq.1.and.incy.eq.1) then

        m = mod(n,3)
        if ( m .ne. 0 ) then
          do 30 i = 1,m
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
   30     continue
          if( n .lt. 3 ) return
        endif

        mp1 = m + 1
        do 50 i = mp1,n,3
          dtemp = dx(i)
          dx(i) = dy(i)
          dy(i) = dtemp
          dtemp = dx(i + 1)
          dx(i + 1) = dy(i + 1)
          dy(i + 1) = dtemp
          dtemp = dx(i + 2)
          dx(i + 2) = dy(i + 2)
          dy(i + 2) = dtemp
   50   continue

      else
!
!       code for unequal increments or equal increments not equal to 1
!
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do 10 i = 1,n
          dtemp = dx(ix)
          dx(ix) = dy(iy)
          dy(iy) = dtemp
          ix = ix + incx
          iy = iy + incy
   10   continue
        return

      endif

      return
end subroutine PCA_dswap



SUBROUTINE PCA_DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )

!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL(kind=dp)      ::   ALPHA
!     .. Array Arguments ..
      REAL(kind=dp)      ::   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - REAL(kind=dp)      :: array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL(kind=dp)      :: array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      !LOGICAL            :: PCA_LSAME
      !EXTERNAL              PCA_LSAME
!     .. External Subroutines ..
      !EXTERNAL              PCA_XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC             MAX
!     .. Local Scalars ..
      LOGICAL            :: LSIDE, NOUNIT, UPPER
      INTEGER            :: I, INFO, J, K, NROWA
      REAL(kind=dp)          ::   TEMP
!     .. Parameters ..
      REAL(kind=dp)   , PARAMETER  ::  ONE  = 1.0D+0
      REAL(kind=dp)   , PARAMETER  ::  ZERO= 0.0D+0
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = PCA_LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = PCA_LSAME( DIAG  , 'N' )
      UPPER  = PCA_LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND. &
               ( .NOT.PCA_LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND. &
               ( .NOT.PCA_LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.PCA_LSAME( TRANSA, 'N' ) ).AND. &
               ( .NOT.PCA_LSAME( TRANSA, 'T' ) ).AND. &
               ( .NOT.PCA_LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.PCA_LSAME( DIAG  , 'U' ) ).AND. &
               ( .NOT.PCA_LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL PCA_XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( PCA_LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )  B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT ) B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT ) TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( PCA_LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!    End of DTRSM .
! 
END SUBROUTINE PCA_DTRSM

integer function PCA_idamax(n,dx,incx)
!
      implicit none

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision :: dx(*),dmax
      integer          :: i,incx,ix,n
!
      PCA_idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      PCA_idamax = 1
      if(n.eq.1)return

!   code for increment equal to 1

      if(incx.eq.1)then

        dmax = dabs(dx(1))
        do 30 i = 2,n
          if(dabs(dx(i)) .gt. dmax) then
            PCA_idamax = i
            dmax = dabs(dx(i))
          endif
   30   continue

!  code for increment not equal to 1

      else
        ix = 1
        dmax = dabs(dx(1))
        ix = ix + incx
        do 10 i = 2,n
           if(dabs(dx(ix)).gt.dmax) then
             PCA_idamax = i
             dmax = dabs(dx(ix))
           endif
          ix = ix + incx
   10   continue

      endif

      return
end function PCA_idamax

!  finish module

End Module pca_auxiliaries_m
