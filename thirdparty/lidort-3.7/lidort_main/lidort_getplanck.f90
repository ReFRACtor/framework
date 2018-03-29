
! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.7 F90                               #
! #  Release Date :   June 2014                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! #       Surface-leaving, BRDF Albedo-scaling     (3.7)    # 
! #       Taylor series, BBF Jacobians, ThreadSafe (3.7)    #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #          get_planckfunction                                 #
! #          get_planckfunction_plus                            #
! #                                                             #
! ###############################################################

      module LIDORT_getPlanck

      PRIVATE
      PUBLIC :: get_planckfunction, &
                get_planckfunction_plus

      contains

      subroutine get_planckfunction           &
           ( WNUMLO, WNUMHI, TEMPERATURE,     & ! Inputs
             BBFUNC, SMALLV, FAIL, MESSAGE )    ! Outputs

      USE LIDORT_PARS

      implicit none

!  Input arguments
!  ---------------

!     WNUMLO, WNUMHI are the wavenumber integration limits
!     TEMPERATURE is self-evident, must be in units K

      real(fpk),     intent(in)  :: WNUMLO, WNUMHI
      real(fpk),     intent(in)  :: TEMPERATURE

!  output arguments
!  ----------------

!     BBFUNC is the Planck function
!     SMALLV is a debug diagnostic

      real(fpk),     intent(out) :: BBFUNC
      INTEGER,       intent(out) :: SMALLV

!  Exception handling

      LOGICAL,       intent(out) :: FAIL
      CHARACTER*(*), intent(out) :: MESSAGE

!  Local variables
!  ---------------

!  Local parameters for:
!  (1) General use
      real(fpk), parameter :: &
        C2     = 1.438786_fpk,         & !Planck function constant #2
        SIGMA  = 5.67032e-8_fpk,       & !Stefan-Boltzmann constant
        SIGDPI = 1.80491891383e-8_fpk, & !SIGMA/Pi
        POWER4 = 4.0_fpk,              &
        CONC   = 1.5398973382e-01_fpk    !15.0d0/(Pi**POWER4)
!  (2) Method discrimination
      real(fpk), parameter :: &
        EPSIL  = 1.0e-8_fpk, &
        VMAX   = 32.0_fpk
!  (2a) Method #1: Simpson's rule
      integer, parameter :: &
        NSIMPSON = 25
      real(fpk), parameter :: &
        CRITERION = 1.0e-10_fpk
!  (2b) Method #2: Polynomial series
      real(fpk), parameter :: &
        A1 =  3.33333333333e-01_fpk, & !  1.0/3.0
        A2 = -1.25e-01_fpk,          & ! -1.0/8.0
        A3 =  1.66666666667e-02_fpk, & !  1.0/60.0
        A4 = -1.98412698413e-04_fpk, & ! -1.0/5040.0
        A5 =  3.67430922986e-06_fpk, & !  1.0/272160.0
        A6 = -7.51563251563e-08_fpk    ! -1.0/13305600.0
!  (2c) Method #3: Exponential series
      real(fpk), parameter :: &
        VCUT = 1.5_fpk
      real(fpk), parameter, dimension(7) :: &
        VCP  = (/ 10.25_fpk, 5.7_fpk, 3.9_fpk, 2.9_fpk, &
                   2.3_fpk,   1.9_fpk, 0.0_fpk /)

!  Other variables

      INTEGER   :: N, K, M, MMAX, I

      real(fpk) :: XLOW, XHIG, XX, X(2), XCUBE, EXPX, OEXPXM1
      real(fpk) :: RANGE, XSTEP, XSTEP3, FACTOR, GAMMA
      real(fpk) :: PLANCK_LOW, PLANCK_HIG, PLANCK, SCALING
      real(fpk) :: VAL, VAL0, OLDVAL, T

      real(fpk) :: EX, F, MV, M4, PL, XSQ
      real(fpk) :: PLANCK_EXP(2)
      real(fpk) :: PLANCK_POL(2)

!  Initialize output

      FAIL    = .false.
      MESSAGE = ' '
      SMALLV  = 0
      BBFUNC  = 0.0d0

!  Check input

      T = TEMPERATURE
      IF( T.LT.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0. ) THEN
        FAIL = .true.
        MESSAGE = 'Bad input--temperature or wavenums. wrong'
        RETURN
      END IF

!  Limits in x-space

      GAMMA   = C2 / T
      X(1) = GAMMA * WNUMLO
      X(2) = GAMMA * WNUMHI

!  Scaling constants

      SCALING  = SIGDPI * T ** POWER4

!  Wavenumbers are very close.  Get integral
!      by iterating Simpson rule to convergence.

      IF ( X(1).GT.EPSIL .AND. X(2).LT.VMAX .AND. &
         ( WNUMHI - WNUMLO ) / WNUMHI .LT. 1.D-2 ) THEN
!write(*,*) 'doing simpson'
         SMALLV = 3

!  interval

         XLOW  = X(1)
         XHIG  = X(2)
         RANGE = XHIG - XLOW

!  Two end values

         EXPX    = DEXP ( XLOW )
         OEXPXM1 = 1.0d0 / ( EXPX - 1.0d0 )
         XCUBE   = XLOW * XLOW * XLOW
         PLANCK_LOW  = XCUBE * OEXPXM1

         EXPX    = DEXP ( XHIG )
         OEXPXM1 = 1.0d0 / ( EXPX - 1.0d0 )
         XCUBE   = XHIG * XHIG * XHIG
         PLANCK_HIG  = XCUBE * OEXPXM1

!  Integral starting points

         VAL0   =  PLANCK_LOW +  PLANCK_HIG

!  First guess

         OLDVAL   = VAL0 * 0.5d0 * RANGE

!  Simpson's Rule up to 10 steps

         DO N = 1, NSIMPSON

!  Interval

          XSTEP  = 0.5d0 * RANGE / DBLE(N)
          XSTEP3 = XSTEP / 3.0d0

!  Integral

          VAL   = VAL0
          DO K = 1, 2*N - 1
             XX      = X(1) + DBLE(K) * XSTEP
             FACTOR  = DBLE(2*(1+MOD(K,2)))
             EXPX    = DEXP ( XX )
             OEXPXM1  = 1.0d0 / ( EXPX - 1.0d0 )
             XCUBE   = XX * XX * XX
             PLANCK  = XCUBE * OEXPXM1
             VAL  = VAL  + FACTOR * PLANCK
          ENDDO
          VAL  = VAL  * XSTEP3

!  Examine convergence

          IF (DABS((VAL-OLDVAL)/VAL).LE.CRITERION) GOTO  30
          OLDVAL = VAL

!  End integration loop

        ENDDO

!  No convergence, error message and return

        MESSAGE = 'Simpson rule didnt converge'
        FAIL    = .true.
        RETURN

!  Continuation point

   30   CONTINUE

!  Set the final answer

        BBFUNC  =  SCALING *  VAL * CONC

        RETURN

!  Finish

      ENDIF

!  Regular cases

      SMALLV = 0

!  Loop over two values

      DO I = 1, 2

!  Power series

       IF( X( I ).LT.VCUT ) THEN
!write(*,*) 'doing polynomial'
        SMALLV = SMALLV + 1
        XX  = X(I)
        XSQ = XX * XX
        PLANCK_POL(I) = CONC * XSQ * XX*( A1 + &
           XX*( A2 + XX*( A3 + XSQ*( A4 + XSQ*( A5 + XSQ*A6 ) ) ) ) )
       ELSE

!  Use exponential series
!write(*,*) 'doing exponential'
!  .......Find the upper limit of the series

        MMAX  = 0
   40   CONTINUE
        MMAX  = MMAX + 1
        IF( X(I) .LT. VCP( MMAX ) ) GO TO  40

!  .......Exponential series integration

        EX  = DEXP( - X(I) )
        F   = 1.0d0
        PL  = 0.0d0
        DO  M = 1, MMAX
          MV = M*X(I)
          F  = EX * F
          M4 = DBLE(M) ** (-POWER4)
          PL = PL + F*(6.0d0+MV*(6.0d0+MV*(3.0d0+MV)))*M4
        ENDDO
        PLANCK_EXP(I)  = PL * CONC
       ENDIF
      ENDDO

!   ** Handle ill-conditioning
!      SMALLV = 0   ---> ** WNUMLO and WNUMHI both small
!      SMALLV = 1   ---> ** WNUMLO small, WNUMHI large
!      SMALLV = 2   ---> ** WNUMLO and WNUMHI both large

      IF( SMALLV.EQ.2 ) THEN
        VAL  =  PLANCK_POL(2) -  PLANCK_POL(1)
      ELSE IF( SMALLV.EQ.1 ) THEN
        VAL  = 1.0d0 - PLANCK_POL (1) -  PLANCK_EXP(2)
      ELSE
        VAL  =  PLANCK_EXP(1) -  PLANCK_EXP(2)
      END IF

!  Set the final answer

      BBFUNC  =  SCALING *  VAL

!  Finish

      RETURN
      END SUBROUTINE get_planckfunction


      subroutine get_planckfunction_plus               &
           ( WNUMLO, WNUMHI, TEMPERATURE,              & ! Inputs
             BBFUNC, DBBFUNC, SMALLV, FAIL, MESSAGE )    ! Outputs

      USE LIDORT_PARS

      implicit none

!  Input arguments
!  ---------------

!     WNUMLO, WNUMHI are the wavenumber integration limits
!     TEMPERATURE is self-evident, must be in units K

      real(fpk),     intent(in)  :: WNUMLO, WNUMHI
      real(fpk),     intent(in)  :: TEMPERATURE

!  output arguments
!  ----------------

!     BBFUNC is the Planck function
!     DBBFUNC is the Planck function Temperature derivative (absolute)
!     SMALLV is a debug diagnostic

      real(fpk),     intent(out) :: BBFUNC
      real(fpk),     intent(out) :: DBBFUNC
      INTEGER,       intent(out) :: SMALLV

!  Exception handling

      LOGICAL,       intent(out) :: FAIL
      CHARACTER*(*), intent(out) :: MESSAGE

!  Local variables
!  ---------------

!  Local parameters for:
!  (1) General use
      real(fpk), parameter :: &
        C2     = 1.438786_fpk,         & !Planck function constant #2
        SIGMA  = 5.67032e-8_fpk,       & !Stefan-Boltzmann constant
        SIGDPI = 1.80491891383e-8_fpk, & !SIGMA/Pi
        POWER4 = 4.0_fpk,              &
        CONC   = 1.5398973382e-01_fpk    !15.0d0/(Pi**POWER4)
!  (2) Method discrimination
      real(fpk), parameter :: &
        EPSIL  = 1.0e-8_fpk, &
        VMAX   = 32.0_fpk
!  (2a) Method #1: Simpson's rule
      integer, parameter :: &
        NSIMPSON = 25
      real(fpk), parameter :: &
        CRITERION = 1.0e-10_fpk
!  (2b) Method #2: Polynomial series
      real(fpk), parameter :: &
        A1 =  3.33333333333e-01_fpk, & !  1.0/3.0
        A2 = -1.25e-01_fpk,          & ! -1.0/8.0
        A3 =  1.66666666667e-02_fpk, & !  1.0/60.0
        A4 = -1.98412698413e-04_fpk, & ! -1.0/5040.0
        A5 =  3.67430922986e-06_fpk, & !  1.0/272160.0
        A6 = -7.51563251563e-08_fpk    ! -1.0/13305600.0
      real(fpk), parameter :: &
        D1 =  3.33333333333e-01_fpk, & !  A1
        D2 =  0.0d0,                 & !  0.0d0
        D3 = -1.66666666667e-02_fpk, & ! -1.0d0 * A3
        D4 =  5.95238095238e-04_fpk, & ! -3.0d0 * A4
        D5 = -1.83715461493e-05_fpk, & ! -5.0d0 * A5
        D6 =  5.26094276094e-07_fpk    ! -7.0d0 * A6
!  (2c) Method #3: Exponential series
      real(fpk), parameter :: &
        VCUT = 1.5_fpk
      real(fpk), parameter, dimension(7) :: &
        VCP  = (/ 10.25_fpk, 5.7_fpk, 3.9_fpk, 2.9_fpk, &
                   2.3_fpk,   1.9_fpk, 0.0_fpk /)

!  Other variables

      INTEGER   :: N, K, M, MMAX, I
      real(fpk) :: XLOW, XHIG, XX, X(2), XCUBE, EXPX, OEXPXM1
      real(fpk) :: RANGE, XSTEP, XSTEP3, FACTOR, GAMMA
      real(fpk) :: PLANCK_LOW, PLANCK_HIG, PLANCK, SCALING, T
      real(fpk) :: DPLANCK_LOW, DPLANCK_HIG, DPLANCK, DSCALING
      real(fpk) :: VAL, VAL0, OLDVAL, DVAL, DVAL0, DOLDVAL

      real(fpk) :: EX, F, MV, M4, PL, DL, XSQ
      real(fpk) :: PLANCK_EXP(2), DPLANCK_EXP(2)
      real(fpk) :: PLANCK_POL(2), DPLANCK_POL(2)

!  Initialize output

      FAIL    = .false.
      MESSAGE = ' '
      SMALLV  = 0
      BBFUNC  = 0.0d0
      DBBFUNC = 0.0d0

!  Check input

      T = TEMPERATURE
      IF( T.LT.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0. ) THEN
        FAIL = .true.
        MESSAGE = 'Bad input--temperature or wavenums. wrong'
        RETURN
      END IF

!  Limits in x-space

      GAMMA   = C2 / T
      X(1) = GAMMA * WNUMLO
      X(2) = GAMMA * WNUMHI

!  Scaling constants

      SCALING  = SIGDPI * T ** POWER4
      DSCALING = SCALING / T 

!  Wavenumbers are very close.  Get integral
!      by iterating Simpson rule to convergence.

      IF ( X(1).GT.EPSIL .AND. X(2).LT.VMAX .AND. &
          ( WNUMHI - WNUMLO ) / WNUMHI .LT. 1.D-2 ) THEN

         SMALLV = 3

!  interval

         XLOW  = X(1)
         XHIG  = X(2)
         RANGE = XHIG - XLOW

!  Two end values

         EXPX    = DEXP ( XLOW )
         OEXPXM1 = 1.0d0 / ( EXPX - 1.0d0 )
         XCUBE   = XLOW * XLOW * XLOW
         PLANCK_LOW  = XCUBE * OEXPXM1
         DPLANCK_LOW = PLANCK_LOW * XLOW * EXPX * OEXPXM1

         EXPX    = DEXP ( XHIG )
         OEXPXM1 = 1.0d0 / ( EXPX - 1.0d0 )
         XCUBE   = XHIG * XHIG * XHIG
         PLANCK_HIG  = XCUBE * OEXPXM1
         DPLANCK_HIG = PLANCK_HIG * XHIG * EXPX * OEXPXM1

!  Integral starting points

         VAL0   =  PLANCK_LOW +  PLANCK_HIG
         DVAL0  = DPLANCK_LOW + DPLANCK_HIG

!  First guess

         OLDVAL   = VAL0 * 0.5d0 * RANGE
         DOLDVAL  = DVAL0 * 0.5d0 * RANGE

!  Simpson's Rule up to 10 steps

         DO N = 1, NSIMPSON

!  Interval

          XSTEP  = 0.5d0 * RANGE / DBLE(N)
          XSTEP3 = XSTEP / 3.0d0

!  Integral

          VAL   = VAL0
          DVAL  = DVAL0
          DO K = 1, 2*N - 1
             XX      = X(1) + DBLE(K) * XSTEP
             FACTOR  = DBLE(2*(1+MOD(K,2)))
             EXPX    = DEXP ( XX )
             OEXPXM1  = 1.0d0 / ( EXPX - 1.0d0 )
             XCUBE   = XX * XX * XX
             PLANCK  = XCUBE * OEXPXM1
             DPLANCK = PLANCK * XX * EXPX * OEXPXM1
             VAL  = VAL  + FACTOR * PLANCK
             DVAL = DVAL + FACTOR * DPLANCK
          ENDDO
          VAL  = VAL  * XSTEP3
          DVAL = DVAL * XSTEP3

!  Examine convergence

          IF (DABS((VAL-OLDVAL)/VAL).LE.CRITERION) GOTO  30
          OLDVAL = VAL
          DOLDVAL = DVAL

!  End integration loop

        ENDDO

!  No convergence, error message and return

        MESSAGE = 'Simpson rule didnt converge'
        FAIL    = .true.
        RETURN

!  Continuation point

   30   CONTINUE

!  Set the final answer

        BBFUNC  =  SCALING *  VAL * CONC
        DBBFUNC = DSCALING * DVAL * CONC

        RETURN

!  Finish

      ENDIF

!  Regular cases

      SMALLV = 0

!  Loop over two values

      DO I = 1, 2

!  Power series

       IF( X( I ).LT.VCUT ) THEN
        SMALLV = SMALLV + 1
        XX  = X(I)
        XSQ = XX * XX
        PLANCK_POL(I) = CONC * XSQ * XX*( A1 +  &
           XX*( A2 + XX*( A3 + XSQ*( A4 + XSQ*( A5 + XSQ*A6 ) ) ) ) )
        DPLANCK_POL(I) = CONC * XSQ * XX*( D1 + &
           XX*( D2 + XX*( D3 + XSQ*( D4 + XSQ*( D5 + XSQ*D6 ) ) ) ) )
       ELSE

!  Use exponential series

!  .......Find the upper limit of the series

        MMAX  = 0
   40   CONTINUE
        MMAX  = MMAX + 1
        IF( X(I) .LT. VCP( MMAX ) ) GO TO  40

!  .......Exponential series integration

        EX  = DEXP( - X(I) )
        F   = 1.0d0
        PL  = 0.0d0
        DL  = 0.0D0
        DO  M = 1, MMAX
          MV = M*X(I)
          F  = EX * F
          M4 = DBLE(M) ** (-POWER4)
          PL = PL + F*(6.0d0+MV*(6.0d0+MV*(3.0d0+MV)))*M4
          DL = DL + F*(24.0d0+MV*(24.0d0+MV*(12.0d0+MV*(4.0d0+MV))))*M4
        ENDDO
        PLANCK_EXP(I)  = PL * CONC
        DPLANCK_EXP(I) = DL * CONC
       ENDIF
      ENDDO

!   ** Handle ill-conditioning
!      SMALLV = 0   ---> ** WNUMLO and WNUMHI both small
!      SMALLV = 1   ---> ** WNUMLO small, WNUMHI large. Look out for constants
!      SMALLV = 2   ---> ** WNUMLO and WNUMHI both large

      IF( SMALLV.EQ.2 ) THEN
        VAL  =  PLANCK_POL(2) -  PLANCK_POL(1)
        DVAL = DPLANCK_POL(2) - DPLANCK_POL(1)
      ELSE IF( SMALLV.EQ.1 ) THEN
        VAL  = 1.0d0 - PLANCK_POL (1) -  PLANCK_EXP(2)
        DVAL = 4.0d0 - DPLANCK_POL(1) - DPLANCK_EXP(2)
      ELSE
        VAL  =  PLANCK_EXP(1) -  PLANCK_EXP(2)
        DVAL = DPLANCK_EXP(1) - DPLANCK_EXP(2)
      END IF

!  Set the final answer

      BBFUNC  =  SCALING *  VAL
      DBBFUNC = DSCALING * DVAL

!  Finish

      RETURN
      end subroutine get_planckfunction_plus

      end module LIDORT_getPlanck
