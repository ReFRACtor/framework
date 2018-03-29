! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #        --           -            -        -        -    #
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
! #       NEW: THERMAL SUPPLEMENT INCLUDED     (3.2)        #
! #       NEW: OUTGOING SPHERICITY CORRECTION  (3.2)        #
! #       NEW: TOTAL COLUMN JACOBIANS          (3.3)        #
! #       VLIDORT COMPATIBILITY                (3.4)        #
! #       THREADED/OPTIMIZED F90 code          (3.5)        #
! #       EXTERNAL SS / NEW I/O STRUCTURES     (3.6)        #
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
! #            LISPARSE_FUNCTION_PLUS                           #
! #            LIDENSE_FUNCTION_PLUS                            #
! #            HAPKE_FUNCTION_PLUS                              #
! #            RAHMAN_FUNCTION_PLUS                             #
! #            COXMUNK_FUNCTION_PLUS                            #
! #            COXMUNK_FUNCTION_MSR_PLUS                        #
! #                                                             #
! # New BPDF Subroutines in this Module (Version 3.7)           #
! #                                                             #
! #            BPDFVEGN_FUNCTION_PLUS                           #
! #            BPDFSOIL_FUNCTION_PLUS                           #
! #            BPDFNDVI_FUNCTION_PLUS                           #
! #            FRESNEL_SCALAR_PLUS (Private, called by BPDF)    #
! #                                                             #
! # New Cox-Munk Subroutines in this Module (Version 3.7)       #
! #                                                             #
! #            BRDF_Generalized_Glint_plus                      #
! #            BRDF_WhiteCap_Reflectance_plus                   #
! #                                                             #
! ###############################################################

      MODULE brdf_LinSup_kernels_m

!  Rob Extension 12/2/14. BPDF Kernels (replace BREONVEG, BREONSOIL)

      use brdf_sup_aux_m, only : derfc_e, BRDF_Fresnel_Complex

      PRIVATE
      PUBLIC :: LISPARSE_FUNCTION_PLUS, &
                LIDENSE_FUNCTION_PLUS, &
                HAPKE_FUNCTION_PLUS, &
                RAHMAN_FUNCTION_PLUS, &
                COXMUNK_FUNCTION_PLUS, &
                COXMUNK_FUNCTION_MSR_PLUS,    &
                BPDFVEGN_FUNCTION_PLUS,       &
                BPDFSOIL_FUNCTION_PLUS,       &
                BPDFNDVI_FUNCTION_PLUS,       &
                BRDF_Generalized_Glint_plus,  &
                BRDF_WhiteCap_Reflectance_plus

      CONTAINS

      SUBROUTINE LISPARSE_FUNCTION_PLUS              &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,   &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,    &
              LISPARSE_KERNEL, LISPARSE_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, HALF, PIE, TWO

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LISPARSE_KERNEL
      REAL(fpk), intent(out) :: LISPARSE_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(fpk)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(fpk)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(fpk)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(fpk)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(fpk)  :: DX_H, DX_R, DX_Q, DX_P, DX_QR, DY_Q, TERM2       ! Variables A2, R2 removed 10/12/12, TERM2 added
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LISPARSE_KERNEL      = ZERO
!mick fix 8/8/2014 - initialize all elements
      LISPARSE_DERIVATIVES = ZERO
      !DO J = 1, NPARS
      !  LISPARSE_DERIVATIVES(J) = ZERO
      !ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Rob Fix 9/25/14. Hot spot condition  was wrongly inserted here
!      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      TERM2 = SX_INC * SX_REF * CKPHI
      CKSI  = X_INC  * X_REF + TERM2

!  contributions P and R and derivatives (if flagged)

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version
!      Linearization of P has changed in the new version

!  Old   P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
! Old       R2   = R * R
! Old       A2   = A * A
! Old       DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
        DX_P = ( ( SX_INC * SX_INC + SX_REF * SX_REF ) * ( P - ONE ) + &
                  TERM2 * ( X_INC * B + X_REF * A ) ) / PARS(2)
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_KERNEL = HALF * P - QR

!  Set derivatives
!  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LISPARSE_DERIVATIVES(1) = - R * DY_Q
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_QR = DX_R * Q + DX_Q * R
        LISPARSE_DERIVATIVES(2) = HALF * DX_P - DX_QR
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LISPARSE_FUNCTION_PLUS

!

      SUBROUTINE LIDENSE_FUNCTION_PLUS               &
             ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
               XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
               LIDENSE_KERNEL, LIDENSE_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, PIE, TWO

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LIDENSE_KERNEL
      REAL(fpk), intent(out) :: LIDENSE_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(fpk)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(fpk)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(fpk)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(fpk)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      REAL(fpk)  :: DX_H, DX_R, DX_Q, DX_P, DX_P_QR, DY_Q, TERM2       ! Variables A2, R2 removed 10/12/12, TERM2 added
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LIDENSE_KERNEL      = ZERO
!mick fix 8/8/2014 - initialize all elements
      LIDENSE_DERIVATIVES = ZERO
      !DO J = 1, NPARS
      !  LIDENSE_DERIVATIVES(J) = ZERO
      !ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Rob Fix 9/25/14. Hot spot condition  was wrongly inserted here
!      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      TERM2 = SX_INC * SX_REF * CKPHI
      CKSI  = X_INC  * X_REF + TERM2

!  contributions P and R and derivatives (if flagged)

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version
!      Linearization of P has changed in the new version

!  Old   P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
! Old       R2   = R * R
! Old       A2   = A * A
! Old       DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
        DX_P = ( ( SX_INC * SX_INC + SX_REF * SX_REF ) * ( P - ONE ) + &
                  TERM2 * ( X_INC * B + X_REF * A ) ) / PARS(2)
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - TWO

!  Set derivatives
!  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LIDENSE_DERIVATIVES(1) = - P_QR * DY_Q / Q 
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_P_QR = ( DX_P / P ) - ( DX_R / R ) - ( DX_Q / Q )
        LIDENSE_DERIVATIVES(2) = P_QR * DX_P_QR
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDENSE_FUNCTION_PLUS

!

      SUBROUTINE HAPKE_FUNCTION_PLUS                 &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,   &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,    &
              HAPKE_KERNEL, HAPKE_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, TWO, HALF, QUARTER, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: HAPKE_KERNEL
      REAL(fpk), intent(out) :: HAPKE_DERIVATIVES ( MAXPARS )

!  Hapke Kernel function.
!    - New version, Fresh Coding
!    - Old version uses DISORT code; for validation.

!  input variables:

!    XI, SXI  : Cosine/Sine of angle of reflection (positive)
!    XJ, SXJ  : Cosine/Sine of angle of incidence (positive)
!    XPHI     : Difference of azimuth angles of incidence and reflection
!    PARS(1)  : single scattering albedo in Hapke's BDR model
!    PARS(2)  : angular width parameter of opposition effect in Hapke's model
!    PARS(3)  : Empirical hot spot multiplier

!  local variables
!    B0_EMPIR : empirical factor to account for the finite size of
!               particles in Hapke's BDR model
!    B_HOT    : term that accounts for the opposition effect
!               (retroreflectance, hot spot) in Hapke's BDR model
!    CTHETA   : cosine of phase angle in Hapke's BDR model
!    GAMMA    : albedo factor in Hapke's BDR model
!    PHASE    : scattering phase function in Hapke's BDR model
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model

!  Local variables

      INTEGER    :: J
      REAL(fpk)  :: CTHETA, THETA, PHASE
      REAL(fpk)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(fpk)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(fpk)  :: HELP_J, GHELP_J, TERM_J
      REAL(fpk)  :: HELP_I, GHELP_I, TERM_I
      REAL(fpk)  :: TI_TJ, DT1, DT2
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      HAPKE_KERNEL      = ZERO
!mick fix 8/8/2014 - initialize all elements
      HAPKE_DERIVATIVES = ZERO
      !DO J = 1, NPARS
      !  HAPKE_DERIVATIVES(J) = ZERO
      !ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!        CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      THETA  = DACOS( CTHETA )
      PHASE  = ONE + HALF * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( HALF * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( ONE - SSALBEDO )
      HELP_J   = TWO * XJ
      GHELP_J  = ( ONE + HELP_J * GAMMA )
      TERM_J   = ( ONE + HELP_J ) / GHELP_J
      HELP_I   = TWO * XI
      GHELP_I  = ( ONE + HELP_I * GAMMA )
      TERM_I   = ( ONE + HELP_I ) / GHELP_I
      TI_TJ    = TERM_J * TERM_I

!  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION     = ( ONE + B_HOT ) * PHASE + TI_TJ - ONE
      HAPKE_KERNEL = REFLEC * FUNCTION

!  ssalbedo derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        DT1 = HAPKE_KERNEL / SSALBEDO
        DT2 = ( HELP_J / GHELP_J ) + ( HELP_I / GHELP_I )
        DT2 = DT2 * TI_TJ * HALF / GAMMA
        HAPKE_DERIVATIVES(1) = DT1 + DT2 * REFLEC
      ENDIF

!  Hotspot  derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        DT1 = B_HOT * ( B0_EMPIR - B_HOT ) / B0_EMPIR / HOTSPOT
        HAPKE_DERIVATIVES(2) = DT1 * REFLEC * PHASE
      ENDIF

!  empirical factor derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        DT1 = B_HOT / B0_EMPIR 
        HAPKE_DERIVATIVES(3) = DT1 * REFLEC * PHASE
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE HAPKE_FUNCTION_PLUS

!

      SUBROUTINE RAHMAN_FUNCTION_PLUS              &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,  &
              RAHMAN_KERNEL, RAHMAN_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, TWO, ONEP5, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: RAHMAN_KERNEL
      REAL(fpk), intent(out) :: RAHMAN_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(fpk)  :: T_INC, T_REF, DT1, DT2, D_FACT, D_HELPM
      REAL(fpk)  :: CXI, DELTA, K1_SQ, FACT, D_K0, D_K1, D_K2
      REAL(fpk)  :: GEOM, PHASE, RFAC, RFAC1, K0, K1, K2, HELPR
      REAL(fpk)  :: XPHI, CKPHI, HSPOT, UPPER_LIMIT, HELPG, HELPM
      REAL(fpk), PARAMETER  :: SMALL = 1.0d-04

!  Initialise

      RAHMAN_KERNEL      = ZERO
!mick fix 8/8/2014 - initialize all elements
      RAHMAN_DERIVATIVES = ZERO
      !DO J = 1, NPARS
      !  RAHMAN_DERIVATIVES(J) = ZERO
      !ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( 2.0d0 - K0 )
      FACT = FACT * ( 1.0d0 - K1 ) / ( 1.0d0 + K1 ) / ( 1.0d0 + K1 )
      GEOM = ( 2.0d0 * XJ * XJ * XJ ) ** ( K2 - 1.0d0 )
      HSPOT = FACT * GEOM

!  Upper limit ( 5 times hotspot value ). Follwing comments inserted.
!     This function needs more checking; some constraints are 
!     required to avoid albedos larger than 1; in particular,
!     the BDREF is limited to 5 times the hotspot value to
!     avoid extremely large values at low polar angles

      UPPER_LIMIT = 5.0d0 * HSPOT

!  hot spot value

      IF ( DABS(PHI) .LT. SMALL .AND. XI.EQ.XJ ) THEN
        RAHMAN_KERNEL = HSPOT
        RETURN
      ENDIF

!  Use upper limit value at edges (low incidence or reflection)

      IF ( XI.LT.SMALL .OR. XJ.LT.SMALL ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
        RETURN
      ENDIF

!  Main section
!  ------------

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. ONE ) CXI = ONE

!  Phase function

      K1_SQ = K1 * K1
      HELPM = ONE - K1_SQ 
      FACT  = ONE + K1_SQ + TWO * K1 * CXI
      PHASE = HELPM / ( FACT ** ONEP5 )

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      HELPR = ONE / ( ONE + DELTA )
      RFAC  = ( ONE - K0 ) * HELPR
      RFAC1 = ONE + RFAC

!  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * RFAC1 * GEOM

!  K0 derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_K0   = ( ONE / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_DERIVATIVES(1) = RAHMAN_KERNEL * D_K0
      ENDIF

!  Phase function derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        D_FACT  =   TWO * K1 + TWO * CXI
        D_HELPM = - TWO * K1
        D_K1    = ( D_HELPM / HELPM ) - ONEP5 * ( D_FACT / FACT )
        RAHMAN_DERIVATIVES(2) = RAHMAN_KERNEL * D_K1
      ENDIF

!  K2 derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_DERIVATIVES(3) = RAHMAN_KERNEL * D_K2
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAHMAN_FUNCTION_PLUS

!

      SUBROUTINE COXMUNK_FUNCTION_PLUS              &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
              COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, MINUS_ONE, TWO, FOUR, HALF, QUARTER, PIE, PIO2
      USE BRDF_SUP_AUX_M

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: COXMUNK_KERNEL
      REAL(fpk), intent(inout) :: COXMUNK_DERIVATIVES ( MAXPARS )

!  Critical exponent taken out

      REAL(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      INTEGER    :: J
      REAL(fpk)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(fpk)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(fpk)  :: XPHI, CKPHI, T1_R, T2_R, DCOT_R, T1_I, T2_I, DCOT_I
      REAL(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(fpk)  :: SHADOWI, SHADOWR, SHADOW
      REAL(fpk)  :: H1H2, H2Z2, TA_SQ, DFAC2, DH1, DH2, DRP, DRL
      REAL(fpk)  :: D_S1, D_S2, D_T1, D_T2
      REAL(fpk)  :: D_SHADOWI, D_SHADOWR, D_SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL      = ZERO
!mick fix 8/8/2014 - initialize all elements
      COXMUNK_DERIVATIVES = ZERO
      !DO J = 1, NPARS
      !  COXMUNK_DERIVATIVES(J) = ZERO
      !ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI   
! old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI   
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      H1H2 = H1 + H2
      RP = ( H1 - H2 ) / H1H2
      H2Z2 = Z2 + H2
      RL = ( Z2 - H2 ) / H2Z2
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A  = PIO2 - DASIN(B)
      TA = DTAN(A)
      TA_SQ = TA * TA
      ARGUMENT = TA_SQ  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DFAC2 = ( PARS(1) - TA_SQ ) / PARS(1) / PARS(1)
          COXMUNK_DERIVATIVES(1) = - COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  square refractive index derivative
!  --This section of code was formerly at the end of routine
!    -- Now moved here before the shadowing option
!    -- otherwise derivative will not get done
!         Bug found by V. Natraj. 02 February 2007.

      IF ( DO_DERIV_PARS(2) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DH1 = Z2
          DH2 = HALF / H2
          DRP = ( DH1 * ( ONE - RP ) - DH2 * ( ONE + RP ) ) / H1H2
          DRL =  - DH2 * ( ONE + RL ) / H2Z2
          DFAC2 = ( RP*DRP + RL*DRL ) / XMP
          COXMUNK_DERIVATIVES(2) = COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code
!  -----------

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT_I = XI/DSQRT(ONE-XXI)
      T1_I   = DEXP(-DCOT_I*DCOT_I*S2)
!      T2_I   = DERFC(DCOT_I*S3)
      T2_I   = DERFC_E(DCOT_I*S3)
      SHADOWI = HALF * ( S1*T1_I/DCOT_I - T2_I )

      XXJ  = XJ*XJ
      DCOT_R = XJ/DSQRT(ONE-XXJ)
      T1_R   = DEXP(-DCOT_R*DCOT_R*S2)
!      T2_R   = DERFC(DCOT_R*S3)
      T2_R   = DERFC_E(DCOT_R*S3)
      SHADOWR = HALF * ( S1*T1_R/DCOT_R - T2_R )

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

!  Update Scalar derivatives
!  -------------------------

!  add the shadow derivative to inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_S1 = HALF / PIE / S1
        D_S2 = - S2 * S2
        D_T1 = - T1_I * DCOT_I * DCOT_I * D_S2
        D_T2 = S2 * S2 * DCOT_I * S1 * T1_I 
        D_SHADOWI = HALF * ( D_S1*T1_I/DCOT_I + S1*D_T1/DCOT_I - D_T2 )
        D_T1 = - T1_R * DCOT_R * DCOT_R * D_S2
        D_T2 = S2 * S2 * DCOT_R * S1 * T1_R
        D_SHADOWR = HALF * ( D_S1*T1_R/DCOT_R + S1*D_T1/DCOT_R - D_T2 )
        D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
        COXMUNK_DERIVATIVES(1) = COXMUNK_DERIVATIVES(1) * SHADOW + &
                                 COXMUNK_KERNEL * D_SHADOW / SHADOW
      ENDIF

!  Refractive index derivative, update

      IF ( DO_DERIV_PARS(2) ) THEN
        COXMUNK_DERIVATIVES(2) = COXMUNK_DERIVATIVES(2) * SHADOW 
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION_PLUS

!

      SUBROUTINE COXMUNK_FUNCTION_MSR_PLUS               &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,       &
              ORDER, N_MUQUAD, N_PHIQUAD,                &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,        &
              X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, &
              X_PHIQUAD, W_PHIQUAD,                      &
              COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, PIE, DEG_TO_RAD, &
                              MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD
      USE BRDF_SUP_AUX_M

      implicit none

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      INTEGER  , intent(in)  :: ORDER
      INTEGER  , intent(in)  :: N_MUQUAD, N_PHIQUAD
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(in)  :: X_MUQUAD   ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: W_MUQUAD   ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: SX_MUQUAD  ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: WXX_MUQUAD ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: X_PHIQUAD  ( MAX_MSRS_PHIQUAD )
      REAL(fpk), intent(in)  :: W_PHIQUAD  ( MAX_MSRS_PHIQUAD )

      REAL(fpk), intent(out) :: COXMUNK_KERNEL
      REAL(fpk), intent(out) :: COXMUNK_DERIVATIVES ( MAXPARS )

!  local variables
!  ---------------

!  help variables

      integer    :: s, n, k, i, i1, N_phiquad_HALF, ni, ki, nr, kr, q
      REAL(fpk)  :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr, sumr, sum, w_p
      REAL(fpk)  :: reflec_0, reflec_s, R0Q, D_R0Q, reflec
      REAL(fpk)  :: phi_sub1, cphi_sub1, sphi_sub1
      REAL(fpk)  :: phi_sub2, cphi_sub2, sphi_sub2
      REAL(fpk)  :: d_reflec_0(3), d_reflec_s(3), d_reflec(3)

!  arrays

      REAL(fpk)  :: R0_QUAD_IN   ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: R0_OUT_QUAD  ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: D_R0_QUAD_IN ( MAXPARS, MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: D_R0_OUT_QUAD( MAXPARS, MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )

      REAL(fpk)  :: RHOLD   ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: D_RHOLD ( MAXPARS, MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )

      REAL(fpk)  :: R0_MSRS_QUAD   ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: D_R0_MSRS_QUAD ( MAXPARS, MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )

!  Safety first zeroing

      REFLEC_0 = ZERO
      COXMUNK_KERNEL = ZERO

      D_REFLEC_0 = ZERO
      COXMUNK_DERIVATIVES = ZERO

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_FUNCTION_PLUS                   &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,  &
              REFLEC_0, D_REFLEC_0 )

!  Higher orders scattering

      REFLEC   = REFLEC_0
      REFLEC_S = ZERO

      D_REFLEC   = D_REFLEC_0
      D_REFLEC_S = ZERO

!  Quadrature output for first order R/T calculations
!        This will be overwritten as the orders increase

      IF ( ORDER.GE.1 ) THEN
        DO K = 1, N_MUQUAD
          XM  = X_MUQUAD(K)
          SXM = SX_MUQUAD(K)
          DO N = 1, N_PHIQUAD
            PHI_SUB1  = X_PHIQUAD(N)
            CPHI_SUB1 = DCOS(PHI_SUB1)
            SPHI_SUB1 = DSIN(PHI_SUB1)
            PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
            CPHI_SUB2 = DCOS(PHI_SUB2)
            SPHI_SUB2 = DSIN(PHI_SUB2)
            CALL COXMUNK_FUNCTION_PLUS                              &
               ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,               &
                 XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,  &
                 R0_OUT_QUAD(K,N), D_R0_OUT_QUAD(:,K,N) )
            CALL COXMUNK_FUNCTION_PLUS                              &
               ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,               &
                 XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,  &
                 R0_QUAD_IN(K,N),  D_R0_QUAD_IN(:,K,N) )
          ENDDO
        ENDDO
      ENDIF

!  Compute the successive orders of scattering.
!  Compute higher orders.

      DO S = 1, ORDER

!  Compute result for this order

         SUMR = ZERO
         DO K = 1, N_MUQUAD
            SUM_PR = ZERO
            DO N = 1, N_PHIQUAD
               W_P = W_PHIQUAD(N)
               SUM = R0_QUAD_IN(K,N) * R0_OUT_QUAD(K,N)
               SUM_PR = SUM_PR + W_P * SUM
            ENDDO
            SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
         ENDDO
         REFLEC_S = SUMR

!  Derivatives

         DO Q = 1, 2
            IF ( DO_DERIV_PARS(Q) ) THEN
               SUMR = ZERO
               DO K = 1, N_MUQUAD
                  SUM_PR = ZERO
                  DO N = 1, N_PHIQUAD
                     W_P  = W_PHIQUAD(N)
                     SUM_PR = SUM_PR + W_P *   R0_QUAD_IN(K,N)   * D_R0_OUT_QUAD(Q,K,N) &
                                     + W_P * D_R0_QUAD_IN(Q,K,N) *   R0_OUT_QUAD(K,N)
                  ENDDO
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
               ENDDO
               D_REFLEC_S(Q) = SUMR
            ENDIF
         ENDDO

!  Finish if reached the scattering order desired

         IF ( S.EQ.ORDER ) GO TO 67

!  Compute Reflectance for next order and update
!    Quad-Quad results get computed each time: very wasteful.
!    Have to do this, as the memory is a killer.

         DO KR = 1, N_MUQUAD
            XMR  = X_MUQUAD(KR)
            SXMR = SX_MUQUAD(KR)
            DO NR = 1, N_PHIQUAD

!  Quad-quad calculations

               DO KI = 1, N_MUQUAD
                  XMI  = X_MUQUAD(KI)
                  SXMI = SX_MUQUAD(KI)
                  DO NI = 1, N_PHIQUAD
                     PHI_SUB1  = X_PHIQUAD(NR) - X_PHIQUAD(NI)
                     CPHI_SUB1 = DCOS(PHI_SUB1)
                     SPHI_SUB1 = DSIN(PHI_SUB1)
                     CALL COXMUNK_FUNCTION_PLUS &
                       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1, &
                         R0_MSRS_QUAD(KI,NI), D_R0_MSRS_QUAD(:,KI,NI) )
                  ENDDO
               ENDDO

!  Multiple reflection

               SUMR = ZERO
               DO KI = 1, N_MUQUAD
                  SUM_PR = ZERO
                  DO NI = 1, N_PHIQUAD
                     W_P  = W_PHIQUAD(NI)
                     R0Q = R0_MSRS_QUAD(KI,NI)
                     SUM_PR = SUM_PR + W_P * R0_QUAD_IN(KI,NI) * R0Q
                  ENDDO
!                  SUMR = SUMR + SUM_PR * W_MUQUAD(KI)
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(KI)
               ENDDO
               RHOLD(KR,NR) = SUMR

!  Derivatives

               DO Q = 1, 2
                  IF ( DO_DERIV_PARS(Q) ) THEN
                     SUMR = ZERO
                     DO KI = 1, N_MUQUAD
                        SUM_PR = ZERO
                        DO NI = 1, N_PHIQUAD
                           W_P = W_PHIQUAD(NI)
                           R0Q   = R0_MSRS_QUAD(KI,NI)
                           D_R0Q = D_R0_MSRS_QUAD(Q,KI,NI)
                           SUM_PR = SUM_PR + W_P *   R0_QUAD_IN(KI,NI)   * D_R0Q &
                                           + W_P * D_R0_QUAD_IN(Q,KI,NI) *   R0Q
                        ENDDO
                        SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
                     ENDDO
                     D_RHOLD(Q,KR,NR) = SUMR
                  ENDIF
               ENDDO

!  End KR, NR loops

            ENDDO
         ENDDO

!  Update

         DO KR = 1, N_MUQUAD
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(KR,NR) = RHOLD(KR,NR)
            ENDDO
         ENDDO
         DO Q = 1, 2
            IF ( DO_DERIV_PARS(Q) ) THEN
               DO KR = 1, N_MUQUAD
                  DO NR = 1, N_PHIQUAD
                     D_R0_QUAD_IN(Q,KR,NR) = D_RHOLD(Q,KR,NR)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO

!  Continuation point for finishing MSR

 67      CONTINUE

!  Add to total

          REFLEC   = REFLEC + REFLEC_S
          D_REFLEC = D_REFLEC + D_REFLEC_S

!  End scattering order loop

      ENDDO

!  Compute total

      COXMUNK_KERNEL = REFLEC
      DO Q = 1, 2
         IF ( DO_DERIV_PARS(Q) ) THEN
            COXMUNK_DERIVATIVES(Q) = D_REFLEC(Q)
         ENDIF
      ENDDO

!  debug

!      write(34,'(1p4e14.5)') coxmunk_kernel, &
!        dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION_MSR_PLUS

!

      SUBROUTINE BPDFSOIL_FUNCTION_PLUS  &
   ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
     BPDFSOIL_KERNEL, BPDFSOIL_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, TWO, HALF, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: BPDFSOIL_KERNEL
      REAL(fpk), intent(out) :: BPDFSOIL_DERIVATIVES ( MAXPARS )

!  Local variables

      REAL(fpk)  :: REFSQ, Z, Z1, Z2,  FP, XPHI, CKPHI, ATTEN, FP1, DFP
      REAL(fpk)  :: sgamma, cgamma, calpha, calpha_sq, salpha

!  F-.M. Breon BPDF Soil model (2009).
!   This is just the (1,1) component - no polarization

!  Initialise

      BPDFSOIL_KERNEL      = zero
      BPDFSOIL_DERIVATIVES = zero
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_PLUS ( REFSQ, Z2, FP, DFP ) ; DFP = TWO * PARS(1) * DFP

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + DSQRT((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)

!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = HALF * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = sqrt(one - calpha_sq)

      Fp1 = 0.25_fpk / xi / xj

!  BRDF  with attenuation factor

      cgamma = Z2
      sgamma = sqrt ( one - cgamma * cgamma )
      atten  = one - sgamma
      BPDFSOIL_KERNEL = Fp1 * FP * atten

!  Derivatives

      IF ( DO_DERIV_PARS(1) )  BPDFSOIL_DERIVATIVES(1) = Fp1 * DFP * atten

!     Finish

      RETURN
      END SUBROUTINE BPDFSOIL_FUNCTION_PLUS

!

      SUBROUTINE BPDFVEGN_FUNCTION_PLUS  &
   ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
     BPDFVEGN_KERNEL, BPDFVEGN_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, TWO, HALF, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: BPDFVEGN_KERNEL
      REAL(fpk), intent(out) :: BPDFVEGN_DERIVATIVES ( MAXPARS )

!  Local variables

      REAL(fpk)  :: REFSQ, Z, Z1, Z2, FP, FP0, XPHI, CKPHI, ATTEN, PROJECTIONS
      REAL(fpk)  :: sgamma, cgamma, calpha, calpha_sq, salpha
      REAL(fpk)  :: PLEAF, GS, GV, DFP

!  Data coefficients

      REAL(fpk)  :: PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098_fpk,  0.011187479_fpk, &
                               0.043329567_fpk, 0.19262991_fpk/

!  F-.M. Breon BPDF vegetation model (2009).
!   This is just the (1,1) component - no polarization

!  Initialize

      BPDFVEGN_KERNEL      = zero
      BPDFVEGN_DERIVATIVES = zero
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_PLUS ( REFSQ, Z2, FP, DFP ) ; DFP = TWO * PARS(1) * DFP

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + DSQRT((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)

!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)

!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = HALF * (xi + xj) / Z2
      calpha_sq = calpha*calpha
      salpha    = sqrt(one - calpha_sq)

! Projection of leaf surface to outgoing direction

      gv = PLAGIOPHILE_COEFFS(1) + xi * &
          (PLAGIOPHILE_COEFFS(2) + xi * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xi))

! Projection of leaf surface to incident direction

      gs = PLAGIOPHILE_COEFFS(1) + xj * &
          (PLAGIOPHILE_COEFFS(2) + xj * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xj))

! Probability of leaf orientation (plagiophile distr.)

      Pleaf = 16.0_fpk * calpha_sq * salpha  / pie

! Polarization model for vegetation

      PROJECTIONS =  Gv/xi + Gs/xj
      Fp0 = 0.25_fpk * PLEAF / xi / xj / PROJECTIONS

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = sqrt ( one - cgamma * cgamma )
      atten  = one - sgamma
      BPDFVEGN_KERNEL = Fp0 * FP * atten

!  Derivatives

      IF ( DO_DERIV_PARS(1) )  BPDFVEGN_DERIVATIVES(1) = Fp0 * DFP * atten

!     Finish

      RETURN
      END SUBROUTINE BPDFVEGN_FUNCTION_PLUS

!

      SUBROUTINE BPDFNDVI_FUNCTION_PLUS &
   ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
     BPDFNDVI_KERNEL, BPDFNDVI_DERIVATIVES )

!  include file of constants

      USE LIDORT_PARS, only : fpk, ZERO, ONE, TWO, HALF, MINUS_ONE, PIE

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: BPDFNDVI_KERNEL
      REAL(fpk), intent(out) :: BPDFNDVI_DERIVATIVES ( MAXPARS )

!  Local variables

      REAL(fpk)  :: REFSQ, Z, Z1, Z2, FP, XPHI, CKPHI, ATTEN, NDVI, EXPNDVI, C
      REAL(fpk)  :: DFP, DNDVI, DEXPNDVI, KERNEL
      REAL(fpk)  :: sgamma, cgamma

!  F-.M. Breon BPDF NDVI model (2009).
!   This is just the (1,1) component - no polarization

!  Initialise

      BPDFNDVI_KERNEL      = zero
      BPDFNDVI_DERIVATIVES = zero
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_PLUS ( REFSQ, Z2, FP, DFP ) ; DFP = TWO * PARS(1) * DFP

!  PARS(2) = NDVI
!  Exponential of the NDVI ( Out of range values default to zero )

      NDVI = PARS(2) ; DNDVI = ONE
      IF ( NDVI .GT. One .or. NDVI .lt. MINUS_ONE) THEN
        NDVI = zero ; DNDVI = zero
      ENDIF
      EXPNDVI  = EXP ( - NDVI )
      DEXPNDVI = - EXPNDVI * DNDVI

! attenuation factor

      cgamma = Z2
      sgamma = sqrt ( one - cgamma * cgamma )
      atten  = exp ( - sgamma / cgamma )

!  PARS(3) = Scaling Factor

      C = PARS(3) 
      KERNEL = 0.25_fpk * atten / ( xi + xj )
      BPDFNDVI_KERNEL = KERNEL * C * FP * EXPNDVI

!  Derivatives

      IF ( DO_DERIV_PARS(1) ) BPDFNDVI_DERIVATIVES(1) = KERNEL * C * DFP * EXPNDVI
      IF ( DO_DERIV_PARS(2) ) BPDFNDVI_DERIVATIVES(2) = KERNEL * C * FP  * DEXPNDVI
      IF ( DO_DERIV_PARS(3) ) BPDFNDVI_DERIVATIVES(3) = KERNEL * FP * EXPNDVI

!  Finish

      RETURN
      END SUBROUTINE BPDFNDVI_FUNCTION_PLUS

!

      SUBROUTINE FRESNEL_SCALAR_PLUS ( PAR1, Z2, FP, DFP )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, MINUS_ONE, ONE, HALF

      implicit none

!  Subroutine arguments

      REAL(fpk), intent(in)  :: PAR1, Z2
      REAL(fpk), intent(out) :: FP, DFP

!  Local variables

      REAL(fpk)  :: Z2_SQ_M1, H1, H2, RP, RL, DH1, DH2, DRP, DRL

!  Code

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PAR1 * Z2
      H2 = SQRT ( PAR1 + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      FP = HALF * ( RP*RP + RL*RL )

      DH1 = Z2   ; DH2 = HALF / H2
      DRP = ( ( DH1 - DH2 ) - RP * ( DH1 + DH2 ) ) / ( H1 + H2 )
      DRL =  - DH2 * ( ONE + RL ) / ( Z2 + H2 )
      DFP = RP * DRP + RL * DRL

!  End

      RETURN
      END SUBROUTINE FRESNEL_SCALAR_PLUS

!

SUBROUTINE BRDF_Generalized_Glint_plus &
         ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
           REFRAC_R, REFRAC_I, WINDSPEED,         &
           PHI_W, CPHI_W, SPHI_W,                 &
           XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
           SUNGLINT_COEFFS, DSUNGLINT_COEFFS,     &
           SUNGLINT_REFLEC, DSUNGLINT_REFLEC )

      use lidort_pars, only : fpk, zero, one, two, three, four, half, quarter, PIE

      implicit none

!  Subroutine Input arguments
!  --------------------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Real and imaginary parts of refractive index

      real(fpk), intent(in)    :: REFRAC_R
      real(fpk), intent(in)    :: REFRAC_I

!  Windspeed m/s

      real(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      real(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Incident and reflected ddirections: sines/cosines. Relative azimuth (angle in radians)

      real(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SPHI

!  Subroutine output arguments
!  ---------------------------

!   Glitter reflectance

      real(fpk), intent(out)   :: SUNGLINT_REFLEC
      real(fpk), intent(out)   :: dSUNGLINT_REFLEC

!  Cox-Munk Coefficients. Intent(inout).

      real(fpk), intent(inout) :: SUNGLINT_COEFFS(7)
      real(fpk), intent(inout) :: DSUNGLINT_COEFFS(7)

!  Local arguments
!  ---------------

      real(fpk), PARAMETER :: CRITEXP = 88.0_fpk
      real(fpk), PARAMETER :: six = two * three, twentyfour = six * four

!  Local variables

      real(fpk)  :: B, ZX, ZY, Z, Z1, Z2, XMP
      real(fpk)  :: TILT, TANTILT, TANTILT_SQ, COSTILT
      real(fpk)  :: ARGUMENT, PROB, FAC2, COEFF, VAR, WSigC, WSigU
      real(fpk)  :: XE, XN, XE_sq, XN_sq, XE_sq_1, XN_sq_1
      real(fpk)  :: XPHI, CKPHI, SKPHI, XPHI_W, CKPHI_W, SKPHI_W
      real(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      real(fpk)  :: SHADOWI, SHADOWR, SHADOW

      real(fpk)  :: dARGUMENT, dPROB, dCOEFF, dVAR, dWSigC, dWSigU
      real(fpk)  :: AN, AE, dXE, dXN, dXE_sq, dXN_sq, EXPO, dEXPO
      real(fpk)  :: T3, T4, T5, T6, T7, dT3, dT4, dT5, dT6, dT7

!  Initialise output

      SUNGLINT_REFLEC  = ZERO
      dSUNGLINT_REFLEC = ZERO

!  COmpute coefficients, according to 6S formulation

      IF ( DO_COEFFS ) THEN
         SUNGLINT_COEFFS = zero ; DSUNGLINT_COEFFS = zero
         IF ( DO_ISOTROPIC ) THEN
            SUNGLINT_COEFFS(1)  = 0.003_fpk + 0.00512_fpk * WINDSPEED
            DSUNGLINT_COEFFS(1) = 0.00512_fpk
         ELSE
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00192_fpk * WINDSPEED ! sigmaC
            SUNGLINT_COEFFS(2) =             0.00316_fpk * WINDSPEED ! sigmaU
            SUNGLINT_COEFFS(3) = 0.010_fpk - 0.00860_fpk * WINDSPEED ! C21
            SUNGLINT_COEFFS(4) = 0.040_fpk - 0.03300_fpk * WINDSPEED ! C03
            SUNGLINT_COEFFS(5) = 0.400_fpk                           ! C40
            SUNGLINT_COEFFS(6) = 0.230_fpk                           ! C04
            SUNGLINT_COEFFS(7) = 0.120_fpk                           ! C22
            DSUNGLINT_COEFFS(1) = 0.00192_fpk
            DSUNGLINT_COEFFS(2) = 0.00316_fpk 
            DSUNGLINT_COEFFS(3) = - 0.00860_fpk
            DSUNGLINT_COEFFS(4) = - 0.03300_fpk
         ENDIF
         DO_COEFFS = .false.
      ENDIF

!  Local angles

      XPHI   = PIE - PHI       ! Not used
!     CKPHI  = - CPHI          ! Original, not correct.

      CKPHI  = + CPHI
      SKPHI  = + SPHI

      XPHI_W  = PHI_W
      CKPHI_W = CPHI_W
      SKPHI_W = SPHI_W

!  Tilt angle

      B  = ONE / ( XI + XJ )
      ZX = - SXI * SKPHI * B
      ZY = ( SXJ + SXI * CKPHI ) * B
      TANTILT_SQ = ZX * ZX + ZY * ZY
      TANTILT    = SQRT ( TANTILT_SQ )
      TILT       = ATAN(TANTILT)
      COSTILT    = COS(TILT)

!  Scatter angle

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)

!  Fresnel
!  -------

       CALL BRDF_Fresnel_Complex ( REFRAC_R, REFRAC_I, Z2, XMP )

!  Anisotropic
!  -----------

      IF ( .not. DO_ISOTROPIC ) THEN

!  Variance

         WSigC = ONE / Sqrt(SUNGLINT_COEFFS(1)) ; dWSigC = - half * WSigC * WSigC * WSigC * DSUNGLINT_COEFFS(1)
         WSigU = ONE / Sqrt(SUNGLINT_COEFFS(2)) ; dWSigU = - half * WSigU * WSigU * WSigU * DSUNGLINT_COEFFS(2)
         VAR   = WSigC * WSigU * HALF ; dVAR = half * ( dWSigC * WSigU + WSigC * dWSigU )
         
!  angles

         AE = (  CKPHI_W * ZX + SKPHI_W * ZY )
         AN = ( -SKPHI_W * ZX + CKPHI_W * ZY )
         XE = AE * WSigC ; XE_sq = XE * XE ; XE_sq_1 = xe_sq - one ; dXE = AE * dWSigC ; dXE_sq = two * dXE * XE
         XN = AN * WSigU ; XN_sq = XN * XN ; XN_sq_1 = xn_sq - one ; dXN = AN * dWSigU ; dXN_sq = two * dXN * XN

!  GC Coefficient

         T3  = XE_sq_1 * XN * half
         dT3 = ( XE_sq_1 * dXN + dXE_sq * XN ) * half
         T4  = ( XN_sq - three ) * XN / six
         dT4 = ( ( XN_sq - three ) * dXN + dXN_sq * XN ) / six
         T5  = ( XE_sq * XE_sq - six * XE_sq + three ) / twentyfour
         dT5 = ( two * dXE_sq * XE_sq - six * dXE_sq ) / twentyfour
         T6  = ( XN_sq * XN_sq - six * XN_sq + three ) / twentyfour
         dT6 = ( two * dXN_sq * XN_sq - six * dXN_sq ) / twentyfour
         T7  = XE_sq_1 * XN_sq_1 / four
         dT7 = ( dXE_sq * XN_sq_1 + XE_sq_1 * dXN_sq ) / four

         Coeff  = ONE - SUNGLINT_COEFFS(3) * T3 &
                      - SUNGLINT_COEFFS(4) * T4 &
                      + SUNGLINT_COEFFS(5) * T5 &
                      + SUNGLINT_COEFFS(6) * T6 &
                      + SUNGLINT_COEFFS(7) * T7
         dCoeff =  - dSUNGLINT_COEFFS(3) * T3 - SUNGLINT_COEFFS(3) * dT3 &
                   - dSUNGLINT_COEFFS(4) * T4 - SUNGLINT_COEFFS(4) * dT4 &
                                              + SUNGLINT_COEFFS(5) * dT5 &
                                              + SUNGLINT_COEFFS(6) * dT6 &
                                              + SUNGLINT_COEFFS(7) * dT7

!  Probability and finish

         ARGUMENT  = (  XE_sq  +  XN_sq ) * HALF
         dARGUMENT = ( dXE_sq  + dXN_sq ) * HALF
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            EXPO = EXP ( - ARGUMENT ) ; dEXPO = - dARGUMENT * EXPO
            PROB = COEFF * EXPO * VAR ; dPROB =  dCOEFF * EXPO * VAR + COEFF * dEXPO * VAR + COEFF * EXPO * dVAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC  = XMP * PROB  * FAC2
            dSUNGLINT_REFLEC = XMP * dPROB * FAC2
         ENDIF

      ENDIF

!  Isotropic
!  ---------

      IF ( DO_ISOTROPIC ) THEN

!  Compute Probability and finish

         VAR   = SUNGLINT_COEFFS(1) ; dVAR = dSUNGLINT_COEFFS(1)
         ARGUMENT = TANTILT_SQ / VAR
         dARGUMENT = - ARGUMENT * dVAR / VAR
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            EXPO = EXP ( - ARGUMENT ) ; dEXPO = - dARGUMENT * EXPO
            PROB = EXPO / VAR ; dPROB =  ( dEXPO - PROB * dVAR ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC  = XMP * PROB  * FAC2
            dSUNGLINT_REFLEC = XMP * dPROB * FAC2
         ENDIF

      ENDIF

!  No Shadow code if not flagged

      IF ( .not. DO_SHADOW  ) RETURN

!  Shadow code

      S1 = SQRT ( VAR / PIE )
      S3 = ONE / ( SQRT(VAR) )
      S2 = S3 * S3

      XXI  = XI*XI
      DCOT = XI / SQRT ( ONE - XXI )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DERFC_E ( DCOT * S3 )
!      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function usage in SLEAVE code
      SHADOWI = HALF * ( S1 * T1 / DCOT - T2 )

      XXJ  = XJ*XJ
      DCOT = XJ / SQRT ( ONE - XXJ )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DERFC_E ( DCOT * S3 )
!      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function usage in SLEAVE code
      SHADOWR = HALF * ( S1 * T1 / DCOT - T2 )

      SHADOW = ONE / ( ONE + SHADOWI + SHADOWR )
      SUNGLINT_REFLEC = SUNGLINT_REFLEC * SHADOW

!     Finish

      RETURN
END SUBROUTINE BRDF_Generalized_Glint_plus

!

subroutine BRDF_WhiteCap_Reflectance_plus &
    ( WindSpeed, Wavelength, &
      WC_Reflectance, WC_Lambertian, DWC_Reflectance, DWC_Lambertian )

   use lidort_pars, only : fpk, zero

!  Stand-alone routine for computing the WhiteCap Reflectance
!   Based on 6S code, as updated by A. Sayer (2011)

!  Linearization with respect to Wind-speed

!   Made compatible with LIDORT SURFACE LEAVING code
!   renamed for BRDF code (useful if BRDF and SLEAVE are operating together)
!   R. Spurr, 23 April 2014, 28 April 2014

   implicit none

!  Inputs
!    (Wind speed in [m/s], Wavelength in Microns)

   real(fpk), intent(in)  :: WindSpeed
   real(fpk), intent(in)  :: Wavelength

!  output

   real(fpk), intent(out) :: WC_Reflectance
   real(fpk), intent(out) :: WC_Lambertian
   real(fpk), intent(out) :: DWC_Reflectance
   real(fpk), intent(out) :: DWC_Lambertian

!  Data
!  ----

!  Single precision

   real :: Effective_WCRef(39)

! effective reflectance of the whitecaps (Koepke, 1984)
! These are the original values - superseded, A Sayer 05 Jul 2011.
!      data Effective_WCRef/ &
!     0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190,&
!     0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,&
!     0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,&
!     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

! effective reflectance of the whitecaps (Frouin et al, 1996)
! Assume linear trends between the node points they give
! This is the spectral shape

      data Effective_WCRef/ &
     1.000,1.000,1.000,1.000,0.950,0.900,0.700,0.550,0.500,0.450,&
     0.400,0.350,0.300,0.250,0.200,0.150,0.100,0.050,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

!  Local variables
!  ---------------

!  Single precision in the original code

   integer :: iwl, iref
   real    :: Wlb, DWlb, WLP, Ref(39), wspd, wl, Ref_i, Rwc, DRwc

!  Initialize

   WC_Reflectance = zero
   WC_Lambertian  = zero
   DWC_Lambertian = zero

!  Single precision inputs in the original

   wspd = real(WindSpeed)
   wl   = real(Wavelength)

!  Scale data for value of 0.22 in the midvisible.

   DO iref = 1,39
      Ref(iref) = 0.22 * Effective_WCRef(iref)
   ENDDO

!  COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)

   Wlb    = 0.0 ; DWlb = 0.0
   IF (wspd .le. 9.25) THEN
       Wlb  = 0.01*((3.18e-03)*((wspd-3.7)**3.0))
       DWlb = 0.03*((3.18e-03)*((wspd-3.7)**2.0))
   ELSE IF (wspd .gt. 9.25) THEN
       Wlb  = 0.01*((4.82e-04)*((wspd+1.8)**3.0))
       DWlb = 0.03*((4.82e-04)*((wspd+1.8)**2.0))
   END IF

! Original whitecap calculation - superseded, A. Sayer 05 Jul 2011.
!      W=2.95e-06*(wspd**3.52)

!  Find data point, Linearly interpolate

   iwl   = 1+int((wl-0.2)/0.1)
   wlp   = 0.5+(iwl-1)*0.1
   Ref_i = Ref(iwl+1) + ( wl-wlp)/0.1*(Ref(iwl)-Ref(iwl+1))
   Rwc   = Wlb*Ref_i
   DRwc  = DWlb*Ref_i

!  Final values

   WC_Lambertian   = real(Wlb,fpk)
   DWC_Lambertian  = real(DWlb,fpk)
   WC_Reflectance  = real(Rwc,fpk)
   DWC_Reflectance = real(DRwc,fpk)

!  Finish

   return
end subroutine BRDF_WhiteCap_Reflectance_plus

! end Module

      END MODULE brdf_LinSup_kernels_m
