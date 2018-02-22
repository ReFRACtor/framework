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
! # Auxiliary subroutines for Surface-Leaving                   #
! #                                                             #
! #               * Sleave_GAULEG                               #
! #               * HOMEGROWN_ERRFUNC                           #
! #                                                             #
! ###############################################################

      MODULE sleave_sup_aux_m

      PUBLIC
      CONTAINS


SUBROUTINE Sleave_GAULEG(X1,X2,X,W,N)

      IMPLICIT NONE
      INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

      INTEGER  , intent(in)  :: N
      REAL(fpk), intent(in)  :: X1, X2
      REAL(fpk), intent(out) :: X(N),W(N)

      INTEGER     :: I, M, J
      REAL(fpk)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(fpk), PARAMETER :: EPS = 3.0D-14

      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)

      DO I=1,M
            Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z1 = 0.0d0
            DO WHILE (DABS(Z-Z1).GT.EPS)
                  P1=1.D0
                  P2=0.D0
                  DO J=1,N
                        P3=P2
                        P2=P1
                        P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
                  ENDDO
                  PP=N*(Z*P1-P2)/(Z*Z-1.D0)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            X(I)=XM-XL*Z
            X(N+1-I)=XM+XL*Z
            W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
            W(N+1-I)=W(I)
      ENDDO
      RETURN
      END SUBROUTINE Sleave_GAULEG

!

      SUBROUTINE HOMEGROWN_ERRFUNC ( X )

      IMPLICIT NONE

      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

      REAL(fpk), intent(inout) :: x

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7. 

      REAL(fpk) :: t,z,xin

      z = abs(x) ; xin = x
      t = 1.d0/(1.d0+0.5d0*z)
      x = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
              t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t* &
              (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+t* &
              .17087277d0)))))))))
      if ( xin .lt. 0.d0 ) x = 2.d0 - x

      return
      END SUBROUTINE HOMEGROWN_ERRFUNC

!  End module

      END MODULE sleave_sup_aux_m

