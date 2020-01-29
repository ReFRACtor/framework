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
! ###########################################################

! ###########################################################
! #                                                         #
! #                 FIRST-ORDER MODEL                       #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)        #
! #                                                         #
! #  This Version :   1.4 F90                               #
! #  Release Date :   August 2013                           #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3a, 29 October  2012, Obsgeom Multi-geom.   #
! #   Version 1.3b, 24 January  2013, BRDF/SL Supplements   #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_ScalarSS_spherfuncs_m

!  Legendre polynomials

!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)

public

contains

SUBROUTINE FO_ScalarSS_spherfuncs ( STARTER, MAXMOMS, MAXGEOMS, NMOMS, NGEOMS, DF1, DF2, COSSCAT, SS_PLEG )

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  I/O

   LOGICAL  , intent(inout) :: STARTER
   INTEGER  , intent(in)    :: MAXMOMS,  NMOMS
   INTEGER  , intent(in)    :: MAXGEOMS, NGEOMS
   REAL(fpk), intent(in)    :: COSSCAT(MAXGEOMS)
   REAL(fpk), intent(out)   :: SS_PLEG(0:MAXMOMS,MAXGEOMS)
   REAL(fpk), intent(inout) :: DF1(MAXMOMS)
   REAL(fpk), intent(inout) :: DF2(MAXMOMS)

!  Local

   integer  :: L, V ; real(fpk) :: MU
     
!  Help arrays

   IF ( STARTER ) THEN
      DF1(1) = 0.0d0 ; DF2(1) = 0.0d0
      DO L = 2, NMOMS
        DF1(L) = DBLE(2*L-1) / DBLE(L)
        DF2(L) = DBLE(L-1)   / DBLE(L)
      ENDDO
      STARTER = .false.
   ENDIF

!  Legendre

   DO V = 1, NGEOMS
     MU = COSSCAT(V)
     SS_PLEG(0,V) = 1.0d0
     SS_PLEG(1,V) = MU
     DO L = 2, NMOMS
        SS_PLEG(L,V) = DF1(L) * SS_PLEG(L-1,V) * MU - DF2(L) * SS_PLEG(L-2,V)
     ENDDO
   ENDDO

!  Finish

   RETURN
END SUBROUTINE FO_ScalarSS_spherfuncs

!  Finish

end module FO_ScalarSS_spherfuncs_m

