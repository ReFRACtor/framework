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

      module LIDORT_LinOutputs_def

!  Version 3.7, Internal threading removed, 02 May 2014

!  This module contains the following LIDORT output structures,
!  with intents :

!         LIDORT_LinAtmos    nested in LIDORT_LinOutputs
!          LIDORT_LinSurf    nested in LIDORT_LinOutputs
!       LIDORT_LinOutputs    Intent(Out)

      use LIDORT_PARS, only : fpk, MAX_GEOMETRIES, MAX_DIRECTIONS, &
                              MAXBEAMS, MAX_USER_STREAMS, MAX_USER_LEVELS,     &
                              MAXLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinAtmos

!  Column weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_COLUMNWF

!  Mean-intensity and Flux, column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_MINT_COLUMNWF
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_FLUX_COLUMNWF

!  Profile weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_PROFILEWF

!  Mean-intensity and flux, profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_MINT_PROFILEWF
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_FLUX_PROFILEWF

!  Blackbody weighting functions, 18 March 2014, Introduced for Version 3.7

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS) :: TS_ABBWFS_JACOBIANS
      REAL(fpk), dimension ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS)                :: TS_ABBWFS_FLUXES

      END TYPE LIDORT_LinAtmos

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinSurf

!  Surface weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_SURFACEWF

!  Mean-intensity and flux, surface weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_MINT_SURFACEWF
      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_FLUX_SURFACEWF

!  Blackbody weighting functions, 18 March 2014, Introduced for Version 3.7

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS) :: TS_SBBWFS_JACOBIANS
      REAL(fpk), dimension ( MAX_USER_LEVELS, 2, MAX_DIRECTIONS)                :: TS_SBBWFS_FLUXES

      END TYPE LIDORT_LinSurf

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinOutputs

      TYPE(LIDORT_LinAtmos) :: Atmos
      TYPE(LIDORT_LinSurf)  :: Surf

      END TYPE LIDORT_LinOutputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_LinAtmos,&
                LIDORT_LinSurf,&
                LIDORT_LinOutputs

      end module LIDORT_LinOutputs_def
