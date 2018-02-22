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

      module LIDORT_Sup_SS_def

!  Version 3.7, Internal threading removed, 02 May 2014

!  This module contains the following structures,
!  with intents :

!      LIDORT_Sup_SS_def  Intent(InOut)

      USE LIDORT_PARS, only : fpk, MAX_USER_LEVELS, MAX_GEOMETRIES, &
                              MAX_DIRECTIONS

      implicit none

! #####################################################################
! #####################################################################

      type LIDORT_Sup_SS

!  SS Intensity Results at all angles and optical depths

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAX_DIRECTIONS ) :: TS_INTENSITY_SS

!  DB Intensity Results at all angles and optical depths

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_GEOMETRIES ) :: &
        TS_INTENSITY_DB

      end type LIDORT_Sup_SS

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_Sup_SS

      end module LIDORT_Sup_SS_def
