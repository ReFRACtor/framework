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

      module BRDF_Sup_Outputs_def

!  This module contains the following structures:

!  BRDF_Sup_Outputs - Intent(In) for LIDORT,
!                     Intent(Out) for BRDF_Sup

      use LIDORT_PARS, only : fpk, MAXMOMENTS, MAXSTREAMS, &
                              MAXBEAMS, MAX_USER_STREAMS, &
                              MAX_USER_RELAZMS, MAX_MESSAGES

      implicit none

! #####################################################################
! #####################################################################

      type BRDF_Sup_Outputs

!  Exact (direct bounce) BRDF
!  mick fix 12/29/2014 - Name changed from EXACTDB --> DBOUNCE

      REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_USER_RELAZMS, &
        MAXBEAMS ) :: BS_DBOUNCE_BRDFUNC

!  Fourier components of BRDF, in the following order
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )   :: BS_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS ) :: BS_BRDF_F
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )   :: BS_USER_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS ) :: BS_USER_BRDF_F

!  Emissivity

      REAL(fpk), dimension ( MAXSTREAMS )       :: BS_EMISSIVITY
      REAL(fpk), dimension ( MAX_USER_STREAMS ) :: BS_USER_EMISSIVITY

!  Revision 12 august 2014. Add WSA and BSA output (one SZA only for the latter)
!  11/27/14. Added the Kernel output

      REAL(fpk) :: BS_WSA_CALCULATED, BS_WSA_KERNELS(3)
      REAL(fpk) :: BS_BSA_CALCULATED, BS_BSA_KERNELS(3)

      end type BRDF_Sup_Outputs

! #####################################################################
! #####################################################################

      TYPE BRDF_Input_Exception_Handling

!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: BS_STATUS_INPUTREAD
      INTEGER      :: BS_NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_INPUTACTIONS

      END TYPE BRDF_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE BRDF_Output_Exception_Handling

!  Exception handling for Output. New code, 02 April 2014
!     Message Length should be at least 120 Characters

      INTEGER      :: BS_STATUS_OUTPUT
      INTEGER      :: BS_NOUTPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_OUTPUTMESSAGES

      END TYPE BRDF_Output_Exception_Handling

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: BRDF_Sup_Outputs,              &
                BRDF_Input_Exception_Handling, &
                BRDF_Output_Exception_Handling

      end module BRDF_Sup_Outputs_def
