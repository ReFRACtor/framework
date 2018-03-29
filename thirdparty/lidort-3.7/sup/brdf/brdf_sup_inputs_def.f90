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

      module BRDF_Sup_Inputs_def

!  Version 3.6 Notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     New OG inputs are :
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control
!     Added Overall Exact flag for better control

! #####################################################################
! #####################################################################

!  BRDF Upgrades for Version 3.7
!  -----------------------------

!  A. White-sky and Black-sky scaling options
!  ==========================================

!  WSA and BSA scaling options.
!   first introduced 14 April 2014, Revised, 17-22 April 2014
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  These options are mutually exclusive. If either is set, the BRDF code
!  will automatically perform an albedo calculation (either WS or BS) for
!  the (1,1) component of the complete 3-kernel BRDF, then normalize the
!  entire BRDF with this albedo before scaling up with the externally chosen
!  WSA or BSA (from input).

!  Additional exception handling has been introduced to make sure that
!  the spherical or planar albedos for a complete 3-kernel BRDF are in
!  the [0,1] range. Otherwise, the WSA/BSA scaling makes no sense. It is
!  still the case that some of the MODIS-tpe kernels give NEGATIVE albedos.

!  The albedo scaling process has been linearized for all existing surface
!  property Jacobians available to the BRDF linearized supplement. In
!  addition, it is also possible to derive single Jacobians of the BRDFs
!  with respect to the WS or BS albedo - this is separate from (and orthogonal
!  to) the usual kernel derivatives.

!  B. Alternative Cox-Munk Glint Reflectance
!  =========================================

!  In conjunction with new Water-leaving code developed for the SLEAVE
!  supplement, we have given the BRDF supplement a new option to 
!  return the (scalar) Cox-Munk glint reflectance, based on code originally
!  written for the 6S code.

!  Developed and tested by R. Spurr, 21-29  April 2014
!  Based in part on Modified-6S code by A. Sayer (NASA-GSFC).
!  Validated against Modified-6S OCEABRDF.F code, 24-28 April 2014.

!  The new glint option depends on Windspeed/direction, with refractive
!  indices now computed using salinity and wavelength. There is now an 
!  optional correction for (Foam) Whitecaps (Foam). These choices come from
!  the 6S formulation.

!  Need to make sure that the wind input information is the same as that
!  use for glint calculations in the SLEAVE supplement when the Glitter
!  kernels are in use. Also, the Foam correction applied here in the
!  surface-leaving code should also be applied in the SLEAVE system..

!  Choosing this new glint option bypasses the normal kernel inputs (and
!  also the WSA/BSA inputs). Instead, a single-kernel glint reflectance
!  is calculated with amplitude 1.0, based on a separate set of dedicated
!  inputs. This option does not apply for surface emission - the solar
!  sources flag must be turned on. There is only 1 Jacobian available - 
!  with respect to the windspeed. 

!  Note that the use of the facet isotropy flag is recommended for
!  multi-beam runs. This is because the wind-direction is a function of
!  the solar angle, and including this wind-direction in the glint
!  calculations for BRDF will only work if there is just one SZA. This
!  condition is checked for internally.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################
! #####################################################################

!  This module contains the following structures:

!  BRDF_Sup_Inputs - Intent(In) for BRDF_Sup

      use LIDORT_PARS, only : fpk, MAXBEAMS, MAX_USER_RELAZMS,     &
                              MAX_USER_STREAMS, MAX_USER_OBSGEOMS, & !@@
                              MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS

      implicit none

! #####################################################################
! #####################################################################

      type BRDF_Sup_Inputs

!  Stream angle flag

      LOGICAL :: BS_DO_USER_STREAMS

!  BRDF surface flag

      LOGICAL :: BS_DO_BRDF_SURFACE

!  Surface emission

      LOGICAL :: BS_DO_SURFACE_EMISSION

!   !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL :: BS_DO_SOLAR_SOURCES
      LOGICAL :: BS_DO_USER_OBSGEOMS

!  Number of discrete ordinate streams

      INTEGER :: BS_NSTREAMS

!  number of solar beams to be processed

      INTEGER :: BS_NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk), dimension (MAXBEAMS) :: BS_BEAM_SZAS

!  user-defined relative azimuths

      INTEGER                                 :: BS_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: BS_USER_RELAZMS

!  User-defined zenith angle input

      INTEGER                                 :: BS_N_USER_STREAMS
      REAL(fpk), dimension (MAX_USER_STREAMS) :: BS_USER_ANGLES_INPUT

!  !@@ Observational geometry inputs

      INTEGER                                    :: BS_N_USER_OBSGEOMS
      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: BS_USER_OBSGEOMS

!  BRDF-specific inputs
!  --------------------

!   Number and index-list of bidirectional functions

      INTEGER                                            :: BS_N_BRDF_KERNELS
      CHARACTER (LEN=10), dimension ( MAX_BRDF_KERNELS ) :: BS_BRDF_NAMES
      INTEGER, dimension ( MAX_BRDF_KERNELS )            :: BS_WHICH_BRDF

!  Parameters required for Kernel families

      INTEGER  , dimension ( MAX_BRDF_KERNELS ) :: &
        BS_N_BRDF_PARAMETERS
      REAL(fpk), dimension ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS ) :: &
        BS_BRDF_PARAMETERS

!  Lambertian Surface control

      LOGICAL, dimension ( MAX_BRDF_KERNELS )   :: BS_LAMBERTIAN_KERNEL_FLAG

!  Input kernel amplitude factors

      REAL(fpk), dimension ( MAX_BRDF_KERNELS ) :: BS_BRDF_FACTORS

!  Number of azimuth quadrature streams for BRDF

      INTEGER :: BS_NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL :: BS_DO_SHADOW_EFFECT

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!  Rob Fix 9/25/14. Two variables replaced

!      LOGICAL :: BS_DO_EXACT
!      LOGICAL :: BS_DO_EXACTONLY
      LOGICAL :: BS_DO_DIRECTBOUNCE_ONLY
!      LOGICAL :: BS_DO_DIRECTBOUNCE_TRUNCATED

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised, 12 August 2014, added option to output calculated WSA/BSA

      LOGICAL   :: BS_DO_WSABSA_OUTPUT
      LOGICAL   :: BS_DO_WSA_SCALING
      LOGICAL   :: BS_DO_BSA_SCALING
      REAL(fpk) :: BS_WSA_VALUE, BS_BSA_VALUE

!  New Cox-Munk Glint options (bypasses the usual Kernel system)
!  -------------------------------------------------------------

!  Overall flag for this option

      LOGICAL   :: BS_DO_NewCMGLINT

!  Input Salinity in [ppt]

      REAL(fpk) :: BS_SALINITY

!  Input wavelength in [Microns]

      REAL(fpk) :: BS_WAVELENGTH

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: BS_WINDSPEED, BS_WINDDIR ( MAXBEAMS )

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: BS_DO_GlintShadow
      LOGICAL   :: BS_DO_FoamOption
      LOGICAL   :: BS_DO_FacetIsotropy

!  Multiple-scattering Glitter options
!  -----------------------------------

!  Multiple reflectance corrections for GLITTER kernels (All of them!)

      LOGICAL :: BS_DO_GLITTER_MSRCORR

!  Multiple reflectance correction for exact-term Glitter kernels only
!  mick fix 12/29/2014 - Name changed from EXACTONLY --> DBONLY

      LOGICAL :: BS_DO_GLITTER_MSRCORR_DBONLY

!  Correction order for the Multiple reflectance computations
!    ( = 0, no correction, 1, 2, 3   ec.)
!  Warning; using S > 0 can increase CPU dramatically

      INTEGER :: BS_GLITTER_MSRCORR_ORDER

!  Quadrature orders for MSRCORR

      INTEGER :: BS_GLITTER_MSRCORR_NMUQUAD
      INTEGER :: BS_GLITTER_MSRCORR_NPHIQUAD

      end type BRDF_Sup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: BRDF_Sup_Inputs

      end module BRDF_Sup_Inputs_def

