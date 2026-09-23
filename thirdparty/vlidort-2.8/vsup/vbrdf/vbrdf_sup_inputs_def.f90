
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

module VBRDF_Sup_Inputs_def_m

!  Version 2.6 Notes
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

!  VBRDF Upgrades for Version 2.7
!  ------------------------------

!  A. White-sky and Black-sky scaling options
!  ==========================================

!  WSA and BSA scaling options.
!   first introduced 02 April 2014, Revised, 14-15 April 2014
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  These options are mutually exclusive. If either is set, the VBRDF code
!  will automatically perform an albedo calculation (either WS or BS) for
!  the (1,1) component of the complete 3-kernel BRDF, then normalize the
!  entire BRDF with this albedo before scaling up with the externally chosen
!  WSA or BSA (from input).

!  Additional exception handling has been introduced to make sure that
!  the spherical or planar albedos for a complete 3-kernel BRDF are in
!  the [0,1] range. Otherwise, the WSA/BSA scaling makes no sense. It is
!  still the case that some of the MODIS-tpe kernels give NEGATIVE albedos.

!  The albedo scaling process has been linearized for all existing surface
!  property Jacobians available to the VBRDF linearized supplement. In
!  addition, it is also possible to derive single Jacobians of the BRDFs
!  with respect to the WS or BS albedo - this is separate from (and orthogonal
!  to) the usual kernel derivatives.

!  B. Alternative Cox-Munk Glint Reflectance
!  =========================================

!  In conjunction with new Water-leaving code developed for the VSLEAVE
!  supplement, we have given the VBRDF supplement a new option to 
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
!  use for glint calculations in the VSLEAVE supplement when the Glitter
!  kernels are in use. Also, the Foam correction applied here in the
!  surface-leaving code should also be applied in the VSLEAVE system..

!  Choosing this new glint option bypasses the normal kernel inputs (and
!  also the WSA/BSA inputs). Instead, a single-kernel glint reflectance
!  is calculated with amplitude 1.0, based on a separate set of dedicated
!  inputs. This option does not apply for surface emission - the solar
!  sources flag must be turned on. There is only 1 Jacobian available - 
!  with respect to the windspeed. 

!  Note that the use of the facet isotropy flag is recommended for
!  multi-beam runs. This is because the wind-direction is a function of
!  the solar angle, and including this wind-direction in the glint
!  calculations for VBRDF will only work if there is just one SZA. This
!  condition is checked for internally.

! #####################################################################
! #####################################################################

!  Upgrades for Version 2.8
!  ------------------------

!  Presence of 2 new kernels, 

!  RTK_HOTSPOT : This is the Old Ross-Thick kernel with a hot-spot modification
!     Derives from the following references.
!  F. M. Breon, F. Maignan, M. Leroy and I. Grant, 
!    "Analysis of hot spot directional siganatures measured from space",
!      J. Geophys. Res., 107, D16, 4282, (2002)
!  E. Vermote C. Justice, and F. M. Breon,
!    "Towards a generalized approach for correction of the BRDF effect in MODIS reflectances"
!      IEEE Trans. Geo. Rem. Sens., 10.1109/TGRS.2008.2005997 (2008)

!  MODFRESNEL. This is a "Modified Fresnel" kernel developed for polarized reflectances
!   Taken from the following reference.
!  P. Litvinov, O. Hasekamp and B. Cairns,
!    "Models for surface reflection of radiance and polarized radiance: Comparison
!     with airborne multi-angle photopolarimetric measurements and implications for
!     modeling top-of-atmopshere measurements,
!       Rem. Sens. Env., 115, 781-792, (2011).

!  4 parameters now allowed in Kernels

!  1/31/21, Version 2.8.3. Analytical Model for Snow BRDF.
!     -- New VLIDORT BRDF Kernel. First introduced to VLIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  1/31/21, Version 2.8.3. Doublet geometry option added

! #####################################################################
! #####################################################################

!  This module contains the following structures:

!  VBRDF_Sup_Inputs - Intent(In) for VBRDFSup

use vlidort_pars_m

implicit none

! #####################################################################
! #####################################################################

type VBRDF_Sup_Inputs

!  Top level flags (same as VLIDORT)
!  --------------------------------

!  BRDF surface flag

      LOGICAL :: BS_DO_BRDF_SURFACE

!  Stream angle flag

      LOGICAL :: BS_DO_USER_STREAMS

!   !@@ Solar sources + Observational Geometry flag !@@
!  1/31/21, Version 2.8.3. Doublet geometry option added

      LOGICAL :: BS_DO_SOLAR_SOURCES
      LOGICAL :: BS_DO_USER_OBSGEOMS
      LOGICAL :: BS_DO_DOUBLET_GEOMETRY

!  Surface emission

      LOGICAL :: BS_DO_SURFACE_EMISSION

!  Numbers and Geometry (similar to VLIDORT)
!  -----------------------------------------

!  Number of Stokes components

      INTEGER :: BS_NSTOKES

!  Number of discrete ordinate streams

      INTEGER :: BS_NSTREAMS

!  number of solar beams to be processed

      INTEGER :: BS_NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk), dimension (MAXBEAMS) :: BS_BEAM_SZAS

!  User-defined zenith angle input 

      INTEGER                                 :: BS_N_USER_STREAMS
      REAL(fpk), dimension (MAX_USER_STREAMS) :: BS_USER_ANGLES_INPUT

!  User-defined relative azimuths

      INTEGER                                 :: BS_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: BS_USER_RELAZMS

!  Observational geometry inputs

      INTEGER                                    :: BS_N_USER_OBSGEOMS
      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: BS_USER_OBSGEOMS

!  Doublet-Geometry input control and angles
!    -- 1/31/21. Version 2.8.3. (2/16/21). Newly introduced

      INTEGER                                   :: BS_N_USER_DOUBLETS
      REAL(fpk), dimension (MAX_USER_STREAMS,2) :: BS_USER_DOUBLETS

!  BRDF-specific inputs
!  --------------------

!  Number of azimuth quadrature streams for BRDF

      INTEGER :: BS_NSTREAMS_BRDF

!  Number and index-list of bidirectional functions

      INTEGER                                            :: BS_N_BRDF_KERNELS
      CHARACTER (LEN=10), dimension ( MAX_BRDF_KERNELS ) :: BS_BRDF_NAMES
      INTEGER, dimension ( MAX_BRDF_KERNELS )            :: BS_WHICH_BRDF

!  Lambertian surface control

      LOGICAL, dimension ( MAX_BRDF_KERNELS )   :: BS_LAMBERTIAN_KERNEL_FLAG

!  Input kernel amplitude factors

      REAL(fpk), dimension ( MAX_BRDF_KERNELS ) :: BS_BRDF_FACTORS

!  Parameters required for Kernel families

      INTEGER  , dimension ( MAX_BRDF_KERNELS ) :: BS_N_BRDF_PARAMETERS
      REAL(fpk), dimension ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS ) :: BS_BRDF_PARAMETERS

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL   :: BS_DO_EXACT
!      LOGICAL   :: BS_DO_EXACTONLY
      LOGICAL   :: BS_DO_DIRECTBOUNCE_ONLY

!  Shadow effect flag (only for old style Cox-Munk type kernels)

      LOGICAL :: BS_DO_SHADOW_EFFECT

!  WSA and BSA scaling options
!  ---------------------------
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Rob Fix 9/27/14; added output option flag

      LOGICAL   :: BS_DO_WSABSA_OUTPUT
      LOGICAL   :: BS_DO_WSA_SCALING
      LOGICAL   :: BS_DO_BSA_SCALING
      REAL(fpk) :: BS_WSA_VALUE, BS_BSA_VALUE

!  New Cox-Munk Glint options (bypasses the usual Kernel system)
!  -------------------------------------------------------------

!  Overall flags for this option

      LOGICAL   :: BS_DO_NewCMGLINT, BS_DO_NewGCMGLINT

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
!   Rob Fix 9/25/14, name change

      LOGICAL :: BS_DO_GLITTER_MSRCORR_DBONLY

!  Correction order for the Multiple reflectance computations
!    ( = 0, no correction, 1, 2, 3   ec.)
!  Warning; using S > 0 can increase CPU dramatically

      INTEGER :: BS_GLITTER_MSRCORR_ORDER

!  Quadrature orders for MSRCORR

      INTEGER :: BS_GLITTER_MSRCORR_NMUQUAD
      INTEGER :: BS_GLITTER_MSRCORR_NPHIQUAD

end type VBRDF_Sup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VBRDF_Sup_Inputs

end module VBRDF_Sup_Inputs_def_m

