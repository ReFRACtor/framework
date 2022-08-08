
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

module VSLEAVE_Sup_Inputs_def_m

!  Version 2.6 Notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 December 2012. 
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  1/31/21. Version 2.8.3. Several changes.

!    1. Azimuth dependence in direct SLEAVE term, controlled by flag SL_AZIMUTHDEP (configuration-file read)
!       -- Requires a new foQ interpolation routine (Interpolate_fOQ_BS3) to allow dependence

!    2.  Fourier-term dependence controlled by flag SL_DO_FOURIER_OUTPUT (configuration-file read)
!       -- if set, use new foQ interpolation routine (Interpolate_fOQ_BSF) to generate L_w for azimuth quadrature
!       -- if set, new code will generate all Fourier components necessary for diffuse-field L_w
!            ** Fourier generation code is similar to that used to create Fourier BRDF entries
!       -- if not set, code reverts to azimuth averaging of Fourier-zero component only (old way)

!    3. Explicit doublet-geometry inputs and control
!       -- flag SL_DO_DOUBLET_GEOMETRY, inputs SL_N_USER_DOUBLETS, SL_USER_DOUBLETS
!       -- Configuration-file reads

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  This module contains the following structures:

!  VSLEAVE_Sup_Inputs - Intent(In) for VSLEAVE_Sup

use VLIDORT_PARS_m, only : fpk, MAXBEAMS, MAX_USER_RELAZMS, &
                                MAX_USER_STREAMS, MAX_USER_OBSGEOMS

implicit none

! #####################################################################
! #####################################################################

type VSLEAVE_Sup_Inputs

!  General control variables
!  -------------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: SL_DO_SLEAVING

!  Isotropic flag

      LOGICAL :: SL_DO_ISOTROPIC

!  Rough-Surface flag
!    ** New, 05 October 2015. Separation of the Rough-Surface Function

      LOGICAL :: SL_DO_ROUGHSURFACE

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: SL_DO_EXACT
      LOGICAL :: SL_DO_EXACTONLY

!  Fluorescence flag

      LOGICAL :: SL_DO_FLUORESCENCE

!  Solar sources flag

      LOGICAL :: SL_DO_SOLAR_SOURCES

!  Path to VSLEAVE_DATA. New 10/28/15

      CHARACTER*200 :: SL_VSLEAVE_DATAPATH

!  Geometry and integer control
!  ----------------------------

!  Stream angle flag

      LOGICAL :: SL_DO_USER_STREAMS

!  Observational Geometry flag 
!  1/31/21, Version 2.8.3. Doublet geometry option added

      LOGICAL :: SL_DO_USER_OBSGEOMS
      LOGICAL :: SL_DO_DOUBLET_GEOMETRY

!  Number of Stokes components

      INTEGER :: SL_NSTOKES

!  Number of discrete ordinate streams

      INTEGER :: SL_NSTREAMS

!  number of solar beams to be processed

      INTEGER :: SL_NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk), dimension (MAXBEAMS) :: SL_BEAM_SZAS

!  user-defined relative azimuths

      INTEGER                                 :: SL_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: SL_USER_RELAZMS

!  User-defined zenith angle input 

      INTEGER                                 :: SL_N_USER_STREAMS
      REAL(fpk), dimension (MAX_USER_STREAMS) :: SL_USER_ANGLES_INPUT

!  Observational geometry inputs

      INTEGER                                    :: SL_N_USER_OBSGEOMS
      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: SL_USER_OBSGEOMS

!  Doublet geometry inputs
!     -- 1/31/21. Version 2.8.3. Added variables

      INTEGER                                   :: SL_N_USER_DOUBLETS
      REAL(fpk), dimension (MAX_USER_STREAMS,2) :: SL_USER_DOUBLETS

!  Water-leaving variables
!  =======================

!  Basic control
!  -------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SL_SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: SL_CHLORCONC

!  Input wavelength in [Microns]

      REAL(fpk) :: SL_WAVELENGTH

!  1/31/21. Version 2.8.3. Azimuth dependence in the water-leaving radiance
!    - Controlling flag now set in the configuration file read

      LOGICAL   :: SL_AZIMUTHDEP
      
!  1/31/21, Version 2.8.3. Fourier dependence in the water-leaving radiance
!    - Controlling flag now set in the configuration file read

      LOGICAL   :: SL_DO_FOURIER_OUTPUT

!  Rough-Surface Control
!  ---------------------

!  Changed for Version 2.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: SL_WINDSPEED, SL_WINDDIR ( MAXBEAMS )

!  Removed, Version 2.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: SL_NSTREAMS_AZQUAD

!  New for Version 2.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: SL_DO_GlintShadow
      LOGICAL   :: SL_DO_FoamOption
      LOGICAL   :: SL_DO_FacetIsotropy

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: SL_FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: SL_FL_Latitude, SL_FL_Longitude

!  Input Epoch

      INTEGER :: SL_FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: SL_FL_Amplitude755

!  Flag for using Data Gaussians

      LOGICAL :: SL_FL_DO_DataGaussian

!  Input (non-data) Gaussians

      REAL(fpk) :: SL_FL_InputGAUSSIANS(3,2)

end type VSLEAVE_Sup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VSLEAVE_Sup_Inputs

end module VSLEAVE_Sup_Inputs_def_m

