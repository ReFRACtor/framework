
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              VLIDORT_WRITEINPUT                             #
! #              VLIDORT_WRITESCEN                              #
! #              VLIDORT_WRITE_STD_INPUT                        #
! #              VLIDORT_WRITE_SUP_BRDF_INPUT                   #
! #              VLIDORT_WRITE_SUP_SS_INPUT                     #
! #              VLIDORT_WRITE_SUP_SLEAVE_INPUT                 #
! #              VLIDORT_WRITEFOURIER                           #
! #              VLIDORT_WRITERESULTS                           #
! #                                                             #
! ###############################################################

      MODULE vlidort_writemodules_m

!  1/31/21. Version 2.8.3. 
!     -- Add 3 new input flags: DO_CLASSICAL_SOLUTION, DO_MSSTS, DO_DOUBLET_GEOMETRY
!     -- VLIDORT_WRITEINPUT, VLIDORT_WRITE_STD_INPUT, re-order Boolean Variables
!     -- VLIDORT_WRITERESULTS needs upgrading

      PRIVATE
      PUBLIC :: VLIDORT_WRITEINPUT,             &
                VLIDORT_WRITESCEN,              &
                VLIDORT_WRITE_STD_INPUT,        &
                VLIDORT_WRITE_SUP_BRDF_INPUT,   &
                VLIDORT_WRITE_SUP_SS_INPUT,     &
                VLIDORT_WRITE_SUP_SLEAVE_INPUT, &
                VLIDORT_WRITEFOURIER,           &
                VLIDORT_WRITERESULTS
                
      CONTAINS

      SUBROUTINE VLIDORT_WRITEINPUT ( IUNIT, &
        DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY, &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_RAYLEIGH_ONLY,   &
        DO_NO_AZIMUTH, DO_CLASSICAL_SOLUTION, DO_MSSTS,                &
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_LAMBERTIAN_SURFACE,     &
        DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_DIRECT_BEAM,      &
        NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     &
        N_USER_LEVELS, NGREEK_MOMENTS_INPUT, VLIDORT_ACCURACY, FLUX_FACTOR )

!  1/31/21. Version 2.8.3. Add 3 new input flags. Add also DO_OBSERVATION_GEOMETRY
!    -- 3 new flags are: DO_CLASSICAL_SOLUTION, DO_MSSTS, DO_DOUBLET_GEOMETRY
!    -- re-order some of the Boolean input arguments (first 3 lines)

!  Write to file of all control input

      USE VLIDORT_PARS_m, Only : VLIDORT_VERSION_NUMBER, FMT_SECTION, FMT_HEADING, FMT_CHAR, FMT_REAL, FMT_INTEGER

      IMPLICIT NONE

!  INPUTS

      INTEGER, INTENT (IN) ::          IUNIT

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_NO_AZIMUTH

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_RAYLEIGH_ONLY

!  1/31/21. Version 2.8.3. Add 3 new input flags, and the OBSERVATION geometry flag
!    DO_CLASSICAL_SOLUTION, DO_MSSTS, DO_DOUBLET_GEOMETRY

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_MSSTS
      LOGICAL, INTENT (IN) ::          DO_DOUBLET_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

!  Removed for Version 2.8
!      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT

      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE

      LOGICAL, INTENT (IN) ::          DO_THERMAL_EMISSION
      LOGICAL, INTENT (IN) ::          DO_SURFACE_EMISSION
      LOGICAL, INTENT (IN) ::          DO_DIRECT_BEAM

      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION, INTENT (IN) :: VLIDORT_ACCURACY
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR

!  Heading and version number

      WRITE(IUNIT, FMT_SECTION) ' Control input variables for run of VLIDORT'
      WRITE(IUNIT, FMT_CHAR)    ' VLIDORT Version number = ',VLIDORT_VERSION_NUMBER

!  General control

      WRITE(IUNIT, FMT_HEADING) ' Control integers'

      WRITE(IUNIT, FMT_INTEGER) ' Double Gaussian quadrature over [-1,0] & [0,1]'
      WRITE(IUNIT, FMT_INTEGER) '  ** Number of half space streams = ', NSTREAMS
      WRITE(IUNIT, FMT_INTEGER) '  ** Number of atmospheric layers = ', NLAYERS
      WRITE(IUNIT, FMT_INTEGER) '  ** Number of Greek moments (input) = ', NGREEK_MOMENTS_INPUT

!      IF ( DO_THERMAL_EMISSION ) THEN
!        WRITE(IUNIT, FMT_INTEGER) '  ** Number of thermal emission coefficients = ',N_THERMAL_COEFFS
!      ENDIF

      WRITE(IUNIT, FMT_HEADING) ' Flux/accuracy control'
      WRITE(IUNIT, FMT_REAL)    ' Flux constant = ', FLUX_FACTOR

!mick thermal fix - DO_NO_AZIMUTH not defined at this point
!      IF ( .NOT.DO_NO_AZIMUTH ) THEN
!        WRITE(IUNIT, FMT_REAL) ' accuracy criterion (Fourier convergence) = ',VLIDORT_ACCURACY
!      ELSE
!        WRITE(IUNIT, FMT_CHAR) '  ** No Fourier series -- Azimuth-independent term only'
!      ENDIF

!  RTE control

      WRITE(IUNIT, FMT_HEADING) ' RTE solution control'
      IF ( DO_DIRECT_BEAM ) THEN
        WRITE(IUNIT, FMT_CHAR) ' Direct beam will be included in solution'
      ELSE
        WRITE(IUNIT, FMT_CHAR) ' Direct beam will NOT be included in solution'
      ENDIF

!      IF ( DO_THERMAL_EMISSION ) THEN
!        WRITE(IUNIT, FMT_CHAR) ' Thermal Emission will be included in solution'
!      ELSE
!        WRITE(IUNIT, FMT_CHAR) ' NO Thermal emission in the solution'
!      ENDIF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(IUNIT, FMT_CHAR) ' Surface Thermal Emission will be included in solution'
      ELSE
        WRITE(IUNIT, FMT_CHAR) ' NO Surface Thermal emission in the solution'
      ENDIF

!  Surface input write

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        WRITE(IUNIT, FMT_CHAR) ' Surface will be treated as Lambertian'
      ELSE
        WRITE(IUNIT, FMT_CHAR) ' Supplementary Bidirectional surface input will be used'
      ENDIF

!  1/31/21. Version 2.8.3. Added DO_MSSTS flag 

      IF ( DO_MSSTS ) THEN
        WRITE(IUNIT, FMT_CHAR) ' Multiple scattering source terms will be calculated and outputted'
      ELSE
        WRITE(IUNIT, FMT_CHAR) ' Multiple scattering source terms will NOT be done'
      ENDIF

!  Other writes. Isotropic-only option now removed. 01/17/06.

!      IF ( DO_ISOTROPIC_ONLY ) THEN
!        WRITE(IUNIT, FMT_CHAR)' Medium is isotropic'
!      ENDIF

      IF ( DO_RAYLEIGH_ONLY ) THEN
        WRITE(IUNIT, FMT_CHAR) ' Medium has Rayleigh scattering only'
      ELSE
        WRITE(IUNIT, FMT_CHAR) ' Medium has general scattering phase function'
      ENDIF

!      IF ( DO_TRANSMITTANCE_ONLY ) THEN
!        WRITE(IUNIT, FMT_CHAR) ' RTE solution is Transmission-only (no scattering)'
!      ELSE
!        WRITE(IUNIT, FMT_CHAR) ' RTE solution is full multiple-scattering'
!      ENDIF

!  Beam particular integral control

      WRITE(IUNIT, FMT_HEADING) ' Beam solution control'

!  1/31/21. Version 2.8.3. Restored DO_CLASSICAL_SOLUTION 

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        WRITE(IUNIT, FMT_CHAR) ' Beam solution determined by classical (Chandrasekhar) method'
      ELSE
        WRITE(IUNIT, FMT_CHAR) ' Beam solution determined by Greens function method'
      ENDIF

!  Resume

      IF ( DO_PLANE_PARALLEL ) THEN
        WRITE(IUNIT, FMT_CHAR) ' Beam solution in Plane-parallel approximation'
        ELSE
        WRITE(IUNIT, FMT_CHAR) ' Beam solution in Pseudo-spherical approximation'
        WRITE(IUNIT, FMT_CHAR) '  ** attenuation will use average secant estimation'
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          WRITE(IUNIT, FMT_CHAR) '  ** Refractive geometry'
        ELSE
          WRITE(IUNIT, FMT_CHAR) '  ** Straight line geometry'
        ENDIF
      ENDIF

!  Output control

      WRITE(IUNIT, FMT_HEADING) ' Output control'

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        WRITE(IUNIT, FMT_CHAR) ' Output of intensities AND fluxes & mean intensities'
      ELSE
        IF ( DO_MVOUT_ONLY ) THEN
          WRITE(IUNIT, FMT_CHAR) ' Output for fluxes & mean intensities ONLY'
        ELSE
          WRITE(IUNIT, FMT_CHAR) ' Output for intensities ONLY'
        ENDIF
      ENDIF

      WRITE(IUNIT, FMT_INTEGER) ' Number of Solar Zenith angles = ', NBEAMS

!mick thermal fix - DO_NO_AZIMUTH not defined at this point
!      IF ( .NOT.DO_NO_AZIMUTH ) THEN
!        WRITE(IUNIT, FMT_INTEGER) &
!           ' Number of user-defined azimuth angles = ', N_USER_RELAZMS
!      ENDIF

      WRITE(IUNIT, FMT_INTEGER) ' Total number of output levels = ', N_USER_LEVELS

!  Revision Version 2.8, no more QUAD output

      IF ( DO_USER_STREAMS ) THEN
          WRITE(IUNIT, FMT_CHAR)    ' Stream angle output for user-defined angles only'
          WRITE(IUNIT, FMT_INTEGER) '  ** Total Number of output stream angles = ', N_USER_STREAMS
      ENDIF

!  1/31/21. Version 2.8.3. Added DO_DOUBLET_GEOMETRY flag, use OBservation Geometry flagin

      IF (DO_USER_STREAMS ) THEN
        IF ( DO_DOUBLET_GEOMETRY ) THEN
          WRITE(IUNIT, FMT_CHAR) ' Doublet-geometry post-processing will be performed'
        ELSE IF ( DO_OBSERVATION_GEOMETRY ) THEN
          WRITE(IUNIT, FMT_CHAR) ' Observational geometry triplets for the post-processing will be performed'
        ELSE
          WRITE(IUNIT, FMT_CHAR) ' Post-processing will be done for lattice geometry inputs'
        ENDIF
      ENDIF

!  Old code (Pre 2.8) commented out
!      IF ( DO_USER_STREAMS ) THEN
!        IF ( DO_QUAD_OUTPUT ) THEN
!          WRITE(IUNIT, FMT_CHAR) ' Stream angle output will include quadratures'
!          WRITE(IUNIT, FMT_INTEGER) '  ** Total Number of output stream angles = ', NSTREAMS + N_USER_STREAMS
!        ELSE
!          WRITE(IUNIT, FMT_CHAR) ' Stream angle output for user-defined angles only'
!          WRITE(IUNIT, FMT_INTEGER) '  ** Total Number of output stream angles = ', N_USER_STREAMS
!        ENDIF
!      ELSE
!        IF ( .NOT. DO_MVOUT_ONLY ) THEN
!          WRITE(IUNIT, FMT_CHAR)    ' Stream output at Quadrature angles only'
!          WRITE(IUNIT, FMT_INTEGER) '  ** Total Number of output stream angles = ',NSTREAMS
!        ENDIF
!      ENDIF

      WRITE(IUNIT, *)

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_WRITEINPUT

!

      SUBROUTINE VLIDORT_WRITESCEN ( &
        SUNIT, DO_DELTAM_SCALING, NSTREAMS, NLAYERS, NGREEK_MOMENTS_INPUT,  &
        N_SZANGLES, SZANGLES, N_USER_RELAZMS, USER_RELAZMS,                 &
        N_USER_VZANGLES, USER_VZANGLES, N_USER_LEVELS, USER_LEVELS,         &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        DO_THERMAL_EMISSION, N_THERMAL_COEFFS, DO_SURFACE_EMISSION, SURFBB, &
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_ANGLES, DO_NO_AZIMUTH,             &
        NMOMENTS, GREEKMAT_INDEX, TAUGRID_INPUT, OMEGA_TOTAL, GREEKMAT_TOTAL, TAUGRID )

!  Write to file of all geophysical VLIDORT input

      USE VLIDORT_PARS_m, Only : MAXMOMENTS_INPUT, MAXMOMENTS, MAXSTREAMS, MAXLAYERS, MAX_SZANGLES, MAXSTOKES_SQ, &
                                 MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_LEVELS, MAX_USER_OBSGEOMS,         &
                                 VLIDORT_VERSION_NUMBER, FMT_SECTION, FMT_HEADING, FMT_CHAR, FMT_REAL, ONE, DEG_TO_RAD

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          SUNIT
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (IN) ::          N_SZANGLES
      DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: USER_LEVELS ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      LOGICAL, INTENT (IN) ::          DO_THERMAL_EMISSION
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      LOGICAL, INTENT (IN) ::          DO_SURFACE_EMISSION
      DOUBLE PRECISION, INTENT (IN) :: SURFBB
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_ANGLES  ( MAXSTREAMS )
      LOGICAL, INTENT (IN) ::          DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          GREEKMAT_INDEX ( 6 )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID ( 0:MAXLAYERS )

!  Local variables

      INTEGER :: S, N, L

!  Heading and version number

      WRITE(SUNIT,FMT_SECTION) ' Geophysical scenario variables for run of VLIDORT'
      WRITE(SUNIT,FMT_CHAR)    ' VLIDORT Version number = ',VLIDORT_VERSION_NUMBER

!  Basic input
!  -----------

      WRITE(SUNIT,FMT_SECTION) ' Basic atmospheric input for stokes vector calculation'

!  Layer optical depths and single scatter albedos

      IF ( DO_DELTAM_SCALING ) THEN
        WRITE(SUNIT, FMT_HEADING) ' Scaled and unscaled inputs with Delta-M turned ON'
        WRITE(SUNIT, FMT_HEADING) ' Layer optical depths and single-scatter-albedos '
        WRITE(SUNIT,'(a,a25/a,a25)') &
             '      (unscaled) (scaled) ','  (unscaled)   (scaled)  ', &
             'Layer  Op-depth  Op-depth ','  s.s albedo  s.s albedo '
        WRITE(SUNIT,'(a)')' '
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(I3,T6,2(f9.5,1x),2x,2(f9.5,3x))') &
                N,TAUGRID_INPUT(N),TAUGRID(N), &
                OMEGA_TOTAL_INPUT(N),OMEGA_TOTAL(N)
        ENDDO
      ELSE
        WRITE(SUNIT, FMT_HEADING) ' Unscaled inputs only: Delta-M turned OFF'
        WRITE(SUNIT, FMT_HEADING) ' Layer optical depths and single-scatter-albedos '
        WRITE(SUNIT,'(a,a13/a,a13)') &
              '      (unscaled) ','  (unscaled) ', &
              'Layer  Op-depth  ','  s.s albedo '
        WRITE(SUNIT,'(a)')' '
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(I3,T6,f9.5,4x,f9.5)') N,TAUGRID_INPUT(N),OMEGA_TOTAL_INPUT(N)
        ENDDO
      ENDIF

!  Phase matrix elements

      IF ( DO_DELTAM_SCALING ) THEN
        WRITE(SUNIT, FMT_HEADING) ' Phase matrix elements (Greek constants), SCALED'
        DO N = 1, NLAYERS
          IF (N /= 1) WRITE(SUNIT,*)
          WRITE(SUNIT,'(a, I3/)') ' Matrix elements (A-F) for layer ',N
          DO L = 0, NMOMENTS
            WRITE(SUNIT,'(1x,I2,1x,1p6e12.4)') L,(GREEKMAT_TOTAL(L,N,GREEKMAT_INDEX(S)),S = 1, 6)
          ENDDO
        ENDDO
      ELSE
        WRITE(SUNIT, FMT_HEADING) ' Phase matrix elements (Greek constants), UNSCALED'
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(a, I3/)') ' Matrix elements (A-F) for layer ',N
          DO L = 0, NGREEK_MOMENTS_INPUT
            WRITE(SUNIT,'(1x,I2,1x,1p6e12.4)') L,(GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(S)),S = 1, 6)
          ENDDO
        ENDDO
      ENDIF

!  Commented out thermal expansion coefficients
!      IF ( DO_THERMAL_EMISSION ) THEN
!        WRITE(SUNIT,FMT_HEADING)'thermal emission coefficients'
!        WRITE(SUNIT,'(a/)') 'Layer | thermal emission expansion coeffs-->'
!        DO N = 1, NLAYERS
!          WRITE(SUNIT,'(I3,4x,10(f10.5))') N,(THERMAL_COEFFS(N,S),S=1,N_THERMAL_COEFFS)
!        ENDDO
!      ENDIF

!  Surface property

      WRITE(SUNIT,FMT_HEADING)' Surface reflecting property'
      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        WRITE(SUNIT,FMT_REAL) ' (Lambertian) surface albedo is',LAMBERTIAN_ALBEDO
        IF ( DO_SURFACE_EMISSION ) THEN
          WRITE(SUNIT,FMT_REAL) ' (Lambertian) emissivity = ',ONE-LAMBERTIAN_ALBEDO
          WRITE(SUNIT,FMT_REAL) ' Surface blackbody function is', SURFBB
        ENDIF
      ELSE
        WRITE(SUNIT,FMT_CHAR) ' BRDF Supplementary input'
      ENDIF

!  Geometry input
!  --------------

      WRITE(SUNIT,FMT_SECTION)' Geometry input'

      WRITE(SUNIT,FMT_HEADING) ' Solar zenith angles'
      WRITE(SUNIT,'(a/)') ' Number |   Angle'
      DO N = 1, N_SZANGLES
        WRITE(SUNIT,'(I3,5x,1x,f9.5)') N,SZANGLES(N)
      ENDDO

      WRITE(SUNIT,FMT_HEADING) ' Computational (quadrature) angles in the half space'
      WRITE(SUNIT,'(a/)') ' Stream  |  Angle    |   Cosine   |   Weight'
      DO N = 1, NSTREAMS
        WRITE(SUNIT,'(I3,5x,3(1x,f9.5,3x))')N,QUAD_ANGLES(N),QUAD_STREAMS(N),QUAD_WEIGHTS(N)
      ENDDO

      WRITE(SUNIT,FMT_HEADING) ' User-defined viewing zenith angles'
      WRITE(SUNIT,'(a/)')' Number  |  Angle    |   Cosine'
      DO N = 1, N_USER_VZANGLES
        WRITE(SUNIT,'(I3,5x,2(1x,f9.5,3x))') N,USER_VZANGLES(N),DCOS(USER_VZANGLES(N)*DEG_TO_RAD)
      ENDDO

      IF ( DO_NO_AZIMUTH ) THEN
        WRITE(SUNIT,FMT_HEADING) ' No azimuth angles'
      ELSE
        WRITE(SUNIT,FMT_HEADING) ' User-defined relative azimuth angles'
        WRITE(SUNIT,'(a/)')' Number |   Angle'
        DO N = 1, N_USER_RELAZMS
          WRITE(SUNIT,'(I3,5x,1x,f9.5)') N,USER_RELAZMS(N)
        ENDDO
      ENDIF

!  Output levels
!  -------------

      WRITE(SUNIT,FMT_SECTION) ' User-defined levels for output'
      WRITE(SUNIT,'(a/)')' # | Level/Layer of occurrence'
      DO N = 1, N_USER_LEVELS
        WRITE(SUNIT,'(i3,4x,f9.5)')N, USER_LEVELS(N)
      ENDDO

      WRITE(SUNIT,*)

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_WRITESCEN

!

      SUBROUTINE VLIDORT_WRITE_STD_INPUT ( DO_DEBUG_WRITE, &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY, DO_SURFACE_EMISSION,                   & ! Sources Main
        DO_TOA_ILLUMINATION, DO_BOA_ILLUMINATION, DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM,                    & ! Sources Other
        DO_FULLRAD_MODE, DO_UPWELLING, DO_DNWELLING, DO_CLASSICAL_SOLUTION, DO_MSSTS, DO_FOURIER0_NSTOKES2, & ! Main RT control
        DO_RAYLEIGH_ONLY, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION,                   & ! RT Control
        DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_NADIR,  DO_FOCORR_OUTGOING, DO_SSCORR_USEFMAT,             & ! FO (first-order)
        DO_DELTAM_SCALING, DO_DOUBLE_CONVTEST, DO_SOLUTION_SAVING,  DO_BVP_TELESCOPING,                     & ! RT Performance
        DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, & ! RT Post-processing
        DO_TOA_CONTRIBS,  DO_SPECIALIST_OPTION_1,  DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3,          & ! RT Specialist
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_EXTERNAL_WLEAVE,                     & ! Surface Control
        DO_WATER_LEAVING, DO_FLUORESCENCE, DO_TF_ITERATION, DO_WLADJUSTED_OUTPUT, TF_MAXITER, TF_CRITERION, & ! Lw Control
        TAYLOR_ORDER, NSTOKES, NSTREAMS, NLAYERS, NFINELAYERS, N_THERMAL_COEFFS, NGREEK_MOMENTS_INPUT, & ! Main numbers
        NLAYERS_NOMS, NLAYERS_CUTOFF, VLIDORT_ACCURACY, FLUX_FACTOR, EARTH_RADIUS, RFINDEX_PARAMETER,  & ! Flux/Acc/Radius
        ASYMTX_TOLERANCE, N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES, N_USER_OBSGEOMS, N_USER_LEVELS, & ! Geometry/Output control
        SZANGLES,   USER_RELAZMS,   USER_VZANGLES,   USER_OBSGEOMS,   USER_LEVELS,                     & ! Geometry/Output control
        HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, FINEGRID, GEOMETRY_SPECHEIGHT, ATMOS_WAVELENGTH, & ! Atmospheric inputs
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, FMATRIX_UP, FMATRIX_DN,            & ! Optical Inputs
        THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, SURFBB, TOA_ILLUMINATION, BOA_ILLUMINATION,               & ! Optical Inputs
        DO_WRITE_INPUT,    INPUT_WRITE_FILENAME,    DO_WRITE_FOURIER,       DO_WRITE_RESULTS,          & ! Debug control
        DO_WRITE_SCENARIO, SCENARIO_WRITE_FILENAME, FOURIER_WRITE_FILENAME, RESULTS_WRITE_FILENAME )     ! Debug control

!  3/28/14. Changes for Version 2.7. Remove LTE linearization references
!  3/1/17 . Changes for Version 2.8. Complete reorganization of argument lists
!mick fix 9/19/2017 - added TAYLOR_ORDER to argument list

!  Additional Control for SLEAVE (DO_WLADJUSTED_OUTPUT,DO_EXTERNAL_WLEAVE).
!     Introduced 3/18/19 for Version 2.8.1

!   Version 2.8.1, Controls for TOA/BOA isotropic illumination added, 3/23/19
!  --- Introduce flags and values for including TOA/BOA isotropic illumination
        
!  1/31/21. Version 2.8.3. Add 4 new input flags. 
!    -- 3 new flags are: DO_CLASSICAL_SOLUTION, DO_MSSTS, DO_DOUBLET_GEOMETRY
!    -- re-order some of the Boolean input arguments (first 8 lines)
!    -- Drop SSCORR_TRUNCATION (disabled). Add ASYMTX_TOLERANCE
!    -- (RTS 2/16/21). Add DO_FOURIER0_NSTOKES2 flag

      USE VLIDORT_PARS_m, Only : MAXMOMENTS_INPUT, MAXLAYERS, MAX_SZANGLES, MAX_USER_RELAZMS, &
                                 MAX_USER_VZANGLES, MAX_USER_LEVELS, MAX_USER_OBSGEOMS, &
                                 MAX_GEOMETRIES, MAXSTOKES_SQ, &
                                 DWFC, DWFL, DWFI, DWFI1, DWFR, DWFR1, DWFR2, DWFR3, DWFR4, DWFR1_3

      IMPLICIT NONE

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

!  1/31/21. Version 2.8.3. Add 3 new input flags. 
!    -- 3 new flags are: DO_CLASSICAL_SOLUTION, DO_MSSTS, DO_DOUBLET_GEOMETRY
!    -- re-order some of the Boolean input arguments (first 8 lines)
!    -- More commentary and categorization of inputs

!  RT Sources (thermal or solar)

      LOGICAL, INTENT(IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) ::          DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT(IN) ::          DO_SURFACE_EMISSION

!  RT Sources (specialist options, mostly constant illumminations)
!   -- Added for Version 2.8.1 as indicated

      LOGICAL, INTENT (IN) ::         DO_TOA_ILLUMINATION  ! New 3/23/19 
      LOGICAL, INTENT (IN) ::         DO_BOA_ILLUMINATION  ! New 3/23/19
      LOGICAL, INTENT(IN) ::          DO_ALBTRN_MEDIA(2)   ! New 4/26/19 
      LOGICAL, INTENT(IN) ::          DO_PLANETARY_PROBLEM ! New 4/28/19 

!  Top-level RT Control
!    1/31/21. Version 2.8.3. DO_CLASSICAL_SOLUTION and DO_MSSTS. (RTS 2/16/21).  Add DO_FOURIER0_NSTOKES2

      LOGICAL, INTENT(IN) ::          DO_FULLRAD_MODE
      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_DNWELLING
      LOGICAL, INTENT(IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT(IN) ::          DO_MSSTS
      LOGICAL, INTENT(IN) ::          DO_FOURIER0_NSTOKES2

!  Top-level RT Control

      LOGICAL, INTENT(IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT(IN) ::          DO_CHAPMAN_FUNCTION
      LOGICAL, INTENT(IN) ::          DO_RAYLEIGH_ONLY

!  First-Order control (single scattering)
!    -- mick mod 9/19/2017 - DO_FOCORR_ALONE now defined internally
!    -- 1/31/21. Version 2.8.3. Dropped DO_SSCORR_TRUNCATION
!       LOGICAL, INTENT(IN) ::          DO_FOCORR_ALONE
!       LOGICAL, INTENT(IN) ::          DO_SSCORR_TRUNCATION

      LOGICAL, INTENT(IN) ::          DO_FOCORR
      LOGICAL, INTENT(IN) ::          DO_FOCORR_EXTERNAL
      LOGICAL, INTENT(IN) ::          DO_FOCORR_NADIR
      LOGICAL, INTENT(IN) ::          DO_FOCORR_OUTGOING
      LOGICAL, INTENT(IN) ::          DO_SSCORR_USEFMAT

!  RT Performance

      LOGICAL, INTENT(IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) ::          DO_DOUBLE_CONVTEST
      LOGICAL, INTENT(IN) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT(IN) ::          DO_BVP_TELESCOPING

!  RT Geometry and output control
!    -- 1/31/21. Version 2.8.3. Added DO_DOUBLET_GEOMETRY

      LOGICAL, INTENT(IN) ::          DO_USER_VZANGLES
      LOGICAL, INTENT(IN) ::          DO_DOUBLET_GEOMETRY
      LOGICAL, INTENT(IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT(IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN) ::          DO_MVOUT_ONLY

!  Specialist options not available to general user

      !LOGICAL, INTENT(IN) ::          DO_QUAD_OUTPUT   removed, version 2.8
      LOGICAL, INTENT(IN) ::          DO_TOA_CONTRIBS
      LOGICAL, INTENT(IN) ::          DO_SPECIALIST_OPTION_1
      LOGICAL, INTENT(IN) ::          DO_SPECIALIST_OPTION_2
      LOGICAL, INTENT(IN) ::          DO_SPECIALIST_OPTION_3

!  Surface control

      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::          DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) ::          DO_SL_ISOTROPIC
      LOGICAL, INTENT(IN) ::          DO_EXTERNAL_WLEAVE

!  Surface Water-leaving control
!   -- New flag for Version 2.8.1. Introduced 3/18/19.

      LOGICAL, INTENT(IN) ::          DO_WATER_LEAVING
      LOGICAL, INTENT(IN) ::          DO_FLUORESCENCE
      LOGICAL, INTENT(IN) ::          DO_TF_ITERATION
      LOGICAL, INTENT(IN) ::          DO_WLADJUSTED_OUTPUT ! 3/18/19
      INTEGER, INTENT(IN) ::          TF_MAXITER
      DOUBLE PRECISION, INTENT(IN) :: TF_CRITERION

!  Main computational control integers

      INTEGER, INTENT(IN) ::          TAYLOR_ORDER
      INTEGER, INTENT(IN) ::          NSTOKES
      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NLAYERS
      INTEGER, INTENT(IN) ::          NFINELAYERS
      INTEGER, INTENT(IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT(IN) ::          NGREEK_MOMENTS_INPUT

!  other specialist control integers

      INTEGER, INTENT(IN) ::          NLAYERS_NOMS
      INTEGER, INTENT(IN) ::          NLAYERS_CUTOFF

!  Control numbers
!  -- 1/31/21. Version 2.8.3. Include ASYMTX_TOLERANCE output

      DOUBLE PRECISION, INTENT(IN) :: VLIDORT_ACCURACY
      DOUBLE PRECISION, INTENT(IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT(IN) :: RFINDEX_PARAMETER
      DOUBLE PRECISION, INTENT(IN) :: ASYMTX_TOLERANCE
      DOUBLE PRECISION, INTENT(IN) :: EARTH_RADIUS

!  geometry and output control (integers)

      INTEGER, INTENT(IN) ::          N_SZANGLES
      INTEGER, INTENT(IN) ::          N_USER_VZANGLES
      INTEGER, INTENT(IN) ::          N_USER_RELAZMS
      INTEGER, INTENT(IN) ::          N_USER_LEVELS
      INTEGER, INTENT(IN) ::          N_USER_OBSGEOMS

!  geometry and output control (angles/levels)

      DOUBLE PRECISION, INTENT(IN) :: SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT(IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT(IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT(IN) :: USER_LEVELS   ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT(IN) :: USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )


!  Atmospheric inputs

      DOUBLE PRECISION, INTENT(IN) :: HEIGHT_GRID      ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: PRESSURE_GRID    ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER         , INTENT(IN) :: FINEGRID         ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: GEOMETRY_SPECHEIGHT
      DOUBLE PRECISION, INTENT(IN) :: ATMOS_WAVELENGTH

!  Optical inputs (Main atmospheric)

      DOUBLE PRECISION, INTENT(IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) :: FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION, INTENT(IN) :: FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Optical inputs (surface/thermal/uniform)
!    -- TOA/BOA Fluxes (new 3/23/19). Version 2.8.1

      DOUBLE PRECISION, INTENT(IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT(IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: SURFBB
      DOUBLE PRECISION, INTENT(IN) :: TOA_ILLUMINATION
      DOUBLE PRECISION, INTENT(IN) :: BOA_ILLUMINATION

!  Code superseded, Version 2.7. Remove this
      !DOUBLE PRECISION, INTENT(IN) ::   LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
      !DOUBLE PRECISION, INTENT(IN) ::   LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  debug control

      LOGICAL, INTENT(IN) ::          DO_DEBUG_WRITE
      LOGICAL, INTENT(IN) ::          DO_WRITE_INPUT
      LOGICAL, INTENT(IN) ::          DO_WRITE_SCENARIO
      LOGICAL, INTENT(IN) ::          DO_WRITE_FOURIER
      LOGICAL, INTENT(IN) ::          DO_WRITE_RESULTS

      CHARACTER (LEN=60), INTENT(IN) :: INPUT_WRITE_FILENAME
      CHARACTER (LEN=60), INTENT(IN) :: SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60), INTENT(IN) :: FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60), INTENT(IN) :: RESULTS_WRITE_FILENAME

!  Local variables

      INTEGER :: OUTUNIT, NGEOMS
      INTEGER :: GEO,K,LAY,MOM,NGREEK,SS,SZA,ULEV,URA,UVA,UOG
      INTEGER :: GMASK(8)
      DATA GMASK / 1, 2, 5, 6, 11, 12, 15, 16 /

!  Define local variables
!  -- 1/31/21. Version 2.8.3. Expand code to allow for additional DO_DOUBLET_GEOMETRY option
!mick fix 1/5/2021 - adjusted defining of NGEOMS for observational geometry and lattice
 
      !NGEOMS = N_SZANGLES
      IF ( DO_USER_VZANGLES ) THEN
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            !NGEOMS = N_SZANGLES * N_USER_VZANGLES * N_USER_RELAZMS
            NGEOMS = N_SZANGLES
         ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            NGEOMS = N_SZANGLES * N_USER_VZANGLES
         ELSE
            !NGEOMS = N_SZANGLES
            NGEOMS = N_SZANGLES * N_USER_VZANGLES * N_USER_RELAZMS
         ENDIF
      ELSE
         NGEOMS = N_SZANGLES
      ENDIF
      
!  Open output file

      OUTUNIT = 101
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_STD_INPUT.dbg',status = 'replace')

!  Changed 16 oct 14, only output non-zero entries of GreekMat

      NGREEK = 1
      IF ( NSTOKES.eq.3 ) NGREEK = 5
      IF ( NSTOKES.eq.4 ) NGREEK = 8

!  Write all input to file
!   -- 1/31/21. Version 2.8.3. (RTS 2/16/21). Output more spaced out according to natural groupings.

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_FULLRAD_MODE        = ',DO_FULLRAD_MODE
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_EMISSION    = ',DO_THERMAL_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_EMISSION    = ',DO_SURFACE_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_PLANE_PARALLEL      = ',DO_PLANE_PARALLEL
      WRITE(OUTUNIT,DWFL)  'DO_UPWELLING           = ',DO_UPWELLING
      WRITE(OUTUNIT,DWFL)  'DO_DNWELLING           = ',DO_DNWELLING
      WRITE(OUTUNIT,DWFL)  'DO_LAMBERTIAN_SURFACE  = ',DO_LAMBERTIAN_SURFACE

      !WRITE(OUTUNIT,DWFL)  'DO_QUAD_OUTPUT         = ',DO_QUAD_OUTPUT  ! removed Version 2.8

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_MSSTS               = ',DO_MSSTS                ! New 1/31/21. Version 2.8.3
      WRITE(OUTUNIT,DWFL)  'DO_FOURIER0_NSTOKES2   = ',DO_FOURIER0_NSTOKES2    ! New 1/31/21. Version 2.8.3. (RTS 2/16/21). 

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_TOA_CONTRIBS        = ',DO_TOA_CONTRIBS
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_1 = ',DO_SPECIALIST_OPTION_1
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_2 = ',DO_SPECIALIST_OPTION_2
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_3 = ',DO_SPECIALIST_OPTION_3

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LEAVING     = ',DO_SURFACE_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_SL_ISOTROPIC        = ',DO_SL_ISOTROPIC
      WRITE(OUTUNIT,DWFL)  'DO_WATER_LEAVING       = ',DO_WATER_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_FLUORESCENCE        = ',DO_FLUORESCENCE
      WRITE(OUTUNIT,DWFL)  'DO_TF_ITERATION        = ',DO_TF_ITERATION
      WRITE(OUTUNIT,DWFL)  'DO_WLADJUSTED_OUTPUT   = ',DO_WLADJUSTED_OUTPUT   !    Introduced 3/18/19 for Version 2.8.1
      WRITE(OUTUNIT,DWFL)  'DO_TOA_ILLUMINATION    = ',DO_TOA_ILLUMINATION          !    Introduced 3/23/19 for Version 2.8.1
      WRITE(OUTUNIT,DWFL)  'DO_BOA_ILLUMINATION    = ',DO_BOA_ILLUMINATION          !    Introduced 3/23/19 for Version 2.8.1

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'TAYLOR_ORDER     = ',TAYLOR_ORDER
      WRITE(OUTUNIT,DWFI)  'NSTOKES          = ',NSTOKES
      WRITE(OUTUNIT,DWFI)  'NSTREAMS         = ',NSTREAMS
      WRITE(OUTUNIT,DWFI)  'NLAYERS          = ',NLAYERS
      WRITE(OUTUNIT,DWFI)  'NFINELAYERS      = ',NFINELAYERS
      WRITE(OUTUNIT,DWFI)  'N_THERMAL_COEFFS = ',N_THERMAL_COEFFS
      WRITE(OUTUNIT,DWFR)  'VLIDORT_ACCURACY = ',VLIDORT_ACCURACY

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NLAYERS_NOMS     = ',NLAYERS_NOMS
      WRITE(OUTUNIT,DWFI)  'NLAYERS_CUTOFF   = ',NLAYERS_CUTOFF

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'TF_MAXITER       = ',TF_MAXITER
      WRITE(OUTUNIT,DWFR)  'TF_CRITERION     = ',TF_CRITERION

!  -- 1/31/21. Version 2.8.3. Expand code to output ASYMTX_TOLERANCE

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'ASYMTX_TOLERANCE = ',ASYMTX_TOLERANCE
      
!  4/26/19. Record the Media-problem and Planetary-problem inputs

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_ALBTRN_MEDIA(1)   = ',DO_ALBTRN_MEDIA(1)
      WRITE(OUTUNIT,DWFL)  'DO_ALBTRN_MEDIA(2)   = ',DO_ALBTRN_MEDIA(2)
      WRITE(OUTUNIT,DWFL)  'DO_PLANETARY_PROBLEM = ',DO_PLANETARY_PROBLEM

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'FLUX_FACTOR      = ',FLUX_FACTOR
      
      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'TOA_ILLUMINATION = ',TOA_ILLUMINATION          !    Introduced 3/23/19 for Version 2.8.1
      WRITE(OUTUNIT,DWFR)  'BOA_ILLUMINATION = ',BOA_ILLUMINATION          !    Introduced 3/23/19 for Version 2.8.1

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_LEVELS    = ',N_USER_LEVELS

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' HEIGHT_GRID(LAY)      = ',HEIGHT_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' PRESSURE_GRID(LAY)    = ',PRESSURE_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' TEMPERATURE_GRID(LAY) = ',TEMPERATURE_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,' FINEGRID (LAY)        = ',FINEGRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'RFINDEX_PARAMETER = ',RFINDEX_PARAMETER

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' DELTAU_VERT_INPUT(LAY) = ',DELTAU_VERT_INPUT(LAY)
      END DO

!mick mod 1/5/2021 - added blank line before each new K set in GREEKMAT_TOTAL_INPUT
!                    FMATRIX_UP, & FMATRIX_DN

      DO K=1,NGREEK
        SS = GMASK(K) ; WRITE(OUTUNIT,*)
        DO LAY=1,NLAYERS
          DO MOM=0,NGREEK_MOMENTS_INPUT
            WRITE(OUTUNIT,DWFR3)  'SS = ',SS,' LAY = ',LAY,' MOM = ',MOM,&
              ' GREEKMAT_TOTAL_INPUT(MOM,LAY,SS) = ',GREEKMAT_TOTAL_INPUT(MOM,LAY,SS)
          END DO
        END DO
      END DO

!  New section for 2.8 code. F-Matrices

      DO K=1,6
        WRITE(OUTUNIT,*)
        DO GEO=1,NGEOMS
          DO LAY=1,NLAYERS
            WRITE(OUTUNIT,DWFR3)  'K = ',K,' GEO = ',GEO,' LAY = ',LAY,&
              ' FMATRIX_UP(LAY,GEO,K) = ',FMATRIX_UP(LAY,GEO,K)
          END DO
        END DO
      END DO

      DO K=1,6
        WRITE(OUTUNIT,*)
        DO GEO=1,NGEOMS
          DO LAY=1,NLAYERS
            WRITE(OUTUNIT,DWFR3)  'K = ',K,' GEO = ',GEO,' LAY = ',LAY,&
              ' FMATRIX_DN(LAY,GEO,K) = ',FMATRIX_DN(LAY,GEO,K)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'LAMBERTIAN_ALBEDO = ',LAMBERTIAN_ALBEDO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' THERMAL_BB_INPUT(LAY) = ',THERMAL_BB_INPUT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'SURFBB            = ',SURFBB

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'ATMOS_WAVELENGTH  = ',ATMOS_WAVELENGTH

!  Next two loops removed, Version 2.7
 !     DO LAY=1,NLAYERS
 !       DO I=1,2
 !         WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' I = ',I,&
 !           ' LTE_DELTAU_VERT_INPUT(I,LAY) = ',LTE_DELTAU_VERT_INPUT(I,LAY)
 !       END DO
 !     END DO
 !     DO LAY=0,NLAYERS
 !       WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
 !         ' LTE_THERMAL_BB_INPUT(LAY) = ',LTE_THERMAL_BB_INPUT(LAY)
 !     END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_DEBUG_WRITE          = ',DO_DEBUG_WRITE
      WRITE(OUTUNIT,DWFL)  'DO_WRITE_INPUT          = ',DO_WRITE_INPUT
      WRITE(OUTUNIT,DWFL)  'DO_WRITE_SCENARIO       = ',DO_WRITE_SCENARIO
      WRITE(OUTUNIT,DWFL)  'DO_WRITE_FOURIER        = ',DO_WRITE_FOURIER
      WRITE(OUTUNIT,DWFL)  'DO_WRITE_RESULTS        = ',DO_WRITE_RESULTS

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFC)  'INPUT_WRITE_FILENAME    = ',INPUT_WRITE_FILENAME
      WRITE(OUTUNIT,DWFC)  'SCENARIO_WRITE_FILENAME = ',SCENARIO_WRITE_FILENAME
      WRITE(OUTUNIT,DWFC)  'FOURIER_WRITE_FILENAME  = ',FOURIER_WRITE_FILENAME
      WRITE(OUTUNIT,DWFC)  'RESULTS_WRITE_FILENAME  = ',RESULTS_WRITE_FILENAME

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Variable'
      WRITE(OUTUNIT,'(A)') '--------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR               = ',DO_FOCORR
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_EXTERNAL      = ',DO_FOCORR_EXTERNAL
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_NADIR         = ',DO_FOCORR_NADIR
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_OUTGOING      = ',DO_FOCORR_OUTGOING
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_USEFMAT       = ',DO_SSCORR_USEFMAT

      !WRITE(OUTUNIT,DWFL)  'DO_FOCORR_ALONE         = ',DO_FOCORR_ALONE       ! Disabled Verison 2.8.1
      !WRITE(OUTUNIT,DWFL)  'DO_SSCORR_TRUNCATION    = ',DO_SSCORR_TRUNCATION ! Disabled 1/31/21. Version 2.8.3

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_EXTERNAL_WLEAVE      = ',DO_EXTERNAL_WLEAVE   !    Introduced 3/18/19 for Version 2.8.1
      WRITE(OUTUNIT,DWFL)  'DO_DOUBLE_CONVTEST      = ',DO_DOUBLE_CONVTEST
      WRITE(OUTUNIT,DWFL)  'DO_SOLAR_SOURCES        = ',DO_SOLAR_SOURCES
      WRITE(OUTUNIT,DWFL)  'DO_REFRACTIVE_GEOMETRY  = ',DO_REFRACTIVE_GEOMETRY
      WRITE(OUTUNIT,DWFL)  'DO_CHAPMAN_FUNCTION     = ',DO_CHAPMAN_FUNCTION

      WRITE(OUTUNIT,DWFL)  'DO_CLASSICAL_SOLUTION   = ',DO_CLASSICAL_SOLUTION ! New 1/31/21. Version 2.8.3

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_RAYLEIGH_ONLY        = ',DO_RAYLEIGH_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_DELTAM_SCALING       = ',DO_DELTAM_SCALING
      WRITE(OUTUNIT,DWFL)  'DO_SOLUTION_SAVING      = ',DO_SOLUTION_SAVING
      WRITE(OUTUNIT,DWFL)  'DO_BVP_TELESCOPING      = ',DO_BVP_TELESCOPING

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_USER_VZANGLES        = ',DO_USER_VZANGLES
      WRITE(OUTUNIT,DWFL)  'DO_ADDITIONAL_MVOUT     = ',DO_ADDITIONAL_MVOUT
      WRITE(OUTUNIT,DWFL)  'DO_MVOUT_ONLY           = ',DO_MVOUT_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_TRANSONLY    = ',DO_THERMAL_TRANSONLY

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_OBSERVATION_GEOMETRY = ',DO_OBSERVATION_GEOMETRY
      WRITE(OUTUNIT,DWFL)  'DO_DOUBLET_GEOMETRY     = ',DO_DOUBLET_GEOMETRY   ! New 1/31/21. Version 2.8.3

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NGREEK_MOMENTS_INPUT    = ',NGREEK_MOMENTS_INPUT

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_SZANGLES     = ',N_SZANGLES
      DO SZA=1,N_SZANGLES
        WRITE(OUTUNIT,DWFR1)  'SZA = ',SZA,' SZANGLES(SZA) = ',SZANGLES(SZA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_RELAZMS = ',N_USER_RELAZMS
      DO URA=1,N_USER_RELAZMS
        WRITE(OUTUNIT,DWFR1)  'URA = ',URA,' USER_RELAZMS(URA) = ',USER_RELAZMS(URA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_VZANGLES = ',N_USER_VZANGLES
      DO UVA=1,N_USER_VZANGLES
        WRITE(OUTUNIT,DWFR1)  'UVA = ',UVA,' USER_VZANGLES(UVA) = ',USER_VZANGLES(UVA)
      END DO

      WRITE(OUTUNIT,*)
      DO ULEV=1,N_USER_LEVELS
        WRITE(OUTUNIT,DWFR1)  'ULEV = ',ULEV,' USER_LEVELS(ULEV)  = ',USER_LEVELS(ULEV)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'GEOMETRY_SPECHEIGHT = ',GEOMETRY_SPECHEIGHT
      WRITE(OUTUNIT,DWFI)  'N_USER_OBSGEOMS     = ',N_USER_OBSGEOMS

      IF (N_USER_OBSGEOMS > 0) WRITE(OUTUNIT,*)
      DO UOG=1,N_USER_OBSGEOMS
        WRITE(OUTUNIT,DWFR1_3)  'UOG = ',UOG,&
          ' USER_OBSGEOMS(UOG,1:3) = ',USER_OBSGEOMS(UOG,1),&
                                   ',',USER_OBSGEOMS(UOG,2),&
                                   ',',USER_OBSGEOMS(UOG,3)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'EARTH_RADIUS        = ',EARTH_RADIUS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' OMEGA_TOTAL_INPUT(LAY) = ',OMEGA_TOTAL_INPUT(LAY)
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_STD_INPUT

!

      SUBROUTINE VLIDORT_WRITE_SUP_BRDF_INPUT ( &
        DO_USER_STREAMS, DO_SURFACE_EMISSION,                           &
        NSTOKES, NSTREAMS, N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS, &
        EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,      &
        EMISSIVITY,USER_EMISSIVITY)

!mick mod 9/19/2017 - added DO_USER_STREAMS & DO_SURFACE_EMISSION for better output control
!                   - added NMOMENTS for VLIDORT internal consistency

!  1/31/21. Version 2.8.3. (2/16/21). BRDF arrays directly from type structure. Restore MAXMOMENTS

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, & 
                                 MAXSTOKES_SQ, MAXMOMENTS, DWFR2, DWFR4

      IMPLICIT NONE

      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::   DO_SURFACE_EMISSION

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NSTREAMS

!  1/31/21. Version 2.8.3. (2/16/21). BRDF remove NMOMENTS (equal to 2*NSTREAMS - 1)
!      INTEGER, INTENT(IN) ::   NMOMENTS

      DOUBLE PRECISION, INTENT(IN) ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. (2/16/21). BRDF arrays directly from type structure. Restore MAXMOMENTS

      DOUBLE PRECISION, INTENT(IN) ::   BRDF_F_0  ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   BRDF_F    ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT(IN) ::   USER_BRDF_F_0 ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_BRDF_F   ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT(IN) ::   EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SS,MOM,NREFL,K,STRM,S,USTRM,IB,URA,NMOMS,&
                 STRMI,STRMJ, RMASK3(9), RMASK4(16), RMASK(16)

      data RMASK3 / 1, 2, 3, 5, 6, 7, 9, 10, 11 /
      data RMASK4 / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 /

!  Open output file

      OUTUNIT = 102
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_BRDF_INPUT.dbg',status = 'replace')

!  Define local variables
!   Changed 16 October 2014, according to correct entries in 4x4 BRDF matrices

      NREFL = NSTOKES * NSTOKES ; RMASK(1) = 1
      IF ( NSTOKES.eq.3 ) RMASK(1:NREFL) = RMASK3(1:NREFL)
      IF ( NSTOKES.eq.4 ) RMASK(1:NREFL) = RMASK4(1:NREFL)

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '----------------------'
      WRITE(OUTUNIT,'(A)') 'BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '----------------------'

      WRITE(OUTUNIT,*)
      DO K=1,NREFL
        SS = RMASK(K)
        DO IB=1,N_SZANGLES
          DO URA=1,N_USER_RELAZMS
            DO USTRM=1,N_USER_VZANGLES
              WRITE(OUTUNIT,DWFR4) &
                'SS = ',SS,' IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,&
                ' EXACTDB_BRDFUNC(SS,USTRM,URA,IB) = ',&
                  EXACTDB_BRDFUNC(SS,USTRM,URA,IB)
            END DO
          END DO
        END DO
      END DO

!  1/31/21. Version 2.8.3. (2/16/21). BRDF arrays directly from type structure
!    -- Restore full Fourier output, put back MOM loops. Notes on the usage by MC, commented output

!      MOM = 0 ; WRITE(OUTUNIT,*)
!      IF ( .NOT.DO_USER_STREAMS ) THEN
!        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING TWO BRDF INPUT SETS:'
!      ELSE
!        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING FOUR BRDF INPUT SETS:'
!      ENDIF

      NMOMS = 2 * NSTREAMS - 1

      WRITE(OUTUNIT,*)
      DO K=1,NREFL
        SS = RMASK(K)
        DO IB=1,N_SZANGLES ; DO STRM=1,NSTREAMS ; DO MOM = 0, NMOMS
           WRITE(OUTUNIT,DWFR4) &
                'SS = ',SS,' IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,&
                ' BRDF_F_0(0,SS,STRM,IB) = ',BRDF_F_0(MOM,SS,STRM,IB)
        END DO ; END DO ; END DO
      END DO

      WRITE(OUTUNIT,*)
      DO K=1,NREFL
        SS = RMASK(K)
        DO STRMJ=1,NSTREAMS ; DO STRMI=1,NSTREAMS ; DO MOM = 0, NMOMS
          WRITE(OUTUNIT,DWFR4) &
                'SS = ',SS,' STRMJ = ',STRMJ,' STRMI = ',STRMI,' MOM = ',MOM,&
                ' BRDF_F(MOM,SS,STRMI,STRMJ) = ',BRDF_F(MOM,SS,STRMI,STRMJ)
        END DO ; END DO ; END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO K=1,NREFL
          SS = RMASK(K)
          DO IB=1,N_SZANGLES ; DO USTRM=1,N_USER_VZANGLES ; DO MOM = 0, NMOMS
            WRITE(OUTUNIT,DWFR4) &
                  'SS = ',SS,' IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,&
                  ' USER_BRDF_F_0(MOM,SS,USTRM,IB) = ',USER_BRDF_F_0(MOM,SS,USTRM,IB)
          END DO ; END DO ; END DO
        END DO
        WRITE(OUTUNIT,*)
        DO K=1,NREFL
          SS = RMASK(K)
          DO STRM=1,NSTREAMS ; DO USTRM=1,N_USER_VZANGLES ; DO MOM = 0, NMOMS
            WRITE(OUTUNIT,DWFR4) &
              'SS = ',SS,' STRM = ',STRM,' USTRM = ',USTRM,' MOM = ',MOM,&
              ' USER_BRDF_F(MOM,SS,USTRM,STRM) = ',USER_BRDF_F(MOM,SS,USTRM,STRM)
          END DO ; END DO ; END DO
        END DO
      END IF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(OUTUNIT,*)
        DO S=1,NSTOKES
          DO STRM=1,NSTREAMS
            WRITE(OUTUNIT,DWFR2)  'S = ',S,' STRM = ',STRM,&
              ' EMISSIVITY (S,STRM) = ',EMISSIVITY (S,STRM)
          END DO
        END DO

        WRITE(OUTUNIT,*)
        DO S=1,NSTOKES
          DO USTRM=1,N_USER_VZANGLES
            WRITE(OUTUNIT,DWFR2)  'S = ',S,' USTRM = ',USTRM,&
              ' USER_EMISSIVITY (S,USTRM) = ',USER_EMISSIVITY(S,USTRM)
          END DO
        END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_SUP_BRDF_INPUT

!

      SUBROUTINE VLIDORT_WRITE_SUP_SS_INPUT ( &
        NSTOKES,N_USER_LEVELS,STOKES_SS,STOKES_DB)

      USE VLIDORT_PARS_m, Only : MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, DWFR3, DWFR4

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_USER_LEVELS

      DOUBLE PRECISION, INTENT(IN) ::  STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) ::  STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,S,ULEV

!  Open output file

      OUTUNIT = 103
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_SS_INPUT.dbg',status = 'replace')

!  Write all single-scatter input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------'
      WRITE(OUTUNIT,'(A)') 'SS Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '--------------------'

      WRITE(OUTUNIT,*)
      DO DIR=1,MAX_DIRECTIONS
        DO S=1,NSTOKES
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
            WRITE(OUTUNIT,DWFR4) &
              'DIR = ',DIR,' S = ',S,' GEO = ',GEO,' ULEV = ',ULEV,&
              ' STOKES_SS(ULEV,GEO,S,DIR) = ',STOKES_SS(ULEV,GEO,S,DIR)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO S=1,NSTOKES
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
          WRITE(OUTUNIT,DWFR3) &
            'S = ',S,' GEO = ',GEO,' ULEV = ',ULEV,&
            ' STOKES_DB(ULEV,GEO,S) = ',STOKES_DB(ULEV,GEO,S)
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_SUP_SS_INPUT

!

      SUBROUTINE VLIDORT_WRITE_SUP_SLEAVE_INPUT ( &
        DO_USER_STREAMS,NSTOKES,NSTREAMS,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,&
        SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)

!mick mod 9/19/2017 - added DO_USER_STREAMS for better output control
!                   - added NMOMENTS for VLIDORT internal consistency (dropped 2/8/3)

!  1/31/21. Version 2.8.3. (2/16/21). SLEAVE arrays directly from type structure (restore MAXMOMENTS)

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAX_USER_STREAMS, &
                                 MAX_USER_RELAZMS, MAXBEAMS, MAXMOMENTS, DWFR2, DWFR4

      IMPLICIT NONE

      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NSTREAMS

!  1/31/21. Version 2.8.3. (2/16/21). SLEAVE remove NMOMENTS (equal to 2*NSTREAMS - 1)
!      INTEGER, INTENT(IN) ::   NMOMENTS

      DOUBLE PRECISION, INTENT(IN) :: SLTERM_ISOTROPIC  ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) :: SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!  1/31/21. Version 2.8.3. (2/16/21). SLEAVE arrays directly from type structure. Restore MAXMOMENTS

      DOUBLE PRECISION, INTENT(IN) :: SLTERM_F_0      ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) :: USER_SLTERM_F_0 ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: S,USTRM,IB,URA,MOM,STRM,NMOMS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_SLEAVE_INPUT.dbg',status = 'replace')

!  Write all surface-leaving input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '------------------------'
      WRITE(OUTUNIT,'(A)') 'SLEAVE Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '------------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,N_SZANGLES
        DO S=1,NSTOKES
          WRITE(OUTUNIT,DWFR2) &
            'IB = ',IB,' S = ',S,&
            ' SLTERM_ISOTROPIC(S,IB) = ',SLTERM_ISOTROPIC(S,IB)
        END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,N_SZANGLES
          DO URA=1,N_USER_RELAZMS
            DO USTRM=1,N_USER_VZANGLES
              DO S=1,NSTOKES
                WRITE(OUTUNIT,DWFR4) &
                  'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,' S = ',S,&
                  ' SLTERM_USERANGLES(S,USTRM,URA,IB) = ',&
                    SLTERM_USERANGLES(S,USTRM,URA,IB)
              END DO
            END DO
          END DO
        END DO
      END IF

!  1/31/21. Version 2.8.3. SLEAVE arrays directly from type structure
!    -- Restore full Fourier output, put back MOM loops. Notes on the usage by MC, commented output

!      MOM = 0 ; WRITE(OUTUNIT,*)
!      IF ( .NOT.DO_USER_STREAMS ) THEN
!        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING SLEAVE INPUT SET:'
!      ELSE
!        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING TWO SLEAVE INPUT SETS:'
!      ENDIF

      NMOMS = 2 * NSTREAMS - 1

      WRITE(OUTUNIT,*)
      DO IB=1,N_SZANGLES
        DO STRM=1,NSTREAMS ; DO S=1,NSTOKES ; DO MOM = 0, NMOMS
          WRITE(OUTUNIT,DWFR4) &
            'IB = ',IB,' STRM = ',STRM,' S = ',S,' MOM = ',MOM,&
            ' SLTERM_F_0(MOM,S,STRM,IB) = ',SLTERM_F_0(MOM, S,STRM,IB)
        END DO ; END DO ; END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,N_SZANGLES ; DO USTRM=1,N_USER_VZANGLES ; DO MOM = 0, NMOMS
          DO S=1,NSTOKES
             WRITE(OUTUNIT,DWFR4) &
               'IB = ',IB,' USTRM = ',USTRM,' S = ',S,' MOM = ',MOM,&
                  ' USER_SLTERM_F_0(MOM,S,USTRM,IB) = ',USER_SLTERM_F_0(MOM,S,USTRM,IB)
          END DO
        END DO ; END DO ; END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_SUP_SLEAVE_INPUT

!

      SUBROUTINE VLIDORT_WRITEFOURIER ( &
        RUN, FOURIER, N_USER_LEVELS, N_OUT_STREAMS, LOCAL_N_USERAZM, &
        NBEAMS, NSTOKES, N_DIRECTIONS, WHICH_DIRECTIONS, STOKES_F )

      USE VLIDORT_pars_m, Only : MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: RUN, FOURIER, N_USER_LEVELS, N_OUT_STREAMS, LOCAL_N_USERAZM, NBEAMS, &
                             NSTOKES, N_DIRECTIONS
      INTEGER, INTENT(IN) :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: STOKES_F ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS )

!  Local variables
!  ---------------

      INTEGER :: I, IB, IDIR, O, UT, W

!  Write Fourier output

      DO IDIR = 1, N_DIRECTIONS
        W = WHICH_DIRECTIONS(IDIR)
        DO IB = 1, NBEAMS
          DO UT = 1, N_USER_LEVELS
            DO I = 1, N_OUT_STREAMS
              DO O = 1, NSTOKES
                !WRITE(RUN,*)
                !WRITE(RUN,*) 'FOURIER = ',FOURIER,' IBEAM = ',IBEAM,' W = ',W,' UT = ',UT,' I = ',I,' O = ',O
                !WRITE(RUN,*) 'STOKES_F(UT,I,IBEAM,O,W) = ',STOKES_F(UT,I,IBEAM,O,W)
                WRITE(*,*)
                WRITE(*,*) 'FOURIER = ',FOURIER,' W = ',W,' IB = ',IB,' UT = ',UT,' I = ',I,' O = ',O
                WRITE(*,*) 'STOKES_F(UT,I,IB,O,W) = ',STOKES_F(UT,I,IB,O,W)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE VLIDORT_WRITEFOURIER

!

      SUBROUTINE VLIDORT_WRITERESULTS ( &
        RUN, DO_FULLRAD_MODE, DO_FOCORR_NADIR,  DO_FOCORR_OUTGOING, DO_DOUBLE_CONVTEST, &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,  &
        NSTOKES, VLIDORT_ACCURACY, SZANGLES, N_USER_RELAZMS,                            &
        USER_RELAZMS, N_USER_LEVELS, USER_LEVELS, HEIGHT_GRID, DELTAU_VERT_INPUT,       &
        DO_NO_AZIMUTH, NBEAMS, N_DIRECTIONS, WHICH_DIRECTIONS, N_OUT_STREAMS,           &
        OUT_ANGLES, PARTLAYERS_OUTFLAG, VZA_OFFSETS, TAUGRID_INPUT, DO_MULTIBEAM,       &
        VLIDORT_Out, MEAN_STOKES, FLUX_STOKES, FOURIER_SAVED )

!mick mod 1/5/2021 - modified to use VLIDORT output type structure "VLIDORT_Out" directly
!                    for STOKES output

!  Include file of dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS, MAXLAYERS, MAX_SZANGLES,           &
                                 MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_USER_LEVELS, &
                                 MAX_USER_OBSGEOMS, MAXSTOKES_SQ, MAXFOURIER, MAX_GEOMETRIES,            &
                                 FMT_SECTION, FMT_HEADING, FMT_CHAR, FMT_REAL, FMT_INTEGER, ONE, UPIDX, DNIDX

      USE VLIDORT_Outputs_def_m

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          RUN
      LOGICAL, INTENT (IN) ::          DO_FULLRAD_MODE
      LOGICAL, INTENT (IN) ::          DO_FOCORR_NADIR
      LOGICAL, INTENT (IN) ::          DO_FOCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_DOUBLE_CONVTEST
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY
      INTEGER, INTENT (IN) ::          NSTOKES
      DOUBLE PRECISION, INTENT (IN) :: VLIDORT_ACCURACY
      DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN)::  USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: USER_LEVELS ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )

      LOGICAL, INTENT (IN) ::          DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::          N_OUT_STREAMS

      DOUBLE PRECISION, INTENT (IN) :: OUT_ANGLES ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

      DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_MULTIBEAM  ( MAXBEAMS, 0:MAXFOURIER )

      !DOUBLE PRECISION, INTENT (IN) :: STOKES        ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      TYPE(VLIDORT_Main_Outputs), INTENT (IN) :: VLIDORT_Out

      DOUBLE PRECISION, INTENT (IN) :: MEAN_STOKES   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: FLUX_STOKES   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::          FOURIER_SAVED ( MAX_SZANGLES )

!  Local variables
!  ---------------

      INTEGER ::          I, UA, IB, V, UT, S, F, FMAX, NT, UTA
      INTEGER ::          LOCAL_NUSERAZMS, IDIR, WDIR

      DOUBLE PRECISION :: DT, USER_HEIGHTS(MAX_USER_LEVELS)
      DOUBLE PRECISION :: USER_OPDEPS (MAX_USER_LEVELS)

      CHARACTER (LEN=15), DIMENSION (MAXSTOKES) :: &
        STOKES_LABEL = (/ '  I component  ','  Q component  ', &
                          '  U component  ','  V component  ' /)

!  Beam Attenuation summary

      write(RUN,'(/a/a/)')' Results Output summary', &
                          ' ----------------------'

!  Removed for Version 2.8 7/6/16
      !IF ( DO_CLASSICAL_SOLUTION ) THEN
      !  write(RUN,'(a)') 'Classical solution of beam particular integral'
      !ELSE
      !  write(RUN,'(a)') 'Green function solution of beam particular integral'
      !ENDIF

      IF ( DO_PLANE_PARALLEL ) THEN
        write(RUN,'(a)') ' Plane parallel approximation to beam attenuation'
      ELSE
        write(RUN,'(a)') ' Pseudo-spherical approximation to beam attenuation'
        write(RUN,'(a)') '  * Average secant approximation was used'
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          write(RUN,'(a)') '  * Refractive geometry was used (Snells Law)'
        ELSE
          write(RUN,'(a)') '  * Straight line geometry was used (no refraction)'
        ENDIF
      ENDIF

      IF ( DO_FULLRAD_MODE ) THEN
        write(RUN,'(a)') ' Full Stokes calculation has been performed'
        IF ( DO_FOCORR_NADIR ) THEN
          write(RUN,'(a)') '  --> Nakajima-Tanaka TMS single scatter correction (nadir)'
        ENDIF
        IF ( DO_FOCORR_OUTGOING ) THEN
          write(RUN,'(a)') '  --> Nakajima-Tanaka TMS single scatter correction (outgoing)'
        ENDIF
        IF ( .NOT.DO_FOCORR_NADIR .AND. .NOT.DO_FOCORR_OUTGOING ) THEN
          write(RUN,'(a)') '  --> No single scatter correction has been applied'
        ENDIF
      ELSE
        write(RUN,'(a)') ' ONLY Multiple-scatter radiance calculation has been performed'
      ENDIF

      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        IF ( DO_DOUBLE_CONVTEST ) THEN
          write(RUN,'(/a)') ' Double convergence test was used for Fourier Azimuth Series'
        ELSE
          write(RUN,'(/a)') ' Single convergence test was used for Fourier Azimuth Series'
        ENDIF
        write(RUN,'(a,F10.7/)') '  --> Accuracy level was pre-set at : ',VLIDORT_ACCURACY

        write(RUN,'(a,I5)')' - Number Solar zenith angles : ',NBEAMS
        write(RUN,'(/a/)') '   SZA  | # Fourier | Fourier breakdown -->'
        FMAX = -100
        DO IB = 1, NBEAMS
          FMAX = MAX(FMAX,FOURIER_SAVED(IB))
        END DO
        DO IB = 1, NBEAMS
           write(RUN,'(1X,F7.2,4X,I3,6X,50(L1))') SZANGLES(IB),FOURIER_SAVED(IB), &
                                                  (DO_MULTIBEAM(IB,F),F=0,FMAX)
        END DO
      ELSE
        write(RUN,'(/A)') ' Azimuth independent output only (Fourier = 0)'
      ENDIF

!  Fix output

      DO UTA = 1, N_USER_LEVELS
        NT = INT(USER_LEVELS(UTA))
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          DT = USER_LEVELS(UTA) - DBLE(NT)
          USER_HEIGHTS(UTA) = HEIGHT_GRID(NT-1) * (ONE-DT) + HEIGHT_GRID(NT) * DT
          USER_OPDEPS(UTA)  = TAUGRID_INPUT(NT-1) + DELTAU_VERT_INPUT(NT)*DT
        ELSE
          USER_HEIGHTS(UTA) = HEIGHT_GRID(NT)
          USER_OPDEPS(UTA)  = TAUGRID_INPUT(NT)
        ENDIF
      ENDDO

!  Control point for avoiding intensity output

      IF ( DO_MVOUT_ONLY ) GO TO 400

!  Stokes vector output
!  --------------------

!  Local number of azimuths

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_NUSERAZMS = 1
      ELSE
        LOCAL_NUSERAZMS = N_USER_RELAZMS
      ENDIF

!  Overall header

      write(RUN,'(/a/a)')' Stokes vector output', &
                         ' --------------------'
!  Start beam loop

      DO IB = 1, NBEAMS

       write(RUN,*)
       write(RUN,FMT_REAL) ' * * RESULTS FOR SOLAR ZENITH ANGLE (degs) = ',SZANGLES(IB)

!  Start azimuth loop

       DO UA = 1, LOCAL_NUSERAZMS

!  Azimuth angle header

        IF (UA /= 1) write(RUN,*)
        IF ( DO_NO_AZIMUTH ) THEN
          write(RUN,FMT_CHAR) ' * * RESULTS FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
        ELSE
          write(RUN,FMT_REAL) ' * * RESULTS FOR RELATIVE AZIMUTH ANGLE (degs) = ',USER_RELAZMS(UA)
        ENDIF

        write(RUN,'(a)')' '
        write(RUN,FMT_INTEGER) ' Total number of output levels = ',N_USER_LEVELS
        write(RUN,FMT_INTEGER) ' Total number of output angles = ',N_OUT_STREAMS

!  Detailed output

        DO IDIR = 1, N_DIRECTIONS

          WDIR = WHICH_DIRECTIONS(IDIR)

!  Direction header

          IF (WDIR .EQ. UPIDX ) THEN
            write(RUN,'(/A)') '  --> Upwelling intensities all output levels and angles'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            write(RUN,'(/A)') '  --> Downwelling intensities all output levels and angles'
          ENDIF

!  Output loop

          DO UT = 1, N_USER_LEVELS
            write(RUN,'(/a,3(f10.5,2x))') ' Layer #, Height and Optical depth', &
                                           USER_LEVELS(UT), USER_HEIGHTS(UT), USER_OPDEPS(UT)
            write(RUN,'(a,4(a15))') ' Output angle |  ',(STOKES_LABEL(S),S=1,NSTOKES)
            DO I = 1, N_OUT_STREAMS
              V = VZA_OFFSETS(IB,I) + UA
              write(RUN,'(2x,F9.5,2x,5(1PE15.5))') OUT_ANGLES(I), (VLIDORT_Out%TS_STOKES(UT,V,S,WDIR),S=1,NSTOKES)
            ENDDO
          ENDDO

!  Direction and solar/azimuth loops - end

        ENDDO
       ENDDO
      ENDDO

!  Integrated value output
!  -----------------------

400   CONTINUE
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

       write(RUN,'(/a/a)')' Integrated value output', &
                          ' -----------------------'

       DO IB = 1, NBEAMS

        write(RUN,*)
        write(RUN,FMT_REAL) ' * * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)

!  Detailed output

        DO IDIR = 1, N_DIRECTIONS

          WDIR = WHICH_DIRECTIONS(IDIR)

!  Direction header

          IF (WDIR .EQ. UPIDX ) THEN
            write(RUN,'(/A/)') '  --> Upwelling mean values & fluxes, all output levels'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            write(RUN,'(/A/)') '  --> Downwelling mean values & fluxes, all output levels'
          ENDIF

!  Diffuse values

          write(RUN,'(a/)') '      ** Diffuse Mean-values for Stokes components  ----> '
          write(RUN,'(a13,3x,4(a15)/)')' Output levels', (STOKES_LABEL(S),S=1,NSTOKES)
          DO UT = 1, N_USER_LEVELS
            write(RUN,'(2x,F9.5,3x,1p4e15.5)') USER_LEVELS(UT), &
              (VLIDORT_Out%TS_MEANST_DIFFUSE(UT,IB,S,WDIR),S=1,NSTOKES)
          ENDDO

          write(RUN,'(/a/)') '      ** Diffuse Fluxes for Stokes components  ----> '
          write(RUN,'(a13,3x,4(a15)/)')' Output levels', (STOKES_LABEL(S),S=1,NSTOKES)
          DO UT = 1, N_USER_LEVELS
            write(RUN,'(2x,F9.5,3x,1p4e15.5)')USER_LEVELS(UT), &
              (VLIDORT_Out%TS_FLUX_DIFFUSE(UT,IB,S,WDIR),S=1,NSTOKES)
          ENDDO

!  Direct

          IF (WDIR .EQ. DNIDX ) THEN
            write(RUN,'(/a/)') '      ** Direct Downwelling Mean-values for Stokes components  ----> '
            write(RUN,'(a13,3x,4(a15)/)')' Output levels', (STOKES_LABEL(S),S=1,NSTOKES)
            DO UT = 1, N_USER_LEVELS
              write(RUN,'(2x,F9.5,3x,1p4e15.5)') USER_LEVELS(UT), &
                (VLIDORT_Out%TS_DNMEANST_DIRECT(UT,IB,S),S=1,NSTOKES)
            ENDDO

            write(RUN,'(/a/)') '      ** Direct Downwelling Fluxes for Stokes components  ----> '
            write(RUN,'(a13,3x,4(a15)/)')' Output levels', (STOKES_LABEL(S),S=1,NSTOKES)
            DO UT = 1, N_USER_LEVELS
              write(RUN,'(2x,F9.5,3x,1p4e15.5)')USER_LEVELS(UT), &
                (VLIDORT_Out%TS_DNFLUX_DIRECT(UT,IB,S),S=1,NSTOKES)
            ENDDO
          ENDIF

!  End direction loop

        ENDDO

       ENDDO
      ENDIF

      write(RUN,*)

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_WRITERESULTS

      END MODULE vlidort_writemodules_m
