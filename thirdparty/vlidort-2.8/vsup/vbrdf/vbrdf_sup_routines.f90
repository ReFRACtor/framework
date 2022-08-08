
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
! # General Subroutines in this Module                          #
! #                                                             #
! #              VBRDF_MAKER                                    #
! #              VBRDF_GCMCRI_MAKER                             #
! #              SCALING_FOURIER_ZERO (New, Version 2.7)        #
! #              VBRDF_FOURIER                                  #
! #                                                             #
! # New Cox-Munk Subroutine in this Module  (Version 2.7)       #
! #                                                             #
! #              VBRDF_NewCM_MAKER                              #
! #              VBRDF_NewGCM_MAKER                             #
! #                                                             #
! ###############################################################

!  Changes for Version 2.8
!  -----------------------

!  Allow for presence of Two new kernels. Set-up calls for them,

!  Version 2p8p1. 2019 Overhaul.
!  -----------------------------

!  In the Fourier routines--->

!     !@@ Rob Fix, 3/20/19. Bug: Lambertian case wrongly used NSTOKESSQ instead of 1
!     !@@ Rob Fix, 11/1/19. Bug: Cossin_Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)
!     !@@ Rob Fix, 11/1/19. Recoded the BRDF emissivity case
!     !@@ Rob Fix, 11/1/19. Used the Dot-Product command throughout. NSTOKESSQ Dropped.

!  Version 2p8p3.
!  --------------

!  1/31/21, Version 2.8.3. Analytical Model for Snow BRDF.
!     -- New VLIDORT BRDF Kernel. First introduced to VLIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  2/28/21, Version 3.8.3. Doublet geometry option added

!  7/28/21. Version 2.8.3. Some Changes
!    -- Half-Range azimuth integration introduced, now the default option. Controlled by parameter setting
!    -- Adjustment of +/- signs in Fourier-component sine-series settings. Validates against TestBed
!    -- MAXSTHALF_BRDF Dimensioning for CXE/SXE/BAX/EBRDFUNC/USER_EBRDFUNC arrays, replaced by MAXSTREAMS_BRDF

      MODULE vbrdf_sup_routines_m

      PRIVATE
      PUBLIC :: VBRDF_MAKER, &
                VBRDF_GCMCRI_MAKER, VBRDF_NewCM_MAKER, VBRDF_NewGCM_MAKER, &
                VBRDF_FOURIER, SCALING_FOURIER_ZERO

      CONTAINS

      SUBROUTINE VBRDF_MAKER &
         ( BRDF_VFUNCTION, BRDF_VFUNCTION_MSR,                               &
           DO_WSA_SCALING, DO_BSA_SCALING, DO_HALF_RANGE,                    & ! New line, Version 2.7
           DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          & ! Doublet-geometry flag added, Version 2.8.2
           DO_DBONLY, DO_MSRCORR, DO_MSRCORR_DBONLY, MSRCORR_ORDER,          &
           DO_USER_STREAMS, DO_SURFACE_EMISSION, N_MUQUAD, N_PHIQUAD,        &
           NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
           NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
           QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
           SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, BRDF_PARS,          &
           SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
           X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
           X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
           DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                          & ! output
           BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,               & ! output
           SCALING_BRDFUNC, SCALING_BRDFUNC_0  )                               ! output, New line, Version 2.7

!  1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag added to argument list

!  7/28/21. Some changes
!    -- Half-Range azimuth integration control input introduced (DO_HALF_RANGE)
!    -- MAXSTHALF_BRDF Dimensioning for CXE/SXE/EBRDFUNC/USER_EBRDFUNC arrays, replaced by MAXSTREAMS_BRDF

!  include file of dimensions and numbers
!    -- 7/28/21.  MAXSTHALF_BRDF Dimensioning dropped

      USE VLIDORT_PARS_m, only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                 MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS_SCALING, & 
                                 MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF,         &
                                 MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD, ZERO

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012.
!     Observation-Geometry input control.         (DO_USER_OBSGEOMS)
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!  Input arguments
!  ===============

!  BRDF functions (external calls)

      EXTERNAL         BRDF_VFUNCTION
      EXTERNAL         BRDF_VFUNCTION_MSR

!  White-sky and Black-sky albedo scaling flags. New Version 2.7

      LOGICAL ::          DO_WSA_SCALING
      LOGICAL ::          DO_BSA_SCALING

!  7/28/21. Half-range integration variable introduced
!   -- if set, azimuthal range is [0,pi], if not set the range is [pi,pi]

      LOGICAL ::          DO_HALF_RANGE

!  Solar sources + Observational Geometry flag
!    - 1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag 

      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL ::          DO_EXACT
!      LOGICAL ::          DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL ::          DO_MSRCORR
      LOGICAL ::          DO_MSRCORR_DBONLY   !  name change 9/28/14
      INTEGER ::          MSRCORR_ORDER
      INTEGER ::          N_MUQUAD, N_PHIQUAD

!  Local flags

      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_SURFACE_EMISSION

!  Number of Azimuth waudrature streams

      INTEGER ::          NSTREAMS_BRDF
      INTEGER ::          NBRDF_HALF

!  Local number of Stokes component matrix entries
!    value = 1 for most kernels, except GISS Cox-Munk

      INTEGER ::          NSTOKESSQ

!  Local number of Kernel parameters

      INTEGER ::          BRDF_NPARS

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  Discrete ordinates (local, for Albedo scaling). New Version 2.7

      INTEGER          :: SCALING_NSTREAMS
      DOUBLE PRECISION :: SCALING_QUAD_STREAMS(MAXSTREAMS_SCALING)
      DOUBLE PRECISION :: SCALING_QUAD_SINES  (MAXSTREAMS_SCALING)

!  Local parameter array

      DOUBLE PRECISION :: BRDF_PARS ( MAX_BRDF_PARAMETERS )

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )

!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: CXE_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SXE_BRDF ( MAXSTREAMS_BRDF )

!  Local arrays for MSR quadrature

      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  Output BRDF functions
!  =====================

!  At quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  At user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity
!   -- 7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: EBRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_EBRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Output for WSA/BSA scaling options. New, Version 2.7

      DOUBLE PRECISION :: SCALING_BRDFUNC &
          ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_0 &
          ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER ::            I, UI, J, K, KE, IB, ORDER, NSQ, NKE_RANGE
      INTEGER, PARAMETER :: LUM = 1
      INTEGER, PARAMETER :: LUA = 1
      DOUBLE PRECISION   :: KERNEL(16)

!  Local

      ORDER = MSRCORR_ORDER

!  Exact DB calculation
!  --------------------

!    !@@ Observational Geometry choice 12/31/12
!    !@@ Rob Fix 9/25/14. Always calculated for Solar Sources
!mick fix 9/19/2017 - initialize DBKERNEL_BRDFUNC

!  1/31/21. Version 2.8.3. Add doublet geometry option, in which the Exact-DB calculation
!   has the Azimuth angle coupled with the user zenith angle in a doublet.

      DBKERNEL_BRDFUNC = ZERO
      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
            !if(ORDER.eq.2)write(*,*)'Doing SZA/VZA = ',IB,UI
            IF ( DO_MSRCORR .or. DO_MSRCORR_DBONLY ) THEN
              CALL BRDF_VFUNCTION_MSR &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,  &
                 n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB),                     &
                 USER_STREAMS(IB), USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB),  &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,       &
                 DBKERNEL_BRDFUNC(1,LUM,LUA,IB) )
            ELSE
              CALL BRDF_VFUNCTION &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,   &
                 NSTOKESSQ, SZASURCOS(IB), SZASURSIN(IB),      &
                 USER_STREAMS(IB), USER_SINES(IB), PHIANG(IB), &
                 COSPHI(IB), SINPHI(IB),                       &
                 DBKERNEL_BRDFUNC(1,LUM,LUA,IB) )
            ENDIF
          ENDDO
        ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
          DO IB = 1, NBEAMS
             DO UI = 1, N_USER_STREAMS
                IF ( DO_MSRCORR .or. DO_MSRCORR_DBONLY ) THEN
                  CALL BRDF_VFUNCTION_MSR &
                      ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,          &
                        n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB),                     &
                        USER_STREAMS(UI), USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI),  &
                        X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,       &
                        DBKERNEL_BRDFUNC(1,UI,LUA,IB) )
                 ELSE
                   CALL BRDF_VFUNCTION &
                      ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,   &
                        NSTOKESSQ, SZASURCOS(IB), SZASURSIN(IB),      &
                        USER_STREAMS(UI), USER_SINES(UI), PHIANG(UI), &
                        COSPHI(UI), SINPHI(UI),                       &
                        DBKERNEL_BRDFUNC(1,UI,LUA,IB) )
                 ENDIF
              ENDDO
           ENDDO
        ELSE
          DO K = 1, N_USER_RELAZMS
             DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                   !if(ORDER.eq.2)write(*,*)'Doing SZA/VZA = ',IB,UI
                   IF ( DO_MSRCORR .or. DO_MSRCORR_DBONLY ) THEN
                     CALL BRDF_VFUNCTION_MSR &
                      ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,  &
                        n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB),                     &
                        USER_STREAMS(UI), USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),     &
                        X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,       &
                        DBKERNEL_BRDFUNC(1,UI,K,IB) )
                   ELSE
                     CALL BRDF_VFUNCTION &
                      ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,  &
                        NSTOKESSQ, SZASURCOS(IB), SZASURSIN(IB),     &
                        USER_STREAMS(UI), USER_SINES(UI), PHIANG(K), &
                        COSPHI(K), SINPHI(K),                        &
                        DBKERNEL_BRDFUNC(1,UI,K,IB) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
        END IF
      ENDIF

!  SCALING OPTIONS (New Section, Version 2.7)
!  ------------------------------------------

!  White-sky albedo, scaling. Only requires the (1,1) component
!     Use Local "Scaling_streams", both incident and outgoing

      IF ( DO_WSA_SCALING ) THEN
         NSQ = 1
         DO I = 1, SCALING_NSTREAMS
            DO J = 1, SCALING_NSTREAMS
               DO K = 1, NSTREAMS_BRDF
                  IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
                     CALL BRDF_VFUNCTION_MSR &
                        ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,      &
                          ORDER, NSQ, n_muquad, n_phiquad,                 &
                          SCALING_QUAD_STREAMS(J), SCALING_QUAD_SINES(J),  &
                          SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),  &         
                          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),               &
                          X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,      &
                          KERNEL )
                     SCALING_BRDFUNC(I,J,K) = KERNEL(1)
                  ELSE
                     CALL BRDF_VFUNCTION &
                        ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, NSQ, &
                          SCALING_QUAD_STREAMS(J), SCALING_QUAD_SINES(J),  &
                          SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),  &         
                          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),               &
                          KERNEL )
                     SCALING_BRDFUNC(I,J,K) = KERNEL(1)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Black-sky albedo, scaling. Only requires the (1,1) component
!     Use Local "Scaling_streams" for outgoing, solar beam for incoming (IB = 1)

      IF ( DO_BSA_SCALING .and. DO_SOLAR_SOURCES ) THEN
         IB = 1 ; NSQ = 1
         DO I = 1, SCALING_NSTREAMS
            DO K = 1, NSTREAMS_BRDF
               IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
                  CALL BRDF_VFUNCTION_MSR &
                      ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSQ,           &
                        n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB),                &
                        SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),                   &
                        X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                                &
                        X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
                        KERNEL )
                  SCALING_BRDFUNC_0(I,K) = KERNEL(1)
               ELSE
                  CALL BRDF_VFUNCTION &
                      ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,     &
                        NSQ, SZASURCOS(IB), SZASURSIN(IB),              &
                        SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I), &
                        X_BRDF(K), CX_BRDF(K), SX_BRDF(K),              &
                        KERNEL )
                  SCALING_BRDFUNC_0(I,K) = KERNEL(1)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  Return if the Direct-bounce BRDF is all that is required (scaled or not!)

      IF ( DO_DBONLY ) RETURN

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam
!    !@@  Solar Optionality. 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IB = 1, NBEAMS
          DO I = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
                CALL BRDF_VFUNCTION_MSR &
                ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,       &
                  n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I), &
                  QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                   &
                  X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,    &
                  BRDFUNC_0(1,I,IB,K) )
              ELSE
                CALL BRDF_VFUNCTION &
                ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, NSTOKESSQ, &
                  SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),         &
                  QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                  BRDFUNC_0(1,I,IB,K) )
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
           IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
            CALL BRDF_VFUNCTION_MSR &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,         &
                 n_muquad, n_phiquad, QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I), &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                     &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,      &
                 BRDFUNC(1,I,J,K) )
           ELSE
            CALL BRDF_VFUNCTION &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, NSTOKESSQ, &
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),       &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 BRDFUNC(1,I,J,K) )
           ENDIF
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions
!    - 7/28/21. Use the HALF_RANGE flag to set the number of azimuth points NKE_RANGE.

      IF ( DO_SURFACE_EMISSION ) THEN
        NKE_RANGE = NSTREAMS_BRDF ; if (.not. DO_HALF_RANGE ) NKE_RANGE = NBRDF_HALF
        DO I = 1, NSTREAMS
          DO KE = 1, NKE_RANGE
            DO K = 1, NSTREAMS_BRDF
             IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
              CALL BRDF_VFUNCTION_MSR &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,      &
                   n_muquad, n_phiquad, CXE_BRDF(KE), SXE_BRDF(KE),                   &
                   QUAD_STREAMS(I), QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,   &
                   EBRDFUNC(1,I,KE,K) )
             ELSE
              CALL BRDF_VFUNCTION &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, NSTOKESSQ,      &
                   CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I), QUAD_SINES(I), &
                   X_BRDF(K),  CX_BRDF(K), SX_BRDF(K),                         &
                   EBRDFUNC(1,I,KE,K) )
             ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam, Outgoing User-stream
!    !@@ Observational Geometry choice + Solar Optionality. 12/31/12

        IF (DO_SOLAR_SOURCES ) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            DO IB = 1, NBEAMS
              DO K = 1, NSTREAMS_BRDF
               IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
                 CALL BRDF_VFUNCTION_MSR &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,        &
                   n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB), &
                   USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                   &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                   USER_BRDFUNC_0(1,LUM,IB,K) )
               ELSE
                 CALL BRDF_VFUNCTION &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, NSTOKESSQ, &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),        &
                   USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),     &
                   USER_BRDFUNC_0(1,LUM,IB,K) )
               ENDIF
              ENDDO
            ENDDO
          ELSE
            DO IB = 1, NBEAMS
             DO UI = 1, N_USER_STREAMS
              DO K = 1, NSTREAMS_BRDF
               IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
                CALL BRDF_VFUNCTION_MSR &
                   ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,        &
                     n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI), &
                     USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                   &
                     X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                     USER_BRDFUNC_0(1,UI,IB,K) )
               ELSE
                CALL BRDF_VFUNCTION &
                   ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, NSTOKESSQ, &
                     SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),        &
                     USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),     &
                     USER_BRDFUNC_0(1,UI,IB,K) )
               ENDIF
              ENDDO
             ENDDO
            ENDDO
          ENDIF
        ENDIF

!  incident quadrature directions, Outgoing User-stream

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
             IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
              CALL BRDF_VFUNCTION_MSR &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,        &
                   n_muquad, n_phiquad, QUAD_STREAMS(J), QUAD_SINES(J),                 &
                   USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                   USER_BRDFUNC(1,UI,J,K) )
             ELSE
              CALL BRDF_VFUNCTION &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,  &
                   NSTOKESSQ, QUAD_STREAMS(J), QUAD_SINES(J),   &
                   USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), &
                   CX_BRDF(K), SX_BRDF(K),                      &
                   USER_BRDFUNC(1,UI,J,K) )
             ENDIF
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions
!  7/28/21. Use the HALF_RANGE flag to set the number of azimth points NKE_RANGE.

        IF ( DO_SURFACE_EMISSION ) THEN
          NKE_RANGE = NSTREAMS_BRDF ; if (.not. DO_HALF_RANGE ) NKE_RANGE = NBRDF_HALF
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NKE_RANGE
              DO K = 1, NSTREAMS_BRDF
               IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
                CALL BRDF_VFUNCTION_MSR &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, ORDER, NSTOKESSQ,        &
                   n_muquad, n_phiquad, CXE_BRDF(KE), SXE_BRDF(KE),                     &
                   USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                   USER_EBRDFUNC(1,UI,KE,K) )
               ELSE
                CALL BRDF_VFUNCTION &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,  &
                   NSTOKESSQ, CXE_BRDF(KE), SXE_BRDF(KE),       &
                   USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), &
                   CX_BRDF(K), SX_BRDF(K),                      &
                   USER_EBRDFUNC(1,UI,KE,K) )
               ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_MAKER

! 

      SUBROUTINE VBRDF_GCMCRI_MAKER &
         ( DO_WSA_SCALING, DO_BSA_SCALING, DO_HALF_RANGE,                     & ! New line, Version 2.7
           DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,           &
           DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION, DO_SHADOW_EFFECT, &
           DO_MSRCORR, DO_MSRCORR_DBONLY, MSRCORR_ORDER,                      &
           NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, NPARS,                       &
           NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                  &
           QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                &
           SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,                &
           SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,        & ! New line, Version 2.7
           X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, n_muquad, n_phiquad, &
           X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,   &
           DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                           & ! output
           BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,                & ! Output
           SCALING_BRDFUNC, SCALING_BRDFUNC_0  )                                ! output, New line, Version 2.7

!  1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag added to argument list

!  7/28/21. Some changes
!    -- Half-Range azimuth integration control input introduced (DO_HALF_RANGE)
!    -- MAXSTHALF_BRDF Dimensioning for CXE/SXE/EBRDFUNC/USER_EBRDFUNC arrays, replaced by MAXSTREAMS_BRDF

!  include file of dimensions and numbers
!    -- 7/28/21.  MAXSTHALF_BRDF Dimensioning dropped

      USE VLIDORT_PARS_m, only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                 MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS_SCALING, & 
                                 MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF,         &
                                 MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD, ZERO

      USE vbrdf_sup_kernels_m, only : GCMCRI_VFUNCTION_MSR, GCMCRI_VFUNCTION

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012.
!     Observation-Geometry input control.         (DO_USER_OBSGEOMS)
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!  Input arguments
!  ===============

!  White-sky and Black-sky albedo scaling flags. New Version 2.7

      LOGICAL ::          DO_WSA_SCALING
      LOGICAL ::          DO_BSA_SCALING

!  7/28/21. Half-range azimuth integration variable introduced
!   -- if set, azimuthal range is [0,pi], if not set the range is [pi,pi]

      LOGICAL ::          DO_HALF_RANGE

!   !@@ Solar sources + Observational Geometry flag !@@
!  1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag 

      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL ::          DO_EXACT
!      LOGICAL ::          DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  Local flags

      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_SURFACE_EMISSION
      LOGICAL ::          DO_SHADOW_EFFECT

!  Multiple reflectance correction for Glitter kernels

      LOGICAL ::          DO_MSRCORR
      INTEGER ::          MSRCORR_ORDER
      LOGICAL ::          DO_MSRCORR_DBONLY

!  Number of Azimuth waudrature streams

      INTEGER ::          NSTREAMS_BRDF
      INTEGER ::          NBRDF_HALF

!  Local number of Stokes component matrix entries
!    value = 1 for most kernels, except GISS Cox-Munk

      INTEGER ::          NSTOKESSQ

!  Local number of Kernel parameters

      INTEGER ::          NPARS

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  Discrete ordinates (local, for Albedo scaling). New Version 2.7

      INTEGER          :: SCALING_NSTREAMS
      DOUBLE PRECISION :: SCALING_QUAD_STREAMS(MAXSTREAMS_SCALING)
      DOUBLE PRECISION :: SCALING_QUAD_SINES  (MAXSTREAMS_SCALING)

!  Local parameter array

      DOUBLE PRECISION :: PARS ( MAX_BRDF_PARAMETERS )

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )

!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: CXE_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SXE_BRDF ( MAXSTREAMS_BRDF )

!  Local MSR quadrature control

      INTEGER          :: n_muquad, n_phiquad

!  Local arrays for MSR quadrature

      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity
!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: EBRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_EBRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Output for WSA/BSA scaling options. New, Version 2.7

      DOUBLE PRECISION :: SCALING_BRDFUNC &
          ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_0 &
          ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER ::          I, UI, J, K, KE, IB, ORDER, NSQ, NKE_RANGE
      LOGICAL ::          DOSHADOW
      INTEGER, PARAMETER :: LUM = 1
      INTEGER, PARAMETER :: LUA = 1
      DOUBLE PRECISION   :: KERNEL(16)

!  Local
!  -----

      DOSHADOW = DO_SHADOW_EFFECT
      ORDER = MSRCORR_ORDER

!  Direct-bounce calculation
!  -------------------------

!  9/27/14, remove DO_EXACT flag
!mick fix 9/19/2017 - initialize DBKERNEL_BRDFUNC

!  1/31/21. Version 2.8.3. Add doublet geometry option, in which the Exact-DB calculation
!   has the Azimuth angle coupled with the user zenith angle in a doublet.

      DBKERNEL_BRDFUNC(:,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = ZERO

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( .NOT. DO_USER_OBSGEOMS ) THEN
          IF ( DO_MSRCORR .or. DO_MSRCORR_DBONLY ) THEN
            DO K = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  CALL GCMCRI_VFUNCTION_MSR &
                  ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,        &
                    n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI), &
                    USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),                     &
                    X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                    DBKERNEL_BRDFUNC(1,UI,K,IB) )
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO K = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  CALL GCMCRI_VFUNCTION &
                  ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),        &
                    USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),       &
                    DBKERNEL_BRDFUNC(1,UI,K,IB) )
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
          IF ( DO_MSRCORR .or. DO_MSRCORR_DBONLY ) THEN
            DO IB = 1, NBEAMS
              DO UI = 1, N_USER_STREAMS
                CALL GCMCRI_VFUNCTION_MSR &
                  ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,        &
                    n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI), &
                    USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI),                  &
                    X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                    DBKERNEL_BRDFUNC(1,UI,LUA,IB) )
              ENDDO
            ENDDO
          ELSE
            DO IB = 1, NBEAMS
              DO UI = 1, N_USER_STREAMS
                CALL GCMCRI_VFUNCTION &
                  ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),        &
                    USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI),    &
                    DBKERNEL_BRDFUNC(1,UI,LUA,IB) )
              ENDDO
            ENDDO
          ENDIF
        ELSE
          IF ( DO_MSRCORR .or. DO_MSRCORR_DBONLY ) THEN
            DO IB = 1, NBEAMS
              CALL GCMCRI_VFUNCTION_MSR &
              ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,        &
                n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB), &
                USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB),                     &
                X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                DBKERNEL_BRDFUNC(1,LUM,LUA,IB) )
            ENDDO
          ELSE
            DO IB = 1, NBEAMS
              CALL GCMCRI_VFUNCTION &
              ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),        &
                USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB),       &
                DBKERNEL_BRDFUNC(1,LUM,LUA,IB) )
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!  SCALING OPTIONS (New Section, Version 2.7)
!  ------------------------------------------

!  White-sky albedo, scaling
!     Use Local "Scaling_streams", both incident and outgoing

      IF ( DO_WSA_SCALING ) THEN
         NSQ = 1
         DO I = 1, SCALING_NSTREAMS
            DO J = 1, SCALING_NSTREAMS
               DO K = 1, NSTREAMS_BRDF
                  IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
                     CALL GCMCRI_VFUNCTION_MSR &
                        ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER,         &
                          NSQ, DOSHADOW, n_muquad, n_phiquad,              &
                          SCALING_QUAD_STREAMS(J), SCALING_QUAD_SINES(J),  &
                          SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),  &         
                          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),               &
                          X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,      &
                          KERNEL )
                     SCALING_BRDFUNC(I,J,K)  = KERNEL(1)
                  ELSE
                     CALL GCMCRI_VFUNCTION &
                        ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                          SCALING_QUAD_STREAMS(J), SCALING_QUAD_SINES(J),        &
                          SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),        &         
                          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                     &
                          KERNEL )
                     SCALING_BRDFUNC(I,J,K)  = KERNEL(1)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Black-sky albedo, scaling
!     Use Local "Scaling_streams" for outgoing, solar beam for incoming (IB = 1)

      IF ( DO_BSA_SCALING .and. DO_SOLAR_SOURCES ) THEN
         IB = 1
         DO I = 1, SCALING_NSTREAMS
            DO K = 1, NSTREAMS_BRDF
               IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
                  CALL GCMCRI_VFUNCTION_MSR &
                     ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,     &
                       n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB),                &
                       SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),                   &
                       X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                                &
                       X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
                        KERNEL )
                     SCALING_BRDFUNC_0(I,K)  = KERNEL(1)
               ELSE
                  CALL GCMCRI_VFUNCTION &
                      ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ,    &
                        DOSHADOW, SZASURCOS(IB), SZASURSIN(IB),         &
                        SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I), &
                        X_BRDF(K), CX_BRDF(K), SX_BRDF(K),              &
                        KERNEL )
                     SCALING_BRDFUNC_0(I,K)  = KERNEL(1)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  Return if the Direct-bounce BRDF is all that is required (scaled or not!)

      IF ( DO_DBONLY ) RETURN

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam
!    !@@  Solar Optionality. 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
          DO IB = 1, NBEAMS
            DO I = 1, NSTREAMS
              DO K = 1, NSTREAMS_BRDF
                CALL GCMCRI_VFUNCTION_MSR &
                ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,       &
                  n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I), &
                  QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                   &
                  X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,    &
                  BRDFUNC_0(1,I,IB,K) )
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            DO I = 1, NSTREAMS
              DO K = 1, NSTREAMS_BRDF
                CALL GCMCRI_VFUNCTION &
                ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                  SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),         &
                  QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                  BRDFUNC_0(1,I,IB,K) )
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  incident quadrature directions

    IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL GCMCRI_VFUNCTION_MSR &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,         &
                 n_muquad, n_phiquad, QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I), &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                     &
                 X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,      &
                 BRDFUNC(1,I,J,K) )
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL GCMCRI_VFUNCTION &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),       &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 BRDFUNC(1,I,J,K) )
          ENDDO
        ENDDO
      ENDDO
    ENDIF

!  Emissivity (optional) - BRDF quadrature input directions
!  7/28/21. Use the HALF_RANGE flag to set the number of azimuth points NKE_RANGE.

    IF ( DO_SURFACE_EMISSION ) THEN
      NKE_RANGE = NSTREAMS_BRDF ; if (.not. DO_HALF_RANGE ) NKE_RANGE = NBRDF_HALF
      IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
        DO I = 1, NSTREAMS
          DO KE = 1, NKE_RANGE
            DO K = 1, NSTREAMS_BRDF
              CALL GCMCRI_VFUNCTION_MSR &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,     &
                   n_muquad, n_phiquad, CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I), &
                   QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                 &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
                   EBRDFUNC(1,I,KE,K) )
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          DO KE = 1, NKE_RANGE
            DO K = 1, NSTREAMS_BRDF
              CALL GCMCRI_VFUNCTION &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                   CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I),           &
                   QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                   EBRDFUNC(1,I,KE,K) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam, Outgoing User-stream
!    !@@ Observational Geometry choice + Solar Optionality. 12/31/12

        IF (DO_SOLAR_SOURCES ) THEN
          IF (.NOT. DO_USER_OBSGEOMS ) THEN
            IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
             DO IB = 1, NBEAMS
              DO UI = 1, N_USER_STREAMS
               DO K = 1, NSTREAMS_BRDF
                 CALL GCMCRI_VFUNCTION_MSR &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,        &
                   n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI), &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                   &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                   USER_BRDFUNC_0(1,UI,IB,K) )
               ENDDO
              ENDDO
             ENDDO
            ELSE
             DO IB = 1, NBEAMS
              DO UI = 1, N_USER_STREAMS
               DO K = 1, NSTREAMS_BRDF
                  CALL GCMCRI_VFUNCTION &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),        &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),     &
                   USER_BRDFUNC_0(1,UI,IB,K) )
               ENDDO
              ENDDO
             ENDDO
            ENDIF
          ELSE
            IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
             DO IB = 1, NBEAMS
               DO K = 1, NSTREAMS_BRDF
                 CALL GCMCRI_VFUNCTION_MSR &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,        &
                   n_muquad, n_phiquad, SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB), &
                   USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                   &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,     &
                   USER_BRDFUNC_0(1,LUM,IB,K) )
               ENDDO
             ENDDO
            ELSE
             DO IB = 1, NBEAMS
               DO K = 1, NSTREAMS_BRDF
                 CALL GCMCRI_VFUNCTION &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),        &
                   USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),     &
                   USER_BRDFUNC_0(1,LUM,IB,K) )
               ENDDO
             ENDDO
            ENDIF
          ENDIF
        ENDIF

!  incident quadrature directions

       IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL GCMCRI_VFUNCTION_MSR &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,          &
                   n_muquad, n_phiquad, QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI), &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                     &
                   X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,       &
                   USER_BRDFUNC(1,UI,J,K) )
            ENDDO
          ENDDO
        ENDDO
       ELSE
        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL GCMCRI_VFUNCTION &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                   QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),      &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),     &
                   USER_BRDFUNC(1,UI,J,K) )
            ENDDO
          ENDDO
        ENDDO
       ENDIF

!  Emissivity (optional) - BRDF quadrature input directions
!  7/28/21. Use the HALF_RANGE flag to set the number of azimuth points NKE_RANGE.

       IF ( DO_SURFACE_EMISSION ) THEN
         NKE_RANGE = NSTREAMS_BRDF ; if (.not. DO_HALF_RANGE ) NKE_RANGE = NBRDF_HALF
         IF ( DO_MSRCORR .and..not.DO_MSRCORR_DBONLY ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NKE_RANGE
              DO K = 1, NSTREAMS_BRDF
                CALL GCMCRI_VFUNCTION_MSR &
                   ( MAX_BRDF_PARAMETERS, NPARS, PARS, ORDER, NSTOKESSQ, DOSHADOW,      &
                     n_muquad, n_phiquad, CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI), &
                     USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                 &
                     X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,   &
                     USER_EBRDFUNC(1,UI,KE,K) )
              ENDDO
            ENDDO
          ENDDO
         ELSE
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NKE_RANGE
              DO K = 1, NSTREAMS_BRDF
                CALL GCMCRI_VFUNCTION &
                   ( MAX_BRDF_PARAMETERS, NPARS, PARS, NSTOKESSQ, DOSHADOW, &
                     CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI),          &
                     USER_SINES(UI), X_BRDF(K),  CX_BRDF(K), SX_BRDF(K),    &
                     USER_EBRDFUNC(1,UI,KE,K) )
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_GCMCRI_MAKER

!

      SUBROUTINE VBRDF_NewCM_MAKER &
             ( DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR,              &
               Refrac_R, Refrac_I, WC_Reflectance, WC_Lambertian,                 &
               DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, DO_DBONLY, &
               NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,   &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                      &
               X_BRDF, CX_BRDF, SX_BRDF,                                          &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                           & ! output
               BRDFUNC_0, USER_BRDFUNC_0  )                                         ! output

!  1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag added to argument list

!  include file of dimensions and numbers

      USE VLIDORT_PARS_m,      only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                      MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS_BRDF, &
                                      ZERO, ONE, DEG_TO_RAD

      USE vbrdf_sup_kernels_m, only : VBRDF_Generalized_Glint

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Input arguments
!  ===============

!  NewCM Glitter options (bypasses the usual Kernel system)
!  -------------------------------------------------------

!  Flags for glint shadowing, Facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FacetIsotropy

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      DOUBLE PRECISION:: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Refractive Index

      DOUBLE PRECISION :: Refrac_R, Refrac_I

!  Whitecap correction (Zero if not flagged)

      DOUBLE PRECISION :: WC_Reflectance, WC_Lambertian

!  Local flags
!    1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag 

      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL ::          DO_EXACT
!      LOGICAL ::          DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  Number of Azimuth quadrature streams

      INTEGER ::          NSTREAMS_BRDF

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  local variables
!  ---------------

      LOGICAL          :: DO_COEFFS, Local_Isotropy
      INTEGER          :: I, UI, J, K, IB
      DOUBLE PRECISION :: PHI_W(MAXBEAMS), CPHI_W(MAXBEAMS), SPHI_W(MAXBEAMS)
      DOUBLE PRECISION :: SUNGLINT_COEFFS(7), WC_correction, KERNEL

      INTEGER, PARAMETER :: LUM = 1
      INTEGER, PARAMETER :: LUA = 1

!   Wind-direction and coefficient set-up

      DO_COEFFS = .true.
      PHI_W = zero ; CPHI_W = one ; SPHI_W = zero
      Local_Isotropy = DO_FacetIsotropy 
      if ( .not.Local_Isotropy ) then
         DO IB = 1, nbeams
            PHI_W(IB)  = WINDDIR(IB)
            CPHI_W(IB) = cos(WINDDIR(IB) * deg_to_rad) 
            SPHI_W(IB) = sin(WINDDIR(IB) * deg_to_rad)
         ENDDO
      endif

!  Whitecap correction to glint

      WC_correction = one - WC_Lambertian

!  Direct Bounce calculation
!  -------------------------
!mick fix 9/19/2017  - initialize DBKERNEL_BRDFUNC
!mick mod 11/14/2019 - trimmed dimensions initialized

!  1/31/21. Version 2.8.3. Add doublet geometry option, in which the Exact-DB calculation
!   has the Azimuth angle coupled with the user zenith angle in a doublet.

      DBKERNEL_BRDFUNC(:,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = ZERO

      IF ( .NOT. DO_USER_OBSGEOMS ) THEN
         DO K = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO  UI = 1, N_USER_STREAMS
                  CALL VBRDF_Generalized_Glint &
                    ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,       &
                      REFRAC_R, REFRAC_I, WINDSPEED,                   &
                      PHI_W(IB), CPHI_W(IB), SPHI_W(IB),               &
                      SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),  &
                      USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K), &
                      SUNGLINT_COEFFS, KERNEL )
                  DBKERNEL_BRDFUNC(1,UI,K,IB) = WC_Reflectance + WC_correction * KERNEL
               ENDDO
            ENDDO
         ENDDO
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
         DO IB = 1, NBEAMS
            DO  UI = 1, N_USER_STREAMS
                  CALL VBRDF_Generalized_Glint &
                    ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
                      REFRAC_R, REFRAC_I, WINDSPEED,                      &
                      PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                  &
                      SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),     &
                      USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI), &
                      SUNGLINT_COEFFS, KERNEL )
               DBKERNEL_BRDFUNC(1,UI,LUA,IB) = WC_Reflectance + WC_correction * KERNEL
            ENDDO
         ENDDO
      ELSE
         DO IB = 1, NBEAMS
            CALL VBRDF_Generalized_Glint &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
               REFRAC_R, REFRAC_I, WINDSPEED,                      &
               PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                  &
               SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),     &
               USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB), &
               SUNGLINT_COEFFS, KERNEL )
            DBKERNEL_BRDFUNC(1,LUM,LUA,IB) = WC_Reflectance + WC_correction * KERNEL
         ENDDO
      ENDIF

!      pause'after direct bounce'

!  Return if this is all you require

      IF ( DO_DBONLY ) RETURN

!  Incident Solar beam
!  ===================
!mick fix 11/14/2019  - initialize BRDFUNC_0, USER_BRDFUNC_0

!  Quadrature outgoing directions

      BRDFUNC_0(:,1:NSTREAMS,1:NBEAMS,1:NSTREAMS_BRDF) = ZERO

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL VBRDF_Generalized_Glint &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               REFRAC_R, REFRAC_I, WINDSPEED,                    &
               PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                &
               SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),    &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, KERNEL )
            BRDFUNC_0(1,I,IB,K) = WC_Reflectance + WC_correction * KERNEL
          ENDDO
        ENDDO
      ENDDO

!  User-streams outgoing directions
!   This is the "Truncated" Direct Bounce calculation

      IF ( DO_USER_STREAMS ) THEN

        USER_BRDFUNC_0(:,1:N_USER_STREAMS,1:NBEAMS,1:NSTREAMS_BRDF) = ZERO

        IF (.NOT. DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              DO K = 1, NSTREAMS_BRDF
                CALL VBRDF_Generalized_Glint &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                   REFRAC_R, REFRAC_I, WINDSPEED,                     &
                   PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                 &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),    &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   SUNGLINT_COEFFS, KERNEL )
                USER_BRDFUNC_0(1,UI,IB,K) = WC_Reflectance + WC_correction * KERNEL
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            DO K = 1, NSTREAMS_BRDF
              CALL VBRDF_Generalized_Glint &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 REFRAC_R, REFRAC_I, WINDSPEED,                     &
                 PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                 &
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),    &
                 USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, KERNEL )
              USER_BRDFUNC_0(1,LUM,IB,K) = WC_Reflectance + WC_correction * KERNEL
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Incident quadrature directions (MULTIPLE SCATTERING)
!  ==============================
!mick fix 11/14/2019  - initialize BRDFUNC, USER_BRDFUNC

!   Can only be treated with 1 Wind direction.....
!     if ( NBEAMS > 1) MUST assume local Facet Isotropy
!            --> set up Local Wind-direction and re-set coefficients flag.
!     if ( NBAMS = 1 ) use the first wind direction, no need to re-calculate coefficients

      if ( NBEAMS .gt. 1 ) then
         local_Isotropy = .true.
!         local_Isotropy = .false.  ! Bug 9/27/14, fixed
         PHI_W      = zero 
         CPHI_W     = one
         SPHI_W     = zero
         DO_COEFFS  = .true.
      endif
 
!  Outgoing quadrature directions

      BRDFUNC(:,1:NSTREAMS,1:NSTREAMS,1:NSTREAMS_BRDF) = ZERO

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL VBRDF_Generalized_Glint &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               REFRAC_R, REFRAC_I, WINDSPEED,                    &
               PHI_W(1), CPHI_W(1), SPHI_W(1),                   &
               QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),  &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, KERNEL )
            BRDFUNC(1,I,J,K) = WC_Reflectance + WC_correction * KERNEL
          ENDDO
        ENDDO
      ENDDO

!  User stream outgoing directions

      IF ( DO_USER_STREAMS ) THEN

        USER_BRDFUNC(:,1:N_USER_STREAMS,1:NSTREAMS,1:NSTREAMS_BRDF) = ZERO

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL VBRDF_Generalized_Glint &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 REFRAC_R, REFRAC_I, WINDSPEED,                     &
                 PHI_W(1), CPHI_W(1), SPHI_W(1),                    &
                 QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),  &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, KERNEL )
              USER_BRDFUNC(1,UI,J,K) = WC_Reflectance + WC_correction * KERNEL
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_NewCM_MAKER

!

      SUBROUTINE VBRDF_NewGCM_MAKER &
             ( DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR,                    &
               Refrac_R, Refrac_I, WC_Reflectance, WC_Lambertian,                       &
               DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, DO_DBONLY, NSSQ, &
               NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,              &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                    &
               X_BRDF, CX_BRDF, SX_BRDF,                                        &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                         & ! output
               BRDFUNC_0, USER_BRDFUNC_0  )                                       ! output

!  1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag added to argument list

!  include file of dimensions and numbers

      USE VLIDORT_PARS_m,      only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                      MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS_BRDF,    &
                                      zero, one, DEG_TO_RAD

      USE vbrdf_sup_kernels_m, only : VBRDF_Generalized_Glint_GCM

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Input arguments
!  ===============

!  NewCM Glitter options (bypasses the usual Kernel system)
!  -------------------------------------------------------

!  Flags for glint shadowing, Facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FacetIsotropy

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      DOUBLE PRECISION:: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Refractive Index

      DOUBLE PRECISION :: Refrac_R, Refrac_I

!  Whitecap correction (Zero if not flagged)

      DOUBLE PRECISION :: WC_Reflectance, WC_Lambertian

!  Local flags
!    1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag 

      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL ::          DO_EXACT
!      LOGICAL ::          DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  Number of stokes_Sq

      INTEGER ::          NSSQ
 
!  Number of Azimuth quadrature streams

      INTEGER ::          NSTREAMS_BRDF

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0 &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  local variables
!  ---------------

      LOGICAL          :: DO_COEFFS, Local_Isotropy
      INTEGER          :: I, UI, J, K, IB, NSSQ_MASK(4), NSTOKES, NSM
      DOUBLE PRECISION :: PHI_W(MAXBEAMS), CPHI_W(MAXBEAMS), SPHI_W(MAXBEAMS)
      DOUBLE PRECISION :: SUNGLINT_COEFFS(7), WC_correction, VKERNEL(16)

      INTEGER, PARAMETER :: LUM = 1
      INTEGER, PARAMETER :: LUA = 1

!   Define "NSTOKES squared" mask and associated upper limit NSM (redefines number
!   of elements from NewGCM kernel output arrays to pass to local output arrays
!   depending on the value of NSTOKES to fix some element-passing bugs)
!mick fix 11/14/2019 - added this section

      NSSQ_MASK = (/ 1,6,11,16 /)
      NSTOKES = INT(SQRT(REAL(NSSQ)))
      NSM = NSSQ_MASK(NSTOKES)

!   Wind-direction and coefficient set-up

      DO_COEFFS = .true.
      PHI_W = zero ; CPHI_W = one ; SPHI_W = zero
      Local_Isotropy = DO_FacetIsotropy 
      if ( .not.Local_Isotropy ) then
         DO IB = 1, nbeams
            PHI_W(IB)  = WINDDIR(IB)
            CPHI_W(IB) = cos(WINDDIR(IB) * deg_to_rad) 
            SPHI_W(IB) = sin(WINDDIR(IB) * deg_to_rad)
         ENDDO
      endif

!  Whitecap correction to glint

      WC_correction = one - WC_Lambertian

!  Direct Bounce calculation
!  -------------------------
!mick fix 11/14/2019  - initialize DBKERNEL_BRDFUNC
!                     - changed dimensional extent of elements being passed from subroutine
!                       output arrays to calling routine output arrays from NSSQ to NSM

!  1/31/21. Version 2.8.3. Add doublet geometry option, in which the Exact-DB calculation
!   has the Azimuth angle coupled with the user zenith angle in a doublet.

      DBKERNEL_BRDFUNC(:,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = ZERO

      IF ( .NOT. DO_USER_OBSGEOMS ) THEN
         DO K = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO  UI = 1, N_USER_STREAMS
                  CALL VBRDF_Generalized_Glint_GCM &
                    ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,       &
                      REFRAC_R, REFRAC_I, WINDSPEED,                   &
                      PHI_W(IB), CPHI_W(IB), SPHI_W(IB),               &
                      SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),  &
                      USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K), &
                      SUNGLINT_COEFFS, VKERNEL )
                  DBKERNEL_BRDFUNC(1:NSM,UI,K,IB) = WC_correction * VKERNEL(1:NSM)
                  DBKERNEL_BRDFUNC(1,UI,K,IB)      = WC_Reflectance + DBKERNEL_BRDFUNC(1,UI,K,IB)
               ENDDO
            ENDDO
         ENDDO
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
         DO IB = 1, NBEAMS
            DO  UI = 1, N_USER_STREAMS
               CALL VBRDF_Generalized_Glint_GCM &
                    ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
                      REFRAC_R, REFRAC_I, WINDSPEED,                      &
                      PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                  &
                      SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),     &
                      USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI), &
                      SUNGLINT_COEFFS, VKERNEL )
               DBKERNEL_BRDFUNC(1:NSM,UI,LUA,IB) = WC_correction * VKERNEL(1:NSM)
               DBKERNEL_BRDFUNC(1,UI,LUA,IB)     = WC_Reflectance + DBKERNEL_BRDFUNC(1,UI,LUA,IB)
            ENDDO
         ENDDO
      ELSE
         DO IB = 1, NBEAMS
            CALL VBRDF_Generalized_Glint_GCM &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
               REFRAC_R, REFRAC_I, WINDSPEED,                      &
               PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                  &
               SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),     &
               USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB), &
               SUNGLINT_COEFFS, VKERNEL )
            DBKERNEL_BRDFUNC(1:NSM,LUM,LUA,IB) = WC_correction * VKERNEL(1:NSM)
            DBKERNEL_BRDFUNC(1,LUM,LUA,IB)      = WC_Reflectance + DBKERNEL_BRDFUNC(1,LUM,LUA,IB)
         ENDDO
      ENDIF

!      pause'after direct bounce'

!  Return if this is all you require

      IF ( DO_DBONLY ) RETURN

!  Incident Solar beam
!  ===================
!mick fix 11/14/2019  - initialize BRDFUNC_0, USER_BRDFUNC_0
!                     - changed dimensional extent of elements being passed from subroutine
!                       output arrays to calling routine output arrays from NSSQ to NSM

!  Quadrature outgoing directions

      BRDFUNC_0(:,1:NSTREAMS,1:NBEAMS,1:NSTREAMS_BRDF) = ZERO

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL VBRDF_Generalized_Glint_GCM &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               REFRAC_R, REFRAC_I, WINDSPEED,                    &
               PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                &
               SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),    &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, VKERNEL )
            BRDFUNC_0(1:NSM,I,IB,K) = WC_correction * VKERNEL(1:NSM)
            BRDFUNC_0(1,I,IB,K)      = WC_Reflectance + BRDFUNC_0(1,I,IB,K)
          ENDDO
        ENDDO
      ENDDO

!  User-streams outgoing directions
!   This is the "Truncated" Direct Bounce calculation

      IF ( DO_USER_STREAMS ) THEN

        USER_BRDFUNC_0(:,1:N_USER_STREAMS,1:NBEAMS,1:NSTREAMS_BRDF) = ZERO

        IF (.NOT. DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              DO K = 1, NSTREAMS_BRDF
                CALL VBRDF_Generalized_Glint_GCM &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                   REFRAC_R, REFRAC_I, WINDSPEED,                     &
                   PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                 &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),    &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   SUNGLINT_COEFFS, VKERNEL )
                USER_BRDFUNC_0(1:NSM,UI,IB,K) = WC_correction * VKERNEL(1:NSM)
                USER_BRDFUNC_0(1,UI,IB,K)      = WC_Reflectance + USER_BRDFUNC_0(1,UI,IB,K)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            DO K = 1, NSTREAMS_BRDF
              CALL VBRDF_Generalized_Glint_GCM &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 REFRAC_R, REFRAC_I, WINDSPEED,                     &
                 PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                 &
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),    &
                 USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, VKERNEL )
              USER_BRDFUNC_0(1:NSM,LUM,IB,K) = WC_correction * VKERNEL(1:NSM)
              USER_BRDFUNC_0(1,LUM,IB,K)      = WC_Reflectance + USER_BRDFUNC_0(1,LUM,IB,K)
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Incident quadrature directions (MULTIPLE SCATTERING)
!  ==============================
!mick fix 11/14/2019  - initialize BRDFUNC, USER_BRDFUNC
!                     - changed dimensional extent of elements being passed from subroutine
!                       output arrays to calling routine output arrays from NSSQ to NSM

!   Can only be treated with 1 Wind direction.....
!     if ( NBEAMS > 1) MUST assume local Facet Isotropy
!            --> set up Local Wind-direction and re-set coefficients flag.
!     if ( NBAMS = 1 ) use the first wind direction, no need to re-calculate coefficients

      if ( NBEAMS .gt. 1 ) then
         local_Isotropy = .true.
!         local_Isotropy = .false.  ! Bug 9/27/14, fixed
         PHI_W      = zero 
         CPHI_W     = one
         SPHI_W     = zero
         DO_COEFFS  = .true.
      endif
 
!  Outgoing quadrature directions

      BRDFUNC(:,1:NSTREAMS,1:NSTREAMS,1:NSTREAMS_BRDF) = ZERO

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL VBRDF_Generalized_Glint_GCM &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               REFRAC_R, REFRAC_I, WINDSPEED,                    &
               PHI_W(1), CPHI_W(1), SPHI_W(1),                   &
               QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),  &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, VKERNEL )
            BRDFUNC(1:NSM,I,J,K) = WC_correction * VKERNEL(1:NSM)
            BRDFUNC(1,I,J,K)      = WC_Reflectance + BRDFUNC(1,I,J,K)
          ENDDO
        ENDDO
      ENDDO

!  User stream outgoing directions

      IF ( DO_USER_STREAMS ) THEN

        USER_BRDFUNC(:,1:N_USER_STREAMS,1:NSTREAMS,1:NSTREAMS_BRDF) = ZERO

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL VBRDF_Generalized_Glint_GCM &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 REFRAC_R, REFRAC_I, WINDSPEED,                     &
                 PHI_W(1), CPHI_W(1), SPHI_W(1),                    &
                 QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),  &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, VKERNEL )
              USER_BRDFUNC(1:NSM,UI,J,K) = WC_correction * VKERNEL(1:NSM)
              USER_BRDFUNC(1,UI,J,K)      = WC_Reflectance + USER_BRDFUNC(1,UI,J,K)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_NewGCM_MAKER

!

      SUBROUTINE SCALING_FOURIER_ZERO &
            ( DO_LOCAL_WSA, DO_LOCAL_BSA, LAMBERTIAN_FLAG,            &
              DO_HALF_RANGE, SCALING_NSTREAMS, NSTREAMS_BRDF, A_BRDF, &
              SCALING_BRDFUNC, SCALING_BRDFUNC_0,                     &
              SCALING_BRDF_F, SCALING_BRDF_F_0 )

!  Version 2p8p1. Patch Overhaul. Code is more compact.
!     !@@ Rob Fix, 11/1/19. Used the Dot-Product command throughout.

!   7/28/21. Have to introduce the Half-range flag here.

!  include file of dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXSTREAMS_BRDF, MAXSTREAMS_SCALING, ZERO, ONE, TWO, HALF

      IMPLICIT NONE

!  Prepares Fourier component of the bidirectional reflectance functions, (1,1) element for aling albedos

!  Input arguments
!  ===============

!  Local flags

      LOGICAL ::          DO_LOCAL_WSA, DO_LOCAL_BSA

!  Control

      LOGICAL ::          LAMBERTIAN_FLAG

!  7/28/21. Introduce Half-range flag

      LOGICAL ::          DO_HALF_RANGE

!  Local numbers

      INTEGER ::          SCALING_NSTREAMS, NSTREAMS_BRDF

!  Azimuth weights

      DOUBLE PRECISION :: A_BRDF ( MAXSTREAMS_BRDF )

!  Input for WSA/BSA scaling options. New, Version 2.7

      DOUBLE PRECISION :: SCALING_BRDFUNC   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: SCALING_BRDF_F   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING )
      DOUBLE PRECISION :: SCALING_BRDF_F_0 ( MAXSTREAMS_SCALING   )

!  local variables
!  ===============

      INTEGER ::          I, J 
      DOUBLE PRECISION :: HELP

!  Zeroing

      SCALING_BRDF_F        = ZERO
      SCALING_BRDF_F_0      = ZERO

!  surface factor
!   -- 7/28/21. set to half the value if you are using the full-range azimuth integration

      if ( DO_HALF_RANGE ) THEN
         HELP = ONE
      else
         HELP = HALF
      endif

!  Quadrature outgoing directions
!  ------------------------------

!  BSA: Incident Solar beam

      IF ( DO_LOCAL_BSA ) THEN
         IF ( .NOT. LAMBERTIAN_FLAG ) THEN
            DO I = 1, SCALING_NSTREAMS
               SCALING_BRDF_F_0(I) = HELP * DOT_PRODUCT(SCALING_BRDFUNC_0(I,1:NSTREAMS_BRDF),A_BRDF(1:NSTREAMS_BRDF))
            ENDDO
         ELSE
            SCALING_BRDF_F_0 = ONE
         ENDIF
      ENDIF

!  WSA: incident quadrature directions

      if ( DO_LOCAL_WSA ) THEN
         IF ( .NOT. LAMBERTIAN_FLAG ) THEN
            DO I = 1, SCALING_NSTREAMS
               DO J = 1, SCALING_NSTREAMS
                  SCALING_BRDF_F(I,J) = HELP * DOT_PRODUCT(SCALING_BRDFUNC(I,J,1:NSTREAMS_BRDF),A_BRDF(1:NSTREAMS_BRDF))
               ENDDO
            ENDDO
         ELSE 
            SCALING_BRDF_F = ONE
         ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SCALING_FOURIER_ZERO

!  

      SUBROUTINE VBRDF_FOURIER &
         ( DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_USER_STREAMS, DO_SURFEMISS, LAMBERTIAN_FLAG,     & ! Control flags
           DO_HALF_RANGE, M, NSTOKES, NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Numbers
           DELFAC, FACTOR, BRDF_COSAZMFAC, BRDF_SINAZMFAC, A_BRDF, BAX_BRDF,                       & ! Surface/Azimuth factors
           BRDFUNC, USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,              & ! Input BRDF matrices
           L_BRDF_F, L_BRDF_F_0, L_USER_BRDF_F, L_USER_BRDF_F_0, L_EMISSIVITY, L_USER_EMISSIVITY )   ! Output Fourier Components

!  Version 2p8p1. Overhaul. 2019. Code is more compact.
!     !@@ Rob Fix, 3/20/19. Bug: Lambertian case wrongly used NSTOKESSQ instead of 1
!     !@@ Rob Fix, 11/1/19. Bug: Cossin_Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)
!     !@@ Rob Fix, 11/1/19. Recoded the BRDF emissivity case
!     !@@ Rob Fix, 11/1/19. Used the Dot-Product command throughout. NSTOKESSQ Dropped.

!  7/28/21. Version 2.8.3. Some Changes
!    -- Half-Range azimuth integration introduced, now the default option. Controlled by parameter setting
!    -- Adjustment of +/- signs in Fourier-component sine-series settings. Validates against TestBed
!    -- MAXSTHALF_BRDF Dimensioning for CXE/SXE/EBRDFUNC/USER_EBRDFUNC/D_EBRDFUNC/D_USER_EBRDFUNC arrays, replaced by MAXSTREAMS_BRDF

!  include file of dimensions and numbers
!    -- MAXSTHALF_BRDF Dimensioning dropped

      USE VLIDORT_PARS_m, Only : MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXBEAMS, &
                                 MAXSTOKES, MAXSTREAMS_BRDF, ZERO, HALF, ONE, TWO

      IMPLICIT NONE

!  Prepares Fourier component of the bidirectional reflectance functions

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012.
!     Observation-Geometry input control.         (DO_USER_OBSGEOMS)
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)

!  Input arguments
!  ===============

!  Control flags
!   !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          LAMBERTIAN_FLAG
      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_SURFEMISS

!  7/28/21. Introduce Half-Range flag for control of azimuth integration

      LOGICAL ::          DO_HALF_RANGE

!  Local numbers

      INTEGER ::          M, NSTOKES, NSTREAMS
      INTEGER ::          NBEAMS, N_USER_STREAMS
      INTEGER ::          NSTREAMS_BRDF, NBRDF_HALF

!  Surface factors

      DOUBLE PRECISION :: DELFAC, FACTOR

!  Azimuth cosines/sines and weights

      DOUBLE PRECISION :: BRDF_COSAZMFAC ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDF_SINAZMFAC ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: A_BRDF         ( MAXSTREAMS_BRDF )

!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: BAX_BRDF       ( MAXSTREAMS_BRDF  )

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Values for Emissivity
!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: EBRDFUNC      ( MAXSTOKES_SQ, MAXSTREAMS,       MAXSTREAMS_BRDF, MAXSTREAMS_BRDF)
      DOUBLE PRECISION :: USER_EBRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF)

!  Output: Local kernel Fourier components
!  =======================================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: L_BRDF_F   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: L_BRDF_F_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      DOUBLE PRECISION :: L_USER_BRDF_F   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: L_USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      DOUBLE PRECISION :: L_EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: L_USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  local variables
!  ===============

!  11/1/19. Some new local variables

      INTEGER ::          I, UI, J, KPHI, IB, Q, O1, O2, OM, OOFF(4), NK, NK2
      DOUBLE PRECISION :: REFL, HELP, HELPC, HELPS, HELP_EMISS, EMISS(4) ! EMISS(16)
      INTEGER ::          COSSIN_MASK(16)

      INTEGER, PARAMETER :: LUM = 1 

      INTEGER ::          K
!      DOUBLE PRECISION :: BSUM

!  define cos/sin mask
!   -- 7/28/21. IMPORTANT. Define mask with 3 entries now.

  !    COSSIN_MASK = (/ 1,1,2,0,1,1,2,0,2,2,1,0,0,0,0,1 /)
      COSSIN_MASK = (/ 1,1,2,0,1,1,2,0,3,3,1,0,0,0,0,1 /)

!  surface factor
!   -- 7/28/21. set to half the value if you are using the full-range azimuth integration
!   -- 7/28/21. Set the emission help-factor ( HELP_EMISS = 2 for the half-range, 1 for the full-range)

      if ( DO_HALF_RANGE ) THEN
         HELP       = DELFAC   ! Alternative
         HELP_EMISS = TWO * DELFAC
      else
         HELP       = HALF * DELFAC
         HELP_EMISS = DELFAC
      endif

!   -- 7/28/21. Introduce +/- signs for cosine/sine

      HELPC = + HELP
      HELPS = - HELP

!  11/1/19. Proxies and offset. New Code.
!   --7/28/21. use NK2 = NBRDF_HALF for full-range, NK2 = NSTREAMS_BRDF for half-range azimuth integration

      NK = NSTREAMS_BRDF ; NK2 = NSTREAMS_BRDF ; if (.not. DO_HALF_RANGE ) NK2 = NBRDF_HALF
      OOFF(1) = 0 ; OOFF(2) = 4 ; OOFF(3) = 8 ; OOFF(4) = 12

!  Zeroing

      L_BRDF_F        = ZERO
      L_BRDF_F_0      = ZERO
      L_USER_BRDF_F   = ZERO
      L_USER_BRDF_F_0 = ZERO

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)
!     !@@ Solar Optionality, added 12/31/12
!     !@@ Rob Fix, 3/20/19. Bug: Lambertian case wrongly used NSTOKESSQ instead of 1
!     !@@ Rob Fix, 11/1/19. Bug: Cosine Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)

!   --7/28/21. IMPORTANT. Cosine Mask now has 3 possibilities.

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO IB = 1, NBEAMS
            DO I = 1, NSTREAMS
!  Start Replacement code --------------------------------
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  OM = OOFF(O1) + O2
                  IF ( COSSIN_MASK(OM) .EQ. 1 ) THEN
                    L_BRDF_F_0(OM,I,IB) = HELPC * DOT_PRODUCT(BRDFUNC_0(OM,I,IB,1:NK),BRDF_COSAZMFAC(1:NK))
                  ELSE IF ( COSSIN_MASK(OM) .EQ. 2 ) THEN
                    L_BRDF_F_0(OM,I,IB) = HELPS * DOT_PRODUCT(BRDFUNC_0(OM,I,IB,1:NK),BRDF_SINAZMFAC(1:NK))
                  ELSE IF ( COSSIN_MASK(OM) .EQ. 3 ) THEN
                    L_BRDF_F_0(OM,I,IB) = - HELPS * DOT_PRODUCT(BRDFUNC_0(OM,I,IB,1:NK),BRDF_SINAZMFAC(1:NK))
                  ENDIF
                ENDDO
              ENDDO
!  End Replacement code ---------------------------------
            ENDDO
          ENDDO
        ELSE IF ( M .EQ. 0 ) THEN
          DO IB = 1, NBEAMS
            DO I = 1, NSTREAMS
              Q = 1 ; L_BRDF_F_0(Q,I,IB) = ONE ! 3/20/19 --> L_BRDF_F_0(1:NSTOKESSQ,I,IB) = ONE
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  incident quadrature directions (surface multiple reflections)
!     !@@ Rob Fix, 3/20/19. Bug: Lambertian case wrongly used NSTOKESSQ instead of 1
!     !@@ Rob Fix, 11/1/19. Bug: Cosine Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)

!   --7/28/21. IMPORTANT. Cosine Mask now has 3 possibilities.

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
         DO I = 1, NSTREAMS
            DO J = 1, NSTREAMS
!  Start Replacement code --------------------------------
               DO O1 = 1, NSTOKES
                  DO O2 = 1, NSTOKES
                     OM = OOFF(O1) + O2
                     IF ( COSSIN_MASK(OM) .EQ. 1 ) THEN
                        L_BRDF_F(OM,I,J) = HELPC * DOT_PRODUCT(BRDFUNC(OM,I,J,1:NK),BRDF_COSAZMFAC(1:NK))
                     ELSE IF ( COSSIN_MASK(OM) .EQ. 2 ) THEN
                        L_BRDF_F(OM,I,J) = HELPS * DOT_PRODUCT(BRDFUNC(OM,I,J,1:NK),BRDF_SINAZMFAC(1:NK))
                     ELSE IF ( COSSIN_MASK(OM) .EQ. 3 ) THEN
                        L_BRDF_F(OM,I,J) = - HELPS * DOT_PRODUCT(BRDFUNC(OM,I,J,1:NK),BRDF_SINAZMFAC(1:NK))
                     ENDIF
                  ENDDO
               ENDDO
!  End Replacement code ---------------------------------
            ENDDO
         ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
         DO I = 1, NSTREAMS
            DO J = 1, NSTREAMS
               Q = 1 ; L_BRDF_F(Q,I,J) = ONE ! 3/20/19 --> L_BRDF_F(1:NSTOKESSQ,I,J) = ONE
            ENDDO
         ENDDO
      ENDIF

!  debug information

!      IF ( DO_DEBUG_WRITE ) THEN
!        WRITE(555,'(A)')'BRDF_1 Fourier 0 quad values'
!        IF ( FOURIER .EQ. 0 ) THEN
!          DO I = 1, NSTREAMS
!          WRITE(555,'(1PE12.5,3x,1P10E12.5)') BIREFLEC_0(1,I,1),(BIREFLEC(1,I,J),J=1,NSTREAMS)
!         ENDDO
!        ENDIF
!      ENDIF

!  albedo check, always calculate the spherical albedo.
!   (Plane albedo calculations are commented out)

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)
!     !@@ Observational Geometry option and Solar Optionality, Installed 12/31/12
!     !@@ Rob Fix, 3/20/19. Bug: Lambertian case wrongly used NSTOKESSQ instead of 1
!     !@@ Rob Fix, 11/1/19. Bug: Cosine Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)

!   --7/28/21. IMPORTANT. Cosine Mask now has 3 possibilities.

        IF ( DO_SOLAR_SOURCES ) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            IF ( .NOT. LAMBERTIAN_FLAG ) THEN
              DO IB = 1, NBEAMS
!  Start Replacement code --------------------------------
                DO O1 = 1, NSTOKES ; DO O2 = 1, NSTOKES
                  OM = OOFF(O1) + O2
                  IF ( COSSIN_MASK(OM) .EQ. 1 ) THEN
                    L_USER_BRDF_F_0(OM,LUM,IB) = HELPC * DOT_PRODUCT(USER_BRDFUNC_0(OM,LUM,IB,1:NK),BRDF_COSAZMFAC(1:NK))
                  ELSE IF ( COSSIN_MASK(OM) .EQ. 2 ) THEN
                    L_USER_BRDF_F_0(OM,LUM,IB) = HELPS * DOT_PRODUCT(USER_BRDFUNC_0(OM,LUM,IB,1:NK),BRDF_SINAZMFAC(1:NK))
                  ELSE IF ( COSSIN_MASK(OM) .EQ. 3 ) THEN
                    L_USER_BRDF_F_0(OM,LUM,IB) = - HELPS * DOT_PRODUCT(USER_BRDFUNC_0(OM,LUM,IB,1:NK),BRDF_SINAZMFAC(1:NK))
                  ENDIF
                  ENDDO ; ENDDO
!  End Replacement code ---------------------------------
              ENDDO
            ELSE IF ( M .EQ. 0 ) THEN
              DO IB = 1, NBEAMS
                Q = 1 ; L_USER_BRDF_F_0(Q,LUM,IB) = ONE ! 3/20/19 --> L_USER_BRDF_F_0(1:NSTOKESSQ,LUM,IB) = ONE 
              ENDDO
            ENDIF
          ELSE
            IF ( .NOT. LAMBERTIAN_FLAG ) THEN
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
!  Start Replacement code --------------------------------
                  DO O1 = 1, NSTOKES ; DO O2 = 1, NSTOKES
                    OM = OOFF(O1) + O2
                    IF ( COSSIN_MASK(OM) .EQ. 1 ) THEN
                      L_USER_BRDF_F_0(OM,UI,IB) = HELPC * DOT_PRODUCT(USER_BRDFUNC_0(OM,UI,IB,1:NK),BRDF_COSAZMFAC(1:NK))
                    ELSE IF ( COSSIN_MASK(OM) .EQ. 2 ) THEN
                      L_USER_BRDF_F_0(OM,UI,IB) = HELPS * DOT_PRODUCT(USER_BRDFUNC_0(OM,UI,IB,1:NK),BRDF_SINAZMFAC(1:NK))
                    ELSE IF ( COSSIN_MASK(OM) .EQ. 3 ) THEN
                      L_USER_BRDF_F_0(OM,UI,IB) = - HELPS * DOT_PRODUCT(USER_BRDFUNC_0(OM,UI,IB,1:NK),BRDF_SINAZMFAC(1:NK))
                    ENDIF
                  ENDDO ; ENDDO
!  End Replacement code ---------------------------------
                ENDDO
              ENDDO
            ELSE IF ( M .EQ. 0 ) THEN
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  Q = 1 ; L_USER_BRDF_F_0(Q,UI,IB) = ONE !  3/20/19 --> L_USER_BRDF_F_0(1:NSTOKESSQ,UI,IB) = ONE
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF

!  incident quadrature directions (surface multiple reflections)
!     !@@ Rob Fix, 3/20/19. Bug: Lambertian case wrongly used NSTOKESSQ instead of 1
!     !@@ Rob Fix, 11/1/19. Bug: Cosine Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)

!   --7/28/21. IMPORTANT. Cosine Mask now has 3 possibilities.

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO UI = 1, N_USER_STREAMS
            DO J = 1, NSTREAMS
!  Start Replacement code --------------------------------
              DO O1 = 1, NSTOKES ; DO O2 = 1, NSTOKES
                OM = OOFF(O1) + O2
                IF ( COSSIN_MASK(OM) .EQ. 1 ) THEN
                  L_USER_BRDF_F(OM,UI,J) = HELPC * DOT_PRODUCT(USER_BRDFUNC(OM,UI,J,1:NK),BRDF_COSAZMFAC(1:NK))
                ELSE IF ( COSSIN_MASK(OM) .EQ. 2 ) THEN
                  L_USER_BRDF_F(OM,UI,J) = HELPS * DOT_PRODUCT(USER_BRDFUNC(OM,UI,J,1:NK),BRDF_SINAZMFAC(1:NK))
                ELSE IF ( COSSIN_MASK(OM) .EQ. 3 ) THEN
                  L_USER_BRDF_F(OM,UI,J) = - HELPS * DOT_PRODUCT(USER_BRDFUNC(OM,UI,J,1:NK),BRDF_SINAZMFAC(1:NK))
                ENDIF
              ENDDO ; ENDDO
!  End Replacement code ---------------------------------
            ENDDO
          ENDDO
        ELSE IF ( M .EQ. 0 ) THEN
          DO UI = 1, N_USER_STREAMS
            DO J = 1, NSTREAMS
              Q = 1 ; L_USER_BRDF_F(Q,UI,J) = ONE ! 3/20/19 --> L_USER_BRDF_F(1:NSTOKESSQ,UI,J) = ONE
            ENDDO
          ENDDO
        ENDIF

!  End user-stream clause

      ENDIF

!  Emissivity
!  ----------

!  Azimuth independent contribution, from Kirchhoff's law

      IF ( DO_SURFEMISS .AND. M .EQ. 0 ) THEN

!  Lambertian case

        IF ( LAMBERTIAN_FLAG ) THEN
          DO I = 1, NSTREAMS
            L_EMISSIVITY(1,I) = FACTOR
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              L_USER_EMISSIVITY(1,UI) = FACTOR
            ENDDO
          ENDIF
        ENDIF

!  bidirectional reflectance

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN

!  Inserted Polarization sum here.   Still to be checked.....!!!!!!!
!     !@@ Rob Fix, 11/1/19. Recoded the BRDF emissivity case

!     !@@ Rob Fix, 11/1/19. Bug: Cosine Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)

!   --7/28/21. IMPORTANT. NK2 = NSTREAMS_BRDF for the HALF_SPACE case.

!  Quadrature polar directions

!mick fix
          EMISS = ZERO

          DO I = 1, NSTREAMS
!  Start Replacement code --------------------------------
            DO O1 = 1, NSTOKES
              DO O2 = 1, NSTOKES
                OM = OOFF(O1) + O2
                REFL = ZERO
                DO KPHI = 1, NSTREAMS_BRDF
                  REFL = REFL + A_BRDF(KPHI) * DOT_PRODUCT(EBRDFUNC(OM,I,1:NK2,KPHI),BAX_BRDF(1:NK2))
                ENDDO
                EMISS(O2) = REFL
              ENDDO
              L_EMISSIVITY(O1,I) = HELP_EMISS * FACTOR * SUM(EMISS(1:NSTOKES))
            ENDDO
!  End Replacement code ---------------------------------

          ENDDO

!  User-defined polar directions

          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
!  Start Replacement code --------------------------------
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  OM = OOFF(O1) + O2
                  REFL = ZERO
                  DO KPHI = 1, NSTREAMS_BRDF
                    REFL = REFL + A_BRDF(KPHI) * DOT_PRODUCT(USER_EBRDFUNC(OM,UI,1:NK2,KPHI),BAX_BRDF(1:NK2))
                  ENDDO
                  EMISS(O2) = REFL
                ENDDO
                L_USER_EMISSIVITY(O1,UI) = HELP_EMISS * FACTOR * SUM(EMISS(1:NSTOKES))
              ENDDO
!  End Replacement code ---------------------------------
            ENDDO
          ENDIF

        ENDIF

!  End emissivity clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_FOURIER

!  End module

      END MODULE vbrdf_sup_routines_m

