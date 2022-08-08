
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
! #            VFO_LPS_MASTER_INTERFACE                         #
! #                                                             #
! ###############################################################

!  Version 2.0 - 2.7. Internal SSCORR/DBCORR routines (Old)
!  Version 2.7. VFO interface to FO code version 1.4 was written
!  Version 2.8. VFO interface to FO code version 1.5 (upgrade)
!  Version 2.8. Internal SSCORR/DBCORR routines removed.

!   -- 4/9/19. Need to add the FO Surface-leaving assignation and saved cumulative transmittance
!              Also their linearizations

!  1/31/21. Version 2.8.3. DO_MSSTS option final installation
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_UP, THETA_ALL, ALPHA (upwelling)
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_DN, THETA_ALL, ALPHA (downwelling)

!  1/31/21. Version 2.8.3. Add DO_DOUBLET Geometry flag, plus offset input arguments

      MODULE vlidort_vfo_lps_interface_m

      USE VFO_LinMasters_m, Only : VFO_LPS_MASTER

      PUBLIC ::  VFO_LPS_MASTER_INTERFACE

      CONTAINS

      SUBROUTINE VFO_LPS_MASTER_INTERFACE ( &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                           & ! Input Sources flags
        DO_PLANE_PARALLEL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, DO_FMATRIX,                   & ! Input SS control flags
        DO_DELTAM, DO_UPWELLING, DO_DNWELLING, DO_PARTLAYERS, DO_OBSGEOM, DO_DOUBLET,         & ! Input RT Control flags
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,         & ! Input Optical/Surface flags
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS,                    & ! Input Jacobian flags
        NSTOKES, NLAYERS, NFINELAYERS, NGREEK_MOMENTS_INPUT, ND_OFFSET, NV_OFFSET, NA_OFFSET, & ! Input numbers/offsets
        N_SZANGLES, SZANGLES, N_USER_VZANGLES, USER_VZANGLES, N_USER_RELAZMS, USER_RELAZMS,   & ! Input geometry
        N_USER_LEVELS, LEVELMASK_UP, LEVELMASK_DN, N_PARTLAYERS,                              & ! Input level-control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,      & ! Input partial-control
        LVARY_FLAG, LVARY_NUMBER, N_SLEAVE_WFS, N_REFLEC_WFS, N_SURFACE_WFS,                & ! Input numbers (Jacobians)
        EARTH_RADIUS, HEIGHT_GRID, SS_FLUX_MULTIPLIER, FLUXVEC,                             & ! Input Flux/Heights
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                         & ! Inputs (Optical - Regular)
        DELTAU_VERT, FMATRIX_UP, FMATRIX_DN, TRUNC_FACTOR, BB_INPUT,                        & ! Inputs (Optical - Regular)
        LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,                        & ! Inputs (Optical - Surface)
        SLTERM_ISOTROPIC, SLTERM_USERANGLES,                                                & ! Inputs (Optical - Surface)
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,                   & ! Inputs (Optical - Lin Atmos)
        L_DELTAU_VERT, L_TRUNC_FACTOR, L_FMATRIX_UP, L_FMATRIX_DN,                          & ! Inputs (Optical - Lin Atmos)
        LS_EXACTDB_BRDFUNC, LS_USER_EMISSIVITY,                                             & ! Inputs (Optical - Lin Surf)
        LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES,                                      & ! Inputs (Optical - Lin Surf)
        FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,                           & ! Output - Stokes vectors
        FO_PROFILEWF_SS,  FO_PROFILEWF_DB, FO_PROFILEWF_DTA, FO_PROFILEWF_DTS,              & ! Output - Profile Jacobians
        FO_SURFACEWF_DB,  FO_SURFACEWF_DTS,                                                 & ! Output - Surface Jacobians
        FO_STOKES_ATMOS,  FO_STOKES_SURF, FO_STOKES,                                        & ! Output - Stokes composites
        FO_PROFILEWF_ATMOS, FO_PROFILEWF_SURF, FO_PROFILEWF, FO_SURFACEWF,                  & ! Output - Jacobian composites
        CUMTRANS, LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA, SLTERM,                       & ! Output - Auxiliary
        LP_CUMTRANS, LSSL_SLTERM, LP_LOSTRANS_UP, LP_LOSTRANS_DN,                           & ! Output - Auxiliary
        FAIL, MESSAGE, TRACE_1, TRACE_2 )                                                     ! Output (Exception handling)

!  1/31/21. Version 2.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA to output list
!    ==> Last 3 lines reorganized output list.

!  1/31/21. Version 2.8.3. Add Doublet geometry flag and include offsets

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_LEVELS, MAXLAYERS, &
                                 MAXSTOKES, MAXMOMENTS_INPUT, MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,      &
                                 MAX_USER_STREAMS, MAXBEAMS, MAX_GEOMETRIES, MAXFINELAYERS, MAXMOMENTS,         &
                                 MAXSTOKES_SQ, MAX_PARTLAYERS, MAX_DIRECTIONS, ZERO, ONE, PIE, DEG_TO_RAD, SMALLNUM

      IMPLICIT NONE

!  Parameter argument
!  ------------------

      INTEGER, PARAMETER  :: FFP = SELECTED_REAL_KIND(15)

!  Inputs
!  ======

!  Flags. Fmatrix flag (7/7/16), Surface leaving inputs (8/3/16), Introduced for Version 2.8
!     Water-leaving flag introduced for Version 2.8.1. 4/9/19.


      LOGICAL, INTENT(IN) ::            DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) ::            DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) ::            DO_SURFACE_EMISSION

      LOGICAL, INTENT(IN) ::            DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) ::            DO_SSCORR_NADIR
      LOGICAL, INTENT(IN) ::            DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::            DO_FMATRIX

      LOGICAL, INTENT(IN) ::            DO_UPWELLING
      LOGICAL, INTENT(IN) ::            DO_DNWELLING
      LOGICAL, INTENT(IN) ::            DO_DELTAM

!  1/31/21. Version 2.8.3. Add Doublet Geometry option

      LOGICAL, INTENT(IN) ::            DO_OBSGEOM
      LOGICAL, INTENT(IN) ::            DO_DOUBLET

      LOGICAL, INTENT(IN) ::            DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::            DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) ::            DO_WATER_LEAVING    ! 4/9/19 Added
      LOGICAL, INTENT(IN) ::            DO_SL_ISOTROPIC

!  Linearization flags. Sleave flag 8/3/16.

      LOGICAL, INTENT(IN) ::            DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT(IN) ::            DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT(IN) ::            DO_SLEAVE_WFS

!  Control integers

      INTEGER, INTENT(IN) ::            NSTOKES
      INTEGER, INTENT(IN) ::            NLAYERS
      INTEGER, INTENT(IN) ::            NFINELAYERS
      INTEGER, INTENT(IN) ::            NGREEK_MOMENTS_INPUT

!  1/31/21. Version 2.8.3. Add OFFSETS inputs (Doublet ND and lattice NV/NA)

      INTEGER, INTENT(IN) ::            ND_OFFSET ( MAX_SZANGLES )
      INTEGER, INTENT(IN) ::            NV_OFFSET ( MAX_SZANGLES )
      INTEGER, INTENT(IN) ::            NA_OFFSET ( MAX_SZANGLES, MAX_USER_VZANGLES )

!  Linearization control
!   Sleave control 8/3/16. Note that N_SURFACE_WFS = N_REFLEC_WFS + N_SLEAVE_WFS

      LOGICAL, INTENT(IN) ::            LVARY_FLAG ( MAXLAYERS )
      INTEGER, INTENT(IN) ::            LVARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT(IN) ::            N_SURFACE_WFS
      INTEGER, INTENT(IN) ::            N_REFLEC_WFS
      INTEGER, INTENT(IN) ::            N_SLEAVE_WFS

!  Geometry and Levels

      INTEGER, INTENT(IN) ::            N_SZANGLES
      DOUBLE PRECISION, INTENT(IN) ::   SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT(IN) ::            N_USER_VZANGLES
      DOUBLE PRECISION, INTENT(IN) ::   USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER, INTENT(IN) ::            N_USER_RELAZMS
      DOUBLE PRECISION, INTENT(IN) ::   USER_RELAZMS  ( MAX_USER_RELAZMS )

!  Levels. Now require the Level Mask inputs, 9/17/16.

      INTEGER, INTENT(IN) ::            N_USER_LEVELS
      INTEGER, INTENT (IN) ::           LEVELMASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           LEVELMASK_DN  ( MAX_USER_LEVELS )

!      DOUBLE PRECISION, INTENT(IN) ::   USER_LEVELS ( MAX_USER_LEVELS )

!  PARTIAL-Layer inputs, added 9/17/16

      LOGICAL, INTENT (IN) ::           DO_PARTLAYERS
      INTEGER, INTENT (IN) ::           N_PARTLAYERS
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

!  other inputs

      DOUBLE PRECISION, INTENT(IN) ::   EARTH_RADIUS
      DOUBLE PRECISION, INTENT(IN) ::   HEIGHT_GRID ( 0:MAXLAYERS )

      DOUBLE PRECISION, INTENT(IN) ::   SS_FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT(IN) ::   FLUXVEC ( MAXSTOKES )

!  Optical unscaled

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) ::   BB_INPUT ( 0:MAXLAYERS )

!  optical scaled

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   TRUNC_FACTOR ( MAXLAYERS )

!  F-matrices, introduced for Version 2.8 of VLIDORT, and Version 1.5 of FO
!mick fix 9/19/2017 - swapped layer & geo dimensions

      !DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_UP ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      !DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_DN ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION, INTENT(IN) ::   FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  surface BRDF and EMiss.

      DOUBLE PRECISION, INTENT(IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT(IN) ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT(IN) ::   SURFBB
      DOUBLE PRECISION, INTENT(IN) ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  Surface leaving inputs, 8/3/16, for Version 2.8

      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_ISOTROPIC  ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized optical properties

      DOUBLE PRECISION, INTENT(IN) ::   L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT(IN) ::   L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )

!mick fix 9/19/2017 - swapped layer & geo indices

      !DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_UP ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 6 )
      !DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_DN ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_UP ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES,  6 )
      DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_DN ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES,  6 )

!  Linearized surface properties

      DOUBLE PRECISION, INTENT(IN) ::   LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )

!  Surface leaving linearizations. New for Version 2.8, 8/3/16

      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_ISOTROPIC ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_USERANGLES &
          ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!  Outputs
!  =======
!mick mod 9/19/2017 - added FO_STOKES_ATMOS, FO_STOKES_SURF, FO_PROFILEWF_ATMOS, FO_PROFILEWF_SURF (new output from vector FO code)

!  Standard
!  --------

!  Solar

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Thermal

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DTA ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DTS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Composite

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_ATMOS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_SURF  ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES       ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  4/9/19. Additional output for the sleave correction

      DOUBLE PRECISION, INTENT (INOUT)  :: CUMTRANS ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  1/31/21. Version 2.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA to output list

      DOUBLE PRECISION, Intent(INOUT)   :: LOSTRANS_UP ( max_geometries, maxlayers )
      DOUBLE PRECISION, Intent(INOUT)   :: LOSTRANS_DN ( max_geometries, maxlayers )

      DOUBLE PRECISION, Intent(INOUT)   :: theta_all ( 0:maxlayers, max_geometries )
      DOUBLE PRECISION, Intent(INOUT)   :: alpha     ( 0:maxlayers, max_geometries )

!  4/9/19. Surface leaving FO assignation

      DOUBLE PRECISION, INTENT (INOUT)  :: SLTERM (  MAXSTOKES, MAX_GEOMETRIES )

!  Linearized
!  ----------

!  Solar

      DOUBLE PRECISION, INTENT (INOUT) :: FO_PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Thermal

      DOUBLE PRECISION, INTENT (INOUT) :: FO_PROFILEWF_DTA &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_PROFILEWF_DTS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF_DTS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Composite

      DOUBLE PRECISION, INTENT (INOUT) :: FO_PROFILEWF_ATMOS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_PROFILEWF_SURF &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_PROFILEWF &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  4/9/19. Additional linearized output for the sleave correction. 5/31/21. Dimensioning re-ordered

!      DOUBLE PRECISION, INTENT (INOUT) :: LP_CUMTRANS ( MAX_ATMOSWFS, maxlayers, max_user_levels, max_geometries )
      DOUBLE PRECISION, INTENT (INOUT) :: LP_CUMTRANS ( max_user_levels, max_geometries, maxlayers, MAX_ATMOSWFS )

!  4/9/19. Surface leaving FO assignation
!mick fix 5/19/2021 - interchanged max dimensions

      DOUBLE PRECISION, INTENT (INOUT) :: LSSL_SLTERM ( MAXSTOKES, MAX_GEOMETRIES, MAX_SLEAVEWFS )

!  1/31/21. Version 2.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LP_LOSTRANS_UP, LP_LOSTRANS_DN to output list

      DOUBLE PRECISION, Intent(INOUT)  :: LP_LOSTRANS_UP ( max_geometries, maxlayers, Max_atmoswfs )
      DOUBLE PRECISION, Intent(INOUT)  :: LP_LOSTRANS_DN ( max_geometries, maxlayers, Max_atmoswfs )

!  Exception Handling

      LOGICAL, INTENT (OUT)             :: FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE_1, TRACE_2

!  FO VARIABLES TAKEN DIRECTLY FROM THE VLIDORT INPUT
!  ==================================================

!  Dimensioning integers

       !INTEGER   :  MAXLAYERS, MAXMOMENTS_INPUT, MAX_USER_LEVELS, max_atmoswfs, max_surfacewfs
       !INTEGER   :: N_USER_LEVELS, NSTOKES, N_USER_LEVELS

!  Sources control, including thermal, deltam_scaling, directional flags

      !LOGICAL   :: DO_SOLAR_SOURCES
      !LOGICAL   :: DO_THERMAL_EMISSION
      !LOGICAL   :: DO_SURFACE_EMISSION
      !LOGICAL   :: DO_DELTAM
      !LOGICAL   :: DO_UPWELLING, DO_DNWELLING

!  FluxVector

      !REAL(FFP) :: FLUXVEC ( MAXSTOKES )

!  surface property

      !REAL(FFP) :: SURFBB

!  Local variables
!  ===============

!  Max dimensions
!  --------------

      INTEGER   :: MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, MAXFINE

!  Dimensions
!  ----------

!  Layer and geometry control. Finelayer divisions may be changed

      INTEGER   :: NGEOMS, NSZAS, NVZAS, NAZMS, NFINE, NMOMENTS_INPUT

!  Number of surface weighting functions

      INTEGER   :: N_REFLECWFS
      INTEGER   :: N_SLEAVEWFS
      INTEGER   :: N_SURFACEWFS

!  Configuration inputs
!  --------------------

!  Flags (sphericity flags should be mutually exclusive)

      LOGICAL   :: DO_PLANPAR
      LOGICAL   :: DO_ENHANCED_PS

!  Lambertian surface flag

      LOGICAL   :: DO_LAMBERTIAN

!  Vector sunlight flag

      LOGICAL   :: DO_SUNLIGHT

!  Linearization flags

      LOGICAL   :: DO_PROFILEWFS, DO_SURFACEWFS, DO_SLEAVEWFS

!  General inputs
!  --------------

!  DTR = degrees-to-Radians

      REAL(FFP) :: DTR

!  Critical adjustment for cloud layers

      LOGICAL   :: DoCrit
      REAL(FFP) :: Acrit

!  Earth radius + heights. Partial heights added 9/17/16

      REAL(FFP) :: ERADIUS
      REAL(FFP) :: HEIGHTS ( 0:MAXLAYERS )
      REAL(FFP) :: PARTIAL_HEIGHTS ( MAX_PARTLAYERS )

!  Output levels. Not required now, 9/17/16.
!      INTEGER   :: FO_USER_LEVELS ( MAX_USER_LEVELS )

!  Geometry inputs
!  ---------------

!  Input angles (Degrees). Note the dimensioning now,

      REAL(FFP) :: obsgeom_boa  ( MAX_GEOMETRIES, 3 )
      REAL(FFP) :: theta_boa ( MAX_SZANGLES )        !SZA
      REAL(FFP) :: alpha_boa ( MAX_USER_VZANGLES )   !UZA
      REAL(FFP) :: phi_boa   ( MAX_USER_RELAZMS )    !RAA

!  Optical inputs
!  --------------

!  Solar flux

      REAL(FFP) :: FLUX

!  Atmosphere

      REAL(FFP) :: EXTINCTION  ( MAXLAYERS )
      REAL(FFP) :: DELTAUS     ( MAXLAYERS )
      REAL(FFP) :: OMEGA       ( MAXLAYERS )
      REAL(FFP) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, MAXSTOKES, MAXSTOKES )
!mick fix 9/19/2017 - added this set of F-matrix arrays for proper passing of linearized
!                     FO F-matrix input
      REAL(FFP) :: FO_L_FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6, MAX_ATMOSWFS )
      REAL(FFP) :: FO_L_FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6, MAX_ATMOSWFS )

!  For TMS correction

      REAL(FFP) :: TRUNCFAC    ( MAXLAYERS )

!  Thermal inputs. Use directly, don't need to copy
!      REAL(FFP) :: BB_INPUT ( 0:MAXLAYERS )

!  Surface properties - reflective (could be the albedo)

      REAL(FFP) :: REFLEC ( MAXSTOKES, MAXSTOKES, MAX_GEOMETRIES )

!  Surface properties - emissivity
!    - polarized surface emissivity. 12/11/17 Rob Add.

!      REAL(FFP) :: EMISS ( MAX_USER_VZANGLES )
      REAL(FFP) :: EMISS ( MAXSTOKES, MAX_USER_VZANGLES )

!  Linearized inputs
!  -----------------

!  Linearization control

      LOGICAL   :: LVARYFLAGS ( MAXLAYERS )
      INTEGER   :: LVARYNUMS  ( MAXLAYERS )
      LOGICAL   :: LVARYMOMS  ( MAXLAYERS, MAX_ATMOSWFS )

!  Linearized optical inputs

      REAL(FFP) :: L_EXTINCTION ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_DELTAUS    ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_OMEGA      ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_GREEKMAT   ( MAXLAYERS, 0:MAXMOMENTS_INPUT, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Linearized TMS correction

      REAL(FFP) :: L_TRUNCFAC ( MAXLAYERS, MAX_ATMOSWFS )

!  Surface properties - reflective

      REAL(FFP) :: LS_REFLEC ( MAXSTOKES, MAXSTOKES, MAX_GEOMETRIES, MAX_SURFACEWFS )

!  Surface properties - emissivity
!    - polarized surface emissivity. 12/11/17 Rob Add.

!      REAL(FFP) :: LS_EMISS  ( MAX_USER_VZANGLES, MAX_SURFACEWFS )
      REAL(FFP) :: LS_EMISS  ( MAXSTOKES, MAX_USER_VZANGLES, MAX_SURFACEWFS )

!  Local variables
!  ---------------

!    -- 1/31/21. Version 2.8.3. do_Polarized_Emissivity added, Offsets now input

      LOGICAL   :: do_Polarized_Emissivity
      INTEGER   :: ns, nv, na, G, GK, I, L, N, O1, O2, PAR, SPAR, UT
      INTEGER   :: LUM=1, LUA=1

!  Define VFO_LPS_MASTER inputs
!  ============================

!  Note: argument passing for variables defined elsewhere commented out here

!  Max dimensions

      maxgeoms          = MAX_GEOMETRIES
      maxfine           = MAXFINELAYERS
      maxszas           = MAX_SZANGLES
      maxvzas           = MAX_USER_VZANGLES
      maxazms           = MAX_USER_RELAZMS

!  Proxies
!     -- 1/31/21. Version 2.8.3. Add Double toption

      if ( do_ObsGeom ) then
         ngeoms  = N_SZANGLES
      else if ( do_Doublet ) then
         ngeoms  = N_SZANGLES * N_USER_VZANGLES
      else
         ngeoms  = N_SZANGLES * N_USER_VZANGLES * N_USER_RELAZMS
      endif
      nszas   = N_SZANGLES
      nvzas   = N_USER_VZANGLES
      nazms   = N_USER_RELAZMS
      nfine           = NFINELAYERS
      nmoments_input  = NGREEK_MOMENTS_INPUT

      n_reflecwfs     = N_REFLEC_WFS
      n_sleavewfs     = N_SLEAVE_WFS
      n_surfacewfs    = N_SURFACE_WFS

!  Offsets. -- 1/31/21. Version 2.8.3. Removed from here. Now inputs
!   na_offset = 0 ; nv_offset = 0
!   if ( .not. do_obsgeom ) then
!     if ( do_doublet ) then
!       do ns = 1, nszas
!         nd_offset(ns) = nvzas * (ns - 1) 
!       enddo
!     else
!       do ns = 1, nszas ;  do nv = 1, nvzas
!           na_offset(ns,nv) = nvzas * nazms * (ns - 1) + nazms * (nv - 1)
!       enddo ; enddo
!     endif
!   endif

!  Configuration inputs. Removed "regular_ps" flag 9/17/16 (Redundant)

      if (DO_PLANE_PARALLEL) then
        do_planpar     = .TRUE.
        do_enhanced_ps = .FALSE.
      else
        do_planpar = .FALSE.
        if (DO_SSCORR_NADIR) then
          do_enhanced_ps = .FALSE.
        else if (DO_SSCORR_OUTGOING) then
          do_enhanced_ps = .TRUE.
        end if
      end if

      do_lambertian      = DO_LAMBERTIAN_SURFACE
      do_sunlight        = .TRUE. !set for now

!  set local FO flags for linearization

      do_profilewfs      = DO_PROFILE_LINEARIZATION
      do_surfacewfs      = DO_SURFACE_LINEARIZATION
      do_sleavewfs       = DO_SLEAVE_WFS

!  General inputs

      dtr    = DEG_TO_RAD

!      DoCrit = .FALSE. !set for now
!      Acrit  = ZERO    !set for now

!  Rob Fix 3/16/15. Default values changed 
!   [ Note for future - these should be inputs ]

      DoCrit = .TRUE.     
      Acrit  = 1.0d-10 

!  Earth Radius and heights. Partials added, 9/17/16

      eradius             = EARTH_RADIUS
      heights(0:nlayers)  = HEIGHT_GRID(0:NLAYERS)
      do ut = 1, n_partlayers
        n = partlayers_layeridx(ut)
        partial_heights(ut) = heights(n-1) - partlayers_values(ut) * ( heights(n-1) - heights(n) )
      enddo

!  This code no longer required, now using same arrays as VLIDORT
!      fo_user_levels(1:n_user_levels) = INT(user_levels(1:n_user_levels))

!  Geometry inputs. all angles in Degrees.
!  These are good for Lattice/Doublet and ObsGeom (1/31/21. Version 2.8.3)

      obsgeom_boa = zero
      if ( do_ObsGeom) then
         obsgeom_boa(1:ngeoms,1) = SZANGLES      (1:ngeoms)
         obsgeom_boa(1:ngeoms,2) = USER_VZANGLES (1:ngeoms)
         obsgeom_boa(1:ngeoms,3) = USER_RELAZMS  (1:ngeoms)
      endif
      theta_boa(1:N_SZANGLES)       = SZANGLES(1:N_SZANGLES)
      alpha_boa(1:N_USER_VZANGLES)  = USER_VZANGLES(1:N_USER_VZANGLES)
      phi_boa  (1:N_USER_RELAZMS)   = USER_RELAZMS (1:N_USER_RELAZMS)

!  Optical inputs

      flux     = SS_FLUX_MULTIPLIER !=FLUX_FACTOR/(4.0d0*PIE)

! @@ Rob 7/30/13
!   If deltam-scaling, must use scaled deltau !!!

      extinction = 0.0d0
      if (do_deltam) then
        extinction(1:nlayers) = DELTAU_VERT(1:NLAYERS) / (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT(1:NLAYERS)
      else
        extinction(1:nlayers) = DELTAU_VERT_INPUT(1:NLAYERS) / (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT_INPUT(1:NLAYERS)
      endif

!  Single scattering albedo is unscaled

      omega(1:nlayers)      = OMEGA_TOTAL_INPUT(1:NLAYERS)

!  Version 2.8 (FO 1.5). Optional use of Fmatrix input
!    Greekmat will not be needed in this case

      greekmat = zero 
      if ( .not. do_fmatrix ) then
        do o1=1,nstokes
          do o2=1,nstokes
            gk = 4*(o1-1) + o2
            do l=0,nmoments_input
              do n=1,nlayers
                greekmat(n,l,o1,o2) = GREEKMAT_TOTAL_INPUT(L,N,GK)
              end do
            end do
          end do
        end do
      endif

! For TMS correction only

      truncfac = zero
      if (DO_DELTAM) then
        do n=1,nlayers
          truncfac(n) = TRUNC_FACTOR(N)
        end do
      end if

!  Version 2.8. Copying not needed
!      bb_input(0:nlayers) = THERMAL_BB_INPUT(0:NLAYERS)

!  BRDF reflection now also for Lattice.
!     -- 1/31/21. Version 2.8.3. (Upgrade 4/15/20, Version 2.8.2). Add Do_Doublet option for BRDF reflection.

      reflec = ZERO
      if (do_lambertian_surface) then
        reflec(1,1,1:ngeoms) = LAMBERTIAN_ALBEDO
      else
        if ( do_ObsGeom ) then
          do o1=1,nstokes ; do o2=1,nstokes
            gk = 4*(o1-1) + o2
            reflec(o1,o2,1:ngeoms) = EXACTDB_BRDFUNC(gk,lum,lua,1:ngeoms)
          end do ; end do
        else  if ( do_doublet ) then
          do ns = 1, nszas ; do nv = 1, nvzas
            g = nd_offset(ns) + nv
            do o1=1,nstokes ; do o2=1,nstokes
              gk = 4*(o1-1) + o2
              reflec(o1,o2,g) = EXACTDB_BRDFUNC(gk,nv,LUA,ns)
            end do ; end do
          end do ; end do
        else
          do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
            g = na_offset(ns,nv) + na
            do o1=1,nstokes ; do o2=1,nstokes
              gk = 4*(o1-1) + o2
              reflec(o1,o2,g) = EXACTDB_BRDFUNC(gk,nv,na,ns)
            end do ; end do
          enddo ; enddo ; enddo
        end if
      end if

!  Surface leaving, introduced 8/3/16. Proper Initializing, 8/3/16
! 4/9/19 ==> This is now an output.
!     -- 1/31/21. Version 2.8.3. (Upgrade 4/15/20, Version 2.8.2). Add Do_Doublet option for surface leaving.

      SLTERM = zero
      if ( do_surface_leaving ) then
        if ( do_ObsGeom ) then
          IF ( DO_SL_ISOTROPIC ) THEN
            SLTERM(1,1:ngeoms) = SLTERM_ISOTROPIC(1,1:ngeoms)
          ELSE
            DO O1 = 1, NSTOKES
              SLTERM(O1,1:ngeoms) = SLTERM_USERANGLES(O1,LUM,LUA,1:ngeoms)
            ENDDO
          ENDIF
        else if ( do_Doublet ) then
          IF ( DO_SL_ISOTROPIC ) THEN
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               SLTERM(1,g) = SLTERM_ISOTROPIC(1,ns)
            end do ; end do 
          ELSE
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               do o1 = 1,nstokes
                  SLTERM(O1,g) = SLTERM_USERANGLES(O1,nv,LUA,ns)
               end do
            end do ; end do
          ENDIF
        else
          IF ( DO_SL_ISOTROPIC ) THEN
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               SLTERM(1,g) = SLTERM_ISOTROPIC(1,ns)
            end do ; end do ; end do 
         ELSE
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               do o1 = 1,nstokes
                  SLTERM(O1,g) = SLTERM_USERANGLES(O1,nv,na,ns)
               end do
            end do ; end do ; end do 
          ENDIF
        endif
      endif

!  Emissivity. Proper Initializing, 8/3/16
!    - polarized surface emissivity. 12/11/17 Rob Add.
!    -- 1/31/21. Version 2.8.3. Flag set here.

      emiss = zero ; do_Polarized_Emissivity = .false.
      if ( do_surface_emission ) then
        if ( do_ObsGeom ) then
          do g=1,ngeoms
            emiss(1:nstokes,g) = USER_EMISSIVITY(1:nstokes,g)
            if ( nstokes.gt.1 .and. (SUM(USER_EMISSIVITY(2:nstokes,g)).ne.zero) ) do_Polarized_Emissivity = .true.
          end do
        else
          do nv=1,nvzas
            emiss(1:nstokes,nv) = USER_EMISSIVITY(1:nstokes,nv)
            if ( nstokes.gt.1 .and. (SUM(USER_EMISSIVITY(2:nstokes,nv)).ne.zero) ) do_Polarized_Emissivity = .true.
          end do
        end if
      endif

!  Linearized control
!  ==================

      Lvaryflags(1:nlayers) = LVARY_FLAG(1:NLAYERS)
      Lvarynums(1:nlayers)  = LVARY_NUMBER(1:NLAYERS)
      !Default setting: see below for modification when
      !Greek moments vary
      Lvarymoms(1:nlayers,1:max_atmoswfs) = .FALSE.

!  Linearized optical

      !Atmospheric quantities
      L_extinction = ZERO
      L_deltaus    = ZERO
      L_omega      = ZERO
      L_truncfac   = ZERO
      L_greekmat   = ZERO

      !Recall that VLIDORT takes FULLY NORMALIZED linearized
      !atmospheric inputs
      !                      (x/y)*(dy/dx)
      !whereas the FO code only takes PARTIALLY NORMALIZED linearized
      !atmospheric inputs
      !                        x*(dy/dx)
      !thus, we must compensate for this difference when formulating the
      !inputs here!

      do n=1,nlayers
        if (Lvaryflags(n)) then
          do par=1,Lvarynums(n)

            if (do_deltam) then
              L_extinction(n,par) = (DELTAU_VERT(N) / (HEIGHT_GRID(N-1)-HEIGHT_GRID(N))) * L_DELTAU_VERT(PAR,N)
              L_deltaus(n,par)    = DELTAU_VERT(N) * L_DELTAU_VERT(PAR,N)
            else
              L_extinction(n,par) = (DELTAU_VERT_INPUT(N) / (HEIGHT_GRID(N-1)-HEIGHT_GRID(N))) * L_DELTAU_VERT_INPUT(PAR,N)
              L_deltaus(n,par)    = DELTAU_VERT_INPUT(N)*L_DELTAU_VERT_INPUT(PAR,N)
            endif

            L_omega(n,par)      = OMEGA_TOTAL_INPUT(N)*L_OMEGA_TOTAL_INPUT(PAR,N)

            if (do_deltam) L_truncfac(n,par) = L_TRUNC_FACTOR(PAR,N)

!  Version 2.8 (FO 1.5). Optional use of Fmatrix input
!    Linearized Greekmat will not be needed in this case

!  Check for variation of Greek moments associated with
!    Jacobian wrt current atmospheric parameter
!mick fix 9/19/2017 - added .not. to the block IF condition (i.e. using traditional GMat
!                     moments, NOT the new FMat facility) 
!                   - changed linearized GMat moment IF condition
!                   - added ELSE section for proper passing of linearized FO F-matrix input

            if ( .not. do_Fmatrix ) then
              do o1=1,nstokes
                do o2=1,nstokes
                  gk = 4*(o1-1) + o2
                  do l=0,nmoments_input
                    L_greekmat(n,l,o1,o2,par) = GREEKMAT_TOTAL_INPUT(L,N,GK) * L_GREEKMAT_TOTAL_INPUT(PAR,L,N,GK)
                    !if ( (l==1).and.(gk==1).and.(L_GREEKMAT_TOTAL_INPUT(PAR,L,N,GK) > 1.0d-8) ) Lvarymoms(n,par) = .TRUE.
                    if ( (gk==1).and.(ABS(L_GREEKMAT_TOTAL_INPUT(PAR,L,N,GK)) >= 1000.0d0*SMALLNUM) ) Lvarymoms(n,par) = .TRUE.
                  end do
                end do
              end do
            else
              do g=1,ngeoms
                do i=1,6
                  FO_L_FMATRIX_UP(n,g,i,par) = L_FMATRIX_UP(par,n,g,i)
                  FO_L_FMATRIX_DN(n,g,i,par) = L_FMATRIX_DN(par,n,g,i)
                end do
              end do
            endif

          end do
        end if
      end do

!  Surface quantities
!  ------------------

!  Initialize

      LS_reflec   = ZERO
      LS_emiss    = ZERO
      LSSL_slterm = ZERO

!  reflectance
!     -- 1/31/21. Version 2.8.3. (Upgrade 4/15/20, Version 2.8.2). Add Do_Doublet option 

      if (do_surfacewfs) then
        if (do_lambertian) then
          do spar=1,n_reflecwfs
            do g=1,ngeoms
              LS_reflec(1,1,1:ngeoms,1:n_reflecwfs) = ONE
            end do
          end do
        else
          do spar=1,n_reflecwfs
            if ( do_ObsGeom ) then
              do g=1,ngeoms
                do o1=1,nstokes ; do o2=1,nstokes
                  gk = 4*(o1-1) + o2
                  LS_reflec(o1,o2,g,spar) = LS_EXACTDB_BRDFUNC(SPAR,GK,LUM,LUA,G)
                end do ; end do
              end do
            else if ( do_Doublet ) then
              do ns = 1, nszas ; do nv = 1, nvzas
                g = nd_offset(ns) + nv
                do o1=1,nstokes ; do o2=1,nstokes
                  gk = 4*(o1-1) + o2
                        LS_reflec(o1,o2,g,spar) = LS_EXACTDB_BRDFUNC(SPAR,GK,NV,LUA,NS)
                end do ; end do
              end do ; end do
            else
              do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                g = na_offset(ns,nv) + na
                do o1=1,nstokes ; do o2=1,nstokes
                  gk = 4*(o1-1) + o2
                  LS_reflec(o1,o2,g,spar) = LS_EXACTDB_BRDFUNC(SPAR,GK,NV,NA,NS)
                end do ; end do
              end do ; end do ; end do
            endif
          end do
        end if
      endif

!  Emissivity
!mick fix 9/19/2017 - added "do_surface_emission" to if condition
!    - polarized surface emissivity. 12/11/17 Rob Add.

      if ( do_surface_emission .and. do_surfacewfs ) then
        do spar=1,n_reflecwfs
          if ( do_ObsGeom ) then
            do g=1,ngeoms
!             LS_emiss(g,spar) = LS_USER_EMISSIVITY(SPAR,1,G)
             LS_emiss(1:nstokes,g,spar) = LS_USER_EMISSIVITY(spar,1:nstokes,g)
            end do
          else
            do nv=1,nvzas
!             LS_emiss(nv,spar) = LS_USER_EMISSIVITY(SPAR,Nstokes,,nv)
             LS_emiss(1:nstokes,nv,spar) = LS_USER_EMISSIVITY(spar,1:nstokes,nv)
            end do
          endif
        end do
      end if

!  Surface leaving, New section 8/3/16 for Version 2.8.  4/9/19 ==> This is now an output.
!     -- 1/31/21. Version 2.8.3. (Upgrade 4/15/20, Version 2.8.2). Add Do_Doublet option.
!mick fix 1/5/2021 - added "nv" to "nd_offset" for "g" under doublet (2 places)

      if ( do_surface_leaving .and. do_sleave_wfs ) then
         do spar = 1, n_sleavewfs
            if ( do_ObsGeom ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  LSSL_SLTERM(1,1:ngeoms,spar) = LSSL_SLTERM_ISOTROPIC(spar,1,1:ngeoms)
               ELSE
                  DO O1 = 1, NSTOKES
                     LSSL_SLTERM(O1,1:ngeoms,spar) = LSSL_SLTERM_USERANGLES(spar,O1,LUM,LUA,1:ngeoms)
                  ENDDO
               ENDIF
            else if ( do_Doublet ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  do ns = 1, nszas ; do nv = 1, nvzas 
                     g = nd_offset(ns) + nv
                     LSSL_SLTERM(1,g,spar) = LSSL_SLTERM_ISOTROPIC(spar,1,ns)
                  end do ; end do
               ELSE
                  do ns = 1, nszas ; do nv = 1, nvzas
                     g = nd_offset(ns) + nv
                     do o1 = 1,nstokes
                        LSSL_SLTERM(O1,g,spar) = LSSL_SLTERM_USERANGLES(spar,O1,nv,LUA,ns)
                     end do
                  end do ; end do
               ENDIF
            else
               IF ( DO_SL_ISOTROPIC ) THEN
                  do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                     g = na_offset(ns,nv) + na
                     LSSL_SLTERM(1,g,spar) = LSSL_SLTERM_ISOTROPIC(spar,1,ns)
                  end do ; end do ; end do 
               ELSE
                  do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                     g = na_offset(ns,nv) + na
                     do o1 = 1,nstokes
                        LSSL_SLTERM(O1,g,spar) = LSSL_SLTERM_USERANGLES(spar,O1,nv,na,ns)
                     end do
                  end do ; end do ; end do 
               ENDIF
            endif
         enddo
      endif

!  Call VFO_LPS_MASTER
!  ===================

!  Upgraded to FO Version 1.5, 7/11/16
!   4/9/19. CUMTRANS and LP_CUMTRANS now additional output. add water-leaving control.

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> MSST situations: Add LOSTRANS_UP/DN, LC_LOSTRANS_UP/DN, THETA_ALL, ALPHA to output list
!    ==> Last 3 lines reorganized output list.

!  1/31/21. Version 2.8.3. Other changes
!     -- Add Doublet flag to input list, add geometry offsets as inputs
!     -- Add do_Polarized_Emissivity flag as input.

      CALL VFO_LPS_MASTER &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, max_partlayers, maxfine,        & ! Input max dims
         maxmoments_input, max_user_levels, max_atmoswfs, max_surfacewfs, max_sleavewfs, & ! Input max dims
         do_solar_sources, do_sunlight, do_thermal_emission, do_surface_emission,        & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_fmatrix, do_obsgeom, do_doublet, do_deltam,      & ! Input flags (general)
         do_lambertian, do_surface_leaving, do_water_leaving, do_Polarized_Emissivity,   & ! Input flags (surface)
         do_Partlayers, do_planpar, do_enhanced_ps, do_profilewfs, do_surfacewfs, do_sleavewfs,     & ! Input flags (geoms/Lin)
         nstokes, ngeoms, nszas, nvzas, nazms, nlayers, nfine, nmoments_input,                      & ! Input Numbers
         nd_offset, nv_offset, na_offset, n_user_levels, LevelMask_up, LevelMask_dn,                & ! Inputs (offsets & cont-levs)
         n_partlayers, partlayers_outindex, partlayers_outflag, partlayers_layeridx,                & ! Inputs (control-partial)
         n_reflecwfs, n_sleavewfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,                  & ! Input Lin control
         dtr, Pie, doCrit, Acrit, eradius, heights, partial_heights,                                & ! Input general
         obsgeom_boa, theta_boa, alpha_boa, phi_boa, flux, fluxvec,                                 & ! Input geometry/flux
         extinction, deltaus, omega, truncfac, greekmat, fmatrix_up, fmatrix_dn,                    & ! Input atmos optical
         bb_input, surfbb, emiss, LS_emiss, reflec, slterm, LS_reflec, LSSL_slterm,                 & ! Input thermal/surf optical
         L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat, fo_L_fmatrix_up, fo_L_fmatrix_dn,& ! Input Lin atmos optical
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,                                  & ! Output - Stokes
         fo_profilewf_ss, fo_profilewf_db, fo_profilewf_dta, fo_profilewf_dts,                      & ! Output - Profile Jacobians
         fo_surfacewf_db, fo_surfacewf_dts,                                                         & ! Output - Surface Jacobians
         fo_stokes_atmos, fo_stokes_surf, fo_stokes,                                                & ! Output - Stokes composites
         fo_profilewf_atmos,  fo_profilewf_surf, fo_profilewf, fo_surfacewf,                        & ! Output - Jacobian composites
         cumtrans, lostrans_up, lostrans_dn, theta_all, alpha,                                      & ! Output - Auxiliary
         LP_cumtrans, LP_lostrans_up, LP_lostrans_dn,                                               & ! Output - Auxiliary
         Fail, message, trace_1, trace_2 )                                                            ! Output - Exception handling

!  Done

      RETURN
END SUBROUTINE VFO_LPS_MASTER_INTERFACE

!  finish module

END MODULE vlidort_vfo_lps_interface_m


