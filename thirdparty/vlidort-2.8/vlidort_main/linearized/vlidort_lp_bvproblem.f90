
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
! #  For linearizations with regular BVP                        #
! #                                                             #
! #            LP_BVP_SOLUTION_MASTER                           #
! #            LP_BVP_COLUMN_SETUP (re-named, Version 2.4)      #
! #            LP_BEAMSOLUTION_NEQK (re-named)                  #
! #            LP_BEAMSOLUTION_NNEK (re-named)                  #
! #                                                             #
! #  For linearizations with telescoped boundary value problem  #
! #                                                             #
! #            LP_BVPTEL_SOLUTION_MASTER                        #
! #            LP_BVPTEL_COLUMN_SETUP                           #
! #            LP_BVPTEL_SURFACE_SETUP  (New, Version 2.8)      #
! #                                                             #
! ###############################################################

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO constructions with IF blocks
!     * Generalized surface treatment with telescoped BVP

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LP_BVP_SOLUTION_MASTER, LP_BVPTEL_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component

      MODULE vlidort_lp_bvproblem_m


private :: LP_BVP_COLUMN_SETUP,    LP_BEAMSOLUTION_NEQK, LP_BEAMSOLUTION_NNEK, &
           LP_BVPTEL_COLUMN_SETUP, LP_BVPTEL_SURFACE_SETUP
public  :: LP_BVP_SOLUTION_MASTER, LP_BVPTEL_SOLUTION_MASTER

      CONTAINS

      SUBROUTINE LP_BVP_SOLUTION_MASTER ( &
         DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTRF,             & ! Flags
         DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                  & ! Flags
         DO_CLASSICAL_SOLUTION, DO_PLANE_PARALLEL, DO_LAYER_SCATTERING, FOURIER,     & ! Flags/Fourier
         IBEAM, NSTOKES, NSTREAMS, NLAYERS, LAYER_TO_VARY, N_LAYERWFS, TAYLOR_ORDER, & ! Numbers (Basic)
         NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,     & ! Numbers (derived)
         MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR, QUAD_STRMWTS,             & ! Bookkeeping
         DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, RF_DIRECT_BEAM, SLTERM,           & ! Optical and direct beam
         BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Input, Beam parameterization
         LP_TRANS_ATMOS, LP_INITIAL_TRANS,  LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,      & ! Input, Linearized Beam parameterization
         CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
         SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC,             & ! Input, Homogeneous/Classical
         BANDMAT2, IPIVOT, SMAT2, SIPIVOT, ALBEDO, BRDF_F,                 & ! Input, BVP matrices, surface
         L_T_DELT_EIGEN, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG,               & ! Input, Linearized Homog solutions
         L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER, LP_BVEC,      & ! Input, Linearized Greens/Thermal/BVEC
         NCON, PCON, L_WLOWER, L_WUPPER,                                   & ! output - Linearized Constants + PI
         STATUS, MESSAGE, TRACE )                                            ! Exception handling

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LC_BVP_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,   &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE VLIDORT_LPC_BVPROBLEM_m, Only : L_BVP_BACKSUB

      IMPLICIT NONE

!  INPUTS
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTRF

      LOGICAL, INTENT (IN) ::          DO_WATER_LEAVING     
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  1/31/21. Version 2.8.3.  Classical_Solution flag (control Greens)

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Basic control numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYERWFS

!  1/31/21. Version 2.8.3. Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)  ::          TAYLOR_ORDER

!  other numbers

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG

!  Layer scattering

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  4/9/19. RF_DIRECT_BEAM is the reflected beam (excludes SLTERM)

      DOUBLE PRECISION, INTENT (IN) :: RF_DIRECT_BEAM   ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SLTERM           ( MAXSTREAMS, MAXSTOKES )

!  optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Solution inputs
!  ---------------


!  Beam parameterization and particular integral vector

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Classical solution

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  homogeneous RTE solutions and inegration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Introduced Green's function solution arrays.

      DOUBLE PRECISION, INTENT (IN) :: ATERM_SAVE (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: BTERM_SAVE (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: CFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: DFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_DN   (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_UP   (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: GAMMA_M    (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GAMMA_P    (MAXEVALUES,MAXLAYERS)

!  BVP matrices

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2   ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

!  Surface inputs
!    -- 1/31/21. Version 2.8.3. BRDF_F defined locally (remove MAXMOMENTS dimension)

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  Linearized solution inputs
!  --------------------------

!  Linearized PI vector, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!  4/9/19 Add linearization of TRANS_ATMOS for the waterleaving contribution

      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_TRANS_ATMOS    ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous solutions
!    -- 1/31/21. Version 2.8.3. Introduced L_KEIGEN

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN       ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3. Introduced Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      DOUBLE PRECISION, INTENT (IN) :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized thermal particular integerals

      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)

!  Linearized Classical solution vector

      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Output
!  ------

!  output linearized integration constants

      DOUBLE PRECISION, INTENT (OUT) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  output particular integrals

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  boundary condition flags

      LOGICAL :: MODIFIED_BCL3, MODIFIED_BCL4

!  error tracing variables

      INTEGER :: STATUS_SUB

!  helper variables

      DOUBLE PRECISION :: COL2_WF  ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  module status and message initialization

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Linearization of the regular BVP case
!  =====================================

!  Set up the column vectors for profile linearizations
!  ----------------------------------------------------

!  Bulk: Compute the main column B' where AX = B'
! 4/9/19. Add water-leaving control, reflected direct beam, surface leaving linearization contribution

!  Profile: Boundary condition flags for special cases
!  Profile: Compute the main column B' where AX = B'

!  Boundary condition flags for special cases

      MODIFIED_BCL3 = ( LAYER_TO_VARY .EQ. 1 )
      MODIFIED_BCL4 = ( LAYER_TO_VARY .EQ. NLAYERS )

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LP_BVP_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component

      CALL LP_BVP_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTRF,               & ! Input, Flags
        DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_PLANE_PARALLEL, & ! Input, Flags
        DO_CLASSICAL_SOLUTION, DO_LAYER_SCATTERING, MODIFIED_BCL3, MODIFIED_BCL4,     & ! Input, Flags/indices
        FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, LAYER_TO_VARY, N_LAYERWFS,        & ! Input, Basic numbers
        TAYLOR_ORDER, NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,               & ! Input, Other numbers
        DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, MUELLER_INDEX, K_REAL, K_COMPLEX,   & ! Input, Optical/Bookkeeping
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM,         & ! Input, Surface Stuff
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                     & ! Beam parameterization
        LP_TRANS_ATMOS, LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Linearized Beam parameterization
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC, LP_BVEC,                & ! Input, Homogeneous/Classical
        CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,   & ! Input, Greens Function
        L_SOLA_XPOS, L_SOLB_XNEG, L_KEIGEN, L_T_DELT_EIGEN,                           & ! Input, Linearized Homogeneous
        L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,                           & ! Input, Linearized Greens/thermal
        L_WLOWER, L_WUPPER, COL2_WF, SCOL2_WF )                                         ! PI Solutions and Column vectors

!  Back-substitution

      CALL L_BVP_BACKSUB ( &
        LAYER_TO_VARY, N_LAYERWFS, NLAYERS, NTOTAL,  & ! Numbers
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS,          & ! Numbers
        NSTKS_NSTRMS_2, K_REAL, K_COMPLEX,           & ! Numbers
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Input BVP matrices
        COL2_WF, SCOL2_WF,                           & ! Input BVP vectors
        NCON, PCON, STATUS_SUB, MESSAGE, TRACE )       ! Output and excpetion handling

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_BVP_SOLUTION_MASTER

!

      SUBROUTINE LP_BVP_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM,             & ! Input, Flags
        DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_PLANE_PARALLEL, & ! Input, Flags
        DO_CLASSICAL_SOLUTION, DO_LAYER_SCATTERING, MODIFIED_BCL3, MODIFIED_BCL4,     & ! Input, Flags/indices
        FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, LAYER_TO_VARY, N_LAYERWFS,        & ! Input, Basic numbers
        TAYLOR_ORDER, NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,               & ! Input, Other numbers
        DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, MUELLER_INDEX, K_REAL, K_COMPLEX,   & ! Input, Optical/Bookkeeping
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM,         & ! Input, Surface Stuff
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                     & ! Beam parameterization
        LP_TRANS_ATMOS, LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Linearized Beam parameterization
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC, LP_BVEC,                & ! Input, Homogeneous/Classical
        CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,   & ! Input, Greens Function
        L_SOLA_XPOS, L_SOLB_XNEG, L_KEIGEN, L_T_DELT_EIGEN,                           & ! Input, Linearized Homogeneous
        L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,                           & ! Input, Linearized Greens/thermal
        L_WLOWER, L_WUPPER, COL2_WF, SCOL2_WF )                                         ! PI Solutions and Column vectors

!  Linearized column-vector setup (profile linearization)
! 4/9/19. Add water-leaving control, reflected direct beam, surface leaving linearization contribution
  
!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LC_BVP_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component
  
!  module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,   &
                                 MAXBANDTOTAL, MAXTOTAL, ZERO

      USE VLIDORT_LPC_BVPROBLEM_m, Only : L_BVP_SURFACE_SETUP

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_WATER_LEAVING
   
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  1/31/21. Version 2.8.3.  New flag for controlling solutions

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

! BVP control and  Layer scattering

      LOGICAL, INTENT (IN) ::          MODIFIED_BCL3, MODIFIED_BCL4
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)  ::          TAYLOR_ORDER

!  Basic control numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYERWFS

!  other numbers

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  optical
!    -- 1/31/21. Version 2.8.3. Add DELTAU_VERT (needed for Green's function Taylor calls)

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  Solution inputs
!  ---------------

!  Beam parameterization
!    -- 1/31/21. Version 2.8.3. Add AVERAGE_SECANT

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Classical solution

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  homogeneous RTE solutions and integration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Introduced Green's function solution arrays.

      DOUBLE PRECISION, INTENT (IN) :: ATERM_SAVE (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: BTERM_SAVE (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: CFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: DFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_DN   (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_UP   (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: GAMMA_M    (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GAMMA_P    (MAXEVALUES,MAXLAYERS)

!  Surface inputs
!    1/31/21. Version 2.8.3. BRDF Fourier inputs are defined locally (remove MAXMOMENTS)

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  4/9/19. RF_DIRECT_BEAM is the reflected beam (excludes SLTERM)

      DOUBLE PRECISION, INTENT (IN) :: RF_DIRECT_BEAM   ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SLTERM           ( MAXSTREAMS, MAXSTOKES )

!  Linearized solution inputs
!  --------------------------

!  Linearized PI Vector, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!  4/9/19 Add linearization of TRANS_ATMOS for the waterleaving contribution
!    -- 1/31/21. Version 2.8.3. Added LP_AVERAGE_SECANT

      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_TRANS_ATMOS    ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN       ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized thermal particular integerals

      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)

!  1/31/21. Version 2.8.3. Introduced Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      DOUBLE PRECISION, INTENT (IN) :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Classical solution vector

      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

!  output particular integrals

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output Column vectors

      DOUBLE PRECISION, INTENT (OUT) :: COL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  local variables
!  ---------------

      LOGICAL ::          REGULAR_BCL3, REGULAR_BCL4
      INTEGER ::          Q, N, N1, I, I1, IR, IROW, CM, C0, O1, NT, NV
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: CPOS, CNEG, L_HOM_R, L_HOM_CR, L_BEAM, FAC
      DOUBLE PRECISION :: T1, T2, T1R, T1I, T2R, T2I, LPTERM

!  reflected arrays

      DOUBLE PRECISION :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  initialize
!  ----------

!  zero the results vectors

!  Enhancement # 1, 6/27/16
      COL2_WF(1:NTOTAL,1:MAX_ATMOSWFS) = ZERO

!  Layer to vary

      NV = LAYER_TO_VARY

!  Copy already existing thermal linearizations
!    This is a very important zeroing.................!!!!!

!  Enhancement # 2, 6/27/16
      L_WUPPER(1:NSTREAMS_2, 1:NSTOKES, 1:NLAYERS, 1:N_LAYERWFS) = ZERO
      L_WLOWER(1:NSTREAMS_2, 1:NSTOKES, 1:NLAYERS, 1:N_LAYERWFS) = ZERO
      IF ( DO_INCLUDE_THERMEMISS ) THEN
        L_WUPPER(1:NSTREAMS_2, 1, 1:NLAYERS, 1:N_LAYERWFS) = L_T_WUPPER(1:NSTREAMS_2, 1:NLAYERS, 1:N_LAYERWFS)
        L_WLOWER(1:NSTREAMS_2, 1, 1:NLAYERS, 1:N_LAYERWFS) = L_T_WLOWER(1:NSTREAMS_2, 1:NLAYERS, 1:N_LAYERWFS)
        L_WUPPER(1:NSTREAMS_2, 2:NSTOKES, 1:NLAYERS, 1:N_LAYERWFS) = ZERO
        L_WLOWER(1:NSTREAMS_2, 2:NSTOKES, 1:NLAYERS, 1:N_LAYERWFS) = ZERO
      ENDIF

!  Get the linearized beam solution for the varying layer
!    -- 1/31/21. Version 2.8.3. Additional arguments to handle Green's function treatment :--
!         DO_CLASSICAL_SOLUTION, TAYLOR_ORDER, AVERAGE_SECANT, LP_AVERAGE_SECANT, NSTKS_NSTRMS, DELTAU_VERT, L_DELTAU_VERT
!         CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL LP_BEAMSOLUTION_NEQK ( &
             DO_LAYER_SCATTERING, DO_CLASSICAL_SOLUTION, TAYLOR_ORDER,      & ! Input, Flags
             FOURIER, IBEAM, NSTOKES, NSTREAMS, LAYER_TO_VARY, N_LAYERWFS,  & ! Input, Numbers
             NSTREAMS_2, NSTKS_NSTRMS, DELTAU_VERT, L_DELTAU_VERT,          & ! Input, optical and bookkeeping
             BEAM_CUTOFF,INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,       & ! Input, Beam Parameterization
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Beam Parameterization (Linearized)
             K_REAL, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,     & ! Input, Homogeneous/Classical
             L_KEIGEN, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,            & ! Input, Homogeneous (Linearizied)
             CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,            & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,            & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                             ! Output solutions

!write(*,*)
!write(*,*) 'L_WUPPER(1,1,1,1) = ',L_WUPPER(1,1,1,1)

      ENDIF

!  complete boundary condition flags

      REGULAR_BCL3 = .NOT. MODIFIED_BCL3
      REGULAR_BCL4 = .NOT. MODIFIED_BCL4

!  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
!  ------------------------------------------------------------------

      N = 1

!    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)

      IF ( MODIFIED_BCL3 ) THEN

!  Start loops

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYERWFS

!  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 3, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R = sum ( ( LCON(1:K_REAL(N),N) * L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) ) + &
!               ( MCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) + &
!                 L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N),N) ) ) )
            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                   L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

          ENDDO
         ENDDO
        ENDDO

!  No variation case (BCL1)

      ELSE

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
!  Enhancement # 4, 6/27/16
          COL2_WF(IROW,1:N_LAYERWFS) = ZERO
         ENDDO
        ENDDO

      ENDIF

!  BCL2 Intermediate levels between top layer and varying layer
!  ------------------------------------------------------------

!  [not required if top layer is varying, case MODIFIED_BCL3 above]

      IF ( REGULAR_BCL3 ) THEN

!  .. nothing varying in these layers

        DO N = 2, LAYER_TO_VARY - 1
          N1 = N - 1
          C0 = N1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
!  Enhancement # 5, 6/27/16
              COL2_WF(CM,1:N_LAYERWFS) = ZERO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  BCL3 - regular upper boundary condition for layer that is varying
!  -----------------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN

!  offsets

        N = LAYER_TO_VARY
        N1  = N - 1
        C0  = N1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYERWFS

!  Beam contribution

            L_BEAM  = + L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 6, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R =sum( ( LCON(1:K_REAL(N),N) *  L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) ) + &
!                ( MCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) + &
!                 L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N),N) ) ) )
            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                   L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  end loops

          ENDDO
         ENDDO
        ENDDO

!  End BCL3 condition

      ENDIF

!  BCL4 - LOWER boundary condition for varying layer
!  -------------------------------------------------

!   special case when layer-to-vary = last (albedo) layer is treated
!   separately below under MODIFIED BCL4.

      IF ( REGULAR_BCL4 ) THEN

!  offsets

        N = LAYER_TO_VARY
        N1  = N + 1
        C0  = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL LP_BEAMSOLUTION_NNEK ( &
             DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION, & ! Input, Flags and order
             NSTOKES, NSTREAMS, N1, NV, N_LAYERWFS, FOURIER, IBEAM,         & ! Input, optical and control
             TAYLOR_ORDER, NSTREAMS_2, NSTKS_NSTRMS, K_REAL, DELTAU_VERT,   & ! Bookkeeping
             BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,      & ! Input, Beam Param.
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Beam Param (Linearized)
             T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,             & ! Input, PIs and linearized  
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,  & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                             ! Output
        ENDIF

!  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYERWFS

!  Beam contributions

            L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 7, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R =sum( ( LCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) + &
!                 L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLA_XPOS(I,O1,1:K_REAL(N),N) ) ) + &
!                 ( MCON(1:K_REAL(N),N) * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) ) )
            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                   L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  End BCL4 condition

      ENDIF

!  BCL5 - Intermediate boundary conditions between varying layer & final
!  ---------------------------------------------------------------------
!mick fix 1/5/2021 - changed N to N1 in call to LP_BEAMSOLUTION_NNEK below

      IF ( REGULAR_BCL4 ) THEN

        DO N = LAYER_TO_VARY + 1, NLAYERS - 1

!  offsets

          N1  = N + 1
          C0  = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Get the linearized beam solution for the next layer

          IF ( DO_SOLAR_SOURCES ) THEN
            CALL LP_BEAMSOLUTION_NNEK ( &
               DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION, & ! Input, Flags and order
               NSTOKES, NSTREAMS, N1, NV, N_LAYERWFS, FOURIER, IBEAM,         & ! Input, optical and control
               TAYLOR_ORDER, NSTREAMS_2, NSTKS_NSTRMS, K_REAL, DELTAU_VERT,   & ! Bookkeeping
               BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,      & ! Input, Beam Param.
               LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Beam Param (Linearized)
               T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,             & ! Input, PIs and linearized  
               GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,  & ! Input, Greens Function
               L_WUPPER, L_WLOWER )                                             ! Output
          ENDIF

!  .. contributions from beam solution (direct assign). No homog. variat

          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
!  Enhancement # 8, 6/27/16
              COL2_WF(CM,1:N_LAYERWFS) = L_WUPPER(I,O1,N1,1:N_LAYERWFS) - L_WLOWER(I,O1,N,1:N_LAYERWFS)
            ENDDO
          ENDDO

!  end layer loop

        ENDDO

!  end BCL5 boundary conditions

      ENDIF

!  Final layer - use BCL6 or BCL4M (last layer is varying)
!  -------------------------------------------------------

      N = NLAYERS

!  Modified BCL4M Component loop

      IF ( MODIFIED_BCL4 ) THEN

!  get the linearized downward-reflected term
!    -- 1/31/21. Version 2.8.3. BRDF_F array defined locally for each Fourier.

        CALL L_BVP_SURFACE_SETUP ( &
          DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BCL4, & ! Flags
          IBEAM, FOURIER, N_LAYERWFS, NSTOKES, NSTREAMS, NLAYERS,   & ! Numbers
          SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,             & ! Surface input
          NSTKS_NSTRMS, MUELLER_INDEX, K_REAL, K_COMPLEX,           & ! bookkeeping
          SOLA_XPOS, T_DELT_EIGEN, L_T_DELT_EIGEN,                  & ! RT Solutions
          L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,                       & ! Linearized RT Solutions
          R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )                           ! Output reflected solutions

!  offsets

        C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  start loops

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYERWFS

!  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 9, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R = sum( ( LCON(1:K_REAL(N),N) * ( ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLA_XPOS(I1,O1,1:K_REAL(N),N,Q) &
!                + L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLA_XPOS(I1,O1,1:K_REAL(N),N) ) - R2_L_HOMP(I,O1,1:K_REAL(N),Q) ) ) + &
!                ( MCON(1:K_REAL(N),N) * ( L_SOLB_XNEG(I1,O1,1:K_REAL(N),N,Q) - R2_L_HOMM(I,O1,1:K_REAL(N),Q) ) ) )
            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                  + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                  - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                  + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contributions

            COL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  ordinary BCL6 Component loop

      ELSE

!  get the linearized downward-reflected term
!    -- 1/31/21. Version 2.8.3. BRDF_F array defined locally for each Fourier.

        CALL L_BVP_SURFACE_SETUP ( &
          DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BCL4, & ! Flags
          IBEAM, FOURIER, N_LAYERWFS, NSTOKES, NSTREAMS, NLAYERS,   & ! Numbers
          SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,             & ! Surface input
          NSTKS_NSTRMS, MUELLER_INDEX, K_REAL, K_COMPLEX,           & ! bookkeeping
          SOLA_XPOS, T_DELT_EIGEN, L_T_DELT_EIGEN,                  & ! RT Solutions
          L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,                       & ! Linearized RT Solutions
          R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )                           ! Output reflected solutions

!  offsets

        C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  start loops : beam contributions only

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
!  Enhancement # 10, 6/27/16
            COL2_WF(CM,1:N_LAYERWFS) = - ( L_WLOWER(I1,O1,N,1:N_LAYERWFS) - R2_L_BEAM(I,O1,1:N_LAYERWFS) )
          ENDDO
        ENDDO

!  End lowest layer boundary value linearization

      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
         DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - RF_DIRECT_BEAM(I,IBEAM,O1) * DELTAU_SLANT(N,LAYER_TO_VARY,IBEAM)
!  Enhancement # 11, 6/27/16
              COL2_WF(CM,1:N_LAYERWFS) = COL2_WF(CM,1:N_LAYERWFS) + ( L_DELTAU_VERT(1:N_LAYERWFS,LAYER_TO_VARY) * FAC )
            ENDDO
         ENDDO
      ENDIF

!  4/9/19 Add the linearization due to adjusted waterleaving term

      IF ( DO_WATER_LEAVING .and. FOURIER .eq. 0 ) THEN
         DO Q = 1, N_LAYERWFS
            LPTERM = LP_TRANS_ATMOS(IBEAM,LAYER_TO_VARY,Q) ; O1 = 1
            DO I = 1, NSTREAMS
               IR = NSTOKES*(I-1) ; I1 = I + NSTREAMS
               IROW = IR + O1     ; CM = C0 + IROW
               COL2_WF(CM,Q) = COL2_WF(CM,Q) + LPTERM * SLTERM(I,O1)
            ENDDO
         ENDDO
      ENDIF

!  copy the single layer vector

      IF ( NLAYERS .EQ. 1 ) THEN
!  Enhancement # 12, 6/27/16
        SCOL2_WF(1:NTOTAL,1:N_LAYERWFS) = COL2_WF(1:NTOTAL,1:N_LAYERWFS)
      ENDIF

!  debug
        !DO NT = 1, NTOTAL
        !  write(95,'(4i4,1p4e17.9)') FOURIER, IBEAM, LAYER_TO_VARY, NT, COL2_WF(NT,1)
        !ENDDO

!  finish

      RETURN
      END SUBROUTINE LP_BVP_COLUMN_SETUP

!

      SUBROUTINE LP_BEAMSOLUTION_NEQK ( &
             DO_LAYER_SCATTERING, DO_CLASSICAL_SOLUTION, TAYLOR_ORDER,      & ! Input, Flags
             FOURIER, IB, NSTOKES, NSTREAMS, N, N_LAYERWFS,                 & ! Input, Numbers
             NSTREAMS_2, NSTKS_NSTRMS, DELTAU_VERT, L_DELTAU_VERT,          & ! Input, optical and bookkeeping
             BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,      & ! Input, Beam Parameterization
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Beam Parameterization (Linearized)
             K_REAL, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,     & ! Input, Homogeneous/Classical
             L_KEIGEN, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,            & ! Input, Homogeneous (Linearizied)
             CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,            & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,            & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                             ! Output solutions

!  Linearization of beam particular integral in one layer.
!  In this module, this is the Layer that contains the variation. N = LAYER_TO_VARY

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops, Rearrange argument lists.

!  1/31/21. Version 2.8.3. Additional arguments to handle Green's function treatment :--
!         DO_CLASSICAL_SOLUTION, TAYLOR_ORDER, AVERAGE_SECANT, LP_AVERAGE_SECANT, NSTKS_NSTRMS, DELTAU_VERT, L_DELTAU_VERT
!         CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE

!  1/31/21. Version 2.8.3. Add 2 more parameters to this list...

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXSTRMSTKS, MAXEVALUES, TAYLOR_SMALL, ONE
       
!  1/31/21. Version 2.8.3. Add Taylor module usage

      USE VLIDORT_Taylor_m, Only : TAYLOR_SERIES_L_1

!  Implicit none
                        
      IMPLICIT NONE

!  Local flags for the solution saving option

      LOGICAL, intent(in)  ::          DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  1/31/21. Version 2.8.3.  New flag for controlling solutions

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  1/31/21. Version 2.8.3. Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)  ::          TAYLOR_ORDER

!  basic control

      INTEGER, INTENT (IN) ::          FOURIER, IB
      INTEGER, INTENT (IN) ::          N, N_LAYERWFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS, NSTREAMS_2, NSTKS_NSTRMS

!  optical 
!    -- 1/31/21. Version 2.8.3. Add DELTAU_VERT (needed for Green's function Taylor calls)

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  bookkeeping
!    -- 1/31/21. Version 2.8.3. Only need real eigensolutions for the Green's function solution

      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )

!  Solution inputs
!  ---------------

!  Beam parameterization
!    -- 1/31/21. Version 2.8.3. Add AVERAGE_SECANT

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Classical solution

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  homogeneous RTE solutions and integration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  1/31/21. Version 2.8.3. Introduced Green's function solution arrays.

      DOUBLE PRECISION, INTENT (IN) :: ATERM_SAVE (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: BTERM_SAVE (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: CFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: DFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_DN   (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_UP   (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: GAMMA_M    (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GAMMA_P    (MAXEVALUES,MAXLAYERS)

!  Linearized solution inputs
!  --------------------------

!  Linearized PI Vector, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!    -- 1/31/21. Version 2.8.3. Added LP_AVERAGE_SECANT

      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized homogeneous solutions

      DOUBLE PRECISION, intent(in)  :: L_KEIGEN       ( MAXEVALUES,   MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES,   MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized thermal particular integerals

!      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
!      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)

!  1/31/21. Version 2.8.3. Introduced Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      DOUBLE PRECISION, INTENT (IN) :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Classical solution vector

      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

!  output particular integrals

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

!  1/31/21. Version 2.8.3. Many new help variables associated with Green's function linearizations

!  Local linearized Green's function multipliers

      DOUBLE PRECISION :: L_GFUNC_UP(MAXEVALUES,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_GFUNC_DN(MAXEVALUES,MAX_ATMOSWFS)

!  Local linearized particular integrals at the layer boundaries

      DOUBLE PRECISION :: LBSOL_LOWER(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION :: LBSOL_UPPER(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS)

!  Help variables

      INTEGER           :: I, I1, K, O1, Q
      DOUBLE PRECISION  :: CONST, WDEL, ZDEL, ZWDEL, AST, BST, EPS, DELTA, LAM, MULT
      DOUBLE PRECISION  :: L_WDEL, L_ZDEL, L_LAM, L_KEG, L_DELTA, L_AST, L_BST
      DOUBLE PRECISION  :: L_MULT(MAX_ATMOSWFS)

!  No linearized particular solution beyond the cutoff layer. ALSO--
!  Nothing if the solution saving mode is on and layer is inactive
!    Exit (Solutions have already been zeroed).

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N) .OR. (N .GT.BEAM_CUTOFF(IB))) THEN
        RETURN
      ENDIF

!  Classical solution
!  ==================

!  1/31/21. Version 2.8.3. Classical solution now an option

      IF ( DO_CLASSICAL_SOLUTION ) THEN

!  Very simple, same code for all situations. Enhancement # 13, 6/27/16

        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)
        DO Q = 1, N_LAYERWFS
          LBSOL_UPPER(1:NSTREAMS_2,1:NSTOKES,Q) = CONST * LP_BVEC(1:NSTREAMS_2,1:NSTOKES,N,N,Q)
          !LBSOL_LOWER(1:NSTREAMS_2,1:NSTOKES,Q) = WDEL  * LBSOL_LOWER(1:NSTREAMS_2,1:NSTOKES,Q) + & 
          !        LP_T_DELT_MUBAR(N,N,IB,Q) * CONST * BVEC(1:NSTREAMS_2,1:NSTOKES,N)
          LBSOL_LOWER(1:NSTREAMS_2,1:NSTOKES,Q) = WDEL  * LBSOL_UPPER(1:NSTREAMS_2,1:NSTOKES,Q) + & 
                  LP_T_DELT_MUBAR(N,N,IB,Q) * CONST * BVEC(1:NSTREAMS_2,1:NSTOKES,N)
        ENDDO

! Add to existing solution. Enhancement # 14, 6/27/16
!mick fix 1/5/2021 - commented out these two extra lines; this accumulation previously moved
!                    to after the Green's function section to apply to both methods 

        !L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,1:N_LAYERWFS) = &
        !  L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,1:N_LAYERWFS) + LBSOL_UPPER(1:NSTREAMS_2,1:NSTOKES,1:N_LAYERWFS)
        !L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,1:N_LAYERWFS) = &
        !  L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,1:N_LAYERWFS) + LBSOL_LOWER(1:NSTREAMS_2,1:NSTOKES,1:N_LAYERWFS)

      ENDIF

!  Green's function solution
!  =========================

!  1/31/21. Version 2.8.3. Completely new, modeled after LIDORT code.

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Set up linearizations of GAMMA constants
!  ----------------------------------------

!  Distinguish two cases:
!  ..(a) quasi-spherical for n > 1
!  ..(b) plane-parallel or QS for n=1

!  Linearizations of optical depth integrations
!  Linearized Green function multipliers

        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)

!  Start Loop over eigenfunctions

        DO K = 1, K_REAL(N)

          ZDEL  = T_DELT_EIGEN(K,N)
          ZWDEL = ZDEL * WDEL
          AST   = CONST * ATERM_SAVE(K,N) 
          BST   = CONST * BTERM_SAVE(K,N) 

!  Downwelling
!  -----------

! .... Make allowances for Taylor series

          IF ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
            EPS   = GAMMA_M(K,N)
            DELTA = DELTAU_VERT(N)
            LAM   = AVERAGE_SECANT(N,IB)
            DO Q = 1, N_LAYERWFS
              L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
              L_KEG  = L_KEIGEN(K,N,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KEG, L_LAM, WDEL, LAM, L_MULT(Q) )
            ENDDO
          ELSE
            MULT = ( ZDEL - WDEL ) / GAMMA_M(K,N)
            DO Q = 1, N_LAYERWFS
              L_ZDEL =  L_T_DELT_EIGEN  (K,N,Q) ; L_KEG  = L_KEIGEN(K,N,Q)
              L_WDEL =  LP_T_DELT_MUBAR (N,N,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
              L_MULT(Q) = ( ( L_ZDEL - L_WDEL ) - MULT * (L_LAM - L_KEG) ) / GAMMA_M(K,N)
            ENDDO
          ENDIF

!  Rob Change 1/15/21. Replace this do-loop with non-logarithmic L_ATERM_SAVE derivatives
!    -- Necessary to use CFUNC = CONST * MULT, as MULT is not available from Taylor clause

!          DO Q = 1, N_LAYERWFS
!            L_AST =  LP_INITIAL_TRANS(N,N,IB,Q)  + L_ATERM_SAVE(K,N,Q)
!            L_GFUNC_DN(K,Q) = GFUNC_DN(K,N) * L_AST + L_MULT(Q) * AST
!          ENDDO

          DO Q = 1, N_LAYERWFS
            L_AST = GFUNC_DN(K,N) * LP_INITIAL_TRANS(N,N,IB,Q)
            L_GFUNC_DN(K,Q) = L_AST + L_MULT(Q) * AST + L_ATERM_SAVE(K,N,Q) * CFUNC(K,N)
          ENDDO

!  Upwelling
!  ---------

          MULT = ( ONE - ZWDEL ) / GAMMA_P(K,N)
          DO Q = 1, N_LAYERWFS
            L_ZDEL =  L_T_DELT_EIGEN  (K,N,Q)   ; L_KEG  = L_KEIGEN(K,N,Q)
            L_WDEL =  LP_T_DELT_MUBAR (N,N,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
            L_MULT(Q) = - ( L_ZDEL*WDEL + L_WDEL*ZDEL + MULT*(L_LAM+L_KEG) ) / GAMMA_P(K,N)
          ENDDO

!  Rob Change 1/15/21. Replace this do-loop with non-logarithmic L_BTERM_SAVE derivatives
!    -- Better to use DFUNC = CONST * MULT
!          DO Q = 1, N_LAYERWFS
!            L_BST =  LP_INITIAL_TRANS(N,N,IB,Q)  + L_BTERM_SAVE(K,N,Q)
!            L_GFUNC_UP(K,Q) = GFUNC_UP(K,N) * L_BST + L_MULT(Q) * BST
!          ENDDO

          DO Q = 1, N_LAYERWFS
            L_BST = GFUNC_UP(K,N) * LP_INITIAL_TRANS(N,N,IB,Q)
            L_GFUNC_UP(K,Q) = L_BST + L_MULT(Q) * BST + L_BTERM_SAVE(K,N,Q) * DFUNC(K,N)
          ENDDO

!  End Eigensolution loop

        ENDDO

!  Now set linearized particular integrals at lower and upper boundaries.
!   Stokes 3/4 with the negative sign.

        DO O1 = 1, NSTOKES
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, N_LAYERWFS
              LBSOL_UPPER(I, O1,Q) = DOT_PRODUCT ( L_GFUNC_UP(1:NSTKS_NSTRMS,Q),   SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N) )  &
                                   + DOT_PRODUCT (   GFUNC_UP(1:NSTKS_NSTRMS,N), L_SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N,Q) )
              LBSOL_UPPER(I1,O1,Q) = DOT_PRODUCT ( L_GFUNC_UP(1:NSTKS_NSTRMS,Q),   SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N) )  &
                                   + DOT_PRODUCT (   GFUNC_UP(1:NSTKS_NSTRMS,N), L_SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N,Q) )
              LBSOL_LOWER(I1,O1,Q) = DOT_PRODUCT ( L_GFUNC_DN(1:NSTKS_NSTRMS,Q),   SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N) )  &
                                   + DOT_PRODUCT (   GFUNC_DN(1:NSTKS_NSTRMS,N), L_SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N,Q) )
              LBSOL_LOWER(I, O1,Q) = DOT_PRODUCT ( L_GFUNC_DN(1:NSTKS_NSTRMS,Q),   SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N) )  &
                                   + DOT_PRODUCT (   GFUNC_DN(1:NSTKS_NSTRMS,N), L_SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N,Q) )
              IF ( O1.gt.2) then
                LBSOL_LOWER(I1,O1,Q) = - LBSOL_LOWER(I1,O1,Q)
                LBSOL_UPPER(I1,O1,Q) = - LBSOL_UPPER(I1,O1,Q)
              ENDIF
            ENDDO
          ENDDO
        ENDDO

!  End Green's function clause

      ENDIF

!  Add to existing solution (both solution methods)
!mick mod 1/5/2021 - replaced original loops with the enhanced lines from the classical section 

      !DO O1 = 1, NSTOKES
      !  DO Q = 1, N_LAYERWFS
      !    L_WUPPER(1:NSTREAMS_2,O1,N,Q) = L_WUPPER(1:NSTREAMS_2,O1,N,Q) + LBSOL_UPPER(1:NSTREAMS_2,O1,Q)
      !    L_WLOWER(1:NSTREAMS_2,O1,N,Q) = L_WLOWER(1:NSTREAMS_2,O1,N,Q) + LBSOL_LOWER(1:NSTREAMS_2,O1,Q)
      !  ENDDO
      !ENDDO
      L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,1:N_LAYERWFS) = L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,1:N_LAYERWFS) &
                                                      + LBSOL_UPPER(1:NSTREAMS_2,1:NSTOKES,1:N_LAYERWFS)
      L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,1:N_LAYERWFS) = L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,1:N_LAYERWFS) &
                                                      + LBSOL_LOWER(1:NSTREAMS_2,1:NSTOKES,1:N_LAYERWFS)

!write(*,*) 'LBSOL_UPPER(1,1,1)       = ',LBSOL_UPPER(1,1,1)
!write(*,*) 'L_WUPPER(1,1,N,1) after  = ',L_WUPPER(1,1,N,1)

!  Finish

      RETURN
      END SUBROUTINE LP_BEAMSOLUTION_NEQK

!

      SUBROUTINE LP_BEAMSOLUTION_NNEK ( &
             DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION, & ! Input, Flags and order
             NSTOKES, NSTREAMS, N, NV, NV_PARAMETERS, FOURIER, IB,          & ! Input, optical and control
             TAYLOR_ORDER, NSTREAMS_2, NSTKS_NSTRMS, K_REAL, DELTAU_VERT,   & ! Bookkeeping
             BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,      & ! Input, Beam Param.
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Beam Param (Linearized)
             T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,             & ! Input, PIs and linearized  
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,  & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                             ! Output

!  Linearization of beam particular integral in one layer.
!  In this module, This is for a layer N not equal to varying layer NV.

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops, Rearrange argument lists.

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, MAXSTREAMS_2
              
!  1/31/21. Version 2.8.3. Additional arguments to handle Green's function treatment :--
!         DO_CLASSICAL_SOLUTION, TAYLOR_ORDER, AVERAGE_SECANT, LP_AVERAGE_SECANT, NSTKS_NSTRMS, DELTAU_VERT, L_DELTAU_VERT
!         CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE

!  1/31/21. Version 2.8.3. Add several  more parameters to this list...

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXSTRMSTKS, MAXEVALUES, TAYLOR_SMALL, ZERO, ONE
       
!  1/31/21. Version 2.8.3. Add Taylor module usage

      USE VLIDORT_Taylor_m, Only : TAYLOR_SERIES_L_1

!  Implicit none
            
      IMPLICIT NONE

!  basic control

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  Local flags for the solution saving option

      LOGICAL, intent(in)  ::          DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  1/31/21. Version 2.8.3.  New flag for controlling solutions

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  1/31/21. Version 2.8.3. Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)  ::          TAYLOR_ORDER

!  basic control

      INTEGER, INTENT (IN) ::          FOURIER, IB
      INTEGER, INTENT (IN) ::          N, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS, NSTREAMS_2, NSTKS_NSTRMS

!  optical 
!    -- 1/31/21. Version 2.8.3. Add DELTAU_VERT (needed for Green's function Taylor calls)

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )

!  bookkeeping
!    -- 1/31/21. Version 2.8.3. Only need real eigensolutions for the Green's function solution

      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )

!  Solution inputs
!  ---------------

!  Beam parameterization
!    -- 1/31/21. Version 2.8.3. Add AVERAGE_SECANT

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Classical solution

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  homogeneous RTE solutions and integration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
!      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
!      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Introduced Green's function solution arrays.

      DOUBLE PRECISION, INTENT (IN) :: ATERM_SAVE (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: BTERM_SAVE (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: GFUNC_DN   (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_UP   (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: GAMMA_M    (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GAMMA_P    (MAXEVALUES,MAXLAYERS)

!  Linearized solution inputs
!  --------------------------

!  Linearized PI Vector, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!    -- 1/31/21. Version 2.8.3. Added LP_AVERAGE_SECANT

      DOUBLE PRECISION, intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized thermal particular integerals
!      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
!      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)

!  Linearized Classical solution vector

      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

!  output particular integrals

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

!  Local linearized Green's function multipliers

      DOUBLE PRECISION  :: L_GFUNC_UP(MAXEVALUES,MAX_ATMOSWFS)
      DOUBLE PRECISION  :: L_GFUNC_DN(MAXEVALUES,MAX_ATMOSWFS)

!  Local linearized particular integrals at the layer boundaries

      DOUBLE PRECISION :: LBSOL_LOWER(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION :: LBSOL_UPPER(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS)

!  Help variables

      INTEGER ::          I, I1, K, O1, Q
      DOUBLE PRECISION :: CONST, ZDEL, WDEL, ZWDEL, TRANS2, AST, BST, T1
      DOUBLE PRECISION :: EPS, DELTA, LAM, MULT, L_WDEL, L_LAM, L_AST, L_BST
      DOUBLE PRECISION :: L_MULT(MAX_ATMOSWFS)
      DOUBLE PRECISION :: VAR_U(MAXSTREAMS_2,MAXSTOKES)

!  No linearized particular solution beyond the cutoff layer. ALSO--
!  Nothing if the solution saving mode is on and layer is inactive
!   -- Solutions are pre-initialized

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N) .OR. (N .GT.BEAM_CUTOFF(IB))) THEN
        RETURN
      ENDIF

!  Classical solution
!  ==================

      IF ( DO_CLASSICAL_SOLUTION ) THEN

!  initialise

        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)
        TRANS2  = CONST * WDEL

!  Distinguish two cases
!  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)
!  ..(b) plane-parallel and qs for n = 1

        IF ( .NOT. DO_PLANE_PARALLEL  ) THEN
!  Enhancement # 15, 6/27/16
          DO Q = 1, NV_PARAMETERS
            VAR_U(1:NSTREAMS_2,1:NSTOKES) = LP_INITIAL_TRANS(N,NV,IB,Q) * BVEC(1:NSTREAMS_2,1:NSTOKES,N) &
                                                + LP_BVEC(1:NSTREAMS_2,1:NSTOKES,N,NV,Q)
            L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,Q) = CONST  * VAR_U(1:NSTREAMS_2,1:NSTOKES)
            L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,Q) = TRANS2 * VAR_U(1:NSTREAMS_2,1:NSTOKES) + & 
                LP_T_DELT_MUBAR(N,NV,IB,Q) * CONST * BVEC(1:NSTREAMS_2,1:NSTOKES,N)
          ENDDO
        ELSE
!  Enhancement # 16, 6/27/16
          DO Q = 1, NV_PARAMETERS
            L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,Q) = LP_INITIAL_TRANS(N,NV,IB,Q) * CONST * BVEC(1:NSTREAMS_2,1:NSTOKES,N)
            L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,Q) = L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,Q) * WDEL
          ENDDO
        ENDIF

!  End classical solution

      ENDIF

!  Green's function solution
!  =========================

      IF ( .not. DO_CLASSICAL_SOLUTION ) THEN

!  Linearizations of optical depth integrations (Linearized Green function multipliers)
!  ------------------------------------------------------------------------------------

!  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)

        IF ( .NOT.DO_PLANE_PARALLEL ) THEN

!  Set up linearizations of GAMMA constants

!  Mick note 9/25/2013 - LP_GAMMA_P & LP_GAMMA_M were previously
!  normalized by GAMMA_P & GAMMA_M, respectively.  With the new
!  definitions of GAMMA_P & GAMMA_M (now straight sums and differences
!  instead of their corresponding reciprocals), we define LP_GAMMA_M as
!  follows:

          CONST   = INITIAL_TRANS(N,IB)
          WDEL    = T_DELT_MUBAR(N,IB)

!  Start discrete ordinate loop

          DO K = 1, K_REAL(N)

            ZDEL  = T_DELT_EIGEN(K,N)
            ZWDEL = ZDEL * WDEL
            AST   = CONST * ATERM_SAVE(K,N) 
            BST   = CONST * BTERM_SAVE(K,N) 

!mick fix 9/4/2013 - small numbers analysis added

!  Downwelling, Make allowances for Taylor series

            IF ( ABS(GAMMA_M(K,N)) .LT. TAYLOR_SMALL ) THEN
              EPS   = GAMMA_M(K,N)
              DELTA = DELTAU_VERT(N)
              LAM   = AVERAGE_SECANT(N,IB)
              DO Q = 1, NV_PARAMETERS
                L_LAM  = LP_AVERAGE_SECANT(N,NV,IB,Q)
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, ZERO, ZERO, L_LAM, WDEL, LAM, L_MULT(Q) )
              ENDDO
            ELSE
              MULT = ( ZDEL - WDEL ) / GAMMA_M(K,N)
              DO Q = 1, NV_PARAMETERS
                L_WDEL =  LP_T_DELT_MUBAR (N,NV,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,NV,IB,Q)
                L_MULT(Q) = ( - L_WDEL - MULT * L_LAM ) / GAMMA_M(K,N)
              ENDDO
            ENDIF
            DO Q = 1, NV_PARAMETERS
              L_AST =  LP_INITIAL_TRANS(N,NV,IB,Q) 
              L_GFUNC_DN(K,Q) = GFUNC_DN(K,N) * L_AST + L_MULT(Q) * AST
            ENDDO

!  Upwelling

            MULT = ( ONE - ZWDEL ) / GAMMA_P(K,N)
            DO Q = 1, NV_PARAMETERS
              L_WDEL =  LP_T_DELT_MUBAR (N,NV,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,NV,IB,Q)
              L_MULT(Q) = - ( L_WDEL*ZDEL + MULT*L_LAM ) / GAMMA_P(K,N)
            ENDDO
            DO Q = 1, NV_PARAMETERS
              L_BST =  LP_INITIAL_TRANS(N,NV,IB,Q)
              L_GFUNC_UP(K,Q) = GFUNC_UP(K,N) * L_BST + L_MULT(Q) * BST
            ENDDO

!  End eigenloop

          ENDDO

!  ..(b) plane-parallel and qs for n = 1

        ELSE

!mick eff 3/22/2017
!  Rob 3/1/21. Bug in code here, only filled up to NSTREAMS (scalar LIDORT hangover!!!)

          DO Q = 1, NV_PARAMETERS
            T1 = LP_INITIAL_TRANS(N,NV,IB,Q)
            !L_GFUNC_DN(1:NSTREAMS,Q) = GFUNC_DN(1:NSTREAMS,N) * T1
            !L_GFUNC_UP(1:NSTREAMS,Q) = GFUNC_UP(1:NSTREAMS,N) * T1
            DO K = 1, K_REAL(N)
              L_GFUNC_DN(K,Q) = GFUNC_DN(K,N) * T1
              L_GFUNC_UP(K,Q) = GFUNC_UP(K,N) * T1
            ENDDO
          ENDDO

!  End clause (b)

        ENDIF

!  Now set linearized particular integrals at lower and upper boundaries.

        DO O1 = 1, NSTOKES
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, NV_PARAMETERS
              LBSOL_UPPER(I, O1,Q) = DOT_PRODUCT ( L_GFUNC_UP(1:NSTKS_NSTRMS,Q),   SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N) )
              LBSOL_UPPER(I1,O1,Q) = DOT_PRODUCT ( L_GFUNC_UP(1:NSTKS_NSTRMS,Q),   SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N) )
              LBSOL_LOWER(I1,O1,Q) = DOT_PRODUCT ( L_GFUNC_DN(1:NSTKS_NSTRMS,Q),   SOLB_XNEG(I,O1,1:NSTKS_NSTRMS,N) )
              LBSOL_LOWER(I, O1,Q) = DOT_PRODUCT ( L_GFUNC_DN(1:NSTKS_NSTRMS,Q),   SOLA_XPOS(I,O1,1:NSTKS_NSTRMS,N) )
              IF ( O1.gt.2) then
                LBSOL_LOWER(I1,O1,Q) = - LBSOL_LOWER(I1,O1,Q)
                LBSOL_UPPER(I1,O1,Q) = - LBSOL_UPPER(I1,O1,Q)
              ENDIF
            ENDDO
          ENDDO
        ENDDO

!  Add to existing solution

        DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
            L_WUPPER(1:NSTREAMS_2,O1,N,Q) = L_WUPPER(1:NSTREAMS_2,O1,N,Q) + LBSOL_UPPER(1:NSTREAMS_2,O1,Q)
            L_WLOWER(1:NSTREAMS_2,O1,N,Q) = L_WLOWER(1:NSTREAMS_2,O1,N,Q) + LBSOL_LOWER(1:NSTREAMS_2,O1,Q)
          ENDDO
        ENDDO

!  End Green's function clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_BEAMSOLUTION_NNEK

!

      SUBROUTINE LP_BVPTEL_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM,              & ! Flags
        DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION, DO_LAYER_SCATTERING, FOURIER, IBEAM, & ! Flags/Indices
        TAYLOR_ORDER, NSTOKES, NSTREAMS, NLAYERS, LAYER_TO_VARY, N_LAYERWFS,           & ! Basic Control Numbers
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NLAYERS_TEL, ACTIVE_LAYERS,          & ! Other Numbers
        N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,              & ! BVPTel Control
        MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR, DIRECT_BEAM,                 & ! Bookkeeping/Surface
        DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Optical.Direct/Disords
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                      & ! Beam parameterization
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC, WLOWER,          & ! Homogeneous Solutions
        QUAD_STRMWTS, ALBEDO, BRDF_F,                                          & ! Surface inputs
        CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                    & ! Input, Greens Function
        ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,                    & ! Input, Greens Function 
        L_T_DELT_EIGEN, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG, LP_BVEC,           & ! Linearized Homogeneous/Classicak
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,                  & ! Linearized Beam parameterization
        BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,                                & ! BVPTel matrices
        L_WLOWER, L_WUPPER, NCON, PCON,                                        & ! output solutions
        STATUS, MESSAGE, TRACE )                                                 ! Exception handling

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LP_BVPTEL_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,   &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO

      USE LAPACK_TOOLS_m, Only : DGBTRS, DGETRS

      USE VLIDORT_LPC_BVPROBLEM_m, Only : L_BVP_BACKSUB

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  1/31/21. Version 2.8.3.  Classical_Solution flag (control Greens)

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  1/31/21. Version 2.8.3. Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)  ::          TAYLOR_ORDER

!  Basic control numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYERWFS

!  other numbers


      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  Layer scattering

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Telescoped BVP Control

      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (IN) ::          NLAYERS_TEL
      INTEGER, INTENT (IN) ::          ACTIVE_LAYERS ( MAXLAYERS )

!  optical and direct beam

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM   ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  discrete ordinate stuff

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR

!  homogeneous RTE solutions and inegration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Beam parameterization

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Particular integral

      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Classical solution

      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  1/31/21. Version 2.8.3. Introduced Green's function solution arrays.

      DOUBLE PRECISION, INTENT (IN) :: ATERM_SAVE (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: BTERM_SAVE (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: CFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: DFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_DN   (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_UP   (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: GAMMA_M    (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GAMMA_P    (MAXEVALUES,MAXLAYERS)

!  Surface inputs
!    -- 1/31/21. Version 2.8.3. BRDF_F defined locally (remove MAXMOMENTS dimension)

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Linearized inputs
!  -----------------

!  Linearized homogeneous solutions
!    -- 1/31/21. Version 2.8.3. Introduced L_KEIGEN

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN       ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3. Introduced Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      DOUBLE PRECISION, intent (IN) :: L_ATERM_SAVE (MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent (IN) :: L_BTERM_SAVE (MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Beam and PI vector
!    -- 1/31/21. Version 2.8.3. Introduced LP_AVERAGE_SECANT

      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  BVP matrices

      DOUBLE PRECISION, INTENT (IN) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOTTEL   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2   ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

!  output linearized integration constants

      DOUBLE PRECISION, INTENT (OUT) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  output particular integrals

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  error tracing variables

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI

!  Other local help variables

      INTEGER ::          I, I1, N, N1, NAL, NAL1, O1, IC, ICOW, Q
      INTEGER ::          K, KO1, K0, K1, K2, C0, NS
      INTEGER ::          IR, IROW, IROW1, IROW_S, IROW1_S
      DOUBLE PRECISION :: SPAR, SHOM, L_HOM1, L_HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, L_HOM1CR, L_HOM2CR
      DOUBLE PRECISION :: LXR, MXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION :: LXR1, MXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION :: LXR2, MXR2, NXR2, PXR2, LLXR2, MLXR2

!  Arrays for linearized column vectors

      DOUBLE PRECISION :: COLTEL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  Initialise

      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  What is this ?
      Q = 0

!  Profile: Boundary condition flags for special cases
!  Profile: Compute the main column B' where AX = B'

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LC_BVP_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component

      CALL LP_BVPTEL_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,     & ! Input, Flags
        DO_CLASSICAL_SOLUTION, DO_LAYER_SCATTERING, FOURIER, IBEAM, NSTOKES,    & ! Input, Flags/indices
        NSTREAMS, NLAYERS, LAYER_TO_VARY, N_LAYERWFS, TAYLOR_ORDER, NSTREAMS_2,  & ! Input, Basic numbers
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,         & ! Input, Tel Control
        DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, MUELLER_INDEX,                & ! Input, Optical/Bookkeeping
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F, DIRECT_BEAM,              & ! Input, Surface Stuff
        T_DELT_DISORDS, L_T_DELT_DISORDS, K_REAL, K_COMPLEX,                    & ! Discrete ordinate transmittances
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,               & ! Beam parameterization
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,                   & ! Linearized Beam parameterization
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC, WLOWER,           & ! Input, Homogeneous/Classical
        CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                     & ! Input, Greens Function
        ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,                     & ! Input, Greens Function
        L_SOLA_XPOS, L_SOLB_XNEG, L_KEIGEN, L_T_DELT_EIGEN, LP_BVEC,            & ! Input, Linearized Homogeneous/Classical
        L_WLOWER, L_WUPPER, COLTEL2_WF, SCOL2_WF )                                ! PI Solutions and Column vectors

!  Solve linearized BVP: Several Active layers
!  ===========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  BV solution for linearized integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUBDIAG, &
              N_BVTELMATRIX_SUPDIAG, N_LAYERWFS, &
              BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, &
              COLTEL2_WF, MAXTOTAL, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (NLAYERS>1)in LP_BVPTEL_SOLUTION_MASTER'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set linearized integration constants, active layers

        C0 = - NSTKS_NSTRMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTKS_NSTRMS_2

!  set real constants from the solution vector

!  Enhancement # 17, 6/27/16
          NCON(1:K_REAL(N),N,1:N_LAYERWFS) = COLTEL2_WF(C0+1:C0+K_REAL(N), 1:N_LAYERWFS)
          PCON(1:K_REAL(N),N,1:N_LAYERWFS) = COLTEL2_WF(C0+NSTKS_NSTRMS+1:C0+NSTKS_NSTRMS+K_REAL(N), 1:N_LAYERWFS)

!  set complex constants from the solution vector

          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = K + K_REAL(N) + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
!  Enhancement # 18, 6/27/16
            NCON(K1,N,1:N_LAYERWFS) = COLTEL2_WF(C0+IROW,    1:N_LAYERWFS)
            NCON(K2,N,1:N_LAYERWFS) = COLTEL2_WF(C0+IROW_S,  1:N_LAYERWFS)
            PCON(K1,N,1:N_LAYERWFS) = COLTEL2_WF(C0+IROW1,   1:N_LAYERWFS)
            PCON(K2,N,1:N_LAYERWFS) = COLTEL2_WF(C0+IROW1_S, 1:N_LAYERWFS)
          ENDDO

!  End number of telescoped layers

        ENDDO

!  Solve linearized BVP: Single Layer only
!  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

        NAL = ACTIVE_LAYERS(1)

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS &
           ( 'N', NSTKS_NSTRMS_2, N_LAYERWFS, &
              SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              SCOL2_WF, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (NLAYERS=1)in LP_BVPTEL_SOLUTION_MASTER'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  set real constants from the solution vector

!  Enhancement # 19, 6/27/16
        NCON(1:K_REAL(NAL),NAL,1:N_LAYERWFS) = SCOL2_WF(1:K_REAL(NAL),1:N_LAYERWFS)
        PCON(1:K_REAL(NAL),NAL,1:N_LAYERWFS) = SCOL2_WF(NSTKS_NSTRMS+1:NSTKS_NSTRMS+K_REAL(NAL),1:N_LAYERWFS)

!  set complex constants from the solution vector

        KO1 = K_REAL(NAL) + 1
        DO K = 1, K_COMPLEX(NAL)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(NAL)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = IROW + K_COMPLEX(NAL)
          IROW1_S = IROW_S + NSTKS_NSTRMS
!  Enhancement # 20, 6/27/16
          NCON(K1,NAL,1:N_LAYERWFS) = SCOL2_WF(IROW,    1:N_LAYERWFS)
          NCON(K2,NAL,1:N_LAYERWFS) = SCOL2_WF(IROW_S,  1:N_LAYERWFS)
          PCON(K1,NAL,1:N_LAYERWFS) = SCOL2_WF(IROW1,   1:N_LAYERWFS)
          PCON(K2,NAL,1:N_LAYERWFS) = SCOL2_WF(IROW1_S, 1:N_LAYERWFS)
        ENDDO

!  End solution of telescoped BVP

      ENDIF

!  Set linearized integration constants for non-active layers
!  ==========================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required if the first active
!       layer is the varying layer

      NAL = ACTIVE_LAYERS(1)
      IF ( NAL .GT. 1 ) THEN

        N1 = NAL - 1

!  Case 1. If active layer is also the varying layer
!  -------------------------------------------------

        IF ( LAYER_TO_VARY.EQ.NAL ) THEN

!  start stream, stokes loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = ( I1 - 1 ) * NSTOKES
            IC = ( I - 1  ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              ICOW = IC + O1

!  start parameter loop

              DO Q = 1, N_LAYERWFS

!  real homogeneous solutions
!  Enhancement # 21, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!             SHOM_R = sum( ( NCON(1:K_REAL(NAL),NAL,Q) * SOLA_XPOS(I1,O1,1:K_REAL(NAL),NAL) + & 
!                       LCON(1:K_REAL(NAL),NAL)   * L_SOLA_XPOS(I1,O1,1:K_REAL(NAL),NAL,Q) ) + &
!                     ( T_DELT_EIGEN(1:K_REAL(NAL),NAL) * & 
!                       ( PCON(1:K_REAL(NAL),NAL,Q) * SOLB_XNEG(I1,O1,1:K_REAL(NAL),NAL) + & 
!                       MCON(1:K_REAL(NAL),NAL) * L_SOLB_XNEG(I1,O1,1:K_REAL(NAL),NAL,Q) ) + &
!                       L_T_DELT_EIGEN(1:K_REAL(NAL),NAL,Q) * MCON(1:K_REAL(NAL),NAL) * & 
!                       SOLB_XNEG(I1,O1,1:K_REAL(NAL),NAL) ) )
                SHOM_R = ZERO
                DO K = 1, K_REAL(NAL)
                  NXR  = NCON(K,NAL,Q) *   SOLA_XPOS(I1,O1,K,NAL)
                  PXR  = PCON(K,NAL,Q) *   SOLB_XNEG(I1,O1,K,NAL)
                  MXR  = MCON(K,NAL)   *   SOLB_XNEG(I1,O1,K,NAL)
                  LLXR = LCON(K,NAL)   * L_SOLA_XPOS(I1,O1,K,NAL,Q)
                  MLXR = MCON(K,NAL)   * L_SOLB_XNEG(I1,O1,K,NAL,Q)
                  L_HOM1 = NXR + LLXR
                  L_HOM2 =   T_DELT_EIGEN(K,NAL)   * ( PXR + MLXR ) + &
                           L_T_DELT_EIGEN(K,NAL,Q) *   MXR
                  SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                ENDDO

!  complex homogeneous solutions

                SHOM_CR = ZERO
                KO1 = K_REAL(NAL) + 1
                DO K = 1, K_COMPLEX(NAL)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  NXR1  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I1,O1,K1,NAL) &
                          - NCON(K2,NAL,Q) *   SOLA_XPOS(I1,O1,K2,NAL)
                  PXR1  =   PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL) &
                          - PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL)
                  PXR2  =   PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL) &
                          + PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL)
                  MXR1  =   MCON(K1,NAL) *   SOLB_XNEG(I1,O1,K1,NAL) &
                          - MCON(K2,NAL) *   SOLB_XNEG(I1,O1,K2,NAL)
                  MXR2  =   MCON(K1,NAL) *   SOLB_XNEG(I1,O1,K2,NAL) &
                          + MCON(K2,NAL) *   SOLB_XNEG(I1,O1,K1,NAL)
                  LLXR1  =   LCON(K1,NAL) * L_SOLA_XPOS(I1,O1,K1,NAL,Q) &
                           - LCON(K2,NAL) * L_SOLA_XPOS(I1,O1,K2,NAL,Q)
                  MLXR1  =   MCON(K1,NAL) * L_SOLB_XNEG(I1,O1,K1,NAL,Q) &
                           - MCON(K2,NAL) * L_SOLB_XNEG(I1,O1,K2,NAL,Q)
                  MLXR2  =   MCON(K1,NAL) * L_SOLB_XNEG(I1,O1,K2,NAL,Q) &
                           + MCON(K2,NAL) * L_SOLB_XNEG(I1,O1,K1,NAL,Q)
                  L_HOM1CR = NXR1 + LLXR1
                  L_HOM2CR = + T_DELT_EIGEN(K1,NAL)     * ( PXR1 + MLXR1 ) &
                             - T_DELT_EIGEN(K2,NAL)     * ( PXR2 + MLXR2 ) &
                             + L_T_DELT_EIGEN(K1,NAL,Q) *   MXR1 &
                             - L_T_DELT_EIGEN(K2,NAL,Q) *   MXR2
                  SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

                SHOM = SHOM_R + SHOM_CR
                SPAR = L_WUPPER(I1,O1,NAL,Q)
                PCON(ICOW,N1,Q) = SPAR + SHOM
                NCON(ICOW,N1,Q) = ZERO

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  Case 2. Active layer is below varying layer
!  -------------------------------------------

        ELSE IF ( LAYER_TO_VARY.LT.NAL) THEN

!  start stream, stokes loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = ( I1 - 1 ) * NSTOKES
            IC = ( I - 1  ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              ICOW = IC + O1

!  start parameter loop

              DO Q = 1, N_LAYERWFS

!  real homogeneous solutions
!  Enhancement # 22, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!             SHOM_R = sum( NCON(1:K_REAL(NAL),NAL,Q)*SOLA_XPOS(I1,O1,1:K_REAL(NAL),NAL) + &
!                     ( PCON(1:K_REAL(NAL),NAL,Q)*SOLB_XNEG(I1,O1,1:K_REAL(NAL),NAL) * T_DELT_EIGEN(1:K_REAL(NAL),NAL) ) )
                SHOM_R = ZERO
                DO K = 1, K_REAL(NAL)
                  NXR = NCON(K,NAL,Q)*SOLA_XPOS(I1,O1,K,NAL)
                  PXR = PCON(K,NAL,Q)*SOLB_XNEG(I1,O1,K,NAL)
                  L_HOM1 = NXR
                  L_HOM2 = PXR * T_DELT_EIGEN(K,NAL)
                  SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                ENDDO

!  complex homogeneous solutions

                SHOM_CR = ZERO
                KO1 = K_REAL(NAL) + 1
                DO K = 1, K_COMPLEX(NAL)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  NXR1  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I1,O1,K1,NAL) &
                          - NCON(K2,NAL,Q) *   SOLA_XPOS(I1,O1,K2,NAL)
                  PXR1  =   PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL) &
                          - PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL)
                  PXR2  =   PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL) &
                          + PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL)
                  L_HOM1CR = NXR1
                  L_HOM2CR =  + T_DELT_EIGEN(K1,NAL) * PXR1 &
                              - T_DELT_EIGEN(K2,NAL) * PXR2
                  SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

                SHOM = SHOM_R + SHOM_CR
                SPAR = L_WUPPER(I1,O1,NAL,Q)
                PCON(ICOW,N1,Q) = SPAR + SHOM
                NCON(ICOW,N1,Q) = ZERO

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  Case 3. Active layer is above varying layer
!  -------------------------------------------

        ELSE IF ( LAYER_TO_VARY.GT.NAL) THEN

!  No solutions

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = ( I1 - 1 ) * NSTOKES
            IC = ( I - 1  ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              ICOW = IC + O1
!  Enhancement # 23, 6/27/16
              PCON(ICOW,N1,1:N_LAYERWFS) = ZERO
              NCON(ICOW,N1,1:N_LAYERWFS) = ZERO
            ENDDO
          ENDDO

!  End case 3

        ENDIF

!  End Clause NAL on NAL

      ENDIF

!  For remaining non-active atmospheric layers to TOA, propagate upwards
!   Additional linearizations if you are passing through the varying lay

      DO N = NAL - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
!  Enhancement # 24, 6/27/16
            NCON(IROW,N,1:N_LAYERWFS) = ZERO
            PCON(IROW,N,1:N_LAYERWFS) = T_DELT_DISORDS(I,N1) * PCON(IROW,N1,1:N_LAYERWFS)
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
!  Enhancement # 25, 6/27/16
            PCON(IROW,N,1:N_LAYERWFS) = PCON(IROW,N,1:N_LAYERWFS) + L_T_DELT_DISORDS(I,N1,1:N_LAYERWFS) * MCON(IROW,N1)
          ENDDO
         ENDDO
        ENDIF
      ENDDO

!  Transmittance layers below active layer(s)
!  -----------------------------------------

!       ** Only do this if active scattering is above (not adjacent to) the surface layer

!   -- NCON values are  propagated downwards from bottom of last active layer
!   -- PCON values also propagated downwards, BUT only present if surface condition
!  1.   Require linearized solutions at bottom of last active layer
!  2.   Set values for layer immediately below last active layer
!  3.   Remaining layers to bottom, just propagate using discrete-ordinate transmittances

      NAL = ACTIVE_LAYERS(NLAYERS_TEL) ; NAL1 = NAL + 1
      IF ( NAL .LT. NLAYERS ) THEN

!  N-CONSTANTS always required
!  ^^^^^^^^^^^^^^^^^^^^^^^^^^^

!  Case 1. If active layer equals the varying layer
!  ------------------------------------------------

        IF ( LAYER_TO_VARY .EQ. NAL ) THEN

!  start stream, stokes and parameter loops

          DO I = 1, NSTREAMS
            IR = ( I - 1 ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              DO Q = 1, N_LAYERWFS

!  real homogeneous solutions

                SHOM_R = ZERO
                DO K = 1, K_REAL(NAL)
                  NXR  = NCON(K,NAL,Q) *   SOLA_XPOS(I,O1,K,NAL)
                  PXR  = PCON(K,NAL,Q) *   SOLB_XNEG(I,O1,K,NAL)
                  LXR  = LCON(K,NAL)   *   SOLA_XPOS(I,O1,K,NAL)
                  LLXR = LCON(K,NAL)   * L_SOLA_XPOS(I,O1,K,NAL,Q)
                  MLXR = MCON(K,NAL)   * L_SOLB_XNEG(I,O1,K,NAL,Q)
                  L_HOM2 = PXR + MLXR
                  L_HOM1 =   T_DELT_EIGEN(K,NAL)   * ( NXR + LLXR ) + &
                           L_T_DELT_EIGEN(K,NAL,Q) *   LXR
                  SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                ENDDO
    
!  complex homogeneous solutions

                SHOM_CR = ZERO
                KO1 = K_REAL(NAL) + 1
                DO K = 1, K_COMPLEX(NAL)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  NXR1  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I,O1,K1,NAL) &
                          - NCON(K2,NAL,Q) *   SOLA_XPOS(I,O1,K2,NAL)
                  NXR2  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I,O1,K2,NAL) &
                          + NCON(K2,NAL,Q) *   SOLA_XPOS(I,O1,K1,NAL)
                  PXR1  =   PCON(K1,NAL,Q) *   SOLB_XNEG(I,O1,K1,NAL) &
                          - PCON(K2,NAL,Q) *   SOLB_XNEG(I,O1,K2,NAL)
                  LXR1  =   LCON(K1,NAL) *   SOLA_XPOS(I,O1,K1,NAL) &
                          - LCON(K2,NAL) *   SOLA_XPOS(I,O1,K2,NAL)
                  LXR2  =   LCON(K1,NAL) *   SOLA_XPOS(I,O1,K2,NAL) &
                          + LCON(K2,NAL) *   SOLA_XPOS(I,O1,K1,NAL)
                  LLXR1  =   LCON(K1,NAL) * L_SOLA_XPOS(I,O1,K1,NAL,Q) &
                           - LCON(K2,NAL) * L_SOLA_XPOS(I,O1,K2,NAL,Q)
                  LLXR2  =   LCON(K1,NAL) * L_SOLA_XPOS(I,O1,K2,NAL,Q) &
                           + LCON(K2,NAL) * L_SOLA_XPOS(I,O1,K1,NAL,Q)
                  MLXR1  =   MCON(K1,NAL) * L_SOLB_XNEG(I,O1,K1,NAL,Q) &
                           - MCON(K2,NAL) * L_SOLB_XNEG(I,O1,K2,NAL,Q)
                  L_HOM2CR = PXR1 + MLXR1
                  L_HOM1CR = + T_DELT_EIGEN(K1,NAL)     * ( NXR1 + LLXR1 ) &
                             - T_DELT_EIGEN(K2,NAL)     * ( NXR2 + LLXR2 ) &
                             + L_T_DELT_EIGEN(K1,NAL,Q) *   LXR1 &
                             - L_T_DELT_EIGEN(K2,NAL,Q) *   LXR2
                  SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

                SHOM = SHOM_R + SHOM_CR
                SPAR = L_WLOWER(I,O1,NAL,Q)
                NCON(IROW,NAL1,Q) = SPAR + SHOM

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  Case 2. Active layer is below varying layer
!  -------------------------------------------

        ELSE IF ( LAYER_TO_VARY.LT.NAL) THEN

!  start stream, stokes  and parameter loops

          DO I = 1, NSTREAMS
            IR = ( I - 1 ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              DO Q = 1, N_LAYERWFS

!  real homogeneous solutions

!  Enhancement # 26, 6/27/16.   Not Implemented, lacks clarity
!             SHOM_R = sum( ( NCON(1:K_REAL(NAL),NAL,Q) * SOLA_XPOS(I,O1,1:K_REAL(NAL),NAL) * T_DELT_EIGEN(1:K_REAL(NAL),NAL) ) + &
!                            PCON(1:K_REAL(NAL),NAL,Q) * SOLB_XNEG(I,O1,1:K_REAL(NAL),NAL) )
                SHOM_R = ZERO
                DO K = 1, K_REAL(NAL)
                  NXR = NCON(K,NAL,Q) * SOLA_XPOS(I,O1,K,NAL)
                  PXR = PCON(K,NAL,Q) * SOLB_XNEG(I,O1,K,NAL)
                  L_HOM2 = PXR
                  L_HOM1 = NXR * T_DELT_EIGEN(K,NAL)
                  SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                ENDDO

!  complex homogeneous solutions

                SHOM_CR = ZERO
                KO1 = K_REAL(NAL) + 1
                DO K = 1, K_COMPLEX(NAL)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  NXR1  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I,O1,K1,NAL) &
                          - NCON(K2,NAL,Q) *   SOLA_XPOS(I,O1,K2,NAL)
                  NXR2  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I,O1,K2,NAL) &
                          + NCON(K2,NAL,Q) *   SOLA_XPOS(I,O1,K1,NAL)
                  PXR1  =   PCON(K1,NAL,Q) *   SOLB_XNEG(I,O1,K1,NAL) &
                          - PCON(K2,NAL,Q) *   SOLB_XNEG(I,O1,K2,NAL)
                  L_HOM2CR = PXR1
                  L_HOM1CR =  + T_DELT_EIGEN(K1,NAL) * NXR1 &
                              - T_DELT_EIGEN(K2,NAL) * NXR2
                  SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

                SHOM = SHOM_R + SHOM_CR
                SPAR = L_WLOWER(I,O1,NAL,Q)
                NCON(IROW,NAL1,Q) = SPAR + SHOM

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  Case 3. Active layer is above varying layer
!  -------------------------------------------

        ELSE IF ( LAYER_TO_VARY.GT.NAL) THEN

!  start stream, stokes  and parameter loops

          DO I = 1, NSTREAMS
            IR = ( I - 1 ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              DO Q = 1, N_LAYERWFS

!  real homogeneous solutions

                SHOM_R = ZERO
                DO K = 1, K_REAL(NAL)
                  NXR = NCON(K,NAL,Q) * SOLA_XPOS(I,O1,K,NAL)
                  PXR = PCON(K,NAL,Q) * SOLB_XNEG(I,O1,K,NAL)
                  L_HOM2 = PXR
                  L_HOM1 = NXR * T_DELT_EIGEN(K,NAL)
                  SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                ENDDO

!  complex homogeneous solutions

                SHOM_CR = ZERO
                KO1 = K_REAL(NAL) + 1
                DO K = 1, K_COMPLEX(NAL)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  NXR1  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I,O1,K1,NAL) &
                          - NCON(K2,NAL,Q) *   SOLA_XPOS(I,O1,K2,NAL)
                  NXR2  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I,O1,K2,NAL) &
                          + NCON(K2,NAL,Q) *   SOLA_XPOS(I,O1,K1,NAL)
                  PXR1  =   PCON(K1,NAL,Q) *   SOLB_XNEG(I,O1,K1,NAL) &
                          - PCON(K2,NAL,Q) *   SOLB_XNEG(I,O1,K2,NAL)
                  L_HOM2CR = PXR1
                  L_HOM1CR =  + T_DELT_EIGEN(K1,NAL) * NXR1 &
                              - T_DELT_EIGEN(K2,NAL) * NXR2
                  SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                ENDDO

!  Final

                NCON(IROW,NAL1,Q) = SHOM_R + SHOM_CR

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  End case 1,2,3

        ENDIF

!  other layers to bottom of medium: propagate downwards.
!   Additional variation if you are passing through the varying layer.

        DO N = NAL + 2, NLAYERS
          N1 = N - 1
          IF ( N1 .EQ. LAYER_TO_VARY ) THEN
            DO I = 1, NSTREAMS
              IR = ( I - 1 ) * NSTOKES
              DO O1 = 1, NSTOKES
                IROW = IR + O1
!  Enhancement # 29, 6/27/16
                NCON(IROW,N,1:N_LAYERWFS) = T_DELT_DISORDS(I,N1) * NCON(IROW,N1,1:N_LAYERWFS) &
                                        + L_T_DELT_DISORDS(I,N1,1:N_LAYERWFS) * LCON(IROW,N1)
              ENDDO
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              IR = ( I - 1 ) * NSTOKES
              DO O1 = 1, NSTOKES
                IROW = IR + O1
!  Enhancement # 28, 6/27/16
                NCON(IROW,N,1:N_LAYERWFS) = T_DELT_DISORDS(I,N1) * NCON(IROW,N1,1:N_LAYERWFS)
              ENDDO
            ENDDO
          ENDIF
        ENDDO

!  P-Constants need to be determined if there is a surface condition. Otherwise zero.
!  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        IF ( DO_INCLUDE_SURFACE ) THEN

!  Case 1. Varying layer same as Last telescoped layer
!  ---------------------------------------------------

          IF ( LAYER_TO_VARY .EQ. NAL ) THEN

!  start stream, stokes and parameter loops

            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                DO Q = 1, N_LAYERWFS

!  real homogeneous solutions

                  SHOM_R = ZERO
                  DO K = 1, K_REAL(NAL)
                    NXR  = NCON(K,NAL,Q) *   SOLA_XPOS(I1,O1,K,NAL)
                    PXR  = PCON(K,NAL,Q) *   SOLB_XNEG(I1,O1,K,NAL)
                    MXR  = MCON(K,NAL)   *   SOLB_XNEG(I1,O1,K,NAL)
                    LLXR = LCON(K,NAL)   * L_SOLA_XPOS(I1,O1,K,NAL,Q)
                    MLXR = MCON(K,NAL)   * L_SOLB_XNEG(I1,O1,K,NAL,Q)
                    L_HOM1 = NXR + LLXR
                    L_HOM2 =   T_DELT_EIGEN(K,NAL) * ( PXR + MLXR ) + L_T_DELT_EIGEN(K,NAL,Q) * MXR
                    SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                  ENDDO

!  complex homogeneous solutions

                  SHOM_CR = ZERO
                  KO1 = K_REAL(NAL) + 1
                  DO K = 1, K_COMPLEX(NAL)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    NXR1  = NCON(K1,NAL,Q) *   SOLA_XPOS(I1,O1,K1,NAL) - NCON(K2,NAL,Q) *   SOLA_XPOS(I1,O1,K2,NAL)
                    PXR1  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL) - PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL)
                    PXR2  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL) + PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL)
                    MXR1  = MCON(K1,NAL) *   SOLB_XNEG(I1,O1,K1,NAL) - MCON(K2,NAL) *   SOLB_XNEG(I1,O1,K2,NAL)
                    MXR2  = MCON(K1,NAL) *   SOLB_XNEG(I1,O1,K2,NAL) + MCON(K2,NAL) *   SOLB_XNEG(I1,O1,K1,NAL)
                    LLXR1 = LCON(K1,NAL) * L_SOLA_XPOS(I1,O1,K1,NAL,Q) - LCON(K2,NAL) * L_SOLA_XPOS(I1,O1,K2,NAL,Q)
                    MLXR1 = MCON(K1,NAL) * L_SOLB_XNEG(I1,O1,K1,NAL,Q) - MCON(K2,NAL) * L_SOLB_XNEG(I1,O1,K2,NAL,Q)
                    MLXR2 = MCON(K1,NAL) * L_SOLB_XNEG(I1,O1,K2,NAL,Q) + MCON(K2,NAL) * L_SOLB_XNEG(I1,O1,K1,NAL,Q)
                    L_HOM1CR = NXR1 + LLXR1
                    L_HOM2CR =   T_DELT_EIGEN(K1,NAL)     * ( PXR1 + MLXR1 ) -   T_DELT_EIGEN(K2,NAL)   * ( PXR2 + MLXR2 ) &
                               + L_T_DELT_EIGEN(K1,NAL,Q) *   MXR1           - L_T_DELT_EIGEN(K2,NAL,Q) * MXR2
                    SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                  ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

                  SHOM = SHOM_R + SHOM_CR
                  SPAR = L_WUPPER(I1,O1,NAL,Q)
                  PCON(ICOW,NAL1,Q) = SPAR + SHOM

!  End loops

                ENDDO
              ENDDO
            ENDDO

!  Case 2. Varying layer below Last telescoped layer
!  -------------------------------------------------

          ELSE IF ( LAYER_TO_VARY .LT. NAL ) THEN

!  start stream, stokes and parameter loops

            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                DO Q = 1, N_LAYERWFS

!  real homogeneous solutions

                  SHOM_R = ZERO
                  DO K = 1, K_REAL(NAL)
                    NXR  = NCON(K,NAL,Q) * SOLA_XPOS(I1,O1,K,NAL)
                    PXR  = PCON(K,NAL,Q) * SOLB_XNEG(I1,O1,K,NAL)
                    L_HOM1 = NXR
                    L_HOM2 = T_DELT_EIGEN(K,NAL) * PXR 
                    SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                  ENDDO

!  complex homogeneous solutions

                  SHOM_CR = ZERO
                  KO1 = K_REAL(NAL) + 1
                  DO K = 1, K_COMPLEX(NAL)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    NXR1  = NCON(K1,NAL,Q) *   SOLA_XPOS(I1,O1,K1,NAL) - NCON(K2,NAL,Q) *   SOLA_XPOS(I1,O1,K2,NAL)
                    PXR1  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL) - PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL)
                    PXR2  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL) + PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL)
                    L_HOM1CR = NXR1
                    L_HOM2CR = T_DELT_EIGEN(K1,NAL) * PXR1 - T_DELT_EIGEN(K2,NAL) * PXR2
                    SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                  ENDDO

!  Final

                  SHOM = SHOM_R + SHOM_CR
                  PCON(ICOW,NAL1,Q) = SHOM_R + SHOM_CR

!  End loops

                ENDDO
              ENDDO
            ENDDO

!  Case 3. Varying layer above Last telescoped layer
!  -------------------------------------------------

          ELSE

!  start stream, stokes and parameter loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IC = ( I - 1  ) * NSTOKES
            DO O1 = 1, NSTOKES
              ICOW = IC + O1
              DO Q = 1, N_LAYERWFS

!  real homogeneous solutions

                  SHOM_R = ZERO
                  DO K = 1, K_REAL(NAL)
                    NXR  = NCON(K,NAL,Q) * SOLA_XPOS(I1,O1,K,NAL)
                    PXR  = PCON(K,NAL,Q) * SOLB_XNEG(I1,O1,K,NAL)
                    L_HOM1 = NXR
                    L_HOM2 = T_DELT_EIGEN(K,NAL) * PXR 
                    SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                  ENDDO

!  complex homogeneous solutions

                  SHOM_CR = ZERO
                  KO1 = K_REAL(NAL) + 1
                  DO K = 1, K_COMPLEX(NAL)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    NXR1  = NCON(K1,NAL,Q) *   SOLA_XPOS(I1,O1,K1,NAL) - NCON(K2,NAL,Q) *   SOLA_XPOS(I1,O1,K2,NAL)
                    PXR1  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL) - PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL)
                    PXR2  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL) + PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL)
                    L_HOM1CR = NXR1
                    L_HOM2CR = T_DELT_EIGEN(K1,NAL) * PXR1 - T_DELT_EIGEN(K2,NAL) * PXR2
                    SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                  ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

                  SHOM = SHOM_R + SHOM_CR
                  SPAR = L_WUPPER(I1,O1,NAL,Q)
                  PCON(ICOW,NAL1,Q) = SPAR + SHOM

!  End loops

                ENDDO
              ENDDO
            ENDDO

!  End 3 cases

          ENDIF

!  other constants propagated

          DO N = NAL + 2, NLAYERS
            N1 = N - 1
            IF ( N.EQ.LAYER_TO_VARY ) THEN
              DO I = 1, NSTREAMS
                I1  = I + NSTREAMS
                IC = ( I - 1  ) * NSTOKES
                DO O1 = 1, NSTOKES
                  ICOW = IC + O1
                  DO Q = 1, N_LAYERWFS
                    PCON(ICOW,N,Q) = ( PCON(ICOW,N1,Q) - L_T_DELT_DISORDS(I,N,Q) * MCON(ICOW,N) ) / T_DELT_DISORDS(I,N) 
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO I = 1, NSTREAMS
                I1  = I + NSTREAMS
                IC = ( I - 1  ) * NSTOKES
                DO O1 = 1, NSTOKES
                  ICOW = IC + O1
                  PCON(ICOW,N,1:N_LAYERWFS) = PCON(ICOW,N1,1:N_LAYERWFS) / T_DELT_DISORDS(I,N) 
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  Otherwise all P-constants are zero

        ELSE

          DO N = NAL1, NLAYERS
            DO I = 1, NSTREAMS
              I1  = I + NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                PCON(ICOW,N,1:N_LAYERWFS)   = ZERO
              ENDDO
            ENDDO
          ENDDO

!  End general surface treatment

        ENDIF

!  End clause for non-active layers below telescoped problem

      ENDIF

!  finish

      RETURN
      END SUBROUTINE LP_BVPTEL_SOLUTION_MASTER

!

      SUBROUTINE LP_BVPTEL_SURFACE_SETUP ( &
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, IBEAM, FOURIER,         & ! Flags
            NSTOKES, NSTREAMS, N, NLAYERS, LAYER_TO_VARY, N_LAYERWFS,          & ! Numbers
            SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,                      & ! Surface input
            MUELLER_INDEX, K_REAL, K_COMPLEX,                                  & ! bookkeeping
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input transmittances
            SOLA_XPOS, SOLB_XNEG, WLOWER, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,  & ! RT Solutions
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output reflected solutions

!  Linearized surface reflectance terms, Telescoping

!  1/31/21. Version 2.8.3. BRDF_F array defined locally for each Fourier component. Drop MAXMOMENTS dimension

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAX_ATMOSWFS, MAXLAYERS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTOKES_SQ, ZERO, ONE

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE

!  Numbers

      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS, N
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYERWFS

!  Surface inputs
!  1/31/21. Version 2.8.3. BRDF_F array defined locally for each Fourier component. Drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F (MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  Bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  Solutions

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Linearized solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Discrete ordinate transmittances

      DOUBLE PRECISION, intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output - linearized surface-reflected solutions

      DOUBLE PRECISION, INTENT (OUT) :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  cumulative tranmsittance (and linearization) from bottom of lowest active layer to surface
!     Calculated here, but  could be done earlier and passed in

      DOUBLE PRECISION, intent (out) :: CUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION, intent (out) :: L_CUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      DOUBLE PRECISION  :: QCUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION  :: L_QCUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  help arrays

      DOUBLE PRECISION :: PV_W ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: PS_W ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: HS_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: HS_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  other variables

      LOGICAL ::          MODIFIED_BOUNDARY
      INTEGER ::          I, J, O1, IB, Q, N1, M
      INTEGER ::          K, KO1, K0, K1, K2, MO2(MAXSTOKES)

      DOUBLE PRECISION :: H1R, H1I, H2R, H2I, KMULT
      DOUBLE PRECISION :: H1, H2, HP, HM, L_H1, L_H2, L_H1R, L_H2R, L_H1I, L_H2I
      DOUBLE PRECISION :: BEAM, L_BEAM, HOMP, HOMM, L_HOMP, L_HOMM
      DOUBLE PRECISION :: HOMPR, HOMMR, L_HOMPR, L_HOMMR, HOMPI, HOMMI, L_HOMPI, L_HOMMI
      DOUBLE PRECISION :: H1_S_CR, H2_S_CR, H1_S_CI, H2_S_CI

!  Initial section
!  ---------------

!  Always zero the result to start

      R2_L_BEAM = ZERO
      R2_L_HOMP = ZERO
      R2_L_HOMM = ZERO

      CUMTRANS   = ONE
      L_CUMTRANS = ZERO

!  Beam index and Fourier component

      IB = IBEAM
      M = FOURIER

!  Boundary flag

      MODIFIED_BOUNDARY = .true.

!  Return if no surface contributions

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if Fourier component > 0 (Lambertian)

      IF ( DO_LAMBERTIAN_SURFACE .and. M.gt.0 ) RETURN

!  Cumulative transmittance and its linearization

      CUMTRANS(1:NSTREAMS) = ONE ; L_CUMTRANS(1:NSTREAMS,:) = ZERO
      DO N1 = NLAYERS, N+1, -1
         IF ( N1.eq.LAYER_TO_VARY ) THEN
            DO Q = 1, N_LAYERWFS
               DO J = 1, NSTREAMS
                  L_CUMTRANS(J,Q) = L_CUMTRANS(J,Q) * T_DELT_DISORDS(J,N1) + CUMTRANS(J) * L_T_DELT_DISORDS(J,N1,Q)
               ENDDO
            ENDDO
         ELSE
            DO Q = 1, N_LAYERWFS
               DO J = 1, NSTREAMS
                  L_CUMTRANS(J,Q) = L_CUMTRANS(J,Q) * T_DELT_DISORDS(J,N1)
               ENDDO
            ENDDO
         ENDIF
         DO J = 1, NSTREAMS
            CUMTRANS(J) = CUMTRANS(J) * T_DELT_DISORDS(J,N1)
         ENDDO
      ENDDO

!  stored variables

      QCUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * QUAD_STRMWTS(1:NSTREAMS)
      DO Q = 1, N_LAYERWFS
         L_QCUMTRANS(1:NSTREAMS,Q) = L_CUMTRANS(1:NSTREAMS,Q) * QUAD_STRMWTS(1:NSTREAMS)
      ENDDO

!  Set up Auxiliary arrays
!  -----------------------

!  Particular integral

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          PS_W(J,O1) = WLOWER(J,O1,N) * QCUMTRANS(J)
          DO Q = 1, N_LAYERWFS
            PV_W(J,O1,Q) = L_WLOWER(J,O1,N,Q) * QCUMTRANS(J) + WLOWER(J,O1,N) * L_QCUMTRANS(J,Q)
          ENDDO
        ENDDO
      ENDDO

!    Modified boundary condition: homogeneous parts

      IF ( MODIFIED_BOUNDARY ) THEN

!  If N is the varying layer
!  -------------------------

        if ( N.EQ.LAYER_TO_VARY ) THEN

!  start loops

          DO J = 1, NSTREAMS
            DO O1 = 1, NSTOKES

!  real homogeneous solution contributions

              DO K = 1, K_REAL(N)
                H1 = SOLA_XPOS(J,O1,K,N) * T_DELT_EIGEN(K,N) 
                H2 = SOLB_XNEG(J,O1,K,N)
                HS_P(J,O1,K) = QCUMTRANS(J)*H1
                HS_M(J,O1,K) = QCUMTRANS(J)*H2
                DO Q = 1, N_LAYERWFS
                  L_H1 = L_SOLA_XPOS(J,O1,K,N,Q) *   T_DELT_EIGEN(K,N) + &
                           SOLA_XPOS(J,O1,K,N)   * L_T_DELT_EIGEN(K,N,Q)
                  L_H2 = L_SOLB_XNEG(J,O1,K,N,Q)
                  HV_P(J,O1,K,Q) = QCUMTRANS(J) * L_H1 + L_QCUMTRANS(J,Q) * H1
                  HV_M(J,O1,K,Q) = QCUMTRANS(J) * L_H2 + L_QCUMTRANS(J,Q) * H2
                ENDDO
              ENDDO

!  Complex homogeneous solution contributions

              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                H1R = SOLA_XPOS(J,O1,K1,N) * T_DELT_EIGEN(K1,N) - SOLA_XPOS(J,O1,K2,N) * T_DELT_EIGEN(K2,N)
                H1I = SOLA_XPOS(J,O1,K1,N) * T_DELT_EIGEN(K2,N) + SOLA_XPOS(J,O1,K2,N) * T_DELT_EIGEN(K1,N)
                H2R = SOLB_XNEG(J,O1,K1,N)
                H2I = SOLB_XNEG(J,O1,K2,N)
                HS_P(J,O1,K1) = QCUMTRANS(J) * H1R
                HS_P(J,O1,K2) = QCUMTRANS(J) * H1I
                HS_M(J,O1,K1) = QCUMTRANS(J) * H2R
                HS_M(J,O1,K2) = QCUMTRANS(J) * H1R
                DO Q = 1, N_LAYERWFS
                  L_H1R = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K1,N)   &
                        - L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K2,N)   &
                        +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K1,N,Q) &
                        -   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K2,N,Q)
                  L_H1I = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K2,N)   &
                        + L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K1,N)   &
                        +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K2,N,Q) &
                        +   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K1,N,Q)
                  L_H2R = L_SOLB_XNEG(J,O1,K1,N,Q)
                  L_H2I = L_SOLB_XNEG(J,O1,K2,N,Q)
                  HV_P(J,O1,K1,Q) = QCUMTRANS(J) * L_H1R + L_QCUMTRANS(J,Q) * H1R
                  HV_P(J,O1,K2,Q) = QCUMTRANS(J) * L_H1I + L_QCUMTRANS(J,Q) * H1I
                  HV_M(J,O1,K1,Q) = QCUMTRANS(J) * L_H2R + L_QCUMTRANS(J,Q) * H2R
                  HV_M(J,O1,K2,Q) = QCUMTRANS(J) * L_H2I + L_QCUMTRANS(J,Q) * H2I
                ENDDO
              ENDDO

!  End loops

            ENDDO
          ENDDO

!  If N is NOT the varying layer
!  -----------------------------

        ELSE

!  start loops

          DO J = 1, NSTREAMS
            DO O1 = 1, NSTOKES

!  real homogeneous solution contributions

              DO K = 1, K_REAL(N)
                H1 = SOLA_XPOS(J,O1,K,N) * T_DELT_EIGEN(K,N) 
                H2 = SOLB_XNEG(J,O1,K,N)
                HS_P(J,O1,K) = QCUMTRANS(J)*H1
                HS_M(J,O1,K) = QCUMTRANS(J)*H2
                DO Q = 1, N_LAYERWFS
                  HV_P(J,O1,K,Q) = L_QCUMTRANS(J,Q) * H1
                  HV_M(J,O1,K,Q) =  L_QCUMTRANS(J,Q) * H2
                ENDDO
              ENDDO

!  Complex homogeneous solution contributions

              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                H1R = SOLA_XPOS(J,O1,K1,N) * T_DELT_EIGEN(K1,N) - SOLA_XPOS(J,O1,K2,N) * T_DELT_EIGEN(K2,N)
                H1I = SOLA_XPOS(J,O1,K1,N) * T_DELT_EIGEN(K2,N) + SOLA_XPOS(J,O1,K2,N) * T_DELT_EIGEN(K1,N)
                H2R = SOLB_XNEG(J,O1,K1,N)
                H2I = SOLB_XNEG(J,O1,K2,N)
                HS_P(J,O1,K1) = QCUMTRANS(J) * H1R
                HS_P(J,O1,K2) = QCUMTRANS(J) * H1I
                HS_M(J,O1,K1) = QCUMTRANS(J) * H2R
                HS_M(J,O1,K2) = QCUMTRANS(J) * H1R
                DO Q = 1, N_LAYERWFS
                  HV_P(J,O1,K1,Q) = L_QCUMTRANS(J,Q) * H1R
                  HV_P(J,O1,K2,Q) = L_QCUMTRANS(J,Q) * H1I
                  HV_M(J,O1,K1,Q) = L_QCUMTRANS(J,Q) * H2R
                  HV_M(J,O1,K2,Q) = L_QCUMTRANS(J,Q) * H2I
                ENDDO
              ENDDO

!  End loops

            ENDDO
          ENDDO

!  End clause N is or isnt Varying layer

        ENDIF

!  End modified boundary condition

      ENDIF

!  Integrated Downward reflection (Lambertian case)
!  ================================================

      if ( DO_LAMBERTIAN_SURFACE ) THEN

!  reflection

        KMULT = SURFACE_FACTOR * ALBEDO

!  only 1 component

        O1 = 1

!  Particular solution (only for the first Stokes component)

        BEAM = SUM (PS_W(1:NSTREAMS,O1) )
        DO Q = 1, N_LAYERWFS
          L_BEAM = SUM(PV_W(1:NSTREAMS,O1,Q))
          R2_L_BEAM(1:NSTREAMS,O1,Q) = KMULT * ( L_BEAM * CUMTRANS(1:NSTREAMS) + BEAM * L_CUMTRANS(1:NSTREAMS,Q) )
        ENDDO

!  Homogeneous solutions for the modified condition

        IF ( MODIFIED_BOUNDARY ) THEN

!  Homogeneous real solutions

          DO K = 1, K_REAL(N)
            HOMP = SUM(HS_P(1:NSTREAMS,O1,K))
            HOMM = SUM(HS_M(1:NSTREAMS,O1,K))
            DO Q = 1, N_LAYERWFS
              L_HOMP = SUM(HV_P(1:NSTREAMS,O1,K,Q))
              L_HOMM = SUM(HV_M(1:NSTREAMS,O1,K,Q))
              R2_L_HOMP(1:NSTREAMS,O1,K,Q) = KMULT * ( L_HOMP * CUMTRANS(1:NSTREAMS) + HOMP * L_CUMTRANS(1:NSTREAMS,Q) )
              R2_L_HOMM(1:NSTREAMS,O1,K,Q) = KMULT * ( L_HOMM * CUMTRANS(1:NSTREAMS) + HOMM * L_CUMTRANS(1:NSTREAMS,Q) )
            ENDDO
          ENDDO

!  Homogeneous Complex solutions

          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            HOMPR = SUM(HV_P(1:NSTREAMS,O1,K1,Q))
            HOMPI = SUM(HV_P(1:NSTREAMS,O1,K2,Q))
            HOMMR = SUM(HV_M(1:NSTREAMS,O1,K1,Q))
            HOMMI = SUM(HV_M(1:NSTREAMS,O1,K2,Q))
            DO Q = 1, N_LAYERWFS
              L_HOMPR = SUM(HV_P(1:NSTREAMS,O1,K1,Q))
              L_HOMPI = SUM(HV_P(1:NSTREAMS,O1,K2,Q))
              L_HOMMR = SUM(HV_M(1:NSTREAMS,O1,K1,Q))
              L_HOMMI = SUM(HV_M(1:NSTREAMS,O1,K2,Q))
              R2_L_HOMP(1:NSTREAMS,O1,K1,Q) = KMULT * ( L_HOMPR * CUMTRANS(1:NSTREAMS) + HOMPR * L_CUMTRANS(1:NSTREAMS,Q) )
              R2_L_HOMP(1:NSTREAMS,O1,K2,Q) = KMULT * ( L_HOMPI * CUMTRANS(1:NSTREAMS) + HOMPI * L_CUMTRANS(1:NSTREAMS,Q) )
              R2_L_HOMM(1:NSTREAMS,O1,K1,Q) = KMULT * ( L_HOMMR * CUMTRANS(1:NSTREAMS) + HOMMR * L_CUMTRANS(1:NSTREAMS,Q) )
              R2_L_HOMM(1:NSTREAMS,O1,K2,Q) = KMULT * ( L_HOMMI * CUMTRANS(1:NSTREAMS) + HOMMI * L_CUMTRANS(1:NSTREAMS,Q) )
            ENDDO
          ENDDO

!  End modified boundary condition clause

        ENDIF

!  BRDF surface condition
!  ======================

!  1/31/21. Version 2.8.3. BRDF_F array defined locally for each Fourier component. 
!    -- Drop FOURIER = M index

      ELSE

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            MO2(1:NSTOKES) = MUELLER_INDEX(O1,1:NSTOKES)

!  Particular solution

            BEAM = ZERO
            DO J = 1, NSTREAMS
              H1 = DOT_PRODUCT(PS_W(J,1:NSTOKES),BRDF_F(MO2(1:NSTOKES),J,I))
              BEAM = BEAM + H1
            ENDDO
            DO Q = 1, N_LAYERWFS
              L_BEAM = ZERO
              DO J = 1, NSTREAMS
                L_H1 = DOT_PRODUCT(PV_W(J,1:NSTOKES,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                L_BEAM = L_BEAM + L_H1
              ENDDO
              R2_L_BEAM(I,O1,Q) = SURFACE_FACTOR * ( L_BEAM * CUMTRANS(I) + BEAM * L_CUMTRANS(I,Q) )
            ENDDO

!  Homogeneous solutions for the modified condition

            IF ( MODIFIED_BOUNDARY ) THEN

!  Homogeneous real solutions

              DO K = 1, K_REAL(N)
                HOMP = ZERO ; HOMM = ZERO
                DO J = 1, NSTREAMS
                  HP = DOT_PRODUCT(HS_P(J,1:NSTOKES,K),BRDF_F(MO2(1:NSTOKES),J,I))
                  HM = DOT_PRODUCT(HS_M(J,1:NSTOKES,K),BRDF_F(MO2(1:NSTOKES),J,I))
                  HOMP = HOMP + HP ; HOMM = HOMM + HM
                ENDDO
                DO Q = 1, N_LAYERWFS
                  L_HOMP = ZERO ; L_HOMM = ZERO
                  DO J = 1, NSTREAMS
                    HP = DOT_PRODUCT(HV_P(J,1:NSTOKES,K,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    HM = DOT_PRODUCT(HV_M(J,1:NSTOKES,K,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    L_HOMP = L_HOMP + HP ; L_HOMM = L_HOMM + HM
                  ENDDO
                  R2_L_HOMP(I,O1,K,Q) = SURFACE_FACTOR * ( L_HOMP * CUMTRANS(I) + HOMP * L_CUMTRANS(I,Q) )
                  R2_L_HOMM(I,O1,K,Q) = SURFACE_FACTOR * ( L_HOMM * CUMTRANS(I) + HOMM * L_CUMTRANS(I,Q) )
                ENDDO
              ENDDO

!  homogeneous complex solutions

              KO1 = K_REAL(NLAYERS) + 1
              DO K = 1, K_COMPLEX(NLAYERS)
                K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                HOMPR = ZERO ; HOMMR = ZERO ; HOMPI = ZERO ; HOMMI = ZERO
                DO J = 1, NSTREAMS
                  H1_S_CR = DOT_PRODUCT(HS_P(J,1:NSTOKES,K1),BRDF_F(MO2(1:NSTOKES),J,I))
                  H2_S_CR = DOT_PRODUCT(HS_M(J,1:NSTOKES,K1),BRDF_F(MO2(1:NSTOKES),J,I))
                  H1_S_CI = DOT_PRODUCT(HS_P(J,1:NSTOKES,K2),BRDF_F(MO2(1:NSTOKES),J,I))
                  H2_S_CI = DOT_PRODUCT(HS_M(J,1:NSTOKES,K2),BRDF_F(MO2(1:NSTOKES),J,I))
                  HOMPR = HOMPR + H1_S_CR
                  HOMMR = HOMMR + H2_S_CR
                  HOMPI = HOMPI + H1_S_CI
                  HOMMI = HOMMI + H2_S_CI
                ENDDO
                DO Q = 1, N_LAYERWFS
                  L_HOMPR = ZERO ; L_HOMMR = ZERO ; L_HOMPI = ZERO ; L_HOMMI = ZERO
                  H1_S_CR = DOT_PRODUCT(HV_P(J,1:NSTOKES,K1,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                  H2_S_CR = DOT_PRODUCT(HV_M(J,1:NSTOKES,K1,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                  H1_S_CI = DOT_PRODUCT(HV_P(J,1:NSTOKES,K2,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                  H2_S_CI = DOT_PRODUCT(HV_M(J,1:NSTOKES,K2,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                  L_HOMPR = L_HOMPR + H1_S_CR
                  L_HOMMR = L_HOMMR + H2_S_CR
                  L_HOMPI = L_HOMPI + H1_S_CI
                  L_HOMMI = L_HOMMI + H2_S_CI
                  R2_L_HOMP(I,O1,K1,Q) = SURFACE_FACTOR * ( L_HOMPR * CUMTRANS(I) + HOMPR * L_CUMTRANS(I,Q) )
                  R2_L_HOMM(I,O1,K1,Q) = SURFACE_FACTOR * ( L_HOMMR * CUMTRANS(I) + HOMMR * L_CUMTRANS(I,Q) )
                  R2_L_HOMP(I,O1,K2,Q) = SURFACE_FACTOR * ( L_HOMPI * CUMTRANS(I) + HOMPI * L_CUMTRANS(I,Q) )
                  R2_L_HOMM(I,O1,K2,Q) = SURFACE_FACTOR * ( L_HOMMI * CUMTRANS(I) + HOMMI * L_CUMTRANS(I,Q) )
                ENDDO
              ENDDO

!  End modified condition

            ENDIF

!  end stream and Stokes loops

          ENDDO
        ENDDO

!  End clause BRDF vs. Lambertian

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_BVPTEL_SURFACE_SETUP

!

      SUBROUTINE LP_BVPTEL_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,     & ! Input, Flags
        DO_CLASSICAL_SOLUTION, DO_LAYER_SCATTERING, FOURIER, IBEAM, NSTOKES,          & ! Input, Flags/indices
        NSTREAMS, NLAYERS, LAYER_TO_VARY, N_LAYERWFS, TAYLOR_ORDER, NSTREAMS_2,       & ! Input, Basic numbers
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS, & ! Input, Tel Control
        DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, MUELLER_INDEX,                & ! Input, Optical/Bookkeeping
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F, DIRECT_BEAM,              & ! Input, Surface Stuff
        T_DELT_DISORDS, L_T_DELT_DISORDS, K_REAL, K_COMPLEX,                    & ! Discrete ordinate transmittances
        BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,               & ! Beam parameterization
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,                   & ! Linearized Beam parameterization
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC, WLOWER,           & ! Input, Homogeneous/Classical
        CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                     & ! Input, Greens Function
        ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,                     & ! Input, Greens Function
        L_SOLA_XPOS, L_SOLB_XNEG, L_KEIGEN, L_T_DELT_EIGEN, LP_BVEC,            & ! Input, Linearized Homogeneous/Classical
        L_WLOWER, L_WUPPER, COLTEL2_WF, SCOL2_WF )                                ! PI Solutions and Column vectors

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LC_BVP_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LP_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, MAXMOMENTS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ, MAXTOTAL, ZERO

      IMPLICIT NONE

!  Flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  1/31/21. Version 2.8.3.  New flag for controlling solutions

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Layer scattering flag

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)  ::          TAYLOR_ORDER

!  Numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYERWFS

!  More numbers

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  BVP Tel control

      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::          NLAYERS_TEL
      INTEGER, INTENT (IN) ::          ACTIVE_LAYERS ( MAXLAYERS )

!  optical 
!    -- 1/31/21. Version 2.8.3. Add DELTAU_VERT (needed for Green's function Taylor calls)

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR

!  Surface inputs
!    1/31/21. Version 2.8.3. BRDF Fourier inputs are defined locally (remove MAXMOMENTS)

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  4/9/19. DIRECT_BEAM is the reflected beam 

      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM   ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  Solution inputs
!  ---------------

!  Beam parameterization
!    -- 1/31/21. Version 2.8.3. Add AVERAGE_SECANT

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  discrete ordinate stuff

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  homogeneous RTE solutions and inegration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  1/31/21. Version 2.8.3. Introduced Green's function solution arrays.

      DOUBLE PRECISION, INTENT (IN) :: ATERM_SAVE (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: BTERM_SAVE (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: CFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: DFUNC      (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_DN   (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GFUNC_UP   (MAXEVALUES,MAXLAYERS)

      DOUBLE PRECISION, INTENT (IN) :: GAMMA_M    (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: GAMMA_P    (MAXEVALUES,MAXLAYERS)

!  Classical beam-solution vector, and particular integral vector

      DOUBLE PRECISION, INTENT (IN) :: BVEC   ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Linearized solution inputs
!  --------------------------

!  Linearized Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!    -- 1/31/21. Version 2.8.3. Added LP_AVERAGE_SECANT

      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN       ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  1/31/21. Version 2.8.3. Introduced Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      DOUBLE PRECISION, INTENT (IN) :: L_ATERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_BTERM_SAVE(MAXEVALUES,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Classical solution PI vector

      DOUBLE PRECISION, INTENT (IN) :: LP_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

!  output particular integrals

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output Column vectos

      DOUBLE PRECISION, INTENT (OUT) :: COLTEL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          Q, N, I, I1, IR, CM, C0, IROW, O1, N1, NS
      INTEGER ::          K, KO1, K0, K1, K2, STATUS
      DOUBLE PRECISION :: CPOS,CNEG,L_HOM_R,L_HOM_CR, L_BEAM
      DOUBLE PRECISION :: T1,T2,T1R,T1I,T2R,T2I, FAC, BEAM

!  Output arguments from the Surface setup (reflectances and cumulative transmittances)
!  cumulative tranmsittance (and linearization) from bottom of lowest active layer to surface

      DOUBLE PRECISION :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: CUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION :: L_CUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  Try this safety-first zeroing

!  Enhancement # 30, 6/27/16
      L_WUPPER(1:NSTREAMS_2,1:NSTOKES,LAYER_TO_VARY,1:N_LAYERWFS) = ZERO
      L_WLOWER(1:NSTREAMS_2,1:NSTOKES,LAYER_TO_VARY,1:N_LAYERWFS) = ZERO

!  status

      status = 0

!  Get the linearized solutions for the layer that is varying
!    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        IF ( N.EQ.LAYER_TO_VARY ) THEN
          CALL LP_BEAMSOLUTION_NEQK ( &
             DO_LAYER_SCATTERING, DO_CLASSICAL_SOLUTION, TAYLOR_ORDER,      & ! Input, Flags
             FOURIER, IBEAM, NSTOKES, NSTREAMS, LAYER_TO_VARY, N_LAYERWFS,  & ! Input, Numbers
             NSTREAMS_2, NSTKS_NSTRMS, DELTAU_VERT, L_DELTAU_VERT,          & ! Input, optical and bookkeeping
             BEAM_CUTOFF,INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,       & ! Input, Beam Parameterization
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,          & ! Input, Beam Parameterization (Linearized)
             K_REAL, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,     & ! Input, Homogeneous/Classical
             L_KEIGEN, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,            & ! Input, Homogeneous (Linearizied)
             CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,            & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,            & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                             ! Output solutions
        ENDIF
      ENDDO

!  General case. NLAYERS_TEL > 1
!  =============================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  zero column vector
!  Enhancement # 31, 6/27/16
        COLTEL2_WF(1:N_BVTELMATRIX_SIZE, 1:MAX_ATMOSWFS) = ZERO

!  top of first active layer, first boundary condition
!  ---------------------------------------------------

        NS = 1
        N = ACTIVE_LAYERS(NS)
        C0 = 0

!  If this active layer = layer that is varying,
!       then require homogeneous and beam solution linearizations

        IF ( LAYER_TO_VARY .EQ. N ) THEN

!  start loops

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_LAYERWFS

!  Beam contribution

                L_BEAM  = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 32, 6/27/16. NOT IMPLEMENTED, LACKS  CLARITY
!            L_HOM_R = sum( LCON(1:K_REAL(N),N) * L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) + &
!                      MCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N) * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) + &
!                      L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N),N) )  )
                L_HOM_R  = ZERO
                DO K = 1, K_REAL(N)
                  CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
                  CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                       L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                  T1 = LCON(K,N) * CPOS
                  T2 = MCON(K,N) * CNEG
                  L_HOM_R = L_HOM_R + T1 + T2
                ENDDO

!  Linearized Complex homogeneous solution contributions

                L_HOM_CR  = ZERO
                KO1 = K_REAL(N) + 1
                DO K = 1, K_COMPLEX(N)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                       L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
                  T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                       - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                       + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                       - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                  T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                       + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                       + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                       + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                  T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                  L_HOM_CR = L_HOM_CR + T1 + T2
                ENDDO

!  Final contribution

                COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

              ENDDO
            ENDDO
          ENDDO

!  otherwise if varying layer is above first active layer, there are bea
!  solution contributions propagated downwards - find these by calling
!  the appropriate solution module = L_BEAMSOLUTION_NNEK

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

          CALL LP_BEAMSOLUTION_NNEK ( &
               DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION,    & ! Input, Flags and order
               NSTOKES, NSTREAMS, N, LAYER_TO_VARY, N_LAYERWFS, FOURIER, IBEAM,  & ! Input, optical and control
               TAYLOR_ORDER, NSTREAMS_2, NSTKS_NSTRMS, K_REAL, DELTAU_VERT,      & ! Bookkeeping
               BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,         & ! Input, Beam Param.
               LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,             & ! Input, Beam Param (Linearized)
               T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,                & ! Input, PIs and linearized  
               GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,     & ! Input, Greens Function
               L_WUPPER, L_WLOWER )                                                ! Output

!  start loops

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
!  Enhancement # 33, 6/27/16
              COLTEL2_WF(CM,1:N_LAYERWFS) = - L_WUPPER(I,O1,N,1:N_LAYERWFS)
            ENDDO
          ENDDO

        ENDIF

!  Intermediate boundaries between active layers
!  ---------------------------------------------

        DO NS = 1, NLAYERS_TEL - 1

!  offsets

          N  = ACTIVE_LAYERS(NS)
          N1 = N + 1
          C0 = NS*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  if N is the varying layer, immediately above boundary

          IF ( N .EQ. LAYER_TO_VARY ) THEN

!  Get the linearized beam solution for the next layer N1

            CALL LP_BEAMSOLUTION_NNEK ( &
               DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION,    & ! Input, Flags and order
               NSTOKES, NSTREAMS, N1, LAYER_TO_VARY, N_LAYERWFS, FOURIER, IBEAM, & ! Input, optical and control
               TAYLOR_ORDER, NSTREAMS_2, NSTKS_NSTRMS, K_REAL, DELTAU_VERT,      & ! Bookkeeping
               BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,         & ! Input, Beam Param.
               LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,             & ! Input, Beam Param (Linearized)
               T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,                & ! Input, PIs and linearized  
               GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,     & ! Input, Greens Function
               L_WUPPER, L_WLOWER )                                                ! Output

!  start loops

            DO I = 1, NSTREAMS_2
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS

!  Beam contributions

                  L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 34, 6/27/16. NOT IMPLEMENTED, LACKS  CLARITY
!            L_HOM_R = sum( ( LCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) + &
!                   L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLA_XPOS(I,O1,1:K_REAL(N),N) ) ) + & 
!                        ( MCON(1:K_REAL(N),N) * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) ) )
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
                    CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                         L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                    T1 = LCON(K,N) * CPOS
                    T2 = MCON(K,N) * CNEG
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO

!  Linearized Complex homogeneous solution contributions

                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                         L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
                    T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                         - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                         + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                         - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                    T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                         + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                         + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                         + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO

!  Final contribution

                  COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

                ENDDO
              ENDDO
            ENDDO

!  If N1 is the varying layer, immediately below boundary
!    Only require contributions from this layer

          ELSE IF ( N1 .EQ. LAYER_TO_VARY ) THEN

!  start loops

            DO I = 1, NSTREAMS_2
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS

!  Beam contribution

                  L_BEAM  = + L_WUPPER(I,O1,N1,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 35, 6/27/16. NOT IMPLEMENTED, LACKS  CLARITY
!            L_HOM_R = sum( LCON(1:K_REAL(N1),N1) * L_SOLA_XPOS(I,O1,1:K_REAL(N1),N1,Q) + &
!                   MCON(1:K_REAL(N1),N1) *  ( T_DELT_EIGEN(1:K_REAL(N1),N1)   * L_SOLB_XNEG(I,O1,1:K_REAL(N1),N1,Q) + &
!                   L_T_DELT_EIGEN(1:K_REAL(N1),N1,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N1),N1) ) )
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N1)
                    CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
                    CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + &
                         L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
                    T1 = LCON(K,N1) * CPOS
                    T2 = MCON(K,N1) * CNEG
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO

!  Linearized Complex homogeneous solution contributions

                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N1) + 1
                  DO K = 1, K_COMPLEX(N1)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) - &
                         L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
                    T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                         - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                         + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1) &
                         - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
                    T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                         + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                         + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1) &
                         + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
                    T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO

!  Final contribution

                  COLTEL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  end loops

                ENDDO
              ENDDO
            ENDDO

!  non-zero variations if LAYER_TO_VARY is an active layer above N
!    Get the linearized beam solution for the next layer
!  .. contributions from beam solutions on both sides.

          ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

!  Get the linearized beam solution for the next layer

            CALL LP_BEAMSOLUTION_NNEK ( &
               DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION,    & ! Input, Flags and order
               NSTOKES, NSTREAMS, N1, LAYER_TO_VARY, N_LAYERWFS, FOURIER, IBEAM, & ! Input, optical and control
               TAYLOR_ORDER, NSTREAMS_2, NSTKS_NSTRMS, K_REAL, DELTAU_VERT,      & ! Bookkeeping
               BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,         & ! Input, Beam Param.
               LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,             & ! Input, Beam Param (Linearized)
               T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,                & ! Input, PIs and linearized  
               GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,     & ! Input, Greens Function
               L_WUPPER, L_WLOWER )                                                ! Output

!  .. contributions from beam solution (direct assign). No homog. variat

            DO I = 1, NSTREAMS_2
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
!  Enhancement # 36, 6/27/16
                COLTEL2_WF(CM,1:N_LAYERWFS) = L_WUPPER(I,O1,N1,1:N_LAYERWFS) - L_WLOWER(I,O1,N,1:N_LAYERWFS)
              ENDDO
            ENDDO

!  Finish different types of boundary condition linearizations

          ENDIF

!  End loop over intermediate active layer boundaries

        ENDDO

!  Final boundary, bottom of lowest active layer
!  ---------------------------------------------

        NS = NLAYERS_TEL
        N  = ACTIVE_LAYERS(NS)
        C0 = (NS-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  Old code. Condition is now completely general.
!  If this is the surface and Specialist option #2 is in place
!       if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 &
!             .AND.  N.EQ.NLAYERS ) THEN

!  get the linearized downward-reflected term. New 6/29/16 Rob Fix

        CALL LP_BVPTEL_SURFACE_SETUP &
          ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, IBEAM, FOURIER,         & ! Flags
            NSTOKES, NSTREAMS, N, NLAYERS, LAYER_TO_VARY, N_LAYERWFS,          & ! Numbers
            SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,           & ! Surface input
            MUELLER_INDEX, K_REAL, K_COMPLEX,                                  & ! bookkeeping
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input transmittances
            SOLA_XPOS, SOLB_XNEG, WLOWER, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,  & ! RT Solutions
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output reflected solutions

!   If varying layer is same as last active telescoped layer, need full solutions
!   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        IF ( LAYER_TO_VARY .EQ. N ) THEN

!  With-surface calculation
!  ^^^^^^^^^^^^^^^^^^^^^^^^

          IF ( DO_INCLUDE_SURFACE ) THEN

!  start loops

            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS

!  Beam contributions

                  L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions
!  Enhancement # 37, 6/27/16. NOT IMPLEMENTED, LACKS  CLARITY
!            L_HOM_R =sum( LCON(1:K_REAL(N),N) * ( ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLA_XPOS(I1,O1,1:K_REAL(N),N,Q) &
!                                      + L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLA_XPOS(I1,O1,1:K_REAL(N),N) ) &
!                                      - R2_L_HOMP(I,O1,1:K_REAL(N),Q) ) + &
!                          MCON(1:K_REAL(N),N) * ( L_SOLB_XNEG(I1,O1,1:K_REAL(N),N,Q) - R2_L_HOMM(I,O1,1:K_REAL(N),Q) ) )
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                        + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                    CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
                    CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                    CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
                    T1 = LCON(K,N) * CPOS
                    T2 = MCON(K,N) * CNEG
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                        - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                        + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                    T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                        + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                        + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                        + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                    T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
                    T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
                    T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO

!  Final contributions

                  COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

                 ENDDO
               ENDDO
             ENDDO

!  No-surface calculation
!  ^^^^^^^^^^^^^^^^^^^^^^

           ELSE

!  start loops

            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS

!  Beam contributions

                  L_BEAM = L_WLOWER(I1,O1,N,Q)

!  Linearized Real homogeneous solution contributions

                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                        + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                    CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                    T1 = LCON(K,N) * CPOS
                    T2 = MCON(K,N) * CNEG
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                        - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                        + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                    T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                        + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                        + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                        + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = L_SOLB_XNEG(I1,O1,K1,N,Q)
                    T2I = L_SOLB_XNEG(I1,O1,K2,N,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO

!  Final contributions

                  COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

                ENDDO
              ENDDO
            ENDDO

!  End surface vs. no-surface clause

          ENDIF 

!  If varying layer is below last active telescoped layer, need only the reflected solutions
!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        ELSE IF ( N .LT. LAYER_TO_VARY ) THEN

          IF ( DO_INCLUDE_SURFACE ) THEN
!  start loops
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS
!  Beam and Real-Hom
                  L_BEAM = R2_L_BEAM(I,O1,Q)
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    T1 = LCON(K,N) * R2_L_HOMP(I,O1,K,Q)
                    T2 = MCON(K,N) * R2_L_HOMM(I,O1,K,Q)
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO
!  Complex Hom
                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                    T1R = R2_L_HOMP(I,O1,K1,Q) ; T1I = R2_L_HOMP(I,O1,K2,Q)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = R2_L_HOMM(I,O1,K1,Q) ; T2I = R2_L_HOMM(I,O1,K2,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO
                  COLTEL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!   If varying layer is above last active layer, need Beam solutions propagated downwards, and reflected
!   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

          IF ( DO_INCLUDE_SURFACE ) THEN
!  start loops
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS
!  Beam and Real-Hom
                  L_BEAM = R2_L_BEAM(I,O1,Q)
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    T1 = LCON(K,N) * R2_L_HOMP(I,O1,K,Q)
                    T2 = MCON(K,N) * R2_L_HOMM(I,O1,K,Q)
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO
!  Complex Hom
                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                    T1R = R2_L_HOMP(I,O1,K1,Q) ; T1I = R2_L_HOMP(I,O1,K2,Q)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = R2_L_HOMM(I,O1,K1,Q) ; T2I = R2_L_HOMM(I,O1,K2,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO
                  COLTEL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
!  Enhancement # 38, 6/27/16
                COLTEL2_WF(CM,1:N_LAYERWFS) = - L_WLOWER(I1,O1,N,1:N_LAYERWFS)
              ENDDO
            ENDDO
          ENDIF 

!  End of variational settings by layer

        ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

!  Add direct beam solution. This is new code, R. Spurr 06/29/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.
!    --- Formerly: Only for Specialist option 2, Fourier = 0 Lambertian, lowest active layer at surface

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              BEAM = DIRECT_BEAM(I,IBEAM,O1) 
              FAC = - BEAM * DELTAU_SLANT(NLAYERS,LAYER_TO_VARY,IBEAM)
              DO Q = 1, N_LAYERWFS
                L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
                COLTEL2_WF(CM,Q) = COLTEL2_WF(CM,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  single-active-layer case
!  ========================

      ELSE

!  zero column vector. Extremely important.
!  Enhancement # 42, 6/27/16
        SCOL2_WF(1:NSTKS_NSTRMS_2, 1:MAX_ATMOSWFS) = ZERO

!  top of active layer
!  -------------------

!  If active layer = layer that is varying,
!       then require homogeneous and beam solution linearizations

        IF ( LAYER_TO_VARY .EQ. N ) THEN

!  Start  loops

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              DO Q = 1, N_LAYERWFS

!  beam solution linearization at top of layer

                L_BEAM = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 43, 6/27/16. NOT IMPLEMENTED, LACKS  CLARITY
!            L_HOM_R = sum( LCON(1:K_REAL(N),N) * L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) + &
!                           MCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) + &
!                                                 L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N),N) ) )
                L_HOM_R  = ZERO
                DO K = 1, K_REAL(N)
                  CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
                  CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                       L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                  T1 = LCON(K,N) * CPOS
                  T2 = MCON(K,N) * CNEG
                  L_HOM_R = L_HOM_R + T1 + T2
                ENDDO

!  Linearized Complex homogeneous solution contributions

                L_HOM_CR  = ZERO
                KO1 = K_REAL(N) + 1
                DO K = 1, K_COMPLEX(N)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                       L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
                  T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                       - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                       + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                       - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                  T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                       + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                       + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                       + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                  T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                  L_HOM_CR = L_HOM_CR + T1 + T2
                ENDDO

!  Final contribution

                SCOL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

              ENDDO
            ENDDO
          ENDDO

!  otherwise if varying layer is above active layer, there are beam
!  solution contributions propagated downwards - find these by calling
!  the appropriate solution module = L_BEAMSOLUTION_NNEK

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

          CALL LP_BEAMSOLUTION_NNEK ( &
             DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION,   & ! Input, Flags and order
             NSTOKES, NSTREAMS, N, LAYER_TO_VARY, N_LAYERWFS, FOURIER, IBEAM, & ! Input, optical and control
             TAYLOR_ORDER, NSTREAMS_2, NSTKS_NSTRMS, K_REAL, DELTAU_VERT,     & ! Bookkeeping
             BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,        & ! Input, Beam Param.
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,            & ! Input, Beam Param (Linearized)
             T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, BVEC, LP_BVEC,               & ! Input, PIs and linearized  
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,    & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                               ! Output

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
!  Enhancement # 44, 6/27/16
              SCOL2_WF(IROW,1:N_LAYERWFS) = - L_WUPPER(I,O1,N,1:N_LAYERWFS)
            ENDDO
          ENDDO

        ENDIF

!  Bottom of active layer
!  ----------------------

        C0 = NSTKS_NSTRMS

!  Old code. Condition is now completely general.
!  If this is the surface and Specialist option #2 is in place
!       if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 &
!             .AND.  N.EQ.NLAYERS ) THEN

!  get the linearized downward-reflected term. New 6/29/16 Rob Fix

        CALL LP_BVPTEL_SURFACE_SETUP &
          ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, IBEAM, FOURIER,         & ! Flags
            NSTOKES, NSTREAMS, N, NLAYERS, LAYER_TO_VARY, N_LAYERWFS,          & ! Numbers
            SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,           & ! Surface input
            MUELLER_INDEX, K_REAL, K_COMPLEX,                                  & ! bookkeeping
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input transmittances
            SOLA_XPOS, SOLB_XNEG, WLOWER, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,  & ! RT Solutions
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output reflected solutions

!   If varying layer is same as last active telescoped layer, need full solutions
!   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        IF ( LAYER_TO_VARY .EQ. N ) THEN

!  With-surface calculation
!  ^^^^^^^^^^^^^^^^^^^^^^^^

          IF ( DO_INCLUDE_SURFACE ) THEN

!  start loops

            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS

!  Beam contributions

                  L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions
!  Enhancement # 37, 6/27/16. NOT IMPLEMENTED, LACKS  CLARITY
!            L_HOM_R =sum( LCON(1:K_REAL(N),N) * ( ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLA_XPOS(I1,O1,1:K_REAL(N),N,Q) &
!                                      + L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLA_XPOS(I1,O1,1:K_REAL(N),N) ) &
!                                      - R2_L_HOMP(I,O1,1:K_REAL(N),Q) ) + &
!                          MCON(1:K_REAL(N),N) * ( L_SOLB_XNEG(I1,O1,1:K_REAL(N),N,Q) - R2_L_HOMM(I,O1,1:K_REAL(N),Q) ) )
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                        + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                    CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
                    CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                    CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
                    T1 = LCON(K,N) * CPOS
                    T2 = MCON(K,N) * CNEG
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                        - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                        + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                    T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                        + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                        + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                        + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                    T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
                    T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
                    T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO

!  Final contributions

                  SCOL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

                 ENDDO
               ENDDO
             ENDDO

!  No-surface calculation
!  ^^^^^^^^^^^^^^^^^^^^^^

           ELSE

!  start loops

            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS

!  Beam contributions

                  L_BEAM = L_WLOWER(I1,O1,N,Q)

!  Linearized Real homogeneous solution contributions

                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                        + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                    CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                    T1 = LCON(K,N) * CPOS
                    T2 = MCON(K,N) * CNEG
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2
                    K1 = KO1 + K0
                    K2 = K1  + 1
                    T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                        - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                        + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                    T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                        + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                        + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                        + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = L_SOLB_XNEG(I1,O1,K1,N,Q)
                    T2I = L_SOLB_XNEG(I1,O1,K2,N,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO

!  Final contributions

                  SCOL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

                ENDDO
              ENDDO
            ENDDO

!  End surface vs. no-surface clause

          ENDIF 

!  If varying layer is below last active telescoped layer, need only the reflected solutions
!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        ELSE IF ( N .LT. LAYER_TO_VARY ) THEN

          IF ( DO_INCLUDE_SURFACE ) THEN
!  start loops
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS
!  Beam and Real-Hom
                  L_BEAM = R2_L_BEAM(I,O1,Q)
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    T1 = LCON(K,N) * R2_L_HOMP(I,O1,K,Q)
                    T2 = MCON(K,N) * R2_L_HOMM(I,O1,K,Q)
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO
!  Complex Hom
                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                    T1R = R2_L_HOMP(I,O1,K1,Q) ; T1I = R2_L_HOMP(I,O1,K2,Q)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = R2_L_HOMM(I,O1,K1,Q) ; T2I = R2_L_HOMM(I,O1,K2,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO
                  SCOL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!   If varying layer is above last active layer, need Beam solutions propagated downwards, and reflected
!   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

          IF ( DO_INCLUDE_SURFACE ) THEN
!  start loops
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                DO Q = 1, N_LAYERWFS
!  Beam and Real-Hom
                  L_BEAM = R2_L_BEAM(I,O1,Q)
                  L_HOM_R  = ZERO
                  DO K = 1, K_REAL(N)
                    T1 = LCON(K,N) * R2_L_HOMP(I,O1,K,Q)
                    T2 = MCON(K,N) * R2_L_HOMM(I,O1,K,Q)
                    L_HOM_R = L_HOM_R + T1 + T2
                  ENDDO
!  Complex Hom
                  L_HOM_CR  = ZERO
                  KO1 = K_REAL(N) + 1
                  DO K = 1, K_COMPLEX(N)
                    K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                    T1R = R2_L_HOMP(I,O1,K1,Q) ; T1I = R2_L_HOMP(I,O1,K2,Q)
                    T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                    T2R = R2_L_HOMM(I,O1,K1,Q) ; T2I = R2_L_HOMM(I,O1,K2,Q)
                    T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                    L_HOM_CR = L_HOM_CR + T1 + T2
                  ENDDO
                  SCOL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
!  Enhancement # 38, 6/27/16
                SCOL2_WF(CM,1:N_LAYERWFS) = - L_WLOWER(I1,O1,N,1:N_LAYERWFS)
              ENDDO
            ENDDO
          ENDIF 

!  End of variational settings by layer

        ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

!  Add direct beam solution. This is new code, R. Spurr 06/29/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.
!    --- Formerly: Only for Specialist option 2, Fourier = 0 Lambertian, lowest active layer at surface

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              BEAM = DIRECT_BEAM(I,IBEAM,O1) 
              FAC = - BEAM * DELTAU_SLANT(NLAYERS,LAYER_TO_VARY,IBEAM)
              DO Q = 1, N_LAYERWFS
                L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
                SCOL2_WF(CM,Q) = SCOL2_WF(CM,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  end clause NLAYERS_TEL > 1 vs NLAYERS_TEL = 1
!  =============================================

      ENDIF

!  finish

      RETURN
      END SUBROUTINE LP_BVPTEL_COLUMN_SETUP

      END MODULE vlidort_lp_bvproblem_m

