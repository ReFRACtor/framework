
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
! ###########################################################

! ###############################################################
! #                                                             #
! #              FIRST-ORDER SCALAR/VECTOR MODEL                #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)            #
! #                                                             #
! #  This Version :   1.5.3                                     #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #   Version 1.1,   13 February  2012, First Code              #
! #   Version 1.2,   01 June      2012, Modularization          #
! #   Version 1.3a,  29 October   2012, Obsgeom Multi-geom.     #
! #   Version 1.3b,  24 January   2013, BRDF/SL Supplements     #
! #   Version 1.4,   31 July      2013, Lattice Multi-geom.     #
! #   Version 1.5,   7  July      2016. Use Fmatrix/Phasfunc    #
! #   Version 1.5,   22 August    2016. Partial-layer output.   #
! #   Version 1.5,   30 April     2017. Shakedown completed.    #
! #   Version 1.5.1, 30 September 2019. Revision Thermal DT.    #
! #   Version 1.5.3, 31 March     2021. Doublet geometry.       #
! #                                                             #
! #   FO Version 1.5   coincides (V)LIDORT Version (2.8)3.8     #
! #   FO Version 1.5.1 coincides (V)LIDORT Version (2.8.1)3.8.1 #
! #   FO Version 1.5.3 coincides (V)LIDORT Version (2.8.3)3.8.3 #
! #                                                             #
! ###############################################################

!    ###########################################################
!    #                                                         #
!    # This is Version 1.5.3 of the FO software library.       #
!    # This library comes with the GNU General Public License, #
!    # Version 3.0. Please read this license carefully.        #
!    #                                                         #
!    #      Copyright (c) 2010-2021.                           #
!    #          Robert Spurr, RT Solutions Inc.                #
!    #                                                         #
!    # This file is part of FO CODE Version 1.5.3.             #
!    #                                                         #
!    # FO CODE is free software: you can redistribute it       #
!    # and/or modify it under the terms of the GNU General     #
!    # Public License as published by the Free Software        #
!    # Foundation, either version 3 of the License, or any     #
!    # later version.                                          #
!    #                                                         #
!    # FO CODE is distributed in the hope that it will be      #
!    # useful, but WITHOUT ANY WARRANTY; without even the      #
!    # implied warranty of MERCHANTABILITY or FITNESS FOR A    #
!    # PARTICULAR PURPOSE.  See the GNU General Public License #
!    # for more details.                                       #
!    #                                                         #
!    # You should have received a copy of the GNU General      #
!    # Public License along with FO CODE Version 1.5.3         #
!    # If not, see <http://www.gnu.org/licenses/>.             #
!    #                                                         #
!    ###########################################################

module FO_VectorSS_RTCalcs_I_m

!  For given wavelength, routine calculates First-Order upwelling+downwelling Stokes vectors

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam solar sources,

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional F-matrices, Surface leaving, and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional F-matrix usage
!    Version 5,  02 August   2016. Surface leaving + Jacobians
!    Version 5,  20 August   2016, Partial-layer output

!  For Solar sources, the subroutines are
!       SSV_Integral_I_UP   (Upwelling only)
!       SSV_Integral_I_DN   (Downwelling only)
!       SSV_Integral_I_UPDN (Upwelling and Downwelling)

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SSV_Integral_I_UP, SSV_Integral_I_UPDN
!    ==> Add LOSTRANS_DN to the output lists from SSV_Integral_I_DN, SSV_Integral_I_UPDN

!  All subroutines public

public

contains

subroutine SSV_Integral_I_UP &
   ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,        & ! Inputs (dimension)
     do_sunlight, do_deltam_scaling, do_fmatrix, do_lambertian,                           & ! Inputs (Flags-General/Surface)
     do_surface_leaving, do_water_leaving,                                                & ! Inputs (Flags-General/Surface)
     do_Partials, do_PlanPar, do_enhanced_ps, flux, fvec, do_sources_up, do_sources_up_p, & ! Inputs(Flags/Flux/criticality)
     nstokes, ngeoms, nlayers, nfinedivs, nmomsinp, n_user_levels, user_levels,           & ! Inputs (control)
     npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
     extinction, deltaus, omega, truncfac, greekmat, fmatrix_up, Reflec, Slterm,          & ! Inputs (Optical/surface)
     Mu0, Mu1, GenSpher, Rotations, LosW_paths, LosP_paths,                               & ! Inputs (Geometry)
     xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                      & ! Inputs (Geometry)
     xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,          & ! Inputs (Geometry)
     stokes_up, stokes_db, cumsource_up, cumtrans, lostrans_up )                            ! Outputs

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Stokes-vector. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   - Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using F Matrices directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/20/16
!    - Partial-layer output introduced
  
!  Version 1.5 for VLIDORT 2.8.1
!    - Introduce water-leaving flag, output cumtrans

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SSV_Integral_I_UP

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   INTEGER, Intent(in) :: maxgeoms
   INTEGER, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials

   INTEGER, Intent(in) :: maxfine
   INTEGER, Intent(in) :: maxmoments_input
   INTEGER, Intent(in) :: max_user_levels

!  flags.
!  Version 1.5: --> F-matrix flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/20/16

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING
   LOGICAL, Intent(in) :: DO_FMATRIX

   LOGICAL, Intent(in) :: DO_LAMBERTIAN
   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   logical, Intent(in) :: DO_WATER_LEAVING

   logical, Intent(in) :: DO_Partials
   LOGICAL, Intent(in) :: DO_PLANPAR
   LOGICAL, Intent(in) :: DO_ENHANCED_PS

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX, FVEC(4)

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)

!  Numbers

   INTEGER, Intent(in) :: NSTOKES
   INTEGER, Intent(in) :: NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   INTEGER, Intent(in) :: NMOMSINP
   INTEGER, Intent(in) :: N_USER_LEVELS
   INTEGER, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/20/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere. Fmatrix input added 7/7/16

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )
   REAL(ffp), Intent(in) :: FMATRIX_UP  ( MAXLAYERS, MAXGEOMS, 6 )

!  Surface reflectivity (Could be the albedo)
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( 4, 4, MAXGEOMS )
   real(ffp), Intent(in) :: SLTERM ( 4,    MAXGEOMS )

!  Geometrical inputs
!  ------------------

!  Ray constants, Mu0, Mu1
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only
!   real(ffp), Intent(in)  :: Raycon(maxgeoms)

   real(ffp), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

   REAL(ffp), Intent(in)  :: GenSpher(0:maxmoments_input,4,maxgeoms)
   REAL(ffp), Intent(in)  :: Rotations(4,maxgeoms)

!  Los paths added, 8/20/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/20/16.

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/20/16.

   integer  , Intent(in)  :: ntraverse     (0:maxlayers,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

   REAL(ffp), Intent(Out)  :: stokes_up     ( max_user_levels, 4, maxgeoms )
   REAL(ffp), Intent(Out)  :: stokes_db     ( max_user_levels, 4, maxgeoms )
   REAL(ffp), Intent(Out)  :: cumsource_up  ( 0:maxlayers,     4, maxgeoms )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS ( max_user_levels, maxgeoms )

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SSV_Integral_I_UP

   real(ffp), Intent(out)  :: lostrans_up (maxgeoms,maxlayers)

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/20/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: attenuationsfine   (maxlayers,maxfine)
   real(ffp)  :: attenuations_p     (maxpartials)
   real(ffp)  :: attenuationsfine_p (maxpartials,maxfine)

!  Scattering

   REAL(ffp)  :: tms (maxlayers)
   REAL(ffp)  :: exactscat_up (maxlayers,4,4)

!  Source function integration results
!    -- 1/31/21. Version 2.8.3. removed LOSTRANS_UP from this list, now an output

   real(ffp)  :: sources_up  (maxlayers, 4), sources_up_p  (maxpartials, 4)
   real(ffp)  :: lostrans_up_p (maxpartials)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: multiplier_p     ( maxpartials )

!  Regular_PS or plane-parallel flag

   LOGICAL    :: do_RegPSorPP

!  Help

   LOGICAL    :: do_regular_ps, layermask_up(maxlayers)
   INTEGER    :: n, ns, uta, nstart, nc, nut, nut_prev, j, L, o1, v, nt, np, ut

   REAL(ffp)  :: sumd, help, sum, tran, factor1, factor2, kn, pi4, dj, path_up, Shelp(4)
   REAL(ffp)  :: cumsource_db(4), fmat(6), sum23, dif23, ctrans, CUMSOURCE_DB_START(4)
   REAL(ffp)  :: help3c1, help3s1, help4c1, help4s1, m4
   real(ffp)  :: suntau(0:maxlayers), suntau_p(maxpartials), lostau

   REAL(ffp), parameter  :: cutoff = 88.0_ffp
   REAL(ffp), parameter  :: zero   = 0.0_ffp
   REAL(ffp), parameter  :: one    = 1.0_ffp

!  Number

!mick fix 9/19/2017 - define pi4 as in VLIDORT_PARS
   !pi4 = acos(-one)/4.0_ffp
   pi4 = acos(-one)*4.0_ffp

!  Zero the output. 4/9/19 include CUMTRANS

   CUMSOURCE_UP  = zero ; STOKES_UP  = zero ; STOKES_DB    = zero ; cumtrans = zero

!  1/31/21. Version 2.8.3. lostrans_up, zeroed here because it has an extra maxgeoms dimension

   lostrans_up = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes ; fmat = zero
   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

!  TMS factors

   do n = 1, nlayers
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
      else
         tms(n) = omega(n)
      endif
   enddo

!  Start geometry loop
!  -------------------

   do v = 1, ngeoms

!  Zero the local sources
!    -- 1/31/21. Version 2.8.3. removed zeroing of lostrans_up

      sources_up   = zero  ; sources_up_p = zero 
      lostrans_up_p = zero ; exactscat_up = zero

!  Scattering functions
!  ====================

!  Version 1.5, F-matrix option introduced. 7/7/16
!mick note 9/19/2017 - defining F-matrix using either (1) F-matrix input ("do_fmatrix")
!                                                  or (2) Greek matrix moment input

!  Scalar only
!  -----------

      if ( nstokes .eq. 1 ) then
        do n = 1, nlayers
          if ( layermask_up(n) ) then

!  Get F-matrix

            if ( do_fmatrix ) then
              fmat(1) = fmatrix_up(n,v,1)
            else
              fmat(1) = zero
              do L = 0, nmomsinp
                fmat(1) = fmat(1) + GenSpher(L,1,V) * Greekmat(n,L,1,1)
              enddo
            endif
            exactscat_up(n,1,1) = fmat(1) * tms(n)
          endif
        enddo
      endif

!  Vector with Sunlight
!  --------------------

      if ( nstokes .gt. 1 .and. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_up(n) ) then

!  Get F-matrix

            if ( do_fmatrix ) then
              fmat(1:2) = fmatrix_up(n,v,1:2)
            else
              fmat(1:2) = zero
              do L = 0, nmomsinp
                fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
              enddo
            endif

!   Apply rotations for Z-matrix, multiply by TMS

            exactscat_up(n,1,1) = + fmat(1)
            exactscat_up(n,2,1) = - fmat(2) * Rotations(3,v)
            exactscat_up(n,3,1) = + fmat(2) * Rotations(4,v)
            exactscat_up(n,1:ns,1) = tms(n) * exactscat_up(n,1:ns,1) 

          endif
        enddo
      endif

!  Vector General case
!  -------------------

! USE FULL 4X4 MATRIX; CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_up(n) ) then

!  Get the F-matrix 

            if ( do_fmatrix ) then
              fmat(1:6) = fmatrix_up(n,v,1:6)
            else
              fmat = zero
              do L = 0, nmomsinp
                fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
                sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                fmat(3) = fmat(3) + Genspher(L,3,v) * sum23
                fmat(4) = fmat(4) + Genspher(L,4,v) * dif23
              enddo
              fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_ffp
              fmat(4) = ( fmat(3) - fmat(4) )
              if ( nstokes.eq.4) then
                do L = 0, nmomsinp
                  fmat(5) = fmat(5) + Genspher(L,2,v) * greekmat(n,L,3,4)
                  fmat(6) = fmat(6) + Genspher(L,1,v) * greekmat(n,L,4,4)
                enddo
              endif
            endif

!  Apply rotations to get Zmatrix

            help3c1 = fmat(3) * Rotations(1,v) ; help3s1 = fmat(3) * Rotations(2,v)
            help4c1 = fmat(4) * Rotations(1,v) ; help4s1 = fmat(4) * Rotations(2,v)
            exactscat_up(n,1,1) = + fmat(1)
            exactscat_up(n,2,1) = - fmat(2) * Rotations(3,v)
            exactscat_up(n,3,1) = + fmat(2) * Rotations(4,v)
            exactscat_up(n,1,2) = + fmat(2) * Rotations(1,v)
            exactscat_up(n,1,3) = - fmat(2) * Rotations(2,v)
            exactscat_up(n,2,2) = + help3c1 * Rotations(3,v) - help4s1 * Rotations(4,v)
            exactscat_up(n,2,3) = - help3s1 * Rotations(3,v) - help4c1 * Rotations(4,v)
            exactscat_up(n,3,2) = + help3c1 * Rotations(4,v) + help4s1 * Rotations(3,v)
            exactscat_up(n,3,3) = - help3s1 * Rotations(4,v) + help4c1 * Rotations(3,v)
            if ( nstokes .eq. 4 ) then
               exactscat_up(n,2,4) = - fmat(5) * Rotations(4,v) 
               exactscat_up(n,4,2) = - fmat(5) * Rotations(2,v) 
               exactscat_up(n,3,4) = + fmat(5) * Rotations(3,v) 
               exactscat_up(n,4,3) = - fmat(5) * Rotations(1,v) 
               exactscat_up(n,4,4) = + fmat(6)
            endif

!  Multiply by TMS

            exactscat_up(n,1:ns,1:ns) = tms(n)*exactscat_up(n,1:ns,1:ns)

!  End calculation

          endif
        enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations   = zero ; suntau = zero   ; Attenuationsfine   = zero
      Attenuations_p = zero ; suntau_p = zero ; Attenuationsfine_p = zero
      nstart = nlayers

!  Critical removed TRC
!  Initialize, only to layer Ncrit if applicable
!     if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = ntraverse(n,v) ; sumd = dot_product(extinction(1:nt),sunpaths(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  RobFix 8/20/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = ntraverse_p(ut,v) ; sumd = dot_product(extinction(1:nt),sunpaths_p(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. cutoff ) Attenuations_p(ut) = exp( - sumd )
        enddo
      endif

!  Enhanced-spherical, fine-layer attenuations, Whole-layer integration

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
            do j = 1, nfinedivs(n,v)
              nt = ntraversefine(n,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine(n,1:nt,j,v))
              if (sumd .lt. cutoff ) Attenuationsfine(n,j) = exp( - sumd )
            enddo
          endif
        enddo
      endif

!  RobFix 8/20/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, nfinedivs_p(ut,v)
              nt = ntraversefine_p(ut,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine_p(ut,1:nt,j,v))
              If (sumd .lt. cutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
            enddo
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

!  1/31/21. Version 2.8.3. lostrans_up, has to be given geometry index, as it is now an output

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_up(n) .and. do_sources_up(n,v)  ) then
            if ( Mu1(v) .gt. zero ) then
              lostau = deltaus(n)  / Mu1(v)
              if ( lostau .lt. cutoff ) lostrans_up(v,n) = exp( - lostau )
              factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(v,n)
              factor2 = one + (suntau(n) - suntau(n-1))/lostau
              multiplier(n) = factor1 / factor2
            endif
          endif
        enddo
      endif

!  RobFix 8/20/16. Plane/Parallel or Regular-PS, Partial-layer output

      if ( do_RegPSorPP .and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v) - LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero
            if ( Mu1(v) .gt. zero ) then
              lostau = kn * path_up
              if ( lostau .lt. cutoff ) lostrans_up_p(ut) = exp( - lostau )
              factor1 = Attenuations_p(ut) - Attenuations(np)*lostrans_up_p(ut)
              factor2 = one + (suntau(np) - suntau_p(ut))/lostau
              multiplier_p(ut) = factor1 / factor2
            endif
          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/20/16 streamlined code using distances
!      Quadratures from Bottom of the layer

!  1/31/21. Version 2.8.3. lostrans_up, has to be given geometry index, as it is now an output
!mick fix 1/13/2021 - swapped indexes in "do_sources_up"

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_up(n) .and. do_sources_up(n,v)   ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_up = LosW_paths(n,v)
            lostau = kn * path_up ; if( lostau.lt.cutoff ) lostrans_up(v,n) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * wfine(n,j,v)
            enddo
            multiplier(n) = sum * kn
          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/02/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v)- LosP_paths(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.cutoff ) lostrans_up_p(ut) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs_p(ut,v)
              dj = path_up - xfine_p(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine_p(ut,j) * tran * wfine_p(ut,j,v)
            enddo
            multiplier_p(ut) = sum * kn
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  General case, Whole layers
!mick fix 9/19/2017 - defined "sources_up" for all stokes vector elements

      do n = nlayers, 1, -1
        if ( layermask_up(n) .and. do_sources_up(n,v)  ) then
          if ( do_sunlight ) then
             shelp(1:ns) = exactscat_up(n,1:ns,1) * fvec(1)
          else
            do o1 = 1, nstokes
              shelp(o1) = dot_product(exactscat_up(n,o1,1:ns),fvec(1:ns))
            enddo
          endif
          !sources_up(n,o1) = shelp(o1) * multiplier(n)
          sources_up(n,1:ns) = shelp(1:ns) * multiplier(n)
        endif
      enddo

!  Partials case

      if ( do_partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v)  ) then
            np = partial_layeridx(ut)
            if ( do_sunlight ) then
              shelp(1:ns) = exactscat_up(np,1:ns,1) * fvec(1)
            else
              do o1 = 1, nstokes
                shelp(o1) = dot_product(exactscat_up(np,o1,1:ns),fvec(1:ns))
              enddo
            endif
            sources_up_p(ut,1:ns) = shelp(1:ns)* multiplier_p(ut)
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  initialize recursion ( For Direct Beam, use PI.mu0.R.Atten )
!  4/9/19. Add start of CUMTRANS recursion (CTRANS = 1.0). Rename CUMSOURCE_DB_START

      NC =  0
      CUMSOURCE_UP(NC,:,V) = zero
      CUMSOURCE_DB_START = zero
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Surface term

      M4 = 4.0d0 * MU0(V)
      if ( DO_LAMBERTIAN ) then
         CUMSOURCE_DB_START(1) = M4 * REFLEC(1,1,V) * attenuations(nlayers) * fvec(1)
      else
         do o1 = 1, nstokes
            sum = dot_product(REFLEC(O1,1:nstokes,V),fvec(1:nstokes))
            CUMSOURCE_DB_START(o1) = M4 * sum * attenuations(nlayers)      ! Bug 3/23/15 attenuations was left out
         enddo
      endif

!  surface-leaving term. Added, 8/2/16
!   -- (modeled after the DBCORRECTION code in Version 2.7)
!   -- 4/9/19. Not done for water-leaving, as need to use adjusted values

     IF ( DO_SURFACE_LEAVING .and. .not. DO_WATER_LEAVING ) THEN
        do o1 = 1, nstokes
           CUMSOURCE_DB_START(o1) = CUMSOURCE_DB_START(o1) + PI4 * SLTERM(o1,v)
        enddo
     ENDIF

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

!  1/31/21. Version 2.8.3. LOSTRANS_UP has to be given geometry index, as it is now an output

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            do o1 = 1, nstokes
               CUMSOURCE_DB_START(O1) = LOSTRANS_UP(V,N) * CUMSOURCE_DB_START(O1)
               CUMSOURCE_UP(NC,O1,V)  = LOSTRANS_UP(V,N) * CUMSOURCE_UP(NC-1,O1,V) + SOURCES_UP(N,O1)
            enddo
         ENDDO
         CUMSOURCE_DB(1:nstokes)     = CUMSOURCE_DB_START(1:nstokes)
         IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            do o1 = 1, nstokes
              STOKES_UP(UTA,O1,V) = FLUX * ( CUMSOURCE_UP(NC,O1,V) * LOSTRANS_UP_p(UT) + SOURCES_UP_p(UT,O1) )
              STOKES_DB(UTA,O1,V) = FLUX * CUMSOURCE_DB(O1) * LOSTRANS_UP_p(UT)
            enddo
         ELSE
           do o1 = 1, nstokes
             STOKES_UP(UTA,O1,V) = FLUX * CUMSOURCE_UP(NC,O1,V)
             STOKES_DB(UTA,O1,V) = FLUX * CUMSOURCE_DB(O1)
           enddo
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  4/9/19.  WATERLEAVING CASE. Add CUMTRANS calculation
!    Add start of CUMTRANS recursion (CTRANS = 1.0).

!  5/5/20. Version 2.8.1 Upgrades. CUMTRANS calculation was not properly initialized
!  1/31/21. Version 2.8.3. LOSTRANS_UP has to be given geometry index, as it is now an output

      if ( do_water_leaving ) then
         NSTART = NLAYERS                ! 5/5/20 Needed for initializing
         NUT_PREV = NSTART + 1           ! 5/5/20 Needed for initializing
         ctrans = one
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               CTRANS = CTRANS * LOSTRANS_UP(V,N)
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
               UT = Partial_OUTINDEX(UTA)
               CUMTRANS(UTA,V)     = CTRANS * LOSTRANS_UP_p(UT)
            ELSE
               CUMTRANS(UTA,V)     = CTRANS
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif
      
!  Finish geometry loop

   enddo

!  Finish

   return
end subroutine SSV_Integral_I_UP

!

subroutine SSV_Integral_I_DN &
   ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,  & ! Inputs (dimension)
     do_sunlight, do_deltam_scaling, do_fmatrix, do_Partials,                       & ! Inputs (Flags)
     do_PlanPar, do_enhanced_ps, flux, fvec, do_sources_dn, do_sources_dn_p,        & ! Inputs (Flags/flux)
     nstokes, ngeoms, nlayers, nfinedivs, nmomsinp, n_user_levels, user_levels,     & ! Inputs (control)
     npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,   & ! Inputs (control-partial)
     extinction, deltaus, omega, truncfac, greekmat, fmatrix_dn,                    & ! Inputs (Optical)
     Mu1, GenSpher, Rotations, LosW_paths, LosP_paths,                              & ! Inputs (Geometry)
     xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                & ! Inputs (Geometry)
     xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,    & ! Inputs (Geometry)
     stokes_dn, cumsource_dn, lostrans_dn )                                           ! Outputs

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of Stokes vector. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using F Matrices directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/20/16
!    - Partial-layer output introduced

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_DN to the output lists from SSV_Integral_I_DN

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   INTEGER, Intent(in) :: maxgeoms
   INTEGER, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials
   INTEGER, Intent(in) :: maxfine
   INTEGER, Intent(in) :: maxmoments_input
   INTEGER, Intent(in) :: max_user_levels

!  flags. T-matrix flag added, 7/7/16

   LOGICAL, Intent(in) ::  DO_SUNLIGHT
   LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
   LOGICAL, Intent(in) ::  DO_FMATRIX

   logical, Intent(in) ::  DO_Partials
   LOGICAL, Intent(in) ::  DO_PLANPAR
   LOGICAL, Intent(in) ::  DO_ENHANCED_PS

!  Solar Flux 

   REAL(ffp), Intent(in) :: FLUX, FVEC(4)

!  Existence flags. 8/20/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Numbers

   INTEGER, Intent(in) ::  NSTOKES
   INTEGER, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   INTEGER, Intent(in) ::  NMOMSINP
   INTEGER, Intent(in) ::  N_USER_LEVELS
   INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/20/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere. Fmatrix input added 7/7/16

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )
   REAL(ffp), Intent(in) :: FMATRIX_DN  ( MAXLAYERS, MAXGEOMS, 6 )

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)
!   real(ffp), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
!   real(ffp), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)

!    Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp), Intent(in)  :: Mu1(maxgeoms)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

   REAL(ffp), Intent(in)  :: GenSpher(0:maxmoments_input,4,maxgeoms)
   REAL(ffp), Intent(in)  :: Rotations(4,maxgeoms)

!  Los paths added, 8/20/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/20/16.

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/20/16.

   integer  , Intent(in)  :: ntraverse     (0:maxlayers,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

   REAL(ffp), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, maxgeoms )
   REAL(ffp), Intent(Out)  :: cumsource_dn  ( 0:maxlayers,     4, maxgeoms )

!  1/31/21. Version 2.8.3.  Add LOSTRANS_DN to the output list

   real(ffp), Intent(out)  :: lostrans_dn (maxgeoms,maxlayers)

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/20/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: attenuationsfine (maxlayers,maxfine)
   real(ffp)  :: attenuations_p     (maxpartials)
   real(ffp)  :: attenuationsfine_p (maxpartials,maxfine)

!  Scattering

   REAL(ffp)  :: tms (maxlayers)
   REAL(ffp)  :: exactscat_dn (maxlayers,4,4)

!  Source function integration results
!  1/31/21. Version 2.8.3. Removed LOSTRANS_DN from this list, now an output

   real(ffp)  :: sources_dn  (maxlayers, 4), sources_dn_p  (maxpartials, 4)
   real(ffp)  :: lostrans_dn_p (maxpartials)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: multiplier_p     ( maxpartials )

!  Regular_PS or plane-parallel flag

   LOGICAL    :: do_RegPSorPP

!  Help

   INTEGER    :: n, ns, uta, nstart, nc, nut, nut_prev, j, L, O1, v, nt, ut, np
   logical    :: do_regular_ps, layermask_dn(maxlayers)

   REAL(ffp)  :: sumd, help, sum, tran, factor1, factor2, kn, dj
   real(ffp)  :: suntau(0:maxlayers), suntau_p(maxpartials), lostau, path_dn, Shelp(4)
   REAL(ffp)  :: fmat(6), sum23, dif23, help3c1, help3s1, help4c1, help4s1

   REAL(ffp), parameter  :: cutoff = 88.0_ffp
   REAL(ffp), parameter  :: zero   = 0.0_ffp
   REAL(ffp), parameter  :: one    = 1.0_ffp

!  Zero the output

   CUMSOURCE_DN = zero ; STOKES_DN  = zero

!  1/31/21. Version 2.8.3. lostrans_dn, zeroed here because it has an extra maxgeoms dimension

   lostrans_dn = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes ; fmat = zero
   NUT = USER_LEVELS(N_USER_LEVELS) + 1
   IF ( NUT > NLAYERS ) NUT = NLAYERS
   LAYERMASK_DN = .false.
   LAYERMASK_DN(1:NUT) = .true.

!  TMS factors

   do n = 1, nlayers
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
      else
         tms(n) = omega(n)
      endif
   enddo

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources
!    -- 1/31/21. Version 2.8.3. removed zeroing of lostrans_dn

      sources_dn    = zero ; exactscat_dn = zero
      lostrans_dn_p = zero ; sources_dn_p = zero

!  Scattering functions
!  ====================

!  Version 1.5, F-matrix option introduced. 7/7/16

!  Scalar only
!  -----------

      if ( nstokes .eq. 1 ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            if ( do_fmatrix ) then
              fmat(1) = fmatrix_dn(n,v,1)
            else
              fmat(1) = zero
              do L = 0, nmomsinp
                fmat(1) = fmat(1) + GenSpher(L,1,v) * Greekmat(n,L,1,1)
              enddo
            endif
            exactscat_dn(n,1,1) = fmat(1) * tms(n)
          endif
        enddo
      endif

!  Vector with Sunlight
!  --------------------

      if ( nstokes .gt. 1 .and. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then

!  get F-matrix

            if ( do_fmatrix ) then
              fmat(1:2) = fmatrix_dn(n,v,1:2)
            else
              fmat(1:2) = zero
              do L = 0, nmomsinp
                fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
              enddo
            endif

!   Apply rotations for Z-matrix, multiply by TMS

            exactscat_dn(n,1,1) = + fmat(1)
            exactscat_dn(n,2,1) = - fmat(2) * Rotations(3,v)
            exactscat_dn(n,3,1) = + fmat(2) * Rotations(4,v)
            exactscat_dn(n,1:ns,1) = tms(n) * exactscat_dn(n,1:ns,1) 

          endif
        enddo
      endif

!  Vector General case
!  -------------------

! USE FULL 4X4 MATRIX; CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then

!  Get the F-matrix 

            if ( do_fmatrix ) then
              fmat(1:6) = fmatrix_dn(n,v,1:6)
            else
              fmat = zero
              do L = 0, nmomsinp
                fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
                sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                fmat(3) = fmat(3) + Genspher(L,3,v) * sum23
                fmat(4) = fmat(4) + Genspher(L,4,v) * dif23
              enddo
              fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_ffp
              fmat(4) = ( fmat(3) - fmat(4) )
              if ( nstokes.eq.4) then
                do L = 0, nmomsinp
                  fmat(5) = fmat(5) + Genspher(L,2,v) * greekmat(n,L,3,4)
                  fmat(6) = fmat(6) + Genspher(L,1,v) * greekmat(n,L,4,4)
                enddo
              endif
            endif

!  Apply the rotations to get Z-matrix

            help3c1 = fmat(3) * Rotations(1,v)
            help3s1 = fmat(3) * Rotations(2,v)
            help4c1 = fmat(4) * Rotations(1,v)
            help4s1 = fmat(4) * Rotations(2,v)
            exactscat_dn(n,1,1) = + fmat(1)
            exactscat_dn(n,2,1) = - fmat(2) * Rotations(3,v)
            exactscat_dn(n,3,1) = + fmat(2) * Rotations(4,v)
            exactscat_dn(n,1,2) = + fmat(2) * Rotations(1,v)
            exactscat_dn(n,1,3) = - fmat(2) * Rotations(2,v)
            exactscat_dn(n,2,2) = + help3c1 * Rotations(3,v) - help4s1 * Rotations(4,v)
            exactscat_dn(n,2,3) = - help3s1 * Rotations(3,v) - help4c1 * Rotations(4,v)
            exactscat_dn(n,3,2) = + help3c1 * Rotations(4,v) + help4s1 * Rotations(3,v)
            exactscat_dn(n,3,3) = - help3s1 * Rotations(4,v) + help4c1 * Rotations(3,v)
            if ( nstokes .eq. 4 ) then
              exactscat_dn(n,2,4) = - fmat(5) * Rotations(4,v) 
              exactscat_dn(n,4,2) = - fmat(5) * Rotations(2,v) 
              exactscat_dn(n,3,4) = + fmat(5) * Rotations(3,v) 
              exactscat_dn(n,4,3) = - fmat(5) * Rotations(1,v) 
              exactscat_dn(n,4,4) = + fmat(6)
            endif

!  Multiply by TMS

            exactscat_dn(n,1:ns,1:ns) = tms(n)*exactscat_dn(n,1:ns,1:ns)

          endif
        enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize

      Attenuations   = zero ; suntau = zero   ; Attenuationsfine   = zero
      Attenuations_p = zero ; suntau_p = zero ; Attenuationsfine_p = zero
      nstart = nlayers

!  Critical removed TRC
!  Initialize, only to layer Ncrit if applicable
!     if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = ntraverse(n,v) ; sumd = dot_product(extinction(1:nt),sunpaths(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  RobFix 8/20/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = ntraverse_p(ut,v) ; sumd = dot_product(extinction(1:nt),sunpaths_p(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. cutoff ) Attenuations_p(ut) = exp( - sumd )
        enddo
      endif

!  Adjust nstart

      do n = 1, nlayers
         if ( layermask_dn(n) .and. attenuations(n-1).ne.zero )  nstart = n
      enddo

!  Enhanced-spherical, fine-layer attenuations, Whole-layer integration

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
            do j = 1, nfinedivs(n,v)
              nt = ntraversefine(n,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine(n,1:nt,j,v))
              if (sumd .lt. cutoff ) Attenuationsfine(n,j) = exp( - sumd )
            enddo
          endif
        enddo
      endif

!  RobFix 8/20/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, nfinedivs_p(ut,v)
              nt = ntraversefine_p(ut,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine_p(ut,1:nt,j,v))
              If (sumd .lt. cutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
            enddo
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

!  1/31/21. Version 2.8.3. LOSTRANS_DN has to be given geometry index, as it is now an output

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
            if ( Mu1(v) .gt. zero ) then
              lostau = deltaus(n)  / Mu1(v)
              if ( lostau .lt. cutoff ) lostrans_dn(v,n) = exp( - lostau )
              factor1 = Attenuations(n-1)*lostrans_dn(v,n) - Attenuations(n)
              factor2 = ((suntau(n) - suntau(n-1))/lostau) - one
              multiplier(n) = factor1 / factor2
            endif
          endif
        enddo
      endif


!  RobFix 8/20/16. Plane/Parallel or Regular-PS, Partial-layer output

      if ( do_RegPSorPP .and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero
            if ( Mu1(v) .gt. zero ) then
              lostau = kn * path_dn
              if ( lostau .lt. cutoff ) lostrans_dn_p(ut) = exp( - lostau )
              factor1 = Attenuations(np-1)*lostrans_dn_p(ut) - Attenuations_p(ut)
              factor2 = ((suntau_p(ut) - suntau(np-1))/lostau) - one
              multiplier_p(ut) = factor1 / factor2
            endif
          endif
        enddo
      endif

!  Enhanced PS: General case. RobFix 8/20/16 streamlined code using distances
!   Quadratures measured from the bottom

!  1/31/21. Version 2.8.3. LOSTRANS_DN has to be given geometry index, as it is now an output

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_dn = LosW_paths(n,v)
            lostau = kn * path_dn ; if( lostau.lt.cutoff ) lostrans_dn(v,n) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * wfine(n,j,v)
            enddo 
            multiplier(n) = sum * kn
          endif
        enddo
      endif

!  Enhanced PS:  RobFix 8/20/16 Partials

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.cutoff ) lostrans_dn_p(ut) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs_p(ut,v)
              dj = path_dn - xfine_p(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine_p(ut,j) * tran * wfine_p(ut,j,v)
            enddo
            multiplier_p(ut) = sum * kn
          endif
        enddo
      endif

!  Layer sources
!  -------------

!  General case, Whole layers
!mick fix 9/19/2017 - defined "sources_dn" for all stokes vector elements

      do n = 1, nlayers, 1
        if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
          if ( do_sunlight ) then
             shelp(1:ns) = exactscat_dn(n,1:ns,1) * fvec(1)
          else
            do o1 = 1, nstokes
              shelp(o1) = dot_product(exactscat_dn(n,o1,1:ns),fvec(1:ns))
            enddo
          endif
          !sources_dn(n,o1) = shelp(o1) * multiplier(n)
          sources_dn(n,1:ns) = shelp(1:ns) * multiplier(n)
        endif
      enddo

!  Partials case

      if ( do_partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v)  ) then
            np = partial_layeridx(ut)
            if ( do_sunlight ) then
              shelp(1:ns) = exactscat_dn(np,1:ns,1) * fvec(1)
            else
              do o1 = 1, nstokes
                shelp(o1) = dot_product(exactscat_dn(np,o1,1:ns),fvec(1:ns))
              enddo
            endif
            sources_dn_p(ut,1:ns) = shelp(1:ns)* multiplier_p(ut)
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,:,v) = zero
      NSTART = 1 ; NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working down from NSTART to NUT
!     Check for Updating the recursion

!  1/31/21. Version 2.8.3. LOSTRANS_DN has to be given geometry index, as it is now an output

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            do o1 = 1, nstokes
               CUMSOURCE_DN(NC,O1,v)  = LOSTRANS_DN(V,N) * CUMSOURCE_DN(NC-1,O1,v) + SOURCES_DN(N,O1)
            enddo
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           do o1 = 1, nstokes
             STOKES_DN(UTA,O1,v) = FLUX * ( CUMSOURCE_DN(NC,O1,V) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT,O1) )
           enddo
         ELSE
           do o1 = 1, nstokes
             STOKES_DN(UTA,O1,v) = FLUX * CUMSOURCE_DN(NC,O1,v)
           enddo
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Finish geometry loop

   enddo

!  Finish

   return
end subroutine SSV_Integral_I_DN

!

subroutine SSV_Integral_I_UPDN   &
   ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,                  & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_sunlight, do_deltam_scaling, do_fmatrix, do_Partials,           & ! Inputs (Flags)
     do_lambertian, do_surface_leaving, do_water_leaving, do_PlanPar, do_enhanced_ps,               & ! Inputs (Flags)
     do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,                                & ! Inputs (Flags - Geometry)
     nstokes, ngeoms, nlayers, nfinedivs, nmomsinp, n_user_levels, user_levels,                     & ! Inputs (control output)
     npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,                   & ! Inputs (control-partial)
     flux, fvec, extinction, deltaus, omega, truncfac, Greekmat, fmatrix_up, fmatrix_dn,            & ! Inputs (Flux/Optical)
     reflec, slterm, Mu0, Mu1, GenSpher_up, GenSpher_dn, Rotations_up, Rotations_dn, LosW_paths,    & ! Inputs (Optical)
     LosP_paths, xfine_up, wfine_up, sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,  & ! Inputs (Geometry)
     xfine_dn, wfine_dn, sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,              & ! Inputs (Geometry)
     xfine_up_p, wfine_up_p, sunpaths_up_p, ntraverse_up_p, sunpathsfine_up_p, ntraversefine_up_p,  & ! Inputs (Geometry)
     xfine_dn_p, wfine_dn_p, sunpaths_dn_p, ntraverse_dn_p, sunpathsfine_dn_p, ntraversefine_dn_p,  & ! Inputs (Geometry)
     stokes_up, stokes_db, cumsource_up, cumtrans, lostrans_up, stokes_dn, cumsource_dn, lostrans_dn )  ! Outputs

!  Stand-alone routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of Stokes-vector. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using F Matrices directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/20/16
!    - Partial-layer output introduced

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SSV_Integral_I_UP, SSV_Integral_I_UPDN
!    ==> Add LOSTRANS_DN to the output lists from SSV_Integral_I_DN, SSV_Integral_I_UPDN

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   INTEGER, Intent(in) :: maxgeoms
   INTEGER, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials
   INTEGER, Intent(in) :: maxfine

   INTEGER, Intent(in) :: maxmoments_input
   INTEGER, Intent(in) :: max_user_levels

!  flags
!  Version 1.5: --> F-matrix flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/20/16

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING
   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING
   LOGICAL, Intent(in) :: DO_FMATRIX

   LOGICAL, Intent(in) :: DO_LAMBERTIAN
   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   LOGICAL, Intent(in) :: DO_WATER_LEAVING

   logical, Intent(in) :: DO_Partials
   LOGICAL, Intent(in) :: DO_PLANPAR
   LOGICAL, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Numbers

   INTEGER, Intent(in) :: NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   INTEGER, Intent(in) :: NMOMSINP, NSTOKES
   INTEGER, Intent(in) :: N_USER_LEVELS
   INTEGER, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/20/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX, FVEC(4)

!  Atmosphere. Fmatrix input added 7/7/16

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )
   REAL(ffp), Intent(in) :: FMATRIX_UP  ( MAXLAYERS, MAXGEOMS, 6 )
   REAL(ffp), Intent(in) :: FMATRIX_DN  ( MAXLAYERS, MAXGEOMS, 6 )

!  Surface reflectivity (Could be the albedo)
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( 4, 4, MAXGEOMS )
   real(ffp), Intent(in) :: SLTERM ( 4,    MAXGEOMS )

!  Geometrical inputs
!  ------------------

!  Ray constants, Mu0, Mu1
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only
!   real(ffp), Intent(in)  :: Raycon(maxgeoms)

   real(ffp), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

   REAL(ffp), Intent(in)  :: GenSpher_up(0:maxmoments_input,4,maxgeoms)
   REAL(ffp), Intent(in)  :: GenSpher_dn(0:maxmoments_input,4,maxgeoms)
   REAL(ffp), Intent(in)  :: Rotations_up(4,maxgeoms)
   REAL(ffp), Intent(in)  :: Rotations_dn(4,maxgeoms)

!  Los paths added, 8/20/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/20/16.

   real(ffp), Intent(in)  :: xfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: xfine_dn   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: xfine_dn_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/20/16.

   integer  , Intent(in)  :: ntraverse_up     (0:maxlayers,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_up      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_up (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_up  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_up_p     (maxpartials,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_up_p      (maxpartials,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_up_p (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_up_p  (maxpartials,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_dn     (0:maxlayers,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_dn      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_dn (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_dn  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_dn_p     (maxpartials,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_dn_p      (maxpartials,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_dn_p (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_dn_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

!  4/9/19. Additional output for the sleave correction

   REAL(ffp), Intent(Out)  :: stokes_up     ( max_user_levels, 4, maxgeoms )
   REAL(ffp), Intent(Out)  :: stokes_db     ( max_user_levels, 4, maxgeoms )
   REAL(ffp), Intent(Out)  :: cumsource_up  ( 0:maxlayers, 4, maxgeoms )
   REAL(ffp), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, maxgeoms )
   REAL(ffp), Intent(Out)  :: cumsource_dn  ( 0:maxlayers, 4, maxgeoms )
   
   REAL(ffp), Intent(Out)  :: cumtrans     ( max_user_levels, maxgeoms )

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SSV_Integral_I_UP
!    ==> Add LOSTRANS_DN to the output lists from SSV_Integral_I_DN

   real(ffp), Intent(out) :: lostrans_up (maxgeoms,maxlayers)
   real(ffp), Intent(out) :: lostrans_dn (maxgeoms,maxlayers)

!  Upwelling
!  ---------

   if ( do_upwelling ) then
       call SSV_Integral_I_UP &
        ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,        & ! Inputs (dimension)
          do_sunlight, do_deltam_scaling, do_fmatrix, do_lambertian,                           & ! Inputs (Flags-General/Surface)
          do_surface_leaving, do_water_leaving,                                                & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, flux, fvec, do_sources_up, do_sources_up_p, & ! Inputs(Flags/Flux/criticality)
          nstokes, ngeoms, nlayers, nfinedivs, nmomsinp, n_user_levels, user_levels,           & ! Inputs (control)
          npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
          extinction, deltaus, omega, truncfac, greekmat, fmatrix_up, Reflec, Slterm,          & ! Inputs (Optical/surface)
          Mu0, Mu1, GenSpher_up, Rotations_up, LosW_paths, LosP_paths,                         & ! Inputs (Geometry)
          xfine_up, wfine_up, sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,             & ! Inputs (Geometry)
          xfine_up_p, wfine_up_p, sunpaths_up_p, ntraverse_up_p, sunpathsfine_up_p, ntraversefine_up_p, & ! Inputs (Geometry)
          stokes_up, stokes_db, cumsource_up, cumtrans, lostrans_up )                                     ! Outputs
   endif

   if ( do_dnwelling ) then
       call SSV_Integral_I_DN &
        ( maxgeoms, maxlayers, maxpartials,  maxfine, maxmoments_input, max_user_levels, & ! Inputs (dimension)
          do_sunlight, do_deltam_scaling, do_fmatrix, do_Partials,                       & ! Inputs (Flags)
          do_PlanPar, do_enhanced_ps, flux, fvec, do_sources_dn, do_sources_dn_p,        & ! Inputs (Flags/flux)
          nstokes, ngeoms, nlayers, nfinedivs, nmomsinp, n_user_levels, user_levels,     & ! Inputs (control)
          npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,           & ! Inputs (control-partial)
          extinction, deltaus, omega, truncfac, greekmat, fmatrix_dn,                            & ! Inputs (Optical)
          Mu1, GenSpher_dn, Rotations_dn, LosW_paths, LosP_paths,                                & ! Inputs (Geometry)
          xfine_dn, wfine_dn, sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,             & ! Inputs (Geometry)
          xfine_dn_p, wfine_dn_p, sunpaths_dn_p, ntraverse_dn_p, sunpathsfine_dn_p, ntraversefine_dn_p, & ! Inputs (Geometry)
          stokes_dn, cumsource_dn, lostrans_dn )                                                          ! Outputs
   endif

!  Finish

   return
end subroutine SSV_Integral_I_UPDN

!  End module

end module FO_VectorSS_RTCalcs_I_m
