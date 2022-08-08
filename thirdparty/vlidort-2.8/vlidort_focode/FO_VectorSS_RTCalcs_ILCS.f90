
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

module FO_VectorSS_RTCalcs_ILCS_m

!  For a given wavelength, this routine will calculate upwelling and downwelling
!  First Order Stokes vectors, and any number of LCS Jacobians (column/surface)

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.

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
!    Version 5,  22 August   2016, Partial-layer output

!  For Solar sources, the subroutines are
!       SSV_Integral_ILCS_UP   (Upwelling only)
!       SSV_Integral_ILCS_DN   (Downwelling only)
!       SSV_Integral_ILCS_UPDN (Upwelling and Downwelling)

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP, LC_LOSTRANS_UP to the output lists from SSV_Integral_ILCS_UP, SSV_Integral_ILCS_UPDN
!    ==> Add LOSTRANS_DN, LC_LOSTRANS_DN to the output lists from SSV_Integral_ILCS_DN, SSV_Integral_ILCS_UPDN

!  All subroutines public

public

contains

subroutine SSV_Integral_ILCS_UP &
   ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input,                     & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs, max_sleavewfs,                    & ! Inputs (dimensioning)
     do_sunlight, do_deltam_scaling, do_fmatrix, do_lambertian, do_surface_leaving,   & ! Inputs (Flags - General/Surface)
     do_water_leaving, do_Partials, do_PlanPar, do_enhanced_ps, flux, fvec,           & ! Inputs (Flags - Geometry)
     do_sources_up, do_sources_up_p,do_columnwfs, do_surfacewfs, do_sleavewfs,        & ! Inputs (Flags - Linearz)
     n_reflecwfs, n_sleavewfs, n_surfacewfs, n_columnwfs, Lvarymoms,                  & ! Inputs (Control, Jacobian)
     nstokes, ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels, & ! Inputs (control, output)
     npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,     & ! Inputs (control-partial)
     extinction, deltaus, omega, truncfac, greekmat, fmatrix_up, reflec, slterm,      & ! Inputs (Flux/Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat, L_fmatrix_up,          & ! Inputs (Linearized)
     LS_reflec, LSSL_slterm, Mu0, Mu1, GenSpher, Rotations, LosW_paths, LosP_paths,   & ! Inputs (Geometry)
     xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                  & ! Inputs (Geometry)
     xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,      & ! Inputs (Geometry)
     Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db,         & ! Output
     cumtrans, lostrans_up, LC_cumtrans, LC_Lostrans_up )                               ! Output

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Stokes vectors and LCS Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using F Matrices directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/22/16
!    - Partial-layer output introduced
!    - 4/9/19. Add the CUMTRANS output, add water-leaving control

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP, LC_LosTRANS_UP to the output lists from SSV_Integral_ILCS_UP

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions. Max_sleavewfs added, 8/2/16

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials

   integer, Intent(in) :: maxfine
   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels

   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs
   INTEGER, Intent(in) :: max_sleavewfs

!  flags
!  Version 1.5: --> F-matrix Flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/22/16

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING
   LOGICAL, Intent(in) :: DO_FMATRIX

   LOGICAL, Intent(in) :: DO_LAMBERTIAN
   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   logical, Intent(in) :: DO_WATER_LEAVING    ! 4/9/19 added

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX, fvec(4)

!  Existence flags. 8/22/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)

!  Jacobian Flags. do_sleavewfs added 8/2/16

   LOGICAL, Intent(in) :: do_columnwfs
   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_sleavewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT, NSTOKES
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/22/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  Jacobian control. Reflec and sleave numbers added, 8/2/16
!    Note that n_surfacewfs = n_reflecwfs + n_sleavewfs

   INTEGER, Intent(in) :: n_reflecwfs
   INTEGER, Intent(in) :: n_sleavewfs
   INTEGER, Intent(in) :: n_surfacewfs
   INTEGER, Intent(in) :: n_columnwfs
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere. Fmatrix input added 7/7/16

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )
   REAL(ffp), Intent(in) :: FMATRIX_UP  ( MAXLAYERS, MAXGEOMS, 6 )

!  Surface reflectivity (Could be the albedo) + linearizations
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( 4, 4, MAXGEOMS )
   real(ffp), Intent(in) :: SLTERM ( 4,    MAXGEOMS )
   real(ffp), Intent(in) :: LS_REFLEC   ( 4, 4, MAXGEOMS, max_surfacewfs )
   real(ffp), Intent(in) :: LSSL_SLTERM ( 4,    MAXGEOMS, max_sleavewfs  )

!  Linearized optical inputs. Fmatrix input added 7/7/16

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )
   REAL(ffp), Intent(in) :: L_FMATRIX_UP  ( MAXLAYERS, MAXGEOMS, 6, max_atmoswfs )

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

!  Los paths added, 8/22/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/22/16.

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/22/16.

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
   real(ffp), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, 4, maxgeoms, max_atmoswfs )
   real(ffp), Intent(Out)  :: LC_Jacobians_db  ( max_user_levels, 4, maxgeoms, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, 4, maxgeoms, max_surfacewfs )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS    ( max_user_levels, maxgeoms )
   real(ffp), Intent(out)  :: LC_CUMTRANS ( max_user_levels, maxgeoms, max_atmoswfs )

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP, LC_LosTRANS_UP to the output lists from SSV_Integral_ILCS_UP

   real(ffp), Intent(out)  :: lostrans_up    (maxgeoms,maxlayers)
   real(ffp), Intent(out)  :: LC_Lostrans_up (maxgeoms,maxlayers,max_atmoswfs)

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/22/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: LC_attenuations   (0:maxlayers,max_atmoswfs)

   real(ffp)  :: attenuations_p      (maxpartials)
   real(ffp)  :: LC_attenuations_p   (maxpartials,max_atmoswfs)

   real(ffp)  :: Attenuationsfine    (maxlayers,maxfine)
   real(ffp)  :: LC_Attenuationsfine (maxlayers,maxfine,max_atmoswfs)

   real(ffp)  :: Attenuationsfine_p    (maxpartials,maxfine)
   real(ffp)  :: LC_Attenuationsfine_p (maxpartials,maxfine,max_atmoswfs)

!  Scattering

   real(ffp)  :: tms            (maxlayers)
   real(ffp)  :: exactscat_up   (maxlayers,4,4)
   real(ffp)  :: L_tms          (maxlayers,max_atmoswfs)
   real(ffp)  :: L_exactscat_up (maxlayers,4,4, max_atmoswfs)

!  Source function integration results
!    -- 1/31/21. Version 2.8.3. removed LOSTRANS_UP, LC_Lostrans_up from this list, now outputs

   real(ffp)  :: sources_up  (maxlayers,4), sources_up_p  (maxpartials,4)
   real(ffp)  :: lostrans_up_p (maxpartials)
   real(ffp)  :: LC_sources_up  (maxlayers,4,max_atmoswfs), LC_sources_up_p  (maxpartials,4,max_atmoswfs)
   real(ffp)  :: LC_lostrans_up_p (maxpartials,max_atmoswfs)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: LC_multiplier    ( maxlayers, max_atmoswfs )
   real(ffp)  :: multiplier_p     ( maxpartials )
   real(ffp)  :: LC_multiplier_p  ( maxpartials, max_atmoswfs )

!  Local cumulative source terms

   real(ffp)  :: cumsource_db      ( 0:maxlayers, 4 )
   real(ffp)  :: cumsource_up      ( 0:maxlayers, 4 )

   real(ffp)  :: L_cumsource       ( 4, max_atmoswfs )
   real(ffp)  :: LS_cumsource      ( 4, max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, ns, j, q, q1, L, v, o1, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), nt, np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers), Qvary(maxlayers)

   real(ffp)  :: help, sum, kn, tran, factor1, factor2, m4, m4a, rhelp(4), shelp(4), pi4, term1(4)
   real(ffp)  :: L_help, L_sum(max_atmoswfs), L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd, dj, path_up
   real(ffp)  :: lostau, L_lostau, fmat(6), L_fmat(6,max_atmoswfs), sum23, dif23, L_Shelp

   real(ffp)  :: help3c1, help3s1, help4c1, help4s1
   real(ffp)  :: suntau(0:maxlayers), suntau_p(maxpartials)
   real(ffp)  :: LC_suntau(0:maxlayers,max_atmoswfs),LC_suntau_p(maxpartials,max_atmoswfs)
   real(ffp)  :: ctrans, LC_ctrans(max_atmoswfs)

   real(ffp), parameter  :: cutoff = 88.0_ffp
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Number

!mick fix 9/19/2017 - define pi4 as in VLIDORT_PARS
   !pi4 = acos(-one)/4.0_ffp
   pi4 = acos(-one)*4.0_ffp

!  Zero the output. 4/9/19 include CUMTRANS and LC_CUMTRANS

   STOKES_UP       = zero ; STOKES_DB       = zero ; cumtrans = zero
   LC_JACOBIANS_UP = zero ; LC_JACOBIANS_DB = zero ; LS_JACOBIANS_DB = zero ; LC_cumtrans = zero

!  1/31/21. Version 2.8.3. lostrans_up, LC_Lostrans_up, zeroed here because extra maxgeoms dimension

   lostrans_up = zero ; LC_Lostrans_up = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes ; fmat = zero

!mick fix 3/22/2017  - turned on all values of LAYERMASK_UP
   !NUT = USER_LEVELS(1) + 1
   !LAYERMASK_UP = .false.
   !LAYERMASK_UP(NUT:NLAYERS) = .true.
   LAYERMASK_UP = .true.

!  Linearization bookkeeping

   Qvary = .false. ; QNums = 0
   if ( do_columnwfs ) then
      Qvary(1:nlayers) = .true.
      QNums(1:nlayers) =  n_columnwfs
   endif

!  TMS factors and linearizations

   if ( do_deltam_scaling ) then
      do n = 1, nlayers
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
            enddo
         endif
      enddo
   else
      do n = 1, nlayers
         tms(n) = omega(n)
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_tms(n,q) = L_omega(n,q)
            enddo
         endif
      enddo
   endif

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero local sources
!    -- 1/31/21. Version 2.8.3. removed zeroing of lostrans_up, LC_Lostrans_up

! mick mod 3/22/2017 - turned off (not needed)
! mick fix 9/19/2017 - turned "L_exactscat_up" back on

      !sources_up    = zero ; exactscat_up   = zero ; cumsource_up = zero
      !LC_sources_up = zero

      L_exactscat_up = zero

      !lostrans_up_p    = zero  ; sources_up_p    = zero
      !LC_lostrans_up_p = zero  ; LC_sources_up_p = zero 

!  Scattering functions and Linearization
!  ======================================

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
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) L_fmat(1,q) = L_fmatrix_up(n,v,1,q)
                enddo
              endif
            else
              fmat(1) = zero
              do L = 0, nmoments_input
                fmat(1) = fmat(1) + GenSpher(L,1,V) * Greekmat(n,L,1,1)
              enddo
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                    L_fmat(1,q) = zero
                    do L = 0, nmoments_input
                      L_fmat(1,q) = L_fmat(1,q) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                    enddo
                  endif
                enddo
              endif
            endif

!  Multiply by TMS

            exactscat_up(n,1,1) = fmat(1) * tms(n)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( Lvarymoms(n,q) ) then
                  L_exactscat_up(n,1,1,q) = L_fmat(1,q) * tms(n) + fmat(1) * L_tms(n,q)
                else
                  L_exactscat_up(n,1,1,q) = fmat(1) * L_tms(n,q)
                endif
              enddo
            endif

!  End calculation
              
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
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) L_fmat(1:2,q) = L_fmatrix_up(n,v,1:2,q)
                enddo
              endif
            else
              fmat(1:2) = zero
              do L = 0, nmoments_input
                fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
              enddo
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                    L_fmat(1:2,q) = zero
                    do L = 0, nmoments_input
                       L_fmat(1,q) = L_fmat(1,q) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                       L_fmat(2,q) = L_fmat(2,q) + GenSpher(L,2,V)  * L_Greekmat(n,L,1,2,q)
                    enddo
                  endif
                enddo
              endif
            endif

!   Apply rotations for Z-matrix, multiply by TMS

            exactscat_up(n,1,1) = + fmat(1)
            exactscat_up(n,2,1) = - fmat(2) * Rotations(3,V)
            exactscat_up(n,3,1) = + fmat(2) * Rotations(4,V)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( Lvarymoms(n,q) ) then
                  L_exactscat_up(n,1,1,q) = + L_fmat(1,q)
                  L_exactscat_up(n,2,1,q) = - L_fmat(2,q) * Rotations(3,V)
                  L_exactscat_up(n,3,1,q) = + L_fmat(2,q) * Rotations(4,V)
                endif
                L_exactscat_up(n,1:ns,1,q) = L_exactscat_up(n,1:ns,1,q) *   tms(n) &
                                             + exactscat_up(n,1:ns,1)   * L_tms(n,q)
              enddo
            endif               
            exactscat_up(n,1:ns,1) = tms(n) * exactscat_up(n,1:ns,1) 

!  End calculation

          endif
        enddo
      endif

!  Vector General case
!  -------------------

! USE FULL 4X4 MATRIX; CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010
      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_up(n) ) then

!  Get F-matrix

            if ( do_fmatrix ) then
              fmat(1:6) = fmatrix_up(n,v,1:6)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) L_fmat(1:6,q) = L_fmatrix_up(n,v,1:6,q)
                enddo
              endif
            else
              fmat = zero
              do L = 0, nmoments_input
                fmat(1) = fmat(1) + Genspher(L,1,V) * greekmat(n,L,1,1)
                fmat(2) = fmat(2) + Genspher(L,2,V) * greekmat(n,L,1,2)
                sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                fmat(3) = fmat(3) + Genspher(L,3,V) * sum23
                fmat(4) = fmat(4) + Genspher(L,4,V) * dif23
              enddo
              fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_ffp
              fmat(4) = ( fmat(3) - fmat(4) )
              if ( nstokes.eq.4) then
                do L = 0, nmoments_input
                  fmat(5) = fmat(5) + Genspher(L,2,V) * greekmat(n,L,3,4)
                  fmat(6) = fmat(6) + Genspher(L,1,V) * greekmat(n,L,4,4)
                enddo
              endif
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                    L_fmat(:,q) = zero
                    do L = 0, nmoments_input
                      L_fmat(1,q) = L_fmat(1,q) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                      L_fmat(2,q) = L_fmat(2,q) + GenSpher(L,2,V)  * L_Greekmat(n,L,1,2,q)
                      sum23 = L_greekmat(n,L,2,2,q) + L_greekmat(n,L,3,3,q)
                      dif23 = L_greekmat(n,L,2,2,q) - L_greekmat(n,L,3,3,q)
                      L_fmat(3,q) = L_fmat(3,q) + Genspher(L,3,V) * sum23
                      L_fmat(4,q) = L_fmat(4,q) + Genspher(L,4,V) * dif23
                    enddo
                    L_fmat(3,q) = ( L_fmat(3,q) + L_fmat(4,q) ) * 0.5_ffp
                    L_fmat(4,q) = ( L_fmat(3,q) - L_fmat(4,q) )
                    if ( nstokes.eq.4) then
                      do L = 0, nmoments_input
                        L_fmat(5,q) = L_fmat(5,q) + Genspher(L,2,V) * L_greekmat(n,L,3,4,q)
                        L_fmat(6,q) = L_fmat(6,q) + Genspher(L,1,V) * L_greekmat(n,L,4,4,q)
                      enddo
                    endif
                  endif
                enddo
              endif 
            endif

!   Apply rotations for Z-matrix, multiply by TMS

            help3c1 = fmat(3) * Rotations(1,V)
            help3s1 = fmat(3) * Rotations(2,V)
            help4c1 = fmat(4) * Rotations(1,V)
            help4s1 = fmat(4) * Rotations(2,V)
            exactscat_up(n,1,1) = + fmat(1)
            exactscat_up(n,2,1) = - fmat(2) * Rotations(3,V)
            exactscat_up(n,3,1) = + fmat(2) * Rotations(4,V)
            exactscat_up(n,2,2) = + help3c1 * Rotations(3,V) - help4s1 * Rotations(4,V)
            exactscat_up(n,2,3) = - help3s1 * Rotations(3,V) - help4c1 * Rotations(4,V)
            exactscat_up(n,3,2) = + help3c1 * Rotations(4,V) + help4s1 * Rotations(3,V)
            exactscat_up(n,3,3) = - help3s1 * Rotations(4,V) + help4c1 * Rotations(3,V)
            if ( nstokes .eq. 4 ) then
              exactscat_up(n,2,4) = - fmat(5) * Rotations(4,V) 
              exactscat_up(n,4,2) = - fmat(5) * Rotations(2,V) 
              exactscat_up(n,3,4) = + fmat(5) * Rotations(3,V) 
              exactscat_up(n,4,3) = - fmat(5) * Rotations(1,V) 
              exactscat_up(n,4,4) = + fmat(6)
            endif
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( Lvarymoms(n,q) ) then
                  help3c1 = L_fmat(3,q) * Rotations(1,V)
                  help3s1 = L_fmat(3,q) * Rotations(2,V)
                  help4c1 = L_fmat(4,q) * Rotations(1,V)
                  help4s1 = L_fmat(4,q) * Rotations(2,V)
                  L_exactscat_up(n,1,1,q) = + L_fmat(1,q)
                  L_exactscat_up(n,2,1,q) = - L_fmat(2,q) * Rotations(3,V)
                  L_exactscat_up(n,3,1,q) = + L_fmat(2,q) * Rotations(4,V)
                  L_exactscat_up(n,2,2,q) = + help3c1 * Rotations(3,V) - help4s1 * Rotations(4,V)
                  L_exactscat_up(n,2,3,q) = - help3s1 * Rotations(3,V) - help4c1 * Rotations(4,V)
                  L_exactscat_up(n,3,2,q) = + help3c1 * Rotations(4,V) + help4s1 * Rotations(3,V)
                  L_exactscat_up(n,3,3,q) = - help3s1 * Rotations(4,V) + help4c1 * Rotations(3,V)
                  if ( nstokes .eq. 4 ) then
                    L_exactscat_up(n,2,4,q) = - L_fmat(5,q) * Rotations(4,V) 
                    L_exactscat_up(n,4,2,q) = - L_fmat(5,q) * Rotations(2,V) 
                    L_exactscat_up(n,3,4,q) = + L_fmat(5,q) * Rotations(3,V) 
                    L_exactscat_up(n,4,3,q) = - L_fmat(5,q) * Rotations(1,V) 
                    L_exactscat_up(n,4,4,q) = + L_fmat(6,q)
                  endif
                endif
                L_exactscat_up(n,1:ns,1:ns,q) = L_exactscat_up(n,1:ns,1:ns,q) *   tms(n) &
                                                + exactscat_up(n,1:ns,1:ns)   * L_tms(n,q)
              enddo
            endif               
            exactscat_up(n,1:ns,1:ns) = tms(n)*exactscat_up(n,1:ns,1:ns)

!  End calculation

          endif
        enddo
      endif

!  Attenuations
!  ============

!  Initialize

      Attenuations   = zero ; Attenuationsfine    = zero ; suntau = zero
      Attenuations_p = zero ; Attenuationsfine_p  = zero ; suntau_p = zero

      LC_Attenuations   = zero ; LC_Attenuationsfine    = zero ; LC_suntau = zero
      LC_Attenuations_p = zero ; LC_Attenuationsfine_p  = zero ; LC_suntau_p = zero

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = ntraverse(n,v) ; sumd = dot_product(extinction(1:nt),sunpaths(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_columnwfs ) then
           do q = 1, n_columnwfs
             L_sumd = dot_product(L_extinction(1:nt,q),sunpaths(n,1:nt,v))
             LC_suntau(n,q) = L_sumd ; LC_Attenuations(n,q) = - Attenuations(n) * L_sumd
           enddo
         endif
      enddo

!  RobFix 8/22/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = ntraverse_p(ut,v) ; sumd = dot_product(extinction(1:nt),sunpaths_p(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. cutoff ) Attenuations_p(ut) = exp( - sumd )
          if ( do_columnwfs ) then
            do q = 1, n_columnwfs
              L_sumd = dot_product(L_extinction(1:nt,q),sunpaths_p(ut,1:nt,v))
              LC_suntau_p(ut,q) = L_sumd ; LC_Attenuations_p(ut,q) = - Attenuations_p(ut) * L_sumd
            enddo
          endif
        enddo
      endif

!  Enhanced-spherical, fine-layer attenuations, Whole-layer integration

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
            do j = 1, nfinedivs(n,v)
              nt = ntraversefine(n,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine(n,1:nt,j,v))
              if (sumd .lt. cutoff ) Attenuationsfine(n,j) = exp( - sumd )
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_sumd = dot_product(L_extinction(1:nt,q),sunpathsfine(n,1:nt,j,v))
                  LC_Attenuationsfine(n,j,q) = - Attenuationsfine(n,j) * L_sumd 
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  RobFix 8/22/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, nfinedivs_p(ut,v)
              nt = ntraversefine_p(ut,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine_p(ut,1:nt,j,v))
              If (sumd .lt. cutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_sumd = dot_product(L_extinction(1:nt,q),sunpathsfine_p(ut,1:nt,j,v))
!mick fix 9/19/2017 - replaced line below
                  !LC_Attenuationsfine(ut,j,q) = - Attenuationsfine_p(ut,j) * L_sumd
                  LC_Attenuationsfine_p(ut,j,q) = - Attenuationsfine_p(ut,j) * L_sumd
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

!  1/31/21. Version 2.8.3. lostrans_up, LC_Lostrans_up given geometry index, as now they are outputs

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_up(n) .and. do_sources_up(n,v)  ) then

 !  Sources, general case

            if ( Mu1(v) .gt. zero ) then
              lostau = deltaus(n)  / Mu1(v)
              if ( lostau .lt. cutoff ) lostrans_up(v,n) = exp( - lostau )
              factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(v,n)
              factor2 = one + (suntau(n) - suntau(n-1))/lostau
              multiplier(n) = factor1 / factor2
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_lostau           = L_deltaus(n,q) / Mu1(v)
                  LC_lostrans_up(v,n,q) = - L_lostau * lostrans_up(v,n)
                  L_factor1 = LC_Attenuations(n-1,q) - LC_Attenuations(n,q)*lostrans_up(v,n)
                  L_factor2 = (LC_suntau(n,q) - LC_suntau(n-1,q))/lostau
                  L_factor1 = L_factor1 - Attenuations(n)*LC_lostrans_up(v,n,q)
                  L_factor2 = L_factor2 - (factor2 - one)*L_lostau/lostau 
                  LC_multiplier(n,q) = ( L_factor1 - multiplier(n)*L_factor2 ) / factor2
                enddo
              endif
            endif

!  End whole layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  RobFix 8/22/16. Plane/Parallel or Regular-PS, Partial-layer output
!  ------------------------------------------------------------------

      if ( do_RegPSorPP .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then          
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v) - LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero

 !  Sources, general case

            if ( Mu1(v) .gt. zero ) then
              lostau = kn * path_up
              if ( lostau .lt. cutoff ) lostrans_up_p(ut) = exp( - lostau )
              factor1 = Attenuations_p(ut) - Attenuations(np)*lostrans_up_p(ut)
              factor2 = one + (suntau(np) - suntau_p(ut))/lostau
              multiplier_p(ut) = factor1 / factor2
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_lostau = L_extinction(np,q) * path_up
                  LC_lostrans_up_p(ut,q) = - L_lostau * lostrans_up_p(ut)
                  L_factor1 = LC_Attenuations_p(ut,q) - LC_Attenuations(np,q)*lostrans_up_p(ut)
                  L_factor2 = (LC_suntau(np,q) - LC_suntau_p(ut,q))/lostau
                  L_factor1 = L_factor1 - Attenuations(np)*LC_lostrans_up_p(ut,q)
                  L_factor2 = L_factor2 - (factor2 - one)*L_lostau/lostau 
                  LC_multiplier_p(ut,q) = ( L_factor1 - multiplier_p(ut)*L_factor2 ) / factor2
                enddo
              endif
            endif

!  End partial layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/22/16 streamlined code using distances
!      Quadratures from Bottom of the layer

!  1/31/21. Version 2.8.3. lostrans_up, LC_Lostrans_up given geometry index, as now they are outputs

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_up(n) .and. do_sources_up(n,v)  ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_up = LosW_paths(n,v)
            lostau = kn * path_up ; if( lostau.lt.cutoff ) lostrans_up(v,n) = exp ( - lostau )
            if ( do_columnwfs ) then
              do q = 1, Qnums(n)
                L_lostau        = L_extinction(n,q) * path_up
                LC_lostrans_up(v,n,q) = - L_lostau * lostrans_up(v,n)
              enddo
            endif
            sum = zero ; L_sum = zero 
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * wfine(n,j,v)
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_tran = - dj * L_extinction(n,q)
                  L_func = ( LC_attenuationsfine(n,j,q) + L_tran * attenuationsfine(n,j) ) * tran
                  L_sum(q) = L_sum(q) + L_func * wfine(n,j,v)
                enddo
              endif
            enddo
            multiplier(n) = sum * kn
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                LC_multiplier(n,q) = kn * L_sum(q) + sum * L_extinction(n,q)
              enddo
            endif        
          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/22/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v)- LosP_paths(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.cutoff ) lostrans_up_p(ut) = exp ( - lostau )

            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                L_lostau = L_extinction(np,q) * path_up
                LC_lostrans_up_p(ut,q) = - L_lostau * lostrans_up_p(ut)
              enddo
            endif

            sum = zero ; L_sum = zero 
            do j = 1, nfinedivs_p(ut,v)
              dj = path_up - xfine_p(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine_p(ut,j) * tran * wfine_p(ut,j,v)
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_tran = - dj * L_extinction(np,q)
                  L_func = ( LC_attenuationsfine_p(ut,j,q) + L_tran * attenuationsfine_p(ut,j) ) * tran
                  L_sum(q) = L_sum(q) + L_func * wfine_p(ut,j,v)
                enddo
              endif
            enddo
            multiplier_p(ut) = sum * kn

            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                LC_multiplier_p(ut,q) = kn * L_sum(q) + sum * L_extinction(np,q)
              enddo
            endif

          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  General case, Whole layers
!mick fix 9/19/2017 - defined "sources_up" for all stokes vector elements

      do n = nlayers, 1, -1
        if ( layermask_up(n) .and.do_sources_up(n,v)  ) then
          if ( do_sunlight ) then
             shelp(1:ns) = exactscat_up(n,1:ns,1) * fvec(1)
          else
            do o1 = 1, nstokes
              shelp(o1) = dot_product(exactscat_up(n,o1,1:ns),fvec(1:ns))
            enddo
          endif
          !sources_up(n,o1) = shelp(o1) * multiplier(n)
          sources_up(n,1:ns) = shelp(1:ns) * multiplier(n)

          if ( do_columnwfs ) then
            do q = 1, n_columnwfs
              do o1 = 1, nstokes
                LC_sources_up(n,o1,q) = shelp(o1) * LC_multiplier(n,q)
                if ( do_sunlight ) then
                  L_Shelp = L_exactscat_up(n,o1,1,q)*fvec(1)
                else
                  L_Shelp = dot_product(L_exactscat_up(n,o1,1:ns,q),fvec(1:ns))
                endif
                LC_sources_up(n,o1,q) =  LC_sources_up(n,o1,q) + L_Shelp * multiplier(n)
              enddo
            enddo
          endif
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

            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                do o1 = 1, nstokes !index o2 in SSCOR
                  LC_sources_up_p(ut,o1,q) = shelp(o1) * LC_multiplier_p(ut,q)
                  if ( do_sunlight ) then
                    L_Shelp = L_exactscat_up(np,o1,1,q)*fvec(1)
                  else
                    L_Shelp = dot_product(L_exactscat_up(np,o1,1:ns,q),fvec(1:ns))
                  endif
                  LC_sources_up_p(ut,o1,q) =  LC_sources_up_p(ut,o1,q) + L_Shelp * multiplier_p(ut)
                enddo
              enddo
            endif
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  Stokes-vector Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0 
      CUMSOURCE_UP(NC,:) = zero
      CUMSOURCE_DB(NC,:) = zero
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Surface term

      RHELP = zero; M4 = 4.0_ffp * MU0(v) ; M4A = M4 * attenuations(nlayers)
      if ( DO_LAMBERTIAN ) then
         RHELP(1) = M4 * REFLEC(1,1,V) * fvec(1)
         CUMSOURCE_DB(NC,1) = RHELP(1) * attenuations(nlayers)
      else
         do o1 = 1, nstokes
            RHELP(O1) = M4 * dot_product(REFLEC(O1,1:ns,V),fvec(1:ns))
            CUMSOURCE_DB(NC,o1) = RHELP(O1) * attenuations(nlayers)
         enddo
      endif

!  surface-leaving term. Added, 8/2/16
!   -- (modeled after the DBCORRECTION code in Version 2.7)
!   -- 4/9/19. Not done for water-leaving, as need to use adjusted values

      IF ( DO_SURFACE_LEAVING .and. .not. DO_WATER_LEAVING ) THEN
         do o1 = 1, nstokes
            CUMSOURCE_DB(NC,o1) = CUMSOURCE_DB(NC,o1) + PI4 * SLTERM(o1,v)
         enddo
      ENDIF

!  Main Stokes vector loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. RobFix 8/22/16 Partials.

!  1/31/21. Version 2.8.3. lostrans_up given geometry index, as now output

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            do o1 = 1, nstokes
               CUMSOURCE_DB(NC,O1) = lostrans_up(v,n) * CUMSOURCE_DB(NC-1,O1)
               CUMSOURCE_UP(NC,O1) = lostrans_up(v,n) * CUMSOURCE_UP(NC-1,O1) + SOURCES_UP(N,O1)
            enddo
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           STOKES_UP(UTA,1:NS,V) = FLUX * ( CUMSOURCE_UP(NC,1:NS) * LOSTRANS_UP_p(UT) + SOURCES_UP_p(UT,1:NS) )
           STOKES_DB(UTA,1:NS,V) = FLUX * CUMSOURCE_DB(NC,1:NS) * LOSTRANS_UP_p(UT)
         ELSE
           STOKES_UP(UTA,1:NS,V) = FLUX * CUMSOURCE_UP(NC,1:NS)
           STOKES_DB(UTA,1:NS,V) = FLUX * CUMSOURCE_DB(NC,1:NS)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
      ENDDO

!  Surface WFs. Change to n_reflecwfs, 8/2/16
!mick fix 9/19/2017 - swapped indices for last two dimemsions in LS_REFLEC for
!                     non-Lambertain case

      if ( do_surfacewfs ) then
         LS_CUMSOURCE = zero
         if ( DO_LAMBERTIAN ) then
            LS_cumsource(1,1:n_reflecwfs) = M4A * LS_REFLEC(1,1,v,1:n_reflecwfs)
         else
            do q = 1, n_reflecwfs
               do o1 = 1, nstokes
                  LS_cumsource(o1,q) = M4A * dot_product(LS_REFLEC(o1,1:ns,v,q),fvec(1:ns))
               enddo
            enddo
         endif
       endif

!  Sleave WFs. This section added, 8/2/16
!   -- (modeled after the LSSL_DBCORRECTION code in Version 2.7)

!  5/24/21. Version 2.8.3. Rob Fix - Add ".not.DO_WATER_LEAVING" to the condition.

!      if ( do_surface_leaving .and. do_sleavewfs ) then
      if ( do_surface_leaving .and. .not. do_water_leaving .and. do_sleavewfs ) then
         do q = 1, n_sleavewfs
            q1  = q + n_reflecwfs
            do o1 = 1, nstokes
               LS_cumsource(o1,q1) = pi4 * lssl_slterm(o1,v,q)
            enddo
         enddo
      endif

!  Propagation of surface+sleave weighting functions

!  1/31/21. Version 2.8.3. lostrans_up given geometry index, as now output

      if ( do_surfacewfs .or. ( do_surface_leaving.and.do_sleavewfs) ) then
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               do q = 1, n_surfacewfs
                  LS_cumsource(1:ns,q) = lostrans_up(v,n) * LS_CUMSOURCE(1:ns,q)
               enddo
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_surfacewfs
                LS_JACOBIANS_DB(UTA,1:NS,V,Q) = FLUX * LS_CUMSOURCE(1:NS,Q) * LOSTRANS_UP_P(UT)
              enddo
            ELSE
              do q = 1, n_surfacewfs
                LS_JACOBIANS_DB(UTA,1:NS,V,Q) = FLUX * LS_CUMSOURCE(1:NS,Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
         ENDDO
      endif

!  Column Wfs (Atmospheric term)
!  1/31/21. Version 2.8.3. lostrans_up and LC_Lostrans_up given geometry index, as now output

      if ( do_columnwfs ) then
!mick fix 3/22/2017 - initialize NC
         NC = 0
         L_CUMSOURCE = zero
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) = LC_SOURCES_UP(N,1:ns,Q)         + &
                                LC_lostrans_up(v,n,q) * CUMSOURCE_UP(NC-1,1:ns) + &
                                   lostrans_up(v,n)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_columnwfs
                Term1(1:NS) = CUMSOURCE_UP(NC,1:ns) * LC_LOSTRANS_UP_P(UT,q) + &
                              L_CUMSOURCE(1:ns,Q) * LOSTRANS_UP_p(UT)
                LC_JACOBIANS_UP(UTA,1:ns,V,Q) = FLUX * ( TERM1(1:NS) + LC_SOURCES_UP_p(UT,1:ns,Q) )
              enddo
            ELSE
              do q = 1, n_columnwfs
                LC_JACOBIANS_UP(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Column Wfs (Direct beam term)
!  1/31/21. Version 2.8.3. lostrans_up and LC_Lostrans_up given geometry index, as now output

      if ( do_columnwfs ) then
!mick fix 3/22/2017 - initialize NC
         NC = 0
         do q = 1, n_columnwfs
            L_CUMSOURCE(1:ns,q) = RHELP(1:ns) * LC_attenuations(nlayers,q)
         enddo
         NSTART = NLAYERS  ;  NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) =  LC_lostrans_up(v,n,q) * CUMSOURCE_DB(NC-1,1:ns) + &
                                            lostrans_up(v,n)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_columnwfs
                Term1(1:NS) = CUMSOURCE_DB(NC,1:ns) * LC_LOSTRANS_UP_P(UT,q) + &
                                 L_CUMSOURCE(1:ns,Q) * LOSTRANS_UP_p(UT)
                LC_JACOBIANS_DB(UTA,1:ns,V,Q) = FLUX * TERM1(1:NS) 
              enddo
            ELSE
              do q = 1, n_columnwfs
                LC_JACOBIANS_DB(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
         ENDDO
      endif

!  4/9/19. WATERLEAVING CASE. Add CUMTRANS calculation
!    Add start of CUMTRANS recursion (CTRANS = 1.0).

!  1/31/21. Version 2.8.3. (5/5/20. Version 2.8.1 Upgrade)
!    -- CUMTRANS calculation was NOT properly initialized (no Jacobians)
!             BUT was properly initialized for the Column jacobian output.
   
!  1/31/21. Version 2.8.3. lostrans_up and LC_Lostrans_up given geometry index, as now output

      if ( do_water_leaving ) then
         NSTART = NLAYERS                ! 5/5/20 Needed for initializing
         NUT_PREV = NSTART + 1           ! 5/5/20 Needed for initializing
         ctrans = one
         if ( do_columnwfs ) then
            Q = n_columnwfs ; LC_Ctrans = zero
            DO UTA = N_USER_LEVELS, 1, -1
               NUT    = USER_LEVELS(UTA) + 1
               DO N = NSTART, NUT, -1
                  LC_CTRANS(1:q) = LC_CTRANS(1:q) * lostrans_up(v,n) + CTRANS * LC_LOSTRANS_UP(V,N,1:q)
                  CTRANS = CTRANS * lostrans_up(v,n)
               ENDDO
               IF ( Partial_OUTFLAG(UTA) ) THEN
                  UT = Partial_OUTINDEX(UTA)
                  LC_CUMTRANS(uta,v,1:q) = LC_CTRANS(1:q) * LOSTRANS_UP_p(UT) + CTRANS * LC_LOSTRANS_UP_p(UT,1:q)
                  CUMTRANS(UTA,V)     = CTRANS * LOSTRANS_UP_p(UT)
               ELSE
                  CUMTRANS(UTA,V)        = CTRANS
                  LC_CUMTRANS(uta,v,1:q) = LC_CTRANS(1:q)
               ENDIF
               IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
               NUT_PREV = NUT
            ENDDO
         else  
            DO UTA = N_USER_LEVELS, 1, -1
               NUT    = USER_LEVELS(UTA) + 1
               DO N = NSTART, NUT, -1
                  CTRANS = CTRANS * lostrans_up(v,n)
               ENDDO
               IF ( Partial_OUTFLAG(UTA) ) THEN
                  UT = Partial_OUTINDEX(UTA)
                  CUMTRANS(UTA,V) = CTRANS * LOSTRANS_UP_p(UT)
               ELSE
                  CUMTRANS(UTA,V) = CTRANS
               ENDIF
               IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
               NUT_PREV = NUT
            ENDDO
         endif
      endif
      
!  End Geometry Loop

   enddo

!  Finish

   return
end subroutine SSV_Integral_ILCS_UP

!

subroutine SSV_Integral_ILCS_DN &
      ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels, max_atmoswfs, & ! Inputs (dimension)
        do_sunlight, do_deltam_scaling, do_fmatrix, do_Partials, do_PlanPar, do_enhanced_ps,        & ! Inputs (Flags/flux)
        flux, fvec, do_sources_dn, do_sources_dn_p, do_columnwfs, n_columnwfs, Lvarymoms,           & ! Inputs (control, Jacobian )
        nstokes, ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (control)
        npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,                & ! Inputs (control-partial)
        extinction, deltaus, omega, truncfac, greekmat, fmatrix_dn,                  & ! Inputs (Optical)
        L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat, L_fmatrix_dn,      & ! Inputs (Optical - Linearized)
        Mu1, GenSpher, Rotations, LosW_paths, LosP_paths,                            & ! Inputs (Geometry)
        xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,              & ! Inputs (Geometry)
        xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,  & ! Inputs (Geometry)
        Stokes_dn, LC_Jacobians_dn, lostrans_dn, LC_lostrans_dn )                      ! Output

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of stokes vectors and LCS Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using F Matrices directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/22/16
!    - Partial-layer output introduced

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add lostrans_dn, LC_Lostrans_dn to the output lists from SSV_Integral_ILCS_DN

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials

   integer, Intent(in) :: maxfine
   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels

   INTEGER, Intent(in) :: max_atmoswfs

!  flags. F-matrix flag added, 7/7/16

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING
   LOGICAL, Intent(in) :: DO_FMATRIX

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Solar Flux

   real(ffp), Intent(in) :: FLUX, fvec(4)

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT, NSTOKES
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/20/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  Jacobian Flag and control

   LOGICAL, Intent(in) :: do_columnwfs
   INTEGER, Intent(in) :: n_columnwfs
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere. Fmatrix input added 7/7/16

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )
   REAL(ffp), Intent(in) :: FMATRIX_DN  ( MAXLAYERS, MAXGEOMS, 6 )

!  Linearized optical inputs. Fmatrix input added 7/7/16

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )
   REAL(ffp), Intent(in) :: L_FMATRIX_DN  ( MAXLAYERS, MAXGEOMS, 6, max_atmoswfs )

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

!  Los paths added, 8/22/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/22/16.

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/22/16.

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

   real(ffp), Intent(Out)  :: Stokes_dn        ( max_user_levels, 4, maxgeoms )
   real(ffp), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, 4, maxgeoms, max_atmoswfs )

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add lostrans_dn, LC_Lostrans_dn to the output lists from SSV_Integral_ILCS_DN

   real(ffp), Intent(out)  :: lostrans_dn (maxgeoms,maxlayers)
   real(ffp), Intent(out)  :: LC_Lostrans_dn (maxgeoms,maxlayers,Max_atmoswfs)

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/22/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: LC_attenuations   (0:maxlayers,max_atmoswfs)

   real(ffp)  :: attenuations_p     (maxpartials)
   real(ffp)  :: LC_attenuations_p   (maxpartials,max_atmoswfs)

   real(ffp)  :: Attenuationsfine    (maxlayers,maxfine)
   real(ffp)  :: LC_Attenuationsfine (maxlayers,maxfine,max_atmoswfs)

   real(ffp)  :: Attenuationsfine_p    (maxpartials,maxfine)
   real(ffp)  :: LC_Attenuationsfine_p (maxpartials,maxfine,max_atmoswfs)

!  Scattering

   real(ffp)  :: tms            (maxlayers)
   real(ffp)  :: exactscat_dn   (maxlayers,4,4)
   real(ffp)  :: L_tms          (maxlayers,max_atmoswfs)
   real(ffp)  :: L_exactscat_dn (maxlayers,4,4, max_atmoswfs)

!  Source function integration results
!  1/31/21. Version 2.8.3. Removed LOSTRANS_DN and LC_LosTRANS_DN from this list, now outputs

   real(ffp)  :: sources_dn  (maxlayers,4), sources_dn_p  (maxpartials,4)
   real(ffp)  :: lostrans_dn_p (maxpartials)
   real(ffp)  :: LC_sources_dn  (maxlayers,4,max_atmoswfs), LC_sources_dn_p  (maxpartials,4,max_atmoswfs)
   real(ffp)  :: LC_lostrans_dn_p (maxpartials,max_atmoswfs)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: LC_multiplier    ( maxlayers, max_atmoswfs )
   real(ffp)  :: multiplier_p     ( maxpartials )
   real(ffp)  :: LC_multiplier_p  ( maxpartials, max_atmoswfs )

!  Local cumulative source terms

   real(ffp)  :: cumsource_dn      ( 0:maxlayers, 4 )
   real(ffp)  :: L_cumsource       ( 4, max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, ns, j, q, L, v, o1, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), nt, np, ut
   logical    :: do_regular_ps, layermask_dn(maxlayers), Qvary(maxlayers)

   real(ffp)  :: help, sum, kn, tran, factor1, factor2, shelp(4), term1(4)
   real(ffp)  :: L_help, L_sum(max_atmoswfs), L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd, dj, path_dn
   real(ffp)  :: lostau, L_lostau, fmat(6), L_fmat(6,max_atmoswfs), sum23, dif23, L_Shelp

   real(ffp)  :: help3c1, help3s1, help4c1, help4s1
   real(ffp)  :: suntau(0:maxlayers), suntau_p(maxpartials)
   real(ffp)  :: LC_suntau(0:maxlayers,max_atmoswfs),LC_suntau_p(maxpartials,max_atmoswfs)

   real(ffp), parameter  :: cutoff = 88.0_ffp
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Zero the output 

   STOKES_DN       = zero
   LC_JACOBIANS_DN = zero

!  1/31/21. Version 2.8.3. lostrans_dn, LC_Lostrans_dn, zeroed here because extra maxgeoms dimension

   lostrans_dn = zero ; LC_Lostrans_dn = zero

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

!  Linearization bookkeeping

   Qvary = .false. ; QNums = 0
   if ( do_columnwfs ) then
      Qvary(1:nlayers) = .true.
      QNums(1:nlayers) =  n_columnwfs
   endif

!  TMS factors and linearizations

   if ( do_deltam_scaling ) then
      do n = 1, nlayers
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
            enddo
         endif
      enddo
   else
      do n = 1, nlayers
         tms(n) = omega(n)
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_tms(n,q) = L_omega(n,q)
            enddo
         endif
      enddo
   endif

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero local sources
!    -- 1/31/21. Version 2.8.3. removed zeroing of lostrans_dn, LC_Lostrans_dn

      sources_dn    = zero ; exactscat_dn   = zero ; cumsource_dn = zero
      LC_sources_dn = zero ; L_exactscat_dn = zero 

      lostrans_dn_p   = zero  ; sources_dn_p   = zero
      LC_lostrans_dn_p = zero  ; LC_sources_dn_p = zero 

!  Scattering functions and Linearization
!  ======================================

!  Version 1.5, F-matrix option introduced. 7/7/16

!  Scalar only
!  -----------

      if ( nstokes .eq. 1 ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then

!  Get the F matrix

            if ( do_fmatrix ) then
              fmat(1) = fmatrix_dn(n,v,1)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) L_fmat(1,q) = L_fmatrix_dn(n,v,1,q)
                enddo
              endif
            else
              fmat(1) = zero
              do L = 0, nmoments_input
                fmat(1) = fmat(1) + GenSpher(L,1,V) * Greekmat(n,L,1,1)
              enddo
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                    L_fmat(1,q) = zero
                    do L = 0, nmoments_input
                      L_fmat(1,q) = L_fmat(1,q) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                    enddo
                  endif
                enddo
              endif
            endif

!  Multiply by TMS

            exactscat_dn(n,1,1) = fmat(1) * tms(n)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( Lvarymoms(n,q) ) then
                  L_exactscat_dn(n,1,1,q) = L_fmat(1,q)* tms(n) + fmat(1) * L_tms(n,q)
                else
                  L_exactscat_dn(n,1,1,q) = fmat(1) * L_tms(n,q)
                endif
              enddo
            endif

!  done calculation
              
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
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) L_fmat(1:2,q) = L_fmatrix_dn(n,v,1:2,q)
                enddo
              endif
            else
              fmat(1:2) = zero
              do L = 0, nmoments_input
                fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
              enddo
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                    L_fmat(1:2,q) = zero
                    do L = 0, nmoments_input
                       L_fmat(1,q) = L_fmat(1,q) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                       L_fmat(2,q) = L_fmat(2,q) + GenSpher(L,2,V)  * L_Greekmat(n,L,1,2,q)
                    enddo
                  endif
                enddo
              endif
            endif

!   Apply rotations for Z-matrix, multiply by TMS

            exactscat_dn(n,1,1) = + fmat(1)
            exactscat_dn(n,2,1) = - fmat(2) * Rotations(3,V)
            exactscat_dn(n,3,1) = + fmat(2) * Rotations(4,V)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( Lvarymoms(n,q) ) then
                  L_exactscat_dn(n,1,1,q) = + L_fmat(1,q)
                  L_exactscat_dn(n,2,1,q) = - L_fmat(2,q) * Rotations(3,V)
                  L_exactscat_dn(n,3,1,q) = + L_fmat(2,q) * Rotations(4,V)
                endif
                L_exactscat_dn(n,1:ns,1,q) = L_exactscat_dn(n,1:ns,1,q) *   tms(n) &
                                             + exactscat_dn(n,1:ns,1)   * L_tms(n,q)
              enddo
            endif               
            exactscat_dn(n,1:ns,1) = tms(n) * exactscat_dn(n,1:ns,1) 

!  done calculation

          endif
        enddo
      endif

!  Vector General case
!  -------------------

! USE FULL 4X4 MATRIX; CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010
      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then

!  get F-matrix

            if ( do_fmatrix ) then
              fmat(1:6) = fmatrix_dn(n,v,1:6)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) L_fmat(1:6,q) = L_fmatrix_dn(n,v,1:6,q)
                enddo
              endif
            else
              fmat = zero
              do L = 0, nmoments_input
                fmat(1) = fmat(1) + Genspher(L,1,V) * greekmat(n,L,1,1)
                fmat(2) = fmat(2) + Genspher(L,2,V) * greekmat(n,L,1,2)
                sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                fmat(3) = fmat(3) + Genspher(L,3,V) * sum23
                fmat(4) = fmat(4) + Genspher(L,4,V) * dif23
              enddo
              fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_ffp
              fmat(4) = ( fmat(3) - fmat(4) )
              if ( nstokes.eq.4) then
                do L = 0, nmoments_input
                  fmat(5) = fmat(5) + Genspher(L,2,V) * greekmat(n,L,3,4)
                  fmat(6) = fmat(6) + Genspher(L,1,V) * greekmat(n,L,4,4)
                enddo
              endif
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                    L_fmat(:,q) = zero
                    do L = 0, nmoments_input
                      L_fmat(1,q) = L_fmat(1,q) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                      L_fmat(2,q) = L_fmat(2,q) + GenSpher(L,2,V)  * L_Greekmat(n,L,1,2,q)
                      sum23 = L_greekmat(n,L,2,2,q) + L_greekmat(n,L,3,3,q)
                      dif23 = L_greekmat(n,L,2,2,q) - L_greekmat(n,L,3,3,q)
                      L_fmat(3,q) = L_fmat(3,q) + Genspher(L,3,V) * sum23
                      L_fmat(4,q) = L_fmat(4,q) + Genspher(L,4,V) * dif23
                    enddo
                    L_fmat(3,q) = ( L_fmat(3,q) + L_fmat(4,q) ) * 0.5_ffp
                    L_fmat(4,q) = ( L_fmat(3,q) - L_fmat(4,q) )
                    if ( nstokes.eq.4) then
                      do L = 0, nmoments_input
                        L_fmat(5,q) = L_fmat(5,q) + Genspher(L,2,V) * L_greekmat(n,L,3,4,q)
                        L_fmat(6,q) = L_fmat(6,q) + Genspher(L,1,V) * L_greekmat(n,L,4,4,q)
                      enddo
                    endif
                  endif
                enddo
              endif 
            endif

!   Apply rotations for Z-matrix, multiply by TMS

            help3c1 = fmat(3) * Rotations(1,V)
            help3s1 = fmat(3) * Rotations(2,V)
            help4c1 = fmat(4) * Rotations(1,V)
            help4s1 = fmat(4) * Rotations(2,V)
            exactscat_dn(n,1,1) = + fmat(1)
            exactscat_dn(n,2,1) = - fmat(2) * Rotations(3,V)
            exactscat_dn(n,3,1) = + fmat(2) * Rotations(4,V)
            exactscat_dn(n,2,2) = + help3c1 * Rotations(3,V) - help4s1 * Rotations(4,V)
            exactscat_dn(n,2,3) = - help3s1 * Rotations(3,V) - help4c1 * Rotations(4,V)
            exactscat_dn(n,3,2) = + help3c1 * Rotations(4,V) + help4s1 * Rotations(3,V)
            exactscat_dn(n,3,3) = - help3s1 * Rotations(4,V) + help4c1 * Rotations(3,V)
            if ( nstokes .eq. 4 ) then
              exactscat_dn(n,2,4) = - fmat(5) * Rotations(4,V) 
              exactscat_dn(n,4,2) = - fmat(5) * Rotations(2,V) 
              exactscat_dn(n,3,4) = + fmat(5) * Rotations(3,V) 
              exactscat_dn(n,4,3) = - fmat(5) * Rotations(1,V) 
              exactscat_dn(n,4,4) = + fmat(6)
            endif
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( Lvarymoms(n,q) ) then
                  help3c1 = L_fmat(3,q) * Rotations(1,V)
                  help3s1 = L_fmat(3,q) * Rotations(2,V)
                  help4c1 = L_fmat(4,q) * Rotations(1,V)
                  help4s1 = L_fmat(4,q) * Rotations(2,V)
                  L_exactscat_dn(n,1,1,q) = + L_fmat(1,q)
                  L_exactscat_dn(n,2,1,q) = - L_fmat(2,q) * Rotations(3,V)
                  L_exactscat_dn(n,3,1,q) = + L_fmat(2,q) * Rotations(4,V)
                  L_exactscat_dn(n,2,2,q) = + help3c1 * Rotations(3,V) - help4s1 * Rotations(4,V)
                  L_exactscat_dn(n,2,3,q) = - help3s1 * Rotations(3,V) - help4c1 * Rotations(4,V)
                  L_exactscat_dn(n,3,2,q) = + help3c1 * Rotations(4,V) + help4s1 * Rotations(3,V)
                  L_exactscat_dn(n,3,3,q) = - help3s1 * Rotations(4,V) + help4c1 * Rotations(3,V)
                  if ( nstokes .eq. 4 ) then
                    L_exactscat_dn(n,2,4,q) = - L_fmat(5,q) * Rotations(4,V) 
                    L_exactscat_dn(n,4,2,q) = - L_fmat(5,q) * Rotations(2,V) 
                    L_exactscat_dn(n,3,4,q) = + L_fmat(5,q) * Rotations(3,V) 
                    L_exactscat_dn(n,4,3,q) = - L_fmat(5,q) * Rotations(1,V) 
                    L_exactscat_dn(n,4,4,q) = + L_fmat(6,q)
                  endif
                endif
                L_exactscat_dn(n,1:ns,1:ns,q) = L_exactscat_dn(n,1:ns,1:ns,q) *   tms(n) &
                                                + exactscat_dn(n,1:ns,1:ns)   * L_tms(n,q)
              enddo
            endif               
            exactscat_dn(n,1:ns,1:ns) = tms(n)*exactscat_dn(n,1:ns,1:ns)

!  Done calculation

          endif
        enddo
      endif

!  Attenuations
!  ============

!  Initialize

      Attenuations   = zero ; Attenuationsfine    = zero ; suntau = zero
      Attenuations_p = zero ; Attenuationsfine_p  = zero ; suntau_p = zero
      LC_Attenuations   = zero ; LC_Attenuationsfine    = zero ; LC_suntau = zero
      LC_Attenuations_p = zero ; LC_Attenuationsfine_p  = zero ; LC_suntau_p = zero
      nstart = nlayers

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = ntraverse(n,v) ; sumd = dot_product(extinction(1:nt),sunpaths(n,1:nt,v))
         suntau(n) = sumd    ; if (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_columnwfs ) then
           do q = 1, n_columnwfs
             L_sumd = dot_product(L_extinction(1:nt,q),sunpaths(n,1:nt,v))
             LC_suntau(n,q) = L_sumd ; LC_Attenuations(n,q) = - Attenuations(n) * L_sumd
           enddo
         endif
      enddo

!  RobFix 8/22/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = ntraverse_p(ut,v) ; sumd = dot_product(extinction(1:nt),sunpaths_p(ut,1:nt,v))
          suntau_p(ut) = sumd    ; if (sumd .lt. cutoff ) Attenuations_p(ut) = exp( - sumd )
          if ( do_columnwfs ) then
            do q = 1, n_columnwfs
              L_sumd = dot_product(L_extinction(1:nt,q),sunpaths_p(ut,1:nt,v))
              LC_suntau_p(ut,q) = L_sumd ; LC_Attenuations_p(ut,q) = - Attenuations_p(ut) * L_sumd
            enddo
          endif
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
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_sumd = dot_product(L_extinction(1:nt,q),sunpathsfine(n,1:nt,j,v))
                  LC_Attenuationsfine(n,j,q) = - Attenuationsfine(n,j) * L_sumd 
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  RobFix 8/22/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, nfinedivs_p(ut,v)
              nt = ntraversefine_p(ut,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine_p(ut,1:nt,j,v))
              if (sumd .lt. cutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_sumd = dot_product(L_extinction(1:nt,q),sunpathsfine_p(ut,1:nt,j,v))
!mick fix 9/19/2017 - replaced line below
                  !LC_Attenuationsfine(ut,j,q) = - Attenuationsfine_p(ut,j) * L_sumd
                  LC_Attenuationsfine_p(ut,j,q) = - Attenuationsfine_p(ut,j) * L_sumd
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

!  1/31/21. Version 2.8.3. LOSTRANS_DN and LC_LosTRANS_DN given geometry index, as now outputs

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then

!  Sources, general case

            if ( Mu1(v) .gt. zero ) then
              lostau = deltaus(n) / Mu1(v)
              if ( lostau .lt. cutoff ) lostrans_dn(v,n) = exp( - lostau )
              factor1 = Attenuations(n-1)*lostrans_dn(v,n) - Attenuations(n)
              factor2 = ((suntau(n) - suntau(n-1))/lostau)  - one
              multiplier(n) = factor1 / factor2

              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_lostau  = L_deltaus(n,q) / Mu1(v)
                  LC_lostrans_dn(v,n,q) = - L_lostau * lostrans_dn(v,n)
                  L_factor1 = LC_Attenuations(n-1,q)*lostrans_dn(v,n) - LC_Attenuations(n,q)
                  L_factor2 = ( LC_suntau(n,q) - LC_suntau(n-1,q) ) / lostau
                  L_factor1 = L_factor1 + Attenuations(n-1) * LC_lostrans_dn(v,n,q)
!mick fix 9/19/2017 - replaced following line
                  !L_factor2 = L_factor2 - factor2 *  L_lostau/lostau
                  L_factor2 = L_factor2 - (factor2 + one) * L_lostau/lostau 
                  LC_multiplier(n,q) = ( L_factor1 - multiplier(n)*L_factor2 ) / factor2
                enddo
              endif
            endif

!  End whole layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  RobFix 8/22/16. Plane/Parallel or Regular-PS, Partial-layer output
!  ------------------------------------------------------------------

      if ( do_RegPSorPP .and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero

 !  Sources, general case

            if ( Mu1(v) .gt. zero ) then
              lostau = kn * path_dn
              if ( lostau .lt. cutoff ) lostrans_dn_p(ut) = exp( - lostau )
              factor1  = Attenuations(np-1)*lostrans_dn_p(ut) - Attenuations_p(ut)
              factor2 = ((suntau_p(ut) - suntau(np-1))/lostau) - one
              multiplier_p(ut) = factor1 / factor2
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_lostau = L_extinction(np,q) * path_dn
                  LC_lostrans_dn_p(ut,q) = - L_lostau * lostrans_dn_p(ut)
                  L_factor1 = LC_Attenuations(np-1,q)*lostrans_dn_p(ut) - LC_Attenuations_p(ut,q)
                  L_factor2 = (LC_suntau_p(ut,q) - LC_suntau(np-1,q))/lostau
                  L_factor1 = L_factor1 + Attenuations(np-1)*LC_lostrans_dn_p(ut,q)
                  L_factor2 = L_factor2 - (factor2 + one) * L_lostau/lostau 
                  L_lostau = L_extinction(np,q) * path_dn
                  LC_multiplier_p(ut,q) = ( L_factor1 - multiplier_p(ut)*L_factor2 ) / factor2
                enddo
              endif
            endif

!  End partial layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/22/16 streamlined code using distances
!      Quadratures from Bottom of the layer

!  1/31/21. Version 2.8.3. LOSTRANS_DN and LC_LosTRANS_DN given geometry index, as now outputs

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_dn = LosW_paths(n,v)
            lostau = kn * path_dn ; if( lostau.lt.cutoff ) lostrans_dn(v,n) = exp ( - lostau )
            if ( do_columnwfs ) then
              do q = 1, Qnums(n)
                L_lostau        = L_extinction(n,q) * path_dn
                LC_lostrans_dn(v,n,q) = - L_lostau * lostrans_dn(v,n)
              enddo
            endif
            sum = zero ; L_sum = zero 
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * wfine(n,j,v)
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_tran = - dj * L_extinction(n,q)
                  L_func = ( LC_attenuationsfine(n,j,q) + L_tran * attenuationsfine(n,j) ) * tran
                  L_sum(q) = L_sum(q) + L_func * wfine(n,j,v)
                enddo
              endif
            enddo
            multiplier(n) = sum * kn
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                LC_multiplier(n,q) = kn * L_sum(q) + sum * L_extinction(n,q)
              enddo
            endif
          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/22/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.cutoff ) lostrans_dn_p(ut) = exp ( - lostau )
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                L_lostau        = L_extinction(np,q) * path_dn
                LC_lostrans_dn_p(ut,q) = - L_lostau * lostrans_dn_p(ut)
              enddo
            endif
            sum = zero ; L_sum = zero 
            do j = 1, nfinedivs_p(ut,v)
              dj = path_dn - xfine_p(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine_p(ut,j) * tran * wfine_p(ut,j,v)
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_tran = - dj * L_extinction(np,q)
                  L_func = ( LC_attenuationsfine_p(ut,j,q) + L_tran * attenuationsfine_p(ut,j) ) * tran
                  L_sum(q) = L_sum(q) + L_func * wfine_p(ut,j,v)
                enddo
              endif
            enddo
            multiplier_p(ut) = sum * kn
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                LC_multiplier_p(ut,q) = kn * L_sum(q) + sum * L_extinction(np,q)
              enddo
            endif
          endif
        enddo
      endif

!  Layer sources
!  -------------

!  General case, Whole layers
!mick fix 9/19/2017 - defined "sources_dn" for all stokes vector elements

      do n = 1, nlayers
        if ( layermask_dn(n) .and.do_sources_dn(n,v)  ) then
          if ( do_sunlight ) then
             shelp(1:ns) = exactscat_dn(n,1:ns,1) * fvec(1)
          else
            do o1 = 1, nstokes
              shelp(o1) = dot_product(exactscat_dn(n,o1,1:ns),fvec(1:ns))
            enddo
          endif
          !sources_dn(n,o1) = shelp(o1) * multiplier(n)
          sources_dn(n,1:ns) = shelp(1:ns) * multiplier(n)
          if ( do_columnwfs ) then
            do q = 1, n_columnwfs
              do o1 = 1, nstokes
                LC_sources_dn(n,o1,q) = shelp(o1) * LC_multiplier(n,q)
                if ( do_sunlight ) then
                  L_Shelp = L_exactscat_dn(n,o1,1,q)*fvec(1)
                else
                  L_Shelp = dot_product(L_exactscat_dn(n,o1,1:ns,q),fvec(1:ns))
                endif
                LC_sources_dn(n,o1,q) =  LC_sources_dn(n,o1,q) + L_Shelp * multiplier(n)
              enddo
            enddo
          endif
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
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                do o1 = 1, nstokes
                  LC_sources_dn_p(ut,o1,q) = shelp(o1) * LC_multiplier_p(ut,q)
                  if ( do_sunlight ) then
                    L_Shelp = L_exactscat_dn(np,o1,1,q)*fvec(1)
                  else
                    L_Shelp = dot_product(L_exactscat_dn(np,o1,1:ns,q),fvec(1:ns))
                  endif
                  LC_sources_dn_p(ut,o1,q) = LC_sources_dn_p(ut,o1,q) + L_Shelp * multiplier_p(ut)
                enddo
              enddo
            endif
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,:) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main Stokes-vector loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. Rob Fix Partials 8/22/16.

!  1/31/21. Version 2.8.3. LOSTRANS_DN given geometry index, as now output

      NC = 0
      CUMSOURCE_DN(NC,:) = zero
      NSTART = 1 ; NUT_PREV = NSTART - 1
      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            do o1 = 1, nstokes
               CUMSOURCE_DN(NC,O1)  = lostrans_dn(v,n) * CUMSOURCE_DN(NC-1,O1) + SOURCES_DN(N,O1)
            enddo
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           STOKES_DN(UTA,1:NS,V) = FLUX * ( CUMSOURCE_DN(NC,1:NS) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT,1:NS) )
         ELSE
           STOKES_DN(UTA,1:NS,V) = FLUX * CUMSOURCE_DN(NC,1:NS)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1 ; NUT_PREV = NUT
      ENDDO

!  Column Wfs (atmospheric term)

!  1/31/21. Version 2.8.3. LOSTRANS_DN and LC_LosTRANS_DN given geometry index, as now outputs

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = 1 ; NUT_PREV = NSTART - 1
         DO UTA = 1, N_USER_LEVELS
            NUT = USER_LEVELS(UTA)
            DO N = NSTART, NUT
               NC = N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) = LC_SOURCES_DN(N,1:ns,Q)         + &
                                LC_lostrans_dn(v,n,q) * CUMSOURCE_DN(NC-1,1:ns) + &
                                   lostrans_dn(v,n)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO

            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_columnwfs
                Term1(1:NS) = CUMSOURCE_DN(NC,1:ns) * LC_LOSTRANS_DN_P(UT,q) + &
                                 L_CUMSOURCE(1:ns,Q) * LOSTRANS_DN_p(UT)
                LC_JACOBIANS_DN(UTA,1:ns,V,Q) = FLUX * ( TERM1(1:NS) + LC_SOURCES_DN_p(UT,1:ns,Q) )
              enddo
            ELSE
              do q = 1, n_columnwfs
                LC_JACOBIANS_DN(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1 ;  NUT_PREV = NUT
         ENDDO
!stop 'at FO col stopping point'
      endif

!  End geometry Loop

   enddo

!  Finish

   return
end subroutine SSV_Integral_ILCS_DN

!

subroutine SSV_Integral_ILCS_UPDN &
   ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input,                                    & ! Inputs (Dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs, max_sleavewfs,                                   & ! Inputs (Dimensioning)
     do_upwelling, do_dnwelling, do_sunlight, do_deltam_scaling, do_fmatrix, do_Partials,            & ! Inputs (Flags - General)
     do_lambertian, do_surface_leaving, do_water_leaving, do_columnwfs, do_surfacewfs, do_sleavewfs, & ! Inputs (Flags - Surf/Lin)
     do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,     & ! Inputs (Flags - Geometry)
     n_reflecwfs, n_sleavewfs, n_surfacewfs, n_columnwfs, Lvarymoms,                                 & ! Inputs (Control, Lin)
     nstokes, ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,                & ! Inputs (Control, Output)
     npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,                    & ! Inputs (Control, Partial)
     flux, fvec, extinction, deltaus, omega, truncfac, Greekmat, fmatrix_up, fmatrix_dn,             & ! Inputs (Flux/Optical)
     reflec, slterm, LS_reflec, LSSL_slterm,                                                         & ! Inputs (Surface/Sleave)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat, L_fmatrix_up, L_fmatrix_dn,           & ! Inputs (Optical - Lin)
     Mu0, Mu1, GenSpher_up, GenSpher_dn, Rotations_up, Rotations_dn, LosW_paths, LosP_paths,         & ! Inputs (Geometry)
     xfine_up, wfine_up, sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,               & ! Inputs (Geometry)
     xfine_dn, wfine_dn, sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,               & ! Inputs (Geometry)
     xfine_up_p, wfine_up_p, sunpaths_up_p, ntraverse_up_p, sunpathsfine_up_p, ntraversefine_up_p,   & ! Inputs (Geometry)
     xfine_dn_p, wfine_dn_p, sunpaths_dn_p, ntraverse_dn_p, sunpathsfine_dn_p, ntraversefine_dn_p,   & ! Inputs (Geometry)
     Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db, lostrans_up, LC_Lostrans_up, & ! Output (upwelling)
     Stokes_dn, LC_Jacobians_dn, cumtrans, LC_cumtrans, lostrans_dn, LC_Lostrans_dn  )                 ! Output (cumtrans & dnwell.)

!  Stand-alone routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of Stokes-vectors and Column/Surface Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using F Matrices directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/22/16
!    - Partial-layer output introduced

!  1/31/21. Version 2.8.3. Introduce lostrans_up/dn, LC_Lostrans_up/dn outputs

   Implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions. Max_sleavewfs added, 8/2/16

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials
   integer, Intent(in) :: maxfine

   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels

   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs
   INTEGER, Intent(in) :: max_sleavewfs

!  flags.
!  Version 1.5: --> F-matrix flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/22/16

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING
   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING
   LOGICAL, Intent(in) :: DO_FMATRIX

   LOGICAL, Intent(in) :: DO_LAMBERTIAN
   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   LOGICAL, Intent(in) :: DO_WATER_LEAVING      ! 4/9/19 added

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Jacobian Flags. do_sleavewfs added 8/2/16

   LOGICAL, Intent(in) :: do_columnwfs
   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_sleavewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT, NSTOKES
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/22/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  Jacobian control. Reflec and sleave numbers added, 8/2/16
!    Note that n_surfacewfs = n_reflecwfs + n_sleavewfs

   INTEGER, Intent(in) :: n_reflecwfs
   INTEGER, Intent(in) :: n_sleavewfs
   INTEGER, Intent(in) :: n_surfacewfs
   INTEGER, Intent(in) :: n_columnwfs
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX, fvec(4)

!  Atmosphere. Fmatrix input added 7/7/16

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )
   REAL(ffp), Intent(in) :: FMATRIX_UP  ( MAXLAYERS, MAXGEOMS, 6 )
   REAL(ffp), Intent(in) :: FMATRIX_DN  ( MAXLAYERS, MAXGEOMS, 6 )

!  Surface reflectivity (Could be the albedo) + linearizations
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( 4, 4, MAXGEOMS )
   real(ffp), Intent(in) :: SLTERM ( 4,    MAXGEOMS )
   real(ffp), Intent(in) :: LS_REFLEC   ( 4,4,MAXGEOMS, max_surfacewfs )
   real(ffp), Intent(in) :: LSSL_SLTERM ( 4,  MAXGEOMS, max_sleavewfs  )

!  Linearized optical inputs. Fmatrix input added 7/7/16

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )
   REAL(ffp), Intent(in) :: L_FMATRIX_UP  ( MAXLAYERS, MAXGEOMS, 6, max_atmoswfs )
   REAL(ffp), Intent(in) :: L_FMATRIX_DN  ( MAXLAYERS, MAXGEOMS, 6, max_atmoswfs )

!  Geometrical inputs
!  ------------------

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

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
!  Los paths added, 8/22/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/22/16.

   real(ffp), Intent(in)  :: xfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: xfine_dn   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn   (maxlayers,maxfine,maxgeoms)


   real(ffp), Intent(in)  :: xfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: xfine_dn_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/22/16.

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

   REAL(ffp), Intent(Out)  :: stokes_up     ( max_user_levels, 4, maxgeoms )
   REAL(ffp), Intent(Out)  :: stokes_db     ( max_user_levels, 4, maxgeoms )
   real(ffp), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, 4, maxgeoms, max_atmoswfs )
   real(ffp), Intent(Out)  :: LC_Jacobians_db  ( max_user_levels, 4, maxgeoms, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, 4, maxgeoms, max_surfacewfs )

   REAL(ffp), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, maxgeoms )
   real(ffp), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, 4, maxgeoms, max_atmoswfs )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS    ( max_user_levels, maxgeoms )
   real(ffp), Intent(out)  :: LC_CUMTRANS ( max_user_levels, maxgeoms, max_atmoswfs )
   
!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP, LC_LosTRANS_UP to the output lists from SSV_Integral_ILCS_UP
!    ==> Add LOSTRANS_DN, LC_LosTRANS_DN to the output lists from SSV_Integral_ILCS_UP

   real(ffp), Intent(out)  :: lostrans_up    (maxgeoms,maxlayers)
   real(ffp), Intent(out)  :: LC_Lostrans_up (maxgeoms,maxlayers,max_atmoswfs)
   real(ffp), Intent(out)  :: lostrans_dn    (maxgeoms,maxlayers)
   real(ffp), Intent(out)  :: LC_Lostrans_dn (maxgeoms,maxlayers,max_atmoswfs)

!  upwelling
!  4/9/19. Additional output for the sleave correction

   if ( do_upwelling  ) then
      call SSV_Integral_ILCS_UP &
      ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input,                     & ! Inputs (dimensioning)
        max_user_levels, max_atmoswfs, max_surfacewfs, max_sleavewfs,                    & ! Inputs (dimensioning)
        do_sunlight, do_deltam_scaling, do_fmatrix, do_lambertian, do_surface_leaving,   & ! Inputs (Flags - General/Surface)
        do_water_leaving, do_Partials, do_PlanPar, do_enhanced_ps, flux, fvec,           & ! Inputs (Flags - Geometry)
        do_sources_up, do_sources_up_p,do_columnwfs, do_surfacewfs, do_sleavewfs,        & ! Inputs (Flags - Linearz)
        n_reflecwfs, n_sleavewfs, n_surfacewfs, n_columnwfs, Lvarymoms,                  & ! Inputs (Control, Jacobian)
        nstokes, ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels, & ! Inputs (control, output)
        npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,     & ! Inputs (control-partial)
        extinction, deltaus, omega, truncfac, greekmat, fmatrix_up, reflec, slterm,      & ! Inputs (Flux/Optical)
        L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat, L_fmatrix_up,          & ! Inputs (Linearized)
        LS_reflec, LSSL_slterm, Mu0, Mu1, GenSpher_up, Rotations_up, LosW_paths, LosP_paths,          & ! Inputs (Geometry)
        xfine_up, wfine_up, sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,             & ! Inputs (Geometry)
        xfine_up_p, wfine_up_p, sunpaths_up_p, ntraverse_up_p, sunpathsfine_up_p, ntraversefine_up_p, & ! Inputs (Geometry)
        Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db,                      & ! Output
        cumtrans, lostrans_up, LC_cumtrans, LC_Lostrans_up )                                            ! Output
   endif

!  Downwelling

   if ( do_dnwelling  ) then
      call SSV_Integral_ILCS_DN &
      ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels, max_atmoswfs, & ! Inputs (dimension)
        do_sunlight, do_deltam_scaling, do_fmatrix, do_Partials, do_PlanPar, do_enhanced_ps,        & ! Inputs (Flags/flux)
        flux, fvec, do_sources_dn, do_sources_dn_p, do_columnwfs, n_columnwfs, Lvarymoms,           & ! Inputs (control, Jacobian )
        nstokes, ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (control)
        npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,                & ! Inputs (control-partial)
        extinction, deltaus, omega, truncfac, greekmat, fmatrix_dn,                                 & ! Inputs (Optical)
        L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat, L_fmatrix_dn,                     & ! Inputs (Optical - Lin)
        Mu1, GenSpher_dn, Rotations_dn, LosW_paths, LosP_paths,                                       & ! Inputs (Geometry)
        xfine_dn, wfine_dn, sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,             & ! Inputs (Geometry)
        xfine_dn_p, wfine_dn_p, sunpaths_dn_p, ntraverse_dn_p, sunpathsfine_dn_p, ntraversefine_dn_p, & ! Inputs (Geometry)
        Stokes_dn, LC_Jacobians_dn, lostrans_dn, LC_Lostrans_dn )                                       ! Output
   endif

!  Finish

   return
end subroutine SSV_Integral_ILCS_UPDN

!  End module

end module FO_VectorSS_RTCalcs_ILCS_m

