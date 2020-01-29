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

! ###########################################################
! #                                                         #
! #                 FIRST-ORDER MODEL                       #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)        #
! #                                                         #
! #  This Version :   1.4 F90                               #
! #  Release Date :   August 2013                           #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3a, 29 October  2012, Obsgeom Multi-geom.   #
! #   Version 1.3b, 24 January  2013, BRDF/SL Supplements   #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_ScalarSS_RTCalcs_ILPS_Optimized_m

!  optimized for upwelling single-level output
!    -- R. Spurr, 6/25/18

!  HERE are the Full NOTES
!  =======================

!  For a given wavelength, this routine will calculate upwelling and downwelling
!  First Order Intensities(I), and any number of LPS Jacobians (profile/surface)

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.
!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 02 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2 , 01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3 , 19 December 2012, Extension to Multiple Geometries, LCS separation
!    Version 4,  31 July     2013, Lattice Multi-geometry

!  For Solar sources, the subroutines are
!       SS_Integral_ILPS_UP   (Upwelling only)
!       SS_Integral_ILPS_DN   (Downwelling only)
!       SS_Integral_ILPS_UPDN (Upwelling and Downwelling)

!  This is the version with Direct phase function inputs (No moments)

!  All subroutines public

public

contains

subroutine SS_Integral_ILPS_UP_Optimized &
   ( maxgeoms, maxlayers, maxfine, max_atmoswfs, max_surfacewfs,                     & ! Inputs (dimensioning)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, do_sleave,                     & ! Inputs (Flags - General)
     do_profilewfs, do_reflecwfs, do_sleavewfs,                                         & ! Inputs (control Jacobian )
     Lvaryflags, Lvarynums, n_reflecwfs, n_sleavewfs,                                   & ! Inputs (control Jacobian )
     ngeoms, nlayers, nfinedivs, ACLEVEL,                                               & ! Inputs (control - output)
     reflec, SLterm, extinction, deltaus, exactscat_up, flux,                           & ! Inputs (Optical - Regular)
     LS_reflec, LSSL_SLterm, L_extinction, L_deltaus, L_exactscat_up,                   & ! Inputs (Optical - Linearized)
     Mu0, Mu1, NCrit, xfine, wfine, csqfine, cotfine,                                   & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                  & ! Inputs (Geometry)
     intensity_up, intensity_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db )      ! Output

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to multiple geometries, 19 December 2012

!  Optimized for  AVIRIS use, R. Spurr 6 February 2018
!     -- Mainly, remove user-level control, only compute to given input level 

!  Optimized for AVIRIS use, R. Spurr 6/25/18
!    - re-introduced surface-leaving option and SLTERM value.
!    - Added control for surface-leaving weighting functions.
!    - renamed "Reflec_WFs" and "Sleave_WFs", to distinguish surface weighting functions

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfine
   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs

!  flags

   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  surface-leaving flag, 6/25/18 for AVIRIS.

   logical, Intent(in) :: DO_SLEAVE

!  Jacobian Flags. Renamed surface flag, added sleave flag (6/25/18)

   LOGICAL, Intent(in) :: do_reflecwfs
   LOGICAL, Intent(in) :: do_sleavewfs
!   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_profilewfs

!  Numbers

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, ACLEVEL

!  Jacobian control. Distinguish relfec and sleave numbers. (6/25/18)
!    n_surfacewfs = SUM OF n_reflecwfs + n_sleavewfs. [Not declared here]

   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   INTEGER, Intent(in) :: n_reflecwfs
   INTEGER, Intent(in) :: n_sleavewfs
!   INTEGER, Intent(in) :: n_surfacewfs

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_UP( MAXLAYERS, MAXGEOMS )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  surface-leaving term, 6/25/18 for AVIRIS.

   real(fpk), Intent(in) :: SLTERM(MAXGEOMS)

!  Linearized optical inputs. surface-leaving term, 6/25/18 for AVIRIS.

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_EXACTSCAT_UP( MAXLAYERS, MAXGEOMS, max_atmoswfs )
   real(fpk), Intent(in) :: LS_REFLEC     ( MAXGEOMS, max_surfacewfs)
   real(fpk), Intent(in) :: LSSL_SLTERM   ( MAXGEOMS, max_surfacewfs)

!  Geometrical inputs
!  ------------------

!       Ray constant, Cotangents
!       Mu0 = cos(theta_boa), required for the surface term (both regular and enhanced)
!       Mu1 = cos(alpha_boa), required for the Regular PS only
!       solar paths, Legendre Polynomials

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  solar paths 

   integer  , Intent(in)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraverse_fine(maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfine,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( maxgeoms )
   real(fpk), Intent(Out)  :: LP_Jacobians_up  ( maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LP_Jacobians_db  ( maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LS_Jacobians_db  ( maxgeoms, max_surfacewfs )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuations      (0:maxlayers)
   real(fpk)  :: LP_attenuations   (0:maxlayers,maxlayers,max_atmoswfs)

!  Solutions for Enhaced-PS

   real(fpk)  :: Solutions_fine    (maxlayers,maxfine)
   real(fpk)  :: LP_solutions_fine (maxlayers,maxfine,maxlayers,max_atmoswfs)

!  Source function integration results

   real(fpk) :: sources_up         ( maxlayers )
   real(fpk) :: LP_sources_up      ( maxlayers,maxlayers,max_atmoswfs )

   real(fpk) :: lostrans_up        ( maxlayers )
   real(fpk) :: L_lostrans_up      ( maxlayers,max_atmoswfs )

   real(fpk) :: cumsource_db      ( 0:maxlayers )
   real(fpk) :: cumsource_up      ( 0:maxlayers )
   real(fpk) :: L_cumsource       ( max_atmoswfs )
   real(fpk) :: LS_cumsource      ( max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, k, j, q, q1, v, nstart, nc, na1, Qnums(maxlayers), n_surfacewfs
   logical    :: Qvary(maxlayers)

   real(fpk)  :: argum(maxfine), tran(maxfine), func(maxfine)
   real(fpk)  :: cons, sum, tran_1, kn, ke, factor1, factor2, m4, term1
   real(fpk)  :: L_sum, L_tran, L_func, L_factor1, L_factor2
   real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau
   real(fpk)  :: L_multiplier, LP_suntau(0:maxlayers,maxlayers,max_atmoswfs), L_lostau(max_atmoswfs)
   real(fpk)  :: attenuations_fine, L_attenuations_fine, sumd, L_sumd, L_exactscat

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   INTENSITY_UP    = zero
   LP_JACOBIANS_UP = zero

   INTENSITY_DB    = zero
   LP_JACOBIANS_DB = zero
   LS_JACOBIANS_DB = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)
   na1 = aclevel + 1

!  number of surface Jacobians

   n_surfacewfs = n_reflecwfs
   if ( do_sleave ) n_surfacewfs = n_surfacewfs + n_sleavewfs

!  Linearization bookkeeping

   Qvary = .false. ; QNums = 0
   if ( do_profilewfs ) then
      Qvary(1:nlayers) = Lvaryflags(1:nlayers)
      QNums(1:nlayers) = Lvarynums (1:nlayers)
   endif

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero local sources

      lostrans_up   = zero  ; sources_up    = zero ; cumsource_up = zero
      L_lostrans_up = zero  ; LP_sources_up = zero

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO   ; LP_Attenuations = ZERO
      Suntau       = ZERO   ; LP_suntau       = ZERO
      Solutions_fine = zero ; LP_Solutions_fine = zero
      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.ntraverse(n,v) ) then
                  do q = 1, Qnums(k)
                     LP_suntau(n,k,q) = L_extinction(k,q) * sunpaths(n,k,v)
                     LP_Attenuations(n,k,q) = - Attenuations(n) * LP_suntau(n,k,q)
                  enddo
               endif
            enddo
         endif
      enddo

!  Enhanced-spherical, fine-layer attenuations
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = na1, nstart
            do j = 1, nfinedivs(n,v)
               sumd = ZERO
               do k = 1, ntraverse_fine(n,j,v)
                  sumd = sumd + extinction(k) * sunpaths_fine(n,k,j,v)
               enddo
               if (sumd .lt. cutoff ) Attenuations_fine = exp( - sumd )
               Solutions_fine(n,j) = exactscat_up(n,v) * Attenuations_fine
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k) .and. k.le.ntraverse_fine(n,j,v) ) then
                        do q = 1, Qnums(k)
                            L_exactscat = zero ; if ( k.eq.n) L_exactscat = L_exactscat_up(n,v,q)
                            L_sumd = L_extinction(k,q) * sunpaths_fine(n,k,j,v)
                            L_Attenuations_fine = - Attenuations_fine * L_sumd 
                            LP_Solutions_fine(n,j,k,q) = L_exactscat         *   Attenuations_fine   + &
                                                           exactscat_up(n,v) * L_Attenuations_fine
                        enddo
                     endif
                  enddo
               endif
            enddo
         enddo
      endif

!  Plane-Parallel or Regular-PS (Average secant formulation): Layer integrated Solar sources
!  =========================================================================================

      if ( do_RegPSorPP ) then
         do n = nlayers, na1, -1

!  LOS transmittance (not for the Horizontal case)

            if ( Mu1(v) .gt. zero ) then
               lostau = deltaus(n)  / Mu1(v)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( do_profilewfs ) then
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_lostau(q)        = L_deltaus(n,q) / Mu1(v)
                        L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                     enddo
                  endif
               endif
            endif

!  Sources, general case

            if ( n.le.nstart  ) then
              if ( Mu1(v) .gt. zero ) then
                factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                factor2 = one + (suntau(n) - suntau(n-1))/lostau
                multiplier = factor1 / factor2
                sources_up(n) = exactscat_up(n,v) * multiplier
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = LP_Attenuations(n-1,k,q) - LP_Attenuations(n,k,q)*lostrans_up(n)
                           L_factor2 = (LP_suntau(n,k,q) - LP_suntau(n-1,k,q))/lostau
                           if ( k.eq.n) then
                              L_factor1 = L_factor1 - Attenuations(n)*L_lostrans_up(n,q)
                              L_factor2 = L_factor2 - (factor2 - one)*L_lostau(q)/lostau 
                              L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                              LP_sources_up(n,k,q) = L_exactscat_up(n,v,q) *   multiplier + &
                                                      exactscat_up(n,v)   * L_multiplier
                           else
                              L_multiplier = ( L_factor1 - multiplier * L_factor2 ) / factor2
                              LP_sources_up(n,k,q) = exactscat_up(n,v)   * L_multiplier
                           endif
                        enddo
                     endif
                  enddo
                endif
              endif
            endif

!  Sources, special case (horizonal view)

            if ( n.le.nstart  ) then
              if ( Mu1(v) .eq. zero ) then
                factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                sources_up(n) = exactscat_up(n,v) * factor1
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = LP_Attenuations(n-1,k,q)
                           if ( k.eq.n) then
                              LP_sources_up(n,k,q) = L_exactscat_up(n,v,q) *   factor1  + &
                                                       exactscat_up(n,v)   * L_factor1
                           else
                              LP_sources_up(n,k,q) = exactscat_up(n,v)   * L_factor1
                           endif
                        enddo
                     endif
                  enddo
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS: special case (nadir viewing). Layer integrated Solar sources
!  =========================================================================

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, na1, -1

!  LOS transmittance

            kn     = extinction(n)
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( do_profilewfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_deltaus(n,q)
                     L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                  enddo
               endif
            endif

!  Sources

            if ( n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = xfine(n,j,v)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = solutions_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j,v)
               enddo
               sources_up(n) = sum * kn
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = LP_solutions_fine(n,j,k,q) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_up(n,k,q)  = L_sum * kn + L_extinction(N,q) * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_func = LP_solutions_fine(n,j,k,q)  * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_up(n,k,q)  = L_sum * kn
                           endif
                        enddo
                     endif
                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS: General case. Layer integrated Solar sources
!  =========================================================

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, na1, -1

!  LOS transmittance

            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            kn = extinction(n) ;  ke = raycon(v) * kn ; cons = raycon(v) * ( cot_2 - cot_1 )
            tran_1 = kn * cons
            if ( tran_1 .lt. cutoff ) lostrans_up(n) = exp ( - tran_1 )
            if ( do_profilewfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_extinction(n,q) * cons
                     L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                  enddo
               endif
            endif

!  Sources

            if ( n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = Raycon(v) * ( cot_2 - cotfine(n,j,v) )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = solutions_fine(n,j) * csqfine(n,j,v) * tran(j)
                  sum      = sum + func(j) * wfine(n,j,v)
               enddo
               sources_up(n) = sum * ke 
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = LP_solutions_fine(n,j,k,q) * csqfine(n,j,v) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_up(n,k,q)  = L_sum * ke + L_extinction(N,q) * Raycon(v) * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_func = LP_solutions_fine(n,j,k,q) * csqfine(n,j,v) * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_up(n,k,q)  = L_sum * ke
                           endif
                        enddo
                     endif
                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Source function integration
!  ===========================

!  Cumulative source terms : Loop over layers working upwards from NSTART to level NA1

!  Reflectance contribution to surface source

      NC = 0 ; M4 = 4.0_fpk * Mu0(v)
      CUMSOURCE_UP(NC) = zero
      CUMSOURCE_DB(NC) = M4 * REFLEC(v) * attenuations(nlayers)

!  Add surface-leaving term (if flagged) to surface source. New, 6/25/18.

      IF ( Do_Sleave ) THEN
         CUMSOURCE_DB(NC) = CUMSOURCE_DB(NC) + SLTERM(v)
      ENDIF

!  Cumulative source terms : Loop over layers working upwards from NSTART to NA1

      DO N = NLAYERS, NA1, -1
         NC = NLAYERS + 1 - N
         CUMSOURCE_DB(NC) = LOSTRANS_UP(N) * CUMSOURCE_DB(NC-1)
         CUMSOURCE_UP(NC) = SOURCES_UP(N) + LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1)
      ENDDO
      INTENSITY_UP(V) = FLUX * CUMSOURCE_UP(NC)
      INTENSITY_DB(V) = FLUX * CUMSOURCE_DB(NC)

!  Reflectance contribution to linearized surface source. Renamed and re-coded, 6/25/18.

      if ( do_reflecwfs ) then
         Term1 = M4  * attenuations(nlayers)
         do q = 1, n_reflecwfs
            LS_cumsource(q) = term1 * LS_reflec(v,q)
         enddo
      endif

!  Add surface-leaving contribution to linearized surface source, New, 6/25/18.

      if ( do_sleave.and. do_sleavewfs ) then
         do q = 1, n_sleavewfs
            q1  = q + n_reflecwfs
            LS_cumsource(q1) = lssl_slterm(v,q)
         enddo
      endif

!  Propagation of surface+sleave weighting functions. Recoded and expanded, 6/25/18.

      if ( do_reflecwfs .or. ( do_sleave.and.do_sleavewfs) ) then
         DO N = NLAYERS, NA1, -1
            NC = NLAYERS + 1 - N
            do q = 1, n_surfacewfs
               LS_cumsource(q) = LOSTRANS_UP(N) * LS_CUMSOURCE(Q)
            enddo
         ENDDO
         do q = 1, n_surfacewfs
            LS_JACOBIANS_DB(V,Q) = FLUX * LS_CUMSOURCE(Q)
         enddo
      endif

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = zero
               DO N = NLAYERS, NA1, -1
                  NC = NLAYERS + 1 - N
                  if ( k.eq.n ) then
                     do q = 1, Qnums(k)
                        L_cumsource(q) = LP_SOURCES_UP(N,K,Q) + &
                             L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1) + LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                     enddo
                  else
                     do q = 1, Qnums(k)
                        L_cumsource(q) = LP_SOURCES_UP(N,K,Q) + LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                     enddo
                  endif
               ENDDO
               do q = 1, Qnums(k)
                  LP_JACOBIANS_UP(V,K,Q) = FLUX * L_CUMSOURCE(Q)
               enddo
            endif
         enddo
      endif

!  Profile Wfs (direct beam term)

      if ( do_profilewfs ) then
         term1 = M4 * reflec(v)
         do k = 1, nlayers
            if ( Qvary(k) ) then
               do q = 1, Qnums(k)
                  L_CUMSOURCE(q) = Term1 * LP_attenuations(nlayers,k,q)
               enddo
               DO N = NSTART, NA1, -1
                  NC = NLAYERS + 1 - N
                  if ( k.eq.n ) then
                     do q = 1, Qnums(k)
                        L_cumsource(q) =  L_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1) + &
                                            LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                     enddo
                  else
                     do q = 1, Qnums(k)
                        L_cumsource(q) = LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                     enddo
                  endif
               ENDDO
               do q = 1, Qnums(k)
                  LP_JACOBIANS_DB(V,K,Q) = FLUX * L_CUMSOURCE(Q)
               enddo
            endif
         enddo
      endif

!  End Geometry Loop

   enddo

!  Finish

   return
end subroutine SS_Integral_ILPS_UP_Optimized

!  End module

end module FO_ScalarSS_RTCalcs_ILPS_Optimized_m


