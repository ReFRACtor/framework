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

module FO_ScalarSS_RTCalcs_I_Optimized_m

!  optimized for upwelling single-level output
!    -- R. Spurr, 6/25/18

!  HERE are the Full NOTES
!  =======================

!  For a given wavelength, this routine will calculate First-Order upwelling+downwelling Intensities (I):

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry

!  For Solar sources, the subroutines are
!       SS_Integral_I_UP   (Upwelling only)
!       SS_Integral_I_DN   (Downwelling only)
!       SS_Integral_I_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine SS_Integral_I_UP_Optimized &
   ( maxgeoms, maxlayers, maxfine,                                   & ! Inputs (dimensioning)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, do_sleave,  & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, ACLEVEL,                            & ! Inputs (control output)
     reflec, slterm, extinction, deltaus, exactscat_up, flux,        & ! Inputs (Optical)
     Mu0, Mu1, NCrit, xfine, wfine, csqfine, cotfine,                & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpathsfine, ntraversefine, & ! Inputs (Geometry)
     intensity_up, intensity_db, cumsource_up )                           ! Outputs

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)

!  Optimized for AVIRIS use, R. Spurr 6 February 2018
!    --  Mainly, user-level inputs are gone, plus only do calculation to a given level.

!  Optimized for AVIRIS use, R. Spurr 6/25/18
!    - re-introduced surface-leaving option and SLTERM value.

!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfine

!  flags

   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  surface-leaving flag, 6/25/18 for AVIRIS.

   logical, Intent(in) :: DO_SLEAVE

!  Numbers

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, ACLEVEL

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_UP( MAXLAYERS,MAXGEOMS )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  surface-leaving term, 6/25/18 for AVIRIS.

   real(fpk), Intent(in) :: SLTERM(MAXGEOMS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  solar paths 

   integer  , Intent(in)  :: ntraverse     (0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( maxgeoms )
   real(fpk), Intent(Out)  :: cumsource_up     ( 0:maxlayers,maxgeoms )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuations      (0:maxlayers)
   real(fpk)  :: attenuationsfine (maxlayers,maxfine)

!  Solutions

   real(fpk)  :: Solutionsfine (maxlayers,maxfine)
   real(fpk)  :: Solutions (0:maxlayers)

!  Source function integration results

   real(fpk)  :: sources_up       ( maxlayers )
   real(fpk)  :: lostrans_up      ( maxlayers )

!  Average secant type arrays (REGULAR-PS or Plane-Parallel)

   real(fpk)  :: factor1       ( maxlayers )
   real(fpk)  :: factor2       ( maxlayers )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, nstart, nc, j, v, nt, ntj, na1
   real(fpk)  :: sumd, sum, tran_1, tran, func, kn, ke, xjkn
   real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau
   real(fpk)  :: cumsource_db

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   CUMSOURCE_UP = zero ; INTENSITY_UP = zero ; INTENSITY_DB = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  bookkeeping for the single-level output

   na1 = aclevel + 1

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero the local sources

      lostrans_up = zero  ; sources_up   = zero

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO ; Attenuationsfine = ZERO
      Solutions    = zero ; Solutionsfine    = zero
      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = aclevel, nlayers
         nt = ntraverse(n,v) ;  suntau(n) = DOT_PRODUCT(extinction(1:nt),sunpaths(n,1:nt,v))
         If( suntau(n) .lt. cutoff ) Attenuations(n) = exp( - suntau(n) )
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = na1, nlayers
            do j = 1, nfinedivs(n,v)
               sumd = ZERO
               ntj = ntraversefine(n,j,v) ;  sumd = DOT_PRODUCT(extinction(1:nt),sunpathsfine(n,1:nt,j,v))
               if (sumd .lt. cutoff ) Attenuationsfine(n,j) = exp( - sumd )
               Solutionsfine(n,j) = exactscat_up(n,v) * Attenuationsfine(n,j)
            enddo
         enddo
      endif

!  Plane/Parallel or Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

      if ( do_RegPSorPP ) then
         factor1(na1:nlayers) = zero
         factor2(na1:nlayers) = zero
         if ( Mu1(v) .eq. zero ) then
            do n = na1, nlayers
               Solutions(n) = exactscat_up(n,v) * Attenuations(n-1)
               if ( attenuations(n-1).ne.zero ) then
                  factor1(n) = Attenuations(n)/Attenuations(n-1) ; nstart = n
               endif
            enddo
         else
            do n = na1, nlayers
               lostau = deltaus(n) / Mu1(v)
               Solutions(n) = exactscat_up(n,v) * Attenuations(n-1)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( attenuations(n-1).ne.zero ) then
                  factor1(n) = Attenuations(n)/Attenuations(n-1)
                  factor2(n) = (suntau(n) - suntau(n-1))/lostau
                  nstart = n
               endif
            enddo
         endif
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS Layer integrated source terms

      if ( do_RegPSorPP ) then
         do n = nlayers, na1, -1
            if ( n.le.nstart ) then
               multiplier = ( one - Factor1(n)*lostrans_up(n) ) / (factor2(n) + one)
            else
               multiplier = zero
            endif  
            sources_up(n) = solutions(n) * multiplier
         enddo
      endif

!  Enhanced PS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, na1, -1
            kn = extinction(n)
            lostrans_up(n)  = exp ( - deltaus(n))
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  xjkn = xfine(n,j,v) * kn
                  func = solutionsfine(n,j) * exp ( - xjkn )
                  sum = sum + func * wfine(n,j,v)
               enddo
               sources_up(n) = sum * kn
            endif
         enddo
      endif

!  Enhanced PS: General case

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, na1, -1
            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            kn = extinction(n) ;  ke = raycon(v) * kn
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_up(n) = tran_1
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  tran = exp ( - ke * ( cot_2 - cotfine(n,j,v) ) )
                  func = solutionsfine(n,j) * csqfine(n,j,v) * tran
                  sum  = sum + func * wfine(n,j,v)
               enddo
               sources_up(n) = sum * ke 
            endif        
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion ( FOr Direct Beam, use PI.mu0.R.Atten )

      NC =  0
      CUMSOURCE_UP(NC,V) = zero
      CUMSOURCE_DB       = 4.0_fpk * Mu0(v) * REFLEC(v) * attenuations(nlayers)

!  surface-leaving term. Added, 6/25/18.

      IF ( DO_Sleave ) THEN
         CUMSOURCE_DB = CUMSOURCE_DB + SLTERM(v)
      ENDIF

!   Cumulative source terms : Loop over layers working upwards from NSTART to NA1

      DO N = NLAYERS, NA1, -1
         NC = NLAYERS + 1 - N
         CUMSOURCE_DB       = LOSTRANS_UP(N) * CUMSOURCE_DB
         CUMSOURCE_UP(NC,V) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,V) + SOURCES_UP(N)
      ENDDO
      INTENSITY_UP(V) = FLUX * CUMSOURCE_UP(NC,V)
      INTENSITY_DB(V) = FLUX * CUMSOURCE_DB

!  End Geometry Loop

   enddo

!  Finish

   return
end subroutine SS_Integral_I_UP_Optimized

!  End module

end module FO_ScalarSS_RTCalcs_I_Optimized_m
