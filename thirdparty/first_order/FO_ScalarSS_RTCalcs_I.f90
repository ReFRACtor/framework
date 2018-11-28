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
! #     FIRST-ORDER SCALAR MODEL (EXACT SINGLE SCATTERING)  #
! #                                                         #
! #  This Version :   1.3 F90                               #
! #  Release Date :   March 2013                            #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3,  29 October  2012, Multiple geometries   #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_ScalarSS_RTCalcs_I_m

!  For a given wavelength, this routine will calculate First-Order upwelling+downwelling Intensities (I):

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry

!  For Solar sources, the subroutines are
!       SS_Integral_I_UP   (Upwelling only)
!       SS_Integral_I_DN   (Downwelling only)
!       SS_Integral_I_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine SS_Integral_I_UP &
   ( maxgeoms, maxlayers, maxfine, max_user_levels, & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,  & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     reflec, extinction, deltaus, exactscat_up, flux,        & ! Inputs (Optical)
     Mu0, Mu1, NCrit, xfine, wfine, csqfine, cotfine,        & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpathsfine, ntraversefine,         & ! Inputs (Geometry)
     intensity_up, intensity_db, cumsource_up )                                ! Outputs

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)

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
   integer, Intent(in) :: max_user_levels

!  flags

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_UP( MAXLAYERS, MAXGEOMS )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

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

   real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels,maxgeoms )
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

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, k, L, v
   logical    :: layermask_up(maxlayers)
   real(fpk)  :: sumd, help, sum, tran_1, tran, func, kn, ke, xjkn
   real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau
   real(fpk)  :: cumsource_db

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   CUMSOURCE_UP = zero ; INTENSITY_UP = zero ; INTENSITY_DB = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

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

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_up(n) ) then
               do j = 1, nfinedivs(n,v)
                  sumd = ZERO
                  do k = 1, ntraversefine(n,j,v)
                     sumd = sumd + extinction(k) * sunpathsfine(n,k,j,v)
                  enddo
                  if (sumd .lt. cutoff ) Attenuationsfine(n,j) = exp( - sumd )
                  Solutionsfine(n,j) = exactscat_up(n,v) * Attenuationsfine(n,j)
               enddo
            endif
         enddo
      endif

!  Plane/Parallel or Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

      if ( do_RegPSorPP ) then
         factor1 = zero ; factor2 = zero
         if ( Mu1(v) .eq. zero ) then
            do n = 1, nlayers
               Solutions(n) = exactscat_up(n,v) * Attenuations(n-1)
               if ( layermask_up(n) ) then
                  if ( attenuations(n-1).ne.zero ) then
                     factor1(n) = Attenuations(n)/Attenuations(n-1)
                     nstart = n
                  endif
               endif
            enddo
         else
            do n = 1, nlayers
               lostau = deltaus(n) / Mu1(v)
               Solutions(n) = exactscat_up(n,v) * Attenuations(n-1)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( layermask_up(n) ) then
                  if ( attenuations(n-1).ne.zero ) then
                     factor1(n) = Attenuations(n)/Attenuations(n-1)
                     factor2(n) = (suntau(n) - suntau(n-1))/lostau
                     nstart = n
                  endif
               endif
            enddo
         endif
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS Layer integrated source terms

      if ( do_RegPSorPP ) then
         do n = nlayers, 1, -1
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
         do n = nlayers, 1, -1
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
         do n = nlayers, 1, -1
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
      CUMSOURCE_DB     = 4.0_fpk * Mu0(v) * REFLEC(v) * attenuations(nlayers)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DB       = LOSTRANS_UP(N) * CUMSOURCE_DB
            CUMSOURCE_UP(NC,V) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,V) + SOURCES_UP(N)
         ENDDO
         INTENSITY_UP(UTA,V) = FLUX * CUMSOURCE_UP(NC,V)
         INTENSITY_DB(UTA,V) = FLUX * CUMSOURCE_DB
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  End GEometry Loop

   enddo

!  Finish

   return
end subroutine SS_Integral_I_UP

!

subroutine SS_Integral_I_DN &
   ( maxgeoms, maxlayers, maxfine, max_user_levels, & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,    & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     extinction, deltaus, exactscat_dn, flux,                & ! Inputs (Optical)
     Mu1, NCrit, RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
     Raycon, radii, cota, sunpaths, ntraverse, sunpathsfine, ntraversefine,    & ! Inputs (Geometry)
     intensity_dn, cumsource_dn )                                                ! Outputs

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)

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
   integer, Intent(in) :: max_user_levels

!  flags

   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_REGULAR_PS
   logical, Intent(in) ::  DO_ENHANCED_PS
   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  doNadir(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  NGEOMS

   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_DN( MAXLAYERS,MAXGEOMS )

!  Solar Flux

   real(fpk), Intent(in) ::  FLUX

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu1(maxgeoms)

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

   real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: cumsource_dn     ( 0:maxlayers,maxgeoms )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuations      (0:maxlayers)
   real(fpk)  :: attenuationsfine (maxlayers,maxfine)

!  Solutions

   real(fpk)  :: Solutionsfine (maxlayers,maxfine)
   real(fpk)  :: Solutions (0:maxlayers)

!  Source function integration results

   real(fpk)  :: sources_dn       ( maxlayers )
   real(fpk)  :: lostrans_dn      ( maxlayers )

!  Average secant type arrays (REGULAR-PS only)

   real(fpk)  :: factor1       ( maxlayers )
   real(fpk)  :: factor2       ( maxlayers )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, k, L, v
   logical    :: layermask_dn(maxlayers)
   real(fpk)  :: sumd, help, sum, tran_1, tran, func, kn, ke, xjkn, trand
   real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau, rdiff, cot_c

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zeron the output

      CUMSOURCE_DN = zero ; INTENSITY_DN = zero
      lostrans_dn  = zero ; sources_dn   = zero

!  Regular_PS or plane-parallel flag

      do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

      NUT = USER_LEVELS(N_USER_LEVELS) + 1
      IF ( NUT > NLAYERS ) NUT = NLAYERS
      LAYERMASK_DN = .false.
      LAYERMASK_DN(1:NUT) = .true.

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero the local sources

      lostrans_dn = zero  ; sources_dn  = zero

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO ; Attenuationsfine = ZERO
      Solutions    = zero ; Solutionsfine    = zero
      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_dn(n) ) then
               do j = 1, nfinedivs(n,v)
                  sumd = ZERO
                  do k = 1, ntraversefine(n,j,v)
                     sumd = sumd + extinction(k) * sunpathsfine(n,k,j,v)
                  enddo
                  if (sumd .lt. cutoff ) Attenuationsfine(n,j) = exp( - sumd )
                  Solutionsfine(n,j) = exactscat_dn(n,v) * Attenuationsfine(n,j)
               enddo
            endif
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

      if ( do_RegPSorPP ) then
        factor1 = zero ; factor2 = zero ; Solutions = zero
        do n = nlayers, 1, -1 
          if ( Mu1(v) .gt. zero ) then
            lostau = deltaus(n) / Mu1(v)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
          endif
          if ( layermask_dn(n) .and. n.le.nstart  ) then
            if ( Mu1(v) .gt. zero ) then
               factor1(n) = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
               factor2(n) = ((suntau(n) - suntau(n-1))/lostau) - one
               multiplier = factor1(n) / factor2(n)
               sources_dn(n) = exactscat_dn(n,v) * multiplier
            else if ( Mu1(v) .eq. zero ) then
               factor1(n) = Attenuations(n)
               sources_dn(n) = exactscat_dn(n,v) * factor1(n)
            endif
          endif
        enddo
      endif

!  Enhanced PS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_dn(n)  = exp ( - deltaus(n))
            rdiff = radii(n-1) - radii(n) ; if ( n.eq.NCrit(v)) rdiff = radii(n-1) - RadCrit(v)
            trand = one ; if ( n.eq.NCrit(v)) trand = exp ( -kn * (RadCrit(v) -radii(n) ) )
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  xjkn = ( rdiff - xfine(n,j,v) ) * kn
                  func = solutionsfine(n,j) * exp ( - xjkn )
                  sum = sum + func * wfine(n,j,v)
               enddo
               sources_dn(n) = sum * kn
!  @@ Robfix, add following line
               if ( n.eq.NCrit(v) ) sources_dn(n) = sources_dn(n) * trand
!  @@ End Robfix add line
            endif
         enddo
      endif

!  Enhanced PS: General case

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, 1, -1
            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            cot_c = cot_1  ; if ( n.eq.NCrit(v) ) cot_c = CotCrit(v)
            kn = extinction(n) ;  ke = raycon(v) * kn
            trand = one  ; if ( n.eq.NCrit(v) ) trand = exp ( - ke * ( CotCrit(v) - cot_1 ) )
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_dn(n) = tran_1
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  tran = exp ( - ke * ( cotfine(n,j,v) - cot_c ) )   !  Down
                  func = solutionsfine(n,j) * csqfine(n,j,v) * tran
                  sum  = sum + func * wfine(n,j,v)
               enddo
               sources_dn(n) = sum * ke * trand
            endif        
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,V) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!     Check for dndating the recursion

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC,v) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1,v)
         ENDDO
         INTENSITY_DN(UTA,V) = FLUX * CUMSOURCE_DN(NC,v)
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  End geometry loop

   enddo

!  Finish

   return
end subroutine SS_Integral_I_DN

!

subroutine SS_Integral_I_UPDN   &
   ( maxgeoms, maxlayers, maxfine, max_user_levels, & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_deltam_scaling, & ! Inputs (Flags)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,     & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     reflec, extinction, deltaus, exactscat_up, exactscat_dn, flux,          & ! Inputs (Optical)
     Mu0, Mu1, NCrit, RadCrit, CotCrit,                      & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                    & ! Inputs (Geometry)
     sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,           & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,           & ! Inputs (Geometry)
     intensity_up, intensity_db, cumsource_up, intensity_dn, cumsource_dn )    ! Outputs

!  Stand-alone routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)

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
   integer, Intent(in) :: max_user_levels

!  flags

   logical, Intent(in) ::  DO_UPWELLING
   logical, Intent(in) ::  DO_DNWELLING
   logical, Intent(in) ::  DO_DELTAM_SCALING

   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_REGULAR_PS
   logical, Intent(in) ::  DO_ENHANCED_PS
   logical, Intent(in) ::  DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  NGEOMS

   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_UP( MAXLAYERS,MAXGEOMS )
   real(fpk), Intent(in) :: EXACTSCAT_DN( MAXLAYERS,MAXGEOMS )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu0 = cos(theta_boa), required for the surface term (both regular and enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms), RadCrit(maxgeoms), CotCrit(maxgeoms)

!  solar paths 

   integer  , Intent(in)  :: ntraverse_up     (0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_up      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_up (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: sunpathsfine_up  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_dn     (0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_dn      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_dn (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: sunpathsfine_dn  (maxlayers,maxlayers,maxfine,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfine,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: cumsource_up     ( 0:maxlayers,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: cumsource_dn     ( 0:maxlayers,maxgeoms )

!  Upwelling
!  ---------

   if ( do_upwelling ) then
       call SS_Integral_I_UP &
   ( maxgeoms, maxlayers, maxfine, max_user_levels,          & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,     & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     reflec, extinction, deltaus, exactscat_up, flux,        & ! Inputs (Optical)
     Mu0, Mu1, NCrit, xfine, wfine, csqfine, cotfine, Raycon, cota, & ! Inputs (Geometry)
     sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,              & ! Inputs (Geometry)
     intensity_up, intensity_db, cumsource_up )                                   ! Outputs
   endif

   if ( do_dnwelling ) then
       call SS_Integral_I_DN &
   ( maxgeoms, maxlayers, maxfine, max_user_levels,          & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,    & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     extinction, deltaus, exactscat_dn, flux,                & ! Inputs (Optical)
     Mu1, NCrit, RadCrit, CotCrit,                           & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                      & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,             & ! Inputs (Geometry)
     intensity_dn, cumsource_dn )                                                ! Outputs
   endif

!  Finish

   return
end subroutine SS_Integral_I_UPDN

!  End module

end module FO_ScalarSS_RTCalcs_I_m

