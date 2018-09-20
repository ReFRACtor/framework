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
! #  This Version :   1.4 F90                               #
! #  Release Date :   August 2013                           #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3,  29 October  2012, Observation geometry  #
! #   Version 1.4,  31 July     2013, Lattice geometries    #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_geometry_Routines_m

use FO_geometry_Generic_m

!  Plane-parallel and Regular PS routines
!  --------------------------------------

!    subroutine Obsgeom_PlanPar
!    subroutine Lattice_PlanPar

!    subroutine Obsgeom_RegularPS
!    subroutine Lattice_RegularPS

!  Following routines are for Enhanced PS case
!  -------------------------------------------

!  LOS-Outgoing: Apply equally to Thermal and Solar-scatter Cases

!    subroutine LosOut_EnhancedPS_Initial
!    subroutine LosOut_EnhancedPS_Quadrature
!    subroutine LosOut_EnhancedPS_QUpgrade
!    subroutine LosOut_EnhancedPS_FindCrit

!  Copying for solar

!    subroutine LOSCopy1_EnhancedPS_Lattice

!  Solar-incoming routines (solar scatter only)

!   subroutine SolarIn_EnhancedPS_Obsgeom_FindCrit
!   subroutine SolarIn_EnhancedPS_Lattice_FindCrit
!   subroutine SolarIn_EnhancedPS_Obsgeom_SunPaths
!   subroutine SolarIn_EnhancedPS_Lattice_SunPaths

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains


subroutine ObsGeom_PlanPar &
          ( maxgeoms, maxlayers, dtr, vsign, do_Chapman,                & ! Inputs
            ngeoms, nlayers, obsgeom_boa, heights,                      & ! Inputs
            alpha, sunpaths, ntraverse, chapfacs, cosscat )               ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Plane-parallel choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   integer  , intent(In)    :: maxlayers, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign

   logical  , intent(In)    :: do_Chapman
   integer  , intent(In)    :: nlayers, ngeoms

   real(ffp), intent(InOut) :: obsgeom_boa(maxgeoms,3)
   real(ffp), intent(In)    :: heights (0:maxlayers)

!  Los geometry

   real(ffp), intent(Out)  :: alpha         (0:maxlayers,maxgeoms)

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: cosscat(maxgeoms)

!  Local

   logical        :: Do_OverheadSun
   integer        :: n, v
   real(ffp)      :: alpha_boa_R, theta_boa_R
   real(ffp)      :: salpha_boa, calpha_boa
   real(ffp)      :: stheta_boa, ctheta_boa, utheta_boa, cphi_boa, diffhts(maxlayers)
   real(ffp)      :: term1(maxgeoms), term2(maxgeoms)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Initialise output

   alpha = zero
   ntraverse = 0
   sunpaths  = zero ; chapfacs  = zero
   cosscat   = zero ; term1 = zero ; term2 = zero

!  start geometry loop

   do v = 1, ngeoms

!  BOA angles

      alpha_boa_R    = obsgeom_boa(v,2) * DTR
      if ( obsgeom_boa(v,2).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         calpha_boa     = cos(alpha_boa_R)
         salpha_boa     = sin(alpha_boa_R)
      endif

      theta_boa_R    = obsgeom_boa(v,1) * DTR
      stheta_boa     = sin(theta_boa_R)
      ctheta_boa     = cos(theta_boa_R)

      cphi_boa       = cos(obsgeom_boa(v,3) * dtr)

!  Nominal traverse paths for Full illumination. Difference heights

      do n = 1, nlayers
         diffhts(n) = heights(n-1) - heights(n)
         ntraverse(n,v) = n
      enddo

!  Overhead Sun

      Do_OverheadSun = (obsgeom_boa(v,1).eq.zero) 

!  Set Alpha, scattering angle

      alpha(1:nlayers,v) = alpha_boa_R
      if ( Do_OverheadSun ) then
         term1(v) = zero
         term2(v) = calpha_boa
         cosscat(v) = - vsign * term2(v) ; if (term2(v).eq.zero) cosscat(v) = term2(v)
      else
         term1(v) = salpha_boa * stheta_boa * cphi_boa
         term2(v) = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2(v) + term1(v) 
      endif

!  Sunpath/Chapman factor calculations

      utheta_boa     = one / ctheta_boa
      do n = 1, nlayers
         sunpaths(n,1:n,v) = diffhts(1:n) * utheta_boa
         if ( do_Chapman ) chapfacs(n,1:n,v) = utheta_boa
      enddo

!  End geometry routine

   enddo

!  Finish

   return
end subroutine Obsgeom_PlanPar

!

subroutine Lattice_PlanPar &
          ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, dtr, vsign, & ! Input
            do_Chapman, nszas, nvzas, nazms, nv_offset, na_offset,      & ! Input
            nlayers, alpha_boa, theta_boa, phi_boa, heights,            & ! Inputs
            alpha, sunpaths, ntraverse, chapfacs, cosscat )               ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Plane-parallel choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   integer  , intent(In)    :: maxlayers, maxszas, maxvzas, maxazms, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign
   logical  , intent(In)    :: do_Chapman

   integer  , intent(In)    :: nlayers, nszas, nvzas, nazms
   integer  , intent(in)    :: nv_offset(maxszas)
   integer  , intent(in)    :: na_offset(maxszas,maxvzas)

   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)
   real(ffp), intent(In)    :: heights (0:maxlayers)

!  Los geometry outputs

   real(ffp), intent(Out)  :: alpha  (0:maxlayers,maxgeoms)

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: cosscat    (maxgeoms)

!  Local
!  -----

   logical       :: Do_OverheadSun(maxszas)
   integer       :: n, k, v, ns, nv, na, os1, os2, nviews
   real(ffp)     :: utheta_boa, alpha_boa_R, theta_boa_R, term1, term2, fac
   real(ffp)     :: salpha_boa(maxvzas), calpha_boa(maxvzas)
   real(ffp)     :: stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms), diffhts(maxlayers)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Initialise output

   alpha = zero
   ntraverse = 0 ; sunpaths  = zero ; chapfacs  = zero ; cosscat   = zero

!  start geometry loops

   do nv = 1, nvzas
      alpha_boa_R    = alpha_boa(nv) * DTR
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)    = zero
         salpha_boa(nv)    = one
      else
         calpha_boa(nv)    = cos(alpha_boa_R)
         salpha_boa(nv)    = sin(alpha_boa_R)
      endif
      do ns = 1, nszas
         do na = 1, nazms
            v = na_offset(ns,nv) + na
            alpha(1:nlayers,v) = alpha_boa_R
         enddo
      enddo
   enddo

   do na = 1, nazms
      cphi_boa(na)       = cos(phi_boa(na) * dtr)
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = (theta_boa(ns).eq.zero) 
      theta_boa_R    = theta_boa(ns) * DTR
      stheta_boa(ns) = sin(theta_boa_R)
      ctheta_boa(ns) = cos(theta_boa_R)
   enddo

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo

!  Nviews

   nviews = nvzas * nazms

!  Main loop

   do ns = 1, nszas

!  Set scattering angle

      if ( Do_OverheadSun(ns) ) then
         do nv = 1, nvzas
            term1 = zero ;  term2 = -vsign * calpha_boa(nv)
            if ( term2.ne.zero ) then
               do na = 1, nazms
                  v = na_offset(ns,nv) + na
                  cosscat(v) = term2
               enddo
            endif
         enddo
      else
         do nv = 1, nvzas
            term2 = - vsign * calpha_boa(nv) * ctheta_boa(ns)
            term1 = salpha_boa(nv) * stheta_boa(ns)
            do na = 1, nazms
               v = na_offset(ns,nv) + na
               cosscat(v) = term2 + term1 * cphi_boa(na)
            enddo
         enddo
      endif

!  Sunpath/Chapman factor calculations

      utheta_boa     = one / ctheta_boa(ns)
      os1 = nv_offset(ns) + 1
      os2 = nv_offset(ns) + nviews
      do n = 1, nlayers
         ntraverse(n,os1:os2) = n
         do k = 1, n
            fac = diffhts(k) * utheta_boa
            sunpaths(n,k,os1:os2) = fac
            if ( do_Chapman ) chapfacs(n,k,os1:os2) = utheta_boa
         enddo
      enddo

!  End sza loop

   enddo

!  Finish

   return
end subroutine Lattice_PlanPar

!

subroutine Obsgeom_RegularPS &
      ( maxgeoms, maxlayers, dtr, vsign, do_Chapman,                  & ! Inputs
        ngeoms, nlayers, obsgeom_boa, heights, eradius,               & ! Inputs
        Raycon, radii, alpha, sunpaths, ntraverse, chapfacs, cosscat )  ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   integer  , intent(In)    :: maxlayers, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign
   logical  , intent(In)    :: do_Chapman

   integer  , intent(In)    :: nlayers, ngeoms
   real(ffp), intent(InOut) :: obsgeom_boa(maxgeoms,3)
   real(ffp), intent(In)    :: eradius, heights (0:maxlayers)

!  Los geometry outputs

   real(ffp), intent(Out)  :: alpha  (0:maxlayers, maxgeoms)
   real(ffp), intent(Out)  :: radii  (0:maxlayers)
   real(ffp), intent(Out)  :: Raycon (maxgeoms)

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: cosscat(maxgeoms)

!  Local

   logical        :: Do_OverheadSun
   integer        :: n, k, v
   real(ffp)      :: alpha_boa_R, theta_boa_R
   real(ffp)      :: salpha_boa, calpha_boa
   real(ffp)      :: stheta_boa, ctheta_boa, cphi_boa, sunpaths_local(maxlayers), diffhts(maxlayers)
   real(ffp)      :: term1(maxgeoms), term2(maxgeoms)

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp

!  Initialise output

   radii = zero ; alpha = zero ; Raycon = zero
   ntraverse = 0
   sunpaths  = zero ; chapfacs  = zero
   cosscat   = zero ; term1 = zero ; term2 = zero

!  Radii

   do n = 0, nlayers
      radii(n) = eradius + heights(n)
   enddo

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo

!  Start geometry loop

   do v = 1, ngeoms

!  BOA angles

      alpha_boa_R    = obsgeom_boa(v,2) * DTR
      if ( obsgeom_boa(v,2).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         calpha_boa     = cos(alpha_boa_R)
         salpha_boa     = sin(alpha_boa_R)
      endif

      theta_boa_R    = obsgeom_boa(v,1) * DTR
      if ( obsgeom_boa(v,1).eq.90.0_ffp ) then
         ctheta_boa     = zero
         stheta_boa     = one
      else
         stheta_boa     = sin(theta_boa_R)
         ctheta_boa     = cos(theta_boa_R)
      endif

      cphi_boa       = cos(obsgeom_boa(v,3) * dtr)

!  Nominal traverse paths for Full illumination

      do n = 1, nlayers
         ntraverse(n,v) = n
      enddo

!  Overhead Sun

      Do_OverheadSun = (obsgeom_boa(v,1).eq.zero) 

!  Set Alpha, ray constant, scattering angle

      alpha(1:nlayers,v) = alpha_boa_R
      Raycon           = salpha_boa * radii(nlayers)
      if ( Do_OverheadSun ) then
         term1(v) = zero
         term2(v) = calpha_boa
         cosscat(v) = - vsign * term2(v) ; if (term2(v).eq.zero) cosscat(v) = term2(v)
      else
         term1(v) = salpha_boa * stheta_boa * cphi_boa
         term2(v) = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2(v) + term1(v)
      endif

!  Sunpath/Chapman factor calculations

!mick fix 4/12/12 - adjusted dimension of array "sunpaths" assignments to be
!                   compatible with subroutine "FindSunPaths_D" computations
      !sunpaths(0,1:maxlayers) = zero

      sunpaths(0,1:nlayers,v) = zero
      do n = 1, nlayers
         call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,&
           theta_boa_R,stheta_boa,N,sunpaths_local)
         !sunpaths(n,1:maxlayers) = sunpaths_local(1:maxlayers)
         sunpaths(n,1:n,v) = sunpaths_local(1:n)
         if ( do_Chapman ) then
            do k = 1, n
               chapfacs(n,k,v) = sunpaths(n,k,v)/diffhts(k)
            enddo
         endif
      enddo

!  End geometry loop

   enddo

!  Finish

   return
end subroutine Obsgeom_RegularPS


subroutine Lattice_RegularPS &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, dtr, vsign,     & ! Input
         do_Chapman, nszas, nvzas, nazms, nv_offset, na_offset, nlayers, & ! I  nput
         alpha_boa, theta_boa, phi_boa, heights, eradius,                & ! Inputs
         Raycon, radii, alpha, sunpaths, ntraverse, chapfacs, cosscat )    ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

   integer  , intent(In)    :: maxlayers, maxszas, maxvzas, maxazms, maxgeoms
   integer  , intent(In)    :: nlayers, nszas, nvzas, nazms
   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)
   logical  , intent(In)    :: do_Chapman
   real(ffp), intent(In)    :: dtr, heights (0:maxlayers), vsign, eradius

   integer , intent(in)     :: nv_offset(maxszas)
   integer , intent(in)     :: na_offset(maxszas,maxvzas)

!  Los geometry

   real(ffp), intent(Out)  :: alpha  (0:maxlayers, maxgeoms)
   real(ffp), intent(Out)  :: radii  (0:maxlayers)
   real(ffp), intent(Out)  :: Raycon (maxgeoms)

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: cosscat    (maxgeoms)

!  Local

   real(ffp)     :: sunpaths_local(maxlayers)
   logical       :: Do_OverheadSun(maxszas)
   integer       :: n, k, v, ns, nv, na, os1, os2, nviews
   real(ffp)     :: alpha_boa_R, Raycon_R, theta_boa_R(maxszas), term1, term2
   real(ffp)     :: salpha_boa(maxvzas), calpha_boa(maxvzas)
   real(ffp)     :: stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms), diffhts(maxlayers)

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp

!  Initialise output

   radii = zero ; alpha = zero ; Raycon = zero
   ntraverse = 0 ; sunpaths  = zero ; chapfacs  = zero ; cosscat   = zero

!  Radii

   do n = 0, nlayers
      radii(n) = eradius + heights(n)
   enddo

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo

!  start geometry loops

   do nv = 1, nvzas
      alpha_boa_R    = alpha_boa(nv) * DTR
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)    = zero
         salpha_boa(nv)    = one
      else
         calpha_boa(nv)    = cos(alpha_boa_R)
         salpha_boa(nv)    = sin(alpha_boa_R)
      endif
      Raycon_R = salpha_boa(nv) * radii(nlayers)
      do ns = 1, nszas
         do na = 1, nazms
            v = na_offset(ns,nv) + na
            alpha(1:nlayers,v) = alpha_boa_R
            Raycon(v)          = Raycon_R
         enddo
      enddo
   enddo

   do na = 1, nazms
      cphi_boa(na)       = cos(phi_boa(na) * dtr)
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = (theta_boa(ns).eq.zero) 
      theta_boa_R(ns)    = theta_boa(ns) * DTR
      stheta_boa(ns) = sin(theta_boa_R(ns))
      ctheta_boa(ns) = cos(theta_boa_R(ns))
   enddo

!  Nviews

   nviews = nvzas * nazms

!  Main loop

   do ns = 1, nszas

!  Set scattering angle cosine

      if ( Do_OverheadSun(ns) ) then
         do nv = 1, nvzas
            term1 = zero ;  term2 = -vsign * calpha_boa(nv)
            if ( term2.ne.zero ) then
               do na = 1, nazms
                  v = na_offset(ns,nv) + na
                  cosscat(v) = term2
               enddo
            endif
         enddo
      else
         do nv = 1, nvzas
            term2 = - vsign * calpha_boa(nv) * ctheta_boa(ns)
            term1 = salpha_boa(nv) * stheta_boa(ns)
            do na = 1, nazms
               v = na_offset(ns,nv) + na
               cosscat(v) = term2 + term1 * cphi_boa(na)
            enddo
         enddo
      endif

!  Sunpath/Chapman factor calculations

      os1 = nv_offset(ns) + 1
      os2 = nv_offset(ns) + nviews
      do n = 1, nlayers
         ntraverse(n,os1:os2) = n
         call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii(n),Radii,&
           theta_boa_R(ns),stheta_boa(ns),N,sunpaths_local)
         do k = 1, n
            sunpaths(n,k,os1:os2) = sunpaths_local(k)
         enddo
         if ( do_Chapman ) then
            do k = 1, n
               chapfacs(n,k,os1:os2) = sunpaths(n,k,os1:os2)/diffhts(k)
            enddo
         endif
      enddo

!  End sza loop

   enddo


!  Finish

   return
end subroutine Lattice_RegularPS

!

subroutine LosOut_EnhancedPS_Initial  &
          ( maxgeoms, maxlayers, dtr, ngeoms, nlayers, & ! Input
            heights, eradius, alpha_boa,               & ! Input
            doNadir, radii, Raycon, Lospaths,          & ! Output
            alpha, sina, cosa, cota )                    ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!  Extension to Multiple Geometries, 29 October 2012
!  Extension to Lattice  Geometries, 31 July    2013 
!   Equally valid for the Lattice case, if we understand ngeoms = nvzas in this case.

!  This routine: Initial LOS path setup

!    starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                        - height grid, earth radius

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer  , intent(in)   :: maxlayers, maxgeoms
   real(ffp), intent(in)   :: dtr

!  number of geometries

   integer  , intent(in)   :: ngeoms

!  Layer control

   integer  , intent(in)   :: nlayers
   real(ffp), intent(in)   :: eradius, heights (0:maxlayers)

!  input angles

   real(ffp), intent(in)   :: alpha_boa(maxgeoms)

!  Flag for the Nadir case

   logical  , intent(out)  :: doNadir(maxgeoms)
  
!  Alphas, Radii, Ray constant, Lospaths

   real(ffp), intent(out)  :: radii    (0:maxlayers)
   real(ffp), intent(out)  :: Raycon   (maxgeoms)
   real(ffp), intent(out)  :: Lospaths (maxlayers,maxgeoms)

   real(ffp), intent(out)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cosa     (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: sina     (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cota     (0:maxlayers,maxgeoms)

!  Local
!  -----

   integer      :: n, n1, v
   real(ffp)    :: salpha_boa, difh, alpha_boa_R
   real(ffp)    :: calpha, calpha1

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Zero output

   cota = zero ; cosa = zero ; sina = zero ; alpha = zero
   donadir = .false. ; Raycon = zero ; Lospaths = zero

!  Radii

   do n = 0, nlayers
     radii(n) = eradius + heights(n)
   enddo

!  START LOOP
!  ==========

   do v = 1, ngeoms

!  Special case. Direct nadir viewing. Compute everything and Exit.

      if ( alpha_boa(v).eq.zero ) then
         doNadir(v) = .true.
         do n = nlayers,1,-1
            difh = radii(n-1) - radii(n) ; Lospaths(n,v) = difh
         enddo
         go to 67
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  start at BOA

      alpha_boa_R    = alpha_boa(v) * DTR
      if ( alpha_boa(v) .eq. 90.0_ffp ) then
         salpha_boa     = one
         calpha1        = zero
      else
         salpha_boa     = sin(alpha_boa_R)
         calpha1        = cos(alpha_boa_R)
      endif

      cosa(nlayers,v)  = calpha1
      sina(nlayers,v)  = salpha_boa
      cota(nlayers,v)  = calpha1 / salpha_boa
      alpha(nlayers,v) = alpha_boa_R

!  Ray constant

      Raycon(v) = salpha_boa * radii(nlayers)

!  Whole layer values

      do n = nlayers - 1, 0, -1
         n1 = n + 1
         sina(n,v) = Raycon(v) / radii(n) ; alpha(n,v) = asin(sina(n,v))
         calpha  = cos(alpha(n,v)) ; cosa(n,v) = calpha 
         cota(n,v) = cosa(n,v) / sina(n,v)
         Lospaths(n1,v) = radii(n)*calpha - radii(n1)*calpha1
         calpha1 = calpha
      enddo

!  End loop

67    continue
   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_Initial

!

subroutine LosOut_EnhancedPS_Quadrature  &
       ( maxgeoms, maxlayers, maxfine,     & ! Input 
         ngeoms, nlayers, nfine,           & ! Input
         doNadir, radii, alpha, Raycon,    & ! Input
         nfinedivs, radiifine, alphafine,  & ! Output
         xfine, wfine, csqfine, cotfine )    ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!  Extension to Multiple Geometries, 29 October 2012
!  Extension to Lattice  Geometries, 31 July    2013 
!   Equally valid for the Lattice case, if we understand ngeoms = nvzas in this case.

!    starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                        - height grid, earth radius, Layer control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer, intent(in)       :: maxlayers, maxfine, maxgeoms

!  Layer and geometry numbers control

   integer, intent(in)       :: nlayers, ngeoms, nfine

!  Flag for the Nadir case

   logical  , intent(in)     :: doNadir(maxgeoms)
  
!  Alphas, Radii, Ray constant

   real(ffp), intent(in)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii      (0:maxlayers)
   real(ffp), intent(in)  :: Raycon      (maxgeoms)

!  Outputs
!  =======

!  Finelayer divisions is output here

   integer  , intent(out)  :: nfinedivs (maxlayers,maxgeoms)

!  Fine layering

   real(ffp), intent(out)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Quadratures

   real(ffp), intent(out)  :: xfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: wfine    (maxlayers,maxfine,maxgeoms)

!  Local geoemetry arrays

   real(ffp), intent(out)  :: csqfine  (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: cotfine  (maxlayers,maxfine,maxgeoms)

!  Local
!  -----

   integer            :: n, n1, j, v
   real(ffp)          :: difh, csfine
   real(ffp)          :: tfine(maxfine), afine(maxfine)

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Zero output

   nfinedivs  = 0
   alphafine = zero    ; radiifine = zero
   xfine      = zero    ; wfine      = zero
   cotfine    = zero    ; csqfine    = zero

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )

      if ( doNadir(v) ) then
         do n = nlayers,1,-1
            difh  = radii(n-1) - radii(n)
            nfinedivs(n,v) = nfine
            call gauleg_ng (zero,difh,tfine,afine,nfine,maxfine)
            do j = 1, nfine
               radiifine(n,j,v) = radii(n-1) - tfine(j)
               xfine(n,j,v) = tfine(j)
               wfine(n,j,v) = afine(j)
            enddo
         enddo
         go to 67
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  Whole layer values

      do n = nlayers, 1, -1
         n1 = n - 1
         nfinedivs(n,v) = nfine
         call gauleg_ng (alpha(n1,v),alpha(n,v),tfine,afine,nfine,maxfine)
         do j = 1,  nfine
            csfine = one / sin(tfine(j))
            radiifine(n,j,v) = raycon(v) * csfine
            alphafine(n,j,v) = tfine(j)
            xfine(n,j,v)   = radii(n1) - radiifine(n,j,v)
            wfine(n,j,v)   = afine(j)
            cotfine(n,j,v) = cos(tfine(j)) * csfine
            csqfine(n,j,v) = csfine * csfine
         enddo
      enddo

!  Continuation point

67    continue

!  End geometry loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_Quadrature

!
! @@ REACHED THIS POINT

SUBROUTINE LosOut_EnhancedPS_FindCrit &
       ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs

!  Purpose: Given a list of Maximum extinctions and LOS angles
!     Then find Critical layers (NCrit) and point where LOS attenuation wipe-outs (Acrit) are achieved
!     Then find the LOS angles and Radii (AlphaCrit,RadCrit) for these Critical Points

!  Extension to Multiple Geometries, 29 October 2012
!   Equally valid for the Lattice case, if we understand ngeoms = nvzas in this case.

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer  , intent(in)  :: maxlayers, maxgeoms

!  Layer and geometry numbers control

   integer  , intent(in)  :: nlayers, ngeoms

!  Attenuation and other parameters

   real(ffp), intent(in)  :: Acrit, Cutoff

!  Special case, Nadir viewing

   logical  , intent(in)  :: doNadir(maxgeoms)

!  View angles and Radii at layer boundaries

   real(ffp), intent(in)  :: sina (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: cosa (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii(0:maxlayers)

!  Extinctions

   real(ffp), intent(in)  :: Lospaths(maxlayers,maxgeoms)
   real(ffp), intent(in)  :: extinc(maxlayers)

!  Modified inputs
!  ---------------

!  Number of Fine divisions

   integer, intent(inout) :: nfinedivs(maxlayers,maxgeoms)

!  outputs
!  -------

!  Critical layer, Number of Fine divisions for this layer

   integer  , intent(out)  :: Ncrit(maxgeoms)

!  Critical angle and radius and cotangent

   real(ffp), intent(out)  :: AlphaCrit(maxgeoms)
   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Other variables

   logical    ::  trawl
   integer    ::  N, N1, ntrans, v
   real(ffp)  ::  opdep, cumtrans_0, cumtrans_1,trans,transcrit,dcrit,TanCrit

!  Initialize

   Ncrit     = 0
   AlphaCrit = zero
   RadCrit   = zero
   CotCrit   = zero

   fail    = .false.
   message = ' '

!  Start Geometry loop

   do v = 1, ngeoms

!  Set trawl
!   Tested March 17th

      trawl = .true. ; n = 0 ; cumtrans_0 = one
      do while (trawl .and.n.lt.nlayers) 
         n = n + 1
         opdep = Lospaths(n,v) * extinc(n) ; trans = zero
         if ( opdep .lt. cutoff ) trans = exp ( - opdep )
         cumtrans_1 = cumtrans_0 * trans
         if ( cumtrans_1 .gt. Acrit ) then
            ntrans = int(-log10(trans) + 1)
            nfinedivs(n,v) = max(nfinedivs(n,v),ntrans)
            cumtrans_0 = cumtrans_1
         else
            NCrit(v) = n ; trawl = .false.
            transcrit    = Acrit / cumtrans_0
            ntrans = int(-log10(transcrit) + 1)
            nfinedivs(n,v) = max(nfinedivs(n,v),ntrans)
            dcrit        = - log(transcrit) / extinc(n)
            if ( doNadir(v) ) then
               Radcrit(v)    = radii(n-1) - dcrit      
            else
               n1 = n-1 ; TanCrit = radii(n1)*sina(n1,v)/(radii(n1)*cosa(n1,v)-dcrit)
               Cotcrit(v)    = one / Tancrit
               alphacrit(v)  = atan( TanCrit)
               radcrit(v)    = sina(n,v) * radii(n) / sin(alphacrit(v))    
            endif
         endif
      enddo

!  Zero the rest

      if ( NCrit(v) .ne. 0 ) nfinedivs(NCrit(v)+1:nlayers,v) = 0

!  End geometry loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_FindCrit


subroutine LOSCopy1_EnhancedPS_Lattice &
       ( maxgeoms, maxszas, maxvzas, maxlayers, maxfine,      & ! Input dimensions
         nvzas, nszas, nazms, nlayers, na_offset,             & ! Input   (# angles/layers)
         doNadir_LOS, Raycon_LOS, alpha_LOS, cota_LOS,        & ! Input LOS
         nfinedivs_LOS, radiifine_LOS, alphafine_LOS,         & ! Input LOS
         xfine_LOS, wfine_LOS, csqfine_LOS, cotfine_LOS,      & ! Input LOS
         doNadir, Raycon, alpha, cota,                        & ! Output
         nfinedivs, radiifine, alphafine,                     & ! Output
         xfine, wfine, csqfine, cotfine )                       ! Output

!  Purpose: Given a list LOS-angle quantities, Copy to all geometries

!  Extension to Observational Geometries, 29 October 2012
!  Extension to Lattice       Geometries, 31 July    2013

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensioning

   integer  , intent(in)  :: maxgeoms, maxszas, maxvzas, maxlayers, maxfine

!  Layer and geometry numbers control

   integer  , intent(in)  :: nvzas, nszas, nazms, nlayers,  na_offset (maxszas,maxvzas)

!  LOS-only input 

   logical  , intent(in)  :: doNadir_LOS (maxvzas)
   real(ffp), intent(in)  :: alpha_LOS   (0:maxlayers,maxvzas)
   real(ffp), intent(in)  :: cota_LOS    (0:maxlayers,maxvzas)
   real(ffp), intent(in)  :: Raycon_LOS  (maxvzas)

!  LOS Fine layering, inputs

   integer  , intent(in)  :: nfinedivs_los (maxlayers,maxvzas)
   real(ffp), intent(in)  :: alphafine_los (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: radiifine_los (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: xfine_los     (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: wfine_los     (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: csqfine_los   (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: cotfine_los   (maxlayers,maxfine,maxvzas)

!  Outputs
!  =======

!  All-geometry

   logical  , intent(inout)  :: doNadir (maxgeoms)
   real(ffp), intent(inout)  :: alpha   (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cota    (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: Raycon  (maxgeoms)

!  Fine layering, Quadratures and local geometry, Adjusted values

   integer  , intent(inout)  :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: xfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: csqfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cotfine   (maxlayers,maxfine,maxgeoms)

!  Local variables
!  ---------------

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Other variables

   integer     ::  j, n, nv, ns, na, g

!  Initialize arrays with the LOS values

      do nv = 1, nvzas
         do ns = 1, nszas
            do na = 1, nazms
               g =  na_offset(ns,nv) + na
               alpha(0:nlayers,g) = alpha_LOS(0:nlayers,nv)
               cota(0:nlayers,g)  = cota_LOS (0:nlayers,nv)
               Raycon(g)          = Raycon_LOS (nv)
               doNadir(g)         = doNadir_LOS(nv)
               Nfinedivs(1:nlayers,g) = nfinedivs_LOS(1:nlayers,nv)
            enddo
         enddo
      enddo

!  Fine layer stuff (only if no Criticality)

      do nv = 1, nvzas
        do ns = 1, nszas
          do na = 1, nazms
            g = na_offset(ns,nv) + na
            do n = 1, nlayers
              do j = 1, Nfinedivs(n,g)
                radiifine(n,j,g)  = radiifine_LOS(n,j,nv)
                alphafine(n,j,g)  = alphafine_LOS(n,j,nv)
                xfine(n,j,g)      = xfine_LOS(n,j,nv)
                wfine(n,j,g)      = wfine_LOS(n,j,nv)
                csqfine(n,j,g)    = csqfine_LOS(n,j,nv)
                cotfine(n,j,g)    = cotfine_LOS(n,j,nv)
              enddo
            enddo
          enddo
        enddo
      enddo
! finish

   return

end subroutine LOSCopy1_EnhancedPS_Lattice

!

subroutine SolarIn_EnhancedPS_Obsgeom_FindCrit &
       ( maxgeoms, maxlayers, ngeoms, nlayers, doNadir, dtr, Acrit, & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,           & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit,     & ! Outputs
         fail, message )                                              ! Outputs

!  Purpose: Given a list of Maximum extinctions and solar angles at BOA
!           Then find Critical layers (NCrit) and points where TOA attenuation wipe-outs (Acrit) are achieved
!           Then find the LOS angles and Radii (AlphaCrit,RadCrit) for these Critical Points
!           Nadir case, Alpha = 0.0, find only the radius (RadCrit)

!  Extension to Multiple Geometries, 29 October 2012

!  Find the Critical Radius (or angle) in layer Ncrit_i, Use Bisection 
!      based on the Function F(x) = opdep(x) - Crit_opdep

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer, intent(in) :: maxlayers, maxgeoms

!  Layer and geometry numbers control

   integer, intent(in) :: nlayers, ngeoms

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir(maxgeoms)

!  View angles and Radii at layer boundaries, Ray constant

   real(ffp), intent(in)  :: alpha(0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii(0:maxlayers)
   real(ffp), intent(in)  :: Raycon(maxgeoms)

!  Extinctions

   real(ffp), intent(in)  :: extinc(maxlayers)

!  Solar control and other parameters

   real(ffp), intent(in)  :: Acrit, theta_boa(maxgeoms), dtr, cutoff

!  Modified inputs (outputs)
!  -------------------------

!  Overall control (May be switched off if Critical test is negative for all Geometries)

   logical  , intent(inout)  :: DoCrit

!  Critical layer

   integer  , intent(inout)  :: Ncrit(maxgeoms)

!  Number of Fine divisions. This is updated according to Criticality

   integer  , intent(inout)  :: nfinedivs(maxlayers,maxgeoms)

!  Critical angle and radius and cotangent

   real(ffp), intent(inout)  :: AlphaCrit(maxgeoms)
   real(ffp), intent(inout)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Outputs
!  -------

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  Bisection accuracy (lower for the Nadir case, as using distances)

   real(ffp), parameter  :: BisectionAccuracy_Nadir   = 1.0d-6
   real(ffp), parameter  :: BisectionAccuracy_General = 1.0d-9
   integer  , parameter  :: jmax = 50

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Other variables

   character*2::  c2
   logical    ::  Finding, trawl, do_ZeroSunBOA, doCrit_local(maxgeoms)
   integer    ::  j, n, ntrans, ncrit_i, k, v, cc
   real(ffp)  ::  s0, x1, x2, xmid, rtbis, dx, f, fmid, suncon, radii_x, sx, accuracy
   real(ffp)  ::  dist, theta_0, theta_1, stheta_0, stheta_1, ground_R, sunpaths(maxlayers)
   real(ffp)  ::  opdep, trans, atten_0, atten_1, transcrit, theta_n, stheta_n, theta_boa_R

!  Initialize

   fail    = .false.
   message = ' '

!  Crit Count, initialize

   CC = 0

!  Start Geometry loop

   do v = 1, ngeoms

!  Initial setups

      doCrit_local = .true.
      atten_0 = one ; theta_boa_R = theta_boa(v) * dtr
      if ( doNadir(v) ) then
         s0 = sin(theta_boa_R) ; ground_R = zero ; accuracy = BisectionAccuracy_Nadir
      else
         s0 = zero ; ground_R = alpha(nlayers,v) + theta_boa_R ; accuracy = BisectionAccuracy_General
      endif

!  Trawl through layers until Critical layer is reached. Nfinedivs is updated.
!    Only go down to Initial (LOS-path) Critical layer 
!       Condition on ZerosunBOA changed 27 March 2012

      NCrit_i = 0 ; trawl = .true. ; n = 0
      do while ( trawl .and.( (Ncrit(v).eq.0.and.n.lt.nlayers).or.(NCrit(v).ne.0.and.n.lt.NCrit(v)) ) )
         n = n + 1
         do_ZeroSunBOA = (n.eq.nlayers.and.theta_boa(v).eq.zero).or.(doNadir(v).and.theta_boa(v).eq.zero)
         if ( doNadir(v) ) then
            theta_n = theta_boa_R ; stheta_n = s0
         else
            theta_n = ground_R - alpha(n,v) ; stheta_n = sin(theta_n)
         endif
         call FindSunPaths_D (do_ZeroSunBOA,Maxlayers,radii(n),Radii,&
                              theta_n,stheta_n,n,sunpaths(:))
         opdep = sum(extinc(1:n)*sunpaths(1:n)) ; trans = zero
         atten_1 = zero; if ( opdep .lt. cutoff ) atten_1 = exp ( - opdep )
         if ( atten_1 .gt. Acrit ) then
            trans = atten_1 / atten_0
            ntrans = int(-log10(trans) + 1)
            nfinedivs(n,v) = max(nfinedivs(n,v),ntrans)
            atten_0 = atten_1
         else
            NCrit_i = n ; trawl = .false.
            transcrit    = Acrit / atten_0
            ntrans = int(-log10(transcrit) + 1)
            nfinedivs(n,v) = max(nfinedivs(n,v),ntrans)
         endif
      enddo

!  Nothing to do if No criticality (previous Critical values are unchanged)

      if ( trawl .and. NCrit_i.eq. 0 ) then
         if ( trawl .and. NCrit(v) .eq. 0 ) DoCrit_local = .false.
         go to 67
      endif

!  Bisection: set Highest/Lowest value of Function (layer bottom/top). 

      if ( doNadir(v) ) then
         x1 = zero             ; x2 = radii(NCrit_i-1) - radii(NCrit_i)
      else
         x1 = alpha(NCrit_i-1,v) ; x2 = alpha(NCrit_i,v) 
      endif
      F     = log(Atten_0) - cutoff
      FMID  = opdep        - cutoff

!  Bisection: Check bracketing, if OK, perform bisection

      IF(F*FMID.GE.zero) then
         write(c2,'(I2)')v
         fail = .true. ; message = 'Root must be bracketed for bisection, geometry #'//c2 ; return
      ENDIF
      IF(F.LT.zero)THEN
         RTBIS=X1 ; DX=X2-X1
      ELSE
         RTBIS=X2 ; DX=X1-X2
      ENDIF

!  Bisection: Iterate to find the answer

      Finding = .true. ; J = 0
      DO While (Finding .and. j .lt. JMAX)
         J = J + 1 ; dx = 0.5_ffp * dx ; XMID = RTBIS + DX
         if ( doNadir(v) ) then
            theta_0 = theta_boa_R ; stheta_0 = s0
         radii_x = radii(NCrit_i-1) - xmid 
         else
            theta_0 = ground_R - xmid ; stheta_0 = sin(theta_0)
            radii_x = Raycon(v) / sin(xmid)
         endif
         suncon = radii_x * stheta_0
         stheta_1 = suncon / radii(NCrit_i-1) ;  theta_1 = asin(stheta_1)
         dist = radii(NCrit_i-1) * sin(theta_0-theta_1) / stheta_0
         opdep = dist * extinc(NCrit_i)
         theta_0 = theta_1 ; stheta_1 = stheta_0
         do k = n - 1, 1, -1
            stheta_1 = suncon / radii(k-1) ; theta_1 = asin(stheta_1)
            dist = radii(k-1) * sin(theta_0-theta_1) / stheta_0
            opdep = opdep + dist * extinc(k)
            theta_0 = theta_1 ; stheta_0 = stheta_1
         enddo
         fmid = opdep - cutoff
         IF ( FMID.LE.zero ) RTBIS = XMID
         IF(ABS(DX).LT.Accuracy .OR. FMID.EQ.0.) Finding = .false.
      ENDDO

!  Exception (too many bisections)

      if ( Finding ) Then
         write(c2,'(I2)')v
         fail = .true. ; message = 'Too many Bisections (540); Root not found, geometry #'//c2 ; return
      endif

!  Set final output if successful

      CC = CC + 1
      if ( doNadir(v) ) then
         RadCrit(v)   =  radii(NCrit_i-1) - RTBIS
      else
         AlphaCrit(v) = RTBIS ;  SX = sin(AlphaCrit(v))
         RadCrit(v)   = Raycon(v) / SX
         CotCrit(v)   = sqrt(one-SX*SX)/SX
      endif
      NCrit(v)     = NCrit_i

!  Continuation point

67    continue

!  Zero the rest

      if ( NCrit(v) .ne. 0 ) nfinedivs(NCrit(v)+1:nlayers,v) = 0

!  End geometry loop

   enddo

!  Reset criticality

   doCrit = ( CC.gt.0 )
 
!  Finish

end subroutine SolarIn_EnhancedPS_Obsgeom_FindCrit

!

subroutine SolarIn_EnhancedPS_Lattice_FindCrit &
       ( maxgeoms, maxvzas, maxszas, maxlayers, ACrit, dtr, cutoff, & ! Input   (dimensions/constants)
         nvzas, nszas, nazms, na_offset, nlayers, theta_boa,        & ! Input   (# angles/layers)
         doNadir, alpha, radii, extinc, Raycon, nfinedivs,          & ! Inputs  (LOS and extinction)
         Ncrit_los, AlphaCrit_los, RadCrit_los, CotCrit_los,        & ! Inputs  (LOS criticality)
         doCrit, Ncrit, AlphaCrit, RadCrit, CotCrit,                & ! Outputs (Incoming criticality)
         fail, message )                                              ! Outputs (Exception handling)

!  Purpose: Given a list of Maximum extinctions and solar angles at BOA
!           Then find Critical layers (NCrit) and points where TOA attenuation wipe-outs (Acrit) are achieved
!           Then find the LOS angles and Radii (AlphaCrit,RadCrit) for these Critical Points
!           Nadir case, Alpha = 0.0, find only the radius (RadCrit)

!  Extension to Observational Geometries, 29 October 2012
!  Extension to Lattice       Geometries, 31 July    2013

!  Find the Critical Radius (or angle) in layer Ncrit_i, Use Bisection 
!      based on the Function F(x) = opdep(x) - Crit_opdep

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer  , intent(in)  :: maxgeoms, maxvzas, maxszas, maxlayers 

!  other parameters

   real(ffp), intent(in)  :: Acrit,  dtr, cutoff

!  Layer and geometry numbers control

   integer  , intent(in)  :: nvzas, nszas, nazms, nlayers
   integer  , intent(in)  :: na_offset(maxszas,maxvzas)
   real(ffp), intent(in)  :: theta_boa(maxszas)

!  Special case, Nadir viewing

   logical  , intent(in)  :: doNadir(maxgeoms)

!  View angles and Radii at layer boundaries, Ray constant

   real(ffp), intent(in)  :: alpha(0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii(0:maxlayers)
   real(ffp), intent(in)  :: Raycon(maxgeoms)

!  Extinctions

   real(ffp), intent(in)  :: extinc(maxlayers)

!  Criticality for the outgoing path: Layer #, LOS angle, radius and cotangent

   integer  , intent(in)  :: Ncrit_los     (maxvzas)
   real(ffp), intent(in)  :: AlphaCrit_los (maxvzas)
   real(ffp), intent(in)  :: RadCrit_los   (maxvzas)
   real(ffp), intent(in)  :: CotCrit_los   (maxvzas)

!  Updated Number of Fine divisions. Modified input

   integer  , intent(inout)  :: nfinedivs(maxlayers,maxgeoms)

!  outputs
!  -------

!  Overall control (May be switched off if Critical test is negative for all Geometries)

   logical  , intent(inout)  :: DoCrit

!  Critical layer

   integer  , intent(inout)  :: Ncrit(maxgeoms)

!  Critical angle and radius and cotangent

   real(ffp), intent(inout)  :: AlphaCrit(maxgeoms)
   real(ffp), intent(inout)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Outputs
!  -------

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  Bisection accuracy (lower for the Nadir case, as using distances)

   real(ffp), parameter  :: BisectionAccuracy_Nadir   = 1.0d-6
   real(ffp), parameter  :: BisectionAccuracy_General = 1.0d-9
   integer  , parameter  :: jmax = 50

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Other variables

   character*3 ::  c3
   logical     ::  Finding, trawl, do_ZeroSunBOA, doCrit_local(maxgeoms)
   integer     ::  j, n, ntrans, ncrit_i, k, nv, ns, na, g, g1, cc
   real(ffp)   ::  s0, x1, x2, xmid, rtbis, dx, f, fmid, suncon, radii_x, sx, accuracy
   real(ffp)   ::  dist, theta_0, theta_1, stheta_0, stheta_1, ground_R, sunpaths(maxlayers)
   real(ffp)   ::  opdep, trans, atten_0, atten_1, transcrit, theta_n, stheta_n, theta_boa_R

!  Initialize

   fail    = .false.
   message = ' '

!  Crit Count, initialize

   CC = 0

!  Start LOS, SZA, and AZM geometry loops
!  --------------------------------------

   do nv = 1, nvzas
      do ns = 1, nszas

!  Geometry index for lattice, first AZM value

         g = na_offset(ns,nv) + 1

!  Initialize arrays with the LOS values

         NCrit(g)     = Ncrit_los(nv)
         AlphaCrit(g) = AlphaCrit_los(nv)
         RadCrit(g)   = Radcrit_los(nv)
         CotCrit(g)   = CotCrit_los(nv)

!  Initial setups

         doCrit_local = .true.
         atten_0 = one ; theta_boa_R = theta_boa(g) * dtr
         if ( doNadir(g) ) then
            s0 = sin(theta_boa_R) ; ground_R = zero ; accuracy = BisectionAccuracy_Nadir
         else
            s0 = zero ; ground_R = alpha(nlayers,g) + theta_boa_R ; accuracy = BisectionAccuracy_General
         endif

!  Trawl through layers until Critical layer is reached. Nfinedivs is updated.
!    Only go down to Initial (LOS-path) Critical layer 
!       Condition on ZerosunBOA changed 27 March 2012

         NCrit_i = 0 ; trawl = .true. ; n = 0
         do while ( trawl .and.( (Ncrit(g).eq.0.and.n.lt.nlayers).or.(NCrit(g).ne.0.and.n.lt.NCrit(g)) ) )
            n = n + 1
            do_ZeroSunBOA = (n.eq.nlayers.and.theta_boa(ns).eq.zero).or.(doNadir(g).and.theta_boa(ns).eq.zero)
            if ( doNadir(g) ) then
               theta_n = theta_boa_R ; stheta_n = s0
            else
               theta_n = ground_R - alpha(n,g) ; stheta_n = sin(theta_n)
            endif
            call FindSunPaths_D (do_ZeroSunBOA,Maxlayers,radii(n),Radii,&
                                       theta_n,stheta_n,n,sunpaths(:))
            opdep = sum(extinc(1:n)*sunpaths(1:n)) ; trans = zero
             atten_1 = zero; if ( opdep .lt. cutoff ) atten_1 = exp ( - opdep )
            if ( atten_1 .gt. Acrit ) then
               trans = atten_1 / atten_0
               ntrans = int(-log10(trans) + 1)
               nfinedivs(n,g) = max(nfinedivs(n,g),ntrans)
               atten_0 = atten_1
            else
               NCrit_i = n ; trawl = .false.
               transcrit    = Acrit / atten_0
               ntrans = int(-log10(transcrit) + 1)
               nfinedivs(n,g) = max(nfinedivs(n,g),ntrans)
            endif
         enddo

!  Nothing to do if No criticality (previous Critical values are unchanged)

         if ( trawl .and. NCrit_i.eq. 0 ) then
            if ( trawl .and. NCrit(g) .eq. 0 ) DoCrit_local = .false.
            go to 67
         endif

!  Bisection: set Highest/Lowest value of Function (layer bottom/top). 

         if ( doNadir(g) ) then
            x1 = zero             ; x2 = radii(NCrit_i-1) - radii(NCrit_i)
         else
            x1 = alpha(NCrit_i-1,g) ; x2 = alpha(NCrit_i,g) 
         endif
         F     = log(Atten_0) - cutoff
         FMID  = opdep        - cutoff

!  Bisection: Check bracketing, if OK, perform bisection

         IF(F*FMID.GE.zero) then
            write(c3,'(I3)')g
            fail = .true. ; message = 'Root must be bracketed for bisection, geometry #'//c3 ; return
         ENDIF
         IF(F.LT.zero)THEN
            RTBIS=X1 ; DX=X2-X1
         ELSE
            RTBIS=X2 ; DX=X1-X2
         ENDIF

!  Bisection: Iterate to find the answer

         Finding = .true. ; J = 0
         DO While (Finding .and. j .lt. JMAX)
            J = J + 1 ; dx = 0.5_ffp * dx ; XMID = RTBIS + DX
            if ( doNadir(g) ) then
               theta_0 = theta_boa_R ; stheta_0 = s0
               radii_x = radii(NCrit_i-1) - xmid 
            else
               theta_0 = ground_R - xmid ; stheta_0 = sin(theta_0)
               radii_x = Raycon(g) / sin(xmid)
            endif
            suncon = radii_x * stheta_0
            stheta_1 = suncon / radii(NCrit_i-1) ;  theta_1 = asin(stheta_1)
            dist = radii(NCrit_i-1) * sin(theta_0-theta_1) / stheta_0
            opdep = dist * extinc(NCrit_i)
            theta_0 = theta_1 ; stheta_1 = stheta_0
            do k = n - 1, 1, -1
               stheta_1 = suncon / radii(k-1) ; theta_1 = asin(stheta_1)
               dist = radii(k-1) * sin(theta_0-theta_1) / stheta_0
               opdep = opdep + dist * extinc(k)
               theta_0 = theta_1 ; stheta_0 = stheta_1
            enddo
            fmid = opdep - cutoff
            IF ( FMID.LE.zero ) RTBIS = XMID
            IF(ABS(DX).LT.Accuracy .OR. FMID.EQ.0.) Finding = .false.
         ENDDO

!  Exception (too many bisections)

         if ( Finding ) Then
            write(c3,'(I3)')g
            fail = .true. ; message = 'Too many Bisections (540); Root not found, geometry #'//c3 ; return
         endif

!  Set final output if successful

         CC = CC + 1
         if ( doNadir(g) ) then
            RadCrit(g)   =  radii(NCrit_i-1) - RTBIS
         else
            AlphaCrit(g) = RTBIS ;  SX = sin(AlphaCrit(g))
            RadCrit(g)   = Raycon(g) / SX
            CotCrit(g)   = sqrt(one-SX*SX)/SX
         endif
         NCrit(g)     = NCrit_i

!  Continuation point for next geometry

67       continue

!  Zero the rest

         if ( NCrit(g) .ne. 0 ) nfinedivs(NCrit(g)+1:nlayers,g) = 0

!  Copy results for the rest of the azms

         do na = 2, nazms
            g1 = na_offset(ns,nv) +  na
            NCrit(g1)     = NCrit(g)
            AlphaCrit(g1) = AlphaCrit(g)
            RadCrit(g1)   = RadCrit(g)
            CotCrit(g1)   = CotCrit(g)
            nfinedivs(1:nlayers,g1) = nfinedivs(1:nlayers,g) 
         enddo

!  End geometry loops

      enddo
   enddo

!  Reset criticality, overall flag

   doCrit = ( CC.gt.0 )
 
!  Finish

   return
end subroutine SolarIn_EnhancedPS_Lattice_FindCrit

!

subroutine LosOut_EnhancedPS_QUpgrade &
       ( maxgeoms, maxlayers, maxfine, ngeoms,      & ! Input
         doNadir, radii, alpha, Raycon, nfinedivs,  & ! Input LOS layer quantities
         Ncrit, AlphaCrit, RadCrit,                 & ! Input Criticality variables
         radiifine, alphafine, xfine,               & ! Output, Adjusted fine-layering
         wfine, csqfine, cotfine )                    ! Output, Adjusted fine-layering

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!  Extension to Observational Geometries, 29 October 2012
!  Extension to Lattice       Geometries, 31 July    2013

!    starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                        - height grid, earth radius, Layer control
!                        - Critical Layer control

!  Regular Quadrature need not be done if LOSPATHS is set

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer, intent(in)    :: maxlayers, maxfine, maxgeoms

!  geometry numbers control

   integer, intent(in)    :: ngeoms

!  Flag for the Nadir case

   logical  , intent(in)  :: doNadir(maxgeoms)
 
!  Alphas, Radii, Ray constant

   real(ffp), intent(in)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii      (0:maxlayers)
   real(ffp), intent(in)  :: Raycon     (maxgeoms)

!  Critical stuff

   integer  , intent(in)  :: Ncrit     (maxgeoms)
   real(ffp), intent(in)  :: AlphaCrit (maxgeoms)
   real(ffp), intent(in)  :: RadCrit   (maxgeoms)

!  Fine divisions (these are already criticality-adjusted)

   integer  , intent(inout)  :: nfinedivs (maxlayers,maxgeoms)

!  Outputs
!  =======

!  Fine layering, Quadratures and local geometry, Adjusted values

   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: xfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: csqfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cotfine   (maxlayers,maxfine,maxgeoms)

!  Local
!  -----

   integer            :: n, n1, j, nfine, g
   real(ffp)          :: difh, csfine
   real(ffp)          :: tfine(maxfine), afine(maxfine)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Start LOS geometry loop

   do g = 1, ngeoms

!  Special case. Direct nadir viewing
!  ==================================

!  Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )
!    -- Adjust quadrature for the Critical layer

      if ( doNadir(g) ) then
         n = NCrit(g) ; difh = radii(n-1) - Radcrit(g) ; nfine = nfinedivs(n,g)
         call gauleg_ng (zero,difh,tfine,afine,nfine,maxfine)
         do j = 1, nfine
            radiifine(n,j,g) = radii(n-1) - tfine(j)
            xfine(n,j,g) = tfine(j)
            wfine(n,j,g) = afine(j)
         enddo

!  Outgoing sphericity geometry (General case)
!  ===========================================

      else if ( .not.doNadir(g) ) then
         n = NCrit(g) ; n1 = n - 1 ; nfine = nfinedivs(n,g)
         call gauleg_ng (alpha(n1,g),AlphaCrit(g),tfine,afine,nfine,maxfine)
         do j = 1,  nfine
            csfine = one / sin(tfine(j))
            radiifine(n,j,g) = raycon(g) * csfine
            alphafine(n,j,g) = tfine(j)
            xfine(n,j,g)   = radii(n1) - radiifine(n,j,g)
            wfine(n,j,g)   = afine(j)
            cotfine(n,j,g) = cos(tfine(j)) * csfine
            csqfine(n,j,g) = csfine * csfine
         enddo
       endif

!  End  eometry loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_QUpgrade

!

subroutine SolarIn_EnhancedPS_ObsGeom_SunPaths &
       ( maxgeoms, maxlayers, maxfine, do_Chapman,            & ! Input dimensions
         vsign, dtr, Pie, ngeoms, nlayers,                    & ! Input Control
         obsgeom_boa, doNadir, radii, alpha,                  & ! Input layer/level quantities
         nfinedivs, radiifine, alphafine,                     & ! Input Finelayer variables
         DoCrit, NCrit, RadCrit, AlphaCrit,                   & ! Input Criticality variables
         sunpaths, ntraverse, sunpathsfine, ntraversefine,    & ! Output
         chapfacs, Mu0, Mu1, cosscat, theta_all, phi_all )      ! Output

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the incoming Solar Beams
!     This is applicable to Both Upwelling and Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control
!    need also the complete values of all VZAs along outgoing paths

!  This routine has the fine gridding treatment

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

!  Dimensions and constants and flag

   integer  , intent(In)    :: maxgeoms, maxlayers, maxfine
   logical  , intent(In)    :: do_Chapman
   real(ffp), intent(In)    :: vsign, dtr, pie

!  BOA angles

   integer  , intent(In)    :: ngeoms
   real(ffp), intent(InOut) :: obsgeom_boa(maxgeoms,3)

!  Layer quantities

   integer  , intent(In)    :: nlayers
   logical  , intent(In)    :: doNadir  (maxgeoms)
   real(ffp), intent(In)    :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(In)    :: radii    (0:maxlayers)

!  Finelayer quantities

   integer  , intent(In)    :: nfinedivs  (maxlayers,maxgeoms)
   real(ffp), intent(In)    :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(In)    :: radiifine (maxlayers,maxfine,maxgeoms)

!  Criticality quantities

   Logical  , intent(In)    :: DoCrit
   integer  , intent(In)    :: NCrit     (maxgeoms)
   real(ffp), intent(In)    :: AlphaCrit (maxgeoms)
   real(ffp), intent(In)    :: RadCrit   (maxgeoms)

!  OUTPUTS
!  =======

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)

!  Fine level output (geometry)

   integer  , intent(Out)  :: ntraversefine(maxlayers,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine (maxlayers,maxlayers,maxfine,maxgeoms)

!  scattering angle cosines and associated angles

   real(ffp), intent(Out)  :: Mu0        (maxgeoms)
   real(ffp), intent(Out)  :: Mu1        (maxgeoms)
   real(ffp), intent(Out)  :: cosscat    (maxgeoms)
   real(ffp), intent(Out)  :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: phi_all    (0:maxlayers,maxgeoms)

!  Local

   logical       :: DirectSun, Do_OverheadSun, Do_ZeroSunBOA, Do_Normal
   integer       :: n, j, k, v
   real(ffp)     :: SolarDirection(3), Radstart, term1, term2
   real(ffp)     :: salpha_boa, calpha_boa, phi_boa_R, sphi_boa
   real(ffp)     :: theta_boa_R, stheta_boa, ctheta_boa, cphi_boa
   real(ffp)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle, diffhts(maxlayers)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Check
!   real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

   logical         :: DirectSunf(maxfine)
   real(ffp)       :: thetaf(maxfine)
   real(ffp)       :: sthetaf(maxfine)
   real(ffp)       :: cthetaf(maxfine)

!  Initialise output

   ntraverse = 0    ; ntraversefine = 0
   sunpaths  = zero ; sunpathsfine  = zero
   chapfacs  = zero

   phi_all = zero   ; theta_all = zero ; cosscat = zero

!  precompute

   do k = 1, nlayers
      diffhts(k) = radii(k-1) - radii(k)
   enddo

!  Start geometry loop

   do v = 1, ngeoms

!  Nominal traverse paths for Full illumination

      ntraverse(0,v) = 0
      do n = 1, nlayers
         ntraverse(n,v) = n
         do j = 1, nfinedivs(n,v)
            ntraversefine(n,j,v) = n
         enddo
      enddo

!  check range of inputs, already done.....

!  Special case

      Do_OverheadSun = obsgeom_boa(v,1).eq.zero

!  BOA angles

      if ( obsgeom_boa(v,2).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         salpha_boa  = sin(alpha(nlayers,v))
         calpha_boa  = cos(alpha(nlayers,v))
      endif
      Mu1(v) = calpha_boa

      theta_boa_R    = obsgeom_boa(v,1) * DTR
      if ( obsgeom_boa(v,1).eq.90.0_ffp ) then
         ctheta_boa     = zero
         stheta_boa     = one
      else
         stheta_boa     = sin(theta_boa_R)
         ctheta_boa     = cos(theta_boa_R)
      endif
      Mu0(v) = ctheta_boa

      phi_boa_R   = obsgeom_boa(v,3) * dtr
      cphi_boa    = cos(phi_boa_R)
      sphi_boa    = sin(phi_boa_R)

!  define Unit solar vector at BOA

      if ( Do_OverheadSun ) then
         SolarDirection = 0.0_ffp
      else
         SolarDirection(1) = - stheta_boa * cphi_boa * vsign
         SolarDirection(2) = - stheta_boa * sphi_boa
         SolarDirection(3) = - ctheta_boa
      endif

!  Cosine of scattering angle at boa

      if ( Do_OverheadSun ) then
         term1 = zero
         term2 = calpha_boa
         cosscat(v) = - vsign * term2 ; if (term2.eq.zero) cosscat(v) = term2
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2 + term1 
      endif

!  General case: LOS path in spherical geometry
!  ============================================

!  Start loop over all layers

      do n = nlayers, 1, -1

!  Special cases

        DO_ZeroSunBOA  = Do_OverheadSun.and.(n.eq.nlayers.or.doNadir(v))
        DO_Normal      = .not. doCrit .or. ( doCrit .and. n.le. NCrit(v) )

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( doCrit .and. n .eq. NCrit(v) ) then
              CumAngle = alpha(nlayers,v) - AlphaCrit(v) ; Radstart = RadCrit(v)
              call FindSun(DoNadir(v),Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                           theta_all(n,v),stheta,ctheta,DirectSun)
           else
              Radstart = radii(n)
              if ( n.eq. nlayers ) then
                 theta_all(n,v) = theta_boa_R ; stheta = stheta_boa ; ctheta = ctheta_boa ; DirectSun = .true.
              else
                 CumAngle = alpha(nlayers,v) - alpha(n,v)
                 call FindSun(DoNadir(v),Do_OverheadSun,radii(n),SolarDirection,CumAngle,theta_boa_R,&
                              theta_all(n,v),stheta,ctheta,DirectSun)
              endif
           endif
!        endif                                               !   @@RTSFix 9/5/12 (Comment out line)

!  Fine-layer sun positions

        if ( Do_Normal ) then
           do j = 1, nfinedivs(n,v)
              CumAngle = alpha(nlayers,v) - alphafine(n,j,v)
              call FindSun(DoNadir(v),Do_OverheadSun,radiifine(n,j,v),SolarDirection,CumAngle,theta_boa_R,&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
           enddo
        endif

!  Sun paths in layer

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( DirectSun ) then
              call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                theta_all(n,v),stheta,N,sunpaths(n,:,v))
           else
              call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_all(n,v),stheta,N,sunpaths(n,:,v),ntraverse(n,v))
           endif
        if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
           do j = 1, nfinedivs(n,v) 
              if ( DirectSunf(j) ) then
                 call FindSunPaths_D &
                  (Do_ZeroSunBOA,Maxlayers,Radiifine(n,j,v),Radii,&
                   thetaf(j),sthetaf(j),N,sunpathsfine(n,:,J,v))
              else
                 call FindSunPaths_T &
                  (Maxlayers,Pie,Radiifine(n,j,v),Radii,thetaf(j),sthetaf(j),N,sunpathsfine(n,:,J,v),ntraversefine(n,J,v))
              endif
!             if ( n.eq.14 ) write(*,*)j,n,Radiifine(n,j)-radii(n)
           enddo
        endif

!  debugging

!        if ( n.eq.14) then
!       sumd = SUM(sunpaths(n,1:ntraverse(n),v))
!       sth1 = stheta*RadCrit(v)/radii(0)
!       sume = sin(theta_all(n,v) - asin(sth1))*radii(0)/stheta
!       write(*,*)n,sumd,sume
!       do j = 1, nfinedivs(n)
!         sumd = SUM(sunpathsfine(n,1:ntraversefine(n,j,v),j,v))
!         sth1 = sthetaf(j)*radiifine(n,j,v)/radii(0)
!         sume = sin(thetaf(j) - asin(sth1))*radii(0)/sthetaf(j)
!         write(*,*)j,sumd,sume
!       enddo
!       pause
!      endif

!  Fix phi by using constancy of scatter angle
!     If AZM > 180, Subtract from 360 for consistency. (VLIDORT code, 10 October 2011)

        if (Do_OverheadSun.or.doNadir(v) ) then
           phi_all(n,v)     = obsgeom_boa(v,3) * dtr
        else
!           if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
              if ( doCrit .and. n .eq. NCrit(v) ) then
                 salpha = sin(AlphaCrit(v))
                 calpha = cos(AlphaCrit(v))
              else
                 salpha = sin(alpha(n,v))
                 calpha = cos(alpha(n,v))
              endif
              cphi = (cosscat(v)+vsign*calpha*ctheta)/stheta/salpha
              if ( cphi.gt.one)  cphi = one
              if ( cphi.lt.-one) cphi = -one
              phi_all(n,v)     = acos(cphi)
              if ( obsgeom_boa(v,3).gt.180.0_ffp) phi_all(n,v) = 2.0_ffp * Pie - phi_all(n,v)
!           endif                                               !   @@RTSFix 9/5/12 (Comment out line)
        endif

!  End layer loop

      enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

      DO_ZeroSunBOA  = Do_OverheadSun.and.doNadir(v)
      CumAngle = alpha(nlayers,v) - alpha(0,v) ; Radstart = radii(0)
      call FindSun(DoNadir(v),Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                   theta_all(0,v),stheta,ctheta,DirectSun)
      if (.not.DirectSun ) then
          call FindSunPaths_T (Maxlayers,Pie,Radii(0),Radii,theta_all(0,v),stheta,1,sunpaths(0,:,v),ntraverse(0,v))
      endif
      if ( Do_OverheadSun .or. doNadir(v) ) then
         phi_all(0,v)     = phi_boa_R
      else
         cphi = (cosscat(v)+vsign*calpha*ctheta)/stheta/salpha
         if ( cphi.gt.one)  cphi = one ; if ( cphi.lt.-one) cphi = -one
         phi_all(0,v)     = acos(cphi)
         if ( obsgeom_boa(v,3).gt.180.0_ffp) phi_all(0,v) = 2.0_ffp * Pie - phi_all(0,v)
      endif

!  Chapman factor calculations
!  ---------------------------

      if ( do_Chapman ) then
         do n = 1, nlayers
            call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,&
              theta_boa_R,stheta_boa,N,chapfacs(n,:,v))
            do k = 1, n
               chapfacs(n,k,v) = chapfacs(n,k,v)/(radii(k-1)-radii(k))
            enddo
         enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine SolarIn_EnhancedPS_ObsGeom_SunPaths

!

subroutine SolarIn_EnhancedPS_Lattice_SunPaths &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine, & ! Input dimensions
         do_Chapman, vsign, dtr, Pie,               & ! Input flags constants
         nvzas, nszas, nazms, na_offset, nlayers,   & ! Input # angles
         alpha_boa, theta_boa, phi_boa,             & ! Input BOA angles
         doNadir, radii, alpha,                     & ! Input layer/level quantities
         nfinedivs, radiifine, alphafine,           & ! Input Finelayer variables
         DoCrit, NCrit, RadCrit, AlphaCrit,         & ! Input Criticality variables
         sunpaths, ntraverse, sunpathsfine, ntraversefine,  & ! Output
         chapfacs, Mu0, Mu1, cosscat, theta_all, phi_all )    ! Output

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the incoming Solar Beams
!     This is applicable to Both Upwelling and Downwelling LOS-path geometries
!     No partials, this routine

!  Extension to Observational Geometries, 29 October 2012
!  Extension to Lattice       Geometries, 31 July    2013

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control
!    need also the complete values of all VZAs along outgoing paths

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions and constants and flag

   integer  , intent(In)    :: maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine
   logical  , intent(In)    :: do_Chapman
   real(ffp), intent(In)    :: vsign, dtr, pie

!  BOA angles

   integer  , intent(In)    :: nvzas, nszas, nazms, na_offset(maxszas,maxvzas)
   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  Layer quantities

   integer  , intent(In)    :: nlayers
   logical  , intent(In)    :: doNadir  (maxgeoms)
   real(ffp), intent(In)    :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(In)    :: radii    (0:maxlayers)

!  Finelayer quantities

   integer  , intent(In)    :: nfinedivs  (maxlayers,maxgeoms)
   real(ffp), intent(In)    :: alphafine  (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(In)    :: radiifine  (maxlayers,maxfine,maxgeoms)

!  Criticality quantities

   Logical  , intent(In)    :: DoCrit
   integer  , intent(In)    :: NCrit     (maxgeoms)
   real(ffp), intent(In)    :: AlphaCrit (maxgeoms)
   real(ffp), intent(In)    :: RadCrit   (maxgeoms)

!  Outputs
!  -------

!  main output:  solar paths, chapman factors and number of layers traversed

   integer  , intent(Out)   :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)   :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)   :: chapfacs   (  maxlayers,maxlayers,maxgeoms)

!  Fine level output: solar paths and number of layers traversed

   integer  , intent(Out)   :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(Out)   :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

!  scattering angle cosines and associated angles

   real(ffp), intent(Out)   :: Mu0        (maxgeoms)
   real(ffp), intent(Out)   :: Mu1        (maxgeoms)
   real(ffp), intent(Out)   :: cosscat    (maxgeoms)
   real(ffp), intent(Out)   :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)   :: phi_all    (0:maxlayers,maxgeoms)

!  Local
!  -----

   logical       :: DirectSun, Do_OverheadSun(maxszas), Do_ZeroSunBOA, Do_Normal
   integer       :: n, j, k, g, nv, ns, na
   real(ffp)     :: SolarDirection(3), Radstart, term1, term2, alpha_boa_R
   real(ffp)     :: salpha_boa(maxvzas), calpha_boa(maxvzas), phi_boa_R(maxazms), sphi_boa(maxazms)
   real(ffp)     :: theta_boa_R(maxszas), stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms)
   real(ffp)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle, diffhts(maxlayers)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Check
!   real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

   logical         :: DirectSunf(maxfine)
   real(ffp)       :: thetaf(maxfine)
   real(ffp)       :: sthetaf(maxfine)
   real(ffp)       :: cthetaf(maxfine)

!  Initialise output

   ntraverse = 0    ; ntraversefine = 0
   sunpaths  = zero ; sunpathsfine  = zero
   chapfacs  = zero

   phi_all = zero   ; theta_all = zero ; cosscat = zero

!  Pre-computation, BOA quantities

   do nv = 1, nvzas
      alpha_boa_R = alpha_boa(nv)*dtr
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)     = zero
         salpha_boa(nv)     = one
      else
         salpha_boa(nv)  = sin(alpha_boa_R)
         calpha_boa(nv)  = cos(alpha_boa_R)
      endif
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = theta_boa(ns).eq.zero
      theta_boa_R(ns)    = theta_boa(ns) * DTR
      if ( theta_boa(ns).eq.90.0_ffp ) then
         ctheta_boa(ns)     = zero
         stheta_boa(ns)     = one
      else
         stheta_boa(ns)     = sin(theta_boa_R(ns))
         ctheta_boa(ns)     = cos(theta_boa_R(ns))
      endif
      
   enddo

   do na = 1, nazms
      phi_boa_R(na)   = phi_boa(na) * dtr
      cphi_boa(na)    = cos(phi_boa_R(na))
      sphi_boa(na)    = sin(phi_boa_R(na))
   enddo

   do k = 1, nlayers
      diffhts(k) = radii(k-1) - radii(k)
   enddo

!  Start geometry loops
!  ====================

   do nv = 1, nvzas
      do ns = 1, nszas
         do na = 1, nazms

!  Geometry index for lattice

            g =  na_offset(ns,nv) + na

!  Cosines

            Mu0(g) = ctheta_boa(ns)
            Mu1(g) = calpha_boa(nv)
!            write(*,*)g, Mu0(g), Mu1(g)

!  Nominal number of Solar-path traverses for Full illumination

            ntraverse(0,g)         = 0
            do n = 1, nlayers
                ntraverse(n,g) = n
               do j = 1, nfinedivs(n,g)
                  ntraversefine(n,j,g) = n
               enddo
            enddo

!  define Unit solar vector at BOA

            if ( Do_OverheadSun(ns) ) then
               SolarDirection = 0.0_ffp
            else
               SolarDirection(1) = - stheta_boa(ns) * cphi_boa(na) * vsign
               SolarDirection(2) = - stheta_boa(ns) * sphi_boa(na)
               SolarDirection(3) = - ctheta_boa(ns)
            endif

!  Cosine of scattering angle at boa

            if ( Do_OverheadSun(ns) ) then
               term1 = zero
               term2 = calpha_boa(nv)
               cosscat(g) = - vsign * term2 ; if (term2.eq.zero) cosscat(g) = term2
            else
               term1 = salpha_boa(nv) * stheta_boa(ns) * cphi_boa(na)
               term2 = calpha_boa(nv) * ctheta_boa(ns)
               cosscat(g) = - vsign * term2 + term1 
            endif

!  General case: LOS path in spherical geometry
!  ============================================

!  Start loop over all layers

            do n = nlayers, 1, -1

!  Special cases

               DO_ZeroSunBOA  = Do_OverheadSun(ns).and.(n.eq.nlayers.or.doNadir(g))
               DO_Normal      = .not. doCrit .or. ( doCrit .and. n.le. NCrit(g) )

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

!            if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
               if ( doCrit .and. n .eq. NCrit(g) ) then
                  CumAngle = alpha(nlayers,g) - AlphaCrit(g) ; Radstart = RadCrit(g)
                  call FindSun(DoNadir(g),Do_OverheadSun(ns),Radstart,SolarDirection,CumAngle,theta_boa_R(ns),&
                           theta_all(n,g),stheta,ctheta,DirectSun)
               else
                  Radstart = radii(n)
                  if ( n.eq. nlayers ) then
                     theta_all(n,g) = theta_boa_R(ns) 
                     stheta = stheta_boa(ns) ; ctheta = ctheta_boa(ns) ; DirectSun = .true.
                  else
                     CumAngle = alpha(nlayers,g) - alpha(n,g)
                     call FindSun(DoNadir(g),Do_OverheadSun(ns),radii(n),SolarDirection,CumAngle,theta_boa_R(ns),&
                              theta_all(n,g),stheta,ctheta,DirectSun)
                  endif
               endif
!           endif                                               !   @@RTSFix 9/5/12 (Comment out line)

!  Fine-layer sun positions

               if ( Do_Normal ) then
                  do j = 1, nfinedivs(n,g)
                     CumAngle = alpha(nlayers,g) - alphafine(n,j,g)
                     call FindSun(DoNadir(g),Do_OverheadSun(ns),radiifine(n,j,g),SolarDirection,CumAngle,theta_boa_R(ns),&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
                  enddo
               endif

!  Sun paths in layer

!              if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
               if ( DirectSun ) then
                  call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                      theta_all(n,g),stheta,N,sunpaths(n,:,g))
               else
                  call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_all(n,g),stheta,N,sunpaths(n,:,g),ntraverse(n,g))
               endif
               if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
                  do j = 1, nfinedivs(n,g) 
                     if ( DirectSunf(j) ) then
                        call FindSunPaths_D &
                           (Do_ZeroSunBOA,Maxlayers,Radiifine(n,j,g),Radii,&
                            thetaf(j),sthetaf(j),N,sunpathsfine(n,:,j,g))
                     else
                        call FindSunPaths_T &
                           (Maxlayers,Pie,Radiifine(n,j,g),Radii,&
                                thetaf(j),sthetaf(j),N,sunpathsfine(n,:,j,g),ntraversefine(n,j,g))
                     endif
                  enddo
               endif

!  Fix phi by using constancy of scatter angle
!     If AZM > 180, Subtract from 360 for consistency. (VLIDORT code, 10 October 2011)

               if (Do_OverheadSun(ns).or.doNadir(g) ) then
                  phi_all(n,g)     = phi_boa_R(na)
               else
!              if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
                  if ( doCrit .and. n .eq. NCrit(g) ) then
                     salpha = sin(AlphaCrit(g))
                     calpha = cos(AlphaCrit(g))
                  else
                     salpha = sin(alpha(n,g))
                     calpha = cos(alpha(n,g))
                  endif
                  cphi = (cosscat(g)+vsign*calpha*ctheta)/stheta/salpha
                  if ( cphi.gt.one)  cphi = one
                  if ( cphi.lt.-one) cphi = -one
                  phi_all(n,g)     = acos(cphi)
                  if ( phi_boa(na).gt.180.0_ffp) phi_all(n,g) = 2.0_ffp * Pie - phi_all(n,g)
!              endif                                               !   @@RTSFix 9/5/12 (Comment out line)
               endif

!  End layer loop

            enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

            DO_ZeroSunBOA  = Do_OverheadSun(ns).and.doNadir(g)
            CumAngle = alpha(nlayers,g) - alpha(0,g) ; Radstart = radii(0)
            call FindSun(DoNadir(nv),Do_OverheadSun(ns),Radstart,SolarDirection,CumAngle,theta_boa_R(ns),&
                      theta_all(0,g),stheta,ctheta,DirectSun)
            if (.not.DirectSun ) then
               call FindSunPaths_T (Maxlayers,Pie,Radii(0),Radii,theta_all(0,g),stheta,1,sunpaths(0,:,g),ntraverse(0,g))
            endif
            if ( Do_OverheadSun(ns) .or. doNadir(g) ) then
               phi_all(0,g)     = phi_boa(na) * dtr
            else
               cphi = (cosscat(g)+vsign*calpha*ctheta)/stheta/salpha
               if ( cphi.gt.one)  cphi = one ; if ( cphi.lt.-one) cphi = -one
               phi_all(0,g)     = acos(cphi)
               if ( phi_boa(na).gt.180.0_ffp) phi_all(0,g) = 2.0_ffp * Pie - phi_all(0,g)
            endif

!  Chapman factor calculations
!  ---------------------------

            if ( do_Chapman ) then
               do n = 1, nlayers
                  call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii(n),Radii,&
                         theta_boa_R(ns),stheta_boa(ns),N,chapfacs(n,:,g))
                  do k = 1, n
                     chapfacs(n,k,g) = chapfacs(n,k,g)/diffhts(k)
                  enddo
               enddo
            endif
 
!  End geometry loops

         enddo
      enddo
   enddo

!  Finish

   return
end subroutine SolarIn_EnhancedPS_Lattice_SunPaths

!  Finish Module

end module FO_geometry_Routines_m

