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
! #   Version 1.3,  29 October  2012, Obsgeom Multi-geom.   #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_SSGeometry_Master_m

!  Stand alone geometry for solar scattering only

!   Plane-parallel option, September 2012               (Version 1.2)
!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)

   use FO_geometry_Generic_m
   use FO_geometry_Routines_m

private
public FO_SSGeometry_Master

contains

subroutine FO_SSGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine,            & ! Input dimensions
         do_obsgeom, do_Chapman, do_planpar, do_enhanced_ps,                 & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, dtr, Pie, vsign,       & ! Input control and constants
         eradius, heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,       & ! Input
         doNadir, doCrit, Acrit, extinc, Raycon, radii, alpha, cota,         & ! Input/Output(level)
         nfinedivs, xfine, wfine, csqfine, cotfine, alphafine, radiifine,    & ! Output(Fine)
         NCrit, RadCrit, CotCrit, Mu0, Mu1, cosscat, chapfacs,               & ! Output(Crit/scat)
         sunpaths, ntraverse, sunpathsfine, ntraversefine,                   & ! Output(Sunpaths)
         fail, message, trace )                                                ! Output(Status)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine

!  Flags (sphericity flag is mutually exclusive). Obsgeom flag new 7/31/13

   logical  , intent(in)     :: do_ObsGeom
   logical  , intent(in)     :: do_Chapman
   logical  , intent(in)     :: do_enhanced_ps
   logical  , intent(in)     :: do_planpar

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)     :: ngeoms, nszas, nvzas, nazms, nlayers, nfine

!  dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

   real(ffp), intent(in)     :: dtr, Pie, vsign

!  Radius + heights

   real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees). Enough information for Lattice or Obsgeom.
!   Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)
!    In both cases, the Phi angle may be changed.....

   real(ffp), intent(InOut)  :: Obsgeom_boa(maxgeoms,3)
   real(ffp), intent(InOut)  :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  Critical adjustment control for cloud layers

   logical  , intent(inout)  :: doCrit
   real(ffp), intent(in)     :: extinc(maxlayers)
   real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments
!  ======================

!  LOSPATHS flag has been removed now.  8/1/13

!  Flag for the Nadir case

   logical  , intent(inout)  :: doNadir(maxgeoms)
  
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout), input if DO_LOSpaths set)
!    WARNING: Adjusted geometry will require maxgeoms dimension

   real(ffp), intent(inout)  :: radii    (0:maxlayers)
   real(ffp), intent(inout)  :: Raycon   (maxgeoms)
   real(ffp), intent(inout)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cota     (0:maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS, find layering output
!    WARNING: Adjusted geometry will require maxgeoms dimension

   integer  , intent(inout)  :: nfinedivs (maxlayers,maxgeoms)

   real(ffp), intent(inout)  :: xfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: csqfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cotfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Output arguments
!  ================

!  Critical layer
!    WARNING: Adjusted geometry will require maxgeoms dimension

   integer  , intent(out)  :: Ncrit(maxgeoms)
   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  solar paths and Chapman factors

   integer  , Intent(out)  :: ntraverse     (0:maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)

   integer  , Intent(out)  :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(out)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(out)  :: Chapfacs      (maxlayers,maxlayers,maxgeoms)

!  Cosine scattering angle, other cosines

   real(ffp), Intent(out)  :: cosscat (maxgeoms)
   real(ffp), intent(Out)  :: Mu0     (maxgeoms)
   real(ffp), intent(Out)  :: Mu1     (maxgeoms)

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

   real(ffp)  :: Lospaths (maxlayers,maxgeoms)

!  Other angles and cos/sin values.

   real(ffp)  :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp)  :: phi_all    (0:maxlayers,maxgeoms)
   real(ffp)  :: cosa       (0:maxlayers,maxgeoms)
   real(ffp)  :: sina       (0:maxlayers,maxgeoms)

!  Critical values

   real(ffp)  :: AlphaCrit(maxgeoms)

!  1-d equivalents (Only for the Lattice option)
!  ---------------------------------------------

!  LOS layer/level quantities

   real(ffp)  :: Lospaths_LOS (maxlayers,maxvzas)

   logical    :: doNadir_LOS  (maxvzas)
   real(ffp)  :: Raycon_LOS   (maxvzas)
   real(ffp)  :: alpha_LOS    (0:maxlayers,maxvzas)
   real(ffp)  :: cota_LOS     (0:maxlayers,maxvzas)
   real(ffp)  :: cosa_LOS     (0:maxlayers,maxvzas)
   real(ffp)  :: sina_LOS     (0:maxlayers,maxvzas)

!  LOS Quadratures for Enhanced PS, find layering output

   integer    :: nfinedivs_LOS (maxlayers,maxvzas)
   real(ffp)  :: xfine_LOS     (maxlayers,maxfine,maxvzas)
   real(ffp)  :: wfine_LOS     (maxlayers,maxfine,maxvzas)
   real(ffp)  :: csqfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: cotfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: alphafine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: radiifine_LOS (maxlayers,maxfine,maxvzas)

   real(ffp)  :: AlphaCrit_LOS(maxvzas)
   integer    :: Ncrit_LOS(maxvzas)
   real(ffp)  :: RadCrit_LOS(maxvzas), CotCrit_LOS(maxvzas)

!  Help variables
!  --------------

   logical                :: do_regular_ps
   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp
   real(ffp)              :: cutoff
   character*2            :: c2
   integer                :: v, ns, nv, nv_offset(maxszas),na_offset(maxszas,maxvzas)

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
   cutoff = zero; if (ACrit.gt.zero) cutoff = -log(ACrit)

!  check sphericity control
!  ------------------------

!  Cannot have Plane-parallel and Enhanced PS

   if ( do_planpar .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Enhanced PS options'
      trace   = 'Initial Flag Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Other sphericity flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   
!  Check geometry angles
!  ---------------------

!  1. OBSERVATIONAL GEOMETRIES

   if ( do_obsgeom ) then

      do v = 1, ngeoms

         if ( obsgeom_boa(v,2).gt.90.0_ffp.or.obsgeom_boa(v,2).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Boa LOS angle outside range [0,90]); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1, OBSGEOM mode'
            fail    = .true. ;  return
         endif

         if ( do_planpar ) then
            if ( obsgeom_boa(v,1).ge.90.0_ffp.or.obsgeom_boa(v,1).lt.zero ) then
               write(c2,'(I2)')v
               message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
               trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1, OBSGEOM mode'
               fail    = .true. ;  return
            endif
         else
            if ( obsgeom_boa(v,1).gt.90.0_ffp.or.obsgeom_boa(v,1).lt.zero ) then
               write(c2,'(I2)')v
               message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!'
               trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1, OBSGEOM mode'
               fail    = .true. ;  return
            endif
         endif

!  PHI is not limited to <= 180 degs. Also, not negative.
!     Old Code :     if ( phi_boa.gt.180.0_ffp ) phi_boa = 360.0_ffp - phi_boa

         if ( obsgeom_boa(v,3).lt.zero )   obsgeom_boa(v,3) = - obsgeom_boa(v,3)
         if ( obsgeom_boa(v,3).gt.360.0_ffp ) obsgeom_boa(v,3) = 360.0_ffp - obsgeom_boa(v,3) - 360.0_ffp

!  End loop over Obs geometries

      enddo

!  2. LATTICE GEOMETRIES

   else

!  VZA can be 0-90 degrees inclusive, but not outside this range

      do v = 1, nvzas
         if ( alpha_boa(v).gt.90.0_ffp.or.alpha_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Boa LOS angle outside range [0,90]); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1, LATTICE MODE'
            fail    = .true. ;  return
         endif
      enddo

!  SZA can be 0-90 degrees inclusive, but not outside this range
!    For plane-parallel, 90 degrees is not allowed

   do v = 1, nszas
      if ( do_planpar ) then
         if ( theta_boa(v).ge.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1, LATTICE MODE'
            fail    = .true. ;  return
         endif
      else
         if ( theta_boa(v).gt.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!, LATTICE MODE'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1'
            fail    = .true. ;  return
         endif
      endif
   enddo

!  PHI is not limited to <= 180 degs. Also, not negative.
!     Old Code :     if ( phi_boa.gt.180.0_ffp ) phi_boa = 360.0_ffp - phi_boa

      do v = 1, nazms
         if ( phi_boa(v).lt.zero )   phi_boa(v) = - phi_boa(v)
         if ( phi_boa(v).gt.360.0_ffp ) phi_boa(v) = 360.0_ffp - phi_boa(v) - 360.0_ffp
      enddo

!  Offsets

      do ns = 1, nszas
         nv_offset(ns) = nvzas * nazms * (ns - 1) 
         do nv = 1, nvzas
            na_offset(ns,nv) = nv_offset(ns) + nazms * (nv - 1)
         enddo
      enddo

!  End Checking

   endif

!  Plane-parallel, One routine only
!  --------------------------------

   if ( do_planpar ) then
      if ( do_obsgeom ) then
         CALL Obsgeom_PlanPar &
          ( maxgeoms, maxlayers, dtr, vsign, do_Chapman,                & ! Inputs
            ngeoms, nlayers, obsgeom_boa, heights,                      & ! Inputs
            alpha, sunpaths, ntraverse, chapfacs, cosscat )               ! Outputs
      else
         CALL Lattice_PlanPar &
          ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, dtr, vsign, & ! Input
            do_Chapman, nszas, nvzas, nazms, nv_offset, na_offset,      & ! Input
            nlayers, alpha_boa, theta_boa, phi_boa, heights,            & ! Inputs
            alpha, sunpaths, ntraverse, chapfacs, cosscat )               ! Outputs
      endif
      return
   endif

!  Regular PS, One routine only
!  ----------------------------

   if ( do_regular_ps ) then
      if ( do_obsgeom ) then
         CALL Obsgeom_RegularPS &
          ( maxgeoms, maxlayers, dtr, vsign, do_Chapman,                  & ! Inputs
            ngeoms, nlayers, obsgeom_boa, heights, eradius,               & ! Inputs
            Raycon, radii, alpha, sunpaths, ntraverse, chapfacs, cosscat )  ! Outputs
      else
         CALL Lattice_RegularPS &
          ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, dtr, vsign,     & ! Input
            do_Chapman, nszas, nvzas, nazms, nv_offset, na_offset, nlayers, & ! I  nput
            alpha_boa, theta_boa, phi_boa, heights, eradius,                & ! Inputs
            Raycon, radii, alpha, sunpaths, ntraverse, chapfacs, cosscat )    ! Outputs
      endif
      return
   endif

!  Enhanced PS; proceed in 4 Steps
!  ===============================

   doCrit = .false.

!  Step 1; Initial LOS-path quantities, OUTGOING Beam
!  --------------------------------------------------

!    1a.  Given heights and BOA LOS angles, compute path angles and radii
!    1b.  LOS fine-layer quadratures. Non-adjusted, no Criticality
!    1c.  Find Critical-layer adjustments (Optional)

!  OBSGEOM mode
!  ============

   if ( do_Obsgeom ) then

      CALL LosOut_EnhancedPS_Initial &
          ( maxgeoms, maxlayers, dtr, ngeoms, nlayers, & ! Input
            heights, eradius, obsgeom_boa(:,2),        & ! Input
            doNadir, radii, Raycon, Lospaths,          & ! Output
            alpha, sina, cosa, cota )                    ! Output

      CALL LosOut_EnhancedPS_Quadrature &
       ( maxgeoms, maxlayers, maxfine,     & ! Input 
         ngeoms, nlayers, nfine,           & ! Input
         doNadir, radii, alpha, Raycon,    & ! Input
         nfinedivs, radiifine, alphafine,  & ! Output
         xfine, wfine, csqfine, cotfine )    ! Output

      if ( doCrit) then
         CALL LosOut_EnhancedPS_FindCrit &
          ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
            extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
            Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs
         if ( Fail ) then
            trace = 'Error from LosOut_EnhancedPS_FindCrit in SS_Geometry_1, OBSGEOM mode' ; return
         endif
      endif

   endif

!  LATTICE MODE
!  ============

   if ( .not.do_Obsgeom ) then

      CALL LosOut_EnhancedPS_Initial &
          ( maxvzas, maxlayers, dtr, nvzas, nlayers,      & ! inputs
            heights, eradius, alpha_boa,                  & ! Inputs
            doNadir_LOS, radii, Raycon_LOS, Lospaths_LOS, & ! Output
            alpha_LOS, sina_LOS, cosa_LOS, cota_LOS )       ! Output

      CALL LosOut_EnhancedPS_Quadrature &
          ( maxvzas, maxlayers, maxfine, nvzas, nlayers, nfine,  & ! Inputs
            doNadir_LOS, radii, alpha_LOS, Raycon_LOS,           & ! Inputs
            nfinedivs_LOS, radiifine_LOS, alphafine_LOS,         & ! Output
            xfine_LOS, wfine_LOS, csqfine_LOS, cotfine_LOS )       ! Output

      if ( doCrit) then
         CALL LosOut_EnhancedPS_FindCrit &
          ( maxvzas, maxlayers, nvzas, nlayers, Acrit, Cutoff, doNadir_LOS,      & ! Inputs
            extinc, Lospaths_LOS, sina_LOS, cosa_LOS, radii, nfinedivs_LOS,      & ! Inputs
            Ncrit_LOS, AlphaCrit_LOS, RadCrit_LOS, CotCrit_LOS, fail, message )    ! Outputs
         if ( Fail ) then
            trace = 'Error from LosOut_EnhancedPS_FindCrit in SS_Geometry_1, LATTICE mode' ; return
         endif
      endif

!  Transform LOS output to GEOMS output, but leave criticality stuff until later on....

      call LOSCopy1_EnhancedPS_Lattice &
       ( maxgeoms, maxszas, maxvzas, maxlayers, maxfine,          & ! Input dimensions
         nvzas, nszas, nazms, nlayers, na_offset,                 & ! Input   (# angles/layers)
         doNadir_LOS, Raycon_LOS, alpha_LOS, cota_LOS,            & ! Input LOS
         nfinedivs_LOS, radiifine_LOS, alphafine_LOS,             & ! Input LOS
         xfine_LOS, wfine_LOS, csqfine_LOS, cotfine_LOS,          & ! Input LOS
         doNadir, Raycon, alpha, cota,                            & ! Output
         nfinedivs, radiifine, alphafine,                         & ! Output
         xfine, wfine, csqfine, cotfine )                           ! Output

!  End lattice option

   endif

!  Step 2, INCOMING SOLAR BEAMS.
!  -----------------------------
   
!  2a.  Incoming,  Find Critical-layer adjustments       (Optional)
!  2b.  Critical-layer Upgrade of Quadrature done here. (Optional)
!  2c.  Incoming, solar pathlengths

!  Step 2a
!  =======

   if ( doCrit) then
      if ( do_obsgeom ) then
         call SolarIn_EnhancedPS_Obsgeom_FindCrit &
          ( maxgeoms, maxlayers, ngeoms, nlayers, doNadir, dtr, Acrit, & ! Input
            cutoff, alpha, radii, extinc, Raycon, obsgeom_boa(:,1),    & ! Inputs
            doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit,     & ! Outputs
            fail, message )                                              ! Outputs
         if ( Fail ) then
            trace = 'Error from SolarIn_EnhancedPS_Obsgeom_FindCrit in SS_Geometry_1, OBSGEOM mode' ; return
         endif
      else
         call SolarIn_EnhancedPS_Lattice_FindCrit &
          ( maxgeoms, maxvzas, maxszas, maxlayers, ACrit, dtr, cutoff,    & ! Input   (dimensions/constants)
            nvzas, nszas, nazms, na_offset, nlayers, theta_boa,           & ! Input   (# angles/layers)
            doNadir, alpha, radii, extinc, Raycon, nfinedivs,             & ! Inputs  (LOS and extinction)
            Ncrit_los, AlphaCrit_los, RadCrit_los, CotCrit_los,           & ! Inputs  (LOS criticality)
            doCrit, Ncrit, AlphaCrit, RadCrit, CotCrit,                   & ! Outputs (Incoming criticality)
            fail, message )                                                 ! Outputs (Exception handling)
         if ( Fail ) then
            trace = 'Error from SolarIn_EnhancedPS_Lattice_FindCrit in SS_Geometry_1, LATTICE mode' ; return
         endif
      endif
   endif

!  Step 2b
!  =======

!  valid for both options

   if ( doCrit) then
      call LosOut_EnhancedPS_QUpgrade &
          ( maxgeoms, maxlayers, maxfine, ngeoms,      & ! Input
            doNadir, radii, alpha, Raycon, nfinedivs,  & ! Input LOS layer quantities
            Ncrit, AlphaCrit, RadCrit,                 & ! Input Criticality variables
            radiifine, alphafine, xfine,               & ! Output, Adjusted fine-layering
            wfine, csqfine, cotfine )                    ! Output, Adjusted fine-layering
   endif

!  Step 2c
!  =======

!  Lattice routine is for UNADJUSTED GEOMETRIES, at the moment

   if ( do_obsgeom ) then
      CALL SolarIn_EnhancedPS_ObsGeom_SunPaths &
       ( maxgeoms, maxlayers, maxfine, do_Chapman,            & ! Input dimensions
         vsign, dtr, Pie, ngeoms, nlayers,                    & ! Input Control
         obsgeom_boa, doNadir, radii, alpha,                  & ! Input layer/level quantities
         nfinedivs, radiifine, alphafine,                     & ! Input Finelayer variables
         DoCrit, NCrit, RadCrit, AlphaCrit,                   & ! Input Criticality variables
         sunpaths, ntraverse, sunpathsfine, ntraversefine,    & ! Output
         chapfacs, Mu0, Mu1, cosscat, theta_all, phi_all )      ! Output
   else

      CALL SolarIn_EnhancedPS_Lattice_SunPaths &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine, & ! Input dimensions
         do_Chapman, vsign, dtr, Pie, nvzas, nszas, nazms,        & ! Input Control
         na_offset, nlayers, alpha_boa, theta_boa, phi_boa,       & ! Input BOA angles
         doNadir, radii, alpha,                                   & ! Input layer/level quantities
         nfinedivs, radiifine, alphafine,                         & ! Input Finelayer variables
         DoCrit, NCrit, RadCrit, AlphaCrit,                       & ! Input Criticality variables
         sunpaths, ntraverse, sunpathsfine, ntraversefine,        & ! Output
         chapfacs, Mu0, Mu1, cosscat, theta_all, phi_all )          ! Output
   endif

!  Finish

   return
end subroutine FO_SSGeometry_Master


!  Finish

end module FO_SSGeometry_Master_m


