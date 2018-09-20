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

module FO_DTGeometry_Master_m

!  Stand alone geometry for Direct Thermal only

!   Plane-parallel option, September 2012               (Version 1.2)
!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4). Not relevant here.

   use FO_geometry_Routines_m, only :    LosOut_EnhancedPS_Initial,  LosOut_EnhancedPS_Quadrature, &
                                         LosOut_EnhancedPS_QUpgrade, LosOut_EnhancedPS_FindCrit


private
public FO_DTGeometry_Master

contains

subroutine FO_DTGeometry_Master  &
       ( maxgeoms, maxlayers, maxfine, do_planpar, do_enhanced_ps,        & ! Input dimensions/flags
         ngeoms, nlayers, nfine, dtr, eradius, heights, alpha_boa,        & ! Inputs
         doNadir, doCrit, Acrit, extinc, Raycon, radii, alpha, cota,      & ! Input/Output(level)
         nfinedivs, xfine, wfine, csqfine, cotfine, alphafine, radiifine, & ! Output(Fine)
         NCrit, RadCrit, CotCrit, Mu1,                                    & ! Output(Critical/Cos)
         fail, message, trace )                                             ! Output(Status)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxlayers, maxfine, maxgeoms

!  Flags (sphericity flag is mutually exclusive).

   logical  , intent(in)     :: do_enhanced_ps
   logical  , intent(in)     :: do_planpar

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)     :: ngeoms, nlayers, nfine

!  dtr = degrees-to-Radians

   real(ffp), intent(in)     :: dtr

!  Radius + heights

   real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees). Enough information for Lattice or Obsgeom.

   real(ffp), intent(InOut)  :: alpha_boa(maxgeoms)

!  Critical adjustment for cloud layers

   logical  , intent(inout)  :: doCrit
   real(ffp), intent(in)     :: extinc(maxlayers)
   real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments
!  ======================

!  LOSPATHS flag has been removed now.  8/1/13

!  Flag for the Nadir case.

   logical  , intent(inout)  :: doNadir(maxgeoms)
  
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout)

   real(ffp), intent(inout)  :: radii    (0:maxlayers)
   real(ffp), intent(inout)  :: Raycon   (maxgeoms)
   real(ffp), intent(inout)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cota     (0:maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS, fine layering output

   integer  , intent(inout)  :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: xfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: csqfine  (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cotfine  (maxlayers,maxfine,maxgeoms)

   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Output arguments
!  ================

!  Critical layer

   integer  , intent(out)  :: Ncrit(maxgeoms)
   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Cosine

   real(ffp), intent(Out)  :: Mu1     (maxgeoms)

!  Exception handling

   logical      , intent(out)    :: fail
   character*100, intent(out)    :: message
   character*100, intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

   real(ffp)  :: Lospaths (maxlayers,maxgeoms)

!  Other angles

   real(ffp)  :: cosa       (0:maxlayers,maxgeoms)
   real(ffp)  :: sina       (0:maxlayers,maxgeoms)

!  Critical values

   real(ffp)  :: AlphaCrit(maxgeoms)

!  Help variables
!  --------------

   logical                :: do_regular_ps
   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp)              :: cutoff, alpha_boa_R
   character*2            :: c2
   integer                :: v

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   Mu1 = zero ; NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
   cutoff = zero; if (ACrit.gt.zero) cutoff = -log(ACrit)

!  check sphericity control
!  ------------------------

!  Cannot have Plane-parallel and Enhanced PS

   if ( do_planpar .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Enhanced PS options'
      trace   = 'Initial Flag Check in FO_DTGeometry_Master'
      fail    = .true. ;  return
   endif

!  Other sphericity flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   
!  Check geometry angles
!  ---------------------

!  VZA can be 0-90 degrees inclusive, but not outside this range

   do v = 1, ngeoms
      if ( alpha_boa(v).gt.90.0_ffp.or.alpha_boa(v).lt.zero ) then
         write(c2,'(I2)')v
         message = 'Boa LOS angle outside range [0,90]); Check it!'
         trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_DTGeometry_Master'
         fail    = .true. ;  return
      endif
    enddo

!  set cosines if you are here

   do v = 1, ngeoms
      alpha_boa_R  = alpha_boa(v) * DTR
      Mu1(v)       = cos(alpha_boa_R)
   enddo

!  Plane-parallel or regular-PS, just set alpha and return
!  -------------------------------------------------------

   if ( do_planpar .or. do_regular_ps ) then
      alpha = zero
      do v = 1, ngeoms
         alpha_boa_R      = alpha_boa(v) * DTR
         alpha(1:nlayers,v) = alpha_boa_R
      enddo
      return
   endif

!  Enhanced PS; proceed in 3 Steps
!  ===============================

!  Step 1; Initial LOS-path quantities, OUTGOING Beam
!  --------------------------------------------------

!    1a.  Given heights and BOA LOS angles, compute path angles and radii
!    1b.  LOS fine-layer quadratures. Non-adjusted, no Criticality
!    1c.  Find Critical-layer adjustments (Optional)

   CALL LosOut_EnhancedPS_Initial &
          ( maxgeoms, maxlayers, dtr, ngeoms, nlayers, & ! Input
            heights, eradius, alpha_boa,               & ! Input
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
         trace = 'Error from LosOut_EnhancedPS_FindCrit in FO_DTGeometry_Master' ; return
      endif
   endif

   if ( doCrit) then
      call LosOut_EnhancedPS_QUpgrade &
          ( maxgeoms, maxlayers, maxfine, ngeoms,      & ! Input
            doNadir, radii, alpha, Raycon, nfinedivs,  & ! Input LOS layer quantities
            Ncrit, AlphaCrit, RadCrit,                 & ! Input Criticality variables
            radiifine, alphafine, xfine,               & ! Output, Adjusted fine-layering
            wfine, csqfine, cotfine )                    ! Output, Adjusted fine-layering
   endif
 
!  Finish

   return
end subroutine FO_DTGeometry_Master

!  Finish

end module FO_DTGeometry_Master_m

