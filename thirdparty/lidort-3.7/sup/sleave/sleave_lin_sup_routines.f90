! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -          -         #
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
! #  This Version :   3.7 F90                               #
! #  Release Date :   June 2014                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED     (3.2)        #
! #       NEW: OUTGOING SPHERICITY CORRECTION  (3.2)        #
! #       NEW: TOTAL COLUMN JACOBIANS          (3.3)        #
! #       VLIDORT COMPATIBILITY                (3.4)        #
! #       THREADED/OPTIMIZED F90 code          (3.5)        #
! #       EXTERNAL SS / NEW I/O STRUCTURES     (3.6)        #
! #                                                         #
! #       Surface-leaving, BRDF Albedo-scaling     (3.7)    # 
! #       Taylor series, BBF Jacobians, ThreadSafe (3.7)    #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # WaterLeaving Subroutines in this Module                     #
! #                                                             #
! #         Linearized_WaterLeaving (Top-level)                 #
! #           - Lin_Ocean_Reflectance (ex. MORCASIWAT)          #
! #           - Lin_WhiteCap_Reflectance                        #
! #           - Lin_Water_Transmittance                         #
! #               * Lin_GENERAL_SUNGLINT                        #
! #                                                             #
! ###############################################################

      MODULE sleave_lin_sup_routines_m

!  Version 3.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 3.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 3.6.
!  For Version 3.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the BRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the BRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      use sleave_sup_aux_m
      use sleave_sup_routines_m, ONLY :  Water_RefracIndex,         &
                                         Water_Transmittance_Quads, &
                                         Fresnel_Complex

      PRIVATE
      PUBLIC :: Linearized_WaterLeaving

      CONTAINS

subroutine Linearized_WaterLeaving &
     ( Maxszas, Maxvzas, Maxstreams,                                 &
       do_Isotropy, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy, &
       Wavelength, Salinity, PigmentConc, Windspeed, WindSunAngles,  &
       nszas, nvzas, nstreams, szas, vzas, streams,                  &
       WLeaving_ISO, WLeaving_SD, WLeaving_SV,                       &
       dC_WLeaving_ISO, dC_WLeaving_SD, dC_WLeaving_SV,              &
       dW_WLeaving_ISO, dW_WLeaving_SD, dW_WLeaving_SV )

! mick fix 12/28/2014 - Using normalization correction by A Sayer, 04 Nov 2014
! Apply a normalisation factor of 1/(pi/mu0) to output water-leaving reflectance, to
! bring things in line with results from e.g. 6S simulations and expected behaviour.
! Think this is a subtlety related to reflectance vs. normalised radiance treatment,
! although it is very obvious if you don't do it. Correction applied at end of the
! subroutine.

   IMPLICIT NONE
   INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  Newly constructed Water-Leaving supplement routine. Version 3.7 LIDORT code.
!    23-24 April 2014. R. Spurr, RT Solutions Inc.
!    Validated against the Modified 6S code by A. Sayer.. 

!  This is a Stand-alone subroutine.

!  inputs
!  ======

!  Dimensioning

   integer  , intent(in)   :: Maxszas, Maxvzas, Maxstreams

!  Logical flags
!  -------------

!  Isotropic (Fast Calculation) option, assumes all transmittances = 1

   Logical  , intent(in)   :: do_Isotropy

!  Optional inclusion of Foam term

   Logical  , intent(in)   :: Do_FoamOption

!  Optional inclusion of Shadow term for Glitter (Air-water only?)

   Logical  , intent(in)   :: Do_GlintShadow

!  Flag for using Isotropic Facet distribution

   LOGICAL  , intent(in)   :: Do_FacetIsotropy

!  Physical
!  --------

!  Wavelength in Micrometers

   real(fpk), intent(in)   :: Wavelength

!  Salinity

   real(fpk), intent(in)   :: Salinity

!  Pigment concentration

   real(fpk), intent(in)   :: PigmentConc

!  Windspeed m/s

      REAL(fpk), intent(in)    :: WINDSPEED

!  Azimuth angles in Radians

      REAL(fpk), intent(in)    :: WindSunAngles ( Maxszas )

!  Sun, viewing and stream angles
!  ------------------------------

   integer  , intent(in) :: nszas, nvzas, nstreams
   real(fpk), intent(in) :: szas   (Maxszas)
   real(fpk), intent(in) :: vzas   (Maxvzas)
   real(fpk), intent(in) :: streams(Maxstreams)

!  OUTPUT
!  ======

!  Isotropic value. Fast calculation

   REAL(fpk), intent(out)    :: WLeaving_ISO    ( Maxszas )
   REAL(fpk), intent(out)    :: dC_WLeaving_ISO ( Maxszas )
   REAL(fpk), intent(out)    :: dW_WLeaving_ISO ( Maxszas )

!  Input solar, output stream angles

   REAL(fpk), intent(out)    :: WLeaving_SD    ( Maxszas, Maxstreams )
   REAL(fpk), intent(out)    :: dC_WLeaving_SD ( Maxszas, Maxstreams )
   REAL(fpk), intent(out)    :: dW_WLeaving_SD ( Maxszas, Maxstreams )

!  input solar, output view angles

   REAL(fpk), intent(out)    :: WLeaving_SV    ( Maxszas, Maxvzas )
   REAL(fpk), intent(out)    :: dC_WLeaving_SV ( Maxszas, Maxstreams )
   REAL(fpk), intent(out)    :: dW_WLeaving_SV ( Maxszas, Maxstreams )

!  remark: Still no Azimuth dependence here.....

!  LOCAL
!  =====

!  Transmittance Quadratures
!  -------------------------

   integer, parameter   :: Max_PolarQuads=24, Max_AzimQuads=48
   real(fpk) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk) :: CosPolarQuads (Max_PolarQuads)
   real(fpk) :: SinPolarQuads (Max_PolarQuads)
   real(fpk) :: PolarWeights  (Max_PolarQuads)

   real(fpk) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk) :: CosAzimQuads (Max_AzimQuads)
   real(fpk) :: SinAzimQuads (Max_AzimQuads)
   real(fpk) :: AzimWeights  (Max_AzimQuads)

!  Glitter control
!  ---------------

   logical   :: do_coeffs, local_do_shadow
   REAL(fpk) :: SUNGLINT_COEFFS(7), dSUNGLINT_COEFFS(7), Refrac_R, Refrac_I, Refrac_sq
   real(fpk) :: phi_w, cphi_w, sphi_w, local_refrac_R, local_refrac_I

!  other variables
!  ---------------

!  Help

   logical   :: noWL
   integer   :: J, I
   real(fpk) :: dtr, pi, Albedo, Const0, Const, Rwprime, denom, incident_angle, Local_Sine
   real(fpk) :: dW_Const0, dC_Rwprime, dW_Const, dC_Const
   real(fpk) :: Foam_correction, WC_Reflectance, WC_Lambertian
   real(fpk) :: dW_Foam_correction, dWC_Reflectance, dWC_Lambertian
   real(fpk) :: Ocean_Reflecs(Maxszas), dC_Ocean_Reflecs(Maxszas), SZA_cosines(Maxszas)

!  transmittance help variables

   logical   :: do_transmittances
   real(fpk) :: Trans_Norm, Trans_solar(Maxszas), dW_Trans_solar(Maxszas)
   real(fpk) :: Trans_stream(Maxstreams), Trans_Viewing(Maxvzas)
   real(fpk) :: dW_Trans_stream(Maxstreams), dW_Trans_Viewing(Maxvzas)

!  Parameters

   real(fpk), parameter :: zero = 0.0_fpk
   real(fpk), parameter :: one  = 1.0_fpk

!  Initial Setup
!  -------------

!  conversions

   pi = acos(-1.0_fpk)
   dtr = pi / 180.0_fpk
   SZA_cosines = zero
   do J = 1, nszas
      SZA_cosines(J) = cos(szas(J)*dtr)
   enddo

!  Zero the output

   WLeaving_ISO = zero ; dC_WLeaving_ISO = zero ; dW_WLeaving_ISO = zero
   WLeaving_SD  = zero ; dC_WLeaving_SD  = zero ; dW_WLeaving_SD  = zero
   WLeaving_SV  = zero ; dC_WLeaving_SV  = zero ; dW_WLeaving_SV  = zero

!  Refractive index. Formerly INDWAT

   Call  Water_RefracIndex  ( Wavelength, Salinity, Refrac_R, Refrac_I )
   Refrac_sq = Refrac_R * Refrac_R + Refrac_I * Refrac_I

!  Ocean Leaving ; Exit if no conribution. Formerly MORCASIWAT

   call  Lin_Ocean_Reflectance &
     ( Maxszas, nszas, SZA_cosines, Wavelength, PigmentConc, noWL, &
       Ocean_Reflecs, dC_Ocean_Reflecs )
   if ( noWL ) return

!  Foam-reflectance correction.

   Foam_correction = one ; dW_Foam_correction = zero
   if ( Do_FoamOption ) then
      call Lin_WhiteCap_Reflectance &
         ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian, &
                        DWC_Reflectance, DWC_Lambertian )
      Foam_correction    = one - WC_Reflectance
      dW_Foam_correction =     - dWC_Reflectance
   endif
   
!  Initial debug...........
!   write(*,*)Ocean_Reflecs(1),DC_Ocean_Reflecs(1)
!   write(*,*)Foam_correction, dW_Foam_correction
!   pause'after first routines'

!  set  Coeffs flag, initialize local shadow flag

   do_coeffs       = .true.
   local_do_shadow = do_GlintShadow

!  Initialize Transmittances to 1.0

   trans_norm    = one
   Trans_solar   = one ; dW_Trans_solar   = zero
   Trans_viewing = one ; dW_Trans_viewing = zero
   Trans_stream  = one ; dW_Trans_stream  = zero

!  Set Transmittances flag. Not required fort he Fast calculation option

   do_transmittances = .false.
   if ( .not. Do_Isotropy ) then
      do J = 1, nszas
         if ( Ocean_Reflecs(J).gt.0.0001 ) do_transmittances =.true.
      enddo
   endif

!  Get quadratures

   if ( do_transmittances ) then
      call Water_Transmittance_Quads &
       ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
         PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Output
         AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Output
         TRANS_NORM )
   endif

!   write(*,*) trans_norm ; pause'after quads'

!  Solar angle incident loop

   do J = 1, nszas

!  Skip if no transmittance calculation (values are initialized to 1.0)

      if ( .not. do_transmittances ) go to 67

!  Downward Solar transmittances
!     Set shadow flag for this

      if ( Ocean_Reflecs(J).gt.0.0001 ) then
         phi_w = WindSunAngles(J)
         cphi_w = cos(phi_w*dtr)
         sphi_w = sin(phi_w*dtr)
         local_do_shadow = do_GlintShadow
         call Lin_Water_Transmittance &
         ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
           PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
           AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
           do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
           szas(j), REFRAC_R, REFRAC_I,                            & ! Input
           WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
           PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Inpu
           TRANS_SOLAR(J), dW_TRANS_SOLAR(J) )
      endif

!  Upward transmittances into View directions
!     Passing from water to air, use Snell's Law.  no absorption
!     Local shadow flag turned off here.

      local_do_shadow  = .false.
      local_refrac_R   = one / refrac_R
      local_refrac_I   = zero
      if ( Ocean_Reflecs(J).gt.0.0001 ) then
         do i = 1, nvzas
            incident_angle = asin(sin(vzas(i) * dtr)/refrac_R)/dtr
            call Lin_Water_Transmittance &
            ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
              PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
              AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
              do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
              incident_angle, local_refrac_R, local_refrac_I,         & ! Input
              WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
              PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
              TRANS_VIEWING(I), dW_TRANS_VIEWING(I) )
         enddo
      endif

!  Upward transmittances into stream directions
!     Passing from water to air, use Snell's Law.  no absorption

      local_do_shadow  = .false.
      local_refrac_R   = one / refrac_R
      local_refrac_I   = zero
      if ( Ocean_Reflecs(J).gt.0.0001 ) then
         do i = 1, nstreams
            local_sine = sqrt ( one - streams(i) * streams(i) )
            incident_angle = asin(local_sine/refrac_R)/dtr
            call Lin_Water_Transmittance &
            ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
              PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
              AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
              do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
              incident_angle, local_refrac_R, local_refrac_I,         & ! Input
              WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
              PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
              TRANS_STREAM(I), dW_TRANS_STREAM(I) )
         enddo
      endif

!  continuation point for avoiding transmittance calculations

67    continue

!  Final computation

      Albedo  = 0.485_fpk
      denom   = ( one - Albedo * Ocean_Reflecs(J)  )
      Rwprime =  Ocean_Reflecs(J) / denom

      Const0 = ( one / refrac_sq ) * trans_solar(J) * Rwprime
      Const  = Const0 * Foam_correction
      WLeaving_Iso(J) = Const

      dC_Rwprime = dC_Ocean_Reflecs(J) * ( one + Rwprime * Albedo ) / denom
      dC_const   = Const * dC_Rwprime / Rwprime
      dC_WLeaving_Iso(J) = dC_Const

!      write(*,*)dC_WLeaving_Iso(J),WLeaving_Iso(J)
!      write(*,*)dC_Rwprime, Rwprime
!      write(*,*)dC_Ocean_Reflecs(J), Ocean_Reflecs(J) 

      dW_const0  = Const0 * dW_trans_solar(J) / trans_solar(J)
      dW_const   = Const0 * dW_Foam_correction + dW_Const0 * Foam_correction
      dW_WLeaving_Iso(J) = dW_Const

      do i = 1, nstreams
         WLeaving_SD(J,I)    = Const    * Trans_Stream(I)
         dC_WLeaving_SD(J,I) = dC_Const * Trans_Stream(I)
         dW_WLeaving_SD(J,I) = dW_Const * Trans_Stream(I) + Const * dW_Trans_Stream(I)
      enddo
      do i = 1, nvzas
         WLeaving_SV(J,I)    = Const    * Trans_Viewing(I)
         dC_WLeaving_SV(J,I) = dC_Const * Trans_Viewing(I)
         dW_WLeaving_SV(J,I) = dW_Const * Trans_Viewing(I) + Const * dW_Trans_Viewing(I)
      enddo

      ! Correction of normalisation factor: divide by (pi/mu0). A Sayer 04 Nov 2014.
      ! Have kept this outside the above loops so it is more obvious.
      WLeaving_Iso(J)=WLeaving_Iso(J)/(pi/cos(szas(J)*dtr))
      do i = 1, nstreams
         WLeaving_SD(J,I) = WLeaving_SD(J,I)/(pi/cos(szas(J)*dtr))
      enddo
      do i = 1, nvzas
         WLeaving_SV(J,I) = WLeaving_SV(J,I)/(pi/cos(szas(J)*dtr))
      enddo
      ! And I think we need to do the same for the gradients, too.
      dC_WLeaving_Iso(J)=dC_WLeaving_Iso(J)/(pi/cos(szas(J)*dtr))
      dW_WLeaving_Iso(J)=dW_WLeaving_Iso(J)/(pi/cos(szas(J)*dtr))
      do i = 1, nstreams
         dC_WLeaving_SD(J,I) = dC_WLeaving_SD(J,I)/(pi/cos(szas(J)*dtr))
         dW_WLeaving_SD(J,I) = dW_WLeaving_SD(J,I)/(pi/cos(szas(J)*dtr))
      enddo
      do i = 1, nvzas
         dC_WLeaving_SV(J,I) = dC_WLeaving_SV(J,I)/(pi/cos(szas(J)*dtr))
         dW_WLeaving_SV(J,I) = dW_WLeaving_SV(J,I)/(pi/cos(szas(J)*dtr))
      enddo

!  End solar loop

   enddo

! Finish

   return
end subroutine Linearized_WaterLeaving

!

subroutine Lin_WhiteCap_Reflectance &
    ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian, DWC_Reflectance, DWC_Lambertian )

!  Stand-alone routine for computing the WhiteCap Reflectance
!   Based on 6S code, as updated by A. Sayer (2011)

!  Linearization with respect to Wind-speed

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!    (Wind speed in [m/s], Wavelength in Microns)

   real(fpk), intent(in)  :: WindSpeed
   real(fpk), intent(in)  :: Wavelength

!  output

   real(fpk), intent(out) :: WC_Reflectance
   real(fpk), intent(out) :: WC_Lambertian
   real(fpk), intent(out) :: DWC_Reflectance
   real(fpk), intent(out) :: DWC_Lambertian

!  Data
!  ----

!  Single precision

   real :: Effective_WCRef(39)

! effective reflectance of the whitecaps (Koepke, 1984)
! These are the original values - superseded, A Sayer 05 Jul 2011.
!      data Effective_WCRef/ &
!     0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190,&
!     0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,&
!     0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,&
!     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

! effective reflectance of the whitecaps (Frouin et al, 1996)
! Assume linear trends between the node points they give
! This is the spectral shape

      data Effective_WCRef/ &
     1.000,1.000,1.000,1.000,0.950,0.900,0.700,0.550,0.500,0.450,&
     0.400,0.350,0.300,0.250,0.200,0.150,0.100,0.050,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

!  Local variables
!  ---------------

!  Single precision in the original code

   integer :: iwl, iref
   real    :: Wlb, DWlb, WLP, Ref(39), wspd, wl, Ref_i, Rwc, DRwc

!  Initialize

   WC_Reflectance = 0.0_fpk
   WC_Lambertian  = 0.0_fpk
   DWC_Lambertian = 0.0_fpk

!  Single precision inputs in the original

   wspd = real(WindSpeed)
   wl   = real(Wavelength)

!  Scale data for value of 0.22 in the midvisible.

   DO iref = 1,39
      Ref(iref) = 0.22 * Effective_WCRef(iref)
   ENDDO

!  COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)

   Wlb    = 0.0 ; DWlb = 0.0
   IF (wspd .le. 9.25) THEN
       Wlb  = 0.01*((3.18e-03)*((wspd-3.7)**3.0))
       DWlb = 0.03*((3.18e-03)*((wspd-3.7)**2.0))
   ELSE IF (wspd .gt. 9.25) THEN
       Wlb  = 0.01*((4.82e-04)*((wspd+1.8)**3.0))
       DWlb = 0.03*((4.82e-04)*((wspd+1.8)**2.0))
   END IF

! Original whitecap calculation - superseded, A. Sayer 05 Jul 2011.
!      W=2.95e-06*(wspd**3.52)

!  Find data point, Linearly interpolate

   iwl   = 1+int((wl-0.2)/0.1)
   wlp   = 0.5+(iwl-1)*0.1
   Ref_i = Ref(iwl+1) + ( wl-wlp)/0.1*(Ref(iwl)-Ref(iwl+1))
   Rwc   = Wlb*Ref_i
   DRwc  = DWlb*Ref_i

!  Final values

   WC_Lambertian   = real(Wlb,fpk)
   DWC_Lambertian  = real(DWlb,fpk)
   WC_Reflectance  = real(Rwc,fpk)
   DWC_Reflectance = real(DRwc,fpk)

!  Finish

   return
end subroutine Lin_WhiteCap_Reflectance

subroutine Lin_Ocean_Reflectance &
       ( Maxszas, nszas, SZA_cosines, Wavelength, PigmentConc, noWL, Ocean_Reflecs, DOcean_Reflecs )

!  THIS IS FORMERLY CALLED "MORCASIWAT", as modified by A. Sayer for 6S

! mick fix 12/28/2014 - Using updates by A Sayer November 03 2014:
! Extended functionality down to 200 nm. Achieved by:
! - Extended data arrays down to 200 nm (another 40 elements).
! - Changed logic check for contribution to 0.2-0.9 microns from 0.4-0.9 microns, and started table
!   lookup calculation from 0.2 microns instead of 0.4 microns.
! Note, this is based on a simple extension of the published optical model for vis wavelengths.
!   Possible that other scatterers/absorbers
! which are neglected in this model may be important at UV wavelengths.
! Do linear interpolation of optical property LUTs, rather than nearest neighbour, to remove discontinuities. Achieved by:
! - Replicated final element of LUTs to avoid potential for extrapolation errors.
! - Replace nint() call with floor() call to correctly get lower bound
! - Define variable dwl, fractional distance along the 5 nm LUT grid
! - Implement the interpolation using dwl rather than direct lookup of nearest value.
! Also:
! - Corrected Prieur and Sathyendranath, Limnol. Oceanogr. reference year to 1981 instead of 1983.
! - Corrected typo in water scattering coefficient at 410 nm: 0.0068 was written instead of 0.0061.
!   Removes artificial spike at 410 nm in calculated reflectance.

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input/Output
!  ------------

   integer  , intent(in)   :: Maxszas, nszas
   real(fpk), intent(in)   :: SZA_cosines  (Maxszas)

   real(fpk), intent(in)   :: Wavelength
   real(fpk), intent(in)   :: PigmentConc

   logical  , intent(out)  :: noWL
   real(fpk), intent(out)  :: Ocean_Reflecs(Maxszas)
   real(fpk), intent(out)  :: DOcean_Reflecs(Maxszas) ! Derivatives w.r.t PigmentConc

!  Local
!  -----

!      subroutine morcasiwat(wl,C,R2,mu_sol)
! Rewritten, beginning 07 July 2011, Andrew Sayer
! Now extends underlight calculations out to 400-900 nm.
! and more modern formulations,
! but still based on Case 1 principle.

! Spectral diffuse attenuation coefficient of Case I Waters as Predicted
! by MOREL within the spectral range 400-700nm (1988, Journal of Geophysical
! Research, Vol.93, No C9, pp 10749-10768)
!
! input parameters:	wl wavelength (IN MICROMETERS)
!			C  pigment concentration
!                       mu_sol : cosine of solar zenith angle
! output parameter:	R2  reflectance of water below the surface

! Tabulated absorption coefficient for water, scattering coefficient for water,
! and absorption coefficient for chlorophyll-a, tabulated from 400 nm to 900 nm
! in 5 nm increments
      real water_abs_coeff(142),water_scat_coeff(142),abs_coeff_chl(142)
! Input/output parameters
      real mu_sol,r2,C,wl
! Absorption/scattering terms, and parameters for calculation of f
      real a_wat,b_wat,b_wat_all,a_chl,f,a_tot,b_tot,a_ph,a_cdom,v,bp,bbp,eta,x,z,dwl
! Wavelength index for tables
      integer iwl, J

!  derivatives
      real df,da_tot,db_tot,da_ph,da_cdom,dv,dbp,dbbp,deta,dx,dz,dR2

! Smith and Baker, AO(20) 1981, table 1, out to 800 nm. This has 10 nm increments so
! linearly interpolate between these.
! Hale & Qurry, AO(12) 1973, table 1, for 805-900 nm. This has 25 nm increments so
! linearly interpolate between these. Provided as extinction coefficient
! so convert using a=4*pi*k/lambda (note lambda in m for units m^{-1})
      data water_abs_coeff/ &
       3.0700,2.5300,1.9900,1.6500,1.3100,&
       1.1185,0.9270,0.8235,0.7200,0.6395,&
       0.5590,0.5080,0.4570,0.4150,0.3730,&
       0.3305,0.2880,0.2515,0.2150,0.1780,&
       0.1410,0.1230,0.1050,0.0907,0.0844,&
       0.0761,0.0678,0.0620,0.0561,0.0512,&
       0.0463,0.0421,0.0379,0.0340,0.0300,&
       0.0260,0.0220,0.0206,0.0191,0.0181,&
       0.0171,0.0166,0.0162,0.0158,0.0153,&
       0.0149,0.0144,0.0144,0.0145,0.0145,&
       0.0145,0.0150,0.0156,0.0156,0.0156,&
       0.0166,0.0176,0.0186,0.0196,0.0227,&
       0.0257,0.0307,0.0357,0.0417,0.0477,&
       0.0492,0.0507,0.0532,0.0558,0.0598,&
       0.0638,0.0673,0.0708,0.0753,0.0799,&
       0.0940,0.1080,0.1330,0.1570,0.2005,&
       0.2440,0.2660,0.2890,0.2990,0.3090,&
       0.3145,0.3190,0.3245,0.3290,0.3390,&
       0.3490,0.3740,0.4000,0.4150,0.4300,&
       0.4400,0.4500,0.4750,0.5000,0.5750,&
       0.6500,0.7445,0.8390,1.0040,1.1690,&
       1.4840,1.7990,2.0895,2.3800,2.4250,&
       2.4700,2.5100,2.5500,2.5300,2.5100,&
       2.4350,2.3600,2.2600,2.1600,2.1150,&
       2.0700,2.2104,2.3509,2.4913,2.6318,&
       2.7722,3.0841,3.3960,3.7079,4.0198,&
       4.3317,4.5884,4.8451,5.1019,5.3586,&
       5.6153,5.8495,6.0836,6.3177,6.5518,&
       6.7858,6.7858/

! Smith and Baker, 1981 out to 800 nm.
! This is again at 10 nm increments, interpolated to 5 nm.
! Set to 0.0003 out to 870 nm and 0.0002 after this
! based on Morel's power law (lambda^-4.32)
! and the fact that results not very sensitive to
! the exact value here.
      data water_scat_coeff/&
       0.1510,0.1350,0.1190,0.1093,0.0995,&
       0.0908,0.0820,0.0753,0.0685,0.0630,&
       0.0575,0.0530,0.0485,0.0450,0.0415,&
       0.0384,0.0353,0.0329,0.0305,0.0284,&
       0.0262,0.0246,0.0229,0.0215,0.0200,&
       0.0188,0.0175,0.0164,0.0153,0.0144,&
       0.0134,0.0127,0.0120,0.0113,0.0106,&
       0.0100,0.0094,0.0089,0.0084,0.0080,&
       0.0076,0.0072,0.0068,0.0064,0.0061,&
       0.0058,0.0055,0.0052,0.0049,0.0047,&
       0.0045,0.0043,0.0041,0.0039,0.0037,&
       0.0036,0.0034,0.0033,0.0031,0.0030,&
       0.0029,0.0028,0.0026,0.0025,0.0024,&
       0.0023,0.0022,0.0022,0.0021,0.0020,&
       0.0019,0.0019,0.0018,0.0018,0.0017,&
       0.0017,0.0016,0.0016,0.0015,0.0015,&
       0.0014,0.0014,0.0013,0.0013,0.0012,&
       0.0012,0.0011,0.0011,0.0010,0.0010,&
       0.0010,0.0009,0.0008,0.0008,0.0008,&
       0.0008,0.0007,0.0007,0.0007,0.0007,&
       0.0007,0.0007,0.0007,0.0007,0.0006,&
       0.0006,0.0006,0.0006,0.0006,0.0006,&
       0.0005,0.0005,0.0005,0.0005,0.0005,&
       0.0005,0.0004,0.0004,0.0004,0.0004,&
       0.0004,0.0003,0.0003,0.0003,0.0003,&
       0.0003,0.0003,0.0003,0.0003,0.0003,&
       0.0003,0.0003,0.0003,0.0003,0.0003,&
       0.0002,0.0002,0.0002,0.0002,0.0002,&
       0.0002,0.0002/

! Prieur and Sathyendranath, Limnol. Oceanogr. 26, 1983 table 2
! assumed zero after 700 nm: no data given and water absorption becomes very dominating here
      data abs_coeff_chl/&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.053,0.123,0.195,0.264,&
       0.335,0.405,0.476,0.546,0.617,&
       0.687,0.781,0.828,0.883,0.913,&
       0.939,0.973,1.001,1.000,0.971,&
       0.944,0.928,0.917,0.902,0.870,&
       0.839,0.798,0.773,0.750,0.717,&
       0.668,0.645,0.618,0.582,0.528,&
       0.504,0.474,0.444,0.416,0.384,&
       0.357,0.321,0.294,0.273,0.276,&
       0.268,0.291,0.274,0.282,0.249,&
       0.236,0.279,0.252,0.268,0.276,&
       0.299,0.317,0.333,0.334,0.326,&
       0.356,0.389,0.441,0.534,0.595,&
       0.544,0.502,0.420,0.329,0.262,&
       0.215,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000,0.000,0.000,0.000,&
       0.000,0.000/

!  Initialize

      wl     = real(Wavelength)
      C      = real(PigmentConc)
      noWL = .false.  ; Ocean_Reflecs = 0.0_fpk ; DOcean_Reflecs = 0.0_fpk

! If wavelength out of range, no need to calculate underlight

      if (wl.lt.0.200.or.wl.gt.0.900)then
        noWL = .true. ; goto 60
      endif

! Extract tabulated values for this wavelength
      iwl=1+floor((wl-0.200)/0.005)
! 03 Nov 2014 A Sayer now linear interpolation rather than nearest neighbour
!      a_wat=water_abs_coeff(iwl)
!      b_wat_all=water_scat_coeff(iwl)
!      a_chl=abs_coeff_chl(iwl)
      dwl=(wl-0.200)/0.005-floor((wl-0.200)/0.005)
      a_wat=water_abs_coeff(iwl)+dwl*(water_abs_coeff(iwl+1)-water_abs_coeff(iwl))
      b_wat_all=water_scat_coeff(iwl)+dwl*(water_scat_coeff(iwl+1)-water_scat_coeff(iwl))
      a_chl=abs_coeff_chl(iwl)+dwl*(abs_coeff_chl(iwl+1)-abs_coeff_chl(iwl))

! Morel and Maritorena, 2001
      a_ph  = 0.06*a_chl*(C**0.65)
! Equations 2 and 4a from Morel and Gentili, RSE 113, 2009
      a_cdom  = 0.0524*(C**0.63)*exp(-0.018*(wl*1000.-412.))
      a_tot=a_wat + a_ph + a_cdom

      da_ph = 0.65 * a_ph / C
      da_cdom = 0.63 * a_cdom / C
      da_tot  = da_ph + da_cdom

      b_wat=0.5*b_wat_all

! Morel and Maritorena, 2001 (also earlier work by Loisel and Morel)
! exponent for spectral dependence of scattering

      x=alog10(C) ; dx = 1.0/C/log(10.0) 
      if (C .le. 2) then
        v=0.5*(x-0.3) ; dv = 0.5*dx
      endif
      if (C .gt. 2) then
        v=0. ; dv = 0.
      endif

      bp=0.416*(C**0.766)
      dbp = 0.766*bp/C

      z = (wl/0.55)**v
      dz = z * dv* log(wl/0.55)

      bbp  = 0.002+0.01*(0.5-0.25*x)*z
      dbbp = -0.0025*dx*z+0.01*(0.5-0.25*x)*dz

      b_tot=b_wat + bbp*bp ; db_tot = bbp * dbp + dbbp * bp
      eta=b_wat/b_tot      ; deta = -eta * db_tot / b_tot

! Morel and Gentili, 1991.
     do j = 1, nszas
        mu_sol = real(SZA_cosines(J))
        f=0.6279 - (0.2227*eta) - (0.0513*eta*eta) + (0.2465*eta -0.3119 )*mu_sol
        df = - 0.2227 * deta - 2.0*0.0513*eta*deta + 0.2465*deta *mu_sol
! R=f*b(a+b) : constant multiplied by ratio of backscatter to total extinction + backscatter

        R2=f*b_tot/(a_tot+b_tot)
        dR2 = ( (df*b_tot + f*db_tot) - R2*(da_tot + db_tot) )/(a_tot+b_tot) 

        Ocean_Reflecs(J)  = real(R2,fpk)
        dOcean_Reflecs(J) = real(dR2,fpk)
     enddo

!  continuation point

 60  continue

!  Finish

      return
end subroutine Lin_Ocean_Reflectance


subroutine Lin_Water_Transmittance &
    ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
      PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
      AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
      DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,                     & ! Input
      INCIDENT_ANGLE, REFRAC_R, REFRAC_I,                     & ! Input
      WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
      PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
      TRANS, dTRANS )

      IMPLICIT NONE
      INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)   :: Max_PolarQuads, Max_AzimQuads

!  Quadratures

   real(fpk), intent(in) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk), intent(in) :: CosPolarQuads (Max_PolarQuads)
   real(fpk), intent(in) :: SinPolarQuads (Max_PolarQuads)
   real(fpk), intent(in) :: PolarWeights  (Max_PolarQuads)

   real(fpk), intent(in) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk), intent(in) :: CosAzimQuads (Max_AzimQuads)
   real(fpk), intent(in) :: SinAzimQuads (Max_AzimQuads)
   real(fpk), intent(in) :: AzimWeights  (Max_AzimQuads)

!  Windspeed, coefficients
!  -----------------------

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Windspeed m/s

      REAL(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      REAL(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Cox-Munk Coefficients. Intent(inout).

      REAL(fpk), intent(inout) :: SUNGLINT_COEFFS(7)
      REAL(fpk), intent(inout) :: dSUNGLINT_COEFFS(7)

!  Other inputs
!  ------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  incident angle in degrees

      REAL(fpk), intent(in)    :: INCIDENT_ANGLE

!  Real and imaginary parts of refractive index

      REAL(fpk), intent(in)    :: REFRAC_R
      REAL(fpk), intent(in)    :: REFRAC_I

!  Pre-computed Norm

      REAL(fpk), intent(in)    :: TRANS_NORM

!  Output
!  ======

      REAL(fpk), intent(out)   :: TRANS
      REAL(fpk), intent(out)   :: dTRANS

!  Local
!  =====

      integer   :: i, k
      real(fpk) :: dtr, xj, sxj, xi, sxi, phi, cphi, sphi, weight
      real(fpk) :: SUNGLINT_REFLEC, dSUNGLINT_REFLEC

!  Computation setup

      TRANS  = 0.0_fpk
      dTRANS = 0.0_fpk
      DTR   = ACOS(-1.0d0) / 180.0_fpk
      XJ  = COS ( INCIDENT_ANGLE * DTR )
      SXJ = SQRT ( 1.0_fpk - XJ * XJ )

!  Loops

      do k = 1, Max_AzimQuads
         PHI  = AzimQuads(K)/dtr
         CPHI = CosAzimQuads(K)
         SPHI = SinAzimQuads(K)
         do i = 1, Max_PolarQuads
            XI  = CosPolarQuads(I)
            SXI = SinPolarQuads(I)
            Call Lin_GENERAL_SUNGLINT &
             ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
               REFRAC_R, REFRAC_I, WINDSPEED,         &
               PHI_W, CPHI_W, SPHI_W,                 &
               XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
               SUNGLINT_COEFFS, DSUNGLINT_COEFFS,     &
               SUNGLINT_REFLEC, DSUNGLINT_REFLEC )
            WEIGHT = PolarWeights(I) * AzimWeights(k)
            TRANS  = TRANS  + SUNGLINT_REFLEC  * WEIGHT
            dTRANS = dTRANS + dSUNGLINT_REFLEC * WEIGHT
         enddo
      enddo
      dTRANS = -dTRANS/TRANS_NORM
      TRANS = 1.0_fpk - (TRANS/TRANS_NORM)

!  done

      RETURN
end subroutine Lin_Water_Transmittance

!

SUBROUTINE Lin_GENERAL_SUNGLINT &
         ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
           REFRAC_R, REFRAC_I, WINDSPEED,         &
           PHI_W, CPHI_W, SPHI_W,                 &
           XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
           SUNGLINT_COEFFS, DSUNGLINT_COEFFS,     &
           SUNGLINT_REFLEC, DSUNGLINT_REFLEC )

      implicit none

      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Subroutine Input arguments
!  --------------------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Real and imaginary parts of refractive index

      REAL(fpk), intent(in)    :: REFRAC_R
      REAL(fpk), intent(in)    :: REFRAC_I

!  Windspeed m/s

      REAL(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      REAL(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Incident and reflected ddirections: sines/cosines. Relative azimuth (angle in radians)

      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SPHI

!  Subroutine output arguments
!  ---------------------------

!   Glitter reflectance

      REAL(fpk), intent(out)   :: SUNGLINT_REFLEC
      REAL(fpk), intent(out)   :: dSUNGLINT_REFLEC

!  Cox-Munk Coefficients. Intent(inout).

      REAL(fpk), intent(inout) :: SUNGLINT_COEFFS(7)
      REAL(fpk), intent(inout) :: DSUNGLINT_COEFFS(7)

!  Local arguments
!  ---------------

!  parameters from LIDORT

   real(fpk), PARAMETER :: ONE = 1.0_fpk, ZERO  = 0.0_fpk, ONEP5 = 1.5_fpk
   real(fpk), PARAMETER :: TWO = 2.0_fpk, THREE = 3.0_fpk, FOUR  = 4.0_fpk
   real(fpk), PARAMETER :: six = two * three, twentyfour = six * four
   real(fpk), PARAMETER :: QUARTER = 0.25_fpk, HALF = 0.5_fpk

   real(fpk), PARAMETER :: MINUS_ONE = - ONE
   real(fpk), PARAMETER :: MINUS_TWO = - TWO

   real(fpk), PARAMETER :: PIE = ACOS(MINUS_ONE)
   real(fpk), PARAMETER :: DEG_TO_RAD = PIE/180.0_fpk

   real(fpk), PARAMETER :: PI2  = TWO  * PIE
   real(fpk), PARAMETER :: PI4  = FOUR * PIE
   real(fpk), PARAMETER :: PIO2 = HALF * PIE
   real(fpk), PARAMETER :: PIO4 = QUARTER * PIE

   REAL(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      REAL(fpk)  :: B, ZX, ZY, Z, Z1, Z2, XMP
      REAL(fpk)  :: TILT, TANTILT, TANTILT_SQ, COSTILT
      REAL(fpk)  :: ARGUMENT, PROB, FAC2, COEFF, VAR, WSigC, WSigU
      REAL(fpk)  :: XE, XN, XE_sq, XN_sq, XE_sq_1, XN_sq_1
      REAL(fpk)  :: XPHI, CKPHI, SKPHI, XPHI_W, CKPHI_W, SKPHI_W
      REAL(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(fpk)  :: SHADOWI, SHADOWR, SHADOW

      REAL(fpk)  :: dARGUMENT, dPROB, dCOEFF, dVAR, dWSigC, dWSigU
      REAL(fpk)  :: AN, AE, dXE, dXN, dXE_sq, dXN_sq, EXPO, dEXPO
      REAL(fpk)  :: T3, T4, T5, T6, T7, dT3, dT4, dT5, dT6, dT7

!  Initialise output

      SUNGLINT_REFLEC = ZERO
      dSUNGLINT_REFLEC = ZERO

!  COmpute coefficients, according to 6S formulation

      IF ( DO_COEFFS ) THEN
         SUNGLINT_COEFFS = zero ; DSUNGLINT_COEFFS = zero
         IF ( DO_ISOTROPIC ) THEN
            SUNGLINT_COEFFS(1)  = 0.003_fpk + 0.00512_fpk * WINDSPEED
            DSUNGLINT_COEFFS(1) = 0.00512_fpk
         ELSE
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00192_fpk * WINDSPEED ! sigmaC
            SUNGLINT_COEFFS(2) =             0.00316_fpk * WINDSPEED ! sigmaU
            SUNGLINT_COEFFS(3) = 0.010_fpk - 0.00860_fpk * WINDSPEED ! C21
            SUNGLINT_COEFFS(4) = 0.040_fpk - 0.03300_fpk * WINDSPEED ! C03
            SUNGLINT_COEFFS(5) = 0.400_fpk                           ! C40
            SUNGLINT_COEFFS(6) = 0.230_fpk                           ! C04
            SUNGLINT_COEFFS(7) = 0.120_fpk                           ! C22
            DSUNGLINT_COEFFS(1) = 0.00192_fpk
            DSUNGLINT_COEFFS(2) = 0.00316_fpk 
            DSUNGLINT_COEFFS(3) = - 0.00860_fpk
            DSUNGLINT_COEFFS(4) = - 0.03300_fpk
         ENDIF
         DO_COEFFS = .false.
      ENDIF

!  Local angles

      XPHI   = PIE - PHI       ! Not used
!     CKPHI  = - CPHI          ! Original, not correct.

      CKPHI  = + CPHI
      SKPHI  = + SPHI

      XPHI_W  = PHI_W
      CKPHI_W = CPHI_W
      SKPHI_W = SPHI_W

!  Tilt angle

      B  = ONE / ( XI + XJ )
      ZX = - SXI * SKPHI * B
      ZY = ( SXJ + SXI * CKPHI ) * B
      TANTILT_SQ = ZX * ZX + ZY * ZY
      TANTILT    = SQRT ( TANTILT_SQ )
      TILT       = ATAN(TANTILT)
      COSTILT    = COS(TILT)

!  Scatter angle

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)

!  Fresnel
!  -------

       CALL Fresnel_Complex ( REFRAC_R, REFRAC_I, Z2, XMP )

!  Anisotropic
!  -----------

      IF ( .not. DO_ISOTROPIC ) THEN

!  Variance

         WSigC = ONE / Sqrt(SUNGLINT_COEFFS(1)) ; dWSigC = - half * WSigC * WSigC * WSigC * DSUNGLINT_COEFFS(1)
         WSigU = ONE / Sqrt(SUNGLINT_COEFFS(2)) ; dWSigU = - half * WSigU * WSigU * WSigU * DSUNGLINT_COEFFS(2)
         VAR   = WSigC * WSigU * HALF ; dVAR = half * ( dWSigC * WSigU + WSigC * dWSigU )
         VAR = ONE / VAR ; dVAR =  - VAR * VAR * dVAR

!  angles

         AE = (  CKPHI_W * ZX + SKPHI_W * ZY )
         AN = ( -SKPHI_W * ZX + CKPHI_W * ZY )
         XE = AE * WSigC ; XE_sq = XE * XE ; XE_sq_1 = xe_sq - one ; dXE = AE * dWSigC ; dXE_sq = two * dXE * XE
         XN = AN * WSigU ; XN_sq = XN * XN ; XN_sq_1 = xn_sq - one ; dXN = AN * dWSigU ; dXN_sq = two * dXN * XN

!  GC Coefficient

         T3  = XE_sq_1 * XN * half
         dT3 = ( XE_sq_1 * dXN + dXE_sq * XN ) * half
         T4  = ( XN_sq - three ) * XN / six
         dT4 = ( ( XN_sq - three ) * dXN + dXN_sq * XN ) / six
         T5  = ( XE_sq * XE_sq - six * XE_sq + three ) / twentyfour
         dT5 = ( two * dXE_sq * XE_sq - six * dXE_sq ) / twentyfour
         T6  = ( XN_sq * XN_sq - six * XN_sq + three ) / twentyfour
         dT6 = ( two * dXN_sq * XN_sq - six * dXN_sq ) / twentyfour
         T7  = XE_sq_1 * XN_sq_1 / four
         dT7 = ( dXE_sq * XN_sq_1 + XE_sq_1 * dXN_sq ) / four

         Coeff  = ONE - SUNGLINT_COEFFS(3) * T3 &
                      - SUNGLINT_COEFFS(4) * T4 &
                      + SUNGLINT_COEFFS(5) * T5 &
                      + SUNGLINT_COEFFS(6) * T6 &
                      + SUNGLINT_COEFFS(7) * T7
         dCoeff =  - dSUNGLINT_COEFFS(3) * T3 - SUNGLINT_COEFFS(3) * dT3 &
                   - dSUNGLINT_COEFFS(4) * T4 - SUNGLINT_COEFFS(4) * dT4 &
                                              + SUNGLINT_COEFFS(5) * dT5 &
                                              + SUNGLINT_COEFFS(6) * dT6 &
                                              + SUNGLINT_COEFFS(7) * dT7

!  Probability and finish

         ARGUMENT  = (  XE_sq  +  XN_sq ) * HALF
         dARGUMENT = ( dXE_sq  + dXN_sq ) * HALF
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            EXPO = EXP ( - ARGUMENT ) ; dEXPO = - dARGUMENT * EXPO
            PROB = COEFF * EXPO / VAR ; dPROB =  ( dCOEFF * EXPO + COEFF * dEXPO - PROB * dVAR ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC  = XMP * PROB  * FAC2
            dSUNGLINT_REFLEC = XMP * dPROB * FAC2
         ENDIF

      ENDIF

!  Isotropic
!  ---------

      IF ( DO_ISOTROPIC ) THEN

!  Compute Probability and finish

         VAR   = SUNGLINT_COEFFS(1) ; dVAR = dSUNGLINT_COEFFS(1)
         ARGUMENT = TANTILT_SQ / VAR
         dARGUMENT = - ARGUMENT * dVAR / VAR
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            EXPO = EXP ( - ARGUMENT ) ; dEXPO = - dARGUMENT * EXPO
            PROB = EXPO / VAR ; dPROB =  ( dEXPO - PROB * dVAR ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC  = XMP * PROB  * FAC2
            dSUNGLINT_REFLEC = XMP * dPROB * FAC2
         ENDIF

      ENDIF

!  No Shadow code if not flagged

      IF ( .not. DO_SHADOW  ) RETURN

!  Shadow code

      S1 = SQRT ( VAR / PIE )
      S3 = ONE / ( SQRT(VAR) )
      S2 = S3 * S3

      XXI  = XI*XI
      DCOT = XI / SQRT ( ONE - XXI )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function
      SHADOWI = HALF * ( S1 * T1 / DCOT - T2 )

      XXJ  = XJ*XJ
      DCOT = XJ / SQRT ( ONE - XXJ )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function
      SHADOWR = HALF * ( S1 * T1 / DCOT - T2 )

      SHADOW = ONE / ( ONE + SHADOWI + SHADOWR )
      SUNGLINT_REFLEC = SUNGLINT_REFLEC * SHADOW

!     Finish

      RETURN
END SUBROUTINE Lin_GENERAL_SUNGLINT

!  End module

END MODULE sleave_lin_sup_routines_m

