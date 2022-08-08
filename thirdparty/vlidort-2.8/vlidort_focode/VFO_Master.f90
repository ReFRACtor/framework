
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

!  FO Version history
!  ------------------

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional phase matrices, Surface Leaving and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional F-matrix usage
!    Version 5,  02 August   2016. Sleave + Jacobians
!    Version 5,  19 August   2016, Partial-layer output
!    Version 5,  11 December 2017, Optional code for polarized emissivity

!  VLIDORT Interface history
!  -------------------------

!    FO Version 1.4: This module is interface with VLIDORT V2.7. R.Spurr 3/19/15
!    FO Version 1.5: Interface module upgraded to  VLIDORT V2.8. R.Spurr 7/7/16, 9/17/16

!  1/31/21. Version 2.8.3. Following Direct-Thermal upgrades made for 2.8.1, 5/5/20.
!    ==> Thermal geometry alls now introduced inside the downwelling clause (formerly absent)
!    ==> Add radiifine_LOS + radiifine_p_LOS to the argument list
!    ==> These are vertical distances from layer top, needed for FO Outgoing direct-thermal calculation
!    ==> Geometrical calculation for direct thermal outgoing downwelling and upwelling was corrected.

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP/DN to the output lists from SSV_Integral_I_UP/DN, SSV_Integral_I_UPDN
!    ==> Add Theta_all output from Geometry routine for Solar sources 

!  1/31/21. Version 2.8.3. Other upgrades made for 2.8.1, 5/5/20.
!    ==> Geometry routine for Solar sources now has Doublet option.
!    ==> doublet and lattice offsets defined in main subroutine and passed down when needed
!    ==> Water-leaving cumulative tranmsittance was properly initialized.

module VFO_Master_m

!  All subroutines public

public

contains

subroutine VFO_MASTER &
       ( maxgeoms, maxszas, maxvzas, maxazms,                                        & ! Input max dims
         maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,         & ! Input max dims
         do_solar_sources, do_sunlight, do_thermal_emission, do_surface_emission,    & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_fmatrix, do_obsgeom, do_doublet, do_deltam,  & ! Input flags (general)
         do_lambertian, do_surface_leaving, do_water_leaving, do_Polarized_Emiss,    & ! Input flags (surface)
         do_Partials, do_planpar, do_enhanced_ps,                                    & ! Input flags (geoms)
         nstokes, ngeoms, nszas, nvzas, nazms, nlayers, nfine, nmoments_input,       & ! Inputs (control-numbers)
         nd_offset, nv_offset, na_offset, n_user_levels, LevelMask_up, LevelMask_dn, & ! Inputs (offsets and control-levels)
         npartials, partial_outindex, partial_outflag, partial_layeridx,             & ! Inputs (control-partial)
         dtr, Pie, doCrit, Acrit, eradius, heights, partial_heights,                 & ! Input general
         obsgeom_boa, alpha_boa, theta_boa, phi_boa, flux, fluxvec,                  & ! Input geometry/Flux
         extinction, deltaus, omega, greekmat, fmatrix_up, fmatrix_dn,               & ! Input atmos optical
         truncfac, bb_input, surfbb, emiss, reflec, slterm,                          & ! Input thermal/surf optical
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,                   & ! Output
         fo_stokes_atmos, fo_stokes_surf, fo_stokes,                                 & ! Main FO   Output
         cumtrans, lostrans_up, lostrans_dn, theta_all, alpha,                       & ! Auxiliary Output
         Master_fail, message, trace_1, trace_2 )                                      ! Exception handling

!  4/9/19.  Add CUMTRANS output, and Waterleaving input control

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add lostrans_up, Lostrans_dn, theta_all and alpha to the output. 

!  1/31/21. Version 2.8.3. 
!    -- Add Doublet flag to input list, define offsets locally in main routine
!    -- Set do_Polarized_Emiss in main routine (input here)

!  Use modules

   USE FO_SSWPGeometry_Master_m
   USE FO_DTWPGeometry_Master_m

   USE FO_VectorSS_spherfuncs_m
   USE FO_VectorSS_RTCalcs_I_m

   USE FO_Thermal_RTCalcs_I_m

   implicit none

!  parameter arguments

   integer, parameter :: ffp = selected_real_kind(15),&
                         maxstokes = 4, max_directions = 2,&
                         upidx = 1, dnidx = 2

!  Subroutine inputs
!  =================

!  Max dimensions
!  --------------

   integer, intent(in)  :: maxgeoms, maxszas, maxvzas, maxazms
   integer, intent(in)  :: maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels

!  Configuration inputs
!  --------------------

!  Sources control, including thermal and Vector sunlight flag

   logical, intent(in)  :: do_solar_sources
   logical, intent(in)  :: do_thermal_emission
   logical, intent(in)  :: do_surface_emission
   logical, intent(in)  :: do_sunlight

!  Directional Flags

   logical, intent(in)  :: do_upwelling, do_dnwelling

!  deltam scaling flag

   logical, intent(in)  :: do_deltam

!  F-matrix flag. (FO 1.5, VLIDORT 2.8). Introduced 7/7/16
!    If set, FO will use F-matrix input directly 

   logical, intent(in)  :: do_fmatrix

!  Obsgeom flag
!  1/31/21. Version 2.8.3. Add Doublet flag

   logical, intent(in)  :: do_Obsgeom
   logical, intent(in)  :: do_Doublet

!  Lambertian surface flag. Sleave flag added, 8/2/16, water-leaving 4/9/19
!    -- 1/31/21. Version 2.8.3. Polarized Emissivity now an input

   logical, intent(in)  :: do_lambertian
   logical, intent(in)  :: do_surface_leaving
   logical, intent(in)  :: do_water_leaving
   logical, intent(in)  :: do_Polarized_Emiss

!  flags. Version 1.5:  Partials 8/25/16

   logical, intent(in)  :: do_Partials

!  Flags (sphericity flags are mutually exclusive). Regular PS now removed, version 1.5

   logical, intent(in)  :: do_planpar
   logical, intent(in)  :: do_enhanced_ps

!  Numbers
!  -------

!  Number of Stokes components

   integer, intent(in)    :: nstokes

!  Layer and geometry control. Finelayer divisions may be changed

   integer, intent(in)    :: ngeoms, nszas, nvzas, nazms, nlayers, nfine
   integer, intent(in)    :: nmoments_input

!  1/31/21. Version 2.8.3. Offsets now inputs from main calling routine

   integer, intent(in)    :: nd_offset(maxszas), nv_offset(maxszas), na_offset(maxszas,maxvzas)

!  Output levels. Use masking as in main codes, 9/17/16

   integer, intent(in) :: n_user_levels
   integer, intent(in) :: LevelMask_up ( max_user_levels )
   integer, intent(in) :: LevelMask_dn ( max_user_levels )

!  Control for partial-layer output, added 8/25/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( max_user_levels )
   integer, Intent(in) :: partial_outindex( max_user_levels )

!  General inputs
!  --------------

!  dtr = degrees-to-Radians. Pie = 3.14159...

   real(ffp), intent(in) :: dtr, Pie

!  Critical adjustment for cloud layers. Not enabled. 9/17/16

   logical, intent(inout)  :: doCrit
   real(ffp),   intent(in) :: Acrit

!  Earth radius + heights. Partials added 9/17/16.

   real(ffp), intent(in)   :: eradius, heights (0:maxlayers)
   real(ffp), intent(In)   :: partial_heights (maxpartials)

!  Geometry inputs
!  ---------------

!  input angles (Degrees). Enough information for Lattice or Obsgeom.
!   Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)
!    In both cases, the Phi angle may be changed.....

   real(ffp), intent(inout)  :: Obsgeom_boa(maxgeoms,3)
   real(ffp), intent(inout)  :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  optical inputs
!  --------------

!  Solar flux

   real(ffp), intent(in) :: flux, fluxvec(maxstokes)

!  Atmosphere

   real(ffp), intent(in) :: extinction  ( maxlayers )
   real(ffp), intent(in) :: deltaus     ( maxlayers )
   real(ffp), intent(in) :: omega       ( maxlayers )
   real(ffp), intent(in) :: greekmat    ( maxlayers, 0:maxmoments_input, maxstokes, maxstokes )

!  Fmatrix input. (FO 1.5, VLIDORT 2.8). Introduced 7/7/16

   real(ffp), intent(in) :: fmatrix_up  ( maxlayers, maxgeoms, 6 )
   real(ffp), intent(in) :: fmatrix_dn  ( maxlayers, maxgeoms, 6 )

!  For TMS correction

   real(ffp), intent(in) :: truncfac ( maxlayers )

!  Thermal inputs, surface emissivity.
!   Emissivity is now polarized. 12/11/17 Rob add.

   real(ffp), intent(in) :: bb_input ( 0:maxlayers )
   real(ffp), intent(in) :: surfbb
!   real(ffp), intent(in) :: emiss ( maxvzas )
   real(ffp), intent(in) :: emiss ( maxstokes, maxvzas )

!  Surface properties - reflective (could be the albedo), surface leaving added 8/2/16

   real(ffp), intent(in) :: reflec(maxstokes,maxstokes,maxgeoms)
   real(ffp), intent(in) :: slterm(maxstokes,maxgeoms)

!  Subroutine outputs
!  ==================

!  Solar

   real(ffp), intent(out) :: fo_stokes_ss ( max_user_levels,maxgeoms,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_db ( max_user_levels,maxgeoms,maxstokes )

!  Thermal

   real(ffp), intent(out) :: fo_stokes_dta ( max_user_levels,maxgeoms,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_dts ( max_user_levels,maxgeoms,maxstokes )

!  Composite
!mick mod 9/19/2017 - added "fo_stokes_atmos" & "fo_stokes_surf" to vector FO output to
!                     keep FO scalar & vector subroutines on an equal footing
                     
   real(ffp), intent(out) :: fo_stokes_atmos ( max_user_levels,maxgeoms,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_surf  ( max_user_levels,maxgeoms,maxstokes )
   real(ffp), intent(out) :: fo_stokes       ( max_user_levels,maxgeoms,maxstokes,max_directions )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out) :: CUMTRANS ( max_user_levels, maxgeoms )

!  1/31/21. Version 2.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add lostrans_up, Lostrans_dn, theta_all and alpha to the output. 

   real(ffp), Intent(out) :: lostrans_up (maxgeoms,maxlayers)
   real(ffp), Intent(out) :: lostrans_dn (maxgeoms,maxlayers)

   real(ffp), Intent(out) :: theta_all ( 0:maxlayers, maxgeoms )
   real(ffp), Intent(out) :: alpha     ( 0:maxlayers, maxgeoms )

!  Exception handling

   logical, intent(out)           :: Master_fail
   character (len=*), intent(out) :: message
   character (len=*), intent(out) :: trace_1, trace_2

!  Other variables
!  ===============

!  Geometry routine outputs
!  ------------------------

!  VSIGN = +1 (Up); -1(Down)

   real(ffp)  :: vsign

!  LOSPATHS flag has been removed now.  8/1/13

!  Flag for the Nadir case

   logical    :: doNadir(maxgeoms)
  
!  cosa, sina, Radii, Ray constant. 
!    -- 1/31/21. Version 2.8.3. alpha removed from here, now an argument

   real(ffp)  :: radii    (0:maxlayers)
   real(ffp)  :: Raycon   (maxgeoms)
   real(ffp)  :: cosa     (0:maxlayers,maxgeoms)
   real(ffp)  :: sina     (0:maxlayers,maxgeoms)

   real(ffp)  :: radii_p    (maxpartials)
   real(ffp)  :: alpha_p    (maxpartials,maxgeoms)
   real(ffp)  :: cosa_p     (maxpartials,maxgeoms)
   real(ffp)  :: sina_p     (maxpartials,maxgeoms)

!  Critical layer. Not yet active 9/17/16.
!   integer    :: Ncrit(maxgeoms)
!   real(ffp)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Existence flags. 8/19/16. Criticality enters here

   logical    :: do_sources_up       (maxlayers,maxgeoms)
   logical    :: do_sources_dn       (maxlayers,maxgeoms)
   logical    :: do_sources_up_p     (maxpartials,maxgeoms)
   logical    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  1/31/21. Version 2.8.3. Separate existence flags for the Thermal

   logical    :: do_Tsources_up       (maxlayers,maxvzas)
   logical    :: do_Tsources_dn       (maxlayers,maxvzas)
   logical    :: do_Tsources_up_p     (maxpartials,maxvzas)
   logical    :: do_Tsources_dn_p     (maxpartials,maxvzas)

!  Chapman factors

   real(ffp)  :: Chapfacs      (maxlayers,  maxlayers,maxgeoms)
   real(ffp)  :: chapfacs_p    (maxpartials,maxlayers,maxgeoms)

!  Los paths added, 8/17/16

   real(ffp)  :: LosW_paths(maxlayers  ,maxgeoms)
   real(ffp)  :: LosP_paths(maxpartials,maxgeoms)

!  Cosine scattering angle, other cosines

   real(ffp)  :: cosscat (maxgeoms)
   real(ffp)  :: Mu0     (maxgeoms)
   real(ffp)  :: Mu1     (maxgeoms)

!  LOS Quadratures for Enhanced PS

   integer    :: nfinedivs (maxlayers,maxgeoms)
   real(ffp)  :: xfine     (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: radiifine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: wfine     (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: sinfine   (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: cosfine   (maxlayers,maxfine,maxgeoms)

!  Quadratures for partial 

   integer    :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp)  :: xfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: wfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: radiifine_p (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: sinfine_p   (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: cosfine_p   (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/17/16.

   integer    :: ntraverse     (0:maxlayers,maxgeoms)
   real(ffp)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)
   integer    :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer    :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer    :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  Spherfunc routine outputs
!  -------------------------

!  Spherical functions, rotation angles

   real(ffp)  :: rotations_up(4,maxgeoms)
   real(ffp)  :: genspher_up(0:maxmoments_input,4,maxgeoms)
   real(ffp)  :: rotations_dn(4,maxgeoms)
   real(ffp)  :: genspher_dn(0:maxmoments_input,4,maxgeoms)
   real(ffp)  :: gshelp(7,0:maxmoments_input)

!  RT calculation outputs
!  ----------------------

!  Solar routines

   real(ffp)  :: stokes_up    ( max_user_levels,maxstokes,maxgeoms )
   real(ffp)  :: stokes_dn    ( max_user_levels,maxstokes,maxgeoms )
   real(ffp)  :: stokes_db    ( max_user_levels,maxstokes,maxgeoms )

!  Thermal routines (Scalar, no polarization). SEE LOS VARIABLES (THERMAL), below
!   real(ffp)  :: intensity_dta_up ( max_user_levels,maxgeoms )
!   real(ffp)  :: intensity_dta_dn ( max_user_levels,maxgeoms )
!   real(ffp)  :: intensity_dts    ( max_user_levels,maxgeoms )

!  Intermediate RT products
!  ------------------------

!  Composite
!mick mod 9/19/2017 - added "fo_stokes_atmos" & "fo_stokes_surf" to vector FO output to
!                     keep FO scalar & vector subroutines on an equal footing

   !real(ffp)  :: fo_stokes_atmos ( max_user_levels,maxgeoms,maxstokes,max_directions )
   !real(ffp)  :: fo_stokes_surf  ( max_user_levels,maxgeoms,maxstokes )

!  LOS VARIABLES (THERMAL SOLUTION)
!  --------------------------------

!  Output

   real(ffp)  :: intensity_dta_up_LOS ( max_user_levels,maxvzas )
   real(ffp)  :: intensity_dta_dn_LOS ( max_user_levels,maxvzas )
   real(ffp)  :: intensity_dts_LOS    ( max_user_levels,maxvzas )

   real(ffp)  :: cumsource_up_LOS     ( 0:maxlayers,maxvzas )
   real(ffp)  :: cumsource_dn_LOS     ( 0:maxlayers,maxvzas )

!  1/31/21. Version 2.8.3. Following earlier upgrade.
!    ==> LOSTRANS output might be used again for the (upwelling) thermal-NoScattering 

   real(ffp)  :: lostrans_up_LOS      ( maxlayers  ,maxvzas )
   real(ffp)  :: lostrans_up_p_LOS    ( maxpartials,maxvzas )

!  StokesQUV_dts_LOS, contribution from Polarized emissivity. 12/11/17 Rob add.

   real(ffp)  :: StokesQUV_dts_LOS    ( max_user_levels, 3, maxvzas )

!  Geometry. Los paths added, 8/25/16. Partials added 8/22/16

   real(ffp)  :: Mu1_LOS(maxvzas)

   real(ffp)  :: alpha_LOS    (0:maxlayers,maxvzas)
   real(ffp)  :: cosa_LOS     (0:maxlayers,maxvzas)
   real(ffp)  :: sina_LOS     (0:maxlayers,maxvzas)

   real(ffp)  :: alpha_p_LOS    (maxpartials,maxvzas)
   real(ffp)  :: cosa_p_LOS     (maxpartials,maxvzas)
   real(ffp)  :: sina_p_LOS     (maxpartials,maxvzas)

   real(ffp)  :: LosW_paths_LOS (maxlayers,maxvzas)
   real(ffp)  :: LosP_paths_LOS (maxpartials,maxvzas)

!  LOS Quadratures for Enhanced PS. Partials added 8/25/16.

   integer    :: nfinedivs_LOS (maxlayers,maxvzas)
   real(ffp)  :: xfine_LOS     (maxlayers,maxfine,maxvzas)
   real(ffp)  :: wfine_LOS     (maxlayers,maxfine,maxvzas)
   real(ffp)  :: cosfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: sinfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: alphafine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: radiifine_LOS (maxlayers,maxfine,maxvzas)

   integer    :: nfinedivs_p_LOS (maxpartials,maxvzas)
   real(ffp)  :: xfine_p_LOS     (maxpartials,maxfine,maxvzas)
   real(ffp)  :: wfine_p_LOS     (maxpartials,maxfine,maxvzas)
   real(ffp)  :: cosfine_p_LOS   (maxpartials,maxfine,maxvzas)
   real(ffp)  :: sinfine_p_LOS   (maxpartials,maxfine,maxvzas)
   real(ffp)  :: alphafine_p_LOS (maxpartials,maxfine,maxvzas)
   real(ffp)  :: radiifine_p_LOS (maxpartials,maxfine,maxvzas)

!  No criticality yet. 9/17/16
!   integer    :: Ncrit_LOS(maxvzas)
!   real(ffp)  :: RadCrit_LOS(maxvzas), CotCrit_LOS(maxvzas)

!  Other products
!  --------------

!  Albedo. 1/31/21. Version 2.8.3. No longer required, commented out
!   real(ffp) :: albedo

!  Thermal setup

   real(ffp)  :: tcom1(maxlayers,2)

!  Dummies

   real(ffp)  :: SScumsource_up ( 0:maxlayers,maxstokes,maxgeoms )
   real(ffp)  :: SScumsource_dn ( 0:maxlayers,maxstokes,maxgeoms )

!   real(ffp)  :: DTcumsource_up ( 0:maxlayers,maxgeoms )
!   real(ffp)  :: DTcumsource_dn ( 0:maxlayers,maxgeoms )

!  LOCAL HELP VARIABLES
!  --------------------

!  numbers

   real(ffp), parameter :: zero = 0.0_ffp

!   1/31/21. Version 2.8.3. Offsets/Polarized Emissivity removed from here, now inputs

   integer   :: ns, nv, na, g, lev, o1
   logical   :: STARTER, do_Thermset, fail, do_Chapman, Do_spherfunc

!  Initialize output
!mick mod 9/19/2017 - initialized "fo_stokes_atmos" & "fo_stokes_surf"

   fo_stokes_ss    = zero
   fo_stokes_db    = zero
   fo_stokes_dta   = zero
   fo_stokes_dts   = zero
   fo_stokes_atmos = zero
   fo_stokes_surf  = zero
   fo_stokes       = zero

   Master_fail = .false.
   message = ' '
   trace_1 = ' '
   trace_2 = ' '

!  Flags to be set for each calculation (safety)

   do_Chapman  = .false.
   do_Thermset = .true.
   starter     = .true.
!mick fix 3/25/2015 - added initialization
   doNadir     = .false.

!  No need to calculate spherical function if using F-matrix input
!    Turn off the local "SPHERFUNC" flag in this case

   Do_spherfunc = .not. do_fmatrix

!  Offsets. -- 1/31/21. Version 2.8.3. Removed from here. Now inputs
!   na_offset = 0 ; nv_offset = 0
!   if ( .not. do_obsgeom ) then
!     if ( do_doublet ) then
!       do ns = 1, nszas
!         nd_offset(ns) = nvzas * (ns - 1) 
!       enddo
!     else
!       do ns = 1, nszas ;  do nv = 1, nvzas
!           na_offset(ns,nv) = nvzas * nazms * (ns - 1) + nazms * (nv - 1)
!       enddo ; enddo
!     endif
!   endif

!  Set do_Polarized_Emiss flag. 12/11/17 Rob add.
!    -- 1/31/21. Now an input, flag is set in the calling routine
!   do_Polarized_Emiss = .false.
!   if ( nstokes.gt.1 .and. do_surface_emission ) then
!      do nv = 1, nvzas
!        if ( SUM(emiss(2:4,nv)).ne.zero ) do_Polarized_Emiss = .true.
!      enddo
!   endif

!  temporary fix until criticality realized. 9/17/16
!  -------------------------------------------------

!  1/31/21. These apply only to the solar source terms

   do_sources_up   = .true.
   do_sources_dn   = .true.
   do_sources_up_p = .true.
   do_sources_dn_p = .true.

!  Solar sources run (NO THERMAL)
!  ------------------------------

   if ( do_solar_sources ) then

!  Upwelling

     if ( do_upwelling ) then
       vsign =  1.0_ffp

!  Geometry call. Updated 9/17/16.

!    -- 1/31/21. Version 2.8.3. theta_all added to the output, needed for MSST output later on.
!    -- 1/31/21. Version 2.8.3. do_doublet flag added, Offsets added, some inputs rearranged

       call FO_SSWPGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine,          & ! Input dimensions/constants
         do_obsgeom, do_doublet, do_Chapman, do_planpar, do_enhanced_ps, do_Partials,   & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, npartials, partial_layeridx,      & ! Input numbers             
         dtr, Pie, vsign, eradius, nv_offset, na_offset, nd_offset,                     & ! Input Constants/Bookkeeping
         heights, partial_heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,          & ! Input heights/geometry
         doNadir, Raycon, Mu0, Mu1, cosscat,                                                & ! Outputs geometry
         Radii,   LosW_paths, alpha, sina, cosa, sunpaths, ntraverse, chapfacs, theta_all,  & ! Outputs (Layer boundaries)
         Radii_p, LosP_paths, alpha_p, sina_p, cosa_p, sunpaths_p, ntraverse_p, chapfacs_p, & ! Outputs (partial levels)
         nfinedivs,   xfine,   wfine,   radiifine,   alphafine,                             & ! Output Wholelayer
         sinfine,   cosfine,   sunpathsfine,   ntraversefine,                               & ! Output Wholelayer
         nfinedivs_p, xfine_p, wfine_p, radiifine_p, alphafine_p,                           & ! Output partial up
         sinfine_p, cosfine_p, sunpathsfine_p, ntraversefine_p,                             & ! Output partial up
         fail, message, trace_1 )                                                             ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_SSWPGeometry_Master, Solar Sources, Upwelling calculation'
         Master_fail = .true. ; return
       endif

!  Spherical functions call. Updated 3/19/15 for Lattice option. Updated 9/17/16.
!    -- 1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Include option for Doublet-geometry (Add flag and offsets)
!    -- 1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Rearranged the argument list

       Call FO_VectorSS_spherfuncs &
        ( MAXMOMENTS_INPUT, MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, DTR,     & ! Inputs
          NMOMENTS_INPUT, NGEOMS, NSZAS, NVZAS, NAZMS, NSTOKES, VSIGN,    & ! Inputs
          STARTER, DO_OBSGEOM, DO_DOUBLET, DO_SPHERFUNC, DO_SUNLIGHT,     & ! Inputs
          NA_offset, ND_offset, THETA_BOA, ALPHA_BOA, PHI_BOA, COSSCAT, & ! Inputs
          ROTATIONS_UP, GSHELP, GENSPHER_UP )                           ! Outputs

!  RT Call Solar only
!    -- Updated to include surface leaving, 8/2/16. Updated 9/17/16.
!    -- 4/9/19. water-leaving control, Additional cumtrans output for the sleave correction
!    -- 1/31/21. Version 2.8.3. Add LOSTRANS_UP to the output of this routine. 

       Call SSV_Integral_I_UP &   
        ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,           & ! Inputs (dimension)
          do_sunlight, do_deltam, do_fmatrix, do_lambertian, do_surface_leaving, do_water_leaving,& !Inputs (Flags-Gen/Surf)
          do_Partials, do_PlanPar, do_enhanced_ps, flux, fluxvec, do_sources_up, do_sources_up_p, & ! Inputs(Flags/Flux/criticality)
          nstokes, ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, LevelMask_up,       & ! Inputs (control)
          npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,            & ! Inputs (control-partial)
          extinction, deltaus, omega, truncfac, greekmat, fmatrix_up, Reflec, Slterm,             & ! Inputs (Optical/surface)
          Mu0, Mu1, GenSpher_up, Rotations_up, LosW_paths, LosP_paths,                            & ! Inputs (Geometry)
          xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                         & ! Inputs (Geometry)
          xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,             & ! Inputs (Geometry)
          stokes_up, stokes_db, SScumsource_up, Cumtrans, lostrans_up )                             ! Outputs

!  Save results
!mick fix 9/19/2017 - turned off "fo_stokes" (defined later)

       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,upidx) = stokes_up(lev,o1,g)
             fo_stokes_db(lev,g,o1)       = stokes_db(lev,o1,g)
             !fo_stokes(lev,g,o1,upidx)    = fo_stokes_ss(lev,g,o1,upidx) + fo_stokes_db(lev,g,o1)
           enddo
         enddo
       enddo

!  End upwelling

     endif

!  Donwelling

     if ( do_dnwelling ) then
       vsign =  -1.0_ffp

!  Geometry call. Updated 9/17/16.

!    -- 1/31/21. Version 2.8.3. theta_all added to the output, needed for MSST output later on.
!    -- 1/31/21. Version 2.8.3. do_doublet flag added, Offsets added, some inputs rearranged

       call FO_SSWPGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine,          & ! Input dimensions/constants
         do_obsgeom, do_doublet, do_Chapman, do_planpar, do_enhanced_ps, do_Partials,   & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, npartials, partial_layeridx,      & ! Input numbers             
         dtr, Pie, vsign, eradius, nv_offset, na_offset, nd_offset,                     & ! Input Constants/Bookkeeping
         heights, partial_heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,          & ! Input heights/geometry
         doNadir, Raycon, Mu0, Mu1, cosscat,                                                & ! Outputs geometry
         Radii,   LosW_paths, alpha, sina, cosa, sunpaths, ntraverse, chapfacs, theta_all,  & ! Outputs (Layer boundaries)
         Radii_p, LosP_paths, alpha_p, sina_p, cosa_p, sunpaths_p, ntraverse_p, chapfacs_p, & ! Outputs (partial levels)
         nfinedivs,   xfine,   wfine,   radiifine,   alphafine,                             & ! Output Wholelayer
         sinfine,   cosfine,   sunpathsfine,   ntraversefine,                               & ! Output Wholelayer
         nfinedivs_p, xfine_p, wfine_p, radiifine_p, alphafine_p,                           & ! Output partial up
         sinfine_p, cosfine_p, sunpathsfine_p, ntraversefine_p,                             & ! Output partial up
         fail, message, trace_1 )                                                             ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_SSWPGeometry_Master, Solar Sources, Downwelling calculation'
         Master_fail = .true. ; return
       endif

!  Spherical functions call. Updated 3/19/15 for Lattice option. Updated 9/17/16.
!    -- 1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Include option for Doublet-geometry (Add flag and offsets)
!    -- 1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Rearranged the argument list

       Call FO_VectorSS_spherfuncs &
        ( MAXMOMENTS_INPUT, MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, DTR,     & ! Inputs
          NMOMENTS_INPUT, NGEOMS, NSZAS, NVZAS, NAZMS, NSTOKES, VSIGN,    & ! Inputs
          STARTER, DO_OBSGEOM, DO_DOUBLET, DO_SPHERFUNC, DO_SUNLIGHT,     & ! Inputs
          NA_offset, ND_offset, THETA_BOA, ALPHA_BOA, PHI_BOA, COSSCAT, & ! Inputs
          ROTATIONS_DN, GSHELP, GENSPHER_DN )                               ! Outputs

!  RT Call Solar only. Updated 9/17/16.
!    -- 1/31/21. Version 2.8.3. Add lostrans_dn to the output of this routine. 

       Call SSV_Integral_I_DN &
        ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,     & ! Inputs (dimension)
          do_sunlight, do_deltam, do_fmatrix, do_Partials,                                  & ! Inputs (Flags)
          do_PlanPar, do_enhanced_ps, flux, fluxvec, do_sources_dn, do_sources_dn_p,        & ! Inputs (Flags/flux)
          nstokes, ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, LevelMask_dn, & ! Inputs (control)
          npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,      & ! Inputs (control-partial)
          extinction, deltaus, omega, truncfac, greekmat, fmatrix_dn,                       & ! Inputs (Optical)
          Mu1, GenSpher_dn, Rotations_dn, LosW_paths, LosP_paths,                           & ! Inputs (Geometry)
          xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                   & ! Inputs (Geometry)
          xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,       & ! Inputs (Geometry)
          stokes_dn, SScumsource_dn, lostrans_dn )                                            ! Outputs

!  Save results
!mick fix 9/19/2017 - turned off "fo_stokes" (defined later)

       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,dnidx) = stokes_dn(lev,o1,g)
             !fo_stokes(lev,g,o1,dnidx)    = fo_stokes_ss(lev,g,o1,dnidx)
           enddo
         enddo
       enddo

!  End downwelling

     endif

!  End solar only

   endif

!  Thermal sources run
!  -------------------

   if ( do_thermal_emission.and.do_surface_emission ) then

!  Upwelling
!  ---------

     if ( do_upwelling ) then

!  DT Geometry call. Updated 9/17/16. 

!  1/31/21. Version 2.8.3. Following upgrade made for 2.8.1, 5/5/20.
!        ==> Call moved inside the upwelling clause (formerly outside)
!        ==> Add radiifine_LOS + radiifine_p_LOS to the argument list
!        ==> These are vertical distances from layer top, needed for FO Outgoing direct-thermal calculation

       call FO_DTWPGeometry_Master  &
         ( maxvzas, maxlayers, maxpartials, maxfine, dtr, eradius,        & ! Input dimensions/constants
           .true., do_planpar, do_enhanced_ps, do_Partials,               & ! Input flags
           nvzas, nlayers, npartials, nfine, partial_layeridx,            & ! Input control
           heights, alpha_boa, partial_heights,                           & ! Input heights/geometry     
           Mu1_LOS, Radii, LosW_paths_LOS, alpha_LOS, sina_LOS, cosa_LOS, & ! Outputs (Layer boundaries)
           Radii_p, LosP_paths_LOS, alpha_p_LOS, sina_p_LOS, cosa_p_LOS,  & ! Outputs (partial levels)
           nfinedivs_LOS,   xfine_LOS,   wfine_LOS,   radiifine_LOS,      & ! Output Wholelayer
           alphafine_LOS,   sinfine_LOS,   cosfine_LOS,                   & ! Output Wholelayer
           nfinedivs_p_LOS, xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,    & ! Output partial up
           alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS,                 & ! Output partial up
           fail, message, trace_1 )                                         ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_DTWPGeometry_Master, Upwelling'
         Master_fail = .true. ; return
       endif

!  Direct thermal, calculate. Updated 9/17/16
!    -- If Polarized Emissivity flag present, then use optional call. 12/11/17 Rob add.
!    -- Array "Emiss" now has vector dimension.
 
!  1/31/21. Version 2.8.3. Following upgrade made for 2.8.1, 5/5/20.
!        ==> Add radiifine_LOS + radiifine_p_LOS to the argument list, now required inputs
!        ==> Add lostrans_up_LOS + lostrans_up_p_LOS to argument list
!        ==> Use separately defined source-flag arrays (New)

       do_Tsources_up = .true. ; do_Tsources_up_p = .true.

       if ( do_Polarized_Emiss ) then
         call DTE_Integral_I_UP &
         ( maxvzas, maxlayers, maxpartials, maxfine, max_user_levels,                    & ! Inputs (dimensioning)
           Do_Thermset, do_deltam, do_Partials, do_PlanPar,                              & ! Inputs (Flags)
           do_enhanced_ps, do_Tsources_up, do_Tsources_up_p,                             & ! Inputs (Flags)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, LevelMask_up, npartials,        & ! Inputs (control output)
           nfinedivs_p_LOS, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
           bb_input, surfbb, emiss(1,:), extinction, deltaus, omega, truncfac,           & ! Inputs (Optical)
           Mu1_LOS, LosW_paths_LOS, LosP_paths_LOS, xfine_LOS, wfine_LOS, radiifine_LOS, & ! Inputs (Geometry)
           xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,                                    & ! Inputs (Geometry)
           intensity_dta_up_LOS, intensity_dts_LOS, cumsource_up_LOS,                    & ! Main Outputs
           tcom1, lostrans_up_LOS, lostrans_up_p_LOS,                                    & ! Other Outputs
           do_Polarized_Emiss, nstokes, emiss(2:4,:),                                    & ! Optional Input.  12/11/17 Rob Add.
           StokesQUV_dts_LOS )                                                             ! Optional Output. 12/11/17 Rob Add.
       else
         call DTE_Integral_I_UP &
         ( maxvzas, maxlayers, maxpartials, maxfine, max_user_levels,                    & ! Inputs (dimensioning)
           Do_Thermset, do_deltam, do_Partials, do_PlanPar,                              & ! Inputs (Flags)
           do_enhanced_ps, do_Tsources_up, do_Tsources_up_p,                             & ! Inputs (Flags)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, LevelMask_up, npartials,        & ! Inputs (control output)
           nfinedivs_p_LOS, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
           bb_input, surfbb, emiss(1,:), extinction, deltaus, omega, truncfac,           & ! Inputs (Optical)
           Mu1_LOS, LosW_paths_LOS, LosP_paths_LOS, xfine_LOS, wfine_LOS, radiifine_LOS, & ! Inputs (Geometry)
           xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,                                    & ! Inputs (Geometry)
           intensity_dta_up_LOS, intensity_dts_LOS, cumsource_up_LOS,                    & ! Main Outputs
           tcom1, lostrans_up_LOS, lostrans_up_p_LOS )                                     ! Other Outputs
       endif

!  Save results
!mick fix 9/19/2017 - turned off "fo_stokes" (defined later)
!  1/31/21. Version 2.8.3. Do_Doublet upgrade made for 2.8.2, installed here
!  6/30/21.  Ordering wrong with the doublet/lattice geometry settings. SZA then VZA then AZM

       o1=1
       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_stokes_dta(lev,g,o1,upidx) = intensity_dta_up_LOS(lev,g)
             fo_stokes_dts(lev,g,o1)       = intensity_dts_LOS(lev,g)
             !fo_stokes(lev,g,o1,upidx)     = fo_stokes_dta(lev,g,o1,upidx) + fo_stokes_dts(lev,g,o1)
           enddo
         enddo
       else if ( do_Doublet ) then
          do ns = 1, nszas ; do nv = 1, nvzas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,upidx) = intensity_dta_up_LOS(lev,nv)
                fo_stokes_dts(lev,g,o1)       = intensity_dts_LOS(lev,nv)
             enddo
          enddo ; enddo
       else
          do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,upidx) = intensity_dta_up_LOS(lev,nv)
                fo_stokes_dts(lev,g,o1)       = intensity_dts_LOS(lev,nv)
                !fo_stokes(lev,g,o1,upidx)     = fo_stokes_dta(lev,g,o1,upidx) + fo_stokes_dts(lev,g,o1)
             enddo
          enddo ; enddo ; enddo
       endif

!  Save polarized Emissivity results. 12/11/17  Rob add.
!  1/31/21. Version 2.8.3. Do_Doublet upgrade made for 2.8.2, installed here
!  6/30/21.  Ordering wrong with the doublet/lattice geometry settings. SZA then VZA then AZM

       if ( do_Polarized_Emiss.and.nstokes.gt.1 ) then
         if ( do_obsgeom ) then
           do g = 1, nvzas
             do lev=1,n_user_levels
               fo_stokes_dts(lev,g,2:nstokes) = StokesQUV_dts_LOS(lev,1:nstokes-1,g)
             enddo
           enddo
         else if ( do_doublet ) then
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               do lev=1,n_user_levels
                  fo_stokes_dts(lev,g,2:nstokes) = StokesQUV_dts_LOS(lev,1:nstokes-1,nv)
               enddo
            enddo ; enddo
         else
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                  fo_stokes_dts(lev,g,2:nstokes) = StokesQUV_dts_LOS(lev,1:nstokes-1,nv)
               enddo
            enddo ; enddo ; enddo
         endif
       endif

!  End upwelling

     endif

!  Downwelling
!  -----------

     if ( do_dnwelling ) then

!  DT Geometry call. Updated 9/17/16.

!  1/31/21. Version 2.8.3. Following upgrade made for 2.8.1, 5/5/20.
!        ==> Call now introduced inside the downwelling clause (formerly absent)
!        ==> Add radiifine_LOS + radiifine_p_LOS to the argument list
!        ==> These are vertical distances from layer top, needed for FO Outgoing direct-thermal calculation

       call FO_DTWPGeometry_Master  &
         ( maxvzas, maxlayers, maxpartials, maxfine, dtr, eradius,        & ! Input dimensions/constants
           .false., do_planpar, do_enhanced_ps, do_Partials,              & ! Input flags
           nvzas, nlayers, npartials, nfine, partial_layeridx,            & ! Input control
           heights, alpha_boa, partial_heights,                           & ! Input heights/geometry     
           Mu1_LOS, Radii, LosW_paths_LOS, alpha_LOS, sina_LOS, cosa_LOS, & ! Outputs (Layer boundaries)
           Radii_p, LosP_paths_LOS, alpha_p_LOS, sina_p_LOS, cosa_p_LOS,  & ! Outputs (partial levels)
           nfinedivs_LOS,   xfine_LOS,   wfine_LOS,   radiifine_LOS,      & ! Output Wholelayer
           alphafine_LOS,   sinfine_LOS,   cosfine_LOS,                   & ! Output Wholelayer
           nfinedivs_p_LOS, xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,    & ! Output partial up
           alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS,                 & ! Output partial up
           fail, message, trace_1 )                                         ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_DTWPGeometry_Master, Downwelling'
         Master_fail = .true. ; return
       endif

!  Direct thermal, calculate. Updated, 9/17/16.

!  1/31/21. Version 2.8.3. Following upgrade made for 2.8.1, 5/5/20.
!        ==> Add radiifine_LOS + radiifine_p_LOS to the argument list, now required inputs
!        ==> Use separately defined source-flag arrays (New)

       do_Tsources_dn = .true. ; do_Tsources_dn_p = .true.
       
       call DTE_Integral_I_DN &
         ( maxvzas, maxlayers, maxpartials, maxfine, max_user_levels,                    & ! Inputs (dimensioning)
           Do_Thermset, do_deltam, do_Partials, do_PlanPar,                              & ! Inputs (Flags)
           do_enhanced_ps, do_Tsources_dn, do_Tsources_dn_p,                             & ! Inputs (Flags)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, LevelMask_dn, npartials,        & ! Inputs (control output)
           nfinedivs_p_LOS, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
           BB_input, extinction, deltaus, omega, truncfac,                               & ! Inputs (Optical)
           Mu1_LOS, LosW_paths_LOS, LosP_paths_LOS, xfine_LOS, wfine_LOS, radiifine_LOS, & ! Inputs (Geometry)
           xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,                                    & ! Inputs (Geometry)
           intensity_dta_dn_LOS, cumsource_dn_LOS, tcom1 )                                 ! Outputs

!  Save results
!mick fix 9/19/2017 - turned off "fo_stokes" (defined later)
!  1/31/21. Version 2.8.3. Do_Doublet upgrade made for 2.8.2, installed here
!  6/30/21.  Ordering wrong with the doublet/lattice geometry settings. SZA then VZA then AZM

       o1=1
       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,g)
             !fo_stokes(lev,g,o1,dnidx)     = fo_stokes_dta(lev,g,o1,dnidx)
           enddo
         enddo
       else if ( do_doublet ) then
          do ns = 1, nszas ; do nv = 1, nvzas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,nv)
                !fo_stokes(lev,g,o1,dnidx)     = fo_stokes_dta(lev,g,o1,dnidx)
             enddo
          enddo ; enddo
       else
          do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,nv)
                !fo_stokes(lev,g,o1,dnidx)     = fo_stokes_dta(lev,g,o1,dnidx)
             enddo
          enddo ; enddo ; enddo
       endif

!  end downwelling

     endif

!  End thermal run

   endif

!  Final computation of Composites
!  -------------------------------
!mick fix 9/19/2017 - this section overhauled to avoid foreseen bugs & simplify computations

   if ( do_upwelling ) then
     do o1 = 1, nstokes
       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_stokes_atmos(lev,g,o1,upidx) = fo_stokes_ss(lev,g,o1,upidx)    + fo_stokes_dta(lev,g,o1,upidx)
           fo_stokes_surf(lev,g,o1)        = fo_stokes_db(lev,g,o1)          + fo_stokes_dts(lev,g,o1)
           fo_stokes(lev,g,o1,upidx)       = fo_stokes_atmos(lev,g,o1,upidx) + fo_stokes_surf(lev,g,o1)
         enddo
       enddo
     enddo
   endif

   if ( do_dnwelling ) then
     do o1 = 1, nstokes
       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_stokes_atmos(lev,g,o1,dnidx) = fo_stokes_ss(lev,g,o1,dnidx) + fo_stokes_dta(lev,g,o1,dnidx)
           fo_stokes(lev,g,o1,dnidx)       = fo_stokes_atmos(lev,g,o1,dnidx)
         enddo
       enddo
     enddo
   endif

!  Finish

   return
end subroutine VFO_MASTER

end module VFO_Master_m

