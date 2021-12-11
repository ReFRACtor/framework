#include "first_order_driver.h"
#include "fp_serialize_support.h"
#include "wgs84_constant.h"

// Include to use consistent jacobian sizes
#include "lidort_interface_types.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void FirstOrderDriver::serialize(Archive & ar,
				 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpurrRtDriver)
    & FP_NVP_(num_layers)
    & FP_NVP_(num_moments) & FP_NVP_(num_streams)
    & FP_NVP_(surface_type)
    & FP_NVP(height_diffs)
    & FP_NVP(geometry) & FP_NVP(legendre) & FP_NVP_(solar_interface)
    & FP_NVP(l_brdf_driver);
}

FP_IMPLEMENT(FirstOrderDriver);
#endif

//-----------------------------------------------------------------------
/// Construct FirstOrderDriver
//-----------------------------------------------------------------------

FirstOrderDriver::FirstOrderDriver(int number_layers, int surface_type, int number_streams, int number_moments,
                                   bool do_solar, bool do_thermal)
  : SpurrRtDriver(do_solar, do_thermal),
    num_moments_(number_moments), num_streams_(number_streams)
{
    init_interfaces(number_layers, surface_type);
}

void FirstOrderDriver::notify_update(const RtAtmosphere& atm)
{
  int nlayers = atm.number_layer();
  int stype = atm.ground()->spurr_brdf_type();
  if(nlayers != number_layers() || stype != surface_type()) {
    init_interfaces(nlayers, stype);
  }
}

void FirstOrderDriver::init_interfaces(int nlayers, int surface_type)
{

  // Use LIDORT parameters to have consistent sizes for maximum values
  Lidort_Pars lid_pars = Lidort_Pars::instance();
  int max_geoms = lid_pars.max_geometries;
  int max_szas = lid_pars.maxbeams;
  int max_vzas = lid_pars.max_user_streams;
  int max_azms = lid_pars.max_user_relazms;
  int max_layers = lid_pars.maxlayers;
  int max_partials = lid_pars.max_partlayers;
  int max_fine = lid_pars.maxfinelayers;
  int max_moments_input = lid_pars.maxmoments_input;
  int max_user_levels = lid_pars.max_user_levels;
  int max_atmoswfs = lid_pars.max_atmoswfs;
  int max_surfacewfs = lid_pars.max_surfacewfs; 
  int max_sleavewfs = lid_pars.max_sleavewfs;

  // Store the surface type being used, this value is only used for the surface_type() accessor
  surface_type_ = surface_type;

  // Compute actual sizes used for processing < max values
  num_layers_ = nlayers;
  // We process 1 geometry at a time for now
  int nszas = 1;
  int nvzas = 1;
  int nazms = 1;
  int npartials = 0;
  int ngeoms = nszas * nvzas * nazms;
  int n_sleavewfs = 0; // No surface leaving
  int n_user_levels = 1;

  // Match what is used by LIDORT driver
  int nfine = 4;
  
  geometry.reset(new Fo_Sswpgeometry_Master(max_geoms, max_szas, max_vzas, max_azms, max_layers, max_partials, max_fine, ngeoms, nszas, nvzas, nazms, nlayers, nfine, npartials));

  // Use same definitions as Fortran code
  // Define pi and degrees to radians
  geometry->pie(acos(-1.0));
  geometry->dtr(geometry->pie() / 180.0);
  
  // Do not compute chapman factors
  geometry->do_chapman(false);
  
  // input angles (Degrees), VSIGN = +1 (Up); -1(Down)
  // scattering angle control, Upwelling only
  geometry->vsign(1);
  
  // Earth's radius
  geometry->eradius(OldConstant::wgs84_a.convert(units::km).value);
  
  // Specify angles instead of using observation geometry mode
  geometry->do_obsgeom(false);
  
  ///////
  // Initialize legendre interface
  // Used for computing exact scattering from phase moments
  legendre.reset(new Fo_Scalarss_Spherfuncs(max_moments_input, max_geoms, num_moments_, ngeoms));
  
  // Flags to be set for each calculation (safety)
  legendre->starter(true);

  // Need to calculate spherical function 
  legendre->do_spherfunc(true);
  
  ///////
  // Initialize brdf interface

  // Recommended value by manual of 50 in case we use cox-munk
  int n_brdf_stream = 50;
  l_brdf_driver.reset(new LidortBrdfDriver(n_brdf_stream, num_moments_));
  Brdf_Sup_Inputs& brdf_inputs = l_brdf_driver->brdf_interface()->brdf_sup_in();
  
  // Only use 1 beam meaning only one set of sza, azm
  brdf_inputs.bs_nbeams(1);
  brdf_inputs.bs_n_user_streams(1);
  brdf_inputs.bs_n_user_relazms(1);
  
  // This MUST be consistent with streams used for 
  // LIDORT RT calculation
  brdf_inputs.bs_nstreams(num_streams_);
  
  // Number of quadtrature streams for BRDF calculation
  brdf_inputs.bs_nstreams_brdf(n_brdf_stream);
  
  brdf_inputs.bs_do_directbounce_only(true);
  brdf_inputs.bs_do_brdf_surface(true);
  brdf_inputs.bs_do_user_streams(true);
  brdf_inputs.bs_do_solar_sources(do_solar_sources);
  brdf_inputs.bs_do_surface_emission(do_thermal_emission);
  
  // Copy LidortBrdfDriver pointer into SpurrDriver variable
  brdf_driver_ = l_brdf_driver;
  brdf_driver_->initialize_brdf_inputs(surface_type);

  ///////
  // Initialize radiance interface

  // Dynamically determined from BRDF driver setup
  int nsurfacewfs = brdf_driver()->n_surface_wfs();
  
  // Initialize solar mode interface
  solar_interface_.reset(new Fo_Scalarss_Rtcalcs_Ilps(max_geoms, max_layers, max_partials, max_fine, max_moments_input, max_user_levels, max_atmoswfs, ngeoms, nlayers, num_moments_, n_user_levels, npartials, max_surfacewfs, max_sleavewfs, n_sleavewfs, nsurfacewfs));
  
  // Set solar flux to 1.0 for solar spectrum case
  // Adjust flux value to the same meaning as LIDORT's TS_FLUX_FACTOR
  float lidort_flux_factor = 1.0;
  solar_interface_->flux(0.25 * lidort_flux_factor / geometry->pie());

  // Enable sources for all layers
  Array<bool, 2> do_sources_up(solar_interface_->do_sources_up());
  do_sources_up = true;

  // Supply phase moments instead of phase function
  // This is the same way LIDORT would do its internal call to FO
  solar_interface_->do_phasfunc(false);

  // By default enable most accurate mode to match LIDORT defaults
  set_line_of_sight();

  // Enabled by default to match LIDORT behavior and because its generally a good idea
  // to leave this enabled except for testing
  do_deltam_scaling_ = true;
}

/// Set plane parallel sphericity
void FirstOrderDriver::set_plane_parallel()
{
  // LIDORT should have ts_do_sscorr_nadir(false), ts_do_sscorr_outgoing(false)
  geometry->do_planpar(true);
  geometry->do_enhanced_ps(false);
  
  copy_geometry_flags();
}

/// Set pseudo spherical sphericity
void FirstOrderDriver::set_pseudo_spherical()
{
  // Matches LIDORT ts_do_sscorr_nadir(true)
  geometry->do_planpar(false);
  geometry->do_enhanced_ps(false);
  
  copy_geometry_flags();
}

/// Set pseudo spherical sphericity
void FirstOrderDriver::set_line_of_sight()
{
  // Matches LIDORT ts_do_sscorr_outgoing(true)
  geometry->do_planpar(false);
  geometry->do_enhanced_ps(true);
  
  copy_geometry_flags();
}

/// Copy flags for sphericity calculations from geometry object
void FirstOrderDriver::copy_geometry_flags()
{
  // Flags copied from geometry object
  solar_interface_->do_planpar(geometry->do_planpar());
  solar_interface_->do_enhanced_ps(geometry->do_enhanced_ps());
}

void FirstOrderDriver::setup_height_grid(const blitz::Array<double, 1>& height_grid)
{
  int nlayers = height_grid.rows() - 1;
  
  geometry->heights()(Range(0, nlayers)) = height_grid;
  
  // Used for computing extincts value
  height_diffs.resize(nlayers);
  for (int lev_idx = 1; lev_idx < height_grid.rows(); lev_idx++) {
    height_diffs(lev_idx-1) = height_grid(lev_idx-1) - height_grid(lev_idx);
  }
}

void FirstOrderDriver::setup_geometry(double sza, double azm, double zen)
{
  Array<double, 1> alpha_boa(geometry->alpha_boa());
  alpha_boa(0) = zen;
  
  Array<double, 1> theta_boa(geometry->theta_boa());
  theta_boa(0) = sza;
  
  Array<double, 1> phi_boa(geometry->phi_boa());
  phi_boa(0) = azm;
  
  // Run geometry computation
  geometry->run();
  
  if (geometry->fail()) {
    Exception err;
    err << "Error in first order geometry calculation:\n"
        << "Message: " << geometry->message() << "\n"
        << "Trace: " << geometry->trace() << "\n";
    throw err;
  }
    
  // Rerun legendre interface since it dependes on cosscat
  legendre->cosscat(geometry->cosscat());
  legendre->run();

  // Copy dynamic geometry outputs into FO interface object
  solar_interface_->mu0(geometry->mu0());
  solar_interface_->mu1(geometry->mu1());

  solar_interface_->legpoly_up(legendre->ss_pleg());
   
  // Consider .reference to avoid copying?
  Array<int, 2> nfinedivs(solar_interface_->nfinedivs());
  nfinedivs = geometry->nfinedivs();

  solar_interface_->losw_paths(geometry->losw_paths());
  solar_interface_->losp_paths(geometry->losp_paths());
 
  solar_interface_->xfine(geometry->xfine());
  solar_interface_->wfine(geometry->wfine());

  solar_interface_->sunpaths(geometry->sunpaths());
  solar_interface_->ntraverse(geometry->ntraverse());
  solar_interface_->sunpathsfine(geometry->sunpathsfine());
  solar_interface_->ntraversefine(geometry->ntraversefine());

  solar_interface_->xfine_p(geometry->xfine_p());
  solar_interface_->wfine_p(geometry->wfine_p());
  solar_interface_->sunpaths_p(geometry->sunpaths_p());
  solar_interface_->ntraverse_p(geometry->ntraverse_p());
  solar_interface_->sunpathsfine_p(geometry->sunpathsfine_p());
  solar_interface_->ntraversefine_p(geometry->ntraversefine_p());
}

void FirstOrderDriver::setup_thermal_inputs(double UNUSED(surface_bb), const blitz::Array<double, 1>& UNUSED(atmosphere_bb))
{
    // Nothing for now, in future use DT geometry and DT RT
}

/// Compute truncation factor for use in deltam scaling
const blitz::Array<double, 1> FirstOrderDriver::deltam_trunc_factor(const blitz::Array<double, 2>& pf) const
{
    Array<double, 1> truncfac(pf.cols());

    if(pf.rows() < 2*num_streams_+1) {
        // If pf does not contain enough moments then we are likely in a aerosol free case, in that case
        // even if we had a pf matrix that had the number moments available, those entries would be zero
        // since they would not have been able to be filed in with aerosol moments. Therefore its probably
        // a rayleigh only atmosphere.
        truncfac = 0;
    } else {
        double dnm1 = 4 * (num_streams_) + 1;
        truncfac = pf(2*num_streams_, Range::all()) / dnm1;
    }

    return truncfac;
}

void FirstOrderDriver::setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                            const blitz::Array<double, 1>& ssa,
                                            const blitz::Array<double, 2>& pf) 
{
    int nlay = od.rows();
    Range r_all = Range::all();
    Range r_lay = Range(0, nlay-1);
    Range r_mom = Range(0, min(num_moments_, pf.rows()-1)); 

    // Set delta-m scaling flag for solar interface
    solar_interface_->do_deltam_scaling(do_deltam_scaling_);

    // Total per layer optical depth
    Array<double, 1> optical_depth(solar_interface_->deltaus());
    optical_depth(r_lay) = od;

    // Extinction profile
    Array<double, 1> extinction(solar_interface_->extinction());
    extinction(r_lay) = optical_depth(r_lay) / height_diffs;

    // Single scattering albedo
    Array<double, 1> omega(solar_interface_->omega());
    omega(r_lay) = ssa;

    // Compute truncated values when delta_m truncation is turned on or not
    if(do_deltam_scaling_) {
        // Set trunfaction factor for use in FO when deltam scaling is turned on
        Array<double, 1> truncfac(solar_interface_->truncfac());
        truncfac = deltam_trunc_factor(pf);

        // Are thse needed???
        // Doing this to omega would mess up internal deltam scaling code
        optical_depth *= (1 - truncfac * omega);
        //omega *= (1 - truncfac * ssa);
        extinction *= (1 - truncfac * omega);
    }

    // Set phase moments, FO storage is transpose of ReFRACtor ordering
    Array<double, 2> phasmoms(solar_interface_->phasmoms());
    phasmoms(r_lay, r_mom) = pf.transpose(secondDim, firstDim)(r_all, r_mom);

    // Use direct bounce BRDF from LIDORT BRDF supplement for first order reflection
    Array<double, 1> reflectance(solar_interface_->reflec());
    reflectance(0) = l_brdf_driver->brdf_interface()->brdf_sup_out().bs_dbounce_brdfunc()(0, 0, 0);
}

void FirstOrderDriver::clear_linear_inputs() 
{
    solar_interface_->do_profilewfs(false);
    solar_interface_->do_surfacewfs(false); 
    solar_interface_->n_reflecwfs(0);
}

/// Compute truncation factor for use in deltam scaling
const blitz::Array<double, 2> FirstOrderDriver::deltam_linear_trunc_factor(const ArrayAd<double, 2>& pf) const
{
    Array<double, 2> l_truncfac(pf.cols(), pf.number_variable());

    if(pf.rows() < 2*num_streams_+1) {
        // If pf does not contain enough moments then we are likely in a aerosol free case, in that case
        // even if we had a pf matrix that had the number moments available, those entries would be zero
        // since they would not have been able to be filed in with aerosol moments. Therefore its probably
        // a rayleigh only atmosphere.
        l_truncfac = 0;
    } else {
        double dnm1 = 4 * (num_streams_) + 1;
        firstIndex i1; secondIndex i2;
        Range r_all = Range::all();
        l_truncfac = pf.jacobian()(2*num_streams_, r_all, r_all)(i1, i2) / dnm1;
    }

    return l_truncfac;
}


void FirstOrderDriver::setup_linear_inputs
(const ArrayAd<double, 1>& od, 
 const ArrayAd<double, 1>& ssa,
 const ArrayAd<double, 2>& pf,
 bool do_surface_linearization) 
{
    // Set which profile layer jacobians are computed
    int natm_jac = od.number_variable();
    int nlay = od.rows();

    Range r_all = Range::all();
    Range r_jac = Range(0, natm_jac-1);
    Range r_lay = Range(0, nlay-1);

    // Only use the number of common moments available
    Range r_mom = Range(0, min(num_moments_, pf.rows()-1)); 

    // Enable weighting functions
    solar_interface_->do_profilewfs(true);

    Array<bool, 1> layer_jac_flag( solar_interface_->lvaryflags() );
    layer_jac_flag(r_lay) = true;

    Array<int, 1> layer_jac_number( solar_interface_->lvarynums() );
    layer_jac_number(r_lay) = natm_jac;

    // Set up OD related inputs
    Array<double, 2> l_deltau(solar_interface_->l_deltaus());
    Array<double, 2> l_omega(solar_interface_->l_omega());
    Array<double, 2> l_extinction(solar_interface_->l_extinction());

    l_deltau(r_lay, r_jac) = od.jacobian();
    l_omega(r_lay, r_jac) = ssa.jacobian();

    firstIndex i1; secondIndex i2;
    l_extinction(r_lay, r_jac) = od.jacobian()(i1, i2) / height_diffs(i1);

    // Truncate jacobian weights for deltam scaling
    if (do_deltam_scaling_) {
        Array<double, 1> truncfac = deltam_trunc_factor(pf.value());
        Array<double, 2> l_truncfac(solar_interface_->l_truncfac());
        l_truncfac = deltam_linear_trunc_factor(pf);

        Array<double, 1> correction_fac(truncfac.rows());
        correction_fac = 1 - truncfac * ssa.value();

        l_deltau(r_lay, r_jac) = l_deltau(r_lay, r_jac)(i1, i2) * correction_fac(i1) -
            od.value()(i1) * (l_truncfac(i1, i2) * ssa.value()(i1) + truncfac(i1) * ssa.jacobian()(i1, i2));

        l_extinction(r_lay, r_jac) = l_extinction(r_lay, r_jac)(i1, i2) * correction_fac(i1) -
            (od.value()(i1) / height_diffs(i1)) * (l_truncfac(i1, i2) * ssa.value()(i1) + truncfac(i1) * ssa.jacobian()(i1, i2));
    }

    // Set phasmoms jacobian and set correct ordering in copy
    Array<double, 3> l_phasmoms(solar_interface_->l_phasmoms());
    l_phasmoms(r_lay, r_mom, r_jac) = pf.jacobian().transpose(secondDim, firstDim, thirdDim)(r_all, r_mom, r_all);


    // Check for variation of PHASMOMS associated with Jacobian wrt current atmospheric parameter
    // Duplicates behavior in lidort_sfo_lps_interface.f90
    Lidort_Pars lid_pars = Lidort_Pars::instance();

    Array<bool, 2> layer_jac_moms( solar_interface_->lvarymoms() );
    for(int mom_idx = 0; mom_idx < num_moments_; mom_idx++) {
        Array<double, 2> l_mom_vals(l_phasmoms(r_all, mom_idx, r_all));
        layer_jac_moms = where(abs(l_mom_vals) >= 1000.0 * lid_pars.smallnum, true, layer_jac_moms);
    }

    // Set up solar linear inputs
    if (do_surface_linearization) {
        int n_surf_wfs = brdf_driver()->n_surface_wfs();
        Range r_surf_wfs = Range(0, n_surf_wfs-1);

        solar_interface_->do_surfacewfs(true); 
        solar_interface_->n_reflecwfs(n_surf_wfs);
        
        // Set up reflection linearization from BRDF supplement
        Array<double, 1> reflectance(solar_interface_->reflec());
        Array<double, 2> ls_reflec( solar_interface_->ls_reflec() );

        ls_reflec(0, r_surf_wfs) = l_brdf_driver->brdf_interface()->brdf_linsup_out().bs_ls_dbounce_brdfunc()(r_surf_wfs, 0, 0, 0);
    }
 
}


void FirstOrderDriver::calculate_rt() const
{
    // Run RT calculation
    solar_interface_->ss_integral_ilps_up();
}

double FirstOrderDriver::get_intensity() const
{
    return solar_interface_->intensity_db()(0, 0) + solar_interface_->intensity_up()(0, 0);
}

void FirstOrderDriver::copy_jacobians
(blitz::Array<double, 2>& jac_atm,
 blitz::Array<double, 1>& jac_surf_param,
 double& UNUSED(jac_surf_temp),
 blitz::Array<double, 1>& UNUSED(jac_atm_temp)) const
{
    Range ra(Range::all());

    // Output array size is: max_user_levels, maxgeoms, maxlayers, max_atmoswf
    // We only have 1 user_level and 1 geom
    Array<double, 2> jac_total(solar_interface_->lp_jacobians_up()(0, 0, ra, ra) + solar_interface_->lp_jacobians_db()(0, 0, ra, ra));
    // Need to transpose output to be in the expected order of njac, nlay
    jac_total.transposeSelf(secondDim, firstDim);
    jac_atm.reference(jac_total);

    jac_surf_param.resize(solar_interface_->max_surfacewfs());

    // Output array size is: max_user_levels, maxgeoms, max_surfacewfs
    jac_surf_param = solar_interface_->ls_jacobians_db()(0, 0, ra);
}
