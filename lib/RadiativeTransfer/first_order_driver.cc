#include "first_order_driver.h"
#include "wgs84_constant.h"
#include "lidort_driver.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
///
//-----------------------------------------------------------------------

FirstOrderDriver::FirstOrderDriver(int number_layers, int surface_type, int number_moments, bool do_solar, bool do_thermal) 
: num_moments_(number_moments), surface_type_(surface_type), SpurrRtDriver(do_solar, do_thermal)
{
    init_interfaces(number_layers);
}


void FirstOrderDriver::init_interfaces(int nlayers)
{
    // We process 1 geometry at a time for now
    int nszas = 1;
    int nvzas = 1;
    int nazms = 1;
    int ngeoms = nszas * nvzas * nazms;

    // Match what is used by LIDORT driver
    int nfine = 4;

    geometry.reset(new Fo_Ssgeometry_Master(ngeoms, nszas, nvzas, nazms, nlayers, nfine, ngeoms, nszas, nvzas, nazms, nlayers, nfine));

    // Do not compute chapman factors
    geometry->do_chapman(false);

    // Use enhanced pseudo spherical instead of plane parallel
    geometry->do_planpar(false);
    geometry->do_enhanced_ps(true);

    // Use same definitions as Fortran code
    // Define pi and degrees to radians
    geometry->pie(acos(-1.0));
    geometry->dtr(geometry->pie() / 180.0);

    // input angles (Degrees), VSIGN = +1 (Up); -1(Down)
    // scattering angle control, Upwelling only
    geometry->vsign(1);

    // Earth's radius
    geometry->eradius(OldConstant::wgs84_a.convert(units::km).value);

    // Specify angles instead of using observation geometry mode
    geometry->do_obsgeom(false);
   
    // Attenuation control
    // Criticality not set (fast option) valid for lower optical depth atmospheres
    // Only needed when enchanced_ps is enabled for extremely optically thick atmospheres 
    // where light does not reach surface
    geometry->docrit(false);

    // Commented out values not used if docrit == false
    // Critical attenuation
    //geometry->acrit(1.0e-10);
    // Extincion profile
    //geometry->extinct(...)

    // Initialize legendre object

    // Should always be true  before first Spherfuncs call
    bool starter = true;

    legendre.reset(new Fo_Scalarss_Spherfuncs(starter, num_moments_, ngeoms, num_moments_, ngeoms));

    // Initialize radiance object

    // Top of atmosphere only
    int n_user_levels = 1;
    fo_interface.reset(new Fo_Scalarss_Rtcalcs_I(ngeoms, nlayers, nfine, n_user_levels, ngeoms, nlayers, n_user_levels));

    // Same setting as LIDORT would have
    fo_interface->do_deltam_scaling(true);

    // Flags copied from geometry object
    fo_interface->do_planpar(geometry->do_planpar());
    fo_interface->do_regular_ps(!geometry->do_enhanced_ps());
    fo_interface->do_enhanced_ps(geometry->do_enhanced_ps());
    fo_interface->donadir(geometry->donadir());

    // Turn on upwelling calculation or else we get no results
    fo_interface->do_upwelling(true);

    // Set solar flux to 1.0 for solar spectrum case
    fo_interface->flux(1.0);

    // Recommended value by manual of 50 in case we use cox-munk
    int n_brdf_stream = 50;
    boost::shared_ptr<LidortBrdfDriver> l_brdf_driver(new LidortBrdfDriver(n_brdf_stream, num_moments_));
    l_brdf_driver->brdf_interface()->brdf_sup_in().bs_do_directbounce_only(true);
    l_brdf_driver->brdf_interface()->brdf_sup_in().bs_do_brdf_surface(true);
    l_brdf_driver->brdf_interface()->brdf_sup_in().bs_do_user_streams(true);
    l_brdf_driver->brdf_interface()->brdf_sup_in().bs_do_solar_sources(do_solar_sources);
    //l_brdf_driver->brdf_interface()->brdf_sup_in().bs_do_surface_emission(do_thermal_emission);

    brdf_driver_ = l_brdf_driver;
    brdf_driver_->initialize_brdf_inputs(surface_type_);
}

void FirstOrderDriver::setup_height_grid(const blitz::Array<double, 1>& height_grid) const
{
    int nlayers = height_grid.rows() - 1;

    geometry->heights()(Range(0, nlayers)) = height_grid;

    // Used for computing extincts value
    height_diffs.resize(nlayers);
    for (int lev_idx = 1; lev_idx < height_grid.rows(); lev_idx++) {
        height_diffs(lev_idx-1) = height_grid(lev_idx-1) - height_grid(lev_idx);
    }
}

void FirstOrderDriver::setup_geometry(double sza, double azm, double zen) const
{
    Array<double, 1> alpha_boa(geometry->alpha_boa());
    alpha_boa(0) = zen;

    Array<double, 1> theta_boa(geometry->theta_boa());
    theta_boa(0) = sza;

    Array<double, 1> phi_boa(geometry->phi_boa());
    phi_boa(0) = azm;

    // Run geometry computation
    geometry->run();

    // Copy dynamic geometry outputs into FO interface object
    // Consider .reference to avoid copying?
    fo_interface->mu0(geometry->mu0());
    fo_interface->mu1(geometry->mu1());
    fo_interface->ncrit(geometry->ncrit());
    fo_interface->xfine(geometry->xfine());
    fo_interface->wfine(geometry->wfine());
    fo_interface->csqfine(geometry->csqfine());
    fo_interface->cotfine(geometry->cotfine());
    fo_interface->raycon(geometry->raycon());
    fo_interface->cota(geometry->cota());
    fo_interface->sunpaths(geometry->sunpaths());
    fo_interface->ntraverse(geometry->ntraverse());
    fo_interface->sunpathsfine(geometry->sunpathsfine());
    fo_interface->ntraversefine(geometry->ntraversefine());
}

void FirstOrderDriver::setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1> atmosphere_bb) const
{
    // Nothing for now, in future use DT geometry and DT RT
}


void FirstOrderDriver::setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                           const blitz::Array<double, 1>& ssa,
                                           const blitz::Array<double, 2>& pf) const
{

    // Total per layer optical depth
    Array<double, 1> optical_depth(fo_interface->deltaus());
    optical_depth = od;

    // Extinction profile
    Array<double, 1> extinction(fo_interface->extinction());
    extinction = optical_depth / height_diffs;

    // Compute phase function from fourier moments by summing over moments times general spherical function
    Array<double, 2> phase_function(fo_interface->exactscat_up());

    phase_function = 0;
    for(int lay_idx = 0; lay_idx < od.rows(); lay_idx++) {
        for(int mom_idx = 0; mom_idx < num_moments_; mom_idx++) {
            phase_function(lay_idx, 0) = phase_function(lay_idx, 0) + legendre->ss_pleg()(mom_idx, 0) * phase_function(mom_idx, lay_idx);
        }
    }
    
    // Use direct bounce BRDF from LIDORT BRDF supplement for first order reflection
    Array<double, 1> reflectance(fo_interface->reflec());
    boost::shared_ptr<LidortBrdfDriver> l_brdf_driver = boost::dynamic_pointer_cast<LidortBrdfDriver>(brdf_driver());
    reflectance(0) = l_brdf_driver->brdf_interface()->brdf_sup_out().bs_dbounce_brdfunc()(0, 0, 0);
}

void FirstOrderDriver::clear_linear_inputs() const
{
    // Nothing for now
}

void FirstOrderDriver::setup_linear_inputs(const ArrayAd<double, 1>& od, 
                                           const ArrayAd<double, 1>& ssa,
                                           const ArrayAd<double, 2>& pf,
                                           bool do_surface_linearization) const
{
    // Nothing for now
}


void FirstOrderDriver::calculate_rt() const
{
    // Run RT calculation
    fo_interface->ss_integral_i_up();
}

double FirstOrderDriver::get_intensity() const
{
    return fo_interface->intensity_db()(0, 0) + fo_interface->intensity_up()(0, 0);
}

void FirstOrderDriver::copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf) const
{
    // Nothing for now
}
