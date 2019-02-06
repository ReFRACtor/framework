#include "first_order_driver.h"
#include "wgs84_constant.h"
#include "lidort_driver.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
///
//-----------------------------------------------------------------------

FirstOrderDriver::FirstOrderDriver(int number_layers, int surface_type, int number_streams, int number_moments, bool do_solar, bool do_thermal) 
: num_streams_(number_streams), num_moments_(number_moments), SpurrRtDriver(do_solar, do_thermal)
{
    init_interfaces(number_layers, surface_type);
    // Use pseudo spherical correction by default
    set_pseudo_spherical();
}


void FirstOrderDriver::init_interfaces(int nlayers, int surface_type)
{
    // We process 1 geometry at a time for now
    int nszas = 1;
    int nvzas = 1;
    int nazms = 1;
    int ngeoms = nszas * nvzas * nazms;

    // Match what is used by LIDORT driver
    int nfine = 4;

    geometry.reset(new Fo_Ssgeometry_Master(ngeoms, nszas, nvzas, nazms, nlayers, nfine, ngeoms, nszas, nvzas, nazms, nlayers, nfine));

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

    ///////
    // Initialize legendre interface
    // Used for computing exact scattering from phase moments

    // Should always be true before first Spherfuncs call
    bool starter = true;
    legendre.reset(new Fo_Scalarss_Spherfuncs(starter, num_moments_, ngeoms, num_moments_, ngeoms));

    ///////
    // Initialize radiance interface

    // Top of atmosphere only
    int n_user_levels = 1;
    fo_interface.reset(new Fo_Scalarss_Rtcalcs_I(ngeoms, nlayers, nfine, n_user_levels, ngeoms, nlayers, n_user_levels));
    
    // Same setting as LIDORT would have
    fo_interface->do_deltam_scaling(true);

    // Set solar flux to 1.0 for solar spectrum case
    // Adjust flux value to the same meaning as LIDORT's TS_FLUX_FACTOR
    float lidort_flux_factor = 1.0;
    fo_interface->flux(0.25 * lidort_flux_factor / geometry->pie());

    // Recommended value by manual of 50 in case we use cox-munk
    int n_brdf_stream = 50;
    boost::shared_ptr<LidortBrdfDriver> l_brdf_driver(new LidortBrdfDriver(n_brdf_stream, num_moments_));
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
    //brdf_inputs.bs_do_surface_emission(do_thermal_emission);

    brdf_driver_ = l_brdf_driver;
    brdf_driver_->initialize_brdf_inputs(surface_type);
}

/// Set plane parallel sphericity
void FirstOrderDriver::set_plane_parallel() const
{
    geometry->do_planpar(true);
    geometry->do_enhanced_ps(false);

    copy_geometry_flags();
}

/// Set pseudo spherical sphericity
void FirstOrderDriver::set_pseudo_spherical() const
{
    geometry->do_planpar(false);
    geometry->do_enhanced_ps(true);

    copy_geometry_flags();
}

/// Copy flags for sphericity calculations from geometry object
void FirstOrderDriver::copy_geometry_flags() const
{
    // Flags copied from geometry object
    fo_interface->do_planpar(geometry->do_planpar());
    fo_interface->do_regular_ps(!geometry->do_planpar() && !geometry->do_enhanced_ps());
    fo_interface->do_enhanced_ps(geometry->do_enhanced_ps());
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
    // Consider .reference to avoid copying?
    Array<int, 2> nfinedivs(fo_interface->nfinedivs());
    nfinedivs = geometry->nfinedivs();

    fo_interface->donadir(geometry->donadir());
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
    
    if (geometry->do_planpar()) {
        // Account for a bug in the plane parallel version of the geometry routine
        // where these values are not computed 
        Array<double, 1> mu0(fo_interface->mu0());
        Array<double, 1> mu1(fo_interface->mu1());

        for (int ns = 0; ns < geometry->nszas(); ns++) {
            int nv_offset = geometry->nvzas() * geometry->nazms() * ns;
            for(int nv = 0; nv < geometry->nvzas(); nv++) {
                int na_offset = nv_offset + geometry->nazms() * nv;
                for(int na = 0; na < geometry->nazms(); na++) {
                    // Geometry index for lattice
                    int g =  na_offset + na;
                    mu0(g) = cos(geometry->theta_boa()(ns) * geometry->dtr());
                    mu1(g) = cos(geometry->alpha_boa()(nv) * geometry->dtr());
                }
            }
        }
    } else {
        fo_interface->mu0(geometry->mu0());
        fo_interface->mu1(geometry->mu1());
    }
}

void FirstOrderDriver::setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1> atmosphere_bb) const
{
    // Nothing for now, in future use DT geometry and DT RT
}


void FirstOrderDriver::setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                            const blitz::Array<double, 1>& ssa,
                                            const blitz::Array<double, 2>& pf) const
{
    if (pf.rows() != (num_moments_+1)) {
        // By definition num_moments should be one less than the size of pf array
        // because in the Fortran code phase function arrays are index from 0:num_moments
        Exception err;
        err << "Number of moments specified to First Order: " << num_moments_ 
            << " need to be one less than the size of the pf array: " << pf.rows();
        throw err;
    }

    // Total per layer optical depth
    Array<double, 1> optical_depth(fo_interface->deltaus());
    optical_depth = od;

    // Extinction profile
    Array<double, 1> extinction(fo_interface->extinction());
    extinction = optical_depth / height_diffs;

    // Compute phase function from fourier moments by summing over moments times general spherical function
    Array<double, 2> exactscat(fo_interface->exactscat_up());

    exactscat = 0;

    Array<double, 2> moment_sum(pf.cols(), geometry->ngeoms());
    moment_sum = 0;
    for (int geom_idx = 0; geom_idx < geometry->ngeoms(); geom_idx++) {
        for (int lay_idx = 0; lay_idx < pf.cols(); lay_idx++) {
            for(int mom_idx = 0; mom_idx < pf.rows(); mom_idx++) {
                moment_sum(lay_idx, geom_idx) += legendre->ss_pleg()(mom_idx, geom_idx) * pf(mom_idx, lay_idx);
            }
        }
    }

    // For TMS truncation correction
    double dnm1 = 4 * num_streams_ + 1;
    for (int geom_idx = 0; geom_idx < geometry->ngeoms(); geom_idx++) {
        for (int lay_idx = 0; lay_idx < geometry->nlayers(); lay_idx++) {
            exactscat(lay_idx, geom_idx) = moment_sum(lay_idx, geom_idx);
            double omw = ssa(lay_idx);
            double tms;
            if (fo_interface->do_deltam_scaling()) {
                double truncfac =  pf(2 * num_streams_ -1, lay_idx) / dnm1;
                tms = omw / (1.0 - truncfac * omw);
            } else {
                tms = omw;
            }
            exactscat(lay_idx, geom_idx) *= tms;
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

void FirstOrderDriver::copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf_param, double& jac_surf_temp) const
{
    // Nothing for now
}
