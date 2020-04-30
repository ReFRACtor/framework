#include "first_order_driver.h"
#include "wgs84_constant.h"

// Include to use consistent jacobian sizes
#include "lidort_interface_types.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Construct FirstOrderDriver
//-----------------------------------------------------------------------

FirstOrderDriver::FirstOrderDriver(int number_layers, int surface_type, int number_streams, int number_moments,
                                   bool do_solar, bool do_thermal)
  : SpurrRtDriver(do_solar, do_thermal),
    num_moments_(number_moments), num_streams_(number_streams)
{
    init_interfaces(number_layers, surface_type);

    // By default enable most accurate mode to match LIDORT defaults
    set_line_of_sight();

    // Enabled by default to match LIDORT behavior and because its generally a good idea
    // to leave this enabled except for testing
    do_deltam_scaling_ = true;
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

    // Use LIDORT parameters to have consistent sizes for maximum weighting functions
    Lidort_Pars lid_pars = Lidort_Pars::instance();
    int natmoswfs = lid_pars.max_atmoswfs;
    int nsurfacewfs = lid_pars.max_surfacewfs; 

    // Initialize solar mode interface
    solar_interface_.reset(new Fo_Scalarss_Rtcalcs_Ilps(ngeoms, nlayers, nfine, natmoswfs, nsurfacewfs, ngeoms, nlayers));
    
    // Set solar flux to 1.0 for solar spectrum case
    // Adjust flux value to the same meaning as LIDORT's TS_FLUX_FACTOR
    float lidort_flux_factor = 1.0;
    solar_interface_->flux(0.25 * lidort_flux_factor / geometry->pie());

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
}

/// Set plane parallel sphericity
void FirstOrderDriver::set_plane_parallel() const
{
    // LIDORT should have ts_do_sscorr_nadir(false), ts_do_sscorr_outgoing(false)
    geometry->do_planpar(true);
    geometry->do_enhanced_ps(false);

    copy_geometry_flags();
}

/// Set pseudo spherical sphericity
void FirstOrderDriver::set_pseudo_spherical() const
{
    // Matches LIDORT ts_do_sscorr_nadir(true)
    geometry->do_planpar(false);
    geometry->do_enhanced_ps(false);

    copy_geometry_flags();
}

/// Set pseudo spherical sphericity
void FirstOrderDriver::set_line_of_sight() const
{
    // Matches LIDORT ts_do_sscorr_outgoing(true)
    geometry->do_planpar(false);
    geometry->do_enhanced_ps(true);

    copy_geometry_flags();
}

/// Copy flags for sphericity calculations from geometry object
void FirstOrderDriver::copy_geometry_flags() const
{
    // Flags copied from geometry object
    solar_interface_->do_planpar(geometry->do_planpar());
    solar_interface_->do_regular_ps(!geometry->do_planpar() && !geometry->do_enhanced_ps());
    solar_interface_->do_enhanced_ps(geometry->do_enhanced_ps());
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
    Array<int, 2> nfinedivs(solar_interface_->nfinedivs());
    nfinedivs = geometry->nfinedivs();

    solar_interface_->donadir(geometry->donadir());
    solar_interface_->ncrit(geometry->ncrit());
    solar_interface_->xfine(geometry->xfine());
    solar_interface_->wfine(geometry->wfine());
    solar_interface_->csqfine(geometry->csqfine());
    solar_interface_->cotfine(geometry->cotfine());
    solar_interface_->raycon(geometry->raycon());
    solar_interface_->cota(geometry->cota());
    solar_interface_->sunpaths(geometry->sunpaths());
    solar_interface_->ntraverse(geometry->ntraverse());
    solar_interface_->sunpaths_fine(geometry->sunpathsfine());
    solar_interface_->ntraverse_fine(geometry->ntraversefine());
    
    // Compute these values for plane parallel or regular pseudo spherical modes
    if (geometry->do_planpar() || (!geometry->do_planpar() && !geometry->do_enhanced_ps())) {
        Array<double, 1> mu0(solar_interface_->mu0());
        Array<double, 1> mu1(solar_interface_->mu1());

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
        solar_interface_->mu0(geometry->mu0());
        solar_interface_->mu1(geometry->mu1());
    }
}

void FirstOrderDriver::setup_thermal_inputs(double UNUSED(surface_bb), const blitz::Array<double, 1>& UNUSED(atmosphere_bb)) const
{
    // Nothing for now, in future use DT geometry and DT RT
}

/// Compute truncation factor for use in deltam scaling
const blitz::Array<double, 1> FirstOrderDriver::deltam_trunc_factor(const blitz::Array<double, 2>& pf) const
{
    double dnm1 = 4 * (num_streams_) + 1;
    Array<double, 1> truncfac(pf.cols());
    truncfac = pf(2*num_streams_, Range::all()) / dnm1;
    return truncfac;
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
    Array<double, 1> optical_depth(solar_interface_->deltaus());
    optical_depth = od;

    // Extinction profile
    Array<double, 1> extinction(solar_interface_->extinction());
    extinction = optical_depth / height_diffs;

    // Compute truncated values when delta_m truncation is turned on or not
    Array<double, 1> tms;
    if(do_deltam_scaling_) {
        Array<double, 1> truncfac = deltam_trunc_factor(pf);
        optical_depth *= (1 - truncfac * ssa);
        extinction *= (1 - truncfac * ssa);

        tms.resize(ssa.rows());
        tms = ssa / (1 - truncfac * ssa);
    } else {
        tms.reference(ssa);
    }

    // Compute phase function from fourier moments by summing over moments times general spherical function
    // Sum over moment index
    firstIndex lay_idx; secondIndex geom_idx; thirdIndex mom_idx;
    Array<double, 2> exactscat(solar_interface_->exactscat_up());
    exactscat = sum(legendre->ss_pleg()(mom_idx, geom_idx) * pf(mom_idx, lay_idx), mom_idx) * tms(lay_idx);
    
    // Use direct bounce BRDF from LIDORT BRDF supplement for first order reflection
    Array<double, 1> reflectance(solar_interface_->reflec());
    reflectance(0) = l_brdf_driver->brdf_interface()->brdf_sup_out().bs_dbounce_brdfunc()(0, 0, 0);
}

void FirstOrderDriver::clear_linear_inputs() const
{
    solar_interface_->do_profilewfs(false);
    solar_interface_->do_reflecwfs(false); 
    solar_interface_->n_reflecwfs(0);
}

void FirstOrderDriver::setup_linear_inputs
(const ArrayAd<double, 1>& od, 
 const ArrayAd<double, 1>& ssa,
 const ArrayAd<double, 2>& pf,
 bool do_surface_linearization) const
{
    // Set which profile layer jacobians are computed
    int natm_jac = od.number_variable();
    int nlay = od.rows();
    int ngeom = geometry->ngeoms();

    Range r_all = Range::all();
    Range r_jac = Range(0, natm_jac-1);
    Range r_lay = Range(0, nlay-1);

    // Enable weighting functions
    solar_interface_->do_profilewfs(true);

    Array<bool, 1> layer_jac_flag( solar_interface_->lvaryflags() );
    layer_jac_flag = true;

    Array<int, 1> layer_jac_number( solar_interface_->lvarynums() );
    layer_jac_number = natm_jac;

    // Set up OD related inputs
    Array<double, 2> l_deltau(solar_interface_->l_deltaus());
    Array<double, 2> l_extinction(solar_interface_->l_extinction());

    l_deltau(r_lay, r_jac) = od.jacobian();

    firstIndex i1; secondIndex i2;
    l_extinction(r_lay, r_jac) = od.jacobian()(i1, i2) / height_diffs(i1);

    // Truncate jacobian weights for deltam scaling
    Array<double, 2> l_tms;
    if (do_deltam_scaling_) {
        Array<double, 1> truncfac = deltam_trunc_factor(pf.value());

        Array<double, 1> od_correction(od.value().rows());
        for (int par_idx = 0; par_idx < natm_jac; par_idx++) {
            od_correction = (ssa.value() * truncfac * ssa.jacobian()(r_all, par_idx)) /
                            (1 - ssa.value() * truncfac);
            //l_deltau(r_all, par_idx) -= od_correction;
           // l_extinction(r_all, par_idx) -= od_correction;
        }
      
        l_tms.resize(ssa.jacobian().shape());
        l_tms = ssa.jacobian()(i1, i2) / (1 - truncfac(i1) * ssa.value()(i1));
    } else {
        l_tms.reference(ssa.jacobian());
    }

    // l_exactscat takes into account ssa jacobian contributions
    // Compute legendre * pf function first seperately so that separate placeholder
    // variables can be used when computing l_exactscat to ensure correct ordering of values
    Array<double, 3> l_exactscat(solar_interface_->l_exactscat_up());

    Array<double, 2> legpf(nlay, ngeom);
    firstIndex j1; secondIndex j2; thirdIndex j3;
    legpf = sum(legendre->ss_pleg()(j3, j2) * pf.value()(j3, j1), j3);

    firstIndex k1; secondIndex k2; thirdIndex k3;
    l_exactscat(r_lay, r_all, r_jac) = legpf(k1, k2) * l_tms(k1, k3);

    // Set up solar linear inputs
    if (do_surface_linearization) {
        int n_surf_wfs = brdf_driver()->n_surface_wfs();
        Range r_surf_wfs = Range(0, n_surf_wfs-1);

        solar_interface_->do_reflecwfs(true); 
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
    return solar_interface_->intensity_db()(0) + solar_interface_->intensity_up()(0);
}

void FirstOrderDriver::copy_jacobians
(blitz::Array<double, 2>& jac_atm,
 blitz::Array<double, 1>& jac_surf_param,
 double& UNUSED(jac_surf_temp),
 blitz::Array<double, 1>& UNUSED(jac_atm_temp)) const
{
    Range ra(Range::all());

    // Need to transpose output to be in the expected order of njac, nlay
    Array<double, 2> jac_total(solar_interface_->lp_jacobians_up()(0, ra, ra) + solar_interface_->lp_jacobians_db()(0, ra, ra));
    jac_total.transposeSelf(secondDim, firstDim);
    jac_atm.reference(jac_total);

    jac_surf_param.resize(solar_interface_->max_surfacewfs());
    jac_surf_param = solar_interface_->ls_jacobians_db()(0, ra);
}
