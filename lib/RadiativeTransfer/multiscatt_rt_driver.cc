#include "multiscatt_rt_driver.h"

#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "old_constant.h"
#include "wgs84_constant.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void MultiScattRtDriver::serialize(Archive & ar, const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpurrRtDriver)
    & FP_NVP_(surface_type)
    & FP_NVP_(do_multi_scatt_only) 
    & FP_NVP_(pure_nadir);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void MultiScattRtDriver::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void MultiScattRtDriver::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // Nothing more to do
}

FP_IMPLEMENT(MultiScattRtDriver);

#endif

//=======================================================================
// MultiScattRtDriver
//=======================================================================

int MultiScattRtDriver::number_moment() const 
{
  boost::shared_ptr<Spurr_Lps_Masters_Base> interface(rt_interface());

  if (!interface) {
    throw Exception("RT Interface not initialized");
  }
  return interface->modified_inputs_base().modified_control_base().ts_nmoments_input();
}

int MultiScattRtDriver::number_stream() const
{
  boost::shared_ptr<Spurr_Lps_Masters_Base> interface(rt_interface());

  if (!interface) {
    throw Exception("RT Interface not initialized");
  }

  return interface->fixed_inputs_base().fixed_control_base().ts_nstreams();
}

void MultiScattRtDriver::notify_update(const RtAtmosphere& atm)
{
  int stype = atm.ground()->spurr_brdf_type();
  if(stype != surface_type()) {
    initialize_brdf(stype);
  }
}

void MultiScattRtDriver::initialize_brdf(int surface_type)
{
  // Initialize BRDF data structure
  brdf_driver()->initialize_brdf_inputs(surface_type);
  surface_type_ = surface_type;
}

void MultiScattRtDriver::initialize_rt(int nstream, int nmoment, bool do_solar_sources, bool do_thermal_emission, bool do_thermal_scattering) 
{
  // Set up these references for convienence
  Spurr_Fixed_Boolean_Base& fboolean_inputs = rt_interface()->fixed_inputs_base().fixed_boolean_base();
  Spurr_Modified_Control_Base& mcontrol_inputs = rt_interface()->modified_inputs_base().modified_control_base();
  Spurr_Modified_Boolean_Base& mboolean_inputs = rt_interface()->modified_inputs_base().modified_boolean_base();

  Spurr_Fixed_Control_Base& fcontrol_inputs = rt_interface()->fixed_inputs_base().fixed_control_base();
  Spurr_Fixed_Sunrays_Base& fbeam_inputs = rt_interface()->fixed_inputs_base().fixed_sunrays_base();
  Spurr_Modified_Sunrays_Base& mbeam_inputs = rt_interface()->modified_inputs_base().modified_sunrays_base();
  Spurr_Modified_Chapman_Base& mchapman_inputs = rt_interface()->modified_inputs_base().modified_chapman_base();

  Spurr_Modified_Uservalues_Base& muser_inputs = rt_interface()->modified_inputs_base().modified_uservalues_base();
  Spurr_Fixed_Uservalues_Base& fuser_inputs = rt_interface()->fixed_inputs_base().fixed_uservalues_base();

  Spurr_Fixed_Lincontrol_Base& lincontrol = rt_interface()->fixed_lininputs_base().fixed_lincontrol_base();

  // Set values used for all calculations
  // Number of quadtature streams in the cosine half space
  fcontrol_inputs.ts_nstreams(nstream);

  // Number of Legendre expansion coefficients for the phase function
  mcontrol_inputs.ts_nmoments_input(nmoment);

  // Always use BRDF supplement, don't use specialized lambertian_albedo mode
  brdf_interface()->brdf_sup_inputs_base().bs_do_brdf_surface(true);

  // Flags for viewing mode
  fboolean_inputs.ts_do_upwelling(true);
  fboolean_inputs.ts_do_dnwelling(false);

  // Fourier azimuth series is examined twice for convergence
  mboolean_inputs.ts_do_double_convtest(true);

  // Do internal calculation of slanth path optical depths
  mboolean_inputs.ts_do_chapman_function(true);

  // In most instances this flag for delta-m scaling should be set
  mboolean_inputs.ts_do_deltam_scaling(true);

  // There will be output at a number of off-quadrature
  // zenith angles specified by the user, this is the normal case
  mboolean_inputs.ts_do_user_streams(true);
  brdf_interface()->brdf_sup_inputs_base().bs_do_user_streams(true);

  // Accuracy criterion for convergence of Fourier series in
  // relative azimuth. Set to recommended value.
  // Controls lidort_accuracy variable in LIDORT and vlidort_accuracy variable in VLIDORT
  fcontrol_inputs.ts_fourier_accuracy(1e-8);

  // New value in LIDORT 3.8.3 that if left at zero will cause floating point errors in asymtx
  fcontrol_inputs.ts_asymtx_tolerance(1e-20);

  // Beam source flux, same value used for all solar angles
  // Normally set to 1 for "sun-normalized" output
  fbeam_inputs.ts_flux_factor(1.0);

  // Enable solar sources calculations
  if (do_solar_sources) {
      // Needed for atmospheric scattering of sunlight
      mboolean_inputs.ts_do_solar_sources(true);

      // Do internal calculation of slanth path optical depths
      // Only needed for solar sources
      mboolean_inputs.ts_do_chapman_function(true);

      // Enable solar sources in the BRDF driver
      brdf_interface()->brdf_sup_inputs_base().bs_do_solar_sources(true);
  }

  // Enable thermal emission calculation
  if (do_thermal_emission) {
      fboolean_inputs.ts_do_thermal_emission(true);
      fboolean_inputs.ts_do_surface_emission(true);

      // Number of coefficients used in treatment of blackbody emissions in a layer
      // 1 implies constant within a layer
      // 2 implies a lineaer treatment
      fcontrol_inputs.ts_n_thermal_coeffs(2);

      // Enable solar sources in the BRDF driver
      brdf_interface()->brdf_sup_inputs_base().bs_do_surface_emission(true);

      if(!do_thermal_scattering) {
        mboolean_inputs.ts_do_thermal_transonly(true);
      }
  }

  // Number of solar beams
  mbeam_inputs.ts_nbeams(1);

  // Earth's radius in km
  mchapman_inputs.ts_earth_radius( OldConstant::wgs84_a.
                                  convert(units::km).value);

  // Number of user-defined relative azimuth angles
  muser_inputs.ts_n_user_relazms(1);

  // Number of user-defined viewing zenith angles
  muser_inputs.ts_n_user_streams(1);

  // Number of vertical output levels
  fuser_inputs.ts_n_user_levels(1);

  // Flag for calculating profile Jacobians in layer n
  // Calculate jacobians for all layers
  Array<bool, 1> layer_jac_flag(lincontrol.ts_layer_vary_flag());
  layer_jac_flag = true;
  lincontrol.ts_layer_vary_flag(layer_jac_flag);

  // Needs to match the number set in the BRDF structure
  lincontrol.ts_n_surface_wfs(brdf_interface()->brdf_linsup_inputs_base().bs_n_surface_wfs());

} 

/// Set up/reset sphericity mode which may be affected by
/// the current zenith viewing angle
void MultiScattRtDriver::setup_sphericity(blitz::Array<double, 1> UNUSED(zen), bool do_multi_scatt_only, bool pure_nadir)
{
  do_multi_scatt_only_ = do_multi_scatt_only;
  pure_nadir_ = pure_nadir;

  Spurr_Fixed_Boolean_Base& fboolean_inputs = rt_interface()->fixed_inputs_base().fixed_boolean_base();
  Spurr_Modified_Boolean_Base& mboolean_inputs = rt_interface()->modified_inputs_base().modified_boolean_base();
  Spurr_Fixed_Control_Base& fcontrol_inputs = rt_interface()->fixed_inputs_base().fixed_control_base();

  if(do_multi_scatt_only) {
    // Do not do full reflectance calculation if only multiple scattering needed
    // These flags should already be false, but just in case..

    // Only return diffuse multiple scattering component
    fboolean_inputs.ts_do_fullrad_mode(false);

    fboolean_inputs.ts_do_plane_parallel(false);

    // Also SS corr and plane-parallel flags are false
    mboolean_inputs.ts_do_focorr(false);
    mboolean_inputs.ts_do_focorr_nadir(false);
    mboolean_inputs.ts_do_focorr_outgoing(false);

  } else { // if not do multi_scatt_only
    // Do full SS + MS calculation and use LOS correction

    // Do a full reflectance calculation
    fboolean_inputs.ts_do_fullrad_mode(true);

    fboolean_inputs.ts_do_plane_parallel(false);

    // Pseudo-spherical + Line of Sight correction
    mboolean_inputs.ts_do_focorr(true);
    mboolean_inputs.ts_do_focorr_nadir(false);

    if (do_solar_sources) {
        mboolean_inputs.ts_do_focorr_outgoing(true);
    } else {
        // This must be disabled for thermal only mode
        mboolean_inputs.ts_do_focorr_outgoing(false);
    }

    // Number of fine layers subdividing coarse layering
    // Only used during LOS correction
    fcontrol_inputs.ts_nfinelayers(4);
  }

  // Flag for controlling azimuth dependence in the output
  // LIDORT will complain if user zenith is 0 and this is not
  // set, when not using ss correction mode
  if( pure_nadir ) {

    if (mboolean_inputs.ts_do_focorr_outgoing()) {
      // Use Pseudo-spherical correction instead for these small viewing zenith angles
      mboolean_inputs.ts_do_focorr_outgoing(false);
      mboolean_inputs.ts_do_focorr_nadir(true);
    }
  }
}

/// Set plane parallel sphericity
void MultiScattRtDriver::set_plane_parallel()
{
  Spurr_Modified_Boolean_Base& mboolean_inputs = rt_interface()->modified_inputs_base().modified_boolean_base();
  Spurr_Fixed_Boolean_Base& fboolean_inputs = rt_interface()->fixed_inputs_base().fixed_boolean_base();

  mboolean_inputs.ts_do_focorr(false);
  mboolean_inputs.ts_do_focorr_outgoing(false);
  mboolean_inputs.ts_do_focorr_nadir(false);
  fboolean_inputs.ts_do_plane_parallel(true);

  do_multi_scatt_only_ = true;
  pure_nadir_ = false;
}

/// Set pseudo spherical sphericity
void MultiScattRtDriver::set_pseudo_spherical()
{
  // Lidort may cause floating point exceptions when doing setup. This
  // is because it may copy garbage value, which are never used. By
  // chance the garbage values may cause a overflow. We
  // suspend floating point exceptions when doing setup
  FeDisableException disable_fp;

  Spurr_Modified_Boolean_Base& mboolean_inputs = rt_interface()->modified_inputs_base().modified_boolean_base();
  Spurr_Fixed_Boolean_Base& fboolean_inputs = rt_interface()->fixed_inputs_base().fixed_boolean_base();

  mboolean_inputs.ts_do_focorr(true);
  mboolean_inputs.ts_do_focorr_outgoing(false);
  mboolean_inputs.ts_do_focorr_nadir(true);
  fboolean_inputs.ts_do_plane_parallel(false);

  do_multi_scatt_only_ = false;
  pure_nadir_ = false;
}

/// Set plane parallel plus single scattering correction
void MultiScattRtDriver::set_plane_parallel_plus_ss_correction()
{
  Spurr_Modified_Boolean_Base& mboolean_inputs = rt_interface()->modified_inputs_base().modified_boolean_base();
  Spurr_Fixed_Boolean_Base& fboolean_inputs = rt_interface()->fixed_inputs_base().fixed_boolean_base();

  mboolean_inputs.ts_do_focorr(true);
  mboolean_inputs.ts_do_focorr_outgoing(false);
  mboolean_inputs.ts_do_focorr_nadir(true);
  fboolean_inputs.ts_do_plane_parallel(true);

  do_multi_scatt_only_ = false;
  pure_nadir_ = true;
}

/// Set line of sight mode
void MultiScattRtDriver::set_line_of_sight()
{
  Spurr_Modified_Boolean_Base& mboolean_inputs = rt_interface()->modified_inputs_base().modified_boolean_base();
  Spurr_Fixed_Boolean_Base& fboolean_inputs = rt_interface()->fixed_inputs_base().fixed_boolean_base();

  mboolean_inputs.ts_do_focorr(true);
  mboolean_inputs.ts_do_focorr_outgoing(true);
  mboolean_inputs.ts_do_focorr_nadir(false);
  fboolean_inputs.ts_do_plane_parallel(false);

  do_multi_scatt_only_ = false;
  pure_nadir_ = false;
}

void MultiScattRtDriver::setup_height_grid(const blitz::Array<double, 1>& in_height_grid)
{
  Spurr_Fixed_Chapman_Base& fchapman_inputs = rt_interface()->fixed_inputs_base().fixed_chapman_base();
  Spurr_Modified_Uservalues_Base& muser_inputs = rt_interface()->modified_inputs_base().modified_uservalues_base();

  Array<double, 1> rt_height_grid( fchapman_inputs.ts_height_grid() );
  int nlayer = in_height_grid.extent(firstDim) - 1;
  Range lay_range = Range(0,nlayer);
  rt_height_grid(lay_range) = in_height_grid;

  // Set GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
  muser_inputs.ts_geometry_specheight( rt_height_grid(nlayer) );

  // Tell LIDORT number of layers
  Spurr_Fixed_Control_Base& fcontrol_inputs = rt_interface()->fixed_inputs_base().fixed_control_base();
  fcontrol_inputs.ts_nlayers(nlayer);
}

void MultiScattRtDriver::setup_geometry(double sza, double azm, double zen)
{

  Spurr_Modified_Sunrays_Base& mbeam_inputs = rt_interface()->modified_inputs_base().modified_sunrays_base();
  Spurr_Modified_Uservalues_Base& muser_inputs = rt_interface()->modified_inputs_base().modified_uservalues_base();

  // Solar zenith angles (degrees) [0,90]
  Array<double, 1> ld_sza( mbeam_inputs.ts_beam_szas() );
  ld_sza(0) = sza;

  // User-defined relative angles (in degrees) for
  // off-quadrature output.
  Array<double, 1> ld_azm( muser_inputs.ts_user_relazms() );
  ld_azm(0) = azm;

  // User-defined viewing zenith angles (in degrees) for
  // off quadrature output.
  Array<double, 1> ld_zen( muser_inputs.ts_user_angles_input() );
  ld_zen(0) = zen;
}


