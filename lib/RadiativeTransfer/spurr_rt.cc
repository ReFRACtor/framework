#include "spurr_rt.h"
#include "ostream_pad.h"
#include <cmath>

#include "spurr_brdf_types.h"
#include "ground_lambertian.h"
#include "ground_lambertian_piecewise.h"
#include "ground_emissivity_polynomial.h"
#include "ground_emissivity_piecewise.h"
#include "ground_coxmunk.h"
#include "ground_coxmunk_plus_lambertian.h"
#include "ground_brdf.h"

using namespace FullPhysics;
using namespace blitz;

// For debugging purposes, it can be useful to dump out the input
// used by this class. We'll leave this code in place, in case we
// need it again, but normally this turned off.
const bool dump_data = false;

//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Atm The atmpsphere to use
/// \param Stokes_coef The Stokes coefficients.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
//-----------------------------------------------------------------------

SpurrRt::SpurrRt(const boost::shared_ptr<RtAtmosphere>& Atm,
                 const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
                 const blitz::Array<double, 1>& Sza,
                 const blitz::Array<double, 1>& Zen,
                 const blitz::Array<double, 1>& Azm,
                 const bool do_solar,
                 const bool do_thermal)
: RadiativeTransferSingleWn(Stokes_coef, Atm),
  do_solar_sources(do_solar), do_thermal_emission(do_thermal),
  sza(Sza.copy()), zen(Zen.copy()), azm(Azm.copy()),
  alt_spec_index_cache(-1), geo_spec_index_cache(-1)
{
  if(sza.rows() != number_spectrometer() ||
     zen.rows() != number_spectrometer() ||
     azm.rows() != number_spectrometer())
    throw Exception("Sza, Zen, and Atm all need to be size number_spectrometer()");
  for(int i = 0; i < number_spectrometer(); ++i) {
    range_check(sza(i), 0.0, 90.0);
    range_check(azm(i), 0.0, 360.0);
    range_check(zen(i), 0.0, 90.0);
  }

  // Watch atmosphere for changes, so we clear cache if needed.
  atm->add_observer(*this);

  // Looks at the type of the Ground class to determine the surface
  // type integer for use in the Spurr RT Fortran code
  // Do this in the constructor since dynamic casting is an expensive operation
  if(dynamic_cast<GroundLambertian*>(atm->ground().get())) {
    surface_type_int = LAMBERTIAN;
  } else if(dynamic_cast<GroundLambertianPiecewise*>(atm->ground().get())) {
    surface_type_int = LAMBERTIAN;
  } else if(dynamic_cast<GroundEmissivityPolynomial*>(atm->ground().get())) {
    surface_type_int = EMISSIVITY;
  } else if(dynamic_cast<GroundEmissivityPiecewise*>(atm->ground().get())) {
    surface_type_int = EMISSIVITY;
   } else if(dynamic_cast<GroundCoxmunk*>(atm->ground().get())) {
    surface_type_int = COXMUNK;
  } else if(dynamic_cast<GroundCoxmunkPlusLambertian*>(atm->ground().get())) {
    surface_type_int = COXMUNK;
  } else if(dynamic_cast<GroundBrdfVeg*>(atm->ground().get())) {
    surface_type_int = BPDFVEGN;
  } else if(dynamic_cast<GroundBrdfSoil*>(atm->ground().get())) {
    surface_type_int = BPDFSOIL;
  } else {
    Exception err_msg;
    err_msg << "Spurr RT can not determine surface type integer from ground class: "
            << atm->ground();
    throw(err_msg);
  }
}

//-----------------------------------------------------------------------
/// Update the altitude information. This can change the number of
/// layers if desired.
//-----------------------------------------------------------------------

void SpurrRt::update_altitude(int spec_index) const
{
  if(spec_index == alt_spec_index_cache)
    return;
  alt_spec_index_cache = spec_index;

  // Set layers into LIDORT inputs
  Array<double, 1> atm_alt(atm->altitude(spec_index).convert(units::km).value.value());

  rt_driver_->setup_height_grid(atm_alt);
  if(dump_data)
    std::cout << "# atm_alt:\n" << atm_alt << "\n";
}

//-----------------------------------------------------------------------
/// Update the geometry if necessary, only needs to change when 
/// spectrometer index changes
//-----------------------------------------------------------------------

void SpurrRt::update_geometry(int spec_index) const
{
  if(spec_index == geo_spec_index_cache)
    return;
  geo_spec_index_cache = spec_index;

  rt_driver_->brdf_driver()->setup_geometry(sza(spec_index), azm(spec_index), zen(spec_index));
  rt_driver_->setup_geometry(sza(spec_index), azm(spec_index), zen(spec_index));
  if(dump_data)
    std::cout << "# Geometry:\n" << sza(spec_index) << "\n"
              << azm(spec_index) << "\n"
              << zen(spec_index) << "\n";
}

// See base class for description of this
Array<double,1> SpurrRt::stokes_single_wn(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  // Obtain Wn and Spec_index dependent inputs
  Range ra(Range::all());

  Array<double, 1> od_in, ssa_in;
  Array<double, 2> pf;
  if(Opt_prop) {
    od_in.reference( Opt_prop->total_optical_depth().value() );
    ssa_in.reference( Opt_prop->total_single_scattering_albedo().value() );
    if (number_moment() > 0)
      pf.reference( Opt_prop->total_phase_function_moments(number_moment(), 1)(ra, ra, 0).value() );
  } else {
    od_in.reference( atm->optical_depth_wrt_rt(Wn, Spec_index).value() );
    ssa_in.reference( atm->single_scattering_albedo_wrt_rt(Wn, Spec_index).value() );
    if (number_moment() > 0)
      pf.reference( atm->phase_function_moments_wrt_rt(Wn, Spec_index, number_moment(), 1)(ra, ra, 0).value() );
  }

  // Update user levels if necessary
  update_altitude(Spec_index);

  // Update geometry if necessary
  update_geometry(Spec_index);

  // Updae thermal emission inputs if enabled
  if(do_thermal_emission) {
      rt_driver_->setup_thermal_inputs(atm->surface_blackbody(Wn, Spec_index).value(), atm->atmosphere_blackbody(Wn, Spec_index).value());
  }

  // Set up BRDF inputs, here we throw away the jacobian
  // value of the surface parameters
  ArrayAd<double, 1> surface_parameters(atm->ground()->surface_parameter(Wn, Spec_index));
  ArrayAd<double, 1> lidort_surface = rt_driver_->brdf_driver()->setup_brdf_inputs(surface_type(), surface_parameters);

  // Set up LIDORT inputs and run
  rt_driver_->setup_optical_inputs(od_in, ssa_in, pf);
  rt_driver_->clear_linear_inputs();
  rt_driver_->calculate_rt();

  // Copy values from LIDORT
  Array<double, 1> stokes(number_stokes());
  stokes = 0;
  stokes(0) = rt_driver_->get_intensity();
  // Check for NaN from lidort
  if(std::isnan(stokes(0)))
    throw Exception("SpurrRt encountered a NaN in the radiance");
  return stokes;
}

// See base class for description of this
ArrayAd<double, 1> SpurrRt::stokes_and_jacobian_single_wn(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{

  // Obtain Wn and Spec_index dependent inputs
  Range ra(Range::all());
  ArrayAd<double, 1> od, ssa;
  ArrayAd<double, 2> pf;
  Array<double, 3> jac_iv(0,0,0);
  if(Opt_prop) {
    od.reference(Opt_prop->total_optical_depth());
    ssa.reference(Opt_prop->total_single_scattering_albedo());
    if (number_moment() > 0)
      pf.reference(Opt_prop->total_phase_function_moments(number_moment(), 1)(ra, ra, 0));
    jac_iv.reference(Opt_prop->intermediate_jacobian());
  } else {
    od.reference(atm->optical_depth_wrt_rt(Wn, Spec_index));
    ssa.reference(atm->single_scattering_albedo_wrt_rt(Wn, Spec_index));
    if (number_moment() > 0)
      pf.reference(atm->phase_function_moments_wrt_rt(Wn, Spec_index, number_moment(), 1)(ra, ra, 0));
    jac_iv.reference(atm->intermediate_jacobian(Wn, Spec_index));
  }

  // Update user levels if necessary
  update_altitude(Spec_index);

  // Update geometry if necessary
  update_geometry(Spec_index);

  // Updae thermal emission inputs if enabled
  AutoDerivative<double> surface_bb;
  ArrayAd<double, 1> atmos_bb;
  if(do_thermal_emission) {
      surface_bb = atm->surface_blackbody(Wn, Spec_index);
      atmos_bb = atm->atmosphere_blackbody(Wn, Spec_index);
      rt_driver_->setup_thermal_inputs(surface_bb.value(), atmos_bb.value());
  }

  // Setup surface
  ArrayAd<double, 1> surface_parameters(atm->ground()->surface_parameter(Wn, Spec_index));
  ArrayAd<double, 1> lidort_surface = rt_driver_->brdf_driver()->setup_brdf_inputs(surface_type(), surface_parameters);

  // Set up LIDORT inputs and run
  rt_driver_->setup_optical_inputs(od.value(), ssa.value(), pf.value());

  bool do_surface_pd = !surface_parameters.is_constant();
  if(dump_data &&
    fabs(Wn - 13050.47) < 0.005)
    std::cout << "# Surface type:\n " << surface_type() << "\n"
              << "# Surface paramters:\n" << surface_parameters << "\n"
              << "# Od:\n" << od << "\n"
              << "# SSA:\n" << ssa << "\n"
              << "# PF:\n" << pf << "\n"
              << "# Do surface pd:\n" << do_surface_pd << "\n";

  rt_driver_->setup_linear_inputs(od, ssa, pf, do_surface_pd);
  rt_driver_->calculate_rt();

  // Copy values from LIDORT
  double rad = rt_driver_->get_intensity();

  Array<double, 2> jac_atm;
  Array<double, 1> jac_surf_param;
  double jac_surf_temp;
  Array<double, 1> jac_atm_temp;
  rt_driver_->copy_jacobians(jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);

  //-----------------------------------------------------------------------
  /// To speed up the calculation, the Atmosphere Jacobian was
  /// calculated relative to the RtAtmosphere "intermediate
  /// variables". The Surface Jacobian was calculated relative to the
  /// surface parameters. For both of these, calculate these relative to
  /// the state vector variables. Then sum the Atmosphere Jacobian over
  /// the layers and add in the surface Jacobian to give us the total
  /// Jacobian to the reflectance with respect to the state vector.
  //-----------------------------------------------------------------------
  Array<double, 1> jac(jac_iv.depth());
  for(int st_idx = 0; st_idx < jac.rows(); ++st_idx) {
    double val = 0;

    // dimensions swapped on jac_iv and jac_atm
    // jac_atm is njac x nlayer 
    // jac_iv  is nlayer x njac
    for(int lay_idx = 0; lay_idx < jac_iv.rows(); ++lay_idx) {
      for(int var_idx = 0; var_idx < jac_iv.cols(); ++var_idx) {
        val += jac_atm(var_idx, lay_idx) * jac_iv(lay_idx, var_idx, st_idx);
      }

      if(do_thermal_emission and !atmos_bb.is_constant()) {
        val += jac_atm_temp(lay_idx) * atmos_bb.jacobian()(lay_idx, st_idx);
      }
    }

    if(do_surface_pd) {
      // The min() here ensures that we only loop over the number of parameters that
      // either the RT or source parameters both have
      // LIDORT will have a jac_surf_param larger insize (hardcoded allocation) than 
      // lidort_surface.jacobian()
      // 2stream has a jac_surf_param smaller than lidort_surface because due to the way
      // it is set up no parameter index is available for shadowing when using
      // coxmunk mode
      for(int m = 0; m < min(jac_surf_param.rows(), lidort_surface.jacobian().rows()); ++m) {
        val += jac_surf_param(m) * lidort_surface.jacobian()(m, st_idx);
      }
    }

    if(do_thermal_emission and !surface_bb.is_constant()) {
      val += surface_bb.gradient()(st_idx) * jac_surf_temp;
    }

    jac(st_idx) = val;
  }

  ArrayAd<double, 1> rad_jac(number_stokes(), jac.rows());
  rad_jac = 0;
  rad_jac(0) = AutoDerivative<double>(rad, jac);
  // Check for NaN from lidort
  if(std::isnan(rad))
    throw Exception("SpurrRt encountered a NaN in the radiance");
  if(any(blitz_isnan(jac)))
    throw Exception("SpurrRt encountered a NaN in the jacobian");
  return rad_jac;
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void SpurrRt::print(std::ostream& Os, bool Short_form) const
{
  Os << "SpurrRt\n";
  OstreamPad opad1(Os, "  ");
  RadiativeTransferSingleWn::print(opad1, Short_form);
  opad1.strict_sync();
  OstreamPad opad(Os, "    ");
  Os << "  Solar zenith angle: \n";
  opad << sza << "\n";
  opad.strict_sync();
  Os << "  Zenith angle: \n";
  opad << zen << "\n";
  opad.strict_sync();
  Os << "  Azimuth angle: \n";
  opad << azm << "\n";
  opad.strict_sync();

  opad << "\n";
  opad.strict_sync();
}
