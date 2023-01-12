#include "spurr_rt_driver.h"
#include "fp_serialize_support.h"
#include "ground.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void SpurrRtDriver::serialize(Archive & ar,
                              const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(SpurrRtDriver);
  ar & FP_NVP(do_solar_sources) & FP_NVP(do_thermal_emission)
    & FP_NVP_(brdf_driver);
}

FP_IMPLEMENT(SpurrRtDriver);

#endif

//-----------------------------------------------------------------------
/// Calculates intensity value with the given inputs
//-----------------------------------------------------------------------

const Array<double,1> SpurrRtDriver::reflectance_calculate(const Array<double, 1>& height_grid,
                                                           double sza, double azm, double zen,
                                                           int surface_type,
                                                           const Array<double, 1>& surface_parameters,
                                                           const Array<double, 1>& od, 
                                                           const Array<double, 1>& ssa,
                                                           const Array<double, 3>& pf,
                                                           double surface_bb,
                                                           const Array<double, 1>& atmosphere_bb)
{
  // Initialize scene 
  setup_height_grid(height_grid);
  brdf_driver_->setup_geometry(sza, azm, zen);
  setup_geometry(sza, azm, zen);

  // Set up BRDF inputs, here we throw away the jacobian
  // value of the surface parameters
  ArrayAd<double, 1> surf_param_ad(surface_parameters.rows(), 0);
  surf_param_ad.value() = surface_parameters;
  ArrayAd<double, 1> lidort_surf = brdf_driver_->setup_brdf_inputs(surface_type, surf_param_ad);

  // Set up LIDORT inputs and run
  setup_optical_inputs(od, ssa, pf);

  if (do_thermal_emission)
      setup_thermal_inputs(surface_bb, atmosphere_bb);

  clear_linear_inputs();
  calculate_rt();

  // Return answer
  return get_intensity();
}

//-----------------------------------------------------------------------
/// Calculates intensity, profile and surface weighting factors (jacobians)
/// with the given inputs
//-----------------------------------------------------------------------

void SpurrRtDriver::reflectance_and_jacobian_calculate(const Array<double, 1>& height_grid,
                                                       double sza, double azm, double zen,
                                                       int surface_type,
                                                       ArrayAd<double, 1>& surface_parameters,
                                                       const ArrayAd<double, 1>& od, 
                                                       const ArrayAd<double, 1>& ssa,
                                                       const ArrayAd<double, 3>& pf,
                                                       Array<double, 1>& reflectance,
                                                       Array<double, 2>& jac_atm, 
                                                       Array<double, 1>& jac_surf_param,
                                                       double& jac_surf_temp,
                                                       Array<double, 1>& jac_atm_temp,
                                                       double surface_bb,
                                                       const Array<double, 1>& atmosphere_bb)

{
  // Initialize scene 
  setup_height_grid(height_grid);
  brdf_driver_->setup_geometry(sza, azm, zen);
  setup_geometry(sza, azm, zen);

  // Set up BRDF inputs and run
  surface_parameters = brdf_driver_->setup_brdf_inputs(surface_type, surface_parameters);

  setup_optical_inputs(od.value(), ssa.value(), pf.value());

  if (do_thermal_emission)
      setup_thermal_inputs(surface_bb, atmosphere_bb);

  bool do_surface_pd = true;
  setup_linear_inputs(od, ssa, pf, do_surface_pd);
  calculate_rt();

  // Copy values from LIDORT
  reflectance.reference(get_intensity());
  copy_jacobians(jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);
}
