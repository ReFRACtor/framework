// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "spurr_driver.h"
%}

%import "array_ad.i"
%fp_shared_ptr(FullPhysics::SpurrBrdfDriver);
%fp_shared_ptr(FullPhysics::SpurrRtDriver);

namespace FullPhysics {
class SpurrBrdfDriver {
public:
  virtual void initialize_brdf_inputs(int surface_type);
  virtual void setup_geometry(double sza, double azm, double zen) const = 0;
  virtual ArrayAd<double, 1> setup_brdf_inputs(int surface_type, const ArrayAd<double, 1>& surface_parameters) const;

  %python_attribute_abstract(n_brdf_kernels, int)
  %python_attribute_abstract(n_kernel_factor_wfs, int)
  %python_attribute_abstract(n_kernel_params_wfs, int)
  %python_attribute_abstract(n_surface_wfs, int)
  %python_attribute_abstract(do_shadow_effect, bool)
};

class SpurrRtDriver {
public:

  SpurrRtDriver(bool do_solar = true, bool do_thermal = false);

  virtual double reflectance_calculate(const blitz::Array<double, 1>& height_grid,
                                       double sza, double azm, double zen,
                                       int surface_type,
                                       const blitz::Array<double, 1>& surface_parameters,
                                       const blitz::Array<double, 1>& od, 
                                       const blitz::Array<double, 1>& ssa,
                                       const blitz::Array<double, 2>& pf,
                                       double surface_bb = 0,
                                       const blitz::Array<double, 1>& atmosphere_bb = blitz::Array<double,1>());

  virtual void reflectance_and_jacobian_calculate(const blitz::Array<double, 1>& height_grid,
                                                  double sza, double azm, double zen,
                                                  int surface_type,
                                                  ArrayAd<double, 1>& surface_parameters,
                                                  const ArrayAd<double, 1>& od, 
                                                  const ArrayAd<double, 1>& ssa,
                                                  const ArrayAd<double, 2>& pf,
                                                  double& reflectance,
                                                  blitz::Array<double, 2>& jac_atm, 
                                                  blitz::Array<double, 1>& jac_surf_param,
                                                  double &jac_surf_temp,
                                                  double surface_bb = 0,
                                                  const blitz::Array<double, 1>& atmosphere_bb = blitz::Array<double,1>());

  %python_attribute(brdf_driver, boost::shared_ptr<SpurrBrdfDriver>)
  virtual void setup_height_grid(const blitz::Array<double, 1>& height_grid) const = 0;
  virtual void setup_geometry(double sza, double azm, double zen) const = 0;
  virtual void setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1> atmosphere_bb) const = 0;
  virtual void setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                    const blitz::Array<double, 1>& ssa,
                                    const blitz::Array<double, 2>& pf) const = 0;
  virtual void clear_linear_inputs() const =  0;
  virtual void setup_linear_inputs(const ArrayAd<double, 1>& od,
                                   const ArrayAd<double, 1>& ssa,
                                   const ArrayAd<double, 2>& pf,
                                   bool do_surface_linearization) const = 0;
  virtual void calculate_rt() const = 0;
  virtual double get_intensity() const = 0;
  virtual void copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf_params, double& jac_surf_temp) const = 0;
};
}

