#ifndef SPURR_RT_DRIVER_H
#define SPURR_RT_DRIVER_H

#include <boost/shared_ptr.hpp>
#include <blitz/array.h>
#include "array_ad.h"
#include "generic_object.h"

/****************************************************************//**
  Contains classes to abstract away details in various Spurr
  Radiative Transfer software.
*******************************************************************/

namespace FullPhysics {

/****************************************************************//**
  Abstracts away set up of BRDF kernel interfaces.

  The arrays needing setting for the initialization steps are done
  through methods because it is often the case that boolean
  arrays passed to Fortran will need to be copied and can not be
  set directly into the memory read by Fortran. Since this step
  is only done once this is probably acceptable.

  For the setup of the values for the surface to be calculated,
  three double arrays are exposed so that this interface might
  copy values into the memory used by Fortran once instead of 
  twice if method calls were used to do the setting.

  The public inteface only exposes the two main initialization 
  and setup routines. However, most of the protected interface
  needs to be implemented by implementing classes. These are 
  mostly methods to link a parameter to the location it is stored
  by the specific incarnation of Spurr BRDF code.
*******************************************************************/

class SpurrBrdfDriver : public virtual GenericObject {

public:
  virtual void initialize_brdf_inputs(int surface_type);
  virtual void setup_geometry(double sza, double azm, double zen) const = 0;
  virtual ArrayAd<double, 1> setup_brdf_inputs(int surface_type, const ArrayAd<double, 1>& surface_parameters) const;

  virtual int n_brdf_kernels() const = 0;

  virtual int n_kernel_factor_wfs() const = 0;
  virtual int n_kernel_params_wfs() const = 0;
  virtual int n_surface_wfs() const = 0;
  virtual bool do_kparams_derivs(const int kernel_index) const = 0;
  virtual bool do_shadow_effect() const = 0;

protected:
  virtual void initialize_brdf_kernel(int kernel_type);

  virtual void setup_lambertian_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const blitz::Array<int, 1>& parameter_indexes) const;
  virtual void setup_coxmunk_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const blitz::Array<int, 1>& parameter_indexes) const;
  virtual void setup_rahman_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const blitz::Array<int, 1>& parameter_indexes) const;
  virtual void setup_breon_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const blitz::Array<int, 1>& parameter_indexes) const;

  virtual void calculate_brdf() const = 0;

  virtual void n_brdf_kernels(const int n_kernels) = 0;
  virtual void n_kernel_factor_wfs(const int n_factors) = 0;
  virtual void n_kernel_params_wfs(const int n_params) = 0;
  virtual void n_surface_wfs(const int n_wfs) = 0;
  virtual void do_kparams_derivs(const int kernel_index, const bool do_kparams) = 0;
  virtual void do_shadow_effect(const bool do_shadow) const = 0;

  virtual void initialize_kernel_parameters(const int kernel_index,
                                            const int which_brdf,
                                            const bool lambertian_flag,
                                            const int n_brdf_parameters,
                                            const bool do_factor_wfs,
                                            const blitz::Array<bool, 1>& do_params_wfs) = 0;

  // To speed up brdf setup, ie reduce the number of function calls
  // These are set through attributes linked to a valid array by the implementing class. 
  mutable blitz::Array<double, 1> brdf_factors;
  mutable blitz::Array<double, 2> brdf_params;
};

/****************************************************************//**
  Abstracts away set up of Radiative Transfer software from Rob
  Spurr into a simpler common inteface used by the L2 software

  This interface should be independent of the L2 Atmosphere
  class to make unit testing easier.
*******************************************************************/

class SpurrRtDriver : public virtual GenericObject {

public:

  SpurrRtDriver(bool do_solar = true, bool do_thermal = false) 
      : do_solar_sources(do_solar), do_thermal_emission(do_thermal) {}

  /// Computes reflectance without jacobians
  virtual double reflectance_calculate(const blitz::Array<double, 1>& height_grid,
                                       double sza, double azm, double zen,
                                       int surface_type,
                                       const blitz::Array<double, 1>& surface_parameters,
                                       const blitz::Array<double, 1>& od, 
                                       const blitz::Array<double, 1>& ssa,
                                       const blitz::Array<double, 2>& pf,
                                       double surface_bb = 0,
                                       const blitz::Array<double, 1>& atmosphere_bb = blitz::Array<double,1>());
  
  // Computes reflectance and jacobians for profiles as well as surface
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
                                                  blitz::Array<double, 1>& jac_atm_temp,
                                                  double surface_bb = 0,
                                                  const blitz::Array<double, 1>& atmosphere_bb = blitz::Array<double,1>());

  /// Access to BRDF driver
  const boost::shared_ptr<SpurrBrdfDriver> brdf_driver() const { return brdf_driver_; }

  /// Setup height grid, should only be called once per instance or if
  /// the height grid changes
  virtual void setup_height_grid(const blitz::Array<double, 1>& height_grid) const = 0;

  /// Setup viewing geometry, should only be called once per instance or if
  /// the viewing geometry changes
  virtual void setup_geometry(double sza, double azm, double zen) const = 0;

  /// Set up thermal emission inputs
  virtual void setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1>& atmosphere_bb) const = 0;

  /// Set up optical depth, single scattering albedo and phase function
  /// Should be called per spectral point
  virtual void setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                    const blitz::Array<double, 1>& ssa,
                                    const blitz::Array<double, 2>& pf) const = 0;
  
  /// Mark that we are not retrieving weighting functions
  virtual void clear_linear_inputs() const =  0;

  /// Set up linearization, weighting functions
  virtual void setup_linear_inputs(const ArrayAd<double, 1>& od,
                                   const ArrayAd<double, 1>& ssa,
                                   const ArrayAd<double, 2>& pf,
                                   bool do_surface_linearization) const = 0;

  /// Perform radiative transfer calculation with the values
  /// setup by setup_optical_inputs and setup_linear_inputs
  virtual void calculate_rt() const = 0;

  /// Retrieve the intensity value calculated
  virtual double get_intensity() const = 0;

  /// Copy jacobians out of internal xdata structures
  virtual void copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf_params, double& jac_surf_temp, blitz::Array<double, 1>& jac_atm_temp) const = 0;

protected:

  bool do_solar_sources, do_thermal_emission;

  /// Spurr BRDF class interface class to use
  mutable boost::shared_ptr<SpurrBrdfDriver> brdf_driver_;

};

}

#endif
