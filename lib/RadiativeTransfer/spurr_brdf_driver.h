#ifndef SPURR_BRDF_DRIVER_H
#define SPURR_BRDF_DRIVER_H

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

class SpurrBrdfDriver : public Printable<SpurrBrdfDriver> {

public:
  virtual void initialize_brdf_inputs(int surface_type);
  virtual void setup_geometry(double sza, double azm, double zen) = 0;
  virtual ArrayAd<double, 1> setup_brdf_inputs(int surface_type, const ArrayAd<double, 1>& surface_parameters);

  virtual int n_brdf_kernels() const = 0;

  virtual int n_kernel_factor_wfs() const = 0;
  virtual int n_kernel_params_wfs() const = 0;
  virtual int n_surface_wfs() const = 0;
  virtual bool do_kparams_derivs(const int kernel_index) const = 0;
  virtual bool do_shadow_effect() const = 0;
  virtual void print(std::ostream& Os) const {Os << "SpurrBrdfDriver";}
protected:
  virtual void initialize_brdf_kernel(int kernel_type);

  virtual void setup_lambertian_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const blitz::Array<int, 1>& parameter_indexes);
  virtual void setup_coxmunk_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const blitz::Array<int, 1>& parameter_indexes);
  virtual void setup_rahman_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const blitz::Array<int, 1>& parameter_indexes);
  virtual void setup_breon_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const blitz::Array<int, 1>& parameter_indexes);

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
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(SpurrBrdfDriver);

#endif
