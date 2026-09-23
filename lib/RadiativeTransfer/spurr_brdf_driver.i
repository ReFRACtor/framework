%include "fp_common.i"

%{
#include "spurr_brdf_driver.h"
%}

%base_import(generic_object)

%import "array_ad.i"

%fp_shared_ptr(FullPhysics::SpurrBrdfDriver);

namespace FullPhysics {

class SpurrBrdfDriver : public GenericObject {
public:
  std::string print_to_string() const;
  virtual void initialize_brdf_inputs(int surface_type);
  virtual void setup_geometry(double sza, double azm, double zen) = 0;
  virtual ArrayAd<double, 1> setup_brdf_inputs(int surface_type, const ArrayAd<double, 1>& surface_parameters);

  %python_attribute_abstract(n_brdf_kernels, int)
  %python_attribute_abstract(n_kernel_factor_wfs, int)
  %python_attribute_abstract(n_kernel_params_wfs, int)
  %python_attribute_abstract(n_surface_wfs, int)
  %python_attribute_abstract(do_shadow_effect, bool)
  %pickle_serialization();
};

}
