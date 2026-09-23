%include "fp_common.i"

%{
#include "multiscatt_brdf_driver.h"
#include "spurr_interface_types.h"
%}

%base_import(spurr_brdf_driver)

%import "spurr_interface_masters.i"

%fp_shared_ptr(FullPhysics::MultiScattBrdfDriver);

namespace FullPhysics {

class MultiScattBrdfDriver : public SpurrBrdfDriver {
public:
  %python_attribute_abstract(brdf_interface, boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base>)

  virtual void setup_geometry(double sza, double azm, double zen);

  %python_attribute(n_brdf_kernels, int)
  %python_attribute(n_kernel_factor_wfs, int)
  %python_attribute(n_kernel_params_wfs, int)
  %python_attribute(n_surface_wfs, int)
  %python_attribute(do_shadow_effect, bool)
  virtual bool do_kparams_derivs(const int kernel_index) const;

  %pickle_serialization();
};

}
