%include "fp_common.i"

%{
#include "lidort_brdf_driver.h"
%}

%base_import(spurr_brdf_driver)

%import "lidort_interface_masters.i"
%import "lidort_interface_types.i"

%fp_shared_ptr(FullPhysics::LidortBrdfDriver);

namespace FullPhysics {

class LidortBrdfDriver : public SpurrBrdfDriver {
public:
  LidortBrdfDriver(int nstream, int nmoment);
  virtual ~LidortBrdfDriver();

  %python_attribute(brdf_interface, boost::shared_ptr<Brdf_Lin_Sup_Masters>)

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
