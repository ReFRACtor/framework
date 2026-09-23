%include "fp_common.i"

%{
#include "vlidort_brdf_driver.h"
#include "spurr_interface_types.h"
%}

%base_import(multiscatt_brdf_driver)

%import "spurr_interface_masters.i"
%import "vlidort_interface_masters.i"
%import "vlidort_interface_types.i"

%fp_shared_ptr(FullPhysics::VLidortBrdfDriver);

namespace FullPhysics {

class VLidortBrdfDriver : public MultiScattBrdfDriver {
public:
  VLidortBrdfDriver(int nstream, int nmoment);

  %python_attribute(brdf_interface, boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base>)
  %python_attribute(vlidort_brdf_interface, boost::shared_ptr<VBrdf_Linsup_Masters>)

  %pickle_serialization();
};

}
