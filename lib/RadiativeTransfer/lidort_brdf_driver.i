%include "fp_common.i"

%{
#include "lidort_brdf_driver.h"
#include "spurr_interface_types.h"
%}

%base_import(multiscatt_brdf_driver)

%import "spurr_interface_masters.i"
%import "lidort_interface_masters.i"
%import "lidort_interface_types.i"

%fp_shared_ptr(FullPhysics::LidortBrdfDriver);

namespace FullPhysics {

class LidortBrdfDriver : public MultiScattBrdfDriver {
public:
  LidortBrdfDriver(int nstream, int nmoment);

  %python_attribute(brdf_interface, boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base>)
  %python_attribute(lidort_brdf_interface, boost::shared_ptr<Brdf_Lin_Sup_Masters>)

  %pickle_serialization();
};

}
