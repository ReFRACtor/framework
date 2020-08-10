%include "fp_common.i"

%{
#include "state_mapping.h"
%}

%base_import(generic_object)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMapping);

namespace FullPhysics {

class StateMapping : public GenericObject {
public:
  virtual ~StateMapping();
  virtual ArrayAd<double, 1> fm_view(const ArrayAd<double, 1>& updated_coeff) const;
  virtual ArrayAd<double, 1> retrieval_init
  (const ArrayAd<double, 1>& initial_coeff) const;
  %python_attribute(name, std::string);
  virtual boost::shared_ptr<StateMapping> clone() = 0;
  std::string print_to_string();
  %pickle_serialization();
};
}

%template(vector_state_mapping) std::vector<boost::shared_ptr<FullPhysics::StateMapping> >;
