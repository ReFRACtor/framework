%include "fp_common.i"

%{
#include "state_mapping.h"
%}

%base_import(generic_object)
%import "array_ad.i"
%import "pressure.i"
%template(vector_state_mapping) std::vector<boost::shared_ptr<FullPhysics::StateMapping> >;
%fp_shared_ptr(FullPhysics::StateMapping);

namespace FullPhysics {
%feature("director") StateMapping;

class StateMapping : public GenericObject {
public:
  virtual ~StateMapping();
  virtual ArrayAd<double, 1> mapped_state(const ArrayAd<double, 1>& updated_coeff) const;
  virtual ArrayAd<double, 1> retrieval_state(const ArrayAd<double, 1>& initial_values) const;
  virtual blitz::Array<double, 2> jacobian_retrieval
  (const blitz::Array<double, 1>& x,
   const blitz::Array<double, 2>& jacobian_mapped) const;
  virtual int state_vector_name_index(const int retrieval_state_index) const;
  %python_attribute(name, std::string);
  virtual boost::shared_ptr<StateMapping> clone() const = 0;
  std::string print_to_string();
  std::string print_parent() const;
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(state_mapping, StateMapping)

// List of things "import *" will include
%python_export("StateMapping");
