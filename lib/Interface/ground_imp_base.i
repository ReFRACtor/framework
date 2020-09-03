%include "fp_common.i"

%{
#include "ground_imp_base.h"
%}

%base_import(state_vector)
%base_import(sub_state_vector_array)
%base_import(ground)
%base_import(state_mapping)
%base_import(state_mapping_linear)

%fp_shared_ptr(FullPhysics::GroundImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::Ground>);
namespace FullPhysics {
%template(SubStateVectorArrayGround) 
     FullPhysics::SubStateVectorArray<Ground>;

class GroundImpBase: public SubStateVectorArray<Ground> {
public:
  virtual boost::shared_ptr<Ground> clone() const = 0;
  %python_attribute(state_used, blitz::Array<bool, 1>)
  %sub_state_virtual_func(Ground);
  %pickle_serialization();
protected:
  GroundImpBase(const blitz::Array<double, 1>& Coeff,
                boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
};
}

