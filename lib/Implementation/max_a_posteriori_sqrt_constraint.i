%include "fp_common.i"

%{
#include "max_a_posteriori_sqrt_constraint.h"
%}

%base_import(max_a_posteriori)
%base_import(model_measure_standard)

%import "forward_model.i"
%import "observation.i"
%import "state_vector.i"
%import "state_mapping.i"

%fp_shared_ptr(FullPhysics::MaxAPosterioriSqrtConstraint);

namespace FullPhysics {
class MaxAPosterioriSqrtConstraint : public MaxAPosteriori, public ModelMeasureStandard {
public:
  MaxAPosterioriSqrtConstraint(const boost::shared_ptr<ForwardModel>& fm,
          const boost::shared_ptr<Observation>& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint,
          const boost::shared_ptr<StateMapping>& in_map =
			       boost::make_shared<StateMappingLinear>());
  MaxAPosterioriSqrtConstraint(const std::vector<boost::shared_ptr<ForwardModel> >& fm,
 	  const std::vector<boost::shared_ptr<Observation> >& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint,
          const boost::shared_ptr<StateMapping>& in_map = boost::make_shared<StateMappingLinear>());
  %python_attribute(mapping, boost::shared_ptr<StateMapping>);
  %python_attribute_nonconst(jacobian_fm, blitz::Array<double,2>);
  %python_attribute_nonconst(measurement_jacobian_fm, blitz::Array<double,2>);
  %python_attribute_nonconst(model_measure_diff_jacobian_fm, blitz::Array<double,2>);
  %pickle_serialization();
};
}
