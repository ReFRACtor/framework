%include "fp_common.i"

%{
#include "max_a_posteriori_sqrt_constraint.h"
%}

%base_import(max_a_posteriori)
%base_import(model_measure_standard)

%import "forward_model.i"
%import "observation.i"
%import "state_vector.i"

%fp_shared_ptr(FullPhysics::MaxAPosterioriSqrtConstraint);

namespace FullPhysics {
class MaxAPosterioriSqrtConstraint : public MaxAPosteriori, public ModelMeasureStandard {
public:
  MaxAPosterioriSqrtConstraint(const boost::shared_ptr<ForwardModel>& fm,
          const boost::shared_ptr<Observation>& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint);
  MaxAPosterioriSqrtConstraint(const std::vector<boost::shared_ptr<ForwardModel> >& fm,
 	  const std::vector<boost::shared_ptr<Observation> >& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint);
  %pickle_serialization();
};
}
