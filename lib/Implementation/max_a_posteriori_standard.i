%include "fp_common.i"

%{
#include "max_a_posteriori_standard.h"
%}

%base_import(max_a_posteriori)
%base_import(model_measure_standard)

%import "forward_model.i"
%import "observation.i"
%import "state_vector.i"

%fp_shared_ptr(FullPhysics::MaxAPosterioriStandard);

namespace FullPhysics {
class MaxAPosterioriStandard : public MaxAPosteriori, public ModelMeasureStandard {
public:
  MaxAPosterioriStandard(const boost::shared_ptr<ForwardModel>& fm,
          const boost::shared_ptr<Observation>& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> a_priori_cov);

  virtual ~MaxAPosterioriStandard();

};
}
