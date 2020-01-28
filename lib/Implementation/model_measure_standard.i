// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "model_measure_standard.h"
%}

%base_import(model_measure)
%import "forward_model.i"
%import "observation.i"
%import "state_vector.i"

%fp_shared_ptr(FullPhysics::ModelMeasureStandard);

namespace FullPhysics {
class ModelMeasureStandard : virtual public ModelMeasure {
public:
  virtual ~ModelMeasureStandard();
  virtual void model_eval();
  virtual void jacobian_eval();
  virtual void model_jacobian_eval();
  %python_attribute(expected_parameter_size, int)
  %python_attribute_with_set(parameters,blitz::Array<double, 1>)
  %python_attribute(forward_model, std::vector<boost::shared_ptr<ForwardModel> >)
  %python_attribute(observation, std::vector<boost::shared_ptr<Observation> >)
  %python_attribute(state_vector, boost::shared_ptr<StateVector>)
protected:
  ModelMeasureStandard();
};
}
