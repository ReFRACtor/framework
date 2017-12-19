%include "common.i"

%{
#include "model_measure_standard.h"
%}

%base_import(model_measure)

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
protected:
  ModelMeasureStandard();
};
}
