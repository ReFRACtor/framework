// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "model_measure_meyer.h"
%}
%base_import(model_measure)
%fp_shared_ptr(FullPhysics::ModelMeasureMeyer);

namespace FullPhysics {
class ModelMeasureMeyer : virtual public ModelMeasure {
public:
  virtual void model_eval();
  virtual void jacobian_eval();
};
}
