// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "model_measure_bard.h"
%}
%base_import(model_measure)
%fp_shared_ptr(FullPhysics::ModelMeasureBard);

namespace FullPhysics {
class ModelMeasureBard : virtual public ModelMeasure {
public:
  virtual void model_eval();
  virtual void jacobian_eval();
};
}
