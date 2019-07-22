// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "meyer_ml_problem.h"
%}
%base_import(max_likelihood)
%base_import(model_measure_meyer)
%fp_shared_ptr(FullPhysics::MeyerMLProblem);

namespace FullPhysics {
class MeyerMLProblem : public MaxLikelihood, public ModelMeasureMeyer {
public:
  MeyerMLProblem(const blitz::Array<double, 1>& measurement, 
                 const blitz::Array<double, 1>& measurement_error_cov);
};
}
