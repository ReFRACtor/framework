// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "bard_ml_problem.h"
%}
%base_import(max_likelihood)
%base_import(model_measure_bard)
%fp_shared_ptr(FullPhysics::BardMLProblem);

namespace FullPhysics {
class BardMLProblem : public MaxLikelihood, public ModelMeasureBard {
public:
  BardMLProblem(const blitz::Array<double, 1>& measurement, 
		const blitz::Array<double, 1>& measurement_error_cov);
};
}
