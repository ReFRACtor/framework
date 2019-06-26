// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "bard_nlls_problem.h"
%}
%base_import(nlls_problem)
%base_import(nlls_problem_state)
%fp_shared_ptr(FullPhysics::BardNLLSProblem);

namespace FullPhysics {
class BardNLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  BardNLLSProblem();
  %python_attribute(residual_size, int);
  %python_attribute(expected_parameter_size, int);
  %python_attribute_nonconst(residual, blitz::Array<double, 1>);
  %python_attribute_nonconst(jacobian, blitz::Array<double, 2>);
};
}
