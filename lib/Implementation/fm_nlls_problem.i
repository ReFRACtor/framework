// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "fm_nlls_problem.h"
%}
%import "forward_model.i"
%import "state_vector.i"
%fp_shared_ptr(FullPhysics::FmNLLSProblem);

namespace FullPhysics {
class FmNLLSProblem {
public:
  FmNLLSProblem(const boost::shared_ptr<ForwardModel>& Fm,
   const boost::shared_ptr<StateVector>& Sv,
   const blitz::Array<double, 1>& Rad,
   const blitz::Array<double, 1>& Rad_uncer,
   const blitz::Array<double, 1> X_apriori,
   const blitz::Array<double, 2> Apriori_cov);
  virtual ~FmNLLSProblem();
  %python_attribute(residual_size, int)
  %python_attribute(parameter_size, int)
  %python_attribute_nonconst(residual, blitz::Array<double, 1>)
  %python_attribute_nonconst(jacobian, blitz::Array<double, 2>)
};
}
