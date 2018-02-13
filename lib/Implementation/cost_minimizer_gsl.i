// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "cost_minimizer_gsl.h"
%}
%base_import(cost_minimizer)
%import "cost_func.i"

%fp_shared_ptr(FullPhysics::CostMinimizerGSL);

namespace FullPhysics {
  class CostMinimizerGSL : public CostMinimizer {
  public:
  CostMinimizerGSL(const boost::shared_ptr<CostFunc>& p,
                   int max_cost_function_calls, double size_tol=0.001, 
                   const blitz::Array<double,1>& init_step_size=blitz::Array<double,1>(),
                   bool vrbs=false);
  virtual ~CostMinimizerGSL();
  virtual void solve();
};
}
