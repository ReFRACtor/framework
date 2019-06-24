// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "cost_minimizer.h"
%}
%base_import(iterative_solver);
%import "cost_func.i"
%fp_shared_ptr(FullPhysics::CostMinimizer);

namespace FullPhysics {
class CostMinimizer : public IterativeSolver {
public:
  CostMinimizer(const boost::shared_ptr<CostFunc>& p,
                int max_cost_function_calls, bool vrbs);
  virtual ~CostMinimizer();
};
}
