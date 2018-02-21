// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_solver.h"
%}
%base_import(iterative_solver_der)
%import "nlls_problem.i"
%fp_shared_ptr(FullPhysics::NLLSSolver);

namespace FullPhysics {
class NLLSSolver : public IterativeSolverDer {
public:
  NLLSSolver(const boost::shared_ptr<NLLSProblem>& p, 
             int max_cost_function_calls, bool vrbs);
  virtual ~NLLSSolver();
};
}
