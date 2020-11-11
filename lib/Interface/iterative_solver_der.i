// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "iterative_solver_der.h"
%}

%base_import(iterative_solver)
%fp_shared_ptr(FullPhysics::IterativeSolverDer);

namespace FullPhysics {
class IterativeSolverDer : public IterativeSolver {
public:
  IterativeSolverDer(int max_cost_function_calls, bool vrbs);
  virtual ~IterativeSolverDer();
  %python_attribute(gradient_at_accepted_points, std::vector< blitz::Array<double, 1> >)
  %pickle_serialization();
};
}
