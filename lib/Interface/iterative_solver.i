// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "iterative_solver.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::IterativeSolver);

namespace FullPhysics {
class IterativeSolver : public GenericObject {
public:
  enum status_t {SUCCESS, CONTINUE, STALLED, ERROR, UNTRIED};
  IterativeSolver(int max_cost_function_calls, bool vrbs);
  virtual ~IterativeSolver();
  %python_attribute(num_accepted_steps, int)
  %python_attribute(accepted_points, std::vector< blitz::Array<double, 1> >)
  %python_attribute(cost_at_accepted_points, std::vector<double>)
  virtual void solve() = 0;
  %python_attribute(status, status_t)
  %python_attribute(status_str, std::string)
  std::string print_to_string() const;
};
}
