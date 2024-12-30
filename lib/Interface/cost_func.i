// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "cost_func.h"
%}
%base_import(problem_state)

%fp_shared_ptr(FullPhysics::CostFunc);

namespace FullPhysics {
class CostFunc : virtual public ProblemState {
public:
  enum message_t {NONE, SOLVED, ERROR};
  CostFunc();
  virtual double cost() = 0;
  double cost_x(const blitz::Array<double, 1>& x);
  %python_attribute(num_cost_evaluations, int);
  void zero_num_evaluations();
  %python_attribute(message, message_t)
  %python_attribute(message_str, std::string)
  %pickle_serialization();
};
}
