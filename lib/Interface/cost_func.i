// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
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
  virtual ~CostFunc();
  virtual double cost() = 0;
  virtual double cost_x(const blitz::Array<double, 1>& x);
  %python_attribute(num_cost_evaluations, int);
  virtual void zero_num_evaluations();
  %python_attribute(message, message_t)
  %python_attribute(message_str, std::string)
};
}
