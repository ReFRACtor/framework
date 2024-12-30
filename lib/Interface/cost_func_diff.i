// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "cost_func_diff.h"
%}
%base_import(cost_func)

%fp_shared_ptr(FullPhysics::CostFuncDiff);

namespace FullPhysics {
class CostFuncDiff : public CostFunc {
public:
  CostFuncDiff();
  %python_attribute_nonconst(gradient, blitz::Array<double, 1>);
  blitz::Array<double, 1> gradient_x(const blitz::Array<double, 1>& x);
  void cost_gradient(
    double& OUTPUT, blitz::Array<double, 1>& OUTPUT);
  void cost_gradient_x(const blitz::Array<double, 1>& x,
    double& OUTPUT, blitz::Array<double, 1>& OUTPUT);
  %python_attribute(num_der1_evaluations, int);
  void zero_num_evaluations();
  %python_attribute(gradient_size, int);
  %pickle_serialization();
};
}
