// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "nlls_max_a_posteriori.h"
%}
%base_import(nlls_problem)
%base_import(nlls_problem_state)
%import "max_a_posteriori.i"
%fp_shared_ptr(FullPhysics::NLLSMaxAPosteriori);

namespace FullPhysics {
// Allow these classes to be derived from in Python.
%feature("director") NLLSMaxAPosteriori;
  
class NLLSMaxAPosteriori: public NLLSProblem, public NLLSProblemState {
public:
  NLLSMaxAPosteriori(const boost::shared_ptr<MaxAPosteriori>& map);
  virtual double cost();
  void cost_gradient(
    double& c, blitz::Array<double, 1>& g);
  void cost_gradient_x(const blitz::Array<double, 1>& x,
    double& c, blitz::Array<double, 1>& g);
  %python_attribute_derived_nonconst(residual, blitz::Array<double, 1>);
  virtual blitz::Array<double, 1> residual_x(const blitz::Array<double, 1>& x);
  %python_attribute_derived_nonconst(jacobian, blitz::Array<double, 2>);
  virtual blitz::Array<double, 2> jacobian_x(const blitz::Array<double, 1>& x);
  virtual void residual_jacobian
  (blitz::Array<double, 1>& r, blitz::Array<double, 2>& j);
  virtual void residual_jacobian_x
  (const blitz::Array<double, 1>& x,
   blitz::Array<double, 1>& r, blitz::Array<double, 2>& j);
  virtual bool parameters_different(const blitz::Array<double, 1>& x) const;
  %python_attribute(num_residual_evaluations, int);
  %python_attribute(num_jacobian_evaluations, int);
  %python_attribute(residual_size, int);
  %python_attribute(parameter_size, int);
  %python_attribute_derived(expected_parameter_size, int);
  virtual void assert_parameter_set_correctly() const;
  virtual void assert_parameter_correct(const blitz::Array<double, 1>& x) const;
  %python_attribute_with_set(parameters, blitz::Array<double, 1>);
  %python_attribute_nonconst(max_a_posteriori, boost::shared_ptr<MaxAPosteriori>);
  %pickle_serialization();
protected:
  NLLSMaxAPosteriori();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(nlls_max_a_posteriori, NLLSMaxAPosteriori)

// List of things "import *" will include
%python_export("NLLSMaxAPosteriori");
