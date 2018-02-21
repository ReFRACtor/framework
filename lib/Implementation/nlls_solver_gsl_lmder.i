// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_solver_gsl_lmder.h"
%}
%base_import(nlls_solver_gsl)
%fp_shared_ptr(FullPhysics::NLLSSolverGSLLMDER);
%import "nlls_problem.i"

%fp_shared_ptr(FullPhysics::NLLSSolverGSLLMDER);

namespace FullPhysics {
  class NLLSSolverGSLLMDER : public NLLSSolverGSL {
public:
  NLLSSolverGSLLMDER(const boost::shared_ptr<NLLSProblem>& p, int max_cost_function_calls,
                     double dx_tol_abs=0.000001, double dx_tol_rel=0.000001, double g_tol=6.0555e-06, 
                     bool vrbs=false);
  virtual ~NLLSSolverGSLLMDER();
};
}
