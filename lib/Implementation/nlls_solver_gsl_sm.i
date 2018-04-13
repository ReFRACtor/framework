// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_solver_gsl_sm.h"
%}
%base_import(nlls_solver)
%import "nlls_problem.i"

%fp_shared_ptr(FullPhysics::NLLSSolverGSLSM);

namespace FullPhysics {
  class NLLSSolverGSLSM : public NLLSSolver {
  public:
  NLLSSolverGSLSM( const boost::shared_ptr<NLLSProblem>& p, int32_t max_cost_function_calls,
                   gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters(),
                   double x_tol=1.0e-6, double g_tol=6.0555e-06, double f_tol=0.0, bool vrbs=false );
  virtual ~NLLSSolverGSLSM();
  virtual void solve();
};
}
