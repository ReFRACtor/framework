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
  NLLSSolverGSLSM( const boost::shared_ptr<NLLSProblem>& p, int max_cost_function_calls, 
                   double dx_tol_abs=0.000001, double dx_tol_rel=0.000001, double g_tol=6.0555e-06, 
                   gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters(),
                   double x_tol=1.0e-6, double g_tol=6.0555e-06, double f_tol=0.0, bool vrbs=false );
  virtual ~NLLSSolverGSLSM();
  virtual void solve();
};
}
