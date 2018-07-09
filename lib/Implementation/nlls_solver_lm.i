// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_solver_lm.h"
%}
%base_import(nlls_solver)
%import "nlls_problem.i"

%fp_shared_ptr(FullPhysics::NLLSSolverLM);

namespace FullPhysics {
  class NLLSSolverLM : public NLLSSolver {
  public:
    class Options {
    public:
    Options();
    double min_W;
    double tr_rad_tol;
    double tr_rad;
    double cr_ratio_tol;
    };
  NLLSSolverLM( const boost::shared_ptr<NLLSProblem>& p, int max_cost_function_calls, 
                const NLLSSolverLM::Options& opt=NLLSSolverLM::Options(),
                double dx_tol_abs=0.000001, double dx_tol_rel=0.000001,
                double g_tol_abs=0.000001, double g_tol_rel=6.0555e-06,
                bool vrbs=false );
  virtual ~NLLSSolverLM();
  virtual void solve();
};
}
