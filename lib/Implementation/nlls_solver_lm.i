// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "nlls_solver_lm.h"
%}
%base_import(nlls_solver)
%import "nlls_problem.i"

%fp_shared_ptr(FullPhysics::NLLSSolverLM);
%fp_shared_ptr(FullPhysics::NLLSSolverLMOptions);

namespace FullPhysics {

// Pull Options out of NLLSolver so its not a nested class and SWIG will create wrapper code
class NLLSSolverLMOptions {
    public:
    NLLSSolverLMOptions();
    double min_W;
    double tr_rad_tol;
    double tr_rad;
    double cr_ratio_tol;
};

class NLLSSolverLM : public NLLSSolver {
  public:
  NLLSSolverLM( const boost::shared_ptr<NLLSProblem>& p, int max_cost_function_calls, 
                const NLLSSolverLMOptions& opt=NLLSSolverLMOptions(),
                double dx_tol_abs=0.000001, double dx_tol_rel=0.000001,
                double g_tol_abs=0.000001, double g_tol_rel=6.0555e-06,
                bool vrbs=false );
  virtual ~NLLSSolverLM();
  virtual void solve();
};

}

%{
namespace FullPhysics {
    // Make C++ think that NLLSSolver::Options is a global class
    typedef NLLSSolverLM::Options NLLSSolverLMOptions;
}
%}
