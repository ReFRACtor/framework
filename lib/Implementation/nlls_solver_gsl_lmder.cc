#include <nlls_solver_gsl_lmder.h>


using namespace FullPhysics;



boost::shared_ptr<IterativeSolver> nlls_solver_gsl_lmder_create(
	          const boost::shared_ptr<CostFunc>& NLLS,
                  int max_cost_function_calls,
                  double dx_tol_abs, double dx_tol_rel, double g_tol,
	          bool vrbs)
{
  const boost::shared_ptr<NLLSProblem> nlls(boost::dynamic_pointer_cast<NLLSProblem>(NLLS));
  return boost::shared_ptr<IterativeSolver>(new NLLSSolverGSLLMDER(nlls, max_cost_function_calls, dx_tol_abs, dx_tol_rel, g_tol, vrbs));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSSolverGSLLMDER, IterativeSolver)
.scope
[
 luabind::def("create", &nlls_solver_gsl_lmder_create)
]
REGISTER_LUA_END()
#endif
