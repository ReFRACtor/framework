#include <nlls_solver_gsl_lmder.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void NLLSSolverGSLLMDER::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NLLSSolverGSL);
}

FP_IMPLEMENT(NLLSSolverGSLLMDER);
#endif


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
