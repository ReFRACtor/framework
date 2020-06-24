#include <nlls_solver.h>
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void NLLSSolver::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IterativeSolverDer)
    & FP_NVP(P);
}

FP_IMPLEMENT(NLLSSolver);
#endif


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSSolver, IterativeSolver)
REGISTER_LUA_END()
#endif
