#include <iterative_solver_der.h>
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void IterativeSolverDer::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IterativeSolver)
    & FP_NVP(Gradient_at_accepted_points);
}

FP_IMPLEMENT(IterativeSolverDer);
#endif


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(IterativeSolverDer, IterativeSolver)
REGISTER_LUA_END()
#endif


