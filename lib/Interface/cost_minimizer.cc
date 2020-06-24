#include <cost_minimizer.h>
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CostMinimizer::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IterativeSolver)
    & FP_NVP(P);
}

FP_IMPLEMENT(CostMinimizer);
#endif


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(CostMinimizer, IterativeSolver)
REGISTER_LUA_END()
#endif
