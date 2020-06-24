#include <cost_func_diff.h>
#include <fp_exception.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CostFuncDiff::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CostFunc)
    & FP_NVP(d_count);
}

FP_IMPLEMENT(CostFuncDiff);
#endif


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(CostFuncDiff, CostFunc)
REGISTER_LUA_END()
#endif


void CostFuncDiff::cost_gradient(double& c, Array<double, 1>& g)
{
  c = cost();
  g.reference(gradient());
}
