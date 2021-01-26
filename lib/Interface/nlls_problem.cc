#include "nlls_problem.h"
#include "fp_exception.h"
#include "fp_serialize_support.h"


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void NLLSProblem::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CostFuncDiff);
}

FP_IMPLEMENT(NLLSProblem);
#endif


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSProblem, CostFunc)
REGISTER_LUA_END()
#endif



void NLLSProblem::residual_jacobian(Array<double, 1>& r, Array<double, 2>& j)
{
  r.reference(residual());
  j.reference(jacobian());
}

double NLLSProblem::cost()
{
  return sum(residual()*residual())/2.0;
}

blitz::Array<double, 1> NLLSProblem::gradient()
{
  firstIndex i1; secondIndex i2;
  Array<double, 1> ret(sum(jacobian()(i2,i1) * residual()(i2), i2));
  return ret;
}
