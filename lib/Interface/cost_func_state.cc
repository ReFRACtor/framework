#include <cost_func_state.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CostFuncState::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ProblemState)
    & FP_NVP(C);
}

FP_IMPLEMENT(CostFuncState);
#endif

void CostFuncState::set(const CostFuncState& s)
{
  ProblemState::set(s);
  C.reference(s.C.copy());
}
