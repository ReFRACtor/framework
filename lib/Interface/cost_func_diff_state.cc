#include <cost_func_diff_state.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CostFuncDiffState::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CostFuncState)
    & FP_NVP(G);
}

FP_IMPLEMENT(CostFuncDiffState);
#endif

void CostFuncDiffState::set(const CostFuncDiffState& s)
{
  CostFuncState::set(s);
  G.reference(s.G.copy());
}
