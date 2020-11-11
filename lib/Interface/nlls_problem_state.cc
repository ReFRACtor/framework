#include <nlls_problem_state.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void NLLSProblemState::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ProblemState)
    & FP_NVP(R) & FP_NVP(J);
}

FP_IMPLEMENT(NLLSProblemState);
#endif

void NLLSProblemState::set(const NLLSProblemState& s)
{
  ProblemState::set(s);
  R.reference(s.R.copy());
  J.reference(s.J.copy());
}
