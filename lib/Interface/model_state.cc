#include <model_state.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ModelState::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ProblemState)
    & FP_NVP(K) & FP_NVP(M);
}

FP_IMPLEMENT(ModelState);
#endif

void ModelState::set(const ModelState& s)
{
  ProblemState::set(s);
  M.reference(s.M.copy());
  K.reference(s.K.copy());
}
