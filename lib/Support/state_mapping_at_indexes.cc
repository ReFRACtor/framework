#include "state_mapping_at_indexes.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateMappingAtIndexes::serialize(Archive& ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateMapping) &
      FP_NVP(full_state) & FP_NVP(retrieval_indexes);
}

FP_IMPLEMENT(StateMappingAtIndexes);
#endif

