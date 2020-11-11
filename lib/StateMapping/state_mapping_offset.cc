#include "state_mapping_offset.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateMappingOffset::serialize(Archive& ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateMapping)
    & FP_NVP(map_name) & FP_NVP_(initial_offset) & FP_NVP_(offsetee);
}

FP_IMPLEMENT(StateMappingOffset);
#endif

