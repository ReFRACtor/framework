#include "state_mapping_composite.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateMappingComposite::serialize(Archive& ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateMapping)
    & FP_NVP(mappings);
}

FP_IMPLEMENT(StateMappingComposite);
#endif

