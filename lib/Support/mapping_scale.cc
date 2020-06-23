#include "mapping_scale.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void MappingScale::serialize(Archive& ar,
			      const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Mapping)
    & FP_NVP(map_name) & FP_NVP_(initial_scale_factor) & FP_NVP_(scalee);
}

FP_IMPLEMENT(MappingScale);
#endif

