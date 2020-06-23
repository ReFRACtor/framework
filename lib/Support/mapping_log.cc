#include "mapping_log.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void MappingLog::serialize(Archive& ar,
			      const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Mapping)
    & FP_NVP(map_name);
}

FP_IMPLEMENT(MappingLog);
#endif

