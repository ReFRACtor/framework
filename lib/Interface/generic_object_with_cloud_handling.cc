#include "generic_object_with_cloud_handling.h"
#include "fp_serialize_support.h"
#include <fp_exception.h>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void GenericObjectWithCloudHandling::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(GenericObjectWithCloudHandling);
  ar & FP_NVP_(do_cloud);
}

FP_IMPLEMENT(GenericObjectWithCloudHandling);
#endif
