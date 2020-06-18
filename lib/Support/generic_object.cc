#include "generic_object.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void GenericObject::serialize(Archive & UNUSED(ar),
			      const unsigned int UNUSED(version))
{
    // Nothing to do
}

FP_IMPLEMENT(GenericObject);
#endif
