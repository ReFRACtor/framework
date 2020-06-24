#include "max_likelihood.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void MaxLikelihood::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelMeasure);
}

FP_IMPLEMENT(MaxLikelihood);
#endif


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(MaxLikelihood)
REGISTER_LUA_END()
#endif

