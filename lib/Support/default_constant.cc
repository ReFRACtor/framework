#include "default_constant.h"

using namespace FullPhysics;

#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void DefaultConstant::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Constant);
}

FP_IMPLEMENT(DefaultConstant);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(DefaultConstant, Constant)
.def(luabind::constructor<>())
REGISTER_LUA_END()
#endif
