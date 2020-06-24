#include "cost_function.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CostFunction::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(CostFunction);

  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
}

FP_IMPLEMENT(CostFunction);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(CostFunction)
REGISTER_LUA_END()
#endif
