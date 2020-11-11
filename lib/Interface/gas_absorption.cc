#include "gas_absorption.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void GasAbsorption::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(GasAbsorption);

  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
}

FP_IMPLEMENT(GasAbsorption);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(GasAbsorption)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<GasAbsorption> >::*pbt1)(
        const std::vector<boost::shared_ptr<GasAbsorption> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<GasAbsorption> >, VectorGasAbsorption)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<GasAbsorption> >::push_back))
REGISTER_LUA_END()
#endif
