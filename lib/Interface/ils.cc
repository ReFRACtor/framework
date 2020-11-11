#include "ils.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Ils::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableIls);
    
}

FP_IMPLEMENT(Ils);
FP_OBSERVER_SERIALIZE(Ils);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Ils)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<Ils> >::*pbt1)(
        const std::vector<boost::shared_ptr<Ils> >::value_type&);
REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<Ils> >, VectorIls)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<Ils> >::push_back))
REGISTER_LUA_END()
#endif
