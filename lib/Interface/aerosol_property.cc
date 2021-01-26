#include "aerosol_property.h"
#include "fp_serialize_support.h"
#include <vector>
#include <boost/shared_ptr.hpp>

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AerosolProperty::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableAerosolProperty);
    
}

FP_IMPLEMENT(AerosolProperty);
FP_OBSERVER_SERIALIZE(AerosolProperty);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(AerosolProperty)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<AerosolProperty> >::*pbt1)(
        const std::vector<boost::shared_ptr<AerosolProperty> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<AerosolProperty> >, VectorAerosolProperty)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<AerosolProperty> >::push_back))
REGISTER_LUA_END()

#endif
