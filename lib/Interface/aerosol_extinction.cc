#include "aerosol_extinction.h"
#include "fp_serialize_support.h"
#include <vector>
#include <boost/shared_ptr.hpp>

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AerosolExtinction::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableAerosolExtinction);
    
}

FP_IMPLEMENT(AerosolExtinction);
FP_OBSERVER_SERIALIZE(AerosolExtinction);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(AerosolExtinction)
REGISTER_LUA_END()

typedef std::vector<boost::shared_ptr<AerosolExtinction> >::reference 
(std::vector<boost::shared_ptr<AerosolExtinction> >::*vaevt)(std::vector<boost::shared_ptr<AerosolExtinction> >::size_type);

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<AerosolExtinction> >::*pbt1)(
        const std::vector<boost::shared_ptr<AerosolExtinction> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<AerosolExtinction> >, VectorAerosolExtinction)
.def(luabind::constructor<>())
.def("size", &std::vector<boost::shared_ptr<AerosolExtinction> >::size)
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<AerosolExtinction> >::push_back))
.def("value", ((vaevt) &std::vector<boost::shared_ptr<AerosolExtinction> >::operator[]))
REGISTER_LUA_END()
#endif
