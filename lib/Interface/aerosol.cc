#include "aerosol.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Aerosol::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableAerosol);
    
}

FP_IMPLEMENT(Aerosol);
FP_OBSERVER_SERIALIZE(Aerosol);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Aerosol)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Timer for Aerosol
//-----------------------------------------------------------------------

AccumulatedTimer Aerosol::timer("Aerosol");

