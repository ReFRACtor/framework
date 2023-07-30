#include "observation.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Observation::serialize(Archive & ar,
			const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StackedRadianceMixin);
  // Older version wasn't an observable. If we don't load anything,
  // then none of the observable stuff is set up, which is actually
  // what we want to older serialization
  if(version > 0)
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableObservation);
    
}

FP_IMPLEMENT(Observation);
FP_OBSERVER_SERIALIZE(Observation);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Observation)
REGISTER_LUA_END()
#endif
