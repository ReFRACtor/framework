#include "state_vector_observer.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateVectorObserver::serialize(Archive & ar,
				    const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverStateVector);
}

FP_IMPLEMENT(StateVectorObserver);
#endif

