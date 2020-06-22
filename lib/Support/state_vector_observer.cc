#include "state_vector_observer.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateVectorObserver::serialize(Archive & UNUSED(ar),
				    const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(StateVectorObserver);
  // Leave out the Observer<StateVector> part, I'm not sure if we
  // want to serialize that or not
}

FP_IMPLEMENT(StateVectorObserver);
#endif

