#include "observation_sv.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void ObservationSv::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableObservationSv)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Observation);
}

FP_IMPLEMENT(ObservationSv);
FP_OBSERVER_SERIALIZE(ObservationSv);
#endif

