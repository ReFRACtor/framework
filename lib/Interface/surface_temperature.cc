#include "surface_temperature.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SurfaceTemperature::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableSurfaceTemperature);
    
}

FP_IMPLEMENT(SurfaceTemperature);
FP_OBSERVER_SERIALIZE(SurfaceTemperature);
#endif
