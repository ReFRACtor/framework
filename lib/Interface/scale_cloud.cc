#include "scale_cloud.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ScaleCloud::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableScaleCloud);
    
}

FP_IMPLEMENT(ScaleCloud);
FP_OBSERVER_SERIALIZE(ScaleCloud);
#endif
