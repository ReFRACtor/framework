#include "measured_radiance_field.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void MeasuredRadianceField::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableMeasuredRadianceField);
}

FP_IMPLEMENT(MeasuredRadianceField);
FP_OBSERVER_SERIALIZE(MeasuredRadianceField);
#endif

