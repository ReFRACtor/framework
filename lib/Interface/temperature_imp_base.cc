#include "temperature_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(Temperature, SubStateVectorArrayTemperature);

template<class Archive>
void TemperatureImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayTemperature)
    & FP_NVP(press);
}

FP_IMPLEMENT(TemperatureImpBase);
#endif
