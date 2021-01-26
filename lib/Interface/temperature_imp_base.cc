#include "temperature_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(Temperature, SubStateVectorArrayTemperature);

template<class Archive>
void TemperatureImpBaseCache::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CalculationCache);
}

template<class Archive>
void TemperatureImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayTemperature)
    & FP_NVP(cache)
    & FP_NVP(press);
}

FP_IMPLEMENT(TemperatureImpBase);
FP_IMPLEMENT(TemperatureImpBaseCache);
#endif

void TemperatureImpBaseCache::fill_cache(const TemperatureImpBase& T)
{
  T.calc_temperature_grid();
}
