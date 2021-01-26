#include "pressure_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(Pressure, SubStateVectorArrayPressure);

template<class Archive>
void PressureImpBaseCache::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CalculationCache);
}

template<class Archive>
void PressureImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayPressure)
    & FP_NVP(cache);
}

FP_IMPLEMENT(PressureImpBase);
FP_IMPLEMENT(PressureImpBaseCache);
#endif

void PressureImpBaseCache::fill_cache(const PressureImpBase& P)
{
  pgrid.value.resize_number_variable(P.coeff.number_variable());
  pgrid.value.jacobian() = 0;
  P.calc_pressure_grid();
}
