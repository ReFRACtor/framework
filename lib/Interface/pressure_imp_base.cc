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
				const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayPressure)
    & FP_NVP(cache);
  // Older version didn't have the type_preference_, instead it
  // was always PREFER_INCREASING_PRESSURE
  if(version > 0)
    ar & FP_NVP_(type_preference);
}

FP_IMPLEMENT(PressureImpBase);
FP_IMPLEMENT(PressureImpBaseCache);
#endif

ArrayAdWithUnit<double, 1> PressureImpBase::pressure_grid
(Pressure::PressureGridType Gtype) const
{
  cache.fill_cache_if_needed(*this);
  if(Gtype == Pressure::NATIVE_ORDER ||
     (int) Gtype == (int) type_preference_)
    return cache.pgrid;
  return ArrayAdWithUnit<double, 1>
    (ArrayAd<double, 1>(cache.pgrid.value.value().reverse(blitz::firstDim),
			cache.pgrid.value.jacobian().reverse(blitz::firstDim),
			cache.pgrid.value.is_constant()),
     cache.pgrid.units);
}


void PressureImpBaseCache::fill_cache(const PressureImpBase& P)
{
  pgrid.value.resize_number_variable(P.coeff.number_variable());
  pgrid.value.jacobian() = 0;
  P.calc_pressure_grid();
}
