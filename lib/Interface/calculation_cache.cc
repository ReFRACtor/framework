#include "calculation_cache.h"
#include "auto_derivative.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CalculationCache::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CacheInvalidatedObserver);
}

FP_IMPLEMENT(CalculationCache);
#endif
