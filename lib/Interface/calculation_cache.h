#ifndef CALCULATION_CACHE_H
#define CALCULATION_CACHE_H
#include "observer.h"

namespace FullPhysics {
/****************************************************************//**
  General class for caching a calculation. We do the calculation each
  time a value is needed, if the cache has been marked as invalidated.
*******************************************************************/

template<class T> class CalculationCache : public CacheInvalidatedObserver {
public:
  CalculationCache() {}
  virtual ~CalculationCache() {}
  void fill_cache_if_needed(const T& Owner_Obj)
  {
    if(cache_valid)
      return;
    fill_cache(Owner_Obj);
    cache_valid = true;
  }
  virtual void fill_cache(const T& Owner_Obj) = 0;
private:
  // Derived classes can skip this in serialization, and
  // go straight to CacheInvalidatedObserver
};
}

#endif
