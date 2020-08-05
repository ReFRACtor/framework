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

/****************************************************************//**
  Variation of CalculationCache that also takes an index (e.g., a
  spectral index)
*******************************************************************/

template<class T> class CalculationCacheIndexed :
    public CacheInvalidatedObserver {
public:
  CalculationCacheIndexed() {}
  virtual ~CalculationCacheIndexed() {}
  void fill_cache_if_needed(const T& Owner_Obj, int Index)
  {
    if(cache_valid)
      return;
    fill_cache(Owner_Obj, Index);
    cache_valid = true;
  }
  virtual void fill_cache(const T& Owner_Obj, int Index) = 0;
private:
  // Derived classes can skip this in serialization, and
  // go straight to CacheInvalidatedObserver
};
}

#endif
