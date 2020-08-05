#ifndef CALCULATION_CACHE_H
#define CALCULATION_CACHE_H
#include "observer.h"

namespace FullPhysics {
/****************************************************************//**
  General class for caching a calculation. We do the calculation each
  time a value is needed, if the cache has been marked as invalidated.
*******************************************************************/

class CalculationCache : public CacheInvalidatedObserver {
public:
  CalculationCache() {}
  virtual ~CalculationCache() {}
  void fill_cache_if_needed()
  {
    if(cache_valid)
      return;
    fill_cache();
    cache_valid = true;
  }
  virtual void fill_cache() = 0;
  void print(std::ostream& Os) const
  {
    Os << "CalculationCache";
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(CalculationCache);
#endif
