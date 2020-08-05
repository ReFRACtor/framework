// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "calculation_cache.h"
%}
%base_import(observer)

%fp_shared_ptr(FullPhysics::CalculationCache);

namespace FullPhysics {
class CalculationCache : public CacheInvalidatedObserver {
public:
  CalculationCache();
  void fill_cache_if_needed();
  virtual void fill_cache() = 0;
};
}
