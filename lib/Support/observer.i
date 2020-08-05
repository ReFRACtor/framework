// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "observer.h"
%}
%base_import(generic_object)

%fp_shared_ptr(FullPhysics::CacheInvalidatedObserver);

namespace FullPhysics {
template<class T> class Observable;
template<class T> class Observer : public GenericObject {
public:
  Observer();
  virtual ~Observer();
  virtual void notify_update(const T& Observed_object);
  virtual void notify_add(T& Observed_object);
  virtual void notify_remove(T& Observed_object);
};

class CacheInvalidatedObserver: public GenericObject
{
public:
  CacheInvalidatedObserver();
  virtual void invalidate_cache();
  std::string print_to_string() const;
  %pickle_serialization();
protected:
  bool cache_valid;
};
  
template<class T> class Observable : public GenericObject {
public:
  virtual ~Observable();
  void add_observer_and_keep_reference(boost::shared_ptr<Observer<T> >& Obs);
  void add_cache_invalidated_observer(CacheInvalidatedObserver& Obs);
  void remove_cache_invalidated_observer(CacheInvalidatedObserver& Obs);
  virtual void add_observer(Observer<T>& Obs) = 0;
  virtual void remove_observer(Observer<T>& Obs) = 0;
  void clear_observers();
protected:
  void notify_update_do(const T& Self);
  void add_observer_do(Observer<T>& Obs, T& t);
  void remove_observer_do(Observer<T>& Obs, T& t);
  void clean_dead_ptr();
};  
}
