#ifndef OBSERVER_H
#define OBSERVER_H
#include "null_deleter.h"
#include "printable.h"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/foreach.hpp>
#include <list>
#include <vector>
#include <iostream>
namespace FullPhysics {
template<class T> class Observable;
/****************************************************************//**
  Simple Mixin to be and Observer of another object of class T. We get
  notified when the other object is updated.

  A note on lifetime of objects. When an Observer is destroyed, 
  nothing special needs to be done. The Observable is automatically
  notified that the Object no longer exists, and shouldn't receive 
  notify_update messages anymore. Likewise, nothing special happens
  when the Observable is destroyed, this Observe simple receives no
  more messages from that Observable since it never changes after
  it is dead.

  The relationship between Observer and Observable is m to n, a 
  Observer can be attached to any number of Observables, and likewise
  an Observable can have any number of Observers attached.
*******************************************************************/
template<class T> class Observer : public virtual GenericObject {
public:
  Observer() {
    this_obj = to_ptr(*this);
  }
  virtual ~Observer() {}

//-----------------------------------------------------------------------
/// Called when the Observed object is updated.
//-----------------------------------------------------------------------

  virtual void notify_update(const T& UNUSED(Observed_object)) {};

//-----------------------------------------------------------------------
/// Called when an object is added to an Observable. Default is to 
/// do nothing.
//-----------------------------------------------------------------------

  virtual void notify_add(T& UNUSED(Observed_object)) {}
  virtual void notify_add() {}

//-----------------------------------------------------------------------
/// Called when an object is removed from an Observable. Default is to 
/// do nothing.
//-----------------------------------------------------------------------

  virtual void notify_remove(T& UNUSED(Observed_object)) {}
  virtual void notify_remove() {}

private:
  // Shared ptr that we can use to get weak_ptr to. This object has
  // exactly the same lifetime as the containing object.
  boost::shared_ptr<Observer> this_obj;
  friend class Observable<T>;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  Special case observer where  a cache is invalidated in an Observable
  changes  
*******************************************************************/

class CacheInvalidatedObserver: virtual public GenericObject
{
public:
  CacheInvalidatedObserver()
    : cache_valid(false)
  {
    this_obj_ = to_ptr(*this);
  }
  virtual ~CacheInvalidatedObserver() {}

//-----------------------------------------------------------------------
/// Called when the cache is invalidated.
//-----------------------------------------------------------------------

  virtual void invalidate_cache() { cache_valid = false; }

  boost::shared_ptr<CacheInvalidatedObserver>& this_obj()
  { return this_obj_; }

protected:
  bool cache_valid;
  // Shared ptr that we can use to get weak_ptr to. This object has
  // exactly the same lifetime as the containing object.
  boost::shared_ptr<CacheInvalidatedObserver> this_obj_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};
  
/****************************************************************//**
  Mixin for a class that allows other classes to observe it
  state. When the object state has been updated, the class should
  call notify_update.

  A note on lifetime of objects. When an Observer is destroyed, 
  nothing special needs to be done. The Observable is automatically
  notified that the Object no longer exists, and shouldn't receive 
  notify_update messages anymore. Likewise, nothing special happens
  when the Observable is destroyed, this Observe simple receives no
  more messages from that Observable since it never changes after
  it is dead.

  In some cases, you might actually want the lifetime of the Observer
  to be controlled by the Observable. This would be for a class that
  has no other purpose than to observe another class, say something
  that writes out diagnostic messages whenever a StateVector is
  updated or something like that. In those cases, you can use the
  add_observer_and_keep_reference instead of just add_observer. This
  stashes a copy of the RefPtr so that as long as the Observable
  exists the registered Observer will also.

  The relationship between Observer and Observable is m to n, a 
  Observer can be attached to any number of Observables, and likewise
  an Observable can have any number of Observers attached.

  There are two categories of Observer that we encounter in practice:

  1. Observers that care about what object it is observing, so for
     example a class that does logging for when the state vector
     changes. 
  2. Observers that do some kind of calculation caching, and want to 
     know when the cache is invalid because of a change. These 
     caches often observer multiple objects.
   
  It is useful to separate out the caching observers to:

  1. Make the relationship more explicit
  2. Allow a cache to attach itself to multiple Observables, any
     change will invalidate the cache.

  Another complication is serialization (e.g., calling 
  serializate_write_string). There isn't "one" way to handle
  weak_ptr. The first extreme is to serialize none of them - so
  serialized objects loose all connections. Or we could save all 
  of them, but then serializing a simple object like Pressure might
  pull in the whole world though lots of things pointing to each
  other.

  The way we handle this is to do the serialization twice. The first 
  time, we mark all the Observer<T> that get serialized. The second
  time, we keep all weak_ptr that point to one of these objects. So,
  if you serialize just Pressure you pretty much loose all 
  Observer<Pressure> objects. But if you serialize ForwardModel, you
  get all the Observer<Pressure> that are somehow referenced through
  ForwardModel. 

  This seems like what we will want in general, but we can reevaluate
  this if we change our minds over time.
*******************************************************************/
template<class T> class Observable : public virtual GenericObject {
public:
  virtual ~Observable() {}

//-----------------------------------------------------------------------
/// Add an observer 
//-----------------------------------------------------------------------

  virtual void add_observer(Observer<T>& Obs) = 0;

//-----------------------------------------------------------------------
/// Add an CacheInvalidatedObserver
//-----------------------------------------------------------------------

  void add_cache_invalidated_observer(CacheInvalidatedObserver& Obs)
  {
    CPointerEqual p(Obs.this_obj().get());
    bool already_added =
      std::find_if(colist.begin(), colist.end(), p) != colist.end();
    if (!already_added)
      colist.push_back(boost::weak_ptr<CacheInvalidatedObserver>(Obs.this_obj()));
  }
  void remove_cache_invalidated_observer(CacheInvalidatedObserver& Obs)
  {
    CPointerEqual p(Obs.this_obj().get());
    colist.remove_if(p);
  }
  
//-----------------------------------------------------------------------
/// Add an observer and keep a reference to it. See the discussion in
/// the Observer class description for details.
//-----------------------------------------------------------------------

  void add_observer_and_keep_reference
  (boost::shared_ptr<Observer<T> >& Obs) 
  { ref_list.push_back(Obs); add_observer(*Obs); }

//-----------------------------------------------------------------------
/// Remove an observer 
//-----------------------------------------------------------------------

  virtual void remove_observer(Observer<T>& Obs) = 0;

//-----------------------------------------------------------------------
/// Remove all observers through calls to the remove_observer routine
//-----------------------------------------------------------------------
  void clear_observers()
  {
    // Copy out shared_ptrs into a seperate list since olist will be
    // acted upon and have items removed as we iterate
    std::vector<boost::shared_ptr<Observer<T> > > obs_pointers;
    BOOST_FOREACH(boost::weak_ptr<Observer<T> >& obs_ref, olist) {
      boost::shared_ptr<Observer<T> > obs_ptr = obs_ref.lock();
      if(obs_ptr)
        obs_pointers.push_back(obs_ptr);
    }
    std::vector<boost::shared_ptr<CacheInvalidatedObserver> > cobs_pointers;
    BOOST_FOREACH(boost::weak_ptr<CacheInvalidatedObserver>& cobs_ref, colist) {
      boost::shared_ptr<CacheInvalidatedObserver> cobs_ptr = cobs_ref.lock();
      if(cobs_ptr)
        cobs_pointers.push_back(cobs_ptr);
    }

    // Remove each observer with copied list of pointers, olist will
    // be updated on each call to remove_observer
    BOOST_FOREACH(boost::shared_ptr<Observer<T> >& obs_ptr, obs_pointers) {
      remove_observer(*obs_ptr);
    }
    BOOST_FOREACH(boost::shared_ptr<CacheInvalidatedObserver>& cobs_ptr, cobs_pointers) {
      remove_cache_invalidated_observer(*cobs_ptr);
    }
  }

protected:
//-----------------------------------------------------------------------
/// Function to call to notify Observers of a state change. The object
/// should pass itself to this function, so it can be passed to the
/// Observers. 
//-----------------------------------------------------------------------
  void notify_update_do(const T& Self) const
  {
    const_cast<Observable<T>*>(this)->clean_dead_ptr();
    for(auto t : olist) {
      boost::shared_ptr<Observer<T> > t2 = t.lock();
      if(t2)
        t2->notify_update(Self);
    }
    for(auto t : colist) {
      boost::shared_ptr<CacheInvalidatedObserver> t2 = t.lock();
      if(t2)
        t2->invalidate_cache();
    }
  }

//-----------------------------------------------------------------------
// Helper class for comparing against the weak_ptrs used to keep track
// of observables
//-----------------------------------------------------------------------

  class PointerEqual {
  public:
    PointerEqual(Observer<T>* t) : t_(t) {}
    bool operator()(const boost::weak_ptr<Observer<T> >& w)
    { return w.lock().get() == t_; }
    Observer<T>* t_;
  };
  class CPointerEqual {
  public:
    CPointerEqual(CacheInvalidatedObserver* t) : t_(t) {}
    bool operator()(const boost::weak_ptr<CacheInvalidatedObserver>& w)
    { return w.lock().get() == t_; }
    CacheInvalidatedObserver* t_;
  };

//-----------------------------------------------------------------------
/// Add an observer 
//-----------------------------------------------------------------------

  void add_observer_do(Observer<T>& Obs, T& t)
  {
    PointerEqual p(Obs.this_obj.get());
    bool already_added = std::find_if(olist.begin(), olist.end(), p) != olist.end();
    if (!already_added) {
      olist.push_back(boost::weak_ptr<Observer<T> >(Obs.this_obj));
      Obs.notify_add(t);
      Obs.notify_add();
    }
  }

  void add_observer_do(Observer<T>& Obs)
  {
    PointerEqual p(Obs.this_obj.get());
    bool already_added = std::find_if(olist.begin(), olist.end(), p) != olist.end();
    if (!already_added) {
      olist.push_back(boost::weak_ptr<Observer<T> >(Obs.this_obj));
      Obs.notify_add();
    }
  }

//-----------------------------------------------------------------------
/// Remove an observer 
//-----------------------------------------------------------------------

  void remove_observer_do(Observer<T>& Obs, T& t)
  {
    PointerEqual p(Obs.this_obj.get());
    olist.remove_if(p);
    Obs.notify_remove(t);
    Obs.notify_remove();
  }

  void remove_observer_do(Observer<T>& Obs)
  {
    PointerEqual p(Obs.this_obj.get());
    olist.remove_if(p);
    Obs.notify_remove();
  }

  class PointerDead {
  public:
    bool operator()(const boost::weak_ptr<Observer<T> >& w)
    { return w.expired(); }
    bool operator()(const boost::weak_ptr<CacheInvalidatedObserver>& w)
    { return w.expired(); }
  };
//-----------------------------------------------------------------------
/// Remove any dead pointers.
//-----------------------------------------------------------------------
  void clean_dead_ptr() {
    PointerDead p;
    olist.remove_if(p);
    colist.remove_if(p);
  }
  std::list<boost::weak_ptr<Observer<T> > > olist;
  std::list<boost::weak_ptr<CacheInvalidatedObserver> > colist;
  std::vector<boost::shared_ptr<Observer<T> > > ref_list;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};  
}

#define FP_EXPORT_OBSERVER_KEY(NAME) \
namespace FullPhysics { \
typedef Observer<NAME> Observer ## NAME; \
typedef Observable<NAME> Observable ## NAME; \
} \
FP_EXPORT_KEY(Observer ## NAME); \
FP_EXPORT_KEY(Observable ## NAME);

FP_EXPORT_KEY(CacheInvalidatedObserver);
#endif
