#ifndef OBSERVER_H
#define OBSERVER_H
#include "null_deleter.h"
#include "generic_object.h"
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
*******************************************************************/
template<class T> class Observable : public virtual GenericObject {
public:
  virtual ~Observable() {}

//-----------------------------------------------------------------------
/// Add an observer 
//-----------------------------------------------------------------------

  virtual void add_observer(Observer<T>& Obs) = 0;

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

    // Remove each observer with copied list of pointers, olist will
    // be updated on each call to remove_observer
    BOOST_FOREACH(boost::shared_ptr<Observer<T> >& obs_ptr, obs_pointers) {
      remove_observer(*obs_ptr);
    }
  }

protected:
//-----------------------------------------------------------------------
/// Function to call to notify Observers of a state change. The object
/// should pass itself to this function, so it can be passed to the
/// Observers. 
//-----------------------------------------------------------------------
  void notify_update_do(const T& Self)
  {
    clean_dead_ptr();
    BOOST_FOREACH(boost::weak_ptr<Observer<T> >& t, olist) {
      boost::shared_ptr<Observer<T> > t2 = t.lock();
      if(t2)
        t2->notify_update(Self);
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
  };
//-----------------------------------------------------------------------
/// Remove any dead pointers.
//-----------------------------------------------------------------------
  void clean_dead_ptr() {
    PointerDead p;
    olist.remove_if(p);
  }
  std::list<boost::weak_ptr<Observer<T> > > olist;
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

#endif

