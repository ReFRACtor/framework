#ifndef SUB_STATE_VECTOR_OBSERVER_H
#define SUB_STATE_VECTOR_OBSERVER_H

#include "state_vector_observer.h"

namespace FullPhysics {

/****************************************************************//**
  A common StateVectorObserver just "owns" a subset of the
  StateVector. This class gives the common behavior for this case.

 Additionally defines the interface necessary for interacting with 
 a portion of the state vector that represents an individual 
 retrieval component.

 The intended purpose for some of these methods is is to define the 
 interface expected by the configuration system to use retrieval 
 components to set up which components are added registered into 
 the StateVector as well as obtain the initial guess values from 
 the components.
*******************************************************************/

class SubStateVectorObserver : virtual public StateVectorObserver {
public:
  virtual ~SubStateVectorObserver() {}
  virtual void notify_update(const StateVector& Sv);
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;

//-----------------------------------------------------------------------
/// Starting index of state vector used by this object.
//-----------------------------------------------------------------------

  int state_vector_start_index() const {return pstart; }

//-----------------------------------------------------------------------
/// Length of the sub set of the state vector used by this object.
//-----------------------------------------------------------------------

  int sub_vector_size() const {return plen; }

//-----------------------------------------------------------------------
/// Called by update_state with the subset of the state vector used by
/// this class. 
//-----------------------------------------------------------------------

  virtual void update_sub_state(
    const ArrayAd<double, 1>& Sv_sub,
    const blitz::Array<double, 2>& Cov_sub) = 0;

//-----------------------------------------------------------------------
/// Called by state_vector_name with the subset of the Sv_name used by
/// this class. The default function doesn't change anything, but
/// derived classes can overide this.
//-----------------------------------------------------------------------

  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& UNUSED(Sv_name)) const { ; }

  virtual void notify_add(StateVector& Sv)
  { 
    if(pstart != -1)
      throw Exception("A SubStateVectorObserver can only be attached to one state vector");
    pstart = Sv.observer_claimed_size();
    Sv.observer_claimed_size(pstart + plen);
  }

  virtual void notify_remove(StateVector& Sv)
  {
    pstart = -1;

    // Note that desired semantics of what to do when the state vector
    // goes away isn't clear. A reasonable thing is to zero out the
    // coefficients and update_sub_state to zero size, or we can just
    // leave the existing state alone. i.e., do we consider being
    // detached from a state vector as an update to the state or not?
    //
    // The original behavior was to have the state zeroed out. However
    // for muses-py we want to be able to detach from the state vector
    // of one retrieval step and attach to a new state vector for a
    // new retrieval step. Detaching isn't a state change, we just
    // want to keep our state and detach so we can reattach to a new
    // state vector.
    //
    // We handle this by just directly passing the desired behavior
    // through a "Sv.keep_state_when_removed" variable.

    if(!Sv.keep_state_when_removed()) {
      // Object is no longer part of the state vector there for the sub state arrays
      // are now empty.
      ArrayAd<double, 1> sv_sub(0, 0);
      blitz::Array<double, 2> cov_sub(0, 0);
      update_sub_state(sv_sub, cov_sub);
    }
  }

  //-----------------------------------------------------------------------
  /// Return a string to identify this retrieval component that manages
  /// part of the state, this name should be all lower case and seperate 
  /// parts with a /. For example, an aerosol
  /// named strat would be named as:
  ///   aerosol/strat.
  /// A gas named CO2 would be named like this:
  ///   absorber/co2
  /// The name is intended to be used for looking up retrieval values 
  /// for a configuration system. Classes that have the same type of inputs
  /// should have the same name.
  //-----------------------------------------------------------------------

  virtual std::string sub_state_identifier() const {
      return "unknown/not_set";
  }

  //-----------------------------------------------------------------------
  /// The current portion of values from the state vector managed by this
  /// component. 
  /// Before retrieval starts this would be the initial guess. But they will
  /// change based on updates to the state vector.
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> sub_state_vector_values() const = 0;

protected:
//-----------------------------------------------------------------------
/// Take the given number of state vector parameters. We determine
/// where the starting point to use is when we attach to the state
/// vector. 
///
/// Note that it is perfectly legal for Plen to be 0, that just means
/// we don't have any parameters. This is a useful edge case that we
/// support. 
//-----------------------------------------------------------------------

  SubStateVectorObserver(int Plen)
  : pstart(-1)
  {
    state_vector_observer_initialize(Plen);
  }

//-----------------------------------------------------------------------
/// Default constructor. Derived classes should call initialize before
/// finishing there constructor.
//-----------------------------------------------------------------------

  SubStateVectorObserver() : pstart(-1), plen(0) {}

  void state_vector_observer_initialize(int Plen);

//-----------------------------------------------------------------------
/// The last full state vector we have been updated with, saved for
/// reference by derived class
//-----------------------------------------------------------------------
  ArrayAd<double, 1> sv_full;

//-----------------------------------------------------------------------
/// The last full covariance matrix we have been with, saved for
/// reference by derived class.
//-----------------------------------------------------------------------
  blitz::Array<double, 2> sv_cov_full;

//-----------------------------------------------------------------------
/// The subset of sv_full that is "owned" by this class, what was
/// passed through update_sub_state. Saved for reference by derived
/// class. 
//-----------------------------------------------------------------------
  ArrayAd<double, 1> sv_sub;

//-----------------------------------------------------------------------
/// The subset of cov_full that is "owned" by this class, what was
/// passed through update_sub_state. Saved for reference by derived
/// class. 
//-----------------------------------------------------------------------
  blitz::Array<double, 2> sv_cov_sub;
private:
  int pstart;
  int plen;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(SubStateVectorObserver);

#endif
