#ifndef STATE_VECTOR_H
#define STATE_VECTOR_H
#include "printable.h"
#include "observer.h"
#include "array_ad.h"
#include "fp_exception.h"
#include <blitz/array.h>

namespace FullPhysics {

/****************************************************************//**
  This handles informing a set of interested objects when the state
  vector has updated. Those objects then update their internal state
  to account for the new state vector.
*******************************************************************/
class StateVector : public Printable<StateVector>,
		    public Observable<StateVector> {
public:
  StateVector()
    : keep_state_when_removed_(false),
      pstart(0)
  {
  }
  virtual ~StateVector() { }

  virtual void add_observer(Observer<StateVector>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<StateVector>& Obs) override {
    remove_observer_do(Obs, *this);
  }
  virtual void remove_observer(Observer<StateVector>& Obs, bool Keep_state_when_removed) {
    keep_state_when_removed_ = Keep_state_when_removed;
    remove_observer_do(Obs, *this);
  }

  //-----------------------------------------------------------------------
  /// Note that desired semantics of what to do when the state vector
  /// detaches from an Observer isn't clear. A reasonable thing is to
  /// zero out the coefficients and update_sub_state to zero size, or
  /// we can just leave the existing state alone. i.e., do we consider
  /// being detached from a state vector as an update to the state
  /// or not?
  ///
  /// The original behavior was to have the state zeroed out. However
  /// for muses-py we want to be able to detach from the state vector
  /// of one retrieval step and attach to a new state vector for a
  /// new retrieval step. Detaching isn't a state change, we just
  /// want to keep our state and detach so we can reattach to a new
  /// state vector.
  ///
  /// We handle this by just directly passing the desired behavior
  /// through the Keep_state_when_removed variable. If false, we have
  /// the old default semantics of clearing the state to a zero sized
  /// array, If true, we detach from the Observer but the Observer
  /// doesn't change its state.
  //-----------------------------------------------------------------------
  
  void clear_observers(bool Keep_state_when_removed=false)
  {
    keep_state_when_removed_ = Keep_state_when_removed;
    Observable<StateVector>::clear_observers();
  }
  inline bool keep_state_when_removed() const {return keep_state_when_removed_;}
  virtual void print(std::ostream& Os) const;

//-----------------------------------------------------------------------
/// Current state vector.
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& state() const { return x_.value(); }

//-----------------------------------------------------------------------
/// Return the state vector as state() does, but also make each value
/// a AutoDerivative. The derivative is with respect to the state
/// vector, i.e., we treat the state vector as the independent
/// variables. This means the first value has a gradient all 0's
/// except for 1 in the first index, the second value all zeros except
/// for 1 in the second index, etc.
//-----------------------------------------------------------------------

  const ArrayAd<double, 1>& state_with_derivative() const { return x_;}

  blitz::Array<std::string, 1> state_vector_name() const;

//-----------------------------------------------------------------------
/// Get the default state vector names that are used for unobserved
/// portions of the state vector.
//-----------------------------------------------------------------------

  std::vector<std::string> default_state_vector_name() const { return default_sv_names; };

//-----------------------------------------------------------------------
/// Set the default state vector names that are used for unobserved
/// portions of the state vector. If not set then a generic name
/// will be used for those portions of the state vector.
//-----------------------------------------------------------------------

  void default_state_vector_name(const std::vector<std::string>& default_name) {
    default_sv_names = default_name;
  }

//-----------------------------------------------------------------------
/// Current covariance of the state vector.
//-----------------------------------------------------------------------

  const blitz::Array<double, 2>& state_covariance() const {return cov_;}

  void update_state(const blitz::Array<double, 1>& X);
  void update_state(const blitz::Array<double, 1>& X,
		    const blitz::Array<double, 2>& Cov);
  void update_state(const ArrayAd<double, 1>& X,
		    const blitz::Array<double, 2>& Cov);

//-----------------------------------------------------------------------
/// Total "claimed" size of the state vector. For observers that
/// register an interest in a portion of the state vector, we add all
/// of the portions of interest. Note that an actual state vector
/// isn't constrained to this size, it might be larger (with
/// presumably portions ignored), or if the observers handle it
/// correctly it could be smaller.
//-----------------------------------------------------------------------

  int observer_claimed_size() const {return pstart;}

//-----------------------------------------------------------------------
/// Update claimed size of state vector.
//-----------------------------------------------------------------------

  void observer_claimed_size(int Pstart) { pstart = Pstart; }

private:
  ArrayAd<double, 1> x_;
  blitz::Array<double, 2> cov_;

  // Default state vector names for elements not connected to observers
  std::vector<std::string> default_sv_names;

  // Pass argument to observers when we remove what the semantics
  // should be - do we want to keep the state or update it?
  mutable bool keep_state_when_removed_;

  // Helper value that says what portion of state vector has been
  // claimed by observers. We don't do anything with this value in
  // this class, except make it available to StateVectorObservers when 
  // they are attached.
  int pstart;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateVector);
FP_EXPORT_OBSERVER_KEY(StateVector);
#endif
