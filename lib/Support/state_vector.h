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
class StateVector : public Printable<StateVector>, public Observable<StateVector> {
public:
  StateVector()
    : pstart(0)
  {
  }
  virtual ~StateVector() { }

  virtual void add_observer(Observer<StateVector>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<StateVector>& Obs) 
  { remove_observer_do(Obs, *this);}
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
/// Current covariance of the state vector.
//-----------------------------------------------------------------------

  const blitz::Array<double, 2>& state_covariance() const {return cov_;}

  void update_state(const blitz::Array<double, 1>& X);
  void update_state(const blitz::Array<double, 1>& X, const blitz::Array<double, 2>& Cov);
  void update_state(const ArrayAd<double, 1>& X, const blitz::Array<double, 2>& Cov);
  blitz::Array<bool, 1> used_flag() const;

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
  // Helper value that says what portion of state vector has been
  // claimed by observers. We don't do anything with this value in
  // this class, except make it available to StateVectorObservers when 
  // they are attached.
  int pstart;
};
}

#endif
