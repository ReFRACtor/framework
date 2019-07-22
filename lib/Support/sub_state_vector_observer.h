#ifndef SUB_STATE_VECTOR_OBSERVER_H
#define SUB_STATE_VECTOR_OBSERVER_H

#include "state_vector_observer.h"

namespace FullPhysics {

/****************************************************************//**
  A common StateVectorObserver just "owns" a subset of the
  StateVector. This class gives the common behavior for this case.
*******************************************************************/
class SubStateVectorObserver : virtual public StateVectorObserver {
public:
  virtual ~SubStateVectorObserver() {}
  virtual void notify_update(const StateVector& Sv);
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
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
/// Called by mark_used with the subset of the state vector used by
/// this class. The default marks everything as used, but derived
/// classes can override this.
//-----------------------------------------------------------------------

  virtual void mark_used_sub(blitz::Array<bool, 1>& Used) const
  { Used = true; }

//-----------------------------------------------------------------------
/// Called by state_vector_name with the subset of the Sv_name used by
/// this class. The default function doesn't change anything, but
/// derived classes can ovveride this.
//-----------------------------------------------------------------------

  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& UNUSED(Sv_name)) 
    const {}
  virtual void print(std::ostream& Os) const {Os << "SubStateVectorObserver";}
  virtual void notify_add(StateVector& Sv)
  { 
    if(pstart != -1)
      throw Exception("A SubStateVectorObserver can only be attached to one state vector");
    pstart = Sv.observer_claimed_size();
    Sv.observer_claimed_size(pstart + plen);
  }

  virtual void notify_remove(StateVector& UNUSED(Sv))
  {
    pstart = -1;

    // Object is no longer part of the state vector there for the sub state arrays
    // are now empty.
    ArrayAd<double, 1> sv_sub(0, 0);
    blitz::Array<double, 2> cov_sub(0, 0);
    update_sub_state(sv_sub, cov_sub);
  }
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
};
}

FP_EXPORT_KEY(SubStateVectorObserver);

#endif
