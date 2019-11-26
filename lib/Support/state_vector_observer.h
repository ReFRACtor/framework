#ifndef STATE_VECTOR_OBSERVER
#define STATE_VECTOR_OBSERVER

#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This is an observer of a StateVector. If attached to a StateVector,
  this class gets notified when a state vector is updated.

  It is completely unspecified what an observer does with this
  information, but commonly the class will update its internal state
  based on the state vector update.
*******************************************************************/
class StateVectorObserver : public Printable<StateVectorObserver>,
  virtual public Observer<StateVector> {
public:
  StateVectorObserver() {}
  virtual ~StateVectorObserver() {}

//-----------------------------------------------------------------------
/// Mark elements that we are actively using (i.e., that aren't
/// ignored). You only need to mark the ones that are used as true,
/// everything is already initialized as false. Default is to do nothing.
//-----------------------------------------------------------------------

  virtual void mark_used(const StateVector& UNUSED(Sv), 
			 blitz::Array<bool, 1>& UNUSED(Used)) const {}
			 
//-----------------------------------------------------------------------
/// Update any portion of the list of the state vector names that
/// apply to this object. Default is to do nothing.
//-----------------------------------------------------------------------

  virtual void state_vector_name(const StateVector& UNUSED(Sv), 
		 blitz::Array<std::string, 1>& UNUSED(Sv_name)) const {}
			 
  virtual void print(std::ostream& Os) const { Os << "StateVectorObserver";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateVectorObserver);

#endif
