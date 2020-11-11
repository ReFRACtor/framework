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
class StateVectorObserver : public Observer<StateVector> {
public:
  StateVectorObserver() {}
  virtual ~StateVectorObserver() {}

//-----------------------------------------------------------------------
/// Update any portion of the list of the state vector names that
/// apply to this object. Default is to do nothing.
//-----------------------------------------------------------------------

  virtual void state_vector_name(const StateVector& UNUSED(Sv), 
		 blitz::Array<std::string, 1>& Sv_name) const
  {
    if(sv_name_.size() == 0)
      return;
    if((int) sv_name_.size() != Sv_name.rows())
      throw Exception("sv_name_ and Sv_name need to be the same size");
    for(int i = 0; i < Sv_name.rows(); ++i)
      Sv_name(i) = sv_name_[i];
  }
			 
//-----------------------------------------------------------------------
/// Python doesn't work so well with the state_vector_name
/// interface, so give access to places to set this.
//-----------------------------------------------------------------------
  
  const std::vector<std::string>& sv_name() const { return sv_name_; }
  void sv_name(const std::vector<std::string>& V) 
  { sv_name_ = V; }
private:
  std::vector<std::string> sv_name_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateVectorObserver);

#endif
