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
			 blitz::Array<bool, 1>& Used) const
  {
    if(used_flag_.rows() == 0)
      return;
    if(used_flag_.rows() != Used.rows())
      throw Exception("used_flag_ and Used need to be the same size");
    Used = used_flag_;
  }
			 
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
			 
  virtual void print(std::ostream& Os) const { Os << "StateVectorObserver";}

//-----------------------------------------------------------------------
/// Python doesn't work so well with the mark_used and state_vector_name
/// interface, so give access to places to set this.
//-----------------------------------------------------------------------
  
  const blitz::Array<bool, 1>& used_flag() const { return used_flag_;}
  void used_flag(const blitz::Array<bool, 1>& V)
  {used_flag_.reference(V.copy());}
  const std::vector<std::string>& sv_name() const { return sv_name_; }
  void sv_name(const std::vector<std::string>& V) 
  { sv_name_ = V; }
private:
  blitz::Array<bool, 1> used_flag_;
  std::vector<std::string> sv_name_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateVectorObserver);

#endif
