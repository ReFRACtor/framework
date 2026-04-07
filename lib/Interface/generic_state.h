#ifndef GENERIC_STATE_H
#define GENERIC_STATE_H
#include "state_vector_observer.h"
#include "sub_state_vector_array.h"
#include "observer.h"

namespace FullPhysics {
/****************************************************************//**
  Most of our state element follow the same pattern - some basic
  set of coefficients that possibly map to a StateVector, and some
  basic summary of things calculated from the coefficients (possibly
  with a StateMapping).

  For C++, we just derive from a template class. However, we can't
  do that with python. So we have introduced a "GenericState" and
  "GenericStateImpBase" that we can use in python code. Note you should
  *not* use this for a C++ class - instead just derive from
  SubStateVectorArray. But this can be useful in python for
  1) Initial testing before we create a new C++ class or
  2) Just to use indefinitely if there is never a need to pass this
  to C++ (e.g., a python only class uses the GenericState.
*******************************************************************/
class GenericState : public Printable<GenericState>,
		     virtual public StateVectorObserver,
		     public Observable<GenericState> {
public:
  virtual ~GenericState() {}

  virtual void add_observer(Observer<GenericState>& Obs) 
  { add_observer_do(Obs, *this);}

  virtual void remove_observer(Observer<GenericState>& Obs) 
  { remove_observer_do(Obs, *this);}

  //-----------------------------------------------------------------------
  /// Clone a Pcloud object. Note that the cloned version will *not*
  /// be attached to a StateVector or Observer<Pcloud>, although you
  /// can of course attach them after receiving the cloned object.
  ///
  /// Because this isn't attached to the StateVector, one use of the
  /// clone operator is to create a "frozen" Pcloud object.
  //-----------------------------------------------------------------------

  virtual boost::shared_ptr<GenericState> clone() const = 0;
  virtual void print(std::ostream& Os) const
  { Os << "GenericState";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

class GenericStateImpBase :
    virtual public SubStateVectorArray<GenericState> {
public:
  GenericStateImpBase() {}
  virtual ~GenericStateImpBase() {}

  virtual boost::shared_ptr<GenericState> clone() const
  { return boost::make_shared<GenericStateImpBase>(); }
  virtual std::string sub_state_identifier() const
  { return "generic_state"; }
  std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const
  { Os << "GenericStateImpBase";}

private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<GenericState> SubStateVectorArrayGenericState;
}
  
FP_EXPORT_KEY(GenericState);
FP_EXPORT_OBSERVER_KEY(GenericState);
FP_EXPORT_KEY(GenericStateImpBase)
FP_EXPORT_KEY(SubStateVectorArrayGenericState);

#endif
