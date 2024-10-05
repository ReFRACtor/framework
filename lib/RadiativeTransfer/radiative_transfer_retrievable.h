#ifndef RADIATIVE_TRANSFER_RETRIEVABLE_H
#define RADIATIVE_TRANSFER_RETRIEVABLE_H

#include "state_vector_observer.h"
#include "observer.h"
#include "radiative_transfer.h"

namespace FullPhysics {
/****************************************************************//**
 Interface class for radiative transfer implementations that
 happen to have retrievable parameters.								
*******************************************************************/

class RadiativeTransferRetrievable : public RadiativeTransfer,
				     virtual public StateVectorObserver,
				     public Observable<RadiativeTransferRetrievable> {
public:
  RadiativeTransferRetrievable() {}
  virtual ~RadiativeTransferRetrievable() {}

  virtual void add_observer(Observer<RadiativeTransferRetrievable>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<RadiativeTransferRetrievable>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os, bool UNUSED(Short_form) = false) const 
  { Os << "RadiativeTransferRetrievable";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(RadiativeTransferRetrievable);
FP_EXPORT_OBSERVER_KEY(RadiativeTransferRetrievable);
#endif
