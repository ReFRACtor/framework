#ifndef PCLOUD_H
#define PCLOUD_H
#include "state_vector_observer.h"
#include "observer.h"
#include "auto_derivative_with_unit.h"
#include "array_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
  This class represents the pressure of the cloud center under observation. 
*******************************************************************/
class Pcloud : public Printable<Pcloud>,
	       virtual public StateVectorObserver,
	       public Observable<Pcloud> {
public:
  virtual ~Pcloud() {}

  virtual void add_observer(Observer<Pcloud>& Obs) 
  { add_observer_do(Obs, *this);}

  virtual void remove_observer(Observer<Pcloud>& Obs) 
  { remove_observer_do(Obs, *this);}

  //-----------------------------------------------------------------------
  /// Return pressure of the top of the cloud
  //-----------------------------------------------------------------------

  virtual AutoDerivativeWithUnit<double> pressure_cloud(int sensor_index) const = 0;

  //-----------------------------------------------------------------------
  /// Clone a Pcloud object. Note that the cloned version will *not*
  /// be attached to a StateVector or Observer<Pcloud>, although you
  /// can of course attach them after receiving the cloned object.
  ///
  /// Because this isn't attached to the StateVector, one use of the
  /// clone operator is to create a "frozen" Pcloud object.
  //-----------------------------------------------------------------------

  virtual boost::shared_ptr<Pcloud> clone() const = 0;
  virtual void print(std::ostream& Os) const
  { Os << "Pcloud";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}
FP_EXPORT_KEY(Pcloud);
FP_EXPORT_OBSERVER_KEY(Pcloud);
#endif
