#ifndef SCALE_CLOUD_H
#define SCALE_CLOUD_H
#include "state_vector_observer.h"
#include "observer.h"
#include "auto_derivative.h"

namespace FullPhysics {
/****************************************************************//**
  This class represents the scale of cloud log thickness under 
  observation. 
*******************************************************************/
class ScaleCloud : public Printable<ScaleCloud>,
			   virtual public StateVectorObserver,
			   public Observable<ScaleCloud> {
public:
  virtual ~ScaleCloud() {}

  virtual void add_observer(Observer<ScaleCloud>& Obs) 
  { add_observer_do(Obs, *this);}

  virtual void remove_observer(Observer<ScaleCloud>& Obs) 
  { remove_observer_do(Obs, *this);}

  //-----------------------------------------------------------------------
  /// Return the scale to cloud.
  //-----------------------------------------------------------------------

  virtual AutoDerivative<double> scale_cloud(int sensor_index) const = 0;

  //-----------------------------------------------------------------------
  /// Clone a ScaleCloud object. Note that the cloned version will *not*
  /// be attached to a StateVector or Observer<ScaleCloud>, although you
  /// can of course attach them after receiving the cloned object.
  ///
  /// Because this isn't attached to the StateVector, one use of the
  /// clone operator is to create a "frozen" ScaleCloud object.
  //-----------------------------------------------------------------------

  virtual boost::shared_ptr<ScaleCloud> clone() const = 0;
  virtual void print(std::ostream& Os) const
  { Os << "ScaleCloud";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}
FP_EXPORT_KEY(ScaleCloud);
FP_EXPORT_OBSERVER_KEY(ScaleCloud);
#endif
