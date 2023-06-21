#ifndef OBSERVATION_SV_H
#define OBSERVATION_SV_H
#include "state_vector_observer.h"
#include "observation.h"
#include "observer.h"
#include "spectrum.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains a Observation that changes based on the updates
  to the StateVector. This is similiar to the L1B data we used with
  OCO-2, but Muses-py uses a bit of a different formulation where some
  of the corrections we are doing happens to the input radiance rather
  than at the forward model side.  This actually makes sense, so
  things like correction for calibration occurs with the radiance
  input rather than altering the forward model to match the incorrect
  calibration. So basically this is like L1B but with possible
  corrections applied by the StateVector.

  Note that muses sometimes scales the data by the solar model, so the
  solar model is included on the radiance field side rather than the
  ForwardModel side. This means that the Spectrum is actually
  reflectance rather than radiance. But "radiance_or_reflectance" is
  kind of a long name, so we just have this one function to handle
  both cases. If unsure, you can just check the units of the
  SpectralRange in the returned Spectrum to see if the units are
  radiance or dimensionless like a reflectance.
  
  Other objects may depend on the ObservationSv, and should be updated
  when the ObservationSv is updated. To facilitate that, this class is
  an Oberverable, and objects can add themselves as Observers to be
  notified when the ObservationSv is updated.
 *******************************************************************/
class ObservationSv : virtual public Observation,
		      virtual public StateVectorObserver, public Observable<ObservationSv> {
public:
  virtual ~ObservationSv() {}
  virtual void add_observer(Observer<ObservationSv>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<ObservationSv>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Clone a ObservationSv object. Note that the cloned version
/// will *not* be attached to a StateVector, although you can of
/// course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" ObservationSv object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<ObservationSv> clone() const = 0;
  virtual void print(std::ostream& Os) const
  { Os << "ObservationSv";}
protected:
  ObservationSv() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(ObservationSv);
FP_EXPORT_OBSERVER_KEY(ObservationSv);
#endif
