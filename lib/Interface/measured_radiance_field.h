#ifndef MEASURED_RADIANCE_FIELD_H
#define MEASURED_RADIANCE_FIELD_H
#include "state_vector_observer.h"
#include "observer.h"
#include "spectrum.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the measured radiance field. This is similiar
  to the L1B data we used with OCO-2, but Muses-py uses a bit of a
  different formulation where some of the corrections we are doing
  happens to the input radiance rather than at the forward model side.
  This actually makes sense, so things like correction for calibration
  occurs with the radiance input rather than altering the forward
  model to match the incorrect calibration. So basically this is like
  L1B but with possible corrections applied by the StateVector.

  Note that muses sometimes scales the data by the solar model, so the
  solar model is included on the radiance field side rather than the
  ForwardModel side. This means that the Spectrum is actually
  reflectance rather than radiance. But
  "MeasuredRadianceOrReflectanceField" is kind of a long name, so we
  just have this one class to handle both cases. If unsure, you can
  just check the units of the SpectralRange in the returned Spectrum
  to see if the units are radiance or dimensionless like a
  reflectance.
  
  Other objects may depend on the MeasuredRadianceField, and should be updated
  when the MeasuredRadianceField is updated. To facilitate that, this class is
  an Oberverable, and objects can add themselves as Observers to be
  notified when the MeasuredRadianceField is updated.
 *******************************************************************/
class MeasuredRadianceField : public Printable<MeasuredRadianceField>,
	      virtual public StateVectorObserver, public Observable<MeasuredRadianceField> {
public:
  virtual ~MeasuredRadianceField() {}
  virtual void add_observer(Observer<MeasuredRadianceField>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<MeasuredRadianceField>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Measured radiance for the given spectral index. Note in general
/// we have a jacobian and uncertainty (e.g., nesr) in the generated
/// Spectrum.
///
///  Note that muses sometimes scales the data by the solar model, so
///  the solar model is included on the radiance field side rather
///  than the ForwardModel side. This means that the Spectrum is
///  actually reflectance rather than radiance. But
///  "MeasuredRadianceOrReflectanceField" is kind of a long name, , so
///  we just have this one class to handle both cases. If unsure, you
///  can just check the units of the SpectralRange in the returned
///  Spectrum to see if the units are radiance or dimensionless like a
///  reflectance.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Spectrum> measured_radiance_field(int spec_index) const = 0;

//-----------------------------------------------------------------------
/// Clone a MeasuredRadianceField object. Note that the cloned version
/// will *not* be attached to a StateVector, although you can of
/// course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" MeasuredRadianceField object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<MeasuredRadianceField> clone() const = 0;
  virtual void print(std::ostream& Os) const
  { Os << "MeasuredRadianceField";}
protected:
  MeasuredRadianceField() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MeasuredRadianceField);
FP_EXPORT_OBSERVER_KEY(MeasuredRadianceField);
#endif
