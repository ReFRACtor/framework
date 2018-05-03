#ifndef SURFACE_TEMPERATURE_H
#define SURFACE_TEMPERATURE_H
#include "state_vector_observer.h"
#include "observer.h"
#include "auto_derivative_with_unit.h"
#include "array_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
  This class represents the surface of the temperature under 
  observation. This is different than the temperature of the 
  atmosphere near the surface.

  The surface temperature is normally used for black body 
  calculations in thermal infrared bands. Therefore it is 
  unnecessay when doing near infrared or ultraviolet computations.
*******************************************************************/
class SurfaceTemperature : virtual public StateVectorObserver, public Observable<SurfaceTemperature> {
public:
    virtual ~SurfaceTemperature() {}

    virtual void add_observer(Observer<SurfaceTemperature>& Obs) 
    { add_observer_do(Obs, *this);}

    virtual void remove_observer(Observer<SurfaceTemperature>& Obs) 
    { remove_observer_do(Obs, *this);}

    //-----------------------------------------------------------------------
    /// Return the temperature of the surface. This is different than the
    /// temperature near the surface which would be the lowest level of
    /// the temperature grid.
    //-----------------------------------------------------------------------

    virtual AutoDerivativeWithUnit<double> surface_temperature(int channel_index) const = 0;

    //-----------------------------------------------------------------------
    /// Clone a SurfaceTemperature object. Note that the cloned version will *not*
    /// be attached to a StateVector or Observer<SurfaceTemperature>, although you
    /// can of course attach them after receiving the cloned object.
    ///
    /// Because this isn't attached to the StateVector, one use of the
    /// clone operator is to create a "frozen" SurfaceTemperature object.
    //-----------------------------------------------------------------------

    virtual boost::shared_ptr<SurfaceTemperature> clone() const = 0;
};
}
#endif
