#ifndef SURFACE_TEMPERATURE_DIRECT_H
#define SURFACE_TEMPERATURE_DIRECT_H

#include "surface_temperature.h"
#include "sub_state_vector_array.h"
#include "double_with_unit.h"
#include "auto_derivative_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
 Implements a direct representation of the surface temperature
 in the atmospheric state from the supplied value.
*******************************************************************/
class SurfaceTemperatureDirect : public SubStateVectorArray<SurfaceTemperature> {
public:
    SurfaceTemperatureDirect(DoubleWithUnit surf_temp);
    virtual ~SurfaceTemperatureDirect() {}

    //-----------------------------------------------------------------------
    /// Return the temperature of the surface. This is different than the
    /// temperature near the surface which would be the lowest level of
    /// the temperature grid.
    //-----------------------------------------------------------------------

    virtual AutoDerivativeWithUnit<double> surface_temperature() const;
    virtual boost::shared_ptr<SurfaceTemperature> clone() const;
private:
    Unit units;
};
}
#endif
