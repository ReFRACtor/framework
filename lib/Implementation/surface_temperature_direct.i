
%include "common.i"

%{
#include "surface_temperature_direct.h"
%}

%base_import(surface_temperature)
%import "sub_state_vector_array.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::SurfaceTemperatureDirect);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::SurfaceTemperature>)

%template(SubStateVectorArraySurfaceTemperatureDirect) FullPhysics::SubStateVectorArray<FullPhysics::SurfaceTemperature>;

namespace FullPhysics {
class SurfaceTemperatureDirect : public SubStateVectorArray<SurfaceTemperature> {
public:
    SurfaceTemperatureDirect(const ArrayWithUnit<double, 1>& surf_temp, blitz::Array<bool, 1> flag);
    virtual ~SurfaceTemperatureDirect() {}

    //-----------------------------------------------------------------------
    /// Return the temperature of the surface. This is different than the
    /// temperature near the surface which would be the lowest level of
    /// the temperature grid.
    //-----------------------------------------------------------------------

    virtual AutoDerivativeWithUnit<double> surface_temperature(int channel_index) const;
    virtual boost::shared_ptr<SurfaceTemperature> clone() const;

    std::string state_vector_name_i(int i) const;

private:
    Unit units;
};
}
