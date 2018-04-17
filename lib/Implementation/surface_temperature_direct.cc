#include "surface_temperature_direct.h"

using namespace FullPhysics;

SurfaceTemperatureDirect::SurfaceTemperatureDirect(DoubleWithUnit surf_temp)
{
    coeff.value()(0) = surf_temp.value;

    // Save units seperately so it can not be saved into the state vector coefficients
    units = surf_temp.units;
}

AutoDerivativeWithUnit<double> SurfaceTemperatureDirect::surface_temperature() const
{
    return AutoDerivativeWithUnit<double>(coeff(0), units);
}

boost::shared_ptr<SurfaceTemperature> SurfaceTemperatureDirect::clone() const
{
    return boost::shared_ptr<SurfaceTemperatureDirect>(new SurfaceTemperatureDirect(DoubleWithUnit(coeff(0).value(), units)));
}
