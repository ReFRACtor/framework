#include "surface_temperature_direct.h"

using namespace FullPhysics;

SurfaceTemperatureDirect::SurfaceTemperatureDirect(const ArrayWithUnit<double, 1>& surf_temp, blitz::Array<bool, 1> flag)
{
    if(surf_temp.rows() != flag.rows()) {
        Exception error;
        error << "Number of surface temperature channels: " << surf_temp.rows() 
              << " must match retrieval flag array size: " << flag.rows();
        throw error;
    }

    // Save units seperately so it can not be saved into the state vector coefficients
    units = surf_temp.units;

    init(surf_temp.value, flag);
}

AutoDerivativeWithUnit<double> SurfaceTemperatureDirect::surface_temperature(int channel_index) const
{
    return AutoDerivativeWithUnit<double>(coeff(channel_index), units);
}

boost::shared_ptr<SurfaceTemperature> SurfaceTemperatureDirect::clone() const
{
    return boost::shared_ptr<SurfaceTemperatureDirect>(new SurfaceTemperatureDirect(ArrayWithUnit<double, 1>(coefficient().value(), units), used_flag_value()));
}

std::string SurfaceTemperatureDirect::state_vector_name_i(int i) const
{
    std::stringstream sv_name;
    sv_name << "Surface Temperature for channel " << i;
    return sv_name.str();
}
