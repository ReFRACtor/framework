#ifndef ANGLE_UTIL_H
#define ANGLE_UTIL_H

#include <math.h>

#include "unit.h"
#include "double_with_unit.h"

namespace FullPhysics {

    DoubleWithUnit scattering_angle(const DoubleWithUnit& observation_zenith_angle, const DoubleWithUnit& solar_zenith_angle, const DoubleWithUnit& relative_azimuth_angle)
    {
        double vza_rad = observation_zenith_angle.convert(units::rad).value;
        double sza_rad = solar_zenith_angle.convert(units::rad).value;
        double raz_rad = relative_azimuth_angle.convert(units::rad).value;
        
        double temp1 = -1.0 * cos(vza_rad) * cos(sza_rad);
        double temp2 = sqrt(1.0 - cos(vza_rad) * cos(vza_rad));
        double temp3 = sqrt(1.0 - cos(sza_rad) * cos(sza_rad));
        double temp4 = cos(raz_rad);
        double sca_rad = temp1 + temp2 * temp3 * temp4;
        double sca_deg = 180.0 * (1.0 - acos(sca_rad) / units::pi);

        return DoubleWithUnit(sca_deg, units::deg);
    }
}

#endif
