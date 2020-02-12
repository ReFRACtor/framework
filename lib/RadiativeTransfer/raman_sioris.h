#ifndef RAMAN_SIORIS_H
#define RAMAN_SIORIS_H

#include <boost/shared_ptr.hpp>
#include <blitz/array.h>

#include "spectral_domain.h"

namespace FullPhysics {

blitz::Array<double, 1> compute_raman_sioris(double solar_zenith, double viewing_zenith, double scattering_angle, double albedo, bool do_upwelling, const blitz::Array<double, 1> &temperature_layers, const blitz::Array<double, 1>& air_number_density, const SpectralDomain &grid, const blitz::Array<double, 1> &solar_irradiance, const blitz::Array<double, 2> &total_optical_depth);

}

#endif
