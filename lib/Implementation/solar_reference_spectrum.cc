#include "solar_reference_spectrum.h"

// For to_c_order routine
#include "linear_algebra.h"

#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

SolarReferenceSpectrum::SolarReferenceSpectrum(const boost::shared_ptr<Spectrum>& ref_spectrum_in,
                                               const boost::shared_ptr<SolarDopplerShift>& doppler_shift_in) 
    : doppler_shift_(doppler_shift_in), ref_spec_orig(ref_spectrum_in)
{
    Array<double, 1> x = to_c_order(ref_spec_orig->spectral_domain().data());
    Array<double, 1> y = to_c_order(ref_spec_orig->spectral_range().data());

    if(x.rows() != y.rows()) {
        throw Exception("SpectralDomain and SpectralRange data must have the same size");
    }

    ref_spec_interp = LinearInterpolate<double, double>(x.dataFirst(), x.dataFirst() + x.rows(), y.dataFirst());
}

Spectrum SolarReferenceSpectrum::solar_spectrum(const SpectralDomain& spec_domain) const
{
    Array<double, 1> lookup_grid = spec_domain.convert_wave(ref_spec_orig->spectral_domain().units());
    
    ArrayWithUnit<double, 1> result;
    result.value.resize(lookup_grid.shape());
    result.units = ref_spec_orig->spectral_range().units();

    // Derive doppler shift correction, by default 1.0 if no SolarDopplerShift object available
    double solar_dist = 1.0;
    if (doppler_shift_) {
        solar_dist = doppler_shift_->solar_distance().convert(OldConstant::AU).value;
    }

    // Apply interpolation
    for(int grid_idx = 0; grid_idx < result.value.rows(); grid_idx++) {
        result.value(grid_idx) = ref_spec_interp(lookup_grid(grid_idx)) / (solar_dist * solar_dist);
    }

    return Spectrum(spec_domain, SpectralRange(result.value, result.units));
}

void SolarReferenceSpectrum::print(std::ostream& Os) const
{
    Os << "SolarReferenceSpectrum\n";
}
