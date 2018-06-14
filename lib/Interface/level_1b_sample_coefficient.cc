#include "level_1b_sample_coefficient.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

/* TODO: Implement return of SpectralDomain from spectral_coefficient */
SpectralDomain Level1bSampleCoefficient::sample_spectral_domain(int Spec_index) const {
    // Get number samples
    SpectralRange spec_radiance = this->radiance(Spec_index);
    Array<double, 1> rad_data = spec_radiance.data();
    int number_samples = rad_data.extent(0);
    cout << "Hi! Got number_samples" << number_samples << endl;
    // Allocate ArraywithUnit with that number of samples
    // For each sample, calculate wavelength from polynomial
    // Create spectralDomain from data
    ArrayWithUnit<double, 1> data;
    SpectralDomain _temp = SpectralDomain(data);
    return _temp;
}
