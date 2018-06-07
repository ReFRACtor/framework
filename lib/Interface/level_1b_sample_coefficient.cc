#include "level_1b_sample_coefficient.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

/* TODO: Implement return of SpectralDomain from spectral_coefficient */
SpectralDomain Level1bSampleCoefficient::sample_spectral_domain(int Spec_index) const {
    ArrayWithUnit<double, 1> data;
    SpectralDomain _temp = SpectralDomain(data);
    return _temp;
}
