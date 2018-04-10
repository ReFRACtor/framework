#include "uniform_spectrum_sampling.h"

using namespace FullPhysics;
using namespace blitz;

// See base class for description.

SpectralDomain UniformSpectrumSampling::spectral_domain(int spec_index,
                                                        const SpectralDomain& Lowres_grid, 
                                                        const DoubleWithUnit& Ils_half_width) const
{
    range_check(spec_index, 0, spec_spacing.rows());

    DoubleWithUnit spacing = spec_spacing(spec_index);

    // Determine an integral amount of spacing units to add to the beginning and end of the grid
    DoubleWithUnit offset = round(Ils_half_width.convert_wave(spacing.units) / spacing) * spacing; 
    
    // Convert low res grid values to same units as spacing
    Array<double, 1> lowres_spac_units = Lowres_grid.convert_wave(spacing.units);

    // Determine the bounds of the high resolution grid, account for fact the
    // low res grid points might be in a different order before conversion
    double highres_beg = min(lowres_spac_units) - offset.value;
    double highres_end = max(lowres_spac_units) + offset.value;

    int nsamples = (int) floor((highres_end - highres_beg) / spacing.value) + 1;
    Array<double, 1> highres_grid(nsamples);
    for(int i = 0; i < nsamples; ++i) {
        highres_grid(i) = highres_beg + i * spacing.value;
    }

    return SpectralDomain(highres_grid, spacing.units);
}

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

void UniformSpectrumSampling::print(std::ostream& Os) const 
{ 
  Os << "UniformSpectrumSampling\n"
     << "  Spacing: " << spec_spacing << "\n";
}
