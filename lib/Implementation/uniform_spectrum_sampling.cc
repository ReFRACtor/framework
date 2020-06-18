#include "uniform_spectrum_sampling.h"

using namespace FullPhysics;
using namespace blitz;

// See base class for description.

SpectralDomain UniformSpectrumSampling::spectral_domain(int spec_index,
                                                        const SpectralDomain& Lowres_grid, 
                                                        const DoubleWithUnit& Edge_extension) const
{
    range_check(spec_index, 0, spec_spacing.rows());

    DoubleWithUnit spacing = spec_spacing(spec_index);

    // Units we input and want to output and the units we use for computing the grid
    Unit u_out  = Lowres_grid.units();
    Unit u_grid = spacing.units;

    // Convert low res grid values to same units as spacing, in sorted order
    std::vector<double> lowres_conv;
    BOOST_FOREACH(double v, Lowres_grid.data()) {
        lowres_conv.push_back(DoubleWithUnit(v, Lowres_grid.units()).convert_wave(u_grid).value);
    }
    std::sort(lowres_conv.begin(), lowres_conv.end());

    // Determine an integral amount of spacing units to add to the beginning and end of the grid
    DoubleWithUnit offset = round(Edge_extension.convert_wave(u_grid) / spacing) * spacing; 

    // Determine the bounds of the high resolution grid, account for fact the
    // low res grid points might be in a different order before conversion
    double highres_beg = lowres_conv[0] - offset.value;
    double highres_end = lowres_conv[lowres_conv.size()-1] + offset.value;

    int nsamples = (int) ceil((highres_end - highres_beg) / spacing.value) + 1;

    double extension_val = Edge_extension.convert_wave(u_grid).value;

    std::vector<double> highres_points;
    for(int i = 0; i < nsamples; ++i) {
        DoubleWithUnit hr_point(highres_beg + i * spacing.value, u_grid);
        highres_points.push_back(hr_point.convert_wave(u_out).value);
    }

    // Make sure the computed grid points are
    // in sorted order
    std::sort(highres_points.begin(), highres_points.end());

    // Create result object
    Array<double, 1> uniform_vals((int) highres_points.size());
    std::copy(highres_points.begin(), highres_points.end(), uniform_vals.begin());
    SpectralDomain uniform_grid = SpectralDomain(uniform_vals, u_out);

    return uniform_grid;
}

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

void UniformSpectrumSampling::print(std::ostream& Os) const 
{ 
  Os << "UniformSpectrumSampling\n"
     << "  Spacing: " << spec_spacing << "\n";
}
