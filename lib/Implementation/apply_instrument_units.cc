#include "apply_instrument_units.h"

using namespace FullPhysics;

boost::shared_ptr<InstrumentCorrection> ApplyInstrumentUnits::clone() const
{
    return boost::shared_ptr<InstrumentCorrection>(new ApplyInstrumentUnits(spectral_units, scaling));
}

void ApplyInstrumentUnits::apply_correction(const SpectralDomain& Pixel_grid, const std::vector<int>& Pixel_list, SpectralRange& Radiance) const
{
    // Use the copy operator to modify the units
    Radiance = SpectralRange(Radiance.data_ad(), spectral_units, Radiance.uncertainty());

    // Apply scale factor whether or not there are jacobians
    if(Radiance.data_ad().number_variable() == 0) {
        Radiance.data() = Radiance.data() * scaling;
    } else {
        for(int i = 0; i < Radiance.data_ad().rows(); ++i) {
            Radiance.data_ad()(i) = Radiance.data_ad()(i) * scaling;
        }
    }
}

void ApplyInstrumentUnits::print(std::ostream& Os) const {
    Os << "ApplyInstrumentUnits"
       << "\tunits: " << spectral_units << std::endl
       << "\tscale_factor: " << scaling << std::endl;
}
