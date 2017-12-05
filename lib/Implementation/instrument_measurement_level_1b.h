#ifndef INSTRUMENT_MEASUREMENT_LEVEL_1B_H
#define INSTRUMENT_MEASUREMENT_LEVEL_1B_H

#include "instrument_measurement.h"
#include "level_1b.h"
#include "forward_model_spectral_grid.h"

namespace FullPhysics {

class InstrumentMeasurementLevel1b : public InstrumentMeasurement {
public:
    InstrumentMeasurementLevel1b(const boost::shared_ptr<Level1b>& level_1b, const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grids);

    int num_channels() const;

    const SpectralDomain spectral_grid(int channel_index) const;

    Spectrum radiance(int channel_index) const;
private:
    boost::shared_ptr<Level1b> l1b;
    boost::shared_ptr<ForwardModelSpectralGrid> grids;
};
}

#endif
