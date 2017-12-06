%include "common.i"

%{
#include "instrument_measurement_level_1b.h"
%}

%base_import(instrument_measurement)
%include "level_1b.i"
%include "forward_model_spectral_grid.i"

%fp_shared_ptr(FullPhysics::InstrumentMeasurementLevel1b)

namespace FullPhysics {

class InstrumentMeasurementLevel1b : public InstrumentMeasurement {
public:
    InstrumentMeasurementLevel1b(const boost::shared_ptr<Level1b>& level_1b, const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grids);

    int num_channels() const;

    const SpectralDomain spectral_domain(int channel_index) const;

    Spectrum radiance(int channel_index, bool skip_jacobian = false) const;
};
}
