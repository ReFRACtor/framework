%include "fp_common.i"

%{
#include "observation_level_1b.h"
%}

%base_import(observation)
%import "level_1b.i"
%import "forward_model_spectral_grid.i"

%fp_shared_ptr(FullPhysics::ObservationLevel1b)

namespace FullPhysics {

class ObservationLevel1b : public Observation {
public:
    ObservationLevel1b(const boost::shared_ptr<Level1b>& level_1b, 
            const boost::shared_ptr<Instrument> &instrument,
            const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grids);

    int num_channels() const;

    SpectralDomain spectral_domain(int sensor_index) const;

    Spectrum radiance(int sensor_index, bool skip_jacobian = false) const;

    %pickle_serialization();
};
}
