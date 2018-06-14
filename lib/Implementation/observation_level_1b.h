#ifndef OBSERVATION_LEVEL_1B_H
#define OBSERVATION_LEVEL_1B_H

#include "observation.h"
#include "level_1b_sample_coefficient.h"
#include "instrument.h"
#include "forward_model_spectral_grid.h"

namespace FullPhysics {

class ObservationLevel1b : public Observation {
public:
    ObservationLevel1b(const boost::shared_ptr<Level1bSampleCoefficient>& level_1b,
            const boost::shared_ptr<Instrument> &instrument,
            const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grids);

    int num_channels() const;

    const SpectralDomain spectral_domain(int channel_index) const;

    Spectrum radiance(int channel_index, bool skip_jacobian = false) const;

    virtual void print(std::ostream& Os) const
    {
        Os << "ObservationLevel1b";
    }

private:
    boost::shared_ptr<Level1bSampleCoefficient> l1b;
    boost::shared_ptr<Instrument> inst;
    boost::shared_ptr<ForwardModelSpectralGrid> grids;
};
}

#endif
