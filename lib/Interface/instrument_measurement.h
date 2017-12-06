#ifndef INSTRUMENT_MEASUREMENT_H
#define INSTRUMENT_MEASUREMENT_H

#include "printable.h"
#include "stacked_radiance_mixin.h"

namespace FullPhysics {

class InstrumentMeasurement : public StackedRadianceMixin, public Printable<InstrumentMeasurement> {
public:
    /// Number of spectral channels
    virtual int num_channels() const = 0;

    /// The spectral grid of the radiance values, implemented by inheriting class
    virtual const SpectralDomain spectral_domain(int channel_index) const = 0;

    /// Measured spectrum for the given spectral channel. Note that this may be empty.
    ///
    /// \param channel_index Band to give value for
    /// \return The set of radiances, possibly empty.
    virtual Spectrum radiance(int channel_index, bool skip_jacobian = false) const = 0;

    virtual void print(std::ostream& Os) const
    {
        Os << "InstrumentMeasurement";
    }
};
}

#endif
