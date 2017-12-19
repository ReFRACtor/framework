#ifndef OBSERVATION_H
#define OBSERVATION_H

#include "printable.h"
#include "stacked_radiance_mixin.h"

namespace FullPhysics {

class Observation : public StackedRadianceMixin, public Printable<Observation> {
public:
    /// Number of spectral channels
    virtual int num_channels() const = 0;

    /// The spectral grid of the radiance values, implemented by inheriting class
    virtual const SpectralDomain spectral_domain(int channel_index) const = 0;

    /// Measured spectrum for the given spectral channel. Note that this may be empty.
    ///
    /// \param channel_index Band to give value for
    /// \param skip_jacobian If true, don't do the Jacobian
    /// \return The set of radiances, possibly empty.
    virtual Spectrum radiance(int channel_index, bool skip_jacobian = false) const = 0;

    virtual void print(std::ostream& Os) const
    {
        Os << "Observation";
    }
};
}

#endif
