%include "common.i"

%{
#include "stacked_radiance_mixin.h"
%}

%import "spectral_domain.i"
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::StackedRadianceMixin);

namespace FullPhysics {

class StackedRadianceMixin {
public:
    /// Number of spectral channels
    virtual int num_channels() const = 0;

    /// The spectral grid of the radiance values, implemented by inheriting class
    virtual const SpectralDomain spectral_grid(int channel_index) const = 0;

    /// The range of indicies that corresponds to a particular
    /// band in the stacked radiances.
    ///
    /// The range may well be empty if a band is not used at all.
    /// This is a useful edge case, but unfortunately blitz::Range does not
    /// support empty ranges. As a simple work around, we use the
    /// boost::optional class to return a value only if the range is not empty.
    virtual const boost::optional<blitz::Range> stacked_pixel_range(int Spec_index) const;

    /// Per channel radiance data, implemented by inheriting class
    virtual Spectrum radiance(int channel_index) const = 0;

    /// Radiance data all stacked together as one long
    /// spectrum (so band 0, followed by band 1, etc.).
    virtual Spectrum radiance_all() const;
};
}


