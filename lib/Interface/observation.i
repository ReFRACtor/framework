%include "fp_common.i"

%{
#include "observation.h"
%}

%base_import(stacked_radiance_mixin)

%fp_shared_ptr(FullPhysics::Observation)

namespace FullPhysics {

class Observation : public StackedRadianceMixin {
public:
    virtual int num_channels() const = 0;
    virtual const SpectralDomain spectral_domain(int channel_index) const = 0;
    virtual Spectrum radiance(int channel_index) const = 0;
};
}
