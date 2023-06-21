// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"

%{
#include "stacked_radiance_mixin.h"
%}

%import "spectral_domain.i"
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::StackedRadianceMixin);

namespace FullPhysics {

class StackedRadianceMixin : public GenericObject {
public:
  virtual int num_channels() const = 0;
  virtual SpectralDomain spectral_domain(int sensor_index) const = 0;
  boost::optional<blitz::Range> stacked_pixel_range(int sensor_index) const;
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false)
    const = 0;
  virtual Spectrum radiance_all(bool skip_jacobian = false) const;
  std::string print_to_string() const;
  %pickle_serialization();
};
}


