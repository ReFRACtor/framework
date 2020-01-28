// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
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
  virtual SpectralDomain spectral_domain(int channel_index) const = 0;
  virtual Spectrum radiance(int channel_index) const = 0;
};
}

%template(Vector_Observation) std::vector<boost::shared_ptr<FullPhysics::Observation> >;
