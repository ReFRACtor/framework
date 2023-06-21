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
  %python_attribute_abstract(num_channels, int)
  virtual SpectralDomain spectral_domain(int sensor_index) const = 0;
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false) const = 0;
  %pickle_serialization();
};
}

%template(Vector_Observation) std::vector<boost::shared_ptr<FullPhysics::Observation> >;
