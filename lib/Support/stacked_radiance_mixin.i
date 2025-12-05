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

// Allow these classes to be derived from in Python.
%feature("director") StackedRadianceMixin;
  
class StackedRadianceMixin : public GenericObject {
public:
  %python_attribute_abstract(num_channels, int);
  virtual SpectralDomain spectral_domain(int sensor_index) const = 0;
  boost::optional<blitz::Range> stacked_pixel_range(int sensor_index) const;
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false)
    const = 0;
  virtual SpectralDomain spectral_domain_all() const;
  virtual Spectrum radiance_all(bool skip_jacobian = false) const;
  virtual void print(std::ostream& Os) const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(stacked_radiance_mixin, StackedRadianceMixin)

// List of things "import *" will include
%python_export("StackedRadianceMixin");

