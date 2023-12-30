// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "fp_common.i"

%{
#include "forward_model.h"
%}

%base_import(stacked_radiance_mixin)
%import "spectrum.i"
%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::ForwardModel)


namespace FullPhysics {
  
// Allow these classes to be derived from in Python.
%feature("director") ForwardModel;

class ForwardModel : public StackedRadianceMixin {
public:
  %python_attribute_abstract(num_channels, int);
  virtual std::string desc() const;
  std::string print_to_string() const;
  virtual SpectralDomain spectral_domain(int sensor_index) const = 0;
  boost::optional<blitz::Range> stacked_pixel_range(int sensor_index) const;
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false) const = 0;
  virtual SpectralDomain spectral_domain_all() const;
  virtual Spectrum radiance_all(bool skip_jacobian = false) const;
  virtual void setup_grid() = 0;
  virtual SpectralDomain::TypePreference spectral_domain_type_preference() const = 0;
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >);
  %pickle_serialization();
};
}

%template(Vector_ForwardModel) std::vector<boost::shared_ptr<FullPhysics::ForwardModel> >;

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(forward_model, ForwardModel)

// List of things "import *" will include
%python_export("ForwardModel");
