// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "spectral_window.h"
%}

%base_import(generic_object)
%import "spectral_domain.i"
%import "spectral_bound.i"
%import "spectrum.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::SpectralWindow);

namespace FullPhysics {
// Allow these classes to be derived from in Python.
%feature("director") SpectralWindow;
  
class SpectralWindow  : public GenericObject {
public:
  virtual ~SpectralWindow();
  virtual std::string desc() const;
  std::string print_to_string() const;
  SpectralDomain apply(const SpectralDomain& Grid, int Spec_index) const;
  Spectrum apply(const Spectrum& Spec, int Spec_index) const;
  virtual std::vector<int> grid_indexes(const SpectralDomain& Grid, 
					int Spec_index) const = 0;
  %python_attribute_abstract(number_spectrometer, int);
  %python_attribute_abstract(spectral_bound, SpectralBound);
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(spectral_window, SpectralWindow);

// List of things "import *" will include
%python_export("SpectralWindow");
