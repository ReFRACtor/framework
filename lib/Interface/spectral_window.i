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
class SpectralWindow  : public GenericObject {
public:
  virtual ~SpectralWindow();
  std::string print_to_string() const;
  SpectralDomain apply(const SpectralDomain& Grid, int Spec_index) const;
  Spectrum apply(const Spectrum& Spec, int Spec_index) const;
  virtual std::vector<int> grid_indexes(const SpectralDomain& Grid, 
					int Spec_index) const = 0;
  %python_attribute(number_spectrometer, virtual int);
  %python_attribute(spectral_bound, virtual SpectralBound);
  %pickle_serialization();
};
}
