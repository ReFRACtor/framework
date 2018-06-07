// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "level_1b_sample_coefficient.h"
%}

%base_import(level_1b)
%import "array_with_unit.i"
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::Level1bSampleCoefficient);

namespace FullPhysics {
class Level1bSampleCoefficient : public Level1b {
public:
  virtual ~Level1bSampleCoefficient();
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const = 0;
  virtual SpectralDomain sample_spectral_domain(int Spec_index) const;
};
}
