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
  Level1bSampleCoefficient(const bool One_based = true);
  virtual ~Level1bSampleCoefficient();
  virtual int number_sample(int Spec_index) const = 0;
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const = 0;
  virtual SpectralDomain sample_grid(int Spec_index) const;
};
}
