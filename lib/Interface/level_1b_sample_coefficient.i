%include "fp_common.i"

%{
#include "level_1b_sample_coefficient.h"
%}

%base_import(level_1b)
%import "array_with_unit.i"
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::Level1bSampleCoefficient);

namespace FullPhysics {

%feature("director") Level1bSampleCoefficient;

class Level1bSampleCoefficient : public Level1b {
public:
  virtual ~Level1bSampleCoefficient();
  virtual int number_sample(int channel_index) const = 0;
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int channel_index)
    const = 0;
  virtual blitz::Array<double, 1> spectral_variable(int channel_index)
    const = 0 ;
  virtual SpectralDomain sample_grid(int channel_index) const;
  %pickle_serialization();
};
}
