%include "fp_common.i"
%{
#include "uniform_spectrum_sampling.h"
%}

%base_import(spectrum_sampling)
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::UniformSpectrumSampling);

namespace FullPhysics {
// I have no idea why it thinks spectral_domain below might not be implemented, force it to be non abstract
%feature("notabstract") UniformSpectrumSampling;

class UniformSpectrumSampling : public SpectrumSampling {
public:
  UniformSpectrumSampling(const ArrayWithUnit<double, 1>& Spec_spacing);

  virtual SpectralDomain spectral_domain(int spec_index,
					 const SpectralDomain& Lowres_grid, 
					 const DoubleWithUnit& Edge_extension) const;

  virtual void print(std::ostream& Os) const;
  %pickle_serialization();
};
}
