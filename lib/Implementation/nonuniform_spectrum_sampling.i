// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nonuniform_spectrum_sampling.h"
%}
%base_import(spectrum_sampling)
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::NonuniformSpectrumSampling);
namespace FullPhysics {
class NonuniformSpectrumSampling : public SpectrumSampling {
public:
  NonuniformSpectrumSampling(const SpectralDomain& Wn,
     const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling);
  NonuniformSpectrumSampling(const SpectralDomain& Wn1,
     const SpectralDomain& Wn2,
     const SpectralDomain& Wn3,
     const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling);
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const;
};
}
