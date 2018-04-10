// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "simple_fixed_spectrum_sampling.h"
%}
%base_import(spectrum_sampling)
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::SimpleFixedSpectrumSampling);
namespace FullPhysics {
class SimpleFixedSpectrumSampling : public SpectrumSampling {
public:
  SimpleFixedSpectrumSampling(double wn_start, double wn_end, double wn_step);
  SimpleFixedSpectrumSampling(double wn_start1, double wn_end1, double wn_step1,
			  double wn_start2, double wn_end2, double wn_step2,
			  double wn_start3, double wn_end3, double wn_step3);
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const;
};
}
