%include "fp_common.i"

%{
#include "sample_grid_spectral_domain.h"
%}

%base_import(sample_grid_imp_base)

%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::SampleGridSpectralDomain)

namespace FullPhysics {
class SampleGridSpectralDomain: public SampleGridImpBase {
public:
  SampleGridSpectralDomain(const SpectralDomain& Spec_domain, const std::string& Band_name);
  %python_attribute(sub_state_identifier, std::string);
  %python_attribute(sample_grid, SpectralDomain);
  virtual boost::shared_ptr<SampleGrid> clone() const;
  %pickle_serialization();
};

}
