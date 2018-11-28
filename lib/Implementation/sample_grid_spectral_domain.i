%include "common.i"

%{
#include "sample_grid_spectral_domain.h"
%}

%base_import(sub_state_vector_array)
%base_import(sample_grid)

%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::SampleGridSpectralDomain)

namespace FullPhysics {
class SampleGridSpectralDomain: public SubStateVectorArray<SampleGrid> {
public:
    SampleGridSpectralDomain(const SpectralDomain& Spec_domain, const std::string& Band_name);
    virtual ~SampleGridSpectralDomain();
    %python_attribute(sub_state_identifier, std::string);
    %python_attribute(sample_grid, SpectralDomain)
    virtual boost::shared_ptr<SampleGrid> clone() const;
};

}
