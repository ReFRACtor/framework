%include <std_vector.i>
%include "fp_common.i"

%{
#include "forward_model.h"
%}

%base_import(stacked_radiance_mixin)
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::ForwardModel)

namespace FullPhysics {

class ForwardModel : public StackedRadianceMixin {
public:
    virtual ~ForwardModel();
    std::string print_to_string() const;
    virtual void setup_grid();
    %python_attribute(num_channels, virtual int)
    virtual SpectralDomain spectral_domain(int channel_index) const;
    virtual Spectrum radiance(int channel_index, bool skip_jacobian = false) const = 0;
};
}

%template(Vector_ForwardModel) std::vector<boost::shared_ptr<FullPhysics::ForwardModel> >;
