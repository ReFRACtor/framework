%include "fp_common.i"

%{
#include "cloud_3d_effect.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}

%base_import(spectrum_effect_imp_base)

%fp_shared_ptr(FullPhysics::Cloud3dEffect);

namespace FullPhysics {
class Cloud3dEffect : virtual public SpectrumEffectImpBase {
public:
  Cloud3dEffect(const double Offset, const double Slope, const std::string& Band_name, boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());

  %python_attribute(offset, AutoDerivative<double>);
  %python_attribute(slope, AutoDerivative<double>);

  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const;

  virtual boost::shared_ptr<SpectrumEffect> clone() const;

  %python_attribute(sub_state_identifier, std::string);

  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;

  %python_attribute(name, std::string)

  %pickle_serialization();
};
}
