%include "fp_common.i"
%{
#include "instrument_doppler.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}

%base_import(spectrum_effect_imp_base);

%fp_shared_ptr(FullPhysics::InstrumentDoppler);

namespace FullPhysics {
class InstrumentDoppler : public SpectrumEffectImpBase {
public:
  InstrumentDoppler(const DoubleWithUnit& Relative_velocity);

  InstrumentDoppler(const double Relative_velocity_value,
                    const std::string& Relative_velocity_units);

  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const;

  virtual boost::shared_ptr<SpectrumEffect> clone() const;

  %python_attribute(sub_state_identifier, std::string);

  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;

  virtual std::string name() const;
  %pickle_serialization();
};
}
