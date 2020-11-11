%include "fp_common.i"

%{
#include "ground_emissivity_piecewise.h"
%}

%base_import(ground_piecewise)

%fp_shared_ptr(FullPhysics::GroundEmissivityPiecewise);

namespace FullPhysics {
class GroundEmissivityPiecewise: public GroundPiecewise {
public:
  GroundEmissivityPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                            const blitz::Array<double, 1>& point_values,
                            const boost::shared_ptr<StateMapping>& mapping = boost::make_shared<StateMappingLinear>());
  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const AutoDerivative<double> value_at_point(const DoubleWithUnit wave_point) const;

  virtual boost::shared_ptr<Ground> clone() const;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;
  virtual std::string desc() const;
  %pickle_serialization();
};
}

