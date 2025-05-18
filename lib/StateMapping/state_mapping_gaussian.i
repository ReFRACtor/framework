%include "fp_common.i"

%{
#include "state_mapping_gaussian.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMappingGaussian);


namespace FullPhysics {

class StateMappingGaussian : public StateMapping  {
public:
  StateMappingGaussian(const boost::shared_ptr<Pressure>& in_press, bool Linear_AOD,
                       double Min_Desired = 1e-9);
  virtual boost::shared_ptr<StateMapping> clone() const;
  %python_attribute(is_linear_total, bool);
  %python_attribute(min_desired, double);
  %python_attribute(pressure, boost::shared_ptr<Pressure>);
  virtual AutoDerivative<double>
    total_optical_depth(const ArrayAd<double, 1>& component) const;
  %pickle_serialization();
};
}

