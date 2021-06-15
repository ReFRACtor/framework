// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "ground_piecewise.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground_imp_base)
%import "double_with_unit.i"
%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::GroundPiecewise);

namespace FullPhysics {
class GroundPiecewise: public GroundImpBase {
public:
  GroundPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                  const blitz::Array<double, 1>& point_values,
                  const boost::shared_ptr<StateMapping>& mapping = boost::make_shared<StateMappingLinear>());
  virtual boost::shared_ptr<Ground> clone() const = 0;
  virtual const ArrayWithUnit<double, 1>& spectral_points() const;
  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const AutoDerivative<double> value_at_point(const DoubleWithUnit wave_point) const;
  %pickle_serialization();
};
}

