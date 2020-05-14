%include "fp_common.i"

%{
#include "ground_piecewise.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground)
%base_import(sub_state_vector_array)
%import "double_with_unit.i"
%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::GroundPiecewise);

namespace FullPhysics {
class GroundPiecewise: public SubStateVectorArray<Ground> {
public:
    GroundPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                    const blitz::Array<double, 1>& point_values,
                    const blitz::Array<bool, 1>& retrieval_flag);

    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

    virtual const AutoDerivative<double> value_at_point(const DoubleWithUnit wave_point) const;
};
}

