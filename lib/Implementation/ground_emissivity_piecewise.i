%include "fp_common.i"

%{
#include "ground_emissivity_piecewise.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground)
%base_import(sub_state_vector_array)
%import "double_with_unit.i"
%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::GroundEmissivityPiecewise);

namespace FullPhysics {
class GroundEmissivityPiecewise: public SubStateVectorArray<Ground> {
public:
    GroundEmissivityPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                              const blitz::Array<double, 1>& emissivity_values,
                              const blitz::Array<bool, 1>& retrieval_flag);

    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

    virtual const AutoDerivative<double> emissivity(const DoubleWithUnit wave_point) const;
    virtual const AutoDerivative<double> emissivity(const double wn) const;

    virtual void update_sub_state_hook();

    virtual boost::shared_ptr<Ground> clone() const;

    %python_attribute(sub_state_identifier, std::string);

    virtual std::string state_vector_name_i(int i) const;

    virtual void print(std::ostream& Os) const;

    virtual std::string desc() const;
};
}

