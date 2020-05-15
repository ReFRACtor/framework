#include "ground_piecewise.h"

using namespace FullPhysics;
using namespace blitz;

GroundPiecewise::GroundPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                                 const blitz::Array<double, 1>& point_values,
                                 const blitz::Array<bool, 1>& retrieval_flag)
{
    spectral_points_ = spectral_points;
    init(point_values, retrieval_flag);
    update_sub_state_hook();
}


const ArrayWithUnit<double, 1>& GroundPiecewise::spectral_points() const
{
    return spectral_points_;
}

/// Compute surface parameter array used by radiative transfer
ArrayAd<double, 1> GroundPiecewise::surface_parameter(const double wn, const int spec_index) const
{
    AutoDerivative<double> wn_value = value_at_point(DoubleWithUnit(wn, units::inv_cm));
    ArrayAd<double, 1> spars(1, wn_value.number_variable());
    spars(0) = wn_value;
    return spars;
}

/// Return value by querying with an spectral point in arbitrary units
const AutoDerivative<double> GroundPiecewise::value_at_point(const DoubleWithUnit wave_point) const
{
    double point_val = wave_point.convert_wave(spectral_points_.units).value;
    return (*ground_interp)(point_val);
}

void GroundPiecewise::update_sub_state_hook()
{
    typedef LinearInterpolate<double, AutoDerivative<double> > interp_type;

    std::vector<AutoDerivative<double> > coeff_list;
    for(int idx = 0; idx < coefficient().rows(); idx++) {
        coeff_list.push_back(coefficient()(idx));
    }

    ground_interp.reset( new interp_type(spectral_points_.value.begin(), spectral_points_.value.end(), coeff_list.begin()) );
}
