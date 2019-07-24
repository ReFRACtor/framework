#include "ground_piecewise.h"

using namespace FullPhysics;
using namespace blitz;

GroundPiecewise::GroundPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                                 const blitz::Array<double, 1>& point_values,
                                 const blitz::Array<bool, 1>& retrieval_flag)
{
    wavenumbers.reference(spectral_points.convert_wave(units::inv_cm).value);
    init(point_values, retrieval_flag);
    update_sub_state_hook();
}


/// Compute surface parameter array used by radiative transfer
ArrayAd<double, 1> GroundPiecewise::surface_parameter(const double wn, const int spec_index) const
{
    AutoDerivative<double> wn_value = value_at_wavenumber(wn);
    ArrayAd<double, 1> spars(1, wn_value.number_variable());
    spars(0) = wn_value;
    return spars;
}

/// Return value by querying with an spectral point in arbitrary units
const AutoDerivative<double> GroundPiecewise::value_at_point(const DoubleWithUnit wave_point) const
{
    return value_at_wavenumber(wave_point.convert_wave(units::inv_cm).value);
}

/// Return value by querying by wavenumber
const AutoDerivative<double> GroundPiecewise::value_at_wavenumber(const double wn) const
{
    return (*ground_interp)(wn);
}

void GroundPiecewise::update_sub_state_hook()
{
    typedef LinearInterpolate<double, AutoDerivative<double> > interp_type;

    std::vector<AutoDerivative<double> > coeff_list;
    for(int idx = 0; idx < coefficient().rows(); idx++) {
        coeff_list.push_back(coefficient()(idx));
    }

    ground_interp.reset( new interp_type(wavenumbers.begin(), wavenumbers.end(), coeff_list.begin()) );
}
