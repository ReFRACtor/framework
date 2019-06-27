#include "ground_emissivity_piecewise.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

GroundEmissivityPiecewise::GroundEmissivityPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                                                     const blitz::Array<double, 1>& emissivity_values,
                                                     const blitz::Array<bool, 1>& retrieval_flag)
{
    wavenumbers.reference(spectral_points.convert_wave(units::inv_cm).value);
    init(emissivity_values, retrieval_flag);
    update_sub_state_hook();
}


/// Compute surface parameter array used by radiative transfer
ArrayAd<double, 1> GroundEmissivityPiecewise::surface_parameter(const double wn, const int spec_index) const
{
    AutoDerivative<double> wn_emiss = emissivity(wn);
    ArrayAd<double, 1> spars(1, wn_emiss.number_variable());
    spars(0) = wn_emiss;
    return spars;
}

/// Return emissivity value by querying with an spectral point in arbitrary units
const AutoDerivative<double> GroundEmissivityPiecewise::emissivity(const DoubleWithUnit wave_point) const
{
    return emissivity(wave_point.convert_wave(units::inv_cm).value);
}

/// Return emissivity value by querying by wavenumber
const AutoDerivative<double> GroundEmissivityPiecewise::emissivity(const double wn) const
{
    return (*emiss_interp)(wn);
}

void GroundEmissivityPiecewise::update_sub_state_hook()
{
    typedef LinearInterpolate<double, AutoDerivative<double> > interp_type;

    std::vector<AutoDerivative<double> > emiss_coeff_list;
    for(int idx = 0; idx < coefficient().rows(); idx++) {
        emiss_coeff_list.push_back(coefficient()(idx));
    }

    emiss_interp.reset( new interp_type(wavenumbers.begin(), wavenumbers.end(), emiss_coeff_list.begin()) );
}

boost::shared_ptr<Ground> GroundEmissivityPiecewise::clone() const
{
    return boost::shared_ptr<Ground>(new GroundEmissivityPiecewise(ArrayWithUnit<double, 1>(wavenumbers, units::inv_cm), coefficient().value(), used_flag));
}

std::string GroundEmissivityPiecewise::state_vector_name_i(int i) const
{
    return "Ground Emissivity Piecewise at wavenumber " + boost::lexical_cast<std::string>(wavenumbers(i));
}

void GroundEmissivityPiecewise::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "GroundEmissivityPiecewise:" << std::endl;

    opad << wavenumbers.rows() << " emissivity coefficients (" << wavenumbers(0) << "-" << wavenumbers(wavenumbers.rows()-1) << " cm^-1)" << std::endl;

    opad.strict_sync();
}
