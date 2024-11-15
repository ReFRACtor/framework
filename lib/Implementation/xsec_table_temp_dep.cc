#include "xsec_table_temp_dep.h"
#include "fp_serialize_support.h"
#include "fp_logger.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void XSecTableTempDep::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(XSecTableImpBase);
}

FP_IMPLEMENT(XSecTableTempDep);

#endif

XSecTableTempDep::XSecTableTempDep(const ArrayWithUnit<double, 1>& Spectral_grid, const ArrayWithUnit<double, 2>& XSec_values, double Conversion_factor)
: XSecTableImpBase(Spectral_grid, XSec_values, Conversion_factor)
{
}

XSecTableTempDep::XSecTableTempDep(const ArrayWithUnit<double, 1>& Spectral_grid, const std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > Data_interp, const Unit& Data_units, double Conversion_factor)
: XSecTableImpBase(Spectral_grid, Data_interp, Data_units, Conversion_factor)
{
}

ArrayAdWithUnit<double, 1> XSecTableTempDep::optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAdWithUnit<double, 1> gas_density_levels, ArrayAdWithUnit<double, 1> temperature_levels) const
{
    if (temperature_levels.units != units::K) {
        Exception err;
        err << "Expected temperature levels to be in units of K instead of: " << temperature_levels.units;
        throw err;
    }
    
    DoubleWithUnit xsec_value = cross_section_value(spectral_point);
    ArrayWithUnit<double, 1> coeffs = cross_section_coefficients(spectral_point);

    if (coeffs.rows() < 2) {
        Exception err;
        err << "Only " << coeffs.rows() << " cross section coefficients returned, need at least 2";
        throw err;
    } else if(coeffs.rows() > 2) {
        Logger::warning() << "Temperature dependant cross section data only needs 2 coefficients, but data with " << coeffs.rows() << " passed.";
    }

    double c1 = coeffs.value(0);
    double c2 = coeffs.value(1);

    ArrayAd<double, 1> gas_od(gas_density_levels.rows() - 1, gas_density_levels.number_variable());

    for (int lay_idx = 0; lay_idx < gas_od.rows(); lay_idx++) {
        // Convert K to Celsius
        AutoDerivative<double> temp1 = temperature_levels.value(lay_idx) - 273.15;
        AutoDerivative<double> temp2 = temperature_levels.value(lay_idx+1) - 273.15;

        AutoDerivative<double> xsec1 = xsec_value.value + temp1 * c1 + temp1*temp1 * c2;
        AutoDerivative<double> xsec2 = xsec_value.value + temp2 * c1 + temp2*temp2 * c2;

        gas_od(lay_idx) = (gas_density_levels.value(lay_idx) * xsec1 + gas_density_levels.value(lay_idx+1) * xsec2) * 0.5;
    }

    return ArrayAdWithUnit<double, 1>(gas_od, gas_density_levels.units * xsec_value.units);
}

boost::shared_ptr<XSecTable> XSecTableTempDep::clone() const
{
    return boost::shared_ptr<XSecTable>(new XSecTableTempDep(spectral_grid_values, data_interp, xsec_units, conversion_factor));
}
