#include "xsec_table_simple.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void XSecTableSimple::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(XSecTableImpBase);
}

FP_IMPLEMENT(XSecTableSimple);

#endif

XSecTableSimple::XSecTableSimple(const ArrayWithUnit<double, 1>& Spectral_grid, const Array<double, 2>& XSec_values, double Conversion_factor)
: XSecTableImpBase(Spectral_grid, XSec_values, Conversion_factor)
{
}

XSecTableSimple::XSecTableSimple(const ArrayWithUnit<double, 1>& Spectral_grid, const std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > Data_interp, double Conversion_factor)
: XSecTableImpBase(Spectral_grid, Data_interp, Conversion_factor)
{
}

ArrayAd<double, 1> XSecTableSimple::optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAd<double, 1> gas_density_levels, ArrayAd<double, 1> UNUSED(temperature_levels)) const
{
    double xsec_value = cross_section_value(spectral_point);
    ArrayAd<double, 1> gas_od(gas_density_levels.rows() - 1, gas_density_levels.number_variable());

    for (int lay_idx = 0; lay_idx < gas_od.rows(); lay_idx++) {
        gas_od(lay_idx) = (gas_density_levels(lay_idx) * xsec_value + gas_density_levels(lay_idx-1) * xsec_value) * 0.5;
    }

    return gas_od;
}

boost::shared_ptr<XSecTable> XSecTableSimple::clone() const
{
    return boost::shared_ptr<XSecTable>(new XSecTableSimple(spectral_grid_values, data_interp, conversion_factor));
}
