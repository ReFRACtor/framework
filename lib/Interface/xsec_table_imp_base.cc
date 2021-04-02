#include "xsec_table_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void XSecTableImpBase::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(XSecTable)
    & FP_NVP(spectral_grid_values) & FP_NVP(conversion_factor) & FP_NVP(data_interp);
}

FP_IMPLEMENT(XSecTableImpBase);

#endif

XSecTableImpBase::XSecTableImpBase(const ArrayWithUnit<double, 1>& Spectral_grid, const Array<double, 2>& XSec_values, double Conversion_factor)
: spectral_grid_values(Spectral_grid), conversion_factor(Conversion_factor)
{
    init_interpolation(XSec_values);
}

XSecTableImpBase::XSecTableImpBase(const ArrayWithUnit<double, 1>& Spectral_grid, const std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > Data_interp, double Conversion_factor)
: spectral_grid_values(Spectral_grid), data_interp(Data_interp), conversion_factor(Conversion_factor)
{
}

void XSecTableImpBase::init_interpolation(const blitz::Array<double, 2>& xsec_values)
{
    if(xsec_values.cols() < 1 || xsec_values.rows() < 2) {
        throw Exception("No cross section values defined");
    }

    if(xsec_values.rows() != spectral_grid_values.rows()) {
        Exception err;
        err << "The number of cross section values: " << xsec_values.rows() 
            << " must be the same as the number of spectral grid values: " << spectral_grid_values.rows();
        throw err;
    }

    typedef LinearInterpolate<double, double> lin_interp_type;

    for(int xsec_values_col = 0; xsec_values_col < xsec_values.cols(); xsec_values_col++) {
        Array<double, 1> col_values(xsec_values.rows());
        col_values = xsec_values(Range::all(), xsec_values_col);
        
        boost::shared_ptr<lin_interp_type> col_interp(new lin_interp_type(spectral_grid_values.value.begin(), spectral_grid_values.value.end(), col_values.begin()));
        data_interp.push_back(col_interp);
    }
}

const double XSecTableImpBase::cross_section_value(DoubleWithUnit& spectral_point) const
{
    // Convert to grid units
    double interp_point = spectral_point.convert(spectral_grid_values.units).value;

    return (*data_interp[0])(interp_point) / conversion_factor;
}

const Array<double, 1> XSecTableImpBase::cross_section_coefficients(DoubleWithUnit& spectral_point) const
{
    if (data_interp.size() <= 1) {
        return Array<double, 1>(0);
    }

    // Convert to grid units
    double interp_point = spectral_point.convert(spectral_grid_values.units).value;

    Array<double, 1> coeffs(data_interp.size() - 1);

    for(int coeff_col = 1; coeff_col < data_interp.size(); coeff_col++) {
        coeffs(coeff_col - 1) = (*data_interp[coeff_col])(interp_point) / conversion_factor;
    }

    return coeffs;
}
