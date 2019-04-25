#include "ils_convolution.h"
#include "hdf_file.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(IlsConvolution, Ils)
.def(luabind::constructor<const boost::shared_ptr<SampleGrid>&, const boost::shared_ptr<IlsFunction>&>())
.def(luabind::constructor<const boost::shared_ptr<SampleGrid>&, const boost::shared_ptr<IlsFunction>&, DoubleWithUnit&>())
REGISTER_LUA_END()
#endif

IlsConvolution::IlsConvolution(const boost::shared_ptr<SampleGrid>& Sample_grid,
                               const boost::shared_ptr<IlsFunction>& Ils_func,
                               const DoubleWithUnit& Ils_half_width)
: IlsImpBase(Sample_grid, Ils_half_width), ils_func(Ils_func)
{
}

// See base class for description.
blitz::Array<double, 1> IlsConvolution::apply_ils
(const blitz::Array<double, 1>& Hres_wn,
 const blitz::Array<double, 1>& Hres_rad,
 const std::vector<int>& Pixel_list) const
{
    if(Hres_wn.rows() != Hres_rad.rows()) {
        throw Exception("wave_number and radiance need to be the same size");
    }

    Array<double, 1> disp_wn(sample_grid()->pixel_grid().data());
    ArrayAd<double, 1> response;
    Array<double, 1> res((int) Pixel_list.size());

    for(int i = 0; i < res.rows(); ++i) {
        // Find the range of Hres_wn that lie within
        // wn_center +- ils_half_width

        double wn_center = disp_wn(Pixel_list[i]);
        Array<double, 1>::const_iterator itmin =
            std::lower_bound(Hres_wn.begin(), Hres_wn.end(),
                             wn_center - ils_half_width().value);
        Array<double, 1>::const_iterator itmax =
            std::lower_bound(Hres_wn.begin(), Hres_wn.end(),
                             wn_center + ils_half_width().value);
        int jmin = (itmin == Hres_wn.end() ? Hres_wn.rows() - 1 : itmin.position()(0));
        int jmax = (itmax == Hres_wn.end() ? Hres_wn.rows() - 1 : itmax.position()(0));
        Range r(jmin, jmax);

        // Convolve with response

        ils_func->ils(wn_center, Hres_wn(r), response);
        Array<double, 1> conv(response.value() * Hres_rad(r));

        // And integrate to get pixel value.

        res(i) = integrate(Hres_wn(r), conv) /
                 integrate(Hres_wn(r), response.value());
    }

    return res;
}

//-----------------------------------------------------------------------
/// Simple trapezoid rule for integration.
//-----------------------------------------------------------------------

double IlsConvolution::integrate(const blitz::Array<double, 1>& x,
                                 const blitz::Array<double, 1>& y) const
{
    double res = 0;

    for(int i = 1; i < x.rows(); ++i) {
        res += (x(i) - x(i - 1)) * (y(i) + y(i - 1));
    }

    return res / 2.0;
}

//-----------------------------------------------------------------------
/// Simple trapezoid rule for integration.
//-----------------------------------------------------------------------

AutoDerivative<double> IlsConvolution::integrate(const blitz::Array<double, 1>& x,
        const ArrayAd<double, 1>& y) const
{
    double resv = 0;
    Array<double, 1> resgrad(y.number_variable());
    resgrad = 0;

    for(int i = 1; i < x.rows(); ++i) {
        double t = x(i) - x(i - 1);
        resv +=  t * (y.value()(i) + y.value()(i - 1));
        resgrad += t * (y.jacobian()(i, Range::all()) +
                        y.jacobian()(i - 1, Range::all()));
    }

    resv /= 2.0;
    resgrad /= 2.0;
    return AutoDerivative<double>(resv, resgrad);
}


// See base class for description.
ArrayAd<double, 1> IlsConvolution::apply_ils
(const blitz::Array<double, 1>& Hres_wn,
 const ArrayAd<double, 1>& Hres_rad,
 const std::vector<int>& Pixel_list) const
{
    firstIndex i1;
    secondIndex i2;

    if(Hres_wn.rows() != Hres_rad.rows()) {
        throw Exception("wave_number and radiance need to be the same size");
    }

    ArrayAd<double, 1> disp_wn(sample_grid()->pixel_grid().data_ad());
    ArrayAd<double, 1> res((int) Pixel_list.size(),
                           std::max(disp_wn.number_variable(),
                                    Hres_rad.number_variable()));
    // A few scratch variables, defined outside of loop so we don't keep
    // recreating them.
    Array<double, 1> one(1);
    one = 1;
    Array<double, 1> normfact_grad(disp_wn.number_variable());
    ArrayAd<double, 1> conv(1, res.number_variable());
    ArrayAd<double, 1> response;

    for(int i = 0; i < res.rows(); ++i) {
        // Find the range of Hres_wn that lie within
        // wn_center +- ils_half_width

        AutoDerivative<double> wn_center = disp_wn(Pixel_list[i]);
        Array<double, 1>::const_iterator itmin =
            std::lower_bound(Hres_wn.begin(), Hres_wn.end(),
                             wn_center.value() - ils_half_width().value);
        Array<double, 1>::const_iterator itmax =
            std::lower_bound(Hres_wn.begin(), Hres_wn.end(),
                             wn_center.value() + ils_half_width().value);
        int jmin = (itmin == Hres_wn.end() ? Hres_wn.rows() - 1 : itmin.position()(0));
        int jmax = (itmax == Hres_wn.end() ? Hres_wn.rows() - 1 : itmax.position()(0));
        Range r(jmin, jmax);

        // Convolve with response

        // For speed, we just calculate dresponse/dwn_center. This along
        // with the chain rule is enough to calculate the Jacobian, and is
        // faster.
        AutoDerivative<double> wn_center2(wn_center.value(), one);
        ils_func->ils(wn_center2, Hres_wn(r), response);

        // Response may not be normalized, so calculate normalization
        // factor.
        AutoDerivative<double> normfact_dwn = integrate(Hres_wn(r), response);
        // Convert gradient from dwn_center to dState
        normfact_grad = normfact_dwn.gradient()(0) * wn_center.gradient();
        AutoDerivative<double> normfact(normfact_dwn.value(), normfact_grad);

        conv.resize(response.rows(), res.number_variable());
        conv.value() = response.value() * Hres_rad(r).value();

        if(!Hres_rad.is_constant() && !wn_center.is_constant())
            conv.jacobian() = response.value()(i1) * Hres_rad(r).jacobian()(i1, i2) +
                              response.jacobian()(Range::all(), 0)(i1) * wn_center.gradient()(i2) *
                              Hres_rad(r).value()(i1);
        else if(!wn_center.is_constant())
            conv.jacobian() = response.jacobian()(Range::all(), 0)(i1) *
                              wn_center.gradient()(i2) * Hres_rad(r).value()(i1);
        else if(!Hres_rad.is_constant()) {
            conv.jacobian() = response.value()(i1) * Hres_rad(r).jacobian()(i1, i2);
        }

        res(i) = integrate(Hres_wn(r), conv) / normfact;
    }

    return res;
}

boost::shared_ptr<Ils> IlsConvolution::clone() const
{
    return boost::shared_ptr<Ils>(new IlsConvolution(sample_grid()->clone(),
                                  ils_func, ils_half_width()));
}

void IlsConvolution::print(std::ostream& Os) const
{
    Os << "IlsConvolution:\n";
    OstreamPad opad(Os, "  ");
    opad << *sample_grid() << "\n" << *ils_func << "\n"
         << "Half Width: " << ils_half_width() << "\n";
    opad.strict_sync();
}
