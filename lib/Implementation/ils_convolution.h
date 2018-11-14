#ifndef ILS_CONVOLUTION_H
#define ILS_CONVOLUTION_H
#include "ils_imp_base.h"
#include "ils_function.h"
#include "sample_grid.h"

namespace FullPhysics {
/****************************************************************//**
  This is a ILS where we use a Dispersion object to determine the
  wavenumbers of each pixel, and convolve against a IlsFunction.
*******************************************************************/

class IlsConvolution : public IlsImpBase {
public:

    //-----------------------------------------------------------------------
    /// Constructor.
    //-----------------------------------------------------------------------

    IlsConvolution(const boost::shared_ptr<SampleGrid>& Sample_grid,
                   const boost::shared_ptr<IlsFunction>& Ils_func,
                   const DoubleWithUnit& Ils_half_width = DoubleWithUnit(20, units::inv_cm));

    //-----------------------------------------------------------------------
    /// Constructor.
    //-----------------------------------------------------------------------

    /*  IlsConvolution(const boost::shared_ptr<SampleGrid>& Disp,
             const boost::shared_ptr<IlsFunction>& Ils_func,
             double Ils_half_width)
        : disp(Disp), ils_func(Ils_func),
          ils_half_width_(Ils_half_width, units::inv_cm)
      { disp->add_observer(*this); }*/

    virtual ~IlsConvolution() = default;

    virtual blitz::Array<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const blitz::Array<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const;
    virtual ArrayAd<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const ArrayAd<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const;

    virtual void print(std::ostream& Os) const;
    virtual boost::shared_ptr<Ils> clone() const;

    //-----------------------------------------------------------------------
    /// Underlying IlsFunction
    //-----------------------------------------------------------------------
    boost::shared_ptr<IlsFunction> ils_function() const
    {
        return ils_func;
    }

private:
    boost::shared_ptr<IlsFunction> ils_func;
    double integrate(const blitz::Array<double, 1>& x,
                     const blitz::Array<double, 1>& y) const;
    AutoDerivative<double> integrate(const blitz::Array<double, 1>& x,
                                     const ArrayAd<double, 1>& y) const;
};
}
#endif
