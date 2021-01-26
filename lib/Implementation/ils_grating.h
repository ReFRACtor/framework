#ifndef ILS_GRATING_H
#define ILS_GRATING_H
#include "ils_imp_base.h"
#include "ils_function.h"
#include "sample_grid.h"

namespace FullPhysics {
/****************************************************************//**
  This is a ILS where we use a Dispersion object to determine the
  wavenumbers of each pixel, and convolve against a IlsFunction.

  This class is only suitable for use by grating spectrometers.
*******************************************************************/

class IlsGrating : virtual public IlsImpBase {
public:

  //-----------------------------------------------------------------------
  /// Constructor.
  //-----------------------------------------------------------------------

  IlsGrating(const boost::shared_ptr<SampleGrid>& Sample_grid,
	     const boost::shared_ptr<IlsFunction>& Ils_func,
	     const DoubleWithUnit& Ils_half_width = DoubleWithUnit(20, units::inv_cm));

  virtual ~IlsGrating() = default;

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
  IlsGrating() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(IlsGrating);
#endif
