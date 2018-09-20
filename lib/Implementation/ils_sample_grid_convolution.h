#ifndef ILS_SAMPLE_GRID_CONVOLUTION_H
#define ILS_SAMPLE_GRID_CONVOLUTION_H
#include "ils.h"
#include "ils_function.h"
#include "sample_grid.h"

namespace FullPhysics {
  class HdfFile;
/****************************************************************//**
  This is a ILS where we use a SampleGrid object to determine the
  wavenumbers of each pixel, and convolve against a IlsFunction.
*******************************************************************/

class IlsSampleGridConvolution : public Ils, public Observer<SampleGrid> {
public:
//-----------------------------------------------------------------------
/// Constructors.
//-----------------------------------------------------------------------
    IlsSampleGridConvolution(const boost::shared_ptr<SampleGrid>& Sample_grid,
                             const boost::shared_ptr<IlsFunction>& Ils_func,
                             const DoubleWithUnit&
                             Ils_half_width = DoubleWithUnit(20, units::inv_cm))
    : sample_grid_(Sample_grid), ils_func(Ils_func), ils_half_width_(Ils_half_width)  {
        sample_grid_->add_observer(*this);
  }
    IlsSampleGridConvolution(const boost::shared_ptr<SampleGrid>& Sample_grid,
                             const boost::shared_ptr<IlsFunction>& Ils_func,
                             double Ils_half_width)
    : sample_grid_(Sample_grid), ils_func(Ils_func), ils_half_width_(Ils_half_width, units::inv_cm) {
        sample_grid_->add_observer(*this);
    }
  virtual ~IlsSampleGridConvolution() {}

  virtual void notify_update(const SampleGrid& S)  { notify_update_do(*this); }
  virtual blitz::Array<double, 1> apply_ils(const blitz::Array<double, 1>& High_resolution_wave_number,
                                            const blitz::Array<double, 1>& High_resolution_radiance,
                                            const std::vector<int>& Pixel_list) const;
  virtual ArrayAd<double, 1> apply_ils(const blitz::Array<double, 1>& High_resolution_wave_number,
                                       const ArrayAd<double, 1>& High_resolution_radiance,
                                       const std::vector<int>& Pixel_list) const;
  virtual void print(std::ostream& Os) const;
  virtual std::string band_name() const { return ils_func->band_name(); }
  virtual std::string hdf_band_name() const { return ils_func->hdf_band_name(); }
  virtual SpectralDomain pixel_grid() const { return sample_grid_->sample_grid(); }
  virtual DoubleWithUnit ils_half_width() const  {return ils_half_width_;}
  virtual void ils_half_width(const DoubleWithUnit& half_width) { ils_half_width_ = half_width; }
  virtual boost::shared_ptr<Ils> clone() const;

//-----------------------------------------------------------------------
/// Underlying IlsFunction
//-----------------------------------------------------------------------
  boost::shared_ptr<IlsFunction> ils_function() const {return ils_func; }

//-----------------------------------------------------------------------
/// Underlying SampleGrid.
//-----------------------------------------------------------------------
  boost::shared_ptr<SampleGrid> sample_grid() const {return sample_grid_; }

private:
  boost::shared_ptr<SampleGrid> sample_grid_;
  boost::shared_ptr<IlsFunction> ils_func;
  DoubleWithUnit ils_half_width_;
  double integrate(const blitz::Array<double, 1>& x, 
		           const blitz::Array<double, 1>& y) const;
  AutoDerivative<double> integrate(const blitz::Array<double, 1>& x, 
				                   const ArrayAd<double, 1>& y) const;
};
}
#endif
