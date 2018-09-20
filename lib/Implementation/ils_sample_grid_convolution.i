// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "ils_sample_grid_convolution.h"
#include "sub_state_vector_array.h"
%}
%base_import(ils)
%base_import(observer)
%base_import(sample_grid)
%import "double_with_unit.i"
%import "state_vector.i"
%import "array_ad.i"
%import "ils_function.i"
%fp_shared_ptr(FullPhysics::IlsSampleGridConvolution);

namespace FullPhysics {
class IlsSampleGridConvolution : public Ils, public Observer<SampleGrid> {
public:
  IlsSampleGridConvolution(const boost::shared_ptr<SampleGrid>& Sample_grid,
		 const boost::shared_ptr<IlsFunction>& Ils_func,
		 const DoubleWithUnit& Ils_half_width = DoubleWithUnit(20, units::inv_cm));
  IlsSampleGridConvolution(const boost::shared_ptr<SampleGrid>& Sample_grid,
		 const boost::shared_ptr<IlsFunction>& Ils_func,
		 double Ils_half_width);
  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const;
  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const;
  virtual boost::shared_ptr<Ils> clone() const;
  %python_attribute(ils_function, boost::shared_ptr<IlsFunction>);
  %python_attribute(sample_grid, boost::shared_ptr<SampleGrid>);
};
}
