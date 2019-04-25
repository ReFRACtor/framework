// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"

%{
#include "ils_imp_base.h"
#include "sub_state_vector_array.h"
%}

%import "double_with_unit.i"
%import "array_ad.i"
%import "spectral_domain.i"
%import "double_with_unit.i"

%base_import(ils)
%base_import(observer)
%base_import(sample_grid)

%fp_shared_ptr(FullPhysics::IlsImpBase);

%feature("director") FullPhysics::IlsImpBase;

namespace FullPhysics {
class IlsImpBase : public Ils, public Observer<SampleGrid> {
public:
  IlsImpBase(const boost::shared_ptr<SampleGrid>& Sample_grid, const DoubleWithUnit& Ils_half_width);
  virtual void notify_update(const SampleGrid& D);
  virtual blitz::Array<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const blitz::Array<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const = 0;
  virtual ArrayAd<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const ArrayAd<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const = 0;
  virtual const std::string band_name() const;
  virtual const std::string hdf_band_name() const;
  virtual const SpectralDomain pixel_grid() const;
  virtual const DoubleWithUnit ils_half_width() const;
  virtual void ils_half_width(const DoubleWithUnit& half_width);
  boost::shared_ptr<SampleGrid> sample_grid() const {return sample_grid_; }
    virtual boost::shared_ptr<Ils> clone() const = 0;
};
}
