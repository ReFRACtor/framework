#ifndef ILS_IMP_BASE_H
#define ILS_IMP_BASE_H

#include "double_with_unit.h"
#include "sample_grid.h"
#include "ils.h"

namespace FullPhysics {

/****************************************************************//**
 Base class to inherit from to simplify the task of implementing 
 most of the necessary Ils behavior. Inheriting classes need to 
 define the apply_ils methods.
*******************************************************************/

class IlsImpBase : public Ils, public Observer<SampleGrid> {
public:
  //-----------------------------------------------------------------------
  /// Constructor.
  //-----------------------------------------------------------------------
  IlsImpBase(const boost::shared_ptr<SampleGrid>& Sample_grid, const DoubleWithUnit& High_res_extension)
    : sample_grid_(Sample_grid), high_res_extension_(High_res_extension)
  {
    sample_grid_->add_observer(*this);
  }

  virtual ~IlsImpBase() = default;

  virtual void notify_update(const SampleGrid& UNUSED(D))
  {
    notify_update_do(*this);
  }

  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;

  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;

  virtual void print(std::ostream& Os) const
  {
    Os << "IlsImpBase";
  }

  virtual std::string band_name() const
  {
    return desc_band_name_;
  }
  virtual std::string hdf_band_name() const
  {
    return hdf_band_name_;
  }
  virtual SpectralDomain pixel_grid() const
  {
    return sample_grid_->pixel_grid();
  }
  virtual DoubleWithUnit high_res_extension() const
  {
    return high_res_extension_;
  }
  virtual void high_res_extension(const DoubleWithUnit& extension)
  {
    high_res_extension_ = extension;
  }

  //-----------------------------------------------------------------------
  /// Underlying SampleGrid object
  //-----------------------------------------------------------------------
  boost::shared_ptr<SampleGrid> sample_grid() const {return sample_grid_; }

  virtual boost::shared_ptr<Ils> clone() const = 0;

protected:
  IlsImpBase() {}
private:
  boost::shared_ptr<SampleGrid> sample_grid_;
  DoubleWithUnit high_res_extension_;
  std::string desc_band_name_, hdf_band_name_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  Handle the degenerate case when we don't want to apply a ILS, this
  just returns the high resolution spectrum without change.
*******************************************************************/
  
class IdentityIls: public IlsImpBase {
public:
  IdentityIls(const boost::shared_ptr<SampleGrid>& Sample_grid)
    : IlsImpBase(Sample_grid, DoubleWithUnit(0, units::nm)) {}
  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& UNUSED(High_resolution_wave_number),
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& UNUSED(Pixel_list)) const 
  { return High_resolution_radiance; }
  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& UNUSED(High_resolution_wave_number),
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& UNUSED(Pixel_list)) const
  { return High_resolution_radiance; }
  virtual boost::shared_ptr<Ils> clone() const
  { return boost::make_shared<IdentityIls>(sample_grid()); }
private:
  IdentityIls() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};  
}
FP_EXPORT_KEY(IdentityIls);
FP_EXPORT_KEY(IlsImpBase);
#endif
