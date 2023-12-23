// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"

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
%fp_shared_ptr(FullPhysics::IdentityIls);

%feature("director") FullPhysics::IlsImpBase;

namespace FullPhysics {
class IlsImpBase : public Ils, public Observer<SampleGrid> {
public:
  IlsImpBase(const boost::shared_ptr<SampleGrid>& Sample_grid, const DoubleWithUnit& Edge_extension);
  virtual void notify_update(const SampleGrid& D);
  virtual blitz::Array<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const blitz::Array<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const = 0;
  virtual ArrayAd<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const ArrayAd<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const = 0;
  virtual std::string band_name() const;
  virtual std::string hdf_band_name() const;
  virtual SpectralDomain pixel_grid() const;
  virtual DoubleWithUnit high_res_extension() const;
  virtual void high_res_extension(const DoubleWithUnit& extension);
  boost::shared_ptr<SampleGrid> sample_grid() const;
  virtual boost::shared_ptr<Ils> clone() const = 0;
  virtual std::string desc() const;
  %pickle_serialization();
protected:
  IlsImpBase();
};

class IdentityIls: public IlsImpBase {
public:
  IdentityIls(const boost::shared_ptr<SampleGrid>& Sample_grid);
  virtual blitz::Array<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const blitz::Array<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const;
  virtual ArrayAd<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const ArrayAd<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const;
  virtual boost::shared_ptr<Ils> clone() const;
  %pickle_serialization();
};  
  
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(ils_imp_base, IlsImpBase)

// List of things "import *" will include
%python_export("IlsImpBase", "IdentityIls");
