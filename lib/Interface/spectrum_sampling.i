// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "spectrum_sampling.h"
%}
%base_import(generic_object)
%import "spectral_domain.i"
%import "spectrum.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::SpectrumSampling);
%fp_shared_ptr(FullPhysics::IdentitySpectrumSampling);

namespace FullPhysics {
// Allow these classes to be derived from in Python.
%feature("director") SpectrumSampling;
  
class SpectrumSampling : public GenericObject {
public:
  virtual ~SpectrumSampling();
  std::string print_to_string() const;
  %python_attribute(number_spectrometer, int);
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Edge_extension) const = 0;
  virtual SpectralDomain spectral_domain_interpolated(int Spec_index, 
		 const SpectralDomain& Lowres_grid, 
	         const DoubleWithUnit& Edge_extension) const;
  virtual bool need_interpolation(int Spec_index) const;
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >);
  %pickle_serialization();
protected:
  int nspectrometer;
  SpectrumSampling();
  SpectrumSampling(int num_spectrometer);
};
class IdentitySpectrumSampling: public SpectrumSampling {
public:
  IdentitySpectrumSampling(int nspec);
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Edge_extension) const;
  %pickle_serialization();
};  
}
