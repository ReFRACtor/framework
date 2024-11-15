// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "absco_coeff.h"
%}
%base_import(absco)
%import "spectral_bound.i"

%fp_shared_ptr(FullPhysics::AbscoCoeff);
namespace FullPhysics {

// Force to be not abstract
%feature("notabstract") AbscoCoeff;

class AbscoCoeff : public Absco {
public:
  AbscoCoeff(const std::string& Fname, double Table_scale = 1.0, 
	     int Cache_nline = 5000);
  AbscoCoeff(const std::string& Fname, 
	     const SpectralBound& Spectral_bound,
	     const std::vector<double>& Table_scale,
	     int Cache_nline = 5000);
  void load_file(const std::string& Fname);
  void load_file(const std::string& Fname, double Table_scale,
		 int Cache_nline = 5000);
  void load_file(const std::string& Fname, 
		 const SpectralBound& Spectral_bound,
		 const std::vector<double>& Table_scale,
		 int Cache_nline = 5000);
  %python_attribute_derived(pressure_grid, blitz::Array<double, 1>)
  %python_attribute_derived(temperature_grid, blitz::Array<double, 2>)
  %python_attribute(wavenumber_grid, blitz::Array<double, 1>)
  %python_attribute(file_name, std::string)
  void wn_extent(double Wn_in, double& OUTPUT, double& OUTPUT) const;
  bool have_data(double wn) const;
  void wn_index(double Wn_in, int& OUTPUT, double& OUTPUT) const;
  %pickle_serialization();
};
}

