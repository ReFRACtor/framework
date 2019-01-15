// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "absco_aer.h"
%}
%base_import(absco)
%import "spectral_bound.i"

%fp_shared_ptr(FullPhysics::AbscoAer);
namespace FullPhysics {

// Force to be not abstract
%feature("notabstract") AbscoAer;

class AbscoAer : public Absco {
public:
  enum InterpolationType {THROW_ERROR_IF_NOT_ON_WN_GRID, NEAREST_NEIGHBOR_WN,
			  INTERPOLATE_WN};
  AbscoAer(const std::string& Fname, double Table_scale = 1.0, 
	   int Cache_nline = 5000,
	   InterpolationType Itype = THROW_ERROR_IF_NOT_ON_WN_GRID);
  AbscoAer(const std::string& Fname, 
	   const SpectralBound& Spectral_bound,
	   const std::vector<double>& Table_scale,
	   int Cache_nline = 5000,
	   InterpolationType Itype = THROW_ERROR_IF_NOT_ON_WN_GRID);
  void load_file(const std::string& Fname);
  void load_file(const std::string& Fname, double Table_scale,
		 int Cache_nline = 5000);
  void load_file(const std::string& Fname, 
		 const SpectralBound& Spectral_bound,
		 const std::vector<double>& Table_scale,
		 int Cache_nline = 5000);
  %python_attribute_derived(pressure_grid, blitz::Array<double, 1>);
  %python_attribute_derived(temperature_grid, blitz::Array<double, 2>);
  %python_attribute(wavenumber_grid, blitz::Array<double, 1>);
  %python_attribute(file_name, std::string);
  %python_attribute_with_set(interpolation_type, InterpolationType);
  void wn_extent(double Wn_in, double& OUTPUT, double& OUTPUT) const;
  bool have_data(double wn) const;
};
}

