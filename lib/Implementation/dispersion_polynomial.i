// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "dispersion_polynomial.h"
%}
%base_import(sub_state_vector_array)
%base_import(sample_grid)
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::DispersionPolynomial)
namespace FullPhysics {
class DispersionPolynomial: public SubStateVectorArray<SampleGrid> {
public:
  DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
		       const blitz::Array<bool, 1>& Used_flag,
		       const Unit& Coeff_unit,
		       const std::string& Band_name,
		       int Number_pixel,
		       bool Is_one_based);
  DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
		       const blitz::Array<bool, 1>& Used_flag,
		       const std::string& Coeff_unit_name,
		       const std::string& Band_name,
		       int Number_pixel, bool Is_one_based);
  virtual ~DispersionPolynomial() {}
  %python_attribute(sub_state_identifier, std::string);
  %python_attribute(pixel_grid, SpectralDomain)
  virtual boost::shared_ptr<FullPhysics::SampleGrid> clone() const;
};

}
