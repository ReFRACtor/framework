// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "ground_emissivity_polynomial.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground_imp_base)
%import "double_with_unit.i"
%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::GroundEmissivityPolynomial);
namespace FullPhysics {
class GroundEmissivityPolynomial: public GroundImpBase {
public:
  GroundEmissivityPolynomial(const blitz::Array<double, 2>& Spec_coeffs,
			     const ArrayWithUnit<double, 1>& Ref_points,
			     const std::vector<std::string>& Desc_band_names);

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
  virtual const AutoDerivative<double> emissivity(const DoubleWithUnit wave_point, const int spec_index) const;
  virtual const int number_spectrometer() const;
  virtual const int number_params() const;
  virtual const ArrayAd<double, 1> emiss_coefficients(const int spec_index) const;
  virtual const blitz::Array<double, 2> emiss_covariance(const int spec_index) const;
  virtual const DoubleWithUnit reference_point(const int spec_index) const;
  virtual boost::shared_ptr<Ground> clone() const;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;
  virtual std::string desc() const;
  %pickle_serialization();
};
}
