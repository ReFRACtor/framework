#ifndef GROUND_EMISSIVITY_POLYNOMIAL_H
#define GROUND_EMISSIVITY_POLYNOMIAL_H

#include "ground_imp_base.h"
#include "array_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a emissivity as a ground type using
  a polynomial per band to model the change in emissivity
  over the channels.
*******************************************************************/
class GroundEmissivityPolynomial: virtual public GroundImpBase {
public:

  GroundEmissivityPolynomial(const blitz::Array<double, 2>& Spec_coeffs,
			     const blitz::Array<bool, 2>& Flag,
			     const ArrayWithUnit<double, 1>& Ref_points,
			     const std::vector<std::string>& Desc_band_names);

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const AutoDerivative<double> emissivity(const DoubleWithUnit wave_point, const int spec_index) const;

  virtual int number_spectrometer() const
  {
    return desc_band_names.size();
  }
  virtual int number_params() const
  {
    return coefficient().value().rows() / number_spectrometer();
  }

  virtual const ArrayAd<double, 1> emiss_coefficients(const int spec_index) const;
  virtual const blitz::Array<double, 2> emiss_covariance(const int spec_index) const;

  /// Center wavelength that spectrally dependent parameter is referenced to
  virtual const DoubleWithUnit reference_point(const int spec_index) const
  {
    return reference_points(spec_index);
  }

  virtual boost::shared_ptr<Ground> clone() const;

  virtual std::string sub_state_identifier() const
  {
    return "ground/emissivity_polynomial";
  }

  virtual std::string state_vector_name_i(int i) const;

  virtual void print(std::ostream& Os) const;

  virtual std::string desc() const
  {
    return "GroundEmissivityPolynomial";
  }

protected:

  GroundEmissivityPolynomial(const blitz::Array<double, 1>& Spec_coeffs,
			     const blitz::Array<bool, 1>& Flag,
			     const ArrayWithUnit<double, 1>& Ref_points,
			     const std::vector<std::string>& Desc_band_names);

  GroundEmissivityPolynomial() {}
private:
  ArrayWithUnit<double, 1> reference_points;
  std::vector<std::string> desc_band_names;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(GroundEmissivityPolynomial);
#endif
