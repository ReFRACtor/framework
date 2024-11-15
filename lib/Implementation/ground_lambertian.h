#ifndef GROUND_LAMBERTIAN_H
#define GROUND_LAMBERTIAN_H

#include "ground_imp_base.h"

#include "array_with_unit.h"
#include "state_mapping_linear.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a Lambertian albedo as a ground type. 
*******************************************************************/
class GroundLambertian: virtual public GroundImpBase {
public:
  GroundLambertian(const blitz::Array<double, 2>& Spec_coeffs,
                   const ArrayWithUnit<double, 1>& Ref_points,
                   const std::vector<std::string>& Desc_band_names,
                   boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());
  GroundLambertian(const blitz::Array<double, 2>& Spec_coeffs,
                   const ArrayWithUnit<double, 1>& Ref_points,
		   const Unit& Polynomial_unit,
                   const std::vector<std::string>& Desc_band_names,
                   boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());

  virtual SpurrBrdfType spurr_brdf_type() const
  { return SpurrBrdfType::LAMBERTIAN; }

//-----------------------------------------------------------------------
/// Units that a polynomial is in.
//-----------------------------------------------------------------------
  const Unit& polynomial_unit() const { return polynomial_unit_; }
  
  virtual ArrayAd<double, 1> surface_parameter(const double wn,
					       const int spec_index) const;

  virtual const AutoDerivative<double> albedo(const DoubleWithUnit wave_point,
					      const int spec_index) const;

  virtual int number_spectrometer() const { return desc_band_names.size(); }
  virtual int number_params() const
  { return mapping->mapped_state(coeff).value().rows() / number_spectrometer(); }

  virtual const ArrayAd<double, 1> albedo_coefficients(const int spec_index) const;
  virtual const blitz::Array<double, 2> albedo_covariance(const int spec_index) const;

  /// Center wavelength that spectrally dependent parameter is referenced to
  virtual const DoubleWithUnit reference_point(const int spec_index) const { return reference_points(spec_index); }

  virtual boost::shared_ptr<Ground> clone() const;

  virtual std::string sub_state_identifier() const { return "ground/lambertian"; }

  virtual std::string state_vector_name_i(int i) const;

  virtual void print(std::ostream& Os) const;
protected:

  GroundLambertian(const blitz::Array<double, 1>& Spec_coeffs,
                   const ArrayWithUnit<double, 1>& Ref_points,
		   const Unit& Polynomial_unit,
                   const std::vector<std::string>& Desc_band_names,
                   boost::shared_ptr<StateMapping> Mapping);
  GroundLambertian() : polynomial_unit_(units::inv_cm) {}
private:
  
  ArrayWithUnit<double, 1> reference_points;
  std::vector<std::string> desc_band_names;
  Unit polynomial_unit_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(GroundLambertian);
FP_CLASS_VERSION(GroundLambertian, 1);
#endif
