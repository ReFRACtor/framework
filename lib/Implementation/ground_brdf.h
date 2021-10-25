#ifndef GROUND_BREON_H
#define GROUND_BREON_H

#include "ground_imp_base.h"
#include "auto_derivative.h"
#include "array_with_unit.h"
#include "state_mapping.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements the Breon + Rahman ground type.
  The retrieved parameters are for the Rahman type. Refractive
  index is hardcoded at 1.5 as it is hard coded inside of LIDORT
  and l_rad.

  Base class for both Vegetative and Soil types which share the 
  same mechanisms but are implemented differently in RT codes
*******************************************************************/
class GroundBrdf: virtual public GroundImpBase {
public:
    enum ParamIndex {
        RAHMAN_KERNEL_FACTOR_INDEX = 0,
        RAHMAN_OVERALL_AMPLITUDE_INDEX = 1,
        RAHMAN_ASYMMETRY_FACTOR_INDEX = 2,
        RAHMAN_GEOMETRIC_FACTOR_INDEX = 3,
        BREON_KERNEL_FACTOR_INDEX = 4,
        BRDF_WEIGHT_INTERCEPT_INDEX = 5,
        BRDF_WEIGHT_SLOPE_INDEX = 6
    };

  GroundBrdf(const blitz::Array<double, 2>& Coeffs,
             const ArrayWithUnit<double, 1>& Ref_points,
             const std::vector<std::string>& Desc_band_names,
             boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual int number_spectrometer() const { return desc_band_names.size(); }

  virtual int number_weight_parameters() const { return num_weight_params; }

  virtual const AutoDerivative<double> weight(double wn, int spec_index) const;
  virtual const AutoDerivative<double> weight_intercept(int spec_index) const;
  virtual const AutoDerivative<double> weight_slope(int spec_index) const;
  virtual const ArrayAd<double, 1> weight_parameters(const int spec_index) const;

  // Rahman parameters
  virtual const AutoDerivative<double> rahman_factor(int spec_index) const;
  virtual const AutoDerivative<double> hotspot_parameter(int spec_index) const;
  virtual const AutoDerivative<double> asymmetry_parameter(int spec_index) const;
  virtual const AutoDerivative<double> anisotropy_parameter(int spec_index) const;
  virtual const AutoDerivative<double> breon_factor(int spec_index) const;

  const blitz::Array<double, 2> brdf_covariance(int spec_index) const;
   
  /// Returns hard coded value of 1.5 since that is the value hardcoded into LIDORT
  virtual double refractive_index(int UNUSED(Spec_idx)) const { return 1.5; }

  // Uses LIDORT to compute the black sky albedo from the parameters
  virtual double black_sky_albedo(int Spec_index, double Sza) = 0;

  // Computes kernel value using parameters and specified geometry
  virtual double kernel_value(int Spec_index, double Sza, double Vza, double Azm) = 0;
  
  /// String describing which type of Breon surface type, also makes this class abstract
  virtual const std::string breon_type() const = 0;

  /// Center wavelength that spectrally dependent parameter is referenced to
  virtual DoubleWithUnit reference_point(int spec_index) const { return reference_points(spec_index); } 

  virtual boost::shared_ptr<Ground> clone() const = 0;

  virtual std::string sub_state_identifier() const { return "ground/brdf"; }  

  virtual std::string state_vector_name_i(int i) const;
  
  virtual void print(std::ostream& Os) const;
  
  virtual std::string desc() const { return "GroundBrdf"; }

  GroundBrdf(const GroundBrdf&) = default;
protected:

  GroundBrdf(const blitz::Array<double, 1>& Spec_coeffs,
             const ArrayWithUnit<double, 1>& Ref_points,
             const std::vector<std::string>& Desc_band_names,
             boost::shared_ptr<StateMapping> Mapping);

  GroundBrdf() {}

  virtual ArrayAd<double, 1> brdf_parameters_flat() const;

  ArrayWithUnit<double, 1> reference_points;
  std::vector<std::string> desc_band_names;

  // Helper function for routines that call fortran codes
  blitz::Array<double, 1> black_sky_params(int Spec_index);
  blitz::Array<double, 1> kernel_value_params(int Spec_index);

  int num_weight_params;
  int num_coeff;

private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};


class GroundBrdfVeg: virtual public GroundBrdf {
public:
  GroundBrdfVeg(const blitz::Array<double, 2>& Coeffs,
                const ArrayWithUnit<double, 1>& Ref_points,
                const std::vector<std::string>& Desc_band_names,
                boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>()) :
    GroundBrdf(Coeffs, Ref_points, Desc_band_names, Mapping) {}

  virtual SpurrBrdfType spurr_brdf_type() const
  { return SpurrBrdfType::BREONVEG; }
  virtual double black_sky_albedo(int Spec_index, double Sza);

  static double kernel_value_at_params(const blitz::Array<double, 1>& params, double Sza, double Vza, double Azm);
  virtual double kernel_value(int Spec_index, double Sza, double Vza, double Azm);

  virtual const std::string breon_type() const { return "Vegetative"; }

  virtual boost::shared_ptr<Ground> clone() const {
    return boost::shared_ptr<Ground>(new GroundBrdfVeg(brdf_parameters_flat().value(), reference_points, desc_band_names, mapping));
  }
private:
  GroundBrdfVeg(const blitz::Array<double, 1>& Spec_coeffs,
                const ArrayWithUnit<double, 1>& Ref_points,
                const std::vector<std::string>& Desc_band_names,
                boost::shared_ptr<StateMapping> Mapping) :
    GroundBrdf(Spec_coeffs, Ref_points, Desc_band_names, Mapping) {}
  GroundBrdfVeg() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

class GroundBrdfSoil: virtual public GroundBrdf {
public:
  GroundBrdfSoil(const blitz::Array<double, 2>& Coeffs,
                 const ArrayWithUnit<double, 1>& Ref_points,
                 const std::vector<std::string>& Desc_band_names,
                 boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>()) :
    GroundBrdf(Coeffs, Ref_points, Desc_band_names, Mapping) {}

  virtual SpurrBrdfType spurr_brdf_type() const
  { return SpurrBrdfType::BREONSOIL; }
  virtual double black_sky_albedo(int Spec_index, double Sza);

  static double kernel_value_at_params(const blitz::Array<double, 1>& params, double Sza, double Vza, double Azm);
  virtual double kernel_value(int Spec_index, double Sza, double Vza, double Azm);

  virtual const std::string breon_type() const { return "Soil"; }

  virtual boost::shared_ptr<Ground> clone() const {
    return boost::shared_ptr<Ground>(new GroundBrdfSoil(brdf_parameters_flat().value(), reference_points, desc_band_names, mapping));
  }
private:
  GroundBrdfSoil(const blitz::Array<double, 1>& Spec_coeffs,
                 const ArrayWithUnit<double, 1>& Ref_points,
                 const std::vector<std::string>& Desc_band_names,
                 boost::shared_ptr<StateMapping> Mapping) :
    GroundBrdf(Spec_coeffs, Ref_points, Desc_band_names, Mapping) {}
  GroundBrdfSoil() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

} // End of namespace

FP_EXPORT_KEY(GroundBrdf);
FP_EXPORT_KEY(GroundBrdfVeg);
FP_EXPORT_KEY(GroundBrdfSoil);
 #endif
