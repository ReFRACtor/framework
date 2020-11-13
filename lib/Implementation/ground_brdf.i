%include "fp_common.i"

%{
#include "ground_brdf.h"
#include "array_with_unit.h"
#include "double_with_unit.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground_imp_base)
%import "array_with_unit.i"
%import  "double_with_unit.i"

%fp_shared_ptr(FullPhysics::GroundBrdf);
%fp_shared_ptr(FullPhysics::GroundBrdfVeg);
%fp_shared_ptr(FullPhysics::GroundBrdfSoil);

namespace FullPhysics {
%feature("notabstract") GroundBrdfVeg;
%feature("notabstract") GroundBrdfSoil;
class GroundBrdf: public GroundImpBase {
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
  virtual int number_spectrometer() const;
  virtual const int number_weight_parameters() const;
  virtual const AutoDerivative<double> weight(double wn, int spec_index) const;
  virtual const AutoDerivative<double> weight_intercept(int spec_index) const;
  virtual const AutoDerivative<double> weight_slope(int spec_index) const;
  virtual const ArrayAd<double, 1> weight_parameters(const int spec_index) const;
  virtual const AutoDerivative<double> rahman_factor(int spec_index) const;
  virtual const AutoDerivative<double> hotspot_parameter(int spec_index) const;
  virtual const AutoDerivative<double> asymmetry_parameter(int spec_index) const;
  virtual const AutoDerivative<double> anisotropy_parameter(int spec_index) const;
  virtual const AutoDerivative<double> breon_factor(int spec_index) const;

  const blitz::Array<double, 2> brdf_covariance(int spec_index) const;
  virtual double refractive_index(int Spec_idx) const;
  virtual double black_sky_albedo(int Spec_index, double Sza) = 0;
  virtual double kernel_value(int Spec_index, double Sza, double Vza,
                              double Azm) = 0;
  virtual const std::string breon_type() const = 0;
  virtual DoubleWithUnit reference_point(int spec_index) const;
  virtual boost::shared_ptr<Ground> clone() const = 0;
  virtual std::string sub_state_identifier() const;
  virtual std::string state_vector_name_i(int i) const;
  virtual std::string desc() const;
};

class GroundBrdfVeg: public GroundBrdf {
public:
  GroundBrdfVeg(const blitz::Array<double, 2>& Coeffs,
                const ArrayWithUnit<double, 1>& Ref_points, 
                const std::vector<std::string>& Desc_band_names,
                boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());
  virtual double black_sky_albedo(int Spec_index, double Sza);
  static double kernel_value_at_params(const blitz::Array<double, 1>& params, double Sza, double Vza, double Azm);
  virtual double kernel_value(int Spec_index, double Sza, double Vza, double Azm);
  virtual const std::string breon_type() const;
  virtual boost::shared_ptr<Ground> clone() const;
  %pickle_serialization();
};

class GroundBrdfSoil: public GroundBrdf {
public:
  GroundBrdfSoil(const blitz::Array<double, 2>& Coeffs,
                 const ArrayWithUnit<double, 1>& Ref_points,
                 const std::vector<std::string>& Desc_band_names,
                 boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());
  virtual double black_sky_albedo(int Spec_index, double Sza);
  static double kernel_value_at_params(const blitz::Array<double, 1>& params, double Sza, double Vza, double Azm);
  virtual double kernel_value(int Spec_index, double Sza, double Vza, double Azm);
  virtual const std::string breon_type() const;
  virtual boost::shared_ptr<Ground> clone() const;
  %pickle_serialization();
};

}
