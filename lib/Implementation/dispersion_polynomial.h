#ifndef DISPERSION_POLYNOMIAL_H
#define DISPERSION_POLYNOMIAL_H
#include "sample_grid_imp_base.h"
#include "unit.h"
#include "polynomial_eval.h"
#include "state_mapping.h"

namespace FullPhysics {
/****************************************************************//**
  This is an implementation of Dispersion that uses a polynomial
  expression to calculate the wavenumbers.

  Note that there are two minor variations of the dispersion
  polynomial. The first wavenumber returned can either be the
  polynomial evaluated at the value of "1", or a value of "0". By
  convention, the polynomial is 1 based for GOSAT and OCO, but 0 based
  for FTS.
*******************************************************************/
class DispersionPolynomial: virtual public SampleGridImpBase {
public:
  DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
                       const Unit& Coeff_unit,
                       const blitz::Array<double, 1>& Var_values,
                       const std::string& Band_name,
                       boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());
  DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
                       const std::string& Coeff_unit_name,
                       const blitz::Array<double, 1>& Var_values,
                       const std::string& Band_name,
                       boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());
  DispersionPolynomial(const ArrayWithUnit<double, 1>& Coeff, 
                       const blitz::Array<double, 1>& Var_values,
                       const std::string& Band_name,
                       boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());
  virtual ~DispersionPolynomial() {}

//-----------------------------------------------------------------------
/// Dispersion offset. This is just coeff(0), but we wrap this for use
/// by DispersionPolynomialOutput
//-----------------------------------------------------------------------

  double dispersion_offset() const { return coeff.value()(0); }

//-----------------------------------------------------------------------
/// Dispersion spacing. This is just coeff(1), but we wrap this for use
/// by DispersionPolynomialOutput
//-----------------------------------------------------------------------

  double dispersion_spacing() const { return coeff.value()(1); }

//-----------------------------------------------------------------------
/// Dispersion offset uncertainty. This is just sqrt(Cov(0,0)), but we
/// wrap this for use  by DispersionPolynomialOutput
//-----------------------------------------------------------------------

  double dispersion_offset_uncertainty() const 
  { 
    if(sv_cov_sub.rows() < 1)
      return 0;
    double t = sv_cov_sub(0,0);
    return (t < 0 ? 0 : sqrt(t)); 
  }

//-----------------------------------------------------------------------
/// Dispersion spacing uncertainty. This is just sqrt(Cov(1,1)), but we
/// wrap this for use  by DispersionPolynomialOutput
//-----------------------------------------------------------------------

  double dispersion_spacing_uncertainty() const 
  { 
    if(sv_cov_sub.rows() < 2)
      return 0;
    double t = sv_cov_sub(1,1);
    return (t < 0 ? 0 : sqrt(t)); 
  }

  virtual boost::shared_ptr<SampleGrid> clone() const;

  virtual std::string sub_state_identifier() const { return "dispersion/" + band_name_; } 

  virtual std::string state_vector_name_i(int i) const;
  virtual SpectralDomain pixel_grid() const;
  virtual SpectralDomain sample_grid() const { return this->pixel_grid(); };
  virtual void print(std::ostream& Os) const;
private:
  void initialize();

  Unit coeff_unit;
  std::string band_name_;
  /// This is an array like 0,1,2,3 ... number_sample. This is used by 
  /// sample_grid.
  blitz::Array<double, 1> variable_values_;
  // Very similar to index_array, but always 1 based.
  blitz::Array<int, 1> spectral_index;
  DispersionPolynomial() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(DispersionPolynomial);
#endif
