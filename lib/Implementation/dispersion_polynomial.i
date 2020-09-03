%include "fp_common.i"
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
  %python_attribute(sub_state_identifier, std::string);
  %python_attribute(pixel_grid, SpectralDomain)
  %python_attribute(sample_grid, SpectralDomain);
  virtual boost::shared_ptr<FullPhysics::SampleGrid> clone() const;
  %pickle_serialization();
};

}
