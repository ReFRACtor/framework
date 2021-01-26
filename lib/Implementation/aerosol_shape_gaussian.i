// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "aerosol_shape_gaussian.h"
%}
%base_import(aerosol_extinction_imp_base)
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AerosolShapeGaussian)
namespace FullPhysics {

%feature("notabstract") AerosolShapeGaussian;

class AerosolShapeGaussian : public AerosolExtinctionImpBase {
public:
  AerosolShapeGaussian(const boost::shared_ptr<Pressure>& Press,
		       const blitz::Array<double, 1>& Coeffs,
		       const std::string& Aerosol_name,
		       const bool Linear_AOD);
  virtual boost::shared_ptr<AerosolExtinction> clone() const;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  %pickle_serialization();
protected:
  virtual void calc_aerosol_extinction() const;
};
}
