// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "stokes_coefficient_constant.h"
%}

%base_import(stokes_coefficient)

%fp_shared_ptr(FullPhysics::StokesCoefficientConstant);

namespace FullPhysics {
class StokesCoefficientConstant : public StokesCoefficient {
public:
  StokesCoefficientConstant(const blitz::Array<double, 2>& Stokes_coeff);
  virtual boost::shared_ptr<StokesCoefficient> clone() const;
  %python_attribute(stokes_coefficient, ArrayAd<double, 2>)
  void set_stokes_coefficient(const blitz::Array<double, 2> Stokes_coeff);
  %pickle_serialization();
};
}
