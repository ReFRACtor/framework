// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "ground_coxmunk.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground_imp_base)

%fp_shared_ptr(FullPhysics::GroundCoxmunk);
namespace FullPhysics {
class GroundCoxmunk: public GroundImpBase {
public:
  GroundCoxmunk(const double Windspeed,
                const blitz::Array<double, 1>& Refr_index);
  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
  virtual const AutoDerivative<double> windspeed() const;
  virtual const double refractive_index(const int Spec_idx) const;
  virtual boost::shared_ptr<Ground> clone() const;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;
  virtual void update_sub_state_hook();
  %pickle_serialization();
};
}
