// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "fp_common.i"
%{
#include "relative_humidity.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}
%base_import(generic_object)
%import "absorber.i"
%import "temperature.i"
%import "pressure.i"
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::RelativeHumidity);
namespace FullPhysics {
class RelativeHumidity  : public GenericObject {
public:
  RelativeHumidity(const boost::shared_ptr<Absorber>& Abs, 
		   const boost::shared_ptr<Temperature>& Temp,
		   const boost::shared_ptr<Pressure>& Press);
  virtual boost::shared_ptr<RelativeHumidity> clone() const;
  %python_attribute(relative_humidity_grid, ArrayAd<double, 1>);
  %python_attribute(relative_humidity_layer, ArrayAd<double, 1>);
  %python_attribute(specific_humidity_grid, ArrayAd<double, 1>);
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
}
