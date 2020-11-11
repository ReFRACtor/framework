// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "temperature_offset.h"
%}
%base_import(temperature_imp_base)

%fp_shared_ptr(FullPhysics::TemperatureOffset);

namespace FullPhysics {
class TemperatureOffset: public TemperatureImpBase {
public:
  TemperatureOffset(const boost::shared_ptr<Pressure>& Press,
                    double Temp_offset);
  virtual boost::shared_ptr<Temperature> clone() const = 0;
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(temperature_offset, double);
  %python_attribute(temperature_offset_uncertainty, double);
  %pickle_serialization();
protected:
  virtual void calc_temperature_grid() const = 0;
};
}
