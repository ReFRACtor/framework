// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "pressure_level_input.h"
%}

%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::PressureLevelInput);

namespace FullPhysics {
class PressureLevelInput {
public:
  PressureLevelInput(const blitz::Array<double, 1>& Press_level);
  PressureLevelInput(const HdfFile& Hdf_file, 
		     const std::string& Hdf_group = "Pressure");
  %python_attribute(pressure_level, blitz::Array<double, 1>)
  std::string print_to_string() const;
};
}
