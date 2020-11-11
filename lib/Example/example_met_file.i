// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%include "swig_python_attribute.i"
%{
#include "example_met_file.h"
%}

%base_import(meteorology)
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::ExampleMetFile);

namespace FullPhysics {

class ExampleMetFile : public Meteorology {
public:
  ExampleMetFile(const boost::shared_ptr<HdfFile>& input_file,
		 const std::string& observation_id);
  ExampleMetFile(const std::string& input_filename,
		 const std::string& observation_id);
  %python_attribute(pressure_levels, blitz::Array<double, 1>);
  %python_attribute(specific_humidity, blitz::Array<double, 1>);
  %python_attribute(surface_pressure, double);
  %python_attribute(windspeed_u, double);
  %python_attribute(windspeed_v, double);
  %python_attribute(temperature, blitz::Array<double, 1>);
  %python_attribute(input, boost::shared_ptr<HdfFile>);
  %python_attribute(data_index, int);
  %pickle_serialization()
};
}
