// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "example_level_1b.h"
%}

%base_import(level_1b_sample_coefficient)
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::ExampleLevel1b);

namespace FullPhysics {

class ExampleLevel1b: public Level1bSampleCoefficient {
public:
  ExampleLevel1b(const boost::shared_ptr<HdfFile>& input_file,
		 const std::string& observation_id);
  ExampleLevel1b(const std::string& input_filename,
		 const std::string& observation_id);
  %python_attribute(number_spectrometer, virtual int);
  virtual DoubleWithUnit latitude(int i) const;
  virtual DoubleWithUnit longitude(int i) const;
  virtual DoubleWithUnit solar_zenith(int i) const;;
  virtual DoubleWithUnit solar_azimuth(int i) const;
  virtual DoubleWithUnit altitude(int i) const;
  virtual DoubleWithUnit sounding_zenith(int i) const;
  virtual DoubleWithUnit sounding_azimuth(int i) const;
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const;
  virtual DoubleWithUnit relative_velocity(int i) const;
  virtual int number_sample(int channel_index) const;
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int channel_index)
    const;
  virtual blitz::Array<double, 1> spectral_variable(int channel_index)
    const;
  Time time(int i) const;
  SpectralRange radiance(int Spec_index) const;
  %python_attribute(input, boost::shared_ptr<HdfFile>);
  %python_attribute(data_index, int);
  %pickle_serialization()
};
}
