%include "fp_common.i"
%{
#include "example_level_1b.h"
%}

%base_import(level_1b_sample_coefficient)
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::ExampleLevel1b);

namespace FullPhysics {

%feature("notabstract") ExampleLevel1b;

class ExampleLevel1b: public Level1bSampleCoefficient {
public:
    ExampleLevel1b(const boost::shared_ptr<HdfFile>& input_file, const std::string& observation_id);
    ExampleLevel1b(const std::string& input_filename, const std::string& observation_id);
    int number_spectrometer() const;
    DoubleWithUnit latitude(int i) const;
    DoubleWithUnit longitude(int i) const;
    DoubleWithUnit solar_zenith(int i) const;;
    DoubleWithUnit solar_azimuth(int i) const;
    DoubleWithUnit altitude(int i) const;
    DoubleWithUnit sounding_zenith(int i) const;
    DoubleWithUnit sounding_azimuth(int i) const;
    blitz::Array<double, 1> stokes_coefficient(int i) const;
    DoubleWithUnit relative_velocity(int i) const;
    ArrayWithUnit<double, 1> spectral_coefficient(int i) const;
    Time time(int i) const;
    SpectralRange radiance(int Spec_index) const;
};
}
