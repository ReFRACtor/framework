#ifndef EXAMPLE_MET_FILE_H
#define EXAMPLE_MET_FILE_H
#include "meteorology.h"
#include "hdf_file.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements and example Meteorological reader that
  reads data from an HDF file with datasets of the same name
  as the values.
*******************************************************************/

class ExampleMetFile : public Meteorology {
public:
    ExampleMetFile(const boost::shared_ptr<HdfFile>& input_file, const std::string& observation_id);
    ExampleMetFile(const std::string& input_filename, const std::string& observation_id);
    ~ExampleMetFile() {}

    // Define how to read various items
    using Meteorology::pressure_levels;
    blitz::Array<double, 1> pressure_levels() const
        { return read_array(group_name + "/pressure_levels"); }

    using Meteorology::specific_humidity;
    blitz::Array<double, 1> specific_humidity() const
        { return read_array(group_name + "/specific_humidity" ); }

    double surface_pressure() const
        { return read_scalar(group_name + "/surface_pressure" ); }

    double windspeed_u() const
        { return read_scalar(group_name + "/windspeed_u" ); }

    double windspeed_v() const
        { return read_scalar(group_name + "/windspeed_v" ); }

    using Meteorology::temperature;
    blitz::Array<double, 1> temperature() const
        { return read_array(group_name + "/temperature" ); }

    void print(std::ostream& Os) const { Os << "ExampleMetFile"; }

private:

    //-----------------------------------------------------------------------
    /// File format specific reader routines
    //-----------------------------------------------------------------------

    double read_scalar(const std::string& Field) const;
    blitz::Array<double, 1> read_array(const std::string& Field) const;

    boost::shared_ptr<HdfFile> input;
    int data_index;

    // Name of the HDF group to read from
    const std::string group_name = "Meteorology";
};
}
#endif
