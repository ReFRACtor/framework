%include "fp_common.i"
%{
#include "gas_vmr_apriori.h"
#include "temperature.h"
#include "altitude.h"
%}

%base_import(generic_object)
%import "meteorology.i"
%import "level_1b.i"
%import "pressure.i"
%import "temperature.i"
%import "altitude.i"
%import "hdf_file.i"
%import "reference_vmr_apriori.i"

%fp_shared_ptr(FullPhysics::GasVmrApriori);

namespace FullPhysics {

class GasVmrApriori : public GenericObject {
public:
    GasVmrApriori(const boost::shared_ptr<Meteorology>& met_file,
                  const boost::shared_ptr<Level1b>& l1b_file,
                  const boost::shared_ptr<Altitude>& alt,
                  const HdfFile& hdf_static_input,
                  const std::string& hdf_group,
                  const std::string& gas_name,
                  const int temp_avg_window = 11);

    GasVmrApriori(const blitz::Array<double, 1>& pressure_levels,
                  const blitz::Array<double, 1>& temperature_levels,
                  const double& obs_latitude,
                  const Time& obs_time,
                  const boost::shared_ptr<Altitude>& altitude,
                  const HdfFile& reference_file,
                  const std::string& hdf_group,
                  const std::string& gas_name,
                  const int temp_avg_window = 11);

    const blitz::Array<double, 1> apriori_vmr() const;
    const blitz::Array<double, 1> apriori_vmr(const Pressure& pressure) const;

    %python_attribute(reference, boost::shared_ptr<ReferenceVmrApriori>)
    %python_attribute(tropopause_altitude, DoubleWithUnit)
    %python_attribute(tropopause_pressure, double)

    std::string print_to_string() const;
};
}
