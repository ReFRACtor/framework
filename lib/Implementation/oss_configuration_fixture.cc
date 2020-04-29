#include <boost/algorithm/string.hpp>

#include "oss_configuration_fixture.h"
#include "absorber_vmr_level.h"
#include "temperature_level_offset.h"
#include "pressure_sigma.h"
#include "absorber_absco.h"
#include "gas_absorption.h"
#include "constant.h"

using namespace FullPhysics;
using namespace blitz;


OssConfigurationFixture::OssConfigurationFixture(const std::string& input_file)
{
    input_data = boost::make_shared<HdfFile>(oss_run_dir() + input_file);

    /* Pressure */
    Array<float, 1> pressure_nc = input_data->read_field<float, 1>("/Pressure")(Range::all());
    double surface_pressure_nc = static_cast<double>(pressure_nc(0));
    DoubleWithUnit surface_pressure(surface_pressure_nc, units::mbar);
    surface_pressure = surface_pressure.convert(units::Pa);
    Array<double, 1> pressure_top_bottom(pressure_nc.rows());
    for (int i = 0; i < pressure_nc.rows(); i++) {
        int level_index = (pressure_nc.rows() - 1) - i;
        pressure_top_bottom(i) = static_cast<double>(pressure_nc(level_index));
    }
    ArrayWithUnit<double,1> fm_pressure(pressure_top_bottom, units::mbar);
    fm_pressure = fm_pressure.convert(units::Pa);
    config_pressure = boost::make_shared<PressureSigma>(fm_pressure.value, surface_pressure.value, true);

    /* Temperature */
    float skin_temp_nc = input_data->read_field<float>("/SkinTemperature");
    config_skin_temperature = DoubleWithUnit(static_cast<double>(skin_temp_nc), units::K);

    Array<float, 1> temp_nc = input_data->read_field<float, 1>("/Temperature")(Range::all());
    Array<double, 1> temp(temp_nc.rows());
    for (int i = 0; i < temp_nc.rows(); i++) {
        int level_index = (temp_nc.rows() - 1) - i;
        temp(i) = static_cast<double>(temp_nc(level_index));
    }
    Array<bool, 1> temp_flag(temp_nc.rows());
    temp_flag = true;
    /*
     TemperatureFixedLevel::TemperatureFixedLevel(const blitz::Array<bool, 1>& Flag_temp,
     bool Flag_offset, const blitz::Array<double, 1>& Temp,
     double T_offset,
     const boost::shared_ptr<Pressure>& Press,
     const boost::shared_ptr<PressureLevelInput>& Press_level)

    config_temperature = boost::make_shared<TemperatureFixedLevel>(temp_flag, true, temp, surface_temp, config_pressure, pressure_level);
    */

    /*
    TemperatureLevelOffset(const boost::shared_ptr<Pressure>& Press,
                             const blitz::Array<double, 1>& Temp_levels,
                             double Temp_offset,
                             bool Temp_flag);
    */
    config_temperature = boost::make_shared<TemperatureLevelOffset>(config_pressure, temp, 0.0, true);

    /* AbsorberVmrs */
    config_absorber_calc_jacob = std::vector<bool>();
    Array<float, 2> vmr_gas = Array<float, 2>(input_data->read_field<float, 2>("/vmrGas")(Range::all()));
    std::vector<boost::shared_ptr<AbsorberVmr>> vmr = std::vector<boost::shared_ptr<AbsorberVmr>>();

    /* Retrieve list of gases for which we will retrieve jacobians and right trim */
    Array<std::string, 2> hdf_gas_jacob_names = Array<std::string, 2>(input_data->read_field<std::string, 2>("/nameJacobian")(Range::all()));
    std::vector<std::string> gas_jacob_names = std::vector<std::string>();
    for (int gas_jacob_index = 0; gas_jacob_index < hdf_gas_jacob_names.rows(); gas_jacob_index++) {
        std::string full_gas_jacob_name = std::string();
        for (auto& gas_jacob_name_char : hdf_gas_jacob_names(gas_jacob_index, Range::all())) {
            full_gas_jacob_name += gas_jacob_name_char;
        }
        gas_jacob_names.push_back(boost::trim_right_copy(full_gas_jacob_name));
    }

    /* Retrieve list of all gases and right trim, marking those for which we want jacobians. Create AbsorberVmrs */
    Array<std::string, 2> hdf_gas_names = Array<std::string, 2>(input_data->read_field<std::string, 2>("/nameGas")(Range::all()));
    std::vector<std::string> gas_names = std::vector<std::string>();
    for (int gas_index = 0; gas_index < hdf_gas_names.rows(); gas_index++) {
        std::string full_gas_name = std::string();
        for (auto& gas_name_char : hdf_gas_names(gas_index, Range::all())) {
            full_gas_name += gas_name_char;
        }
        std::string gas_name = boost::trim_right_copy(full_gas_name);
        gas_names.push_back(gas_name);

        Array<double, 1> vmr_curr_gas = Array<double, 1>(cast<double>(vmr_gas(gas_index, Range::all())));
        blitz::Array<bool, 1> vmr_curr_flag = blitz::Array<bool, 1>(vmr_curr_gas.rows());
        if (find(gas_jacob_names.begin(), gas_jacob_names.end(), gas_name) != gas_jacob_names.end()) {
            config_absorber_calc_jacob.push_back(true);
            vmr_curr_flag = true;
        } else {
            config_absorber_calc_jacob.push_back(false);
            vmr_curr_flag = false;
        }
        boost::shared_ptr<AbsorberVmrLevel> current_gas = boost::make_shared<AbsorberVmrLevel>(config_pressure,
                vmr_curr_gas, vmr_curr_flag, gas_name);
        if (current_gas) {
            StateVector sv;
            sv.add_observer(*current_gas);
            Array<double, 1> x(sv.observer_claimed_size());
            x = 1e-20;
            sv.update_state(x);
        }
        vmr.push_back(current_gas);
    }
    config_vmr = vmr;
  /*
  AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                     const boost::shared_ptr<Pressure>& pressurev,
                     const boost::shared_ptr<Temperature>& temperaturev,
                     const boost::shared_ptr<Aerosol>& aerosolv,
                     const boost::shared_ptr<RelativeHumidity>& rhv,
                     const boost::shared_ptr<Ground>& groundv,
                     const boost::shared_ptr<SurfaceTemperature>& surface_tempv,
                     const std::vector<boost::shared_ptr<Altitude> >& altv,
                     const boost::shared_ptr<Constant>& C);
   */
}

OssConfigurationFixture::~OssConfigurationFixture() {};
