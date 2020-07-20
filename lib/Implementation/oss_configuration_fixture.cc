#include <boost/algorithm/string.hpp>

#include "oss_configuration_fixture.h"
#include "absorber_vmr_level.h"
#include "surface_temperature_direct.h"
#include "temperature_level_offset.h"
#include "pressure_sigma.h"
#include "absorber_absco.h"
#include "gas_absorption.h"
#include "constant.h"
#include "ground_emissivity_piecewise.h"

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
    Array<double, 1> skin_temp(1);
    skin_temp = static_cast<double>(skin_temp_nc);
    ArrayWithUnit<double, 1> skin_temp_with_unit(skin_temp, units::K);
    // OSS always returns skin temp jacobians
    blitz::Array<bool, 1> retrieve_skin_temp(1);
    retrieve_skin_temp = true;
    config_skin_temperature = boost::make_shared<SurfaceTemperatureDirect>(skin_temp_with_unit, retrieve_skin_temp);
    retrieval_skin_temperature_flag = true;

    Array<float, 1> temp_nc = input_data->read_field<float, 1>("/Temperature")(Range::all());
    Array<double, 1> temp(temp_nc.rows());
    for (int i = 0; i < temp_nc.rows(); i++) {
        int level_index = (temp_nc.rows() - 1) - i;
        temp(i) = static_cast<double>(temp_nc(level_index));
    }
    Array<bool, 1> temp_flag(temp_nc.rows());
    temp_flag = true;
    config_temperature = boost::make_shared<TemperatureLevelOffset>(config_pressure, temp, 0.0, true);
    retrieval_temperature_levels.resize(1);
    retrieval_temperature_levels(0) = 0; // config to retrieve only 0-th temperature level

    /* AbsorberVmrs */
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
    retrieval_gas_levels = std::vector<Array<int, 1>>();
    for (int gas_index = 0; gas_index < hdf_gas_names.rows(); gas_index++) {
        std::string full_gas_name = std::string();
        for (auto& gas_name_char : hdf_gas_names(gas_index, Range::all())) {
            full_gas_name += gas_name_char;
        }
        std::string gas_name = boost::trim_right_copy(full_gas_name);
        gas_names.push_back(gas_name);

        Array<double, 1> vmr_curr_gas = Array<double, 1>(cast<double>(vmr_gas(gas_index, Range::all())));
        blitz::Array<bool, 1> vmr_curr_flag = blitz::Array<bool, 1>(vmr_curr_gas.rows());
        Array<int, 1> gas_levels;
        if (find(gas_jacob_names.begin(), gas_jacob_names.end(), gas_name) != gas_jacob_names.end()) {
            gas_levels.resize(vmr_curr_gas.rows());
            for (int gas_level_index = 0; gas_level_index < vmr_curr_gas.rows(); gas_level_index++) {
                gas_levels(gas_level_index) = gas_level_index;
            }
            vmr_curr_flag = true;
        } else {
            vmr_curr_flag = false;
        }
        retrieval_gas_levels.push_back(gas_levels);
        boost::shared_ptr<AbsorberVmrLevel> current_gas = boost::make_shared<AbsorberVmrLevel>(config_pressure,
                vmr_curr_gas, vmr_curr_flag, gas_name);
        vmr.push_back(current_gas);
    }
    config_vmr = vmr;

    /* Ground */
    Array<float, 1> spectral_points_nc = input_data->read_field<float, 1>("/SurfaceGrid")(Range::all());
    Array<double, 1> spectral_points_double = Array<double, 1>(cast<double>(spectral_points_nc(Range::all())));
    ArrayWithUnit<double, 1> spectral_points(spectral_points_double, units::inv_cm);

    Array<float, 1> emissivity_val_nc = input_data->read_field<float, 1>("/Emissivity")(Range::all());
    Array<double, 1> emissivity_val = Array<double, 1>(cast<double>(emissivity_val_nc(Range::all())));

    // OSS always returns emissivity jacobians
    blitz::Array<bool, 1> retrieve_emiss(spectral_points_nc.rows());
    retrieve_emiss = true;
    retrieval_emissivity_flags.resize(1);
    retrieval_emissivity_flags(0) = 5; // retrieve only the 5th surface point emissivity jacobian
    retrieval_reflectivity_flags.resize(1);
    retrieval_reflectivity_flags(0) = 10; // retrieve only the 10th surface point reflectivity jacobian

    config_ground = boost::make_shared<GroundEmissivityPiecewise>(spectral_points, emissivity_val, retrieve_emiss);

    config_obs_zen_ang = DoubleWithUnit(1.45646667, units::deg);
    config_sol_zen_ang = DoubleWithUnit(90.0, units::deg);
    config_lat = DoubleWithUnit(45.0, units::deg);
    config_surf_alt = DoubleWithUnit(0.0000639999998, units::m);
    config_lambertian = true;
}

OssConfigurationFixture::~OssConfigurationFixture() {};
