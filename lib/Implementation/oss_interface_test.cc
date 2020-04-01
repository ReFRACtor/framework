#include "oss_interface.h"
#include "unit_test_support.h"
#include "hdf_file.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(oss_interface, GlobalFixture)

BOOST_AUTO_TEST_CASE(oss_interface)
{
    std::string sel_file = oss_data_dir() + "aura-tes-B1B2-unapod-loc-clear-23V-M12.4-v1.0.train.sel";
    std::string od_file = oss_data_dir() + "aura-tes-B1B2-unapod-loc-clear-23V-M12.4-v1.0.train.lut";
    std::string sol_file = oss_data_dir() + "newkur.dat";
    std::string fix_file = oss_data_dir() + "default.dat";
    std::string ch_sel_file = "NULL";

    std::vector<std::string> gas_names = std::vector<std::string>();
    gas_names.push_back("H2O");
    gas_names.push_back("CO2");
    gas_names.push_back("O3");
    gas_names.push_back("N2O");
    gas_names.push_back("CO");
    gas_names.push_back("CH4");
    gas_names.push_back("O2");
    gas_names.push_back("NH3");
    gas_names.push_back("CCL4");
    gas_names.push_back("F11");
    gas_names.push_back("F12");

    std::vector<std::string> gas_jacobian_names = std::vector<std::string>();
    gas_jacobian_names.push_back("H2O");
    gas_jacobian_names.push_back("NH3");

    int num_vert_lev = 65;
    int num_surf_points = 501;
    float min_ext_cld = 999;

    OssFixedInputs fixed_inputs = OssFixedInputs(gas_names, gas_jacobian_names, sel_file,
              od_file, sol_file, fix_file, ch_sel_file, num_vert_lev, num_surf_points, min_ext_cld);

    OssMasters oss_master = OssMasters(fixed_inputs);
    oss_master.init();
    OssFixedOutputs fixed_outputs = oss_master.fixed_outputs;
    ArrayWithUnit<float, 1> center_wavenumbers = fixed_outputs.center_wavenumber.convert(Unit("Wavenumbers"));

    float expected_start_wavelength = 923.0;
    float expected_wavelength_step = 0.06;
    int expected_number_wavelength = 3951;
    for (int i = 0; i < expected_number_wavelength; i++) {
        BOOST_CHECK_CLOSE(center_wavenumbers.value(i), expected_start_wavelength + (i * expected_wavelength_step), 1e-5);
    }

    const std::string test_fname(oss_run_dir() + "/tape5_nc4.nc");
    boost::shared_ptr<HdfFile> test_input_file(new HdfFile(test_fname));


    ArrayWithUnit<float, 1> pres = ArrayWithUnit<float, 1>(
            test_input_file->read_field<float, 1>("/Pressure")(Range::all()),
            units::mbar);
    ArrayWithUnit<float, 1> temp = ArrayWithUnit<float, 1>(
            test_input_file->read_field<float, 1>("/Temperature")(Range::all()),
            units::K);
    FloatWithUnit skin_temp = FloatWithUnit(310.0, units::K);
    ArrayWithUnit<float, 2> vmr_gas = ArrayWithUnit<float, 2>(
            test_input_file->read_field<float, 2>("/vmrGas")(Range::all()),
            units::dimensionless);
    ArrayWithUnit<float, 1> emis = ArrayWithUnit<float, 1>(
            test_input_file->read_field<float, 1>("/Emissivity")(Range::all()),
            units::dimensionless);
    ArrayWithUnit<float, 1> refl = ArrayWithUnit<float, 1>(
            test_input_file->read_field<float, 1>("/Reflectivity")(Range::all()),
            units::dimensionless);
    FloatWithUnit scale_cld = FloatWithUnit(0.0, units::dimensionless);
    FloatWithUnit pres_cld = FloatWithUnit(0.0, units::mbar);
    int num_cld = 2;
    ArrayWithUnit<float, 1> ext_cld = ArrayWithUnit<float, 1>(Array<float, 1>(num_cld), Unit("km^-1"));
    ext_cld.value = 0;
    ArrayWithUnit<float, 1> surf_grid = ArrayWithUnit<float, 1>(
            test_input_file->read_field<float, 1>("/SurfaceGrid")(Range::all()),
            units::inv_cm);
    ArrayWithUnit<float, 1> cld_grid = ArrayWithUnit<float, 1>(Array<float, 1>(num_cld), units::inv_cm);
    cld_grid.value = 0;
    FloatWithUnit obs_zen_ang = FloatWithUnit(1.45646667, units::deg);
    FloatWithUnit sol_zen_ang = FloatWithUnit(90.0, units::deg);
    FloatWithUnit lat = FloatWithUnit(45.0, units::deg);
    FloatWithUnit surf_alt = FloatWithUnit(0.0000639999998, units::m);
    bool lambertian = true;

    OssModifiedInputs modified_inputs = OssModifiedInputs(pres, temp, skin_temp, vmr_gas, emis,
            refl, scale_cld, pres_cld, ext_cld, surf_grid, cld_grid, obs_zen_ang, sol_zen_ang,
            lat, surf_alt, lambertian);

    OssModifiedOutputs modified_outputs = oss_master.run_fwd_model(modified_inputs);
    std::cout << modified_outputs.y;
    float expected_rad_0 = 0.123218864;
    BOOST_CHECK_CLOSE(modified_outputs.y.value(0), expected_rad_0, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
