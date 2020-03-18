#include "oss_interface.h"
#include "unit_test_support.h"
#include "oss_configuration_fixture.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(oss_interface, GlobalFixture)

BOOST_AUTO_TEST_CASE(oss_init)
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
    int num_molecules = gas_names.size();

    std::vector<std::string> gas_jacobian_names = std::vector<std::string>();
    gas_jacobian_names.push_back("H2O");
    gas_jacobian_names.push_back("NH3");

    int num_vert_lev = 65;
    int num_surf_points = 501;
    float min_ext_cld = 999;

    OssFixedInputs fixed_inputs = OssFixedInputs(num_molecules, gas_names, gas_jacobian_names, sel_file,
			  od_file, sol_file, fix_file, ch_sel_file, num_vert_lev, num_surf_points, min_ext_cld);

	OssMasters oss_master = OssMasters(fixed_inputs);
	oss_master.init();
	OssFixedOutputs fixed_outputs = oss_master.FixedOutputs();
}

BOOST_AUTO_TEST_SUITE_END()
