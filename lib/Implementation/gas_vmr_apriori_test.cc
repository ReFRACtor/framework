#include "gas_vmr_apriori.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"
#include "pressure_sigma.h"
#include "temperature_met.h"
#include "altitude_hydrostatic.h"
#include "example_met_file.h"
#include "example_level_1b.h"
#
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(gas_vmr_apriori, GlobalFixture)

BOOST_AUTO_TEST_CASE(apriori_calc)
{
    std::vector<std::string> test_sounding_id = { "2014090915251774", "2014120112331638", "2015020202201332", "2015020301255031", "2015020400304333" };
    IfstreamCs expt_vmr_input(test_data_dir() + "expected/gas_vmr_apriori/apriori_vmr");
    IfstreamCs expt_trop_press_input(test_data_dir() + "expected/gas_vmr_apriori/tropopause_pressure");

    for (auto curr_sid : test_sounding_id) {
        // Match setup used in TcconApriori test
        boost::shared_ptr<HdfFile> l1b_file(new HdfFile(test_data_dir() + "in/common/l1b_example_data.h5"));
        boost::shared_ptr<HdfFile> met_file(new HdfFile(test_data_dir() + "in/common/met_example_data.h5"));

        boost::shared_ptr<Level1b> l1b(new ExampleLevel1b(l1b_file, curr_sid));
        boost::shared_ptr<Meteorology> met(new ExampleMetFile(met_file, curr_sid));

        HdfFile hdf_static_input(test_data_dir() + "in/gas_vmr_apriori/reference_atmosphere.h5");

        // Create pressure/temp and hence altitude just like done in the production software
        blitz::Array<double, 1> sigma_a = hdf_static_input.read_field<double, 1>("Pressure/Pressure_sigma_a");
        blitz::Array<double, 1> sigma_b = hdf_static_input.read_field<double, 1>("Pressure/Pressure_sigma_b");

        double surface_pressure = met->surface_pressure();
        boost::shared_ptr<PressureSigma> press
          (new PressureSigma(sigma_a, sigma_b, surface_pressure, false));

        boost::shared_ptr<TemperatureMet> temp
          (new TemperatureMet(met, press, 0, false));

        boost::shared_ptr<Altitude> alt(new AltitudeHydrostatic(press, temp, l1b->latitude(0), l1b->altitude(0)));

        GasVmrApriori gas_vmr(met, l1b, alt, hdf_static_input, "/Reference_Atmosphere/", "CO2");

        Array<double, 1> expt_vmr;
        expt_vmr_input >> expt_vmr;

        Array<double, 1> calc_ap_vmr = gas_vmr.apriori_vmr();
        BOOST_CHECK_MATRIX_CLOSE_TOL(expt_vmr, calc_ap_vmr, 1e-7);

        double expt_trop_pres;
        expt_trop_press_input >> expt_trop_pres;

        BOOST_CHECK_CLOSE(expt_trop_pres, gas_vmr.tropopause_pressure(), 1e-7);
    }
}

/*  From Brendan's presentation on Reichler method
    Expected tropopause pressure:
    2014090915251774 - 114.6 hPa
    2014120112331638 - 126.6 hPa
    2015020202201332 - 96.0 hPa
    2015020301255031 - 95.3 hPa
    2015020400304333 - 316.7 hPa
*/

BOOST_AUTO_TEST_SUITE_END()
