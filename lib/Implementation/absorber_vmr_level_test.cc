#include "absorber_vmr_level.h"
#include "met_data_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_vmr_level, MetDataFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<double, 1> vmr(19);
  vmr = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19;
  AbsorberVmrLevel avmr(pressure, vmr, "CO2");
  AbsorberVmrLevel avmr2(pressure_reverse, vmr.reverse(blitz::firstDim), "CO2");
  for(int i = 0; i < pressure->pressure_grid().rows(); ++i) {
    BOOST_CHECK_CLOSE(avmr.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), vmr(i), 1e-8);
    BOOST_CHECK_CLOSE(avmr2.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), vmr(i), 1e-8);
  }
  StateVector sv;
  StateVector sv2;
  sv.add_observer(avmr);
  sv2.add_observer(avmr2);
  sv.update_state(vmr);
  sv2.update_state(vmr.reverse(blitz::firstDim));
  for(int i = 0; i < pressure->pressure_grid().rows(); ++i) {
    BOOST_CHECK_CLOSE(avmr.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), vmr(i), 1e-8);
    BOOST_CHECK_CLOSE(avmr2.volume_mixing_ratio(pressure->pressure_grid()(i).value).value(), vmr(i), 1e-8);
  }
  std::cerr << avmr << "\n"
	    << avmr2 << "\n"
	    << sv << "\n"
	    << sv2 << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
