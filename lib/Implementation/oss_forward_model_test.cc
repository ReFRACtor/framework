#include "oss_forward_model.h"
#include "unit_test_support.h"
#include "atmosphere_fixture.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(oss_forward_constant, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(radiance)
{
  // TODO: Add cmake variable for OSS test data dir
  std::string sel_file = "/export/refractor_data/OSS/aura-tes-B1B2-unapod-loc-clear-23V-M12.4-v1.0.train.sel";
  std::string od_file = "/export/refractor_data/OSS/aura-tes-B1B2-unapod-loc-clear-23V-M12.4-v1.0.train.lut";
  std::string sol_file = "/export/refractor_data/OSS//newkur.dat";
  std::string fix_file = "/export/refractor_data/OSS/default.dat";
  std::string ch_sel_file = "NULL";
  OssForwardModel fm(atm, sel_file, od_file, sol_file, fix_file, ch_sel_file);
  fm.setup_grid();
  Spectrum radiance = fm.radiance(0);
  // TODO: Add checks
}

BOOST_AUTO_TEST_SUITE_END()
