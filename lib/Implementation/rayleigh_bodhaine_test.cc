#include "rayleigh_bodhaine.h"

#include "pressure.h"
#include "configuration_fixture.h"
#include "altitude_hydrostatic.h"
#include "atmosphere_standard.h"
#include "unit_test_support.h"
#include "default_constant.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(rayleigh_bodhaine, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(cross_section)
{
  boost::shared_ptr<Constant> constant(new DefaultConstant());
  boost::shared_ptr<Pressure> p = config_pressure;
  boost::shared_ptr<Temperature> t = config_temperature;
  DoubleWithUnit lat(11.099972, units::deg);
  DoubleWithUnit height(0, units::m);
  std::vector<boost::shared_ptr<Altitude> > alt;
  alt.push_back(boost::shared_ptr<Altitude>(new AltitudeHydrostatic(p, t, lat, height)));
  RayleighBodhaine r(p, alt, constant);
 
  // Expected values from MUSES version Bodhaine implementation (Fortran)
  IfstreamCs expt_data(test_data_dir() + "expected/rayleigh_bodhaine/cross_section");
  Array<double, 2> xsect_expect;
  expt_data >> xsect_expect;

  for(int wl_idx = 0; wl_idx < xsect_expect.rows(); wl_idx++) {
      DoubleWithUnit c = r.cross_section(DoubleWithUnit(xsect_expect(wl_idx, 0), units::nm));
      BOOST_CHECK_CLOSE(c.convert(Unit("cm^2")).value, xsect_expect(wl_idx, 1), 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()
