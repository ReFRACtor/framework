#include "rayleigh_young.h"

#include "pressure.h"
#include "serialized_configuration_fixture.h"
#include "altitude_hydrostatic.h"
#include "atmosphere_standard.h"
#include "unit_test_support.h"
#include "hdf_file.h"
#include "default_constant.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(rayleigh_young, LambertianConfigurationFixture)

BOOST_AUTO_TEST_CASE(cross_section)
{
  boost::shared_ptr<Constant> constant(new DefaultConstant());
  boost::shared_ptr<Pressure> p = config_pressure;
  boost::shared_ptr<Temperature> t = config_temperature;
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit height(416, units::m);
  std::vector<boost::shared_ptr<Altitude> > alt;
  alt.push_back(boost::shared_ptr<Altitude>(new AltitudeHydrostatic(p, t, lat, height)));
  RayleighYoung r(p, alt, constant);
 
  DoubleWithUnit c = r.cross_section(DoubleWithUnit(532, units::nm));
  BOOST_CHECK_CLOSE(c.convert(Unit("m^2")).value, 5.14127e-31, 1e-4);
}

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<Constant> constant(new DefaultConstant());
  boost::shared_ptr<Pressure> p = config_pressure;
  boost::shared_ptr<Temperature> t = config_temperature;
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit height(416, units::m);
  std::vector<boost::shared_ptr<Altitude> > alt;
  alt.push_back(boost::shared_ptr<Altitude>(new AltitudeHydrostatic(p, t, lat, height)));
  RayleighYoung r(p, alt, constant);
 
  IfstreamCs expt_data(test_data_dir() + "expected/rayleigh_young/optical_depth");

  Array<double, 1> od_expect(18);
  expt_data >> od_expect;

  Array<double, 1> od_calc(r.optical_depth_each_layer(12930.30, 0).value());

  BOOST_CHECK_MATRIX_CLOSE_TOL(od_expect, od_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(jacobian)
{
  Rayleigh& r = *(dynamic_cast<const AtmosphereStandard&>(*config_atmosphere).rayleigh_ptr());
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());

  ArrayAd<double, 1> rgrid = r.optical_depth_each_layer(12929.94, 0);
  Array<double, 1> rgrid0(rgrid.shape());
  rgrid0 = rgrid.value();
  Array<double, 2> jac = rgrid.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(rgrid0.shape());
    jacfd = (r.optical_depth_each_layer(12929.94, 0).value() - rgrid0)
      / epsilon(i);
    if(false) {                        // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
        std::cerr << i << ": " << jac(Range::all(), i) << "\n"
                  << jacfd << "\n"
                  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 1e-11);
  }
}

BOOST_AUTO_TEST_SUITE_END()
