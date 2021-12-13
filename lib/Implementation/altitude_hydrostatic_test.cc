#include "altitude_hydrostatic.h"
#include "lua_configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(altitude_hydrostatic, LuaConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<Pressure> p = config_pressure;
  boost::shared_ptr<Temperature> t = config_temperature;
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit height(416, units::m);
  AltitudeHydrostatic a(p, t, lat, height);
  boost::shared_ptr<Altitude> aclone = a.clone();
  blitz::Array<double, 1> expect(19);
  
  expect =
    62.697223023642123962, 20.577660190341759971, 16.290429242058728221, 13.751577611346435859, 
    11.911469060807121423, 10.463462920915098664, 9.2490807295643264752, 8.1854184180007099059,
    7.2316545873232964681, 6.3651697590650346825, 5.5706413492396471554, 4.8370244868878096156,
    4.1562384410142669822, 3.5211085530775263486, 2.9252184890059860223, 2.3638440396079425376,
    1.8374680402759244746, 1.3420279053836954297, 0.86947271645329127221;

  for(int i = 0; i < expect.rows(); ++i) {
    BOOST_CHECK_CLOSE(a.altitude(p->pressure_grid()(i)).value.value(), expect(i),1e-4);
    BOOST_CHECK_CLOSE(aclone->altitude(p->pressure_grid()(i)).value.value(), expect(i),1e-4);
  }

  expect = 
    9.6561320517985702594, 9.7667372553938367474, 9.7798821359282435139, 9.7876624450940497013,
    9.7933047717039976021, 9.7977476621188479555, 9.8014756913664928817, 9.8047425431084676006,
    9.807673145833788908, 9.8103366689554896141, 9.8127799312168964008, 9.8150366779751543334,
    9.8171316000262880408, 9.8190866304136168452, 9.8209213943675575109, 9.8226503148015709144,
    9.8242718996219515759, 9.8257985771051217228, 9.8272550378505822977;

  for(int i = 0; i < expect.rows(); ++i) {
    BOOST_CHECK_CLOSE(a.gravity(p->pressure_grid()(i)).value.value(), 
                      expect(i),1e-4);
    BOOST_CHECK_CLOSE(aclone->gravity(p->pressure_grid()(i)).value.value(), 
                      expect(i),1e-4);
  }

  // Turn on for debugging
  if (false) {
    std::cerr << setprecision(20);
    std::cerr << "Altitude: " << expect.rows() << std::endl;
    for(int i = 0; i < expect.rows(); ++i) {
      std::cerr << a.altitude(p->pressure_grid()(i)).value.value() << ", ";
    }
    std::cerr << std::endl;
    std::cerr << "Gravity: " << expect.rows() << std::endl;
    for(int i = 0; i < expect.rows(); ++i) {
      std::cerr << a.gravity(p->pressure_grid()(i)).value.value() << ", ";
    }
    std::cerr << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(jacobian)
{
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit height(416, units::m);
  AltitudeHydrostatic a(config_pressure, config_temperature, 
                        lat, height);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());

  ArrayAdWithUnit<double, 1> p(config_pressure->pressure_grid());
  Array<AutoDerivative<double>, 1> alt(p.rows());
  Array<AutoDerivative<double>, 1> grav(p.rows());
  Array<double, 1> alt0(alt.shape()), grav0(grav.shape());
  for(int i = 0; i < p.rows(); ++i) {
    alt(i) = a.altitude(p(i)).value;
    grav(i) = a.gravity(p(i)).value;
    alt0(i) = alt(i).value();
    grav0(i) = grav(i).value();
  }
  Array<double, 2> altjac = FullPhysics::jacobian(alt);
  Array<double, 2> gravjac = FullPhysics::jacobian(grav);
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> altjacfd(alt0.shape()), gravjacfd(grav0.shape());
    for(int j = 0; j < p.rows(); ++j) {
      altjacfd(j) = (a.altitude(p(j)).value.value() - alt0(j)) / epsilon(i);
      gravjacfd(j) = (a.gravity(p(j)).value.value() - grav0(j)) / epsilon(i);
    }
    if(false) {                  // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(altjac(Range::all(), i) - altjacfd));
      std::cerr.precision(20);
      if(diff > 0)
        std::cerr << i << ": " << altjac(Range::all(), i) << "\n"
                  << altjacfd << "\n"
                  << diff << "\n";
      diff = max(abs(gravjac(Range::all(), i) - gravjacfd));
      if(diff > 0)
        std::cerr << i << ": " << gravjac(Range::all(), i) << "\n"
                  << gravjacfd << "\n"
                  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(altjac(Range::all(), i), altjacfd, 1e-5);
    BOOST_CHECK_MATRIX_CLOSE_TOL(gravjac(Range::all(), i), gravjacfd, 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()
