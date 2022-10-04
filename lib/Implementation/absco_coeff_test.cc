#include "absco_coeff.h"
#include "unit_test_support.h"
#include "ifstream_cs.h"
#include "spectral_bound.h"
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absco_coeff, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string NO2_absco_coeff_path = absco_data_dir() + "/coeff/NO2_19840-37879_v0.0_init_new.nc";
  if ( !boost::filesystem::exists(NO2_absco_coeff_path) ) {
    std::cout << "Test Skipped.\n";
    return;
  }
  IfstreamCs expected_data(test_data_dir() + "expected/absco_coeff/basic");
  AbscoCoeff f(NO2_absco_coeff_path);
  // Note scale here is a nonsense value
  double table_scale = 1.2;
  AbscoCoeff fscale(absco_data_dir() + "/coeff/NO2_19840-37879_v0.0_init_new.nc", table_scale);
  BOOST_CHECK_EQUAL(f.number_broadener(), 0);
  if(false) {
    std::cerr << setprecision(20) << std::scientific
	      << "# This is the expected pressure grid.\n"
              << f.pressure_grid() << "\n"
	      << "# This is tsub_expect\n"
	      << f.temperature_grid()(53, Range(0,7)) << "\n"
	      << "# This is readsub_expect\n"
              << f.read<double, 3>(19858.94)(53, Range(0,7), 0) << "\n";
  }
  Array<double, 1> pgrid_expect;
  expected_data >> pgrid_expect;
  BOOST_CHECK_MATRIX_CLOSE(f.pressure_grid(), pgrid_expect);
  // f->temperature_grid is big, so we just select a single pressure
  // and check the temperature grid. No significance to the row picked
  // - I just grabbed a number.
  Array<double, 1> tsub(f.temperature_grid()(53, Range(0,7)));
  Array<double, 1> tsub_expect;
  expected_data >> tsub_expect;
  BOOST_CHECK_MATRIX_CLOSE(tsub, tsub_expect);
  // Same thing with reading the data.
  Array<double, 1> readsub(f.read<double, 3>(19858.94)(53, Range(0,7), 0));
  Array<double, 1> readsub_expect;
  expected_data >> readsub_expect;
  // Numbers are very small, so we have a small tolerance.
  BOOST_CHECK_MATRIX_CLOSE_TOL(readsub, readsub_expect, 1e-35);
  BOOST_CHECK(f.have_data(19858.94));
  BOOST_CHECK(!f.have_data(100));
  ArrayWithUnit<double, 1> pv, tv;
  ArrayWithUnit<double, 2> bv;
  pv.value.resize(3);
  pv.value = 11459.857421875, 12250.0 ,13516.7548828125;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 192.0, 193.2799987792969;
  tv.units = units::K;
  bv.value.resize(1, 3);
  bv.value = 0,0,0;
  bv.units = units::dimensionless;
  Array<double, 1> abs_expect(3);
  abs_expect = 2.0582147030446068e-19, 2.0586911035762673e-19, 2.0590768499973998e-19;
  for(int i = 0; i < 3; ++i) {
    ArrayWithUnit<double, 1> bva;
    bva.value.resize(1);
    bva.value(0) = bv.value(0,i);
    bva.units = bv.units;
    BOOST_CHECK_CLOSE(f.absorption_cross_section
		      (19858.94, pv(i), tv(i), bva).value,
		      abs_expect(i), 1e-4);
  }
  for(int i = 0; i < 3; ++i) {
    ArrayWithUnit<double, 1> bva;
    bva.value.resize(1);
    bva.value(0) = bv.value(0,i);
    bva.units = bv.units;
    BOOST_CHECK_CLOSE(fscale.absorption_cross_section
		      (19858.94, pv(i), tv(i), bva).value,
		      abs_expect(i) * table_scale, 1e-4);
  }
  DoubleWithUnit pvd(pv.value(1), pv.units);
  AutoDerivativeWithUnit<double>
    tvd(AutoDerivative<double>(tv.value(1), 0, 2), tv.units);
  ArrayAdWithUnit<double, 1> bvd;
  bvd.value.resize(1, 2);
  bvd.value(0) = AutoDerivative<double>(bv.value(0,1), 1, 2);
  bvd.units = bv.units;
  AutoDerivative<double> absv = f.absorption_cross_section(19858.94, pvd, tvd, 
							   bvd).value;
  AutoDerivative<double> absvscale = 
    fscale.absorption_cross_section(19858.94, pvd, tvd, 
				    bvd).value;
  BOOST_CHECK_CLOSE(absv.value(), abs_expect(1), 1e-3);
  BOOST_CHECK_CLOSE(absvscale.value(), abs_expect(1) * table_scale, 1e-3);
  double epsilon = 1e-3;
  tvd.value += epsilon;
  double dabs_dt = (f.absorption_cross_section(19858.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  tvd.value -= epsilon;
  bvd.value(0) = bvd.value(0) + epsilon;
  double dabs_db = (f.absorption_cross_section(19858.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  BOOST_CHECK_CLOSE(absv.gradient()(0), dabs_dt, 1e-4);
  BOOST_CHECK_CLOSE(absv.gradient()(1), dabs_db, 1e-4);
  BOOST_CHECK_CLOSE(absvscale.gradient()(0), dabs_dt * table_scale, 1e-4);
  BOOST_CHECK_CLOSE(absvscale.gradient()(1), dabs_db * table_scale, 1e-4);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  std::string NO2_absco_coeff_path = absco_data_dir() + "/coeff/NO2_19840-37879_v0.0_init_new.nc";
  if ( !boost::filesystem::exists(NO2_absco_coeff_path) ) {
    std::cout << "Test Skipped.\n";
    return;
  }
  boost::shared_ptr<AbscoCoeff> a = boost::make_shared<AbscoCoeff>(NO2_absco_coeff_path);
  std::string d = serialize_write_string(a);
  if(false)
    std::cerr << d;
  boost::shared_ptr<AbscoCoeff> ar =
    serialize_read_string<AbscoCoeff>(d);
  BOOST_CHECK_MATRIX_CLOSE(a->pressure_grid(),
			   ar->pressure_grid());
  // This has nans, so we can't check this
  //BOOST_CHECK_MATRIX_CLOSE(a->temperature_grid(),
  //			   ar->temperature_grid());
  ArrayWithUnit<double, 1> pv, tv;
  ArrayWithUnit<double, 2> bv;
  pv.value.resize(3);
  pv.value = 11459.857421875, 12250.0 ,13516.7548828125;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 192.0, 193.2799987792969;
  tv.units = units::K;
  bv.value.resize(1,3);
  bv.value(0,0) = 0;
  bv.value(0,1) = 0.01;
  bv.value(0,2) = 0.05;
  bv.units = units::dimensionless;
  for(int i = 0; i < 3; ++i) {
    ArrayWithUnit<double, 1> bva;
    bva.value.resize(1);
    bva.value(0) = bv.value(0,i);
    bva.units = bv.units;
    BOOST_CHECK_CLOSE(a->absorption_cross_section
		      (19858.94, pv(i), tv(i), bva).value,
		      ar->absorption_cross_section
		      (19858.94, pv(i), tv(i), bva).value, 1e-4);
  }
}
BOOST_AUTO_TEST_SUITE_END()
