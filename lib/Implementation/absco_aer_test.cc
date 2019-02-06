#include "absco_aer.h"
#include "unit_test_support.h"
#include "ifstream_cs.h"
#include "spectral_bound.h"
#include <boost/timer.hpp>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absco_aer, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  IfstreamCs expected_data(test_data_dir() + "expected/absco_aer/basic");
  AbscoAer f(absco_aer_data_dir() + "/CO2_04760-06300_v0.0_init.nc");
  // Note scale here is a nonsense value
  double table_scale = 1.2;
  AbscoAer fscale(absco_aer_data_dir() + "/CO2_04760-06300_v0.0_init.nc",
		  table_scale);
  AbscoAer f2(absco_aer_data_dir() + "/CH4_04760-06300_v0.0_init.nc");
  BOOST_CHECK_EQUAL(f.broadener_name(0), "h2o");
  BOOST_CHECK_EQUAL(f.number_broadener_vmr(0), 2);
  BOOST_CHECK_EQUAL(f.broadener_vmr_grid(0).rows(), 2);
  BOOST_CHECK_EQUAL(f2.number_broadener(), 0);
  if(false) {
    std::cerr << setprecision(20) << std::scientific
	      << "# This is the expected pressure grid.\n"
              << f.pressure_grid() << "\n"
	      << "# This is tsub_expect\n"
	      << f.temperature_grid()(53, Range(0,7)) << "\n"
              << "# This is readsub_expect\n"
              << f.read<double, 3>(4799.9928)(53, Range(0,7), 0) << "\n";
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
  Array<double, 1> readsub(f.read<double, 3>(4799.9928)(53, Range(0,7), 0));
  Array<double, 1> readsub_expect;
  expected_data >> readsub_expect;
  // Numbers are very small, so we have a small tolerance.
  BOOST_CHECK_MATRIX_CLOSE_TOL(readsub, readsub_expect, 1e-30);
  BOOST_CHECK(f.have_data(4799.9928));
  BOOST_CHECK(!f.have_data(100));
  ArrayWithUnit<double, 1> pv, tv;
  ArrayWithUnit<double, 2> bv;
  pv.value.resize(3);
  pv.value = 76000, 66000, 40000;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 190.0, 193.2799987792969;
  tv.units = units::K;
  bv.value.resize(1, 3);
  bv.value = 0,0,0;
  bv.units = units::dimensionless;
  Array<double, 1> abs_expect(3);
  abs_expect = 1.8513911105644267e-24, 1.7466562917995074e-24, 1.114117854087635e-24;
  for(int i = 0; i < 3; ++i) {
    ArrayWithUnit<double, 1> bva;
    bva.value.resize(1);
    bva.value(0) = bv.value(0,i);
    bva.units = bv.units;
    BOOST_CHECK_CLOSE(f.absorption_cross_section
		      (4799.9928, pv(i), tv(i), bva).value,
		      abs_expect(i), 1e-4);
  }
  for(int i = 0; i < 3; ++i) {
    ArrayWithUnit<double, 1> bva;
    bva.value.resize(1);
    bva.value(0) = bv.value(0,i);
    bva.units = bv.units;
    BOOST_CHECK_CLOSE(fscale.absorption_cross_section
		      (4799.9928, pv(i), tv(i), bva).value,
		      abs_expect(i) * table_scale, 1e-4);
  }
  DoubleWithUnit pvd(pv.value(1), pv.units);
  AutoDerivativeWithUnit<double>
    tvd(AutoDerivative<double>(tv.value(1), 0, 2), tv.units);
  ArrayAdWithUnit<double, 1> bvd;
  bvd.value.resize(1, 2);
  bvd.value(0) = AutoDerivative<double>(bv.value(0,1), 1, 2);
  bvd.units = bv.units;
  AutoDerivative<double> absv = f.absorption_cross_section(4799.9928, pvd, tvd, 
							   bvd).value;
  AutoDerivative<double> absvscale = 
    fscale.absorption_cross_section(4799.9928, pvd, tvd, 
				    bvd).value;
  BOOST_CHECK_CLOSE(absv.value(), abs_expect(1), 1e-3);
  BOOST_CHECK_CLOSE(absvscale.value(), abs_expect(1) * table_scale, 1e-3);
  double epsilon = 1e-3;
  tvd.value += epsilon;
  double dabs_dt = (f.absorption_cross_section(4799.9928, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  tvd.value -= epsilon;
  bvd.value(0) = bvd.value(0) + epsilon;
  double dabs_db = (f.absorption_cross_section(4799.9928, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  BOOST_CHECK_CLOSE(absv.gradient()(0), dabs_dt, 1e-3);
  BOOST_CHECK_CLOSE(absv.gradient()(1), dabs_db, 1.0);
  BOOST_CHECK_CLOSE(absvscale.gradient()(0), dabs_dt * table_scale, 1e-3);
  BOOST_CHECK_CLOSE(absvscale.gradient()(1), dabs_db * table_scale, 1.0);
}

BOOST_AUTO_TEST_CASE(interpolation)
{
  IfstreamCs expected_data(test_data_dir() + "expected/absco_aer/basic");
  AbscoAer f(absco_aer_data_dir() + "/CO2_04760-06300_v0.0_init.nc");
  // For the given table, we determined a wavenumber that is not on
  // the grid, and the closest point on the grid.
  double wn_not_on_grid = 6200.0015;
  double wn_closest = 6199.9992;
  DoubleWithUnit press(12250, "Pa");
  DoubleWithUnit temp(190, "K");
  blitz::Array<double, 1> bvalue(1);
  bvalue(0) = 0;
  ArrayWithUnit<double,1> broadener(bvalue, "dimensionless");
  // Default behavior is to throw an error
  BOOST_CHECK_THROW(f.absorption_cross_section(wn_not_on_grid, press, temp,
					       broadener), Exception);
  // Check that nearest neighbor returns closest point.
  f.interpolation_type(AbscoAer::NEAREST_NEIGHBOR_WN);
  double t1 = f.absorption_cross_section(wn_not_on_grid, press, temp,
					 broadener).value;
  double t2 = f.absorption_cross_section(wn_closest, press, temp,
					 broadener).value;
  BOOST_CHECK_CLOSE(t1, t2, 1e-4);
  // Check that interpolation returns something other than closest
  // point.
  f.interpolation_type(AbscoAer::INTERPOLATE_WN);
  t1 = f.absorption_cross_section(wn_not_on_grid, press, temp,
				  broadener).value;
  BOOST_CHECK(fabs(t1-t2) > 1e-26);
}

BOOST_AUTO_TEST_CASE(read_o2)
{
  // O2 has 2 broadners. Want to make sure we can read this.
  AbscoAer f(absco_aer_data_dir() + "//O2_06140-13230_v0.0_init.nc");
  // For the given table, we determined a wavenumber that is not on
  // the grid, and the closest point on the grid.
  double wn = 6200;
  DoubleWithUnit press(12250, "Pa");
  DoubleWithUnit temp(190, "K");
  blitz::Array<double, 1> bvalue(1);
  bvalue(0) = 0;
  ArrayWithUnit<double,1> broadener(bvalue, "dimensionless");
  double t = f.absorption_cross_section(wn, press, temp,
					 broadener).value;
  std::cerr << t << "\n";
}

BOOST_AUTO_TEST_CASE(full_wn_range)
{
  // Read through the full wn range.
  AbscoAer a(absco_aer_data_dir() + "/CO2_04760-06300_v0.0_init.nc");
  DoubleWithUnit press(12250, "Pa");
  DoubleWithUnit temp(190, "K");
  blitz::Array<double, 1> bvalue(1);
  bvalue(0) = 0;
  ArrayWithUnit<double,1> broadener(bvalue, "dimensionless");
  int asize = a.wavenumber_grid().shape()[0];
  int aspace = (int) ceil(asize / 1000.0);
  std::vector<double> res;
  for(int i = 0; i < asize; i += aspace)
    res.push_back(a.absorption_cross_section(a.wavenumber_grid()(i), press,
					     temp, broadener).value);
}

BOOST_AUTO_TEST_SUITE_END()
