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
  AbscoAer f(absco_aer_data_dir() + "/CO2_00600-01200_v0.0_init.chunk.nc");
  // Note scale here is a nonsense value
  double table_scale = 1.2;
  AbscoAer fscale(absco_aer_data_dir() + "/CO2_00600-01200_v0.0_init.chunk.nc",
		  table_scale);
  AbscoAer f2(absco_aer_data_dir() + "/CH4_00600-01200_v0.0_init.nc");
  BOOST_CHECK_EQUAL(f.broadener_name(), "h2o");
  BOOST_CHECK_EQUAL(f.number_broadener_vmr(), 2);
  BOOST_CHECK_EQUAL(f.broadener_vmr_grid().rows(), 2);
  BOOST_CHECK_EQUAL(f2.broadener_name(), "");
  BOOST_CHECK_EQUAL(f2.number_broadener_vmr(), 0);
  BOOST_CHECK_EQUAL(f2.broadener_vmr_grid().rows(), 0);
  if(false) {
    std::cerr << setprecision(20) << std::scientific
	      << "# This is the expected pressure grid.\n"
              << f.pressure_grid() << "\n"
	      << "# This is tsub_expect\n"
	      << f.temperature_grid()(53, Range(0,7)) << "\n"
              << "# This is readsub_expect\n"
              << f.read<double>(1099.9992)(Range(0,7), 53, 0) << "\n";
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
  Array<double, 1> readsub(f.read<double>(1099.9992)(Range(0,7), 53, 0));
  Array<double, 1> readsub_expect;
  expected_data >> readsub_expect;
  // Numbers are very small, so we have a small tolerance.
  BOOST_CHECK_MATRIX_CLOSE_TOL(readsub, readsub_expect, 1e-30);
  BOOST_CHECK(f.have_data(1099.9992));
  BOOST_CHECK(!f.have_data(100));
  ArrayWithUnit<double, 1> pv, tv, bv;
  pv.value.resize(3);
  pv.value = 76000, 66000, 40000;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 190.0, 193.2799987792969;
  tv.units = units::K;
  bv.value.resize(3);
  bv.value = 0,0,0;
  bv.units = units::dimensionless;
  Array<double, 1> abs_expect(3);
  abs_expect = 2.6294165614745598e-26, 2.2473571402739345e-26,
    1.3445900408125961e-26;
  for(int i = 0; i < 3; ++i)
    BOOST_CHECK_CLOSE(f.absorption_cross_section
		      (1099.9992, pv(i), tv(i), bv(i)).value,
		      abs_expect(i), 1e-4);
  for(int i = 0; i < 3; ++i)
    BOOST_CHECK_CLOSE(fscale.absorption_cross_section
		      (1099.9992, pv(i), tv(i), bv(i)).value,
		      abs_expect(i) * table_scale, 1e-4);
  DoubleWithUnit pvd(pv.value(1), pv.units);
  AutoDerivativeWithUnit<double>
    tvd(AutoDerivative<double>(tv.value(1), 0, 2), tv.units);
  AutoDerivativeWithUnit<double>
    bvd(AutoDerivative<double>(bv.value(1), 1, 2), bv.units);
  AutoDerivative<double> absv = f.absorption_cross_section(1099.9992, pvd, tvd, 
							   bvd).value;
  AutoDerivative<double> absvscale = 
    fscale.absorption_cross_section(1099.9992, pvd, tvd, 
				    bvd).value;
  BOOST_CHECK_CLOSE(absv.value(), abs_expect(1), 1e-3);
  BOOST_CHECK_CLOSE(absvscale.value(), abs_expect(1) * table_scale, 1e-3);
  double epsilon = 1e-3;
  tvd.value += epsilon;
  double dabs_dt = (f.absorption_cross_section(1099.9992, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  tvd.value -= epsilon;
  bvd.value += epsilon;
  double dabs_db = (f.absorption_cross_section(1099.9992, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  BOOST_CHECK_CLOSE(absv.gradient()(0), dabs_dt, 1e-3);
  BOOST_CHECK_CLOSE(absv.gradient()(1), dabs_db, 1.0);
  BOOST_CHECK_CLOSE(absvscale.gradient()(0), dabs_dt * table_scale, 1e-3);
  BOOST_CHECK_CLOSE(absvscale.gradient()(1), dabs_db * table_scale, 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
