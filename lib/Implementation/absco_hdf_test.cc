#include "absco_hdf.h"
#include "unit_test_support.h"
#include "ifstream_cs.h"
#include "spectral_bound.h"
#include <boost/timer.hpp>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absco_hdf, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  IfstreamCs expected_data(test_data_dir() + "expected/absco_hdf/basic");
  AbscoHdf f(absco_data_dir() + "/o2_v151005_cia_mlawer_v151005r1_narrow.h5");
  // Note scale here is a nonsense value
  double table_scale = 1.2;
  AbscoHdf fscale(absco_data_dir() + "/o2_v151005_cia_mlawer_v151005r1_narrow.h5", table_scale);
  BOOST_CHECK_EQUAL(f.number_broadener(), 1);
  Array<double, 1> pgrid_expect;
  expected_data >> pgrid_expect;
  BOOST_CHECK_MATRIX_CLOSE(f.pressure_grid(), pgrid_expect);
  // f->temperature_grid is big, so we just select a single pressure
  // and check the temperature grid. No significance to the row picked
  // - I just grabbed a number.
  Array<double, 1> tsub(f.temperature_grid()(53, Range::all()));
  Array<double, 1> tsub_expect;
  expected_data >> tsub_expect;
  BOOST_CHECK_MATRIX_CLOSE(tsub, tsub_expect);
  // Same thing with reading the data.
  Array<double, 1> readsub(f.read<double, 3>(12929.94)(53, Range::all(), 0));
  Array<double, 1> readsub_expect;
  expected_data >> readsub_expect;
  // Numbers are very small, so we have a small tolerance.
  BOOST_CHECK_MATRIX_CLOSE_TOL(readsub, readsub_expect, 1e-35);
  BOOST_CHECK(f.have_data(12929.94));
  BOOST_CHECK(!f.have_data(100));
  ArrayWithUnit<double, 1> pv, tv;
  ArrayWithUnit<double, 2> bv;
  pv.value.resize(3);
  pv.value = 11459.857421875, 12250.0 ,13516.7548828125;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 190.0, 193.2799987792969;
  tv.units = units::K;
  bv.value.resize(1, 3);
  bv.value = 0,0,0;
  bv.units = units::dimensionless;
  Array<double, 1> abs_expect(3);
  abs_expect = 1.0755117626991715212e-29, 1.1256089624500879072e-29, 1.2291417782245949615e-29;
  for(int i = 0; i < 3; ++i) {
    ArrayWithUnit<double, 1> bva;
    bva.value.resize(1);
    bva.value(0) = bv.value(0,i);
    bva.units = bv.units;
    BOOST_CHECK_CLOSE(f.absorption_cross_section
		      (12929.94, pv(i), tv(i), bva).value,
		      abs_expect(i), 1e-4);
  }
  for(int i = 0; i < 3; ++i) {
    ArrayWithUnit<double, 1> bva;
    bva.value.resize(1);
    bva.value(0) = bv.value(0,i);
    bva.units = bv.units;
    BOOST_CHECK_CLOSE(fscale.absorption_cross_section
		      (12929.94, pv(i), tv(i), bva).value,
		      abs_expect(i) * table_scale, 1e-4);
  }
  DoubleWithUnit pvd(pv.value(1), pv.units);
  AutoDerivativeWithUnit<double>
    tvd(AutoDerivative<double>(tv.value(1), 0, 2), tv.units);
  ArrayAdWithUnit<double, 1> bvd;
  bvd.value.resize(1, 2);
  bvd.value(0) = AutoDerivative<double>(bv.value(0,1), 1, 2);
  bvd.units = bv.units;
  AutoDerivative<double> absv = f.absorption_cross_section(12929.94, pvd, tvd, 
							   bvd).value;
  AutoDerivative<double> absvscale = 
    fscale.absorption_cross_section(12929.94, pvd, tvd, 
				    bvd).value;
  BOOST_CHECK_CLOSE(absv.value(), abs_expect(1), 1e-3);
  BOOST_CHECK_CLOSE(absvscale.value(), abs_expect(1) * table_scale, 1e-3);
  double epsilon = 1e-3;
  tvd.value += epsilon;
  double dabs_dt = (f.absorption_cross_section(12929.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  tvd.value -= epsilon;
  bvd.value(0) = bvd.value(0) + epsilon;
  double dabs_db = (f.absorption_cross_section(12929.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  BOOST_CHECK_CLOSE(absv.gradient()(0), dabs_dt, 1e-4);
  BOOST_CHECK_CLOSE(absv.gradient()(1), dabs_db, 1e-4);
  BOOST_CHECK_CLOSE(absvscale.gradient()(0), dabs_dt * table_scale, 1e-4);
  BOOST_CHECK_CLOSE(absvscale.gradient()(1), dabs_db * table_scale, 1e-4);
}

BOOST_AUTO_TEST_CASE(scale_specindex)
{
  ArrayWithUnit<double, 2> sbd;
  sbd.units = units::inv_cm;
  sbd.value.resize(3,2);
  sbd.value = 
    12950.0, 13190.0,
     6166.0,  6286.0,
     4810.0,  4897.0;
  SpectralBound sb(sbd);

  std::vector<double> tscale;
  tscale.push_back(1.0);
  tscale.push_back(1.1);
  tscale.push_back(1.2);
  AbscoHdf f(absco_data_dir() + "/co2_devi2015_wco2scale-nist_sco2scale-unity.h5");
  AbscoHdf fscale(absco_data_dir() + "/co2_devi2015_wco2scale-nist_sco2scale-unity.h5", sb, tscale);
  DoubleWithUnit pv(12250, "Pa");
  DoubleWithUnit tv(190, "K");
  blitz::Array<double, 1> bvalue(1);
  bvalue(0) = 0;
  ArrayWithUnit<double,1> broadener(bvalue, "dimensionless");
  BOOST_CHECK_CLOSE
    (f.absorption_cross_section(6200, pv, tv, broadener).value * 1.1,
     fscale.absorption_cross_section(6200, pv, tv, broadener).value, 1e-4);
  BOOST_CHECK_CLOSE
    (f.absorption_cross_section(4880, pv, tv, broadener).value * 1.2,
     fscale.absorption_cross_section(4880, pv, tv, broadener).value, 1e-4);
}

BOOST_AUTO_TEST_CASE(absco_4d)
{
  IfstreamCs expected_data(test_data_dir() + "expected/absco_hdf/4d");
  AbscoHdf f(absco_data_dir() + "/o2_v151005_cia_mlawer_v151005r1_narrow.h5");
  BOOST_CHECK_EQUAL(f.broadener_name(0), "h2o");
  BOOST_CHECK_EQUAL(f.number_broadener_vmr(0), 3);
  Array<double, 1> bgrid_expect;
  expected_data >> bgrid_expect;
  BOOST_CHECK_MATRIX_CLOSE(f.broadener_vmr_grid(0), bgrid_expect);
  Array<double, 1> pgrid_expect;
  expected_data >> pgrid_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(f.pressure_grid(), pgrid_expect, 1e-6);
  // f->temperature_grid is big, so we just select a single pressure
  // and check the temperature grid. No significance to the row picked
  // - I just grabbed a number.
  Array<double, 1> tsub(f.temperature_grid()(53, Range::all()));
  Array<double, 1> tsub_expect;
  expected_data >> tsub_expect;
  BOOST_CHECK_MATRIX_CLOSE(tsub, tsub_expect);
  // Same thing with reading the data.
  Array<double, 2> readsub(f.read<double, 3>(12929.94)
			   (53, Range::all(), Range::all()));
  Array<double, 2> readsub_expect;
  expected_data >> readsub_expect;
  // Numbers are very small, so we have a small tolerance.
  BOOST_CHECK_MATRIX_CLOSE_TOL(readsub, readsub_expect, 1e-35);
  BOOST_CHECK(f.have_data(12929.94));
  BOOST_CHECK(!f.have_data(100));
  ArrayWithUnit<double, 1> pv, tv;
  ArrayWithUnit<double, 2> bv;
  pv.value.resize(3);
  pv.value = 11459.857421875, 12250.0 ,13516.7548828125;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 190.0, 193.2799987792969;
  tv.units = units::K;
  bv.value.resize(1,3);
  bv.value(0,0) = 0;
  bv.value(0,1) = 0.01;
  bv.value(0,2) = 0.05;
  bv.units = units::dimensionless;
  Array<double, 1> abs_expect(3);
  abs_expect = 1.0755117626991715212e-29, 1.1265069651377409808e-29, 1.234189722667483543e-29;
  for(int i = 0; i < 3; ++i) {
    ArrayWithUnit<double, 1> bva;
    bva.value.resize(1);
    bva.value(0) = bv.value(0,i);
    bva.units = bv.units;
    BOOST_CHECK_CLOSE(f.absorption_cross_section
		      (12929.94, pv(i), tv(i), bva).value,
		      abs_expect(i), 1e-4);
  }
  DoubleWithUnit pvd(pv.value(1), pv.units);
  AutoDerivativeWithUnit<double>
    tvd(AutoDerivative<double>(tv.value(1), 0, 2), tv.units);
  ArrayAdWithUnit<double, 1> bvd;
  bvd.value.resize(1, 2);
  bvd.value(0) = AutoDerivative<double>(bv.value(0,1), 1, 2);
  bvd.units = bv.units;
  AutoDerivative<double> absv = f.absorption_cross_section(12929.94, pvd, tvd, 
							   bvd).value;
  BOOST_CHECK_CLOSE(absv.value(), abs_expect(1), 1e-3);
  double epsilon = 1e-3;
  tvd.value += epsilon;
  double dabs_dt = (f.absorption_cross_section(12929.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  tvd.value -= epsilon;
  bvd.value(0) = bvd.value(0) + epsilon;
  double dabs_db = (f.absorption_cross_section(12929.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  BOOST_CHECK(fabs(absv.gradient()(0) - dabs_dt) < 1e-6);
  BOOST_CHECK(fabs(absv.gradient()(1) - dabs_db) < 1e-6);
}

BOOST_AUTO_TEST_CASE(timing)
{
  is_timing_test();
  // Read through all the data, to make sure it doesn't take to long.
  boost::timer tm;
  std::cerr << "Starting read of all data\n";
  int j = 0;
  AbscoHdf f(absco_data_dir() + "/o2_v151005_cia_mlawer_v151005r1_narrow.h5");
  for(double i = 12929.94; i < 13210.15; i += 0.01, j++) {
    if(j % 1000 == 0)
      std::cerr << "Reading " << j << "\n"
		<< "Total time: " << tm.elapsed() << "\n";
    f.read<float, 3>(i);
  }
  std::cerr << "Done\n"
	    << "Total time: " << tm.elapsed() << "\n";
}

BOOST_AUTO_TEST_CASE(cache_boundary)
{
  AbscoHdf f(absco_data_dir() + "/o2_v151005_cia_mlawer_v151005r1_narrow.h5");
  // This wn corresponds to the value at line 5000, which is the next
  // cache line.
  Array<double,3> data(f.read<double, 3>(12970.0));
  // We determined the expected results by direct inspection of the
  // HDF file.
  BOOST_CHECK_CLOSE(data(0,0, 0), 2.500249851905119712e-33, 1e-8);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  
  boost::shared_ptr<AbscoHdf> a =
    boost::make_shared<AbscoHdf>(absco_data_dir() + "/o2_v151005_cia_mlawer_v151005r1_narrow.h5");
  std::string d = serialize_write_string(a);
  if(false)
    std::cerr << d;
  boost::shared_ptr<AbscoHdf> ar =
    serialize_read_string<AbscoHdf>(d);
  BOOST_CHECK_EQUAL(a->broadener_name(0), ar->broadener_name(0)); 
  BOOST_CHECK_EQUAL(a->number_broadener_vmr(0),ar->number_broadener_vmr(0));
  BOOST_CHECK_MATRIX_CLOSE(a->broadener_vmr_grid(0),
			   ar->broadener_vmr_grid(0));
  BOOST_CHECK_MATRIX_CLOSE(a->pressure_grid(),
			   ar->pressure_grid());
  BOOST_CHECK_MATRIX_CLOSE(a->temperature_grid(),
			   ar->temperature_grid());
  ArrayWithUnit<double, 1> pv, tv;
  ArrayWithUnit<double, 2> bv;
  pv.value.resize(3);
  pv.value = 11459.857421875, 12250.0 ,13516.7548828125;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 190.0, 193.2799987792969;
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
		      (12929.94, pv(i), tv(i), bva).value,
		      ar->absorption_cross_section
		      (12929.94, pv(i), tv(i), bva).value, 1e-4);
  }
}
BOOST_AUTO_TEST_SUITE_END()
