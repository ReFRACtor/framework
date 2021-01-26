#include "unit_test_support.h"
#include "level_1b_cache.h"
#include "example_level_1b.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(level_1b_cache, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<HdfFile> hfile =
    boost::make_shared<HdfFile>(test_data_dir() +
				"in/common/l1b_example_data.h5");
  ExampleLevel1b l1b_ex = ExampleLevel1b(hfile, "2014090915251774");
  Level1bCache l1b_cache(l1b_ex);

  BOOST_CHECK_EQUAL(l1b_cache.number_spectrometer(), 3);

  BOOST_CHECK_CLOSE(l1b_cache.latitude(0).value, -29.802030563354492, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.latitude(1).value, -29.802614212036133, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.latitude(2).value, -29.803081512451172, 1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cache.solar_zenith(0).value, 43.490657806396484, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.solar_zenith(1).value, 43.491470336914062, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.solar_zenith(2).value, 43.492301940917969, 1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cache.solar_azimuth(0).value, 319.0213623046875, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.solar_azimuth(1).value, 319.02096557617188, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.solar_azimuth(2).value, 319.02029418945312, 1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cache.altitude(0).value, 0, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.altitude(1).value, 0, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.altitude(2).value, 0, 1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cache.sounding_zenith(0).value, 34.251796722412109,
		    1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.sounding_zenith(1).value, 34.244930267333984,
		    1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.sounding_zenith(2).value, 34.242877960205078,
		    1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cache.sounding_azimuth(0).value, 138.78759765625, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.sounding_azimuth(1).value, 138.78660583496094,
		    1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.sounding_azimuth(2).value, 138.79510498046875,
		    1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(0)(0), 0.5, 1e-8);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(0)(1), 0.49999275803565979,
		    1e-8);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(0)(2), -0.0026892200112342834,
		    1e-6);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(0)(3), 0.0, 1e-8);
  
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(1)(0), 0.5, 1e-8);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(1)(1), 0.4999927282333374,
		    1e-8);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(1)(2), -0.0026979264803230762,
		    1e-6);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(1)(3), 0.0, 1e-8);
  
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(2)(0), 0.5, 1e-8);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(2)(1), 0.49999341368675232,
		    1e-8);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(2)(2), -0.0025685979053378105,
		    1e-6);
  BOOST_CHECK_CLOSE(l1b_cache.stokes_coefficient(2)(3), 0.0, 1e-8);
  
  BOOST_CHECK_CLOSE(l1b_cache.relative_velocity(0).value, 3350.394775390625,
		    1e-4);
  BOOST_CHECK_CLOSE(l1b_cache.time(0).pgs_time(), 684429925.64504862, 1e-4);
  
  BOOST_CHECK_EQUAL(l1b_cache.radiance(0).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_cache.radiance(0).data()(1005), 6.0255898020264018e+19,
		    1e-4);
  BOOST_CHECK_EQUAL(l1b_cache.radiance(1).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_cache.radiance(1).data()(1005), 2.3494947203480289e+19,
		    1e-4);
  BOOST_CHECK_EQUAL(l1b_cache.radiance(2).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_cache.radiance(2).data()(1005), 8.0934638603730944e+18,
		    1e-4);
  
  BOOST_CHECK_EQUAL(l1b_cache.sample_grid(0).data().extent(blitz::firstDim),
		    1016);
  BOOST_CHECK_CLOSE(l1b_cache.sample_grid(0).data()(1005), 0.7724454693899255,
		    1e-4);
  BOOST_CHECK_EQUAL(l1b_cache.sample_grid(1).data().extent(blitz::firstDim),
		    1016);
  BOOST_CHECK_CLOSE(l1b_cache.sample_grid(1).data()(1005), 1.6214942358960505,
		    1e-4);
  BOOST_CHECK_EQUAL(l1b_cache.sample_grid(2).data().extent(blitz::firstDim),
		    1016);
  BOOST_CHECK_CLOSE(l1b_cache.sample_grid(2).data()(1005), 2.082915352703643,
		    1e-4);
}


BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  boost::shared_ptr<HdfFile> hfile =
    boost::make_shared<HdfFile>(test_data_dir() +
				"in/common/l1b_example_data.h5");
  ExampleLevel1b l1b_ex = ExampleLevel1b(hfile, "2014090915251774");
  boost::shared_ptr<Level1bCache> l1b_cache =
    boost::make_shared<Level1bCache>(l1b_ex);
  std::string d = serialize_write_string(l1b_cache);
  if(true)
    std::cerr << d;
  boost::shared_ptr<Level1bCache> l1b_cacher =
    serialize_read_string<Level1bCache>(d);
  BOOST_CHECK_EQUAL(l1b_cacher->number_spectrometer(), 3);

  BOOST_CHECK_CLOSE(l1b_cacher->latitude(0).value, -29.802030563354492, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->latitude(1).value, -29.802614212036133, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->latitude(2).value, -29.803081512451172, 1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cacher->solar_zenith(0).value, 43.490657806396484, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->solar_zenith(1).value, 43.491470336914062, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->solar_zenith(2).value, 43.492301940917969, 1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cacher->solar_azimuth(0).value, 319.0213623046875, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->solar_azimuth(1).value, 319.02096557617188, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->solar_azimuth(2).value, 319.02029418945312, 1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cacher->altitude(0).value, 0, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->altitude(1).value, 0, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->altitude(2).value, 0, 1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cacher->sounding_zenith(0).value, 34.251796722412109,
		    1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->sounding_zenith(1).value, 34.244930267333984,
		    1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->sounding_zenith(2).value, 34.242877960205078,
		    1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cacher->sounding_azimuth(0).value, 138.78759765625, 1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->sounding_azimuth(1).value, 138.78660583496094,
		    1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->sounding_azimuth(2).value, 138.79510498046875,
		    1e-4);
  
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(0)(0), 0.5, 1e-8);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(0)(1), 0.49999275803565979,
		    1e-8);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(0)(2), -0.0026892200112342834,
		    1e-6);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(0)(3), 0.0, 1e-8);
  
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(1)(0), 0.5, 1e-8);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(1)(1), 0.4999927282333374,
		    1e-8);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(1)(2), -0.0026979264803230762,
		    1e-6);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(1)(3), 0.0, 1e-8);
  
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(2)(0), 0.5, 1e-8);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(2)(1), 0.49999341368675232,
		    1e-8);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(2)(2), -0.0025685979053378105,
		    1e-6);
  BOOST_CHECK_CLOSE(l1b_cacher->stokes_coefficient(2)(3), 0.0, 1e-8);
  
  BOOST_CHECK_CLOSE(l1b_cacher->relative_velocity(0).value, 3350.394775390625,
		    1e-4);
  BOOST_CHECK_CLOSE(l1b_cacher->time(0).pgs_time(), 684429925.64504862, 1e-4);
  
  BOOST_CHECK_EQUAL(l1b_cacher->radiance(0).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_cacher->radiance(0).data()(1005), 6.0255898020264018e+19,
		    1e-4);
  BOOST_CHECK_EQUAL(l1b_cacher->radiance(1).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_cacher->radiance(1).data()(1005), 2.3494947203480289e+19,
		    1e-4);
  BOOST_CHECK_EQUAL(l1b_cacher->radiance(2).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_cacher->radiance(2).data()(1005), 8.0934638603730944e+18,
		    1e-4);
  
  BOOST_CHECK_EQUAL(l1b_cacher->sample_grid(0).data().extent(blitz::firstDim),
		    1016);
  BOOST_CHECK_CLOSE(l1b_cacher->sample_grid(0).data()(1005), 0.7724454693899255,
		    1e-4);
  BOOST_CHECK_EQUAL(l1b_cacher->sample_grid(1).data().extent(blitz::firstDim),
		    1016);
  BOOST_CHECK_CLOSE(l1b_cacher->sample_grid(1).data()(1005), 1.6214942358960505,
		    1e-4);
  BOOST_CHECK_EQUAL(l1b_cacher->sample_grid(2).data().extent(blitz::firstDim),
		    1016);
  BOOST_CHECK_CLOSE(l1b_cacher->sample_grid(2).data()(1005), 2.082915352703643,
		    1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
