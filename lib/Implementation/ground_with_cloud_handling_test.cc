#include "ground_with_cloud_handling.h"
#include "ground_fixture.h"
#include "unit_test_support.h"
#include "state_vector.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_with_cloud_handling, GroundFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  blitz::Array<double, 2> spec_coeff(3,1);
  spec_coeff(0,0) = 0.8;
  spec_coeff(1,0) = 0.8;
  spec_coeff(2,0) = 0.8;
  blitz::Array<double, 1> ref_pointv(3);
  ref_pointv(0) = 1000;
  ref_pointv(1) = 1000;
  ref_pointv(2) = 1000;
  ArrayWithUnit<double, 1> ref_point(ref_pointv, "nm");
  std::vector<std::string> desc_band_name;
  desc_band_name.push_back("dummy");
  desc_band_name.push_back("dummy");
  desc_band_name.push_back("dummy");
  auto ground_cloud = boost::make_shared<GroundLambertian>(spec_coeff,
					   ref_point, desc_band_name);
  GroundWithCloudHandling g(coxmunk, ground_cloud);
  BOOST_CHECK(g.spurr_brdf_type() == SpurrBrdfType::COXMUNK);
  BOOST_CHECK_CLOSE(g.surface_parameter(13000, 0)(1).value(), 1.331, 1e-8);
  BOOST_CHECK_CLOSE(g.surface_parameter(6500, 1)(1).value(), 1.332, 1e-8);
  BOOST_CHECK_CLOSE(g.surface_parameter(5000, 2)(1).value(), 1.334, 1e-8);
  g.do_cloud(true);
  BOOST_CHECK(g.spurr_brdf_type() == SpurrBrdfType::LAMBERTIAN);
  BOOST_CHECK_CLOSE(g.surface_parameter(13000, 0)(0).value(), 0.8,
		    1e-8);
  BOOST_CHECK_CLOSE(g.surface_parameter(6500, 1)(0).value(), 0.8,
		    1e-8);
  BOOST_CHECK_CLOSE(g.surface_parameter(5000, 2)(0).value(), 0.8,
		    1e-8);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  blitz::Array<double, 2> spec_coeff(3,1);
  spec_coeff(0,0) = 0.8;
  spec_coeff(1,0) = 0.8;
  spec_coeff(2,0) = 0.8;
  blitz::Array<double, 1> ref_pointv(3);
  ref_pointv(0) = 1000;
  ref_pointv(1) = 1000;
  ref_pointv(2) = 1000;
  ArrayWithUnit<double, 1> ref_point(ref_pointv, "nm");
  std::vector<std::string> desc_band_name;
  desc_band_name.push_back("dummy");
  desc_band_name.push_back("dummy");
  desc_band_name.push_back("dummy");
  auto ground_cloud = boost::make_shared<GroundLambertian>(spec_coeff,
					   ref_point, desc_band_name);
  auto g = boost::make_shared<GroundWithCloudHandling>(coxmunk, ground_cloud);
  std::string d = serialize_write_string(g);
  if(false)
    std::cerr << d;
  auto gr = serialize_read_string<GroundWithCloudHandling>(d);
  BOOST_CHECK(gr->spurr_brdf_type() == SpurrBrdfType::COXMUNK);
  BOOST_CHECK_CLOSE(gr->surface_parameter(13000, 0)(1).value(), 1.331, 1e-8);
  BOOST_CHECK_CLOSE(gr->surface_parameter(6500, 1)(1).value(), 1.332, 1e-8);
  BOOST_CHECK_CLOSE(gr->surface_parameter(5000, 2)(1).value(), 1.334, 1e-8);
  gr->do_cloud(true);
  BOOST_CHECK(gr->spurr_brdf_type() == SpurrBrdfType::LAMBERTIAN);
  BOOST_CHECK_CLOSE(gr->surface_parameter(13000, 0)(0).value(),
		    0.8,
		    1e-8);
  BOOST_CHECK_CLOSE(gr->surface_parameter(6500, 1)(0).value(),
		    0.8,
		    1e-8);
  BOOST_CHECK_CLOSE(gr->surface_parameter(5000, 2)(0).value(),
		    0.8,
		    1e-8);
}

BOOST_AUTO_TEST_SUITE_END()



