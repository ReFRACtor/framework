#include "ground_with_cloud_handling.h"
#include "ground_fixture.h"
#include "unit_test_support.h"
#include "state_vector.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_with_cloud_handling, GroundFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  GroundWithCloudHandling g(coxmunk, 0.8);
  BOOST_CHECK(g.spurr_brdf_type() == SpurrBrdfType::COXMUNK);
  BOOST_CHECK_CLOSE(g.surface_parameter(13000, 0)(1).value(), 1.331, 1e-8);
  BOOST_CHECK_CLOSE(g.surface_parameter(6500, 1)(1).value(), 1.332, 1e-8);
  BOOST_CHECK_CLOSE(g.surface_parameter(5000, 2)(1).value(), 1.334, 1e-8);
  g.do_cloud(true);
  BOOST_CHECK(g.spurr_brdf_type() == SpurrBrdfType::LAMBERTIAN);
  BOOST_CHECK_CLOSE(g.surface_parameter(13000, 0)(0).value(), g.cloud_albedo(),
		    1e-8);
  BOOST_CHECK_CLOSE(g.surface_parameter(6500, 1)(0).value(), g.cloud_albedo(),
		    1e-8);
  BOOST_CHECK_CLOSE(g.surface_parameter(5000, 2)(0).value(), g.cloud_albedo(),
		    1e-8);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  auto g = boost::make_shared<GroundWithCloudHandling>(coxmunk, 0.8);
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
		    gr->cloud_albedo(),
		    1e-8);
  BOOST_CHECK_CLOSE(gr->surface_parameter(6500, 1)(0).value(),
		    gr->cloud_albedo(),
		    1e-8);
  BOOST_CHECK_CLOSE(gr->surface_parameter(5000, 2)(0).value(),
		    gr->cloud_albedo(),
		    1e-8);
}

BOOST_AUTO_TEST_SUITE_END()



