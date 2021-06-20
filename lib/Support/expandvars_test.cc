#include "unit_test_support.h"
#include "expandvars.h"
#include <cstdlib>

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(expandvars_test, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string h = getenv("HOME");
  std::string s = getenv("abs_top_srcdir");
  BOOST_CHECK_EQUAL(expandvars("${HOME}/hi/${abs_top_srcdir}"),
		    h + "/hi/" + s);
  std::cerr << expandvars("${HOME}/hi/${abs_top_srcdir}") << "\n";
  std::cerr << expandvars("~/hi") << "\n";
  std::cerr << expandvars("$(ls)") << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

