#include "string_vector_to_char.h"
#include "unit_test_support.h"

using namespace FullPhysics;
BOOST_FIXTURE_TEST_SUITE(string_vector_to_char, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic)
{
    std::vector<std::string> gas_names = std::vector<std::string>();
    gas_names.push_back("H2O");
    gas_names.push_back("O2");
    StringVectorToChar str_vec_to_char(gas_names);

    BOOST_CHECK_EQUAL(str_vec_to_char.substrlen, 3);
    BOOST_CHECK_EQUAL(str_vec_to_char.num_substr, 2);
    BOOST_CHECK_EQUAL(str_vec_to_char.c_str, "H2OO2 ");
}

BOOST_AUTO_TEST_SUITE_END()
