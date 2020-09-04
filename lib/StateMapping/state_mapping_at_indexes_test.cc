#include "state_mapping_at_indexes.h"

#include "unit_test_support.h"
#include "state_mapping_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(state_mapping_at_indexes, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    Array<int, 1> retrieval_indexes(5);
    retrieval_indexes = 1, 3, 5, 7, 9;

    ArrayAd<double, 1> initial_state(10, 0);
    initial_state.value() = 10, 11, 12, 13, 14, 15, 16, 17, 18, 19;

    boost::shared_ptr<StateMappingAtIndexes> at_indexes = 
        boost::make_shared<StateMappingAtIndexes>(retrieval_indexes);

    Array<double, 1> subset_expt(5);
    subset_expt = 11, 13, 15, 17, 19;

    ArrayAd<double, 1> subset_calc = at_indexes->retrieval_state(initial_state);

    BOOST_CHECK_MATRIX_CLOSE_TOL(subset_calc.value(), subset_expt, 1e-10);

    subset_calc.value() += 10;

    ArrayAd<double, 1> mapped_calc = at_indexes->mapped_state(subset_calc);

    Array<double, 1> mapped_expt(10);
    mapped_expt = 10, 21, 12, 23, 14, 25, 16, 27, 18, 29;

    BOOST_CHECK_MATRIX_CLOSE_TOL(mapped_calc.value(), mapped_expt, 1e-10);
    

    // Check serialization if supported
    if(!have_serialize_supported())
      return;

    std::string serial_str = serialize_write_string(at_indexes);

    boost::shared_ptr<StateMappingAtIndexes> map_indexes_restore =
        serialize_read_string<StateMappingAtIndexes>(serial_str);
 
    // Full state should have been stored with the retrieval indexes
    ArrayAd<double, 1> mapped_calc_serial = map_indexes_restore->mapped_state(subset_calc);

    BOOST_CHECK_MATRIX_CLOSE_TOL(mapped_calc_serial.value(), mapped_expt, 1e-10);

}

BOOST_AUTO_TEST_CASE(flags)
{
    Array<bool, 1> retrieval_flags(5);
    retrieval_flags = false, true, false, true, false, true, false, true, false, true;

    ArrayAd<double, 1> initial_state(10, 0);
    initial_state.value() = 10, 11, 12, 13, 14, 15, 16, 17, 18, 19;

    boost::shared_ptr<StateMappingAtIndexes> at_indexes = 
        boost::make_shared<StateMappingAtIndexes>(retrieval_flags);

    Array<double, 1> subset_expt(5);
    subset_expt = 11, 13, 15, 17, 19;

    ArrayAd<double, 1> subset_calc = at_indexes->retrieval_state(initial_state);

    BOOST_CHECK_MATRIX_CLOSE_TOL(subset_calc.value(), subset_expt, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
