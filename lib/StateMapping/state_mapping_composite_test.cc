#include "state_mapping_composite.h"
#include "state_mapping_interpolate.h"
#include "state_mapping_log.h"

#include "unit_test_support.h"
#include "state_mapping_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(state_mapping_composite, StateMappingInterpFixture)

BOOST_AUTO_TEST_CASE(interp_plus_log)
{
    boost::shared_ptr<StateMappingInterpolateLogLog> map_interp =
        boost::make_shared<StateMappingInterpolateLogLog>(press_to, press_from);
    boost::shared_ptr<StateMappingLog> map_log =
        boost::make_shared<StateMappingLog>();

    std::vector<boost::shared_ptr<StateMapping> > mappings;
    mappings.push_back(map_log);
    mappings.push_back(map_interp);

    boost::shared_ptr<StateMappingComposite> map_composite =
        boost::make_shared<StateMappingComposite>(mappings);

    Array<double, 1> expt_mapped(10);
    expt_mapped = 1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0;

    firstIndex i1;
    Array<double, 1> expt_init(5);
    expt_init = log(vals_from_ad.value()(i1));

    ArrayAd<double, 1> init_vals = map_composite->retrieval_state(vals_from_ad);

    // Call with values as they would exist in the state vector which is 
    // the results of the retrieval_state calls
    ArrayAd<double, 1> mapped_vals = map_composite->mapped_state(expt_init);

    BOOST_CHECK_MATRIX_CLOSE_TOL(init_vals.value(), expt_init, 1e-7);
    BOOST_CHECK_MATRIX_CLOSE_TOL(mapped_vals.value(), expt_mapped, 1e-7);
}

BOOST_AUTO_TEST_CASE(serialize)
{
  if(!have_serialize_supported())
    return;

  boost::shared_ptr<StateMappingInterpolateLogLog> map_interp =
    boost::make_shared<StateMappingInterpolateLogLog>(press_to, press_from);
  boost::shared_ptr<StateMappingLog> map_log =
    boost::make_shared<StateMappingLog>();
  
  std::vector<boost::shared_ptr<StateMapping> > mappings;
  mappings.push_back(map_log);
  mappings.push_back(map_interp);
  
  boost::shared_ptr<StateMappingComposite> map_comp_orig =
    boost::make_shared<StateMappingComposite>(mappings);
  
  std::string serial_str = serialize_write_string(map_comp_orig);
  
  boost::shared_ptr<StateMappingComposite> map_comp_restore =
    serialize_read_string<StateMappingComposite>(serial_str);
  
  firstIndex i1;
  Array<double, 1> vals_from_log(5);
  vals_from_log = log(vals_from_ad.value()(i1));
  
  ArrayAd<double, 1> vals_expt = map_comp_orig->mapped_state(vals_from_log);
  ArrayAd<double, 1> vals_restore = map_comp_restore->mapped_state(vals_from_log);
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(vals_restore.value(), vals_expt.value(), 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()
