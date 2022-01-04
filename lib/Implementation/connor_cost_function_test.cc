#include <iostream>

#include "unit_test_support.h"
#include "lua_configuration_fixture.h"

#include "connor_cost_function.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(connor_cost_function, LuaConfigurationFixture)

void check_cost_func(const boost::shared_ptr<ConnorCostFunction>& cost_func, const boost::shared_ptr<StateVector>& state_vector, const string& test_data_dir)
{
    Array<double, 1> state_vec(state_vector->state().copy());

    Array<double, 1> residual;
    Array<double, 1> se;
    Array<double, 2> jacobian;

    cost_func->cost_function(state_vec, residual, se, jacobian);

    Array<double, 1> residual_expt;
    Array<double, 1> se_expt;
    Array<double, 2> jacobian_expt;

    // Writes out expected values when updates are needed
    if (false) {
        std::ofstream out(test_data_dir + "expected/connor_cost_function/all_pixels");
        out << std::setprecision(20) << std::scientific;
        out << "# residual" << std::endl << residual << std::endl;
        out << "# se" << std::endl << se << std::endl;
        out << "# jacobian" << std::endl << jacobian << std::endl;
    }

    IfstreamCs fm_cost_expected(test_data_dir + "expected/connor_cost_function/all_pixels");
    fm_cost_expected >> residual_expt >> se_expt >> jacobian_expt;

    BOOST_CHECK_MATRIX_CLOSE(residual_expt, residual);
    BOOST_CHECK_MATRIX_CLOSE(se_expt, se);
    BOOST_CHECK_MATRIX_CLOSE(jacobian_expt, jacobian);
}

BOOST_AUTO_TEST_CASE(all_pixels)
{
    is_long_test();               // Skip unless we are running long tests.
    turn_on_logger();             // Have log output show up.

    boost::shared_ptr<ConnorCostFunction> fm_cost_func(new ConnorCostFunction(config_state_vector, config_forward_model, config_observation));

    check_cost_func(fm_cost_func, config_state_vector, test_data_dir());
}

BOOST_AUTO_TEST_CASE(serialization)
{
  is_long_test();               // Skip unless we are running long tests.
  turn_on_logger();             // Have log output show up.

  if(!have_serialize_supported())
    return;
    
  boost::shared_ptr<ConnorCostFunction> fm_cost_func_orig(new ConnorCostFunction(config_state_vector, config_forward_model, config_observation));

  std::string serial_str = serialize_write_string(fm_cost_func_orig);
  boost::shared_ptr<ConnorCostFunction> fm_cost_func_read = serialize_read_string<ConnorCostFunction>(serial_str);

  check_cost_func(fm_cost_func_read, config_state_vector, test_data_dir());
}

BOOST_AUTO_TEST_SUITE_END()
