#include <boost/make_shared.hpp>

#include "unit_test_support.h"
#include "hdf_file.h"
#include "fp_time.h"
#include "example_level_1b_info.h"


using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(example_level_1b_info, GlobalFixture)


BOOST_AUTO_TEST_CASE(level1b_list)
{
    float expected_lats[5] = { -29.802, -47.2022, 29.7131, 33.6253, 35.1393 };
    int spec_index = 0;
    ExampleLevel1bInfo ex_lev1_info(test_data_dir() + "in/common/l1b_example_data.h5");
    std::vector<boost::shared_ptr<Level1b>> level1b_series = ex_lev1_info.level1b_list();
    for (int i = 0; i < (int) level1b_series.size(); i++) {
        BOOST_CHECK_CLOSE(level1b_series[i].get()->latitude(spec_index).value, expected_lats[i], 1e-3);
    }
}

BOOST_AUTO_TEST_CASE(obs_distance)
{
    DoubleWithUnit expected_dist(4249.6243199759228, "km");
    std::string input_filename(test_data_dir() + "in/common/l1b_example_data.h5");
    boost::shared_ptr<HdfFile> h_file = boost::make_shared<HdfFile>(input_filename);
    ExampleLevel1bInfo example_level_1b_info(input_filename);
    boost::shared_ptr<Level1b> example_level_1b_ptr = boost::make_shared<ExampleLevel1b>(h_file, "2014090915251774");
    DoubleWithUnit dist = example_level_1b_info.obs_distance(std::move(example_level_1b_ptr), 0.0, 0.0).convert(expected_dist.units);
    BOOST_CHECK_CLOSE(dist.value, expected_dist.value, 1e-8);
    h_file->close();
}

BOOST_AUTO_TEST_CASE(close_obs_distance)
{
    std::string input_filename(test_data_dir() + "in/common/l1b_example_data.h5");
    ExampleLevel1bInfo example_level_1b_info(input_filename);
    int spec_index = 0;
    DoubleWithUnit threshold(4250, "km");
    std::vector<boost::shared_ptr<Level1b>> found_close_obs = example_level_1b_info.close_obs(0.0, 0.0, threshold, spec_index);


    std::vector<Time> expected_times = std::vector<Time>();
    expected_times.push_back(Time::time_unix(1410276325.6450486));

    for (int i = 0; i < (int) found_close_obs.size(); i++) {
        BOOST_CHECK_EQUAL(found_close_obs[i].get()->time(spec_index), expected_times[i]);
    }
}

BOOST_AUTO_TEST_CASE(closest_obs_n_distance)
{
    std::string input_filename(test_data_dir() + "in/common/l1b_example_data.h5");
    ExampleLevel1bInfo example_level_1b_info(input_filename);
    int spec_index = 0;

    std::vector<boost::shared_ptr<Level1b>> closest_5 = example_level_1b_info.closest_obs_n(0.0, 0.0, 5, spec_index);


    std::vector<Time> expected_times = std::vector<Time>();
    expected_times.push_back(Time::time_unix(1410276325.6450486));
    expected_times.push_back(Time::time_unix(1417437204.4190402));
    expected_times.push_back(Time::time_unix(1423009851.2847152));
    expected_times.push_back(Time::time_unix(1422843621.2577152));
    expected_times.push_back(Time::time_unix(1422926758.23139));

    for (int i = 0; i < (int) closest_5.size(); i++) {
        BOOST_CHECK_EQUAL(closest_5[i].get()->time(spec_index), expected_times[i]);
    }
}

BOOST_AUTO_TEST_CASE(closest_obs_distance)
{
    std::string input_filename(test_data_dir() + "in/common/l1b_example_data.h5");
    ExampleLevel1bInfo example_level_1b_info(input_filename);
    int spec_index = 0;

    boost::shared_ptr<Level1b> closest = example_level_1b_info.closest_obs(0.0, 0.0, spec_index);

    Time expected_time = Time::time_unix(1410276325.6450486);
    BOOST_CHECK_EQUAL(closest.get()->time(spec_index), expected_time);
}

BOOST_AUTO_TEST_CASE(close_obs_time)
{
    std::string input_filename(test_data_dir() + "in/common/l1b_example_data.h5");
    ExampleLevel1bInfo example_level_1b_info(input_filename);
    int spec_index = 0;
    Time threshold = Time::time_unix(10);

    std::vector<boost::shared_ptr<Level1b>> found_close_obs = example_level_1b_info.close_obs(Time::parse_time("2014-09-09T15:25:34.645049Z"), threshold, spec_index);

    std::vector<Time> expected_times = std::vector<Time>();
    expected_times.push_back(Time::time_unix(1410276325.6450486));

    for (int i = 0; i < (int) found_close_obs.size(); i++) {
        BOOST_CHECK_EQUAL(found_close_obs[i].get()->time(spec_index), expected_times[i]);
    }
}

BOOST_AUTO_TEST_CASE(closest_obs_n_time)
{
    std::string input_filename(test_data_dir() + "in/common/l1b_example_data.h5");
    ExampleLevel1bInfo example_level_1b_info(input_filename);
    int spec_index = 0;

    std::vector<boost::shared_ptr<Level1b>> closest_5 = example_level_1b_info.closest_obs_n(Time::parse_time("2014-09-09T15:30:25.645049Z"), 5, spec_index);

    std::vector<Time> expected_times = std::vector<Time>();

    expected_times.push_back(Time::time_unix(1410276325.6450486));
    expected_times.push_back(Time::time_unix(1417437204.4190402));
    expected_times.push_back(Time::time_unix(1422843621.2577152));
    expected_times.push_back(Time::time_unix(1422926758.23139));
    expected_times.push_back(Time::time_unix(1423009851.2847152));

    for (int i = 0; i < (int) closest_5.size(); i++) {
        BOOST_CHECK_EQUAL(closest_5[i].get()->time(spec_index), expected_times[i]);
    }
}


BOOST_AUTO_TEST_SUITE_END()
