#include <functional>
#include <boost/lexical_cast.hpp>
#include <boost/geometry.hpp>

#include "level_1b_info.h"
#include "default_constant.h"


using namespace FullPhysics;
using namespace blitz;


DoubleWithUnit Level1bInfo::obs_distance(boost::shared_ptr<Level1b> level1b_ptr, double lat, double lon, int spec_index) {
    DefaultConstant cs;
    Level1b* level1b = level1b_ptr.get();
    DoubleWithUnit curr_lat = level1b->latitude(spec_index);
    DoubleWithUnit curr_lon = level1b->longitude(spec_index);
    // TODO: Add typedef for spherical_point someplace globally?
    typedef boost::geometry::model::point<double, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::degree>> spherical_point;
    spherical_point in_point(lon, lat);
    spherical_point curr_point(curr_lon.value, curr_lat.value);
    // Uses strategy::distance::haversine to get distance
    DoubleWithUnit seperation = boost::geometry::distance(in_point, curr_point) * cs.equatorial_radius();
    return seperation;   
}

bool Level1bInfo::compare_level1b_distance(const boost::shared_ptr<Level1b>& level1b_1, const boost::shared_ptr<Level1b>& level1b_2, double lat, double lon, int spec_index) {
    DoubleWithUnit obs1_dist = Level1bInfo::obs_distance(level1b_1, lat, lon, spec_index);
    DoubleWithUnit obs2_dist = Level1bInfo::obs_distance(level1b_2, lat, lon, spec_index).convert(obs1_dist.units);
    return (std::abs(obs1_dist.value) < std::abs(obs2_dist.value));
}

bool Level1bInfo::compare_level1b_time(const boost::shared_ptr<Level1b>& level1b_1, const boost::shared_ptr<Level1b>& level1b_2, Time search_time, int spec_index) {
    double obs1_timediff = std::abs((level1b_1->time(spec_index) - search_time));
    double obs2_timediff = std::abs((level1b_2->time(spec_index) - search_time));
    return (obs1_timediff < obs2_timediff);
}

std::vector<boost::shared_ptr<Level1b>> Level1bInfo::close_obs(double lat, double lon, DoubleWithUnit threshold, int spec_index) {
    std::vector<boost::shared_ptr<Level1b>> all_level1b = level1b_list();

    std::vector<boost::shared_ptr<Level1b>> close_observations = std::vector<boost::shared_ptr<Level1b>>();
    for(auto level1b : all_level1b) {
        DoubleWithUnit dist = Level1bInfo::obs_distance(level1b, lat, lon, spec_index).convert(threshold.units);
        if (std::abs(threshold.value) - std::abs(dist.value) > 0) {
            close_observations.push_back(level1b);
        }
    }
    return close_observations;
}

std::vector<boost::shared_ptr<Level1b>> Level1bInfo::closest_obs_n(double lat, double lon, int n, int spec_index) {
    std::vector<boost::shared_ptr<Level1b>> all_level1b = level1b_list();
    auto compare_level1b_latlon = std::bind(Level1bInfo::compare_level1b_distance, std::placeholders::_1, std::placeholders::_2, lat, lon, spec_index);
    std::sort(all_level1b.begin(), all_level1b.end(), compare_level1b_latlon);

    std::vector<boost::shared_ptr<Level1b>> closest_observations = std::vector<boost::shared_ptr<Level1b>>();
    for (int obs_num = 0; (obs_num < n) && (obs_num < all_level1b.size()); obs_num++) {
        closest_observations.push_back(all_level1b[obs_num]);
    }
    return closest_observations;
}

boost::shared_ptr<Level1b> Level1bInfo::closest_obs(double lat, double lon, int spec_index) {
    std::vector<boost::shared_ptr<Level1b>> closest = closest_obs_n(lat, lon, 1, spec_index);
    if (closest.size()) {
        return closest.front();
    }
    return boost::shared_ptr<Level1b>();
}

std::vector<boost::shared_ptr<Level1b>> Level1bInfo::close_obs(Time search_time, Time threshold, int spec_index) {
    std::vector<boost::shared_ptr<Level1b>> all_level1b = level1b_list();

    std::vector<boost::shared_ptr<Level1b>> close_observations = std::vector<boost::shared_ptr<Level1b>>();
    for(auto level1b : all_level1b) {
        double timediff = std::abs(level1b->time(spec_index) - search_time);
        if (std::abs(threshold.unix_time()) - timediff > 0) {
            close_observations.push_back(level1b);
        }
    }

    return close_observations;
}

std::vector<boost::shared_ptr<Level1b>> Level1bInfo::closest_obs_n(Time search_time, int n, int spec_index) {
    std::vector<boost::shared_ptr<Level1b>> all_level1b = level1b_list();
    auto compare_level1b_timediff = std::bind(Level1bInfo::compare_level1b_time, std::placeholders::_1, std::placeholders::_2, search_time, spec_index);
    std::sort(all_level1b.begin(), all_level1b.end(), compare_level1b_timediff);

    std::vector<boost::shared_ptr<Level1b>> closest_observations = std::vector<boost::shared_ptr<Level1b>>();
    for (int obs_num = 0; (obs_num < n) && (obs_num < all_level1b.size()); obs_num++) {
        closest_observations.push_back(all_level1b[obs_num]);
    }
    return closest_observations;
}

boost::shared_ptr<Level1b> Level1bInfo::closest_obs(Time search_time, int spec_index) {
    std::vector<boost::shared_ptr<Level1b>> closest = closest_obs_n(search_time, 1, spec_index);
    if (closest.size()) {
        return closest.front();
    }
    return boost::shared_ptr<Level1b>();
}



