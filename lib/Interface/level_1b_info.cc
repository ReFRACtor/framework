#include <functional>
#include <boost/lexical_cast.hpp>
#include <boost/geometry.hpp>

#include "level_1b_info.h"


using namespace FullPhysics;
using namespace blitz;


// TODO: DoubleWithUnit is probably a smarter choice here
double Level1bInfo::obs_distance(boost::shared_ptr<Level1b> level1b_ptr, double lat, double lon) {
    Level1b* level1b = level1b_ptr.get();
    // TODO: Add support for varying instrument index
    DoubleWithUnit curr_lat = level1b->latitude(0);
    DoubleWithUnit curr_lon = level1b->longitude(0);        
    // TODO: Add earth radius someplace globally?
    double const earth_radius = 6378.2064; // From Clarke 1866
    // TODO: Add typedef for spherical_point someplace globally?
    typedef boost::geometry::model::point<double, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::degree>> spherical_point;
    spherical_point in_point(lon, lat);
    spherical_point curr_point(curr_lon.value, curr_lat.value);
    // use strategy::distance::haversine to get distance (in km)
    double seperation = boost::geometry::distance(in_point, curr_point) * earth_radius;
    return seperation;   
}

bool Level1bInfo::compare_level1b_distance(boost::shared_ptr<Level1b>& level1b_1, boost::shared_ptr<Level1b>& level1b_2, double lat, double lon) {
    double obs1_dist = std::abs(Level1bInfo::obs_distance(level1b_1, lat, lon));
    double obs2_dist = std::abs(Level1bInfo::obs_distance(level1b_2, lat, lon));
    return (obs1_dist < obs2_dist);
}

bool Level1bInfo::compare_level1b_time(boost::shared_ptr<Level1b>& level1b_1, boost::shared_ptr<Level1b>& level1b_2, Time search_time) {
    // TODO: Add support for varying spec index
    double obs1_timediff = std::abs((level1b_1->time(0) - search_time));
    double obs2_timediff = std::abs((level1b_2->time(0) - search_time));
    return (obs1_timediff < obs2_timediff);
}

std::vector<boost::shared_ptr<Level1b>> Level1bInfo::closest_obs_n(double lat, double lon, int n) {
    std::vector<boost::shared_ptr<Level1b>> all_level1b = level1b_list();
    auto compare_level1b_latlon = std::bind(Level1bInfo::compare_level1b_distance, std::placeholders::_1, std::placeholders::_2, lat, lon);
    std::sort(all_level1b.begin(), all_level1b.end(), compare_level1b_latlon);

    std::vector<boost::shared_ptr<Level1b>> closest_observations = std::vector<boost::shared_ptr<Level1b>>();
    for (int obs_num = 0; (obs_num < n) && (obs_num < all_level1b.size()); obs_num++) {
        closest_observations.push_back(all_level1b[obs_num]);
    }
    return closest_observations;
}

boost::shared_ptr<Level1b> Level1bInfo::closest_obs(double lat, double lon) {
    std::vector<boost::shared_ptr<Level1b>> closest = closest_obs_n(lat, lon, 1);
    // TODO: What do we return if no obs from above? all_observation.size() == 0
    return closest.front();    
}

std::vector<boost::shared_ptr<Level1b>> Level1bInfo::closest_obs_n(Time search_time, int n) {
    std::vector<boost::shared_ptr<Level1b>> all_level1b = level1b_list();
    auto compare_level1b_timediff = std::bind(Level1bInfo::compare_level1b_time, std::placeholders::_1, std::placeholders::_2, search_time);
    std::sort(all_level1b.begin(), all_level1b.end(), compare_level1b_timediff);

    std::vector<boost::shared_ptr<Level1b>> closest_observations = std::vector<boost::shared_ptr<Level1b>>();
    for (int obs_num = 0; (obs_num < n) && (obs_num < all_level1b.size()); obs_num++) {
        closest_observations.push_back(all_level1b[obs_num]);
    }
    return closest_observations;
}

boost::shared_ptr<Level1b> Level1bInfo::closest_obs(Time search_time) {
    std::vector<boost::shared_ptr<Level1b>> closest = closest_obs_n(search_time, 1);
    // TODO: What do we return if no obs from above? all_observation.size() == 0
    return closest.front();    
}


