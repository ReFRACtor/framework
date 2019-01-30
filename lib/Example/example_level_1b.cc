#include <functional>
#include <boost/lexical_cast.hpp>
#include <boost/geometry.hpp>


#include "example_level_1b.h"
#include "fp_exception.h"
#include "example_observation_id.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ExampleLevel1b, Level1bSampleCoefficient)
.def(luabind::constructor<const boost::shared_ptr<HdfFile>&, const std::string&>())
REGISTER_LUA_END()
#endif

ExampleLevel1b::ExampleLevel1b(const boost::shared_ptr<HdfFile>& input_file, const std::string& observation_id)
:input(input_file)
{
    ExampleObservationId<std::string> obs_id(input_file, "observation_ids", observation_id);
    data_index = obs_id.data_index();
}

ExampleLevel1b::ExampleLevel1b(const std::string& input_filename, const std::string& observation_id)
:input(new HdfFile(input_filename))
{
    ExampleObservationId<std::string> obs_id(input, "observation_ids", observation_id);
    data_index = obs_id.data_index();
}

ExampleLevel1b::ExampleLevel1b(const boost::shared_ptr<HdfFile>& input_file, ExampleObservationId<std::string> observation_id)
:input(input_file)
{
    data_index = observation_id.data_index();
}

ExampleLevel1b::ExampleLevel1b(const std::string& input_filename, ExampleObservationId<std::string> observation_id)
:input(new HdfFile(input_filename))
{
    data_index = observation_id.data_index();
}

int ExampleLevel1b::number_spectrometer() const
{
    TinyVector<int, 2> lat_shape = input->read_shape<2>(group_name + "/latitude");
    return lat_shape(1);
}

SpectralRange ExampleLevel1b::radiance(int Spec_index) const
{
    std::string rad_ds_name = group_name + "/radiance_" + boost::lexical_cast<std::string>(Spec_index + 1);
    std::string uncert_ds_name = group_name + "/uncertainty_" + boost::lexical_cast<std::string>(Spec_index + 1);

    TinyVector<int, 2> rad_ds_shape = input->read_shape<2>(rad_ds_name);
    TinyVector<int, 2> rad_start, rad_size;
    rad_start = data_index, 0;
    rad_size  = 1, rad_ds_shape(1);
    ArrayWithUnit<double, 2> radiance = input->read_field_with_unit<double, 2>(rad_ds_name, 
        Unit("Ph sec^{-1} m^{-2} sr^{-1} um^{-1}"),
        rad_start, rad_size);

    if (input->has_object(uncert_ds_name)) {
        TinyVector<int, 2> uncert_ds_shape = input->read_shape<2>(uncert_ds_name);
        TinyVector<int, 2> uncert_start, uncert_size;
        uncert_start = data_index, 0;
        uncert_size  = 1, uncert_ds_shape(1);
        Array<double, 2> uncertainty = input->read_field<double, 2>(uncert_ds_name, uncert_start, uncert_size);

        return SpectralRange(radiance(0, Range::all()).value, radiance.units, uncertainty(0, Range::all()));
    } else {
        return SpectralRange(radiance(0, Range::all()).value, radiance.units);
    }
}

DoubleWithUnit ExampleLevel1b::read_scalar_with_unit(const std::string& dataset_name, int i, const Unit& default_unit) const
{
    range_check(i, 0, number_spectrometer());
    TinyVector<int, 2> start, size;
    start = data_index, i;
    size = 1, 1;
    ArrayWithUnit<double, 2> val_arr = input->read_field_with_unit<double, 2>(dataset_name, default_unit, start, size);
    return val_arr(0, 0);
}

double ExampleLevel1b::read_scalar(const std::string& dataset_name, int i) const
{
    range_check(i, 0, number_spectrometer());
    TinyVector<int, 2> start, size;
    start = data_index, i;
    size = 1, 1;
    Array<double, 2> val_arr = input->read_field<double, 2>(dataset_name, start, size);
    return val_arr(0, 0);
}

ArrayWithUnit<double, 1> ExampleLevel1b::read_array_with_unit(const std::string& dataset_name, int i, const Unit& default_unit) const
{
    range_check(i, 0, number_spectrometer());

    TinyVector<int, 3> ds_shape = input->read_shape<3>(dataset_name);
    TinyVector<int, 3> start, size;
    start = data_index, i, 0;
    size = 1, 1, ds_shape(2);

    ArrayWithUnit<double, 3> val_arr = input->read_field_with_unit<double, 3>(dataset_name, default_unit, start, size);

    ArrayWithUnit<double, 1> ret_val;
    ret_val.units = val_arr.units;
    ret_val.value.resize(ds_shape(2));
    ret_val.value = val_arr.value(0, 0, Range::all());

    return ret_val;
}

Array<double, 1> ExampleLevel1b::read_array(const std::string& dataset_name, int i) const
{
    range_check(i, 0, number_spectrometer());
    TinyVector<int, 3> ds_shape = input->read_shape<3>(dataset_name);
    TinyVector<int, 3> start, size;
    start = data_index, i, 0;
    size = 1, 1, ds_shape(2);
    Array<double, 3> val_arr = input->read_field<double, 3>(dataset_name, start, size);
    return val_arr(0, 0, Range::all());
}

std::vector<ExampleObservationId<std::string>> ExampleLevel1b::obs_list(const std::string& geo_input_filename) {
    std::string dataset_name = "/observation_ids";
    boost::shared_ptr<HdfFile> input_hdf(new HdfFile(geo_input_filename));
    TinyVector<int, 1> ds_shape = input_hdf->read_shape<1>(dataset_name);
    TinyVector<int, 1> start, size;
    start = 0;
    size = ds_shape(0);
    Array<std::string, 1> obs_arr = input_hdf->read_field<std::string, 1>(dataset_name, start, size);

    std::vector<ExampleObservationId<std::string>> observations = std::vector<ExampleObservationId<std::string>>();
    for (auto obs_str : obs_arr) {
        ExampleObservationId<std::string> new_obs(input_hdf, "observation_ids",obs_str);
        observations.push_back(new_obs);
    }

    return observations;

}

bool ExampleLevel1b::compare_obs_distance(ExampleObservationId<std::string> obs1, ExampleObservationId<std::string> obs2, const std::string& geo_input_filename, double lat, double lon)
{
    /*
    double obs1_dist = std::abs(ExampleLevel1b::obs_distance(geo_input_filename, obs1, lat, lon));
    double obs2_dist = std::abs(ExampleLevel1b::obs_distance(geo_input_filename, obs2, lat, lon));
    return (obs1_dist < obs2_dist);
    */
}

bool ExampleLevel1b::compare_obs_time(ExampleObservationId<std::string> obs1, ExampleObservationId<std::string> obs2, const std::string& geo_input_filename, Time search_time)
{
    ExampleLevel1b level1b_1 (geo_input_filename, obs1);
    ExampleLevel1b level1b_2 (geo_input_filename, obs2);
    // TODO: Add support for varying instrument index
    double obs1_timediff = std::abs(level1b_1.time(0) - search_time);
    double obs2_timediff = std::abs(level1b_2.time(0) - search_time);
    return (obs1_timediff < obs2_timediff);
}

std::vector<ExampleObservationId<std::string>> ExampleLevel1b::closest_obs_n(const std::string& geo_input_filename, double lat, double lon, int n) {
    std::vector<ExampleObservationId<std::string>> all_observations = ExampleLevel1b::obs_list(geo_input_filename);

    auto compare_obs_latlon = std::bind(ExampleLevel1b::compare_obs_distance, std::placeholders::_1, std::placeholders::_2, geo_input_filename, lat, lon);
    std::sort(all_observations.begin(), all_observations.end(), compare_obs_latlon);

    std::vector<ExampleObservationId<std::string>> closest_observations = std::vector<ExampleObservationId<std::string>>();
    for (int obs_num = 0; (obs_num < n) && (obs_num < all_observations.size()); obs_num++) {
        closest_observations.push_back(all_observations[obs_num]);
    }
    return closest_observations;
}

ExampleObservationId<std::string> ExampleLevel1b::closest_obs(const std::string& geo_input_filename, double lat, double lon) {
    std::vector<ExampleObservationId<std::string>> closest = closest_obs_n(geo_input_filename, lat, lon, 1);
    // TODO: What do we return if no obs from above? all_observation.size() == 0
    return closest.front();

}

std::vector<ExampleObservationId<std::string>> ExampleLevel1b::closest_obs_n(const std::string& geo_input_filename, Time search_time, int n) {
    std::vector<ExampleObservationId<std::string>> all_observations = ExampleLevel1b::obs_list(geo_input_filename);

    auto compare_obs_timediff = std::bind(ExampleLevel1b::compare_obs_time, std::placeholders::_1, std::placeholders::_2, geo_input_filename, search_time);
    std::sort(all_observations.begin(), all_observations.end(), compare_obs_timediff);

    std::vector<ExampleObservationId<std::string>> closest_observations = std::vector<ExampleObservationId<std::string>>();
    for (int obs_num = 0; (obs_num < n) && (obs_num < all_observations.size()); obs_num++) {
        closest_observations.push_back(all_observations[obs_num]);
    }
    return closest_observations;
}

ExampleObservationId<std::string> ExampleLevel1b::closest_obs(const std::string& geo_input_filename, Time search_time) {
    std::vector<ExampleObservationId<std::string>> closest = closest_obs_n(geo_input_filename, search_time, 1);
    // TODO: What do we return if no obs from above? all_observation.size() == 0
    return closest.front();

}

/*double ExampleLevel1b::obs_distance(const std::string& geo_input_filename, ExampleObservationId<std::string> obs_id, double lat, double lon) {
    ExampleLevel1b curr_level1b (geo_input_filename, obs_id);
    // TODO: Add support for varying instrument index
    DoubleWithUnit curr_lat = curr_level1b.latitude(0);
    DoubleWithUnit curr_lon = curr_level1b.longitude(0);
    // use strategy::distance::haversine to get distance (in km)
    double const earth_radius = 6378.1;
    typedef boost::geometry::model::point<double, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::degree>> spherical_point;
    spherical_point in_point(lon, lat);
    spherical_point curr_point(curr_lon.value, curr_lat.value);
    double seperation = boost::geometry::distance(in_point, curr_point) * earth_radius;
    return seperation;
}*/
