// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "fp_common.i"

%{
#include "level_1b_info.h"
%}

%base_import(generic_object)
%import "fp_time.i"
%import "level_1b.i"
%import "double_with_unit.i"
%fp_shared_ptr(FullPhysics::Level1bInfo);


namespace FullPhysics {
class Level1bInfo {
public:
    virtual std::vector<boost::shared_ptr<Level1b>> level1b_list() = 0;
    virtual std::vector<boost::shared_ptr<Level1b>> close_obs(double lat, double lon, DoubleWithUnit threshold, int spec_index=0);
    virtual std::vector<boost::shared_ptr<Level1b>> closest_obs_n(double lat, double lon, int n, int spec_index=0);
    virtual boost::shared_ptr<Level1b> closest_obs(double lat, double lon, int spec_index=0);
    virtual std::vector<boost::shared_ptr<Level1b>> close_obs(Time search_time, Time threshold, int spec_index=0);
    virtual std::vector<boost::shared_ptr<Level1b>> closest_obs_n(Time search_time, int n, int spec_index=0);
    virtual boost::shared_ptr<Level1b> closest_obs(Time search_time, int spec_index=0);
    static DoubleWithUnit obs_distance(boost::shared_ptr<Level1b>, double lat, double lon, int spec_index=0);
    virtual ~Level1bInfo() { };
private:
    static bool compare_level1b_distance(boost::shared_ptr<Level1b>& lev1b_1, boost::shared_ptr<Level1b>& lev1b_2, double lat, double lon, int spec_index=0);
    static bool compare_level1b_time(boost::shared_ptr<Level1b>& level1b_1, boost::shared_ptr<Level1b>& level1b_2, Time search_time, int spec_index=0);
};
}

%template(vector_level1b) std::vector<boost::shared_ptr<FullPhysics::Level1b>>;