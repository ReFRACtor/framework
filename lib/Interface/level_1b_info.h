#ifndef LEVEL_1B_INFO_H
#define LEVEL_1B_INFO_H

#include <memory>

#include "fp_time.h"
#include "level_1b.h"
#include "observation_id.h"


namespace FullPhysics {
/****************************************************************//**
  Used to get observation information from L1B readers
*******************************************************************/

class Level1bInfo {
public:
    virtual std::vector<std::shared_ptr<Level1b>> level1b_list() = 0;
    
    virtual std::vector<std::shared_ptr<Level1b>> closest_obs_n(double lat, double lon, int n);
    
    virtual std::shared_ptr<Level1b> closest_obs(double lat, double lon);
    
    virtual std::vector<std::shared_ptr<Level1b>> closest_obs_n(Time search_time, int n);
    
    virtual std::shared_ptr<Level1b> closest_obs(Time search_time);
    
    static double obs_distance(std::shared_ptr<Level1b>, double lat, double lon);

    virtual ~Level1bInfo() { };

private:
    static bool compare_level1b_distance(std::shared_ptr<Level1b>& lev1b_1, std::shared_ptr<Level1b>& lev1b_2, double lat, double lon);
    static bool compare_level1b_time(std::shared_ptr<Level1b>& level1b_1, std::shared_ptr<Level1b>& level1b_2, Time search_time);

};

}
#endif
