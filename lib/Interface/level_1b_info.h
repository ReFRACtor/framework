#ifndef LEVEL_1B_INFO_H
#define LEVEL_1B_INFO_H

#include "fp_time.h"
#include "level_1b.h"
#include "double_with_unit.h"


namespace FullPhysics {
/****************************************************************//**
  Used to get observation information from L1B readers
*******************************************************************/

class Level1bInfo : public Printable<GenericObject> {
public:

//-----------------------------------------------------------------------
/// Returns a list containing all valid level1b observations
//-----------------------------------------------------------------------
    virtual std::vector<boost::shared_ptr<Level1b>> level1b_list() = 0;

//-----------------------------------------------------------------------
/// Returns level1b observations that are closer than the threshold
/// to the lat/lon
/// May be empty
//-----------------------------------------------------------------------
    virtual std::vector<boost::shared_ptr<Level1b>> close_obs(double lat, double lon, DoubleWithUnit threshold, int spec_index=0);


//-----------------------------------------------------------------------
/// Returns "n" level1b observations that are closest to the lat/lon
/// May be empty or contain less than n observations
//-----------------------------------------------------------------------
    virtual std::vector<boost::shared_ptr<Level1b>> closest_obs_n(double lat, double lon, int n, int spec_index=0);
    
//-----------------------------------------------------------------------
/// Returns the level1b observation that is closest to the lat/lon
/// May be NULL
//-----------------------------------------------------------------------
    virtual boost::shared_ptr<Level1b> closest_obs(double lat, double lon, int spec_index=0);
    

//-----------------------------------------------------------------------
/// Returns level1b observations that are closer than the threshold
/// to the search_time
/// May be empty
//-----------------------------------------------------------------------
    virtual std::vector<boost::shared_ptr<Level1b>> close_obs(Time search_time, Time threshold, int spec_index=0);

//-----------------------------------------------------------------------
/// Returns "n" level1b observations that are closest to the search_time
/// May be empty or contain less than n observations
//-----------------------------------------------------------------------
    virtual std::vector<boost::shared_ptr<Level1b>> closest_obs_n(Time search_time, int n, int spec_index=0);

//-----------------------------------------------------------------------
/// Returns the level1b observation that is closest to the search_time
/// May be NULL
//-----------------------------------------------------------------------
    virtual boost::shared_ptr<Level1b> closest_obs(Time search_time, int spec_index=0);
    
//-----------------------------------------------------------------------
/// Returns the distance from given level1b observation to the lat/lon
//-----------------------------------------------------------------------
    static DoubleWithUnit obs_distance(boost::shared_ptr<Level1b>, double lat, double lon, int spec_index=0);

    virtual ~Level1bInfo() { };

private:
//-----------------------------------------------------------------------
/// Internal comparator for use in sorting level1b observations by lat/lon
//-----------------------------------------------------------------------
    static bool compare_level1b_distance(const boost::shared_ptr<Level1b>& lev1b_1, const boost::shared_ptr<Level1b>& lev1b_2, double lat, double lon, int spec_index=0);

//-----------------------------------------------------------------------
/// Internal comparator for use in sorting level1b observations by time
//-----------------------------------------------------------------------
    static bool compare_level1b_time(const boost::shared_ptr<Level1b>& level1b_1, const boost::shared_ptr<Level1b>& level1b_2, Time search_time, int spec_index=0);

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "Level1bInfo";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(Level1bInfo);

#endif
