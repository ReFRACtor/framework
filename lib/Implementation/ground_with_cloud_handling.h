#ifndef GROUND_WITH_CLOUD_HANDLING_H
#define GROUND_WITH_CLOUD_HANDLING_H

#include "ground.h"

namespace FullPhysics {
/****************************************************************//**
  This class is like PressureWithCloudHandling, it adds handling
  to alter behavior based on a do_cloud flag. We are either a 
  "clear" Ground, or we are a lambertian surface of a fixed albedo.
*******************************************************************/

class GroundWithCloudHandling: virtual public Ground,
			       public Observer<Ground>{
public:
  GroundWithCloudHandling(const boost::shared_ptr<Ground> Ground_clear,
			  double Cloud_albedo, bool Do_cloud = false);
  virtual ~GroundWithCloudHandling() {}

//-----------------------------------------------------------------------
/// If true, then truncate Press_clear at Cloud_pressure_level,
/// otherwise just return pressure levels from Pres_clear
//-----------------------------------------------------------------------

  bool do_cloud() const { return do_cloud_;}
  void do_cloud(bool F)
  {
    do_cloud_ = F;
    Observable<Ground>::notify_update_do(*this);
  }

//-----------------------------------------------------------------------
/// The clear ground object
//-----------------------------------------------------------------------

  const boost::shared_ptr<Ground> ground_clear() const
  { return ground_clear_; }

//-----------------------------------------------------------------------
/// The cloud albedo.
//-----------------------------------------------------------------------

  double cloud_albedo() const { return cloud_albedo_; }
  void cloud_albedo(double V)
  {
    cloud_albedo_ = V;
    Observable<Ground>::notify_update_do(*this);
  }

  virtual void notify_update(const Ground& UNUSED(G) )
  {
    // Forward notification from ground_clear to anything
    // observing this object.
    Observable<Ground>::notify_update_do(*this);
  }

  virtual SpurrBrdfType spurr_brdf_type() const
  {
    if(do_cloud_)
      return SpurrBrdfType::LAMBERTIAN;
    return ground_clear_->spurr_brdf_type();
  }
  virtual ArrayAd<double, 1> surface_parameter
  (const double wn, const int spec_index) const;
  virtual boost::shared_ptr<Ground> clone() const;
  virtual void print(std::ostream& Os) const;
private:
  boost::shared_ptr<Ground> ground_clear_;
  double cloud_albedo_;
  bool do_cloud_;
  GroundWithCloudHandling() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(GroundWithCloudHandling);
#endif
