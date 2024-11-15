#ifndef GROUND_WITH_CLOUD_HANDLING_H
#define GROUND_WITH_CLOUD_HANDLING_H

#include "ground.h"
#include "generic_object_with_cloud_handling.h"

namespace FullPhysics {
/****************************************************************//**
  This class is like PressureWithCloudHandling, it adds handling
  to alter behavior based on a do_cloud flag. We are either a 
  "clear" Ground, or we are a lambertian surface of a fixed albedo.
*******************************************************************/

class GroundWithCloudHandling: virtual public Ground,
			       public Observer<Ground>,
			       public GenericObjectWithCloudHandling {
public:
  GroundWithCloudHandling(const boost::shared_ptr<Ground>& Ground_clear,
			  const boost::shared_ptr<Ground>& Ground_cloud,
			  bool Do_cloud = false);
  virtual ~GroundWithCloudHandling() {}

  virtual void notify_do_cloud_update()
  {
    Observable<Ground>::notify_update_do(*this);
  }

//-----------------------------------------------------------------------
/// The clear ground object
//-----------------------------------------------------------------------

  const boost::shared_ptr<Ground> ground_clear() const
  { return ground_clear_; }


//-----------------------------------------------------------------------
/// The cloud ground object
//-----------------------------------------------------------------------

  const boost::shared_ptr<Ground> ground_cloud() const
  { return ground_cloud_; }

  virtual void notify_update(const Ground& UNUSED(G) )
  {
    // Forward notification from ground_clear to anything
    // observing this object.
    Observable<Ground>::notify_update_do(*this);
  }

  virtual SpurrBrdfType spurr_brdf_type() const
  {
    if(do_cloud())
      return ground_cloud_->spurr_brdf_type();
    return ground_clear_->spurr_brdf_type();
  }
  virtual ArrayAd<double, 1> surface_parameter
  (const double wn, const int spec_index) const;
  virtual boost::shared_ptr<Ground> clone() const;
  virtual void print(std::ostream& Os) const;
private:
  boost::shared_ptr<Ground> ground_clear_;
  boost::shared_ptr<Ground> ground_cloud_;
  GroundWithCloudHandling() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(GroundWithCloudHandling);
#endif
