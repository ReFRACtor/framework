#ifndef PRESSURE_WITH_CLOUD_HANDLING_H
#define PRESSURE_WITH_CLOUD_HANDLING_H
#include "observer.h"
#include "fp_exception.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
  This class is an adapter, than takes an existing clear Pressure
  object and handles support for working with cloudy data. If the
  variable "do_cloud" is false, then we just act like the clear
  Pressure object. If it is true, then the pressure grid it truncated 
  at the supplied cloud top pressure level.  This is used to support
  the algorithm done in muses-py, where we do 2 forward model calls
  to generate radiance, once for clear and once for cloudy. The two
  calls are then combined with a cloud fraction - see 
  ForwardModelWithCloudFraction.
*******************************************************************/
class PressureWithCloudHandling : virtual public Pressure,
				  public Observer<Pressure> {
public:
  PressureWithCloudHandling(const boost::shared_ptr<Pressure> Press_clear,
			    double Cloud_pressure_level, bool do_cloud = false);
  virtual ~PressureWithCloudHandling() {}
//-----------------------------------------------------------------------
/// If true, then truncate Press_clear at Cloud_pressure_level,
/// otherwise just return pressure levels from Pres_clear
//-----------------------------------------------------------------------

  bool do_cloud() const { return do_cloud_;}
  void do_cloud(bool F)
  {
    do_cloud_ = F;
    Observable<Pressure>::notify_update_do(*this);
  }

//-----------------------------------------------------------------------
/// The clear pressure object
//-----------------------------------------------------------------------

  const boost::shared_ptr<Pressure> pressure_clear() const
  { return pressure_clear_; }

//-----------------------------------------------------------------------
/// Cloud pressure level, is Pascals
//-----------------------------------------------------------------------

  double cloud_pressure_level() const { return cloud_pressure_level_; }
  void cloud_pressure_level(double P)
  {
    cloud_pressure_level_ = P;
    Observable<Pressure>::notify_update_do(*this);
  }

  virtual void notify_update(const Pressure& UNUSED(P) )
  {
    // Forward notification from pressure_clear to anything
    // observing this object.
    Observable<Pressure>::notify_update_do(*this);
  }
  virtual int max_number_level() const
  {
    return pressure_clear_->max_number_level();
  }
  virtual ArrayAdWithUnit<double, 1> pressure_grid() const;
  virtual void print(std::ostream& Os) const;
  virtual boost::shared_ptr<Pressure> clone() const;
private:
  boost::shared_ptr<Pressure> pressure_clear_;
  bool do_cloud_;
  double cloud_pressure_level_;
  PressureWithCloudHandling() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(PressureWithCloudHandling);
#endif
