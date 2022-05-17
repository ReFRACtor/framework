#include "pressure_with_cloud_handling.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void PressureWithCloudHandling::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Pressure)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverPressure)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GenericObjectWithCloudHandling)
    & FP_NVP_(pressure_clear) 
    & FP_NVP_(cloud_pressure_level);
}

FP_IMPLEMENT(PressureWithCloudHandling);
#endif

//-----------------------------------------------------------------------
/// Nominal Constructor.
//-----------------------------------------------------------------------

PressureWithCloudHandling::PressureWithCloudHandling
(const boost::shared_ptr<Pressure> Press_clear,
 double Cloud_pressure_level, bool Do_cloud)
: GenericObjectWithCloudHandling(Do_cloud),
  pressure_clear_(Press_clear),
  cloud_pressure_level_(Cloud_pressure_level)
{
  pressure_clear_->add_observer(*this);
  type_preference_ = pressure_clear_->type_preference();
}

//-----------------------------------------------------------------------
/// Clone a PressureWithCloudHandling object. Note that the cloned
/// version will *not* be attached to a StateVector or
/// Observer<PressureSigma>, although you can of course attach them
/// after receiving the cloned object.
/// -----------------------------------------------------------------------

boost::shared_ptr<Pressure> PressureWithCloudHandling::clone() const
{
  return boost::make_shared<PressureWithCloudHandling>(pressure_clear_->clone(),
                                                       cloud_pressure_level_,
                                                       do_cloud());
}

void PressureWithCloudHandling::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "PressureWithCloudHandling:\n"
     << "  cloud pressure level:    " << cloud_pressure_level() << " Pa\n"
     << "  do_cloud: " << do_cloud() << "\n"
     << "  pressure clear: \n";
  opad << *pressure_clear() << "\n";
  opad.strict_sync();
}

void PressureWithCloudHandling::calc_pressure_grid() const
{
  ArrayAdWithUnit<double, 1> full = pressure_clear_->pressure_grid(Pressure::NATIVE_ORDER);
  cache.pgrid.units = full.units;

  if(!do_cloud()) {
    cache.pgrid.value.reference(full.value);
    return;
  }

  if(type_preference() == Pressure::PREFER_INCREASING_PRESSURE) {
    // increasing pressure
    int i;
    for(i = 0; i < full.rows(); ++i) {
      double v = full(i).convert(units::Pa).value.value();
      if(v > cloud_pressure_level_)
        break;
    }

    // if cloud_pressure_level_ is negative or very low
    //  blitz::Range(0,i-1) will be empty
    if (i == 0) {
      std::stringstream err_msg;
      err_msg << "Cloud pressure level too low: "
              << cloud_pressure_level_
              << " . Pressure grid would be empty.";
      throw Exception(err_msg.str());
    }

    cache.pgrid.value.reference( full.value(blitz::Range(0,i-1)) );

  } else {
    // decreasing pressure
    int i;
    for(i = full.rows() - 1; i >= 0; --i) {
      double v = full(i).convert(units::Pa).value.value();
      if(v > cloud_pressure_level_)
        break;
    }

    // if cloud_pressure_level_ is negative or very low
    //  blitz::Range(i+1,blitz::toEnd) will be empty
    if (i == full.rows() - 1) {
      std::stringstream err_msg;
      err_msg << "Cloud pressure level too low: "
              << cloud_pressure_level_
              << " . Pressure grid would be empty.";
      throw Exception(err_msg.str());
    }

    cache.pgrid.value.reference( full.value(blitz::Range(i+1,blitz::toEnd)) );

  }
}
