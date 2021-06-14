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
    & FP_NVP_(pressure_clear) & FP_NVP_(do_cloud)
    & FP_NVP_(cloud_pressure_level);
}

FP_IMPLEMENT(PressureWithCloudHandling);
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

PressureWithCloudHandling::PressureWithCloudHandling
(const boost::shared_ptr<Pressure> Press_clear,
 double Cloud_pressure_level, bool Do_cloud)
  : pressure_clear_(Press_clear),
    do_cloud_(Do_cloud),
    cloud_pressure_level_(Cloud_pressure_level)
{
  pressure_clear_->add_observer(*this);
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
						       do_cloud_);
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

ArrayAdWithUnit<double, 1>
PressureWithCloudHandling::pressure_grid() const
{
  ArrayAdWithUnit<double, 1> full = pressure_clear_->pressure_grid();
  if(!do_cloud_)
    return full;
  int i;
  for(i = 0; i < full.rows(); ++i) {
    double v = full(i).convert(Unit("Pa")).value.value();
    if(v > cloud_pressure_level_)
      break;
  }
  return full(blitz::Range(0,i-1));
}

