#include "ground_with_cloud_handling.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void GroundWithCloudHandling::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Ground)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverGround)
    & FP_NVP_(ground_clear) & FP_NVP_(do_cloud) & FP_NVP_(cloud_albedo);
}

FP_IMPLEMENT(GroundWithCloudHandling);
#endif

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------
  
GroundWithCloudHandling::GroundWithCloudHandling
(const boost::shared_ptr<Ground> Ground_clear,
 double Cloud_albedo, bool Do_cloud)
  : ground_clear_(Ground_clear), cloud_albedo_(Cloud_albedo),
    do_cloud_(Do_cloud)
{
  ground_clear_->add_observer(*this);
}

ArrayAd<double, 1> GroundWithCloudHandling::surface_parameter
(const double wn, const int spec_index) const
{
  if(do_cloud_) {
    ArrayAd<double, 1> spars(1, 0);
    spars(0) = cloud_albedo_;
    return spars;
  }
  return ground_clear_->surface_parameter(wn, spec_index);
}

boost::shared_ptr<Ground> GroundWithCloudHandling::clone() const
{
  return boost::make_shared<GroundWithCloudHandling>(ground_clear_->clone(),
						     cloud_albedo_, do_cloud_);
}

void GroundWithCloudHandling::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "GroundWithCloudHandling:\n"
     << "  cloud albedo:    " << cloud_albedo() << "\n"
     << "  do_cloud: " << do_cloud() << "\n"
     << "  ground clear: \n";
  opad << *ground_clear() << "\n";
  opad.strict_sync();
}
