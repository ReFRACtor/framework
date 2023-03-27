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
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GenericObjectWithCloudHandling)
    & FP_NVP_(ground_clear)  & FP_NVP_(ground_cloud);
}

FP_IMPLEMENT(GroundWithCloudHandling);
#endif

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------
  
GroundWithCloudHandling::GroundWithCloudHandling
(const boost::shared_ptr<Ground>& Ground_clear,
 const boost::shared_ptr<Ground>& Ground_cloud,
 bool Do_cloud)
: GenericObjectWithCloudHandling(Do_cloud),
  ground_clear_(Ground_clear), ground_cloud_(Ground_cloud)
{
  ground_clear_->add_observer(*this);
  ground_cloud_->add_observer(*this);
}

ArrayAd<double, 1> GroundWithCloudHandling::surface_parameter
(const double wn, const int spec_index) const
{
  if(do_cloud())
    return ground_cloud_->surface_parameter(wn, spec_index);
  return ground_clear_->surface_parameter(wn, spec_index);
}

boost::shared_ptr<Ground> GroundWithCloudHandling::clone() const
{
  return boost::make_shared<GroundWithCloudHandling>(ground_clear_->clone(),
						     ground_cloud_->clone(),
						     do_cloud());
}

void GroundWithCloudHandling::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "GroundWithCloudHandling:\n"
     << "  do_cloud: " << do_cloud() << "\n"
     << "  ground cloud: \n";
  opad << *ground_cloud() << "\n";
  opad.strict_sync();
  Os << "  ground clear: \n";
  opad << *ground_clear() << "\n";
  opad.strict_sync();
}
