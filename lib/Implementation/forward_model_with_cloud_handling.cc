#include "forward_model_with_cloud_handling.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
SUB_STATE_VECTOR_ARRAY_SERIALIZE(CloudFraction, SubStateVectorArrayCloudFraction);

template<class Archive>
void CloudFraction::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableCloudFraction);
}

template<class Archive>
void CloudFractionFromState::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayCloudFraction);
}

template<class Archive>
void ForwardModelWithCloudHandling::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ForwardModel)
    & FP_NVP_(cloud_handling_vector)
    & FP_NVP_(fmodel) & FP_NVP_(cfrac);
}

FP_IMPLEMENT(CloudFraction);
FP_IMPLEMENT(CloudFractionFromState);
FP_OBSERVER_SERIALIZE(CloudFraction);
FP_IMPLEMENT(ForwardModelWithCloudHandling);
#endif

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

ForwardModelWithCloudHandling::ForwardModelWithCloudHandling
(const boost::shared_ptr<ForwardModel>& Fmodel,
 const boost::shared_ptr<CloudFraction>& Cfrac,
 const std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >&
 Cloud_handling_vector)
  : fmodel_(Fmodel),
    cfrac_(Cfrac),
    cloud_handling_vector_(Cloud_handling_vector)
{
}

Spectrum ForwardModelWithCloudHandling::radiance
(int UNUSED(channel_index), bool UNUSED(skip_jacobian)) const
{
  Spectrum dummy;
  return dummy;
}

void ForwardModelWithCloudHandling::print(std::ostream& Os) const	\
{
  OstreamPad opad(Os, "  ");
  Os << "FowardModelWithCloudHandling:\n"
     << "  cloud fraction:    " << cfrac_->cloud_fraction() << " \n"
     << "  underlying forward mode: \n";
  opad << *fmodel_ << "\n";
  opad.strict_sync();
}


void ForwardModelWithCloudHandling::set_do_cloud(bool do_cloud)
{
  for(auto f : cloud_handling_vector_)
    f->do_cloud(do_cloud);
}
