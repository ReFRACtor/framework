#include "forward_model_with_cloud_handling.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"
#include "logger.h"

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
(int sensor_index, bool skip_jacobian) const
{
  Logger::info() << "Computing Forward Model Clear Radiances: Sensor index" << (sensor_index+1) << "\n";
  set_do_cloud(false);
  Spectrum rclear = fmodel_->radiance(sensor_index, skip_jacobian);
  notify_spectrum_update(rclear, "clear", sensor_index);
  auto dclear = rclear.spectral_range().data_ad();

  Logger::info() << "Computing Forward Model Cloudy Radiances: Sensor index" << (sensor_index+1) << "\n";
  set_do_cloud(true);
  Spectrum rcloud = fmodel_->radiance(sensor_index, skip_jacobian);
  notify_spectrum_update(rcloud, "cloud", sensor_index);
  auto dcloud = rcloud.spectral_range().data_ad();
  set_do_cloud(false);

  ArrayAd<double, 1> dcfrac(dclear.rows(), std::max(
                            std::max(cfrac_->cloud_fraction().number_variable(),
                                     dclear.number_variable()),
                            dcloud.number_variable()));

  for(int i = 0; i < dcfrac.rows(); ++i)
    dcfrac(i) = dcloud(i) * cfrac_->cloud_fraction() +
      dclear(i) * (1 - cfrac_->cloud_fraction());

  Spectrum rcfrac(rclear.spectral_domain(),
                  SpectralRange(dcfrac, rclear.spectral_range().units()));
  notify_spectrum_update(rcfrac, "cloud fraction", sensor_index);

  return rcfrac;
}

void ForwardModelWithCloudHandling::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "FowardModelWithCloudHandling:\n"
     << "  cloud fraction:    " << cfrac_->cloud_fraction() << " \n"
     << "  underlying forward mode: \n";
  opad << *fmodel_ << "\n";
  opad.strict_sync();
}


void ForwardModelWithCloudHandling::set_do_cloud(bool do_cloud) const
{
  for(auto f : cloud_handling_vector_)
    f->do_cloud(do_cloud);
}

void ForwardModelWithCloudHandling::notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int sensor_index) const
{
  if (olist.size() > 0)
    notify_update_do(boost::make_shared<NamedSpectrum>(updated_spec, spec_name, sensor_index));
}
