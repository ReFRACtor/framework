#include "radiance_scaling_sv_muses_fit.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void RadianceScalingSvMusesFit::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RadianceScaling)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayInstrumentCorrection);
}

FP_IMPLEMENT(RadianceScalingSvMusesFit);
#endif

void RadianceScalingSvMusesFit::apply_correction
(const SpectralDomain& Pixel_grid,
 const std::vector<int>& Pixel_list,
 SpectralRange& Radiance) const
{
  ArrayAd<double, 1> grid_ad(Pixel_list.size(), Pixel_grid.data_ad().number_variable());
  for (int i = 0; i < (int) Pixel_list.size(); i++) {
    grid_ad(i) = Pixel_grid.data_ad()(Pixel_list[i]);
  }
  grid_ad.resize_number_variable(coeff.number_variable());
  grid_ad.jacobian() = 0;
  SpectralDomain grid_sd(grid_ad, Pixel_grid.units());
  apply_scaling(grid_sd, Radiance);
}

boost::shared_ptr<InstrumentCorrection> RadianceScalingSvMusesFit::clone() const
{
  return boost::shared_ptr<InstrumentCorrection>
    (new RadianceScalingSvMusesFit(coeff.value().copy(), band_ref, band_name));
}

std::string RadianceScalingSvMusesFit::state_vector_name_i(int i) const
{ 
  return "resscale_O" + boost::lexical_cast<std::string>(i) + "_" + band_name;
}
 
void RadianceScalingSvMusesFit::print(std::ostream& Os) const
{
  Os << "RadianceScalingSvMusesFit:" << std::endl;
  OstreamPad opad(Os, "    ");
  RadianceScaling::print(opad);
  opad.strict_sync();
}
