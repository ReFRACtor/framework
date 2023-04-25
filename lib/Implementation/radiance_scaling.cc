#include "radiance_scaling.h"
#include "polynomial_eval.h"
#include "ostream_pad.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void RadianceScaling::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InstrumentCorrection)
    & FP_NVP(scaling_coeff) & FP_NVP(offset) & FP_NVP(band_ref)
    & FP_NVP(band_name);
}

FP_IMPLEMENT(RadianceScaling);
#endif

void RadianceScaling::apply_scaling(const SpectralDomain& Grid, SpectralRange& Radiance) const
{
  double band_ref_conv = band_ref.convert_wave(Grid.units()).value;
  if(Radiance.data_ad().number_variable() == 0) { 
    Poly1d scaling_poly(scaling_coeff.value(), false);

    for(int i = 0; i < Radiance.data().rows(); ++i) {
      double scale_factor = scaling_poly(Grid.data()(i) - band_ref_conv);
      Radiance.data()(i) = Radiance.data()(i) * scale_factor + offset.value();
    }

  } else {
    Poly1d scaling_poly(scaling_coeff, false);

    for(int i = 0; i < Radiance.data_ad().rows(); ++i) {
      AutoDerivative<double> scale_factor(scaling_poly(Grid.data_ad()(i) - band_ref_conv));
      Radiance.data_ad()(i) = Radiance.data_ad()(i) * scale_factor + offset;
    }
  }

}

void RadianceScaling::print(std::ostream& Os) const
{
  Poly1d scaling_poly(scaling_coeff.value(), false);
  Os << "RadianceScaling:" << std::endl;
  OstreamPad opad(Os, "    ");
  opad << "Band: " << band_name << std::endl
       << "Scaling polynomial:" << std::endl
       << "    " << scaling_poly << std::endl
       << "Offset: " << offset.value() << std::endl;
  opad.strict_sync();
}
