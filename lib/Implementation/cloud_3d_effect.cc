#include "cloud_3d_effect.h"
#include "fp_serialize_support.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Cloud3dEffect::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpectrumEffectImpBase)
    & FP_NVP(band_name);
}

FP_IMPLEMENT(Cloud3dEffect);
#endif

Cloud3dEffect::Cloud3dEffect(const double Offset, const double Slope, const std::string& Band_name, boost::shared_ptr<StateMapping> Mapping)
: band_name(Band_name)
{
    Array<double, 1> coeffs(2);
    coeffs(0) = Offset;
    coeffs(1) = Slope;

    init(coeffs, Mapping);
}

const AutoDerivative<double> Cloud3dEffect::offset() const
{
    return coefficient()(0);
}

const AutoDerivative<double> Cloud3dEffect::slope() const
{
    return coefficient()(1);
}

void Cloud3dEffect::apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& UNUSED(Forward_model_grid)) const
{
    ArrayAd<double, 1> spec_range_ad(Spec.spectral_range().data_ad());
    ArrayAd<double, 1> mapped_coeffs = mapping->mapped_state(coefficient());
    AutoDerivative<double> offset(mapped_coeffs(0));
    AutoDerivative<double> slope(mapped_coeffs(1));

    for(int w_idx = 0; w_idx < spec_range_ad.rows(); w_idx++) { 
        spec_range_ad(w_idx) = spec_range_ad(w_idx) * (1 + offset + slope * spec_range_ad(w_idx));
    }
}

boost::shared_ptr<SpectrumEffect> Cloud3dEffect::clone() const
{
    return boost::shared_ptr<SpectrumEffect>(new Cloud3dEffect(offset().value(), slope().value(), band_name, mapping));
}

std::string Cloud3dEffect::state_vector_name_i(int i) const
{
    switch (i) {
        case 0:
            return "Cloud 3D Effect " + band_name + " Offset";
            break;
        case 1:
            return "Cloud 3D Effect " + band_name + " Slope";
            break;
        default:
            return "Cloud 3D Effect " + band_name + " Coefficient " + boost::lexical_cast<std::string>(i + 1);
    }
}

void Cloud3dEffect::print(std::ostream& Os) const
{
    Os << "Cloud3dEffect" << std::endl;
}
