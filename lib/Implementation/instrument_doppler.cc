#include "instrument_doppler.h"
#include "fp_serialize_support.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void InstrumentDoppler::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpectrumEffectImpBase)
    & FP_NVP(vel_units);
}

FP_IMPLEMENT(InstrumentDoppler);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
#include "state_mapping_at_indexes.h"

boost::shared_ptr<SpectrumEffect> instrument_doppler_create(const double Relative_velocity_value,
                                                            const std::string& Relative_velocity_units,
                                                            bool retrieved)
{
    blitz::Array<bool, 1> Flag(1);
    Flag(0) = retrieved;

    boost::shared_ptr<StateMapping> mapping =
        boost::make_shared<StateMappingAtIndexes>(Flag);

    boost::shared_ptr<InstrumentDoppler> inst_doppl =
        boost::make_shared<InstrumentDoppler>(Relative_velocity_value, Relative_velocity_units, mapping);

    return inst_doppl;
}

REGISTER_LUA_DERIVED_CLASS(InstrumentDoppler, SpectrumEffect)
.def(luabind::constructor<const DoubleWithUnit&>())
.def(luabind::constructor<const double, const std::string&>())
.scope
[
 luabind::def("create", &instrument_doppler_create)
]
REGISTER_LUA_END()
#endif
 
InstrumentDoppler::InstrumentDoppler(const DoubleWithUnit& Relative_velocity, boost::shared_ptr<StateMapping> Mapping)
  : vel_units(Relative_velocity.units)
{
  Array<double, 1> c(1);
  c = Relative_velocity.value;
  init(c, Mapping);
}

InstrumentDoppler::InstrumentDoppler(const double Relative_velocity_value,
                                     const std::string& Relative_velocity_units,
                                     boost::shared_ptr<StateMapping> Mapping) 
  : vel_units(Relative_velocity_units)
{
  Array<double, 1> c(1);
  c = Relative_velocity_value;
  init(c, Mapping);
}

void InstrumentDoppler::apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& UNUSED(Forward_model_grid)) const
{
  // Compute the wavelengths of the object in (2) in the spacecraft frame of reference.
  // Let v_ei be the relative velocity between earth and the instrument (=spacecraft).
  //
  //    a) Wavelengths: wl_spacecraft = wl_earth_rest_frame / (1.d0 + v_ei/c)
  //    b) Wavenumbers: wn_spacecraft = wn_earth_rest_frame * (1.d0 + v_ei/c)
  ArrayAd<double, 1> spec_dom_ad(Spec.spectral_domain().data_ad());
  ArrayAd<double, 1> doppler_val = mapping->fm_view(coefficient());
  if (spec_dom_ad.number_variable() == 0) {
    spec_dom_ad.resize_number_variable(doppler_val.number_variable());
    spec_dom_ad.jacobian() = 0;
  }
  for(int w_idx = 0; w_idx < Spec.spectral_domain().data().rows(); w_idx++) {
    if (Spec.spectral_domain().type_preference() == SpectralDomain::PREFER_WAVELENGTH)
      spec_dom_ad(w_idx) = spec_dom_ad(w_idx) 
        / (1.0 + doppler_val(0) / OldConstant::speed_of_light.value);
    else
      spec_dom_ad(w_idx) = spec_dom_ad(w_idx) 
        * (1.0 + doppler_val(0) / OldConstant::speed_of_light.value);
  }
}

boost::shared_ptr<SpectrumEffect> InstrumentDoppler::clone() const
{
  ArrayAd<double, 1> doppler_val = mapping->fm_view(coefficient());
  return boost::shared_ptr<SpectrumEffect>(new InstrumentDoppler(DoubleWithUnit(doppler_val(0).value(), vel_units)));
}

void InstrumentDoppler::print(std::ostream& Os) const
{
  Os << "InstrumentDoppler" << std::endl
     << "   Relative velocity: " << coefficient().value() << " " << vel_units << std::endl;
}
