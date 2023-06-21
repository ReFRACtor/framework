#include "standard_forward_model.h"
#include "ostream_pad.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StandardForwardModel::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ForwardModel)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableNamedSpectrumPtr)
    & FP_NVP(spec_effect)
    & FP_NVP(inst)
    & FP_NVP(swin)
    & FP_NVP(rt)
    & FP_NVP_(spectrum_sampling)
    & FP_NVP(g);
}

FP_IMPLEMENT(StandardForwardModel);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(StandardForwardModel, ForwardModel)
.def(luabind::constructor <
     const boost::shared_ptr<Instrument>&,
     const boost::shared_ptr<SpectralWindow>&,
     const boost::shared_ptr<RadiativeTransfer>&,
     const boost::shared_ptr<SpectrumSampling>&,
     const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& > ())
.def("spectral_grid", &StandardForwardModel::spectral_grid)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

StandardForwardModel::StandardForwardModel(
    const boost::shared_ptr<Instrument>& Inst,
    const boost::shared_ptr<SpectralWindow>& Spectral_window,
    const boost::shared_ptr<RadiativeTransfer>& Rt,
    const boost::shared_ptr<SpectrumSampling>& Spectrum_sampling,
    const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& Spectrum_effect)
    : spec_effect(Spectrum_effect), inst(Inst), swin(Spectral_window), rt(Rt), spectrum_sampling_(Spectrum_sampling)
{
    if(spec_effect.size() == 0) {
        spec_effect.resize(num_channels());
    }

    if(num_channels() != (int) spec_effect.size()) {
        Exception err;
        err << "Spectrum effect size: " << spec_effect.size() 
            << " needs to be the same size as the number of spectrometers: " << num_channels();
        throw err;
                
    }
}

// See bass class for description.
Spectrum StandardForwardModel::radiance
(int sensor_index, bool Skip_jacobian) const
{
    if(!g) {
        throw Exception ("setup_grid needs to be called before calling radiance");
    }

    range_check(sensor_index, 0, num_channels());
    Spectrum highres_spec =
        rt->reflectance(g->high_resolution_grid(sensor_index), sensor_index,
                        Skip_jacobian);
    notify_spectrum_update(highres_spec, "high_res_rt", sensor_index);

    Spectrum convolved_spec = apply_spectrum_corrections(highres_spec, sensor_index);

    notify_spectrum_update(convolved_spec, "convolved", sensor_index);
   
    return convolved_spec;
}

//-----------------------------------------------------------------------
/// Applies corrections and modeling to modeled spectrum:
/// * Interpolation to uniform grid
/// * Application of spectrum effects
/// * Application of instrument model
//-----------------------------------------------------------------------

Spectrum StandardForwardModel::apply_spectrum_corrections(const Spectrum& highres_spec, int sensor_index) const
{
    if(!g) {
        throw Exception ("setup_grid needs to be called before calling apply_spectrum_corrections");
    }

    Spectrum highres_spec_intepolated =
        g->interpolate_spectrum(highres_spec, sensor_index);
    notify_spectrum_update(highres_spec_intepolated, "high_res_interpolated", sensor_index);

    BOOST_FOREACH(const boost::shared_ptr<SpectrumEffect>& i, spec_effect[sensor_index]) {
        if (i) {
            i->apply_effect(highres_spec_intepolated, *g);
            notify_spectrum_update(highres_spec_intepolated, "high_res_spec_effect_" + i->name(), sensor_index);
        }
    }

    return inst->apply_instrument_model(highres_spec_intepolated,
                                        g->pixel_list(sensor_index), sensor_index);
}

void StandardForwardModel::notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int sensor_index) const
{
  if (olist.size() > 0)
    notify_update_do(boost::make_shared<NamedSpectrum>(updated_spec, spec_name, sensor_index));
}

void StandardForwardModel::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "StandardForwardModel:\n";
    Os << "  Spectral Window:\n";
    opad << *swin << "\n";
    opad.strict_sync();
    Os << "  Spectrum Sampling:\n";
    opad << *spectrum_sampling_ << "\n";
    opad.strict_sync();
    Os << "  Instrument:\n";
    opad << *inst << "\n";
    opad.strict_sync();
    Os << "  Radiative Transfer:\n";
    opad << *rt << "\n";
    opad.strict_sync();

    for(int i = 0; i < num_channels(); ++i) {
        Os << "  Spectrum Effect[" << i << "]:\n";
        BOOST_FOREACH(boost::shared_ptr<SpectrumEffect> se, spec_effect[i]) {
          if(se) opad << *se << "\n";
        }
        opad.strict_sync();
    }
}

std::vector<boost::shared_ptr<GenericObject> >
StandardForwardModel::subobject_list() const
{ std::vector<boost::shared_ptr<GenericObject> > res;
  BOOST_FOREACH(auto i, spec_effect) {
    BOOST_FOREACH(auto j, i) {
      res.push_back(j);
    }
  }
  res.push_back(inst);
  res.push_back(swin);
  res.push_back(rt);
  res.push_back(spectrum_sampling_);
  res.push_back(g);
  return res;
}
  
