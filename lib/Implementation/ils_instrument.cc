#include "ils_instrument.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"
#include <boost/foreach.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void IlsInstrument::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Instrument)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverIls)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverInstrumentCorrection)
    & FP_NVP_(ils) & FP_NVP(inst_corr);
}

FP_IMPLEMENT(IlsInstrument);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(IlsInstrument, Instrument)
.def(luabind::constructor<const std::vector<boost::shared_ptr<Ils> >&>())
.def(luabind::constructor<const std::vector<boost::shared_ptr<Ils> >&,
                          const std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >& >())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

IlsInstrument::IlsInstrument(const std::vector<boost::shared_ptr<Ils> >& Ils_list,
                             const std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >&
                             Instrument_correction)
: ils_(Ils_list), inst_corr(Instrument_correction)
{
  if(inst_corr.size() == 0)
    inst_corr.resize(ils_.size());
  if(ils_.size() != inst_corr.size())
    throw Exception("Ils and Instrument_correction need to be the same size");
  BOOST_FOREACH(boost::shared_ptr<Ils> i, ils_)
    i->add_observer(*this);
  BOOST_FOREACH(std::vector<boost::shared_ptr<InstrumentCorrection> >& i, 
                inst_corr) {
    BOOST_FOREACH(boost::shared_ptr<InstrumentCorrection>& j, i)
      j->add_observer(*this);
  }
}

boost::shared_ptr<Instrument> IlsInstrument::clone() const
{
  std::vector<boost::shared_ptr<Ils> > ils_vec;
  BOOST_FOREACH(boost::shared_ptr<Ils> i, ils_)
    ils_vec.push_back(i->clone());
  std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > > 
    inst_corr_vec;
  BOOST_FOREACH(const std::vector<boost::shared_ptr<InstrumentCorrection> >& i, 
                inst_corr) {
    std::vector<boost::shared_ptr<InstrumentCorrection> > t;
    BOOST_FOREACH(const boost::shared_ptr<InstrumentCorrection>& j, i)
      t.push_back(j->clone());
    inst_corr_vec.push_back(t);
  }
  
  return boost::shared_ptr<Instrument>(new IlsInstrument(ils_vec, inst_corr_vec));
}

Spectrum IlsInstrument::apply_instrument_model(
    const Spectrum& High_resolution_spectrum,
    const std::vector<int>& Pixel_list,
    int Spec_index) const 
{
  range_check(Spec_index, 0, number_spectrometer());

  SpectralDomain full = pixel_spectral_domain(Spec_index);
  blitz::Array<double, 1> res_sd((int) Pixel_list.size());
  for(int i = 0; i < res_sd.rows(); ++i)
    res_sd(i) = full.data()(Pixel_list[i]);
  SpectralDomain res_dom(res_sd, full.units());
  // Pixel_list is relative to pixel_spectral_domain, so we prefer those units
  // This avoids the potential re-ordering associated with wavelength / wavenumber conversion
  // Convert High_resolution_spectrum to pixel_spectral_domain units
  // Re-order High_resolution_spectrum if needed
  boost::shared_ptr<SpectralDomain> hres_sd;
  boost::shared_ptr<SpectralRange> hres_sr =
    boost::make_shared<SpectralRange>(High_resolution_spectrum.spectral_range().data_ad(),
				      High_resolution_spectrum.spectral_range().units());
  if (High_resolution_spectrum.spectral_domain().units() != full.units())
    hres_sd = boost::make_shared<SpectralDomain>(High_resolution_spectrum.spectral_domain().convert_wave(full.units()), full.units());
  else
    hres_sd = boost::make_shared<SpectralDomain>(High_resolution_spectrum.spectral_domain().data(), full.units());
  // SpectralDomain expected to be ascending
  if (hres_sd->data()(0) > hres_sd->data()(1)) {
    // We copy the data because ils_grating assumes the data is
    // contiguous and in the same order as the ils
    hres_sd = boost::make_shared<SpectralDomain>
      (hres_sd->data_ad().reverse_and_copy(firstDim), hres_sd->units());
    hres_sr = boost::make_shared<SpectralRange>
      (hres_sr->data_ad().reverse_and_copy(firstDim), hres_sr->units());
  }
    
  // Gain some speed advantages when running things without autoderivatives
  // if they are not needed.
  SpectralRange res_sr;
  if(hres_sr->data_ad().number_variable() > 0) {
    ArrayAd<double, 1> rad_ad =
      ils_[Spec_index]->apply_ils(hres_sd->data(), hres_sr->data_ad(), Pixel_list);
    res_sr = SpectralRange(rad_ad, hres_sr->units());
  } else {
    Array<double, 1> rad =
      ils_[Spec_index]->apply_ils(hres_sd->data(), hres_sr->data(), Pixel_list);
    res_sr = SpectralRange(rad, hres_sr->units());
  }
  BOOST_FOREACH(const boost::shared_ptr<InstrumentCorrection>& i, inst_corr[Spec_index]) {
    if(i) {
        i->apply_correction(ils_[Spec_index]->pixel_grid(), Pixel_list, res_sr);
    }
  }

  return Spectrum(res_dom, res_sr);
}

void IlsInstrument::print(std::ostream& Os) const 
{ 
  Os << "IlsInstrument:\n";
  OstreamPad opad(Os, "    ");
  for(int i = 0; i < number_spectrometer(); ++i) {
    Os << "  Channel " << (i + 1) << ":\n";
    opad << *ils_[i] << "\n";
    BOOST_FOREACH(boost::shared_ptr<InstrumentCorrection> ic, inst_corr[i]) {
      if(ic) opad << *ic << "\n";
    }
    opad.strict_sync();
  }
}

std::vector<boost::shared_ptr<GenericObject> >
IlsInstrument::subobject_list() const
{
  std::vector<boost::shared_ptr<GenericObject> > res;
  BOOST_FOREACH(auto i, ils_)
    res.push_back(i);
  BOOST_FOREACH(auto i, inst_corr) {
    BOOST_FOREACH(auto j, i)
      res.push_back(j);
  }
  return res;
}
