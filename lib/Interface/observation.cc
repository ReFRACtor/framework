#include "observation.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Observation::serialize(Archive & ar,
			const unsigned int version)
{
  FP_GENERIC_BASE(Observation);
  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
  // We use to be a StackedRadianceMixin, I think we can just remove
  // this
  //ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StackedRadianceMixin);
  // Older version wasn't an observable. If we don't load anything,
  // then none of the observable stuff is set up, which is actually
  // what we want to older serialization
  if(version > 0)
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableObservation);
    
}

FP_IMPLEMENT(Observation);
FP_OBSERVER_SERIALIZE(Observation);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Observation)
REGISTER_LUA_END()
#endif

// We basically have a StackedRadianceMixin but with the extra
// include_bad_sample argument.
// The code from StackedRadianceMixin isn't too complicated, so just
// copy here and add the extra argument

boost::optional<blitz::Range> Observation::stacked_pixel_range(int sensor_index, bool include_bad_sample) const
{
  range_check(sensor_index, 0, num_channels());
  int sind = 0;

  for(int i = 0; i < sensor_index; ++i)
    sind += spectral_domain(i,include_bad_sample).data().rows();
  
  int nrow = spectral_domain(sensor_index, include_bad_sample).data().rows();

  if(nrow > 0)
    return boost::optional<blitz::Range>(blitz::Range(sind, sind + nrow - 1));
  else
    return boost::optional<blitz::Range>();
}


Spectrum Observation::radiance_all(bool skip_jacobian,
				   bool include_bad_sample) const
{
  std::vector<Spectrum> sall;
  std::vector<blitz::Range> prall;
  for(int i = 0; i < num_channels(); ++i) {
    boost::optional<blitz::Range> pr = stacked_pixel_range(i, include_bad_sample);

    if(pr) {
      sall.push_back(radiance(i, skip_jacobian, include_bad_sample));
      prall.push_back(*pr);
    }
  }

  if(sall.size() == 0)
    throw Exception("Measured radiance Spectrum empty, pixel ranges must be empty");

  int nrow = 0;
  int nvar = 0;
  bool have_uncertainty = true;
  BOOST_FOREACH(const Spectrum & s, sall) {
    nrow += s.spectral_domain().data().rows();
    nvar = std::max(nvar, s.spectral_range().data_ad().number_variable());

    if(s.spectral_range().uncertainty().rows() == 0)
      have_uncertainty = false;
  }
  Unit ud, ur;
  if(sall.size() > 0) {
    ud = sall[0].spectral_domain().units();
    ur = sall[0].spectral_range().units();
  }

  ArrayAd<double, 1> sr(nrow, nvar);
  blitz::Array<double, 1> sd(nrow);
  blitz::Array<double, 1> uncer;

  if(have_uncertainty)
    uncer.resize(nrow);

  for(int i = 0; i < (int) sall.size(); ++i) {
    sd(prall[i]) = sall[i].spectral_domain().data() *
      FullPhysics::conversion(sall[i].spectral_domain().units(), ud);
    sr.value()(prall[i]) = sall[i].spectral_range().data() *
      FullPhysics::conversion(sall[i].spectral_range().units(), ur);
    if(have_uncertainty)
      uncer(prall[i]) = sall[i].spectral_range().uncertainty() *
	FullPhysics::conversion(sall[i].spectral_range().units(), ur);
    if(nvar > 0)
      sr.jacobian()(prall[i], blitz::Range::all()) =
	sall[i].spectral_range().data_ad().jacobian() *
	FullPhysics::conversion(sall[i].spectral_range().units(), ur);
  }

  return Spectrum(SpectralDomain(sd, ud), SpectralRange(sr, ur, uncer));
}
