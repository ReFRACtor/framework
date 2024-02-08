#include "spectral_window.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SpectralWindow::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(SpectralWindow);
  
  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
}

FP_IMPLEMENT(SpectralWindow);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(SpectralWindow)
.def("number_spectrometer", &SpectralWindow::number_spectrometer)
.def("spectral_bound", &SpectralWindow::spectral_bound)
.def("grid_indexes", &SpectralWindow::grid_indexes)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Apply a spectral window to a SpectralDomain, returning the
/// possibly empty part of the domain that passes through the window.
//-----------------------------------------------------------------------

SpectralDomain SpectralWindow::apply(const SpectralDomain& Grid, 
				     int Spec_index) const
{
  std::vector<int> gi = grid_indexes(Grid, Spec_index);
  Array<double, 1> res((int) gi.size());
  Array<int, 1> sind;
  if(Grid.sample_index().rows() > 0)
    sind.resize((int) gi.size());
  for(int i = 0; i < res.rows(); ++i) {
    res(i) = Grid.data()(gi[i]);
    if(sind.rows() > 0)
      sind(i) = Grid.sample_index()(gi[i]);
  }
  if(sind.rows() > 0)
    return SpectralDomain(res, sind, Grid.units());
  else
    return SpectralDomain(res, Grid.units());
}

//-----------------------------------------------------------------------
/// Apply a spectral window to a Spectrum, returning the
/// possibly empty part of the spectrum that passes through the window.
//-----------------------------------------------------------------------

Spectrum SpectralWindow::apply(const Spectrum& Spec, 
			       int Spec_index) const
{
  std::vector<int> gi = grid_indexes(Spec.spectral_domain(), Spec_index);
  Array<double, 1> res_d((int) gi.size());
  ArrayAd<double, 1> res_r((int) gi.size(), 
			   Spec.spectral_range().data_ad().number_variable());
  Array<double, 1> uncer;
  if(Spec.spectral_range().uncertainty().rows() > 0)
    uncer.resize(res_r.rows());
  Array<int, 1> sind;
  if(Spec.spectral_domain().sample_index().rows() > 0)
    sind.resize((int) gi.size());
  for(int i = 0; i < res_d.rows(); ++i) {
    res_d(i) = Spec.spectral_domain().data()(gi[i]);
    res_r(i) = Spec.spectral_range().data_ad()(gi[i]);
    if(uncer.rows() > 0)
      uncer(i) = Spec.spectral_range().uncertainty()(gi[i]);
    if(sind.rows() > 0)
      sind(i) = Spec.spectral_domain().sample_index()(gi[i]);
  }
  if(sind.rows() > 0)
    return Spectrum(SpectralDomain(res_d, sind, Spec.spectral_domain().units()),
		    SpectralRange(res_r, Spec.spectral_range().units(),
				  uncer));
  else
    return Spectrum(SpectralDomain(res_d, Spec.spectral_domain().units()),
		    SpectralRange(res_r, Spec.spectral_range().units(),
				  uncer));
}

