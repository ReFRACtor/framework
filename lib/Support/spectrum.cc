#include "spectrum.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Spectrum::serialize(Archive& ar,
			      const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(Spectrum);
  ar & FP_NVP_(spec_domain) & FP_NVP_(spec_range);
}

FP_IMPLEMENT(Spectrum);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
typedef const SpectralDomain& (Spectrum::*f1)(void) const;
typedef const SpectralRange& (Spectrum::*f2)(void) const;
REGISTER_LUA_CLASS(Spectrum)
.def("spectral_domain", ((f1) &Spectrum::spectral_domain))
.def("spectral_range", ((f2) &Spectrum::spectral_range))
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

Spectrum::Spectrum(const SpectralDomain& Spec_domain, 
		   const SpectralRange& Spec_range)
: spec_domain_(Spec_domain), spec_range_(Spec_range)
{
  if(spectral_domain().data().rows() !=
     spectral_range().data().rows()) {
    Exception e;
    e << "Spectral domain and range aren't the same size\n"
      << "  Spectral domain size: " << spectral_domain().data().rows() << "\n"
      << "  Spectral raunge size: " << spectral_range().data().rows() << "\n";
    throw e;
  }
}
