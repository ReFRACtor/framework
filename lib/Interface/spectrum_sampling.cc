#include "spectrum_sampling.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SpectrumSampling::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(SpectrumSampling);
  ar & FP_NVP(nspectrometer);
}

FP_IMPLEMENT(SpectrumSampling);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(SpectrumSampling)
REGISTER_LUA_END()
#endif
